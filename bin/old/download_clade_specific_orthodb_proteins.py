#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.
# Adapted from script originally written by Tomas Bruna:
# https://github.com/tomasbruna/orthodb-clades/blob/3163c23f98417e20cddc5da42ce6e92ba0d18dc3/selectClade.py

"""
Get the species IDs corresponding to a species clade, optionnally excluding species / clades from given lists.
"""
                                 
import argparse
import csv
import sys
import logging
import gzip
import io
import httpx

from pathlib import Path
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, Iterator, ClassVar
from dataclasses import dataclass, field

from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_delay,
    wait_exponential,
    before_sleep_log
)

STOP_RETRY_AFTER_DELAY = 120

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logging.getLogger("httpx").setLevel(logging.WARNING)
logger = logging.getLogger(__name__)

PROTEIN_DB_URL = "https://data.orthodb.org/v12/download/odb_data_dump/odb12v2_aa_fasta.gz"

SPECIES_IDS_OUTFILE = "{clade}.species_ids.txt"
CLADE_PROTEIN_OUTFILE = "{clade}.odb12v2.faa"


# ---------------------------------------------------
# FUNCTIONS
# ---------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--clade', dest="clade", type=str, required=True, help='Clade to select proteins from')
    parser.add_argument('--levels', dest="levels_file", type=Path, required=True, help='Path to the levels OrthoDB file')
    parser.add_argument('--levels2species', dest="levels2species_file", type=Path, required=True, help='Path to the levels2species OrthoDB file')
    parser.add_argument('--species', dest="species_file", type=str, required=True, help='Path to the species OrthoDB file')
    parser.add_argument('--outfile', dest="outfile", type=str, required=True, help='Output file path')
    parser.add_argument('--exclude-clades', dest="excluded_clades", type=str, help='A comma separated set of clade names to exclude')
    parser.add_argument(
        '--exclude-species', dest="excluded_species",  type=str,
        help='A comma separated set of species names to exclude. When used, --species must be supplied.'
    )
    return parser.parse_args()


# ---------------------------------------------------
# SPECIES SELECTOR
# ---------------------------------------------------

@dataclass
class SpeciesSelector:

    clade: str
    levels_file: Path
    levels2species_file: Path
    species_file: Path
    excluded_clades: str | None
    excluded_species: str | None

    def get_species_ids(self) -> set[str]:

        clade_id = self.get_clade_id(self.clade)
        species_ids: set = self.parse_species_ids(clade_id)
    
        if self.excluded_clades:
            for clade_to_exclude in self.excluded_clades.split(','):
                excluded_clade_id = self.get_clade_id(clade_to_exclude)
                excluded_species = self.parse_species_ids(excluded_clade_id)
                species_ids.remove(excluded_species)
    
        if self.excluded_species:
            excluded_species = self.get_species_ids_from_names(self.excluded_species)
            species_ids -= excluded_species

        return species_ids
    

    def get_clade_id(self, clade: str):
        with open(self.levels_file) as fin:
            for row in csv.reader(fin, delimiter='\t'):
                if row[1] == clade:
                    return row[0]
        sys.exit(f'Error: Clade {clade} not found in {self.levels_file}')
    
    
    def parse_species_ids(self, clade_id: str) -> set:
        species_ids = set()
        with open(self.levels2species_file) as fin:
            for row in csv.reader(fin, delimiter='\t'):
                if clade_id in row[3].strip("{}").split(","):
                    species_ids.add(row[1])
        return species_ids
    
    
    def get_species_ids_from_names(self, species_names: str) -> set:
        species_names = {s.strip() for s in species_names.split(",")}
        excluded = set()
        with open(self.species_file) as fin:
            for row in csv.reader(fin, delimiter='\t'):
                if row[2] in species_names:
                    excluded.add(row[1])
        return excluded


# ---------------------------------------------------
# PARALLEL + STREAMING DOWNLOAD
# ---------------------------------------------------

class _ChunkStream(io.RawIOBase):
    """Present an iterator of byte chunks as a readable binary file."""

    def __init__(self, chunks: Iterator[bytes]):
        self._it = iter(chunks)
        self._buf = memoryview(b"")

    def readable(self) -> bool:
        return True

    def readinto(self, b) -> int:
        while not len(self._buf):
            try:
                self._buf = memoryview(next(self._it))
            except StopIteration:
                return 0
        n = min(len(b), len(self._buf))
        b[:n] = self._buf[:n]
        self._buf = self._buf[n:]
        return n


@dataclass
class ParallelStreamFastaDownloader:
    """
    Parallel, ordered, streaming download of a .fasta.gz with on-the-fly filtering.
     
    Uses httpx (shared thread-safe client, connection pooling) + threads for
    concurrent byte-range requests. Output is reassembled in order, so gzip
    and the FASTA parser see one continuous stream. Nothing but the filtered
    output touches disk.
    """

    TIMEOUT: ClassVar[int] = 60
    CHUNK_SIZE: ClassVar[int] = 8 * 1024 * 1024
    N_WORKERS: ClassVar[int] = 3 # starts to return responses with status 403 from N_WORKERS = 4 onwards
    
    url: str
    out_path: str
    species_ids: set[str]
    
    client: httpx.Client = field(init=False)
    _last_id: str | None = field(default=None, init=False)

    def __post_init__(self):
        self.client = httpx.Client(
            follow_redirects=True,  # httpx does not follow redirects by default
            headers={"Accept-Encoding": "identity"},  # never let httpx/server re-encode
            timeout=httpx.Timeout(self.TIMEOUT, connect=15),
            limits=httpx.Limits(max_connections=self.N_WORKERS, max_keepalive_connections=self.N_WORKERS),
        )
        
    
    def _probe(self, client: httpx.Client) -> tuple[int | None, str | None]:
        """
        Return (total_size, etag) if the server honors Range requests,
        else (None, None). Uses a streamed request so that a server that
        ignores Range doesn't make us download the whole body.
        """
        with client.stream("GET", self.url, headers={"Range": "bytes=0-0"}) as r:
            if r.status_code == 206:
                total = r.headers.get("Content-Range", "").rsplit("/", 1)[-1]
                if total.isdigit():
                    return int(total), r.headers.get("ETag")
        return None, None


    @retry(
        retry=retry_if_exception_type((httpx.HTTPError, OSError)),
        stop=stop_after_delay(STOP_RETRY_AFTER_DELAY),
        wait=wait_exponential(multiplier=1, min=1, max=30),
        before_sleep=before_sleep_log(logger, logging.WARNING),
    )
    def _fetch_range(self, client: httpx.Client, start: int, end: int, etag: str | None) -> bytes:
        expected = end - start + 1
        headers = {"Range": f"bytes={start}-{end}"}
        if etag and not etag.startswith("W/"):
            # If the file changed on the server, it replies 200 (full body)
            # instead of 206, which we treat as an error below.
            headers["If-Range"] = etag
    
        r = client.get(self.url, headers=headers)
        if r.status_code != 206:
            raise OSError(
                f"expected 206 Partial Content, got {r.status_code} "
                "(Range ignored or file changed)"
            )
        data = r.content
        if len(data) != expected:
            raise OSError(f"short read: got {len(data)} of {expected} bytes")
        return data
            
    
    def _iter_chunks(self, client: httpx.Client, size: int, etag: str | None) -> Iterator[bytes]:
        """Yield chunks strictly in order."""
        ranges = iter(
            (s, min(s + self.CHUNK_SIZE, size) - 1) for s in range(0, size, self.CHUNK_SIZE)
        )
        pending: deque = deque()
        prefetch = 2 * self.N_WORKERS
    
        with ThreadPoolExecutor(max_workers=self.N_WORKERS) as ex:
    
            def submit() -> None:
                r = next(ranges, None)
                if r is not None:
                    pending.append(ex.submit(self._fetch_range, client, *r, etag))
    
            try:
                for _ in range(prefetch):
                    submit()
                while pending:
                    data = pending.popleft().result()  # wait for the *next* chunk
                    submit()                            # keep the window full
                    yield data
            finally:
                for f in pending:
                    f.cancel()
    
    
    def _iter_single_stream(self, client: httpx.Client) -> Iterator[bytes]:
        """Fallback: plain sequential download when Range isn't supported."""
        with client.stream("GET", self.url) as r:
            r.raise_for_status()
            # iter_raw() yields the bytes exactly as sent (no content decoding).
            yield from r.iter_raw(self.CHUNK_SIZE)
    

    @staticmethod
    def _iter_fasta(handle) -> Iterator[tuple[str, list[str]]]:
        """Yield (header, sequence_lines); a record is emitted only once complete."""
        header, lines = None, []
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if header is not None:
                    yield header, lines
                header, lines = line, []
            elif header is not None and line:
                lines.append(line)
        if header is not None:
            yield header, lines


    def _keep(self, header: str):
        """
        Filter out sequences that do not correspond to the provided species IDs
        In OrthoDB v12, sequences have headers like:
            >100_0:000000	100_0
        """
        def _keep(self, header: str) -> bool:
            header_id = header.split()[1]
            if header_id != self._last_id:
                logger.info(f"Now at species ID {header_id}")
                self._last_id = header_id
            return header_id in self.species_ids
                
    
    def stream_filter_fasta_gz(self) -> tuple[int, int]:
        """
        Download a .fasta.gz with parallel range requests, decompress and filter
        it on the fly, and write kept records to `out_path` (gzipped if it ends
        in .gz).
    
        keep(header, sequence) -> bool decides whether a record is written.
        `header` includes the leading '>'; `sequence` is the unwrapped sequence.
    
        prefetch: chunks submitted ahead of the consumer (default workers * 2).
        Memory use is roughly prefetch * chunk_size.
    
        Returns (records_seen, records_kept).
        """
    
        nb_sequences = nb_kept = 0
        with self.client:
            size, etag = self._probe(self.client)
    
            if size is None:
                chunks = self._iter_single_stream(self.client)
            else:
                chunks = self._iter_chunks(self.client, size, etag)
    
            raw = io.BufferedReader(_ChunkStream(chunks), buffer_size=1024 * 1024)
            opener = gzip.open if self.out_path.endswith(".gz") else open
    
            try:
                # GzipFile transparently handles multi-member (concatenated/bgzip) files.
                with gzip.open(raw, "rt") as text, opener(self.out_path, "wt") as out:
                    for header, lines in self._iter_fasta(text):
                        
                        nb_sequences += 1
                        if nb_sequences % int(1e5) == 0:
                            logger.info(f"Parsed {nb_sequences} sequences")
                        
                        if self._keep(header):
                            nb_kept += 1
                            out.write(header + "\n")
                            out.write("\n".join(lines) + "\n")
            finally:
                raw.close()
                if hasattr(chunks, "close"):
                    chunks.close()  # stops the prefetcher if we exit early
    
        return nb_sequences, nb_kept
    


# ---------------------------------------------------
# MAIN
# ---------------------------------------------------

def main():
    args = parse_args()

    # GET SPECIES IDS
    species_selector = SpeciesSelector(
        args.clade, 
        args.levels_file,
        args.levels2species_file,
        args.species_file,
        args.excluded_clades,
        args.excluded_species
    )
    
    species_ids = species_selector.get_species_ids()
    logger.info(f"Found {len(species_ids)} species for clade {args.clade}")
    
    species_ids_outfile = SPECIES_IDS_OUTFILE.format(clade=args.clade)
    logger.info(f"Writing species ids to {species_ids_outfile}")
    with open(species_ids_outfile, "w") as fout:
        fout.writelines([f"{species_id}\n" for species_id in species_ids])

    # DOWNLOAD PART OF THE WHOLE PROTEIN DB THAT CORRESPONDS TO THESE SPECIES IDS
    protein_db_outfile = CLADE_PROTEIN_OUTFILE.format(clade=args.clade)
    logger.info(f"Downloading protein db to {protein_db_outfile}")
    downloader = ParallelStreamFastaDownloader(
        PROTEIN_DB_URL,
        protein_db_outfile,
        species_ids
    )
    
    nb_sequences, nb_kept_sequences = downloader.stream_filter_fasta_gz()
    logger.info(f"Saw {nb_sequences} sequences, kept {nb_kept_sequences}")

if __name__ == '__main__':
    main()
