#!/usr/bin/env python3

from pathlib import Path
from dataclasses import dataclass
import argparse
import subprocess
import logging
import httpx
from tenacity import (
    before_sleep_log,
    retry,
    stop_after_delay,
    wait_exponential,
)

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

AVAILABLE_DBS = ['diamond', 'mmseqs', 'hmmer', 'no_search', 'cache']

BASE_UNVERSIONED_URL = 'https://downloads.eggnogdb.org/emapper/emapperdb-5.0.2/'
HMM_EGGNOG_URL = 'http://eggnog5.embl.de/download/eggnog_5.0/per_tax_level'
EGGNOG_DOWNLOADS_URL = 'http://eggnog5.embl.de/#/app/downloads'

NCBI_API_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy"
NCBI_API_HEADERS = {"accept": "application/json", "content-type": "application/json"}


#####################################################
# FUNCTIONS
#####################################################

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--db', dest="database", required=True, choices=AVAILABLE_DBS,help='Database to download')
    parser.add_argument('--db-version', dest="db_version", required=True, type=str, help='Version of the database to download')
    
    parser.add_argument("--out", dest="data_dir", required=True, type=Path, help='Directory to use for DATA_PATH.')
    parser.add_argument('--taxid', type=str,
                        help=(
                            'Tax ID of eggNOG HMM database to download. '
                            'e.g. "-H -d 2" for Bacteria. Required if "-H". '
                            f'Available tax IDs can be found at {EGGNOG_DOWNLOADS_URL}.'
                        ))
    parser.add_argument("--ncpus", required=True, type=int, help='Number of CPUs to use for downloading.')
    return parser.parse_args()
    

@retry(
    stop=stop_after_delay(600),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    before_sleep=before_sleep_log(logger, logging.WARNING),
)
def send_request_to_ncbi_taxonomy(taxid: str | int) -> dict:
    taxons = [str(taxid)]
    data = {"taxons": taxons}
    headers = dict(NCBI_API_HEADERS)
    response = httpx.post(NCBI_API_URL, headers=headers, json=data)
    response.raise_for_status()
    return response.json()


def get_taxon_name(taxid: str | int) -> str:
    try:
        response = send_request_to_ncbi_taxonomy(taxid)
        return response["taxonomy_nodes"][0]["taxonomy"]["organism_name"]
    except Exception as e:
        logger.error(f'Failed to get taxon name for tax ID {taxid}: {e}')
        return f'hmmer_{taxid}'


def run(cmd: list[str], shell: bool = False):
    str_cmd = " ".join(cmd)
    logger.info(f'Running command: {str_cmd}')
    subprocess.run(cmd, shell=shell, check=True)


def decompress(file: Path):
    cmd = ['pigz', '-df', str(file)]
    run(cmd)
    

def check_level_exists(taxid: str):
    cmd = ["wget", f"{HMM_EGGNOG_URL}/{taxid}/"]
    try:
        run(cmd)
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"Provided taxid {taxid} does not correspond to a valid HMM level. Please check available levels at {HMM_EGGNOG_URL}")


def create_hmm_database(level: str, db_path: Path):
    dbname = db_path.name
    cmd = [
        f'echo {db_path}/* | xargs mv -t ./ && rm -r {db_path} && ',
        f'rm {level}_hmms.tar.gz; ',
        'numf=$(find ./ | grep -c ".hmm$"); ',
        'curr=0; ',
        f'cat /dev/null > {dbname}.hmm_tmp; ',
        'for file in $(find ./ | grep ".hmm$"); do ',
        'curr=$((curr+1)); ',
        'echo "merging HMMs... ${file} (${curr}/${numf})"; ',
        f'cat "${{file}}" | sed -e "s/.faa.final_tree.fa//" -e "s/.faa.final_tree//" >> {dbname}.hmm_tmp; ',
        'rm "${file}"; ',
        'done; ',
        f'mv {dbname}.hmm_tmp {dbname}.hmm; '
        f'(if [ -f {dbname}.hmm.h3i ]; then rm {dbname}.hmm.h3*; fi) && '
        'echo "hmmpress-ing HMMs... " && '
        f'hmmpress {dbname}.hmm && '
        'echo "generating idmap file... " && '
        f'cat {dbname}.hmm | grep "^NAME" | sed -e "s/^NAME *//" | awk \'{{logger.info NR" "$0}}\' > {dbname}.hmm.idmap && '
        'echo "removing single OG hmm files... " && '
        f'echo ./*hmm | xargs rm; '
    ]
    run(cmd, shell=True)


def transform_alignment_to_fasta(level: str, dbpath: Path):
    cmd = [
        f'echo {dbpath}/* | xargs mv -t ./ && rm -rf {dbpath} && ',
        f'rm {level}_raw_algs.tar; ',
        'numf=$(find ./ | grep -c ".faa.gz$"); ',
        'curr=0; ',
        'for file in $(find ./ | grep ".faa.gz$"); do ',
        'curr=$((curr+1)); ',
        'echo "processing FASTAs...  ${file} (${curr}/${numf})"; ',
        'outf=$(echo "$file" | sed "s/\.raw_alg\.faa\.gz/\.fa/"); ',
        'zcat "$file" | awk \'/^>/{logger.info; next}{gsub("-", ""); logger.info}\' > "$outf" && ',
        'rm "$file"; ',
        'done'
    ]
    run(cmd, shell=True)


#####################################################
# DOWNLOADER CLASS
#####################################################

@dataclass
class Downloader:

    data_path: Path
    ncpus: int

    def download(self, url: str):
        cmd = [
            'aria2c', 
            '-s', str(self.ncpus), 
            '-x', str(self.ncpus),  
            '--optimize-concurrent-downloads', 
            '--check-integrity=true', 
            '--dir', 
            str(self.data_path), 
            url
        ]
        run(cmd)

        
    def download_and_decompress(self, url: str):
        self.download(url)
        file = data_path / url.split('/')[-1]
        decompress(file)


    def untar_decompress(self, file: Path):
        cmd = ['tar', '-xzf', str(file), '-C', str(self.data_path)]
        run(cmd)
        file.unlink()
    
    
    def download_annotations(self, base_url: str):
        url = base_url + '/eggnog.db.gz'
        self.download_and_decompress(url)

    
    def download_taxa(self, base_url: str):
        filename = 'eggnog.taxa.tar.gz'
        url = base_url + '/' + filename
        self.download(url)
        self.untar_decompress(data_path / filename)
    
  
    def download_diamond_db(self, base_url: str):
        url = base_url + '/eggnog_proteins.dmnd.gz'
        self.download_and_decompress(url)
    
    
    def download_mmseqs_db(self, base_url: str):
        url = base_url + '/mmseqs.tar.gz'
        self.download_and_decompress(url)
    
    
    def download_pfam_db(self, base_url: str):
        url = base_url + '/pfam.tar.gz'
        self.download_and_decompress(url)
    
    
    def download_hmm_database(self, level: str):
        baseurl = f'{HMM_EGGNOG_URL}/{level}/'
        check_level_exists(level)
        # Create HMMER database
        hmmsurl = f'{baseurl}/{level}_hmms.tar.gz'
        self.download_and_decompress(hmmsurl)
        create_hmm_database(level, db_path)
        # Transform alignment files to fasta files
        seqsurl = f'{baseurl}/{level}_raw_algs.tar'
        self.download_and_decompress(seqsurl)
        transform_alignment_to_fasta(level, db_path)
        

#####################################################
# MAIN
#####################################################

if __name__ == "__main__":
    args = parse_args()

    data_path = args.data_dir
    Path(data_path).mkdir(parents=True, exist_ok=True)

    if args.database == 'hmmer' and not args.taxid:
        raise ValueError('Must specify --taxid when downloading HMMER databases')

    base_url = BASE_UNVERSIONED_URL.format(args.db_version)

    downloader = Downloader(data_path, args.ncpus)
    
    # Annotation DB
    downloader.download_annotations(base_url)

    # NCBI taxa
    downloader.download_taxa(base_url)

    match args.database:
        case 'diamond':
            downloader.download_diamond_db(base_url)
        case 'mmseqs':
            logger.info(f'Downloading MMseqs2 files " at {data_path}...')
            downloader.download_mmseqs_db(base_url)
        case 'hmmer':
            taxon_name = get_taxon_name(args.taxid)
            db_path = data_path / 'hmmer' / taxon_name
            logger.info(f'Downloading HMMER database of tax ID {args.taxid} as "{taxon_name}" to {db_path}')
            logger.info('Note that this can take a long time for large taxonomic levels')
            downloader.download_hmm_database(args.taxid)
        case 'pfam':
            downloader.download_pfam_db(base_url)

    logger.info("Finished")
