#!/usr/bin/env python3

from pathlib import Path
from dataclasses import dataclass
import argparse
import subprocess
import logging

logging.basicConfig(
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.INFO
)
logger = logging.getLogger(__name__)

AVAILABLE_DBS = ['diamond', 'mmseqs', 'pfam']

BASE_URL = 'https://data.cgmlab.org/eggnog-mapper/emapper-3.0/data/'
GO_OBO_FALLBACK_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"

#####################################################
# FUNCTIONS
#####################################################

def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--db', dest="database", required=True, choices=AVAILABLE_DBS,help='Database to download')
    parser.add_argument('--db-version', dest="db_version", required=True, type=str, help='Version of the database to download')
    parser.add_argument("--out", dest="data_dir", required=True, type=Path, help='Directory to use for DATA_PATH.')
    parser.add_argument("--ncpus", required=True, type=int, help='Number of CPUs to use for downloading.')
    return parser.parse_args()

def run(cmd: list[str], shell: bool = False):
    str_cmd = " ".join(cmd)
    logger.info(f'Running command: {str_cmd}')
    subprocess.run(cmd, shell=shell, check=True)


def decompress(file: Path):
    cmd = ['pigz', '-df', str(file)]
    run(cmd)
    
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
        
    def untar_decompress(self, filename: str):
        file = self.data_path / filename
        cmd = ['tar', '-xzf', str(file), '-C', str(self.data_path)]
        run(cmd)
        file.unlink()
    
    def download_annotations(self):
        for filename in ["eggnog.db", "eggnog.db.taxids.bin", "eggnog.db.fieldpresence.bin"]:
            self.download(BASE_URL + filename)

    def download_go_obo(self):
        try:
            self.download(BASE_URL + 'go-basic.obo')
        except:
            self.download(GO_OBO_FALLBACK_URL)
            

    def download_taxa(self):
        for filename in ["eggnog.taxa.db", "eggnog.taxa.db.traverse.pkl"]:
            self.download(BASE_URL + filename)
    
    def download_diamond_db(self):
        self.download(BASE_URL + 'eggnog_proteins.dmnd')
    
    def download_mmseqs_db(self):
        filename = 'mmseqs.tar.gz'
        url = BASE_URL + filename
        self.download(url)
        self.untar_decompress(filename)
    
    def download_pfam_db(self):
        filename = 'mmseqs.tar.gz'
        url = BASE_URL + filename
        self.download(url)
        self.untar_decompress(filename)
        

#####################################################
# MAIN
#####################################################

if __name__ == "__main__":
    args = parse_args()

    data_path = args.data_dir
    Path(data_path).mkdir(parents=True, exist_ok=True)

    downloader = Downloader(data_path, args.ncpus)
    
    # Annotation DB
    downloader.download_annotations()

    # GO OBO
    downloader.download_go_obo()

    # NCBI taxa
    downloader.download_taxa()

    match args.database:
        case 'diamond':
            downloader.download_diamond_db()
        case 'mmseqs':
            downloader.download_mmseqs_db()
        case 'pfam':
            downloader.download_pfam_db()

    logger.info("Finished")
