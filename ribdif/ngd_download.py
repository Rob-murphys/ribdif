#!/usr/bin/env python3
import ncbi_genome_download as ngd
from pathlib import Path
import gzip
import shutil
import logging
# Using Kai's NCBI genome downloader to get all genones of the specified genus
# Avaliable at: https://github.com/kblin/ncbi-genome-download
def genome_download(genus, outdir, threads, frag, sp_ignore, domain, logger):
    genera = genus.replace("_", " ")
    if frag == "On": # if the user wants fragmented genomes
        logger.info(f"Downloading all genome recoreds of {genus} fron NCBI\n")
        ngd.download(section='refseq', 
               file_formats = 'fasta', 
               genera = genera,  
               output = outdir, 
               parallel = threads*2,
               groups = domain)
        
        
    else: # if the user (correctly) only want complete genomes
        logger.info(f"Downloading all complete genome records of {genus} from NCBI\n")
        ngd.download(section='refseq', 
               file_formats = 'fasta', 
               genera = genera,  
               output = outdir, 
               assembly_levels = "complete",
               parallel = threads*2,
               groups = domain)
    
    count = ""    
    try:
        if Path(f"{outdir}/refseq/{domain}").is_dir():
    
            count = len(list(Path(f"{outdir}/refseq/{domain}").glob("*/*.fna.gz")))
            dir_count = len(list(Path(f"{outdir}/refseq/{domain}").glob("*")))
            
            try:
                if count != dir_count:
                    raise FileNotFoundError()
            except FileNotFoundError as err:
                logger.error(str(err), exc_info = True)
                logger.error(f"{genus} is a real genus but some (or no) genomes were not downloaded")
                return 1
            if sp_ignore == "On":
                sp_count = sp_remove(outdir, domain)
                logger.info(f"{count} genomes of {genus} were downloaded and {sp_count} were removed due to being unnamed species\n\n")
            else: 
                logger.info(f"{count} genomes of {genus} were downloaded\n\n")
            return 0, count
        
        else:
            raise NotADirectoryError()
    except NotADirectoryError as err:
        logger.error(str(err), exc_info = True)
        logger.error(f"Download failed because {genus} is invalid or there are no records of the requested type in NCBI")
        return 1, count
    return 0, count

def sp_remove(outdir, domain):
    sp_count = 0
    downloads = list(Path(f"{outdir}/refseq/{domain}").rglob("*.fna.gz"))
    for download in downloads:
        splitline = sp_get(download)
        if splitline[2] == "sp.":
            shutil.rmtree(Path(download).parent)
            sp_count += 1
    return sp_count

def sp_get(download):
    with gzip.open(download, "r") as f_in:
        line = f_in.readline()
        splitline = str(line).strip().split(" ")
    return splitline