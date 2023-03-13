#!/usr/bin/env python3
import ncbi_genome_download as ngd
from pathlib import Path
import gzip
import shutil
import logging
# Using Kai's NCBI genome downloader to get all genones of the specified genus
# Avaliable at: https://github.com/kblin/ncbi-genome-download

def genome_download(genus, outdir, threads, frag, sp_ignore, domain, logger):
    genera = genus.replace("_", " ") # replace "_" with " "
    assembly_level = "all" if frag else "complete" # assign assembly level based on user input
    
    # Download genomes with specific domain and genus/species
    logger.info(f"Downloading all genome records of {genus} from NCBI at {assembly_level} assembly level\n")
    ngd.download(section='refseq', 
           file_formats = 'fasta', 
           genera = genera,  
           output = outdir, 
           assembly_levels = assembly_level,
           parallel = threads*2,
           groups = domain)
    status, count = download_checker(outdir, domain, logger, genus, sp_ignore)
    
    # if non bacteria domain and no complete genomes were downloaded, try again at chromosome level as it is very rare to have "complete" non bacteria genomes
    if domain != "bacteria" and status == 1 and assembly_level == "complete":
        logger.info(f"No complete genomes for {genus} were found so we will try again using 'chromosome' assembly level as this is what non bacteria genomes are usually added as")
        assembly_level = "chromosome"
        ngd.download(section='refseq', 
               file_formats = 'fasta', 
               genera = genera,  
               output = outdir, 
               assembly_levels = assembly_level,
               parallel = threads*2,
               groups = domain)
        
    # Remove genomes of unknown species if genomes were downloaded (status ==0)
    if sp_ignore and status == 0:
        sp_count = sp_remove(outdir, domain) # Remove genomes of unknown species
        logger.info(f"{count} genomes of {genus} were downloaded and {sp_count} were removed due to being unnamed species\n\n")
    elif status == 0: 
        logger.info(f"{count} genomes of {genus} were downloaded\n\n")
    
    
        
    return status

# Checking if the download worked 
def download_checker(outdir, domain, logger, genus, sp_ignore): 
    try:
        # Check if the directory exists, if it does not then nothing was downloaded so raise and catch the exception
        if Path(f"{outdir}/refseq/{domain}").is_dir():
            count = len(list(Path(f"{outdir}/refseq/{domain}").glob("*/*.fna.gz"))) # count the number of genomes downloaded
            dir_count = len(list(Path(f"{outdir}/refseq/{domain}").glob("*"))) # count the number of genomes that should have been downloaded based on dir count
            
            try:
                if count != dir_count: # check if expected and actual download count match and if not raise error
                    raise FileNotFoundError()
            except FileNotFoundError as err:
                logger.error(str(err), exc_info = True)
                logger.error(f"{genus} is a real genus but some genomes were not downloaded, please repeat and hopefully this wont happen again ")
                return 1, count
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