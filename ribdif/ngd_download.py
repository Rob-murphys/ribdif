#!/usr/bin/env python3
import ncbi_genome_download as ngd
from pathlib import Path

# Using Kai's NCBI genome downloader to get all genones of the specified genus
# Avaliable at: https://github.com/kblin/ncbi-genome-download
def genome_download(genus, outdir, threads, frag):
    
    if frag: # if the user wants fragmented genomes
        print(f"Downloading all genome recoreds of {genus} fron NCBI")
        ngd.download(section='genbank', 
               file_formats = 'fasta', 
               genera = genus,  
               output = outdir, 
               parallel = threads*2,
               groups = 'bacteria')
        
        
    else: # if the user (correctly) only want complete genomes
        print(f"Downloading all complete genome records of {genus} from NCBI")
        ngd.download(section='genbank', 
               file_formats = 'fasta', 
               genera = genus,  
               output = outdir, 
               assembly_levels = "complete",
               parallel = threads*2,
               groups = 'bacteria')
    
    if Path(f"{outdir}/genbank/bacteria").is_dir():
        count = len(list(Path(f"{outdir}/genbank/bacteria").glob("*")))
        print(f"\n\n{count} genomes of {genus} were downloaded")
        return count
    else:
        raise Exception(f"Download failed because {genus} is invalid or there are no records of the requested type in NCBI\n\n")
        return