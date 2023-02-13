#!/usr/bin/env python3
import ncbi_genome_download as ngd
from pathlib import Path
import gzip

# Using Kai's NCBI genome downloader to get all genones of the specified genus
# Avaliable at: https://github.com/kblin/ncbi-genome-download
def genome_download(genus, outdir, threads, frag, sp_ignore):
    
    if frag: # if the user wants fragmented genomes
        print(f"Downloading all genome recoreds of {genus} fron NCBI\n")
        ngd.download(section='refseq', 
               file_formats = 'fasta', 
               genera = genus,  
               output = outdir, 
               parallel = threads*2,
               groups = 'bacteria')
        
        
    else: # if the user (correctly) only want complete genomes
        print(f"Downloading all complete genome records of {genus} from NCBI\n")
        ngd.download(section='refseq', 
               file_formats = 'fasta', 
               genera = genus,  
               output = outdir, 
               assembly_levels = "complete",
               parallel = threads*2,
               groups = 'bacteria')
        
    if Path(f"{outdir}/refseq/bacteria").is_dir():

        count = len(list(Path(f"{outdir}/refseq/bacteria").glob("*/*.fna.gz")))
        dir_count = len(list(Path(f"{outdir}/refseq/bacteria").glob("*")))
        
        if count != dir_count:
            print(count)
            print(dir_count)
            raise FileNotFoundError(repr(f"{genus} is a real genus but some (or no) genomes were not downloaded"))
        if sp_ignore:
            sp_count = sp_remove(outdir)
            print(f"\n\n{count} genomes of {genus} were downloaded and {sp_count} were removed due to being unnamed species")
        print(f"\n\n{count} genomes of {genus} were downloaded")
        return count
    
    else:
        raise NotADirectoryError(repr(f"Download failed because {genus} is invalid or there are no records of the requested type in NCBI\n\n"))
        return

def sp_remove(outdir):
    sp_count = 0
    downloads = list(Path(f"{outdir}/genbank/bacteria").rglob("*.fna.gz"))
    for download in downloads:
        with gzip.open(download, "r") as f_in:
            line = f_in.readline()
            splitline = line.strip().split("_")
            if splitline[4] == "sp.":
                Path.remdir(Path(download).parent)
                sp_count += 1
    