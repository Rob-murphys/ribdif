# -*- coding: utf-8 -*-
import ncbi_genome_download as ngd
import os

# Using Kai's NCBI genome downloader to get all genones of the specified genus
# Avaliable at: https://github.com/kblin/ncbi-genome-download
def genome_download(genus, outdir, threads, frag):
    
    if frag: # if the user wants fragmented genomes
        ngd.download(section='genbank', 
               file_formats = 'fasta', 
               genera = genus,  
               output = outdir, 
               parallel = threads*2,
               groups = 'bacteria')
        
    else: # if the user (correctly) only want complete genomes
        ngd.download(section='genbank', 
               file_formats = 'fasta', 
               genera = genus,  
               output = outdir, 
               assembly_levels = "complete",
               parallel = threads*2,
               groups = 'bacteria')


    


