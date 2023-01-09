# -*- coding: utf-8 -*-
"""
A program to analyse the usefulness of amplicon sequences

Created under GPL-3.0 license

@author: Robert Murphy
@email: robert.murphy@bio.ku.dk
"""
# Importing requierd libraries
import argparse
import os
import pandas as pd
from pathlib import Path
import ngd_download
import barrnap_run
import shutil
import gzip
import multiprocessing
import re



# Defining user arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="RibDif2 evaluates the differences in user defined amplicons  within a genus or species to indicate if that amplicon can differentiate between the taxanomic groups")

    parser.add_argument("-g", dest="genus", 
                        help="The genus you want to search within. E.g. 'Staphylococcus' OR 'Staphylococcus aurea' if wanting to use a species", 
                        required = True)
    
    parser.add_argument("-o", dest = "outdir",
                        help = "Output direcotry path. Default is current directory",
                        default = "False")
    
    parser.add_argument("-c", dest = "clobber",
                        help="Delete previous run if present in output directory", 
                        action = "store_true") # Action store true  will default to false when argument is not present
    
    parser.add_argument("-p", dest = "primers", 
                        help = "Path to custom primer-file, must be a tab seperated file with name, forward and reverse primers. See default.primers",
                        default = "False")
    
    parser.add_argument("-a", dest = "ANI", 
                        help = "ANI is off by default, turn on if you care about individual genomes. The $genus-summary.csv-file will only contain a list of genomes when off",
                        action = "store_true")
    
    parser.add_argument("-f", dest="frag", 
                        help = "Allow use of fragmented genomes. Full genomes are recomended/requierd for detecting all 16S-genes, use with caution. Off by default",
                        action = "store_true")
    
    parser.add_argument("-m", dest = "msa",
                        help = "make multiple sequence alignment and trees",
                        action = "store_true")
    
    parser.add_argument("-i", dest = "id",
                        help = "Identity to cluster amplicons at if not using the deafult 100%. e.g. .99. Does not cluster at the genome level, so beware",
                        default = 1)
    
    parser.add_argument("-t", dest = "threads",
                        help = "Number of threads to use. Default is all avaliable",
                        default = os.cpu_count())
    return parser.parse_args()

def arg_handling(args, workingDir):
    print("handeling arguments")
    #Checking if user is running on a species
    if " " in args.genus: # checking is there is a space in genus/species name
        print(f"Detected species {args.genus.split(' ')[1]}.\n\n")
        genus = args.genus.replace(" ", "_") # if so replacing it with an '_'
    else:
        genus = args.genus

    print("looking at primers")
    # Parsing the primers argument
    if args.primers == "False":
        primer_path = workingDir.parent
        primer_file = Path(f"{primer_path}/docs/default.primers")
    else:
        primer_file = args.primers
        
    print("reading primers")
    # Check primer file exists and is not empty
    primerDf = 0
    if Path(f"{primer_file}").is_file():
        try:
            primerDf = pd.read_table(primer_file, header = None)
        except pd.errors.EmptyDataError: # if it exists but is empty they exit with this error message
            print(f"Error: The provided primer file is empty\n{primer_file}")
        #print(primerDf)
    elif Path(f"{primer_file}").is_file() == False: # If it does not exist then raise this exception
        raise Exception(f"Error: {primer_file} does not exist")
    
    print("looking at outdir")
    # Checking if user provided own outdir and if not setting to RibDif directory
    if args.outdir == "False":
        outdir = Path(f"{workingDir.parent}/results/{genus}")
    else:
        outdir = Path(args.outdir)
    
    print("clobber resolution")
    # Resolving clobber argument
    if args.clobber == True:
        print(f"Removing old run of {genus}" )
        try:
            shutil.rmtree(Path(f"{outdir}")) # Remove genus and all subdirectories
        except FileNotFoundError: # catch if directory not found
            print(f"{genus} folder does not exist, ignoring clobber request/n/n")
            pass
    elif Path(f"{outdir}").is_dir(): # catch if genus output already exists and clobber was not used
       raise Exception(f"{outdir} folder already exists. Run again with -c/--clobber or set another output directory/n/n")
    
    return genus, primerDf, outdir

# Un gzip the downloaded NCBI genomes
def decompress(file_path):
    # Open the .gz file and decompress it
    with gzip.open(file_path, "rb") as f_in, open(file_path[:-3], "wb") as f_out:
        shutil.copyfileobj(f_in, f_out) # copy 

# Remove unwanted characters from anywhere is file (should only be in fasta headers)
def modify(file_path):
    # Open and read the contents of the file
    with open(file_path, "r") as f_in:
        contents = f_in.read()
        
    # Preform substitutions for unwanted characters
    contents = re.sub(r"[:,/()=#\x27]", "", contents)
    contents = re.sub(r"[: ]", "_", contents)

    # Write the modified contents back to file
    with open(file_path, "w") as f_out:
        f_out.write(contents)

# =============================================================================
# def modify2(file_path):
#     # Open and read the contents of the file
#     with open(file_path, "rb") as f_in:
#         contentsB = f_in.read()
#         
#     # Preform substitutions for unwanted characters
#     table = bytes().maketrans(b" ", b"_")
#     contentsS = contentsB.translate(None, b":,/()=#").translate(table).decode()
# 
#     # Write the modified contents back to file
#     with open(file_path, "w") as f_out:
#         f_out.write(contentsS)
# =============================================================================



def main():
    workingDir = Path(os.path.realpath(os.path.dirname(__file__))) # getting the path the script is running from

    args = parse_args()
    
    print(f"\n#== RibDif2 is running on: {args.genus} ==#\n\n")
    
    print("preparing to handel arguments")
    genus, primerDf, outdir = arg_handling(args, workingDir) # handeling arguments
    
    print(outdir)
    print(genus)
    Ncpu = os.cpu_count()
    print("downloading genomes")
    ngd_download.genome_download(genus = genus, outdir = outdir, threads = Ncpu, frag = args.frag) # downloading requierd genomes
    
    
    with multiprocessing.Pool(Ncpu) as pool: # Create a multiprocessing pool with Ncpu workers
        all_gz = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.gz'))]# Recursively search the directory for .gz files and convert path to string sotring in a list
        pool.map(decompress, all_gz)
    

    with multiprocessing.Pool(Ncpu) as pool:
        all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.fna'))]
        pool.map(modify, all_fna)
    
    barrnap_run.barnap_call(outdir)


if __name__ == '__main__':
    main()