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
from pathlib import Path
import shutil
import multiprocessing




import ngd_download
import barrnap_run
import pcr_run
import pyani_run
import utils
import msa_run
import summary_16S
import vsearch_run



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
    
    # Mutually exclusive group of rerun and clobber
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-r", dest = "rerun",
                        help = "Rerun on same of different primers. Avoids having to redownload the same genomes",
                        action = "store_true") # Action store true will default to false when argument is not present
    
    group.add_argument("-c", dest = "clobber",
                        help="Delete previous run if present in output directory", 
                        action = "store_true")
    
    
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
                        help = "Identity to cluster amplicons at if not using the deafult 100%%. e.g. .99. Does not cluster at the genome level, so beware",
                        default = 1)
    
    parser.add_argument("-t", dest = "threads",
                        help = "Number of threads to use. Default is all avaliable",
                        default = os.cpu_count())
    return parser.parse_args()

def arg_handling(args, workingDir):
    print("Parsing arguments:\n\n")
    
    #Checking if user is running on a species
    if " " in args.genus: # checking is there is a space in genus/species name
        print(f"Detected species {args.genus.split(' ')[1]}.\n\n")
        genus = args.genus.replace(" ", "_") # if so replacing it with an '_'
    else:
        genus = args.genus


    # Parsing the primers argument
    if args.primers == "False":
        primer_path = workingDir.parent
        primer_file = Path(f"{primer_path}/docs/default.primers")
    else:
        primer_file = args.primers
    # Check primer file exists and is not empty
    if Path(f"{primer_file}").is_file() and os.stat(f"{primer_file}").st_size == 0:
        raise Exception(f"Error: The provided primer file is empty.\n{primer_file}\nPlease provide a populated primer file")

    elif Path(f"{primer_file}").is_file() == False: # If it does not exist then raise this exception
        raise Exception(f"Error: {primer_file} does not exist")
    

    # Checking if user provided own outdir and if not setting to RibDif directory
    if args.outdir == "False":
        outdir = Path(f"{workingDir.parent}/results/{genus}")
    else:
        outdir = Path(args.outdir)
    
    #Resolving rerun argument
    if args.rerun == True:
        if Path(f"{outdir}").is_dir() == False:
            print(f"No records of {genus} exists. Ignoring rerun request and downloading genomes")
            rerun = False
        else:
            print(f"Reusing previously downloaded {genus} genomes")
            rerun = True
    
    # Resolving clobber argument
    if args.clobber == True:
        print(f"Removing old run of {genus}" )
        try:
            shutil.rmtree(Path(f"{outdir}")) # Remove genus and all subdirectories
        except FileNotFoundError: # catch if directory not found
            print(f"{genus} folder does not exist, ignoring clobber request")
            pass
    elif Path(f"{outdir}").is_dir() and args.rerun == False: # catch if genus output already exists and rerun was not specified and clobber was not used
       raise Exception(f"/n/n{outdir} folder already exists. Run again with -c/--clobber, -r/--rerun or set another output directory/n/n")
       
    print("\n\nAll arguments resolved\n\n")
    return genus, primer_file, outdir, rerun


def main():
    workingDir = Path(os.path.realpath(os.path.dirname(__file__))) # getting the path the script is running from

    args = parse_args()
    
    print(f"\n#== RibDif2 is running on: {args.genus} ==#\n\n")
    
    #Argument handeling
    genus, primer_file, outdir, rerun = arg_handling(args, workingDir)
    
    # If rerun is false, download and handle genomes from NCBI
    if rerun == False:
        # Download genomes from NCBI
        ngd_download.genome_download(genus = genus, outdir = outdir, threads = args.threads, frag = args.frag)
    
        # Un gziping fasta files
        with multiprocessing.Pool(args.threads) as pool: # Create a multiprocessing pool with #threads workers
            all_gz = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('**/*.gz'))]# Recursively search the directory for .gz files and convert path to string sotring in a list
            pool.map(utils.decompress, all_gz)
        
        # Remove unwanted characters from anywhere is file (should only be in fasta headers)
        print("\n\nModifying fasta headers.\n\n")
        with multiprocessing.Pool(args.threads) as pool:
            all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('*/*.fna'))]
            pool.map(utils.modify, all_fna)
            
            
    # If using default primers call barrnap and rerun is false - this assume
    if args.primers == "False":
        print("Running barrnap on downloaded sequences\n\n")
        barrnap_run.barnap_call(outdir, threads = args.threads)
        
        # Processing barrnap output > fishing out 16S sequences
        with multiprocessing.Pool(args.threads) as pool:
            all_RNA = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('*/*.rRNA'))]
            pool.map(barrnap_run.barrnap_process, all_RNA)
        
        # Concatinate all 16S to one file
        barrnap_run.barrnap_conc(genus, outdir)
        
        #If ANI is true calculate
        if args.ANI:
            
            # First need to split whole 16S sequences into seperate files
            with multiprocessing.Pool(args.threads) as pool:
                all_16S = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('*/*.16S'))]
                pool.map(barrnap_run.barrnap_split, all_16S)
               
            # Call pyani
            print("Calculating intra-genomic mismatches and ANI for each genome.\n\n")
            pyani_run.pyani_call(outdir, args.threads)
        else:
            print("Skipping detailed intra-genomic analysis and ANI (if needed, use -a/--ANI).\n\n")
        
        # ALignment of full 16S genes recoverd from barrnap
        print("Alligning full-length 16S genes within genomes with muscle and building trees with fastree.\n\n")
        msa_run.muscle_call_multi(outdir, args.threads)
        
        # Genome statistic summary
        with open(f"{outdir}/{genus}_summary.tsv", "w") as f_out:
            f_out.write("GCF\tGenus\tSpecies\t#16S\tMean\tSD\tMin\tMax\tTotalDiv\n")
        summary_16S.multiproc_sumamry(outdir, genus, args.threads)
        
        # Running msa on concatinated 16S sequences
        if args.msa == True:
            print("Alligning all 16S rRNA genes with Muscle and building tree with fasttree.\n")
            infile , outAln, outTree = f"{outdir}/full/{genus}.16S", f"{outdir}/full/{genus}.16sAln", f"{outdir}/full/{genus}.16sAln" # Asigning in and out files
            msa_run.muscle_call_single(infile, outAln, outTree)
        else:
            print("Skipping alignments and tree generation (if needed, use -m/--msa).\n\n")
            
        # PCR for default primers
        infile = f"{outdir}/full/{genus}.16S" # path to concatinated 16S barrnap output
        name =  pcr_run.pcr_call(infile, outdir, genus, primer_file, workingDir)
        
        # Rename amplicon fasta headers to origin contig
        utils.amp_replace(outdir, genus, name)
        

    # PCR for custom primers   
    else:
        # Concatinate all downloaded genomes
        all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('*/*.fna'))]
        with open(f"{outdir}/genbank/bacteria/{genus}_total.fna", "w") as f_out:
            for file in all_fna:
                with open(file, "r") as f_in:
                    f_out.write(f_in.read())
        infile = f"{outdir}/genbank/bacteria/{genus}_total.fna" # path to cocatinate genus genomes
        name = pcr_run.pcr_call(infile, outdir, genus, primer_file, workingDir)
        
        # Rename amplicon fasta headers to origin contig
        utils.amp_replace(outdir, genus, name)
        #pcr_run.pcr_parallel_call(outdir, genus, primer_file, workingDir, threads)
        #pcr_run.pcr_cleaner(outdir, primer_file, genus)
        
    # msa on all amplicons
    if args.msa == True:
        print("Alligning all ampliconss with Muscle and building tree with fasttree.\n")
        infile , outAln, outTree = f"{outdir}/amplicons/{genus}-{name}.amplicons", f"{outdir}/amplicons/{genus}-{name}.aln", f"{outdir}/amplicons/{genus}-{name}.tree" # Asigning in and out files
        msa_run.muscle_call_single(infile, outAln, outTree)
    else:
        print("Skipping alignments and trees generation (if needed, use -m/--msa).\n\n")
    
    print ("Making unique clusters with vsearch.\n\n")
    vsearch_run.vsearch_call(outdir, genus, name)


if __name__ == '__main__':
    main()