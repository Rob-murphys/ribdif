#!/usr/bin/env python3
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
from itertools import repeat
from hatch.cli import command




import ngd_download
import barrnap_run
import pcr_run
import pyani_run
import utils
import msa_run
import summary_16S
import vsearch_run
import overlaps
import make_heatmap



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
    
    parser.add_argument("-s", dest = "speed", 
                        help = "Run in fast mode skipping all but the nessacary steps (No barrnap, msa, ani etc.)",
                        action = "store_true")
    return parser.parse_args()

def arg_handling(args, workingDir):
    print("#= Parsing arguments =#\n\n")
    
    if args.speed:
        print("Running in fast mode skipping all but the nessacary steps (No barrnap, msa, ani etc.)")
    
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
    else:
        rerun = False
    
    # Resolving clobber argument
    if args.clobber == True:
        print(f"Removing old run of {genus}" )
        try:
            shutil.rmtree(Path(f"{outdir}")) # Remove genus and all subdirectories
        except FileNotFoundError: # catch if directory not found
            print(f"{genus} folder does not exist, ignoring clobber request")
            pass
    elif Path(f"{outdir}").is_dir() and args.rerun == False: # catch if genus output already exists and rerun was not specified and clobber was not used
       raise FileExistsError(f"""/n/n{outdir} folder already exists. Run again with -c/--clobber, -r/--rerun or set another output directory/n/n""")
       
    print("\n\n#= All arguments resolved =#\n\n")
    
    return genus, primer_file, outdir, rerun


def main():
    workingDir = Path(os.path.realpath(os.path.dirname(__file__))) # getting the path the script is running from
    
    args = parse_args()
    
    print(f"""\n
          #=========================================#
          #== RibDif2 is running on: {args.genus} ==#
          #=========================================#\n\n""")
    
    #Argument handeling
    genus, primer_file, outdir, rerun = arg_handling(args, workingDir)
    
    log_dir = Path(outdir) / "ribdif_logs"
    Path(log_dir).mkdir(exist_ok = True, parents = True)
    
    # If rerun is false, download and handle genomes from NCBI
    if rerun == False:
        # Download genomes from NCBI
        genome_count = ngd_download.genome_download(genus = genus, outdir = outdir, threads = args.threads, frag = args.frag)
    
        # Un gziping fasta files
        with multiprocessing.Pool(args.threads) as pool: # Create a multiprocessing pool with #threads workers
            all_gz = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('**/*.gz'))]# Recursively search the directory for .gz files and convert path to string sotring in a list
            pool.map(utils.decompress, all_gz)
        
        # Remove unwanted characters from anywhere is file (should only be in fasta headers)
        print("\n\nModifying fasta headers.\n\n")
        with multiprocessing.Pool(args.threads) as pool:
            all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('*/*.fna'))]
            pool.map(utils.modify, all_fna)
    else:
        genome_count = len(list(Path(f"{outdir}/genbank/bacteria").glob("*")))
        print(f"\n\n{genome_count} previously downloaded genomes of {genus} were found")
            
    # If using default primers call barrnap and rerun is false - this assume
    if args.primers == "False" and not args.speed: # Is skipping barrnap actually faster?
        print("#= Running barrnap on downloaded sequences =#\n\n")
        barrnap_run.barnap_call(outdir, threads = args.threads)
        
        # Processing barrnap output > fishing out 16S sequences
        with multiprocessing.Pool(args.threads) as pool:
            all_RNA = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('*/*.rRNA'))]
            #gene_num = [*range(len(all_RNA))] # adding gene num count here (in congruence with v1. Could also just use in silico pcr amp count)
            pool.map(barrnap_run.barrnap_process, all_RNA) # removed zip(gene_num) and was originally starmap
        
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
        
        if args.speed:
            pass
        else:
            # ALignment of full 16S genes recoverd from barrnap
            print("Alligning full-length 16S genes within genomes with muscle and building trees with fastree.\n\n")
            msa_run.muscle_call_multi(outdir, args.threads)
        
        # Genome statistic summary
        with open(f"{outdir}/{genus}_summary.tsv", "w") as f_out:
            f_out.write("GCF\tGenus\tSpecies\t#16S\tMean\tSD\tMin\tMax\tTotalDiv\n")
        summary_16S.multiproc_sumamry(outdir, genus, args.threads, args.speed)
        
        # Running msa on concatinated 16S sequences
        if args.msa == True:
            print(f"Alligning all {genus} 16S rRNA genes with muscle and building tree with fasttree.\n")
            infile , outAln, outTree = f"{outdir}/full/{genus}.16S", f"{outdir}/full/{genus}.16sAln", f"{outdir}/full/{genus}.16sAln" # Asigning in and out files
            msa_run.muscle_call_single(infile, outAln, outTree)
        else:
            print("Skipping alignments and tree generation for 16S rRNA genes (if needed, use -m/--msa).\n\n")
            
        # PCR for default primers
        infile = f"{outdir}/full/{genus}.16S" # path to concatinated 16S barrnap output
        names =  pcr_run.pcr_call(infile, outdir, genus, primer_file, workingDir)
        

        

    # PCR for custom primers   
    else:
        # Concatinate all downloaded genomes
        all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('*/*.fna'))]
        with open(f"{outdir}/genbank/bacteria/{genus}_total.fna", "w") as f_out:
            for file in all_fna:
                with open(file, "r") as f_in:
                    f_out.write(f_in.read())
                    
        
        infile = f"{outdir}/genbank/bacteria/{genus}_total.fna" # path to cocatinate genus genomes
        names = pcr_run.pcr_call(infile, outdir, genus, primer_file, workingDir)
        #pcr_run.pcr_parallel_call(outdir, genus, primer_file, workingDir, threads)
        #pcr_run.pcr_cleaner(outdir, primer_file, genus)
    
    for name in names:
        # Rename amplicon fasta headers to origin contig
        utils.amp_replace(outdir, genus, name)
        
        
    # msa on all amplicons
    if args.msa == True:
        print("Alligning all amplicons with Muscle and building tree with fasttree.\n")
        for name in names:
            infile , outAln, outTree = f"{outdir}/amplicons/{genus}-{name}.amplicons", f"{outdir}/amplicons/{genus}-{name}.aln", f"{outdir}/amplicons/{genus}-{name}.tree" # Asigning in and out files
            msa_run.muscle_call_single(infile, outAln, outTree)
    else:
        print("Skipping alignments and trees generation for amplicons (if needed, use -m/--msa).\n\n")
    
    print ("Making unique clusters with vsearch.\n\n")
    for name in names:
        vsearch_run.vsearch_call(outdir, genus, name, args.id, log_dir, args.threads)
    
    if args.msa == True:
        print("Making amplicon summary file for tree viewer import.\n\n")
        msa_run.format_trees(outdir, genus, name)
    
    # Generating the heatmaps #
    for name in names:
        print(f"Making heatmaps for {name}\n\n")
        # Cleaning vsearch clustering data
        all_gcfs, uc_dict_clean, gcf_species, cluster_count = overlaps.uc_cleaner(outdir, genus, name)
        
        # Generate a dictionary (that will become a matrix) of GCF cluster membership  
        cluster_dict = overlaps.cluster_matrix(all_gcfs, uc_dict_clean, cluster_count)
        
        # Find all species overlap in the cluster_dict
        combinations = overlaps.species_overlap(cluster_dict, cluster_count, gcf_species)
        
        # Find all GCF overlaps in the cluster dictionary
        pairwise_match = overlaps.gcf_overlaps(all_gcfs, uc_dict_clean, gcf_species)
        utils.pairwise_to_csv(pairwise_match, gcf_species, outdir, genus, name)
        
        # Generate metadata for heatmaps
        row_palette, species_series = make_heatmap.heatmap_meta(gcf_species)
        
        # Plot the cluster matrix
        plot_clus, cluster_df = make_heatmap.cluster_heatmap(cluster_dict, row_palette, species_series)
        
        # Plot the GCF overlap matrix
        plot_dendo, pairwise_df = make_heatmap.pairwise_heatmap(pairwise_match, row_palette, species_series)
        
        # Save the heatmaps
        make_heatmap.pdf_save(plot_clus, plot_dendo, outdir, genus, name)
        
        overlaps.overlap_report(combinations, gcf_species, cluster_df, genus, name, outdir)

if __name__ == '__main__':
    main()