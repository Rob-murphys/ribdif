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
import sys
import logging




from ribdif import ngd_download, barrnap_run, pcr_run, pyani_run, utils, msa_run, summary_files, vsearch_run, overlaps, figures, custom_exceptions, logging_config

# =============================================================================
# import barrnap_run
# import pcr_run
# import pyani_run
# import utils
# import msa_run
# import summary_files
# import vsearch_run
# import overlaps
# import figures
# import logging_config
# =============================================================================

 # Initialise the logging
logging_config.cofigure_logging()


# Defining user arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="RibDif2 evaluates the differences in user defined amplicons  within a genus or species to indicate if that amplicon can differentiate between the taxanomic groups")
    
    group1 = parser.add_mutually_exclusive_group()# Mutually exclusive group of rerun and clobber
    group2 = parser.add_mutually_exclusive_group()# Mutually exclusive group of ANI and whole genome mode (for now)
    parser.add_argument("-d", "--domain", dest = "domain",
                        help = "What domain of life does the genus belong to?",
                        choices = ["bacteria", "fungi", "viral", "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other", "invertebrate", "vertebrate_other"],
                        default = "bacteria")
    parser.add_argument("-g", "--genus", dest = "genus", 
                        help="The genus you want to search within. E.g. 'Staphylococcus' OR 'Staphylococcus aurea' if wanting to use a species", 
                        required = True)
    
    parser.add_argument("--ignore-sp", dest = "sp_ignore",
                        help = "Ignore genomes with unspecified species (i.e. their species is 'sp.')",
                        action = "store_true")
    
    parser.add_argument("-o", "--outdir", dest = "outdir",
                        help = "Output direcotry path. Default is current directory",
                        default = "False")
    
    group1.add_argument("-r", "--rerun", dest = "rerun",
                        help = "Rerun on same or different primers. Avoids having to redownload the same genomes. Mutually exclusive with --clobber",
                        action = "store_true") # Action store true will default to false when argument is not present
    
    group1.add_argument("-c", "--clobber", dest = "clobber",
                        help="Delete previous run if present in output directory. Mutually exclusive with --rerun", 
                        action = "store_true") 
    
    parser.add_argument("-p", "--primers", dest = "primers", 
                        help = "Path to custom primer-file, must be a tab seperated file with name, forward and reverse primers. See default.primers",
                        default = "False")
    
    group2.add_argument("-a", "--ani", dest = "ANI", 
                        help = "ANI is off by default, turn on if you care about individual genomes. The $genus-summary.csv-file will only contain a list of genomes when off. Mutually exclusive with --whole-genome",
                        action = "store_true")
    
    parser.add_argument("-f", "--frag", dest = "frag", 
                        help = "Allow use of fragmented genomes. Full genomes are recomended/requierd for detecting all 16S-genes, use with caution. Off by default",
                        action = "store_true")
    
    parser.add_argument("-m", "--msa", dest = "msa",
                        help = "make multiple sequence alignment and trees",
                        action = "store_true")
    
    parser.add_argument("-i", "--id", dest = "id",
                        help = "Identity to cluster amplicons at if not using the deafult 100%%. e.g. .99. Does not cluster at the genome level, so beware",
                        default = 1)
    
    parser.add_argument("-t", "--threads", dest = "threads",
                        help = "Number of threads to use. Default is all avaliable",
                        default = os.cpu_count())
    
    group2.add_argument("-w", "--whole-genome", dest = "whole", 
                        help = "Indicate the primers given are to be run on the whole genome so no barrnap, msa, ani (Required if your primers are non 16S). Mutually exclusive with --ani",
                        action = "store_true")
    return parser.parse_args()

def arg_handling(args, workingDir):
    
    #Checking if user is running on a species
    if " " in args.genus: # checking is there is a space in genus/species name
        logging.info(f"Detected species {args.genus.split(' ')[1]}.\n\n")
        genus = args.genus.replace(" ", "_") # if so replacing it with an '_'
    else:
        genus = args.genus

    # Checking if user provided own outdir and if not setting to root directory
    if args.outdir == "False":
        outdir = Path(f"{Path.cwd()}/results/{genus}")
        
    else:
        outdir = Path(f"{args.outdir}/{genus}")
    # Make the directory
    #Path.mkdir(outdir, parents = True)

    # Resolving clobber argument
    try:
        if args.clobber:
            shutil.rmtree(Path(f"{outdir}")) # Remove genus and all subdirectories
            logging.info(f"Removing old run of {genus}")  
        elif Path(f"{outdir}").is_dir() and args.rerun == False: # catch if genus output already exists and rerun was not specified and clobber was not used
            raise FileExistsError()
    except FileNotFoundError as err: # catch if directory not found
        logging.info(f"{genus} folder does not exist, ignoring clobber request\n")
        pass
    except FileExistsError as err:
        logging.error(str(err), exc_info = True)
        logging.error(f"{outdir} folder already exists. Run again with -c/--clobber, -r/--rerun or set another output directory")
    
    logging_config.replace_log_file(outdir)
    
    #Resolving rerun argument
    if args.rerun:
        if Path(f"{outdir}").is_dir() == False:
            print(f"No records of {genus} exists. Ignoring rerun request and downloading genomes\n")
            rerun = False
        else:
            print(f"Reusing previously downloaded {genus} genomes\n")
            rerun = True
    else:
        rerun = False
    
   
    logging.info("#= Parsing arguments =#\n\n")
    
    if not args.whole and args.primers != "False":
        logging.info("You are using custom primers on only the 16S genes as you didn't enable whole-genome mode. This is not a problem (if they are 16S primers), but we are just letting you know\n")
    elif args.whole and args.primers == "False":
        logging.info("You are running in whole-genome mode but using the default primers. This is not a problem but will 'skip' barrnap and other potentially useful mectrics scrapped from the whole 16S genes\n" )

    # Parsing the primers argument
    if args.primers == "False":
        try:
            if args.domain != "bacteria":
                raise custom_exceptions.IncompatiablityError()
            else:
                primer_file = Path(f"{workingDir}/default.primers")
        except custom_exceptions.IncompatiablityError as err:
            logging.exception(str(err)) 
            logging.error("You can't use the default 16S primers on non bacteria life.")
    else:
        primer_file = args.primers
        
    # Check primer file exists and is not empty
    try:
        if Path(f"{primer_file}").is_file() and os.stat(f"{primer_file}").st_size == 0:
            raise custom_exceptions.EmptyFileError()
        elif Path(f"{primer_file}").is_file() == False: # If it does not exist then raise this exception
            raise FileNotFoundError()          
    except custom_exceptions.EmptyFileError as err:
        logging.error("fThe provided primer file is empty.\n{primer_file}\nPlease provide a populated primer file")
    except FileNotFoundError as err:
        logging.error(f"{primer_file} does not exist")

       
    print("#= All arguments resolved =#\n\n")
    
    return genus, primer_file, outdir, rerun


def main():
    workingDir = Path(os.path.realpath(os.path.dirname(__file__))) # getting the path the script is running from
    
    args = parse_args()
    
    genus_line = f"#== RibDif2 is running on: {args.genus} in the {args.domain} domain ==#"
    block_line = f"#{'=' * (len(genus_line)-2)}#"
    print(f"\n{block_line}\n{genus_line}\n{block_line}\n\n")
    
    #Argument handeling
    genus, primer_file, outdir, rerun = arg_handling(args, workingDir)
    
    log_dir = Path(outdir) / "ribdif_logs"
    Path(log_dir).mkdir(exist_ok = True, parents = True)
    
    # If rerun is false, download and handle genomes from NCBI
    if rerun == False:
        # Download genomes from NCBI
        genome_count = ngd_download.genome_download(genus, outdir, args.threads, args.frag, args.sp_ignore, args.domain)
    
        # Un gziping fasta files
        with multiprocessing.Pool(args.threads) as pool: # Create a multiprocessing pool with #threads workers
            all_gz = [str(i) for i in list(Path(f"{outdir}/refseq/bacteria/").glob('**/*.gz'))]# Recursively search the directory for .gz files and convert path to string sotring in a list
            pool.map(utils.decompress, all_gz)
        
            
        # Remove unwanted characters from anywhere is file (should only be in fasta headers)
        print("Modifying fasta headers.\n\n")
        with multiprocessing.Pool(args.threads) as pool:
            all_fna = [str(i) for i in list(Path(f"{outdir}/refseq/bacteria/").glob('*/*.fna'))]
            pool.map(utils.modify2, all_fna)
    else:
        genome_count = len(list(Path(f"{outdir}/refseq/bacteria").glob("*")))
        print(f"{genome_count} previously downloaded genomes of {genus} were found\n\n")
            
    # If not using whole-genome mode assume the primers being used are 16S (which they are if default)
    if not args.whole:
        print("#= Running barrnap on downloaded sequences =#\n\n")
        barrnap_run.barnap_call(outdir, threads = args.threads)
        
        # Processing barrnap output > fishing out 16S sequences
        with multiprocessing.Pool(args.threads) as pool:
            all_RNA = [str(i) for i in list(Path(f"{outdir}/refseq/bacteria/").glob('*/*.rRNA'))]
            #gene_num = [*range(len(all_RNA))] # adding gene num count here (in congruence with v1. Could also just use in silico pcr amp count)
            pool.map(barrnap_run.barrnap_process, all_RNA) # removed zip(gene_num) and was originally starmap
        
        # Concatinate all 16S to one file
        barrnap_run.barrnap_conc(genus, outdir)
        
        #If ANI is true calculate
        if args.ANI:
            # First need to split whole 16S sequences into seperate files
            with multiprocessing.Pool(args.threads) as pool:
                all_16S = [str(i) for i in list(Path(f"{outdir}/refseq/bacteria/").glob('*/*.16S'))]
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
        #all_aln = [str(i) for i in list(Path(f"{outdir}/refseq/bacteria/").glob('*/*.16s'))]
        #with open(f"{outdir}/{genus}_{summary_type}_summary.tsv", "w") as f_out:
            #f_out.write("GCF\tGenus\tSpecies\t#16S\tMean\tSD\tMin\tMax\tTotalDiv\n")
        summary_type = "16S"
        in_fna = f"{outdir}/full/{genus}.16S"
        summary_files.make_sumamry(in_fna, outdir, genus, args.whole, args.ANI, args.threads, summary_type)
        
        # Running msa on concatinated 16S sequences
        if args.msa == True:
            print(f"Alligning all {genus} 16S rRNA genes with muscle and building tree with fasttree.\n")
            infile , outAln, outTree = f"{outdir}/full/{genus}.16S", f"{outdir}/full/{genus}.16sAln", f"{outdir}/full/{genus}.16sTree" # Asigning in and out files
            msa_run.muscle_call_single(infile, outAln, outTree)
        else:
            print("Skipping alignments and tree generation for 16S rRNA genes (if needed, use -m/--msa).\n\n")
            
        # PCR for default primers
        infile = f"{outdir}/full/{genus}.16S" # path to concatinated 16S barrnap output
        names =  pcr_run.pcr_call(infile, outdir, genus, primer_file, workingDir)
        

        

    # PCR for custom primers   
    elif args.whole:
        # Concatinate all downloaded genomes
# =============================================================================
#         all_fna = [str(i) for i in list(Path(f"{outdir}/refseq/bacteria/").glob('*/*.fna'))]
#         with open(f"{outdir}/refseq/bacteria/{genus}_total.fna", "w") as f_out:
#             for file in all_fna:
#                 with open(file, "r") as f_in:
#                     f_out.write(f_in.read())
# =============================================================================
                    
        
        #infile = f"{outdir}/refseq/bacteria/{genus}_total.fna" # path to cocatinate genus genomes
        #names = pcr_run.pcr_call(infile, outdir, genus, primer_file, workingDir)
        names = pcr_run.pcr_parallel_call(outdir, genus, primer_file, workingDir, args.threads)
        
        for name in names:
            total_sum_dict = pcr_run.multi_cleaner(outdir, name)
            pcr_run.amplicon_cat(outdir, genus, name)
            pcr_run.sum_dict_write(outdir, genus, name, total_sum_dict)

        
    for name in names:
        # Rename amplicon fasta headers to origin contig
        utils.amp_replace(outdir, genus, name)
        
    # Catching if all amplification failed (empty lists evaluate to false)
    if not list(Path(f"{outdir}/amplicons/").glob(f"{genus}-*.amplicons")):
        sys.exit("No amplification for any of the given primers was successfull. Try again with different primers")
    
    for name in names:
        # Make summary file for whole genome mode (has to be after utils.amp_replace so cant have in main args.whole section)
        if args.whole:
            summary_type = f"{name}-amp"
            in_fna = f"{outdir}/amplicons/{genus}-{name}.amplicons"
            summary_files.make_sumamry(in_fna, outdir, genus, args.whole, args.ANI, args.threads, summary_type)

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
        row_palette, species_series, species_palette = figures.heatmap_meta(gcf_species)
        
        # Plot the cluster matrix
        plot_clus, cluster_df = figures.cluster_heatmap(cluster_dict, row_palette, species_series)
        
        # Plot the GCF overlap matrix
        plot_dendo, pairwise_df = figures.pairwise_heatmap(pairwise_match, row_palette, species_series)
        
        plot_clus = figures.figure_fix(plot_clus)
        plot_dendo = figures.figure_fix(plot_dendo)
        
        # Save the heatmaps
        Path.mkdir(f"{outdir}/figures")
        figures.pdf_save(plot_clus, plot_dendo, outdir, genus, name)
        
        # Generate graps from the pairwise dataframe
        adjacency_df = figures.create_adjacency(pairwise_df, cluster_df)
        graph_subs, n_subplots = figures.create_graph(adjacency_df)
        
        # Draw the generated graps into on plot
        figures.draw_graphs(graph_subs, n_subplots, species_palette, row_palette, outdir, genus, name)
        
        overlaps.overlap_report(combinations, gcf_species, cluster_df, genus, name, outdir)

if __name__ == '__main__':
    main()