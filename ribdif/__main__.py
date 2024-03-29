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
import importlib.metadata
from itertools import repeat



import ribdif
from ribdif import ngd_download, barrnap_run, pcr_run, pyani_run, utils, msa_run, summary_files, vsearch_run, logging_config
from ribdif.custom_exceptions import EmptyFileError, IncompatiablityError, StopError, IncorrectFormatError

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


# Defining user arguments
def parse_args():
    parser = argparse.ArgumentParser(
        description="RibDif2 evaluates the differences in user defined amplicons  within a genus or species to indicate if that amplicon can differentiate between the taxanomic groups")
    
    group1 = parser.add_mutually_exclusive_group()# Mutually exclusive group of rerun and clobber
    group2 = parser.add_mutually_exclusive_group()# Mutually exclusive group of ANI and whole genome mode (for now)
    group3 = parser.add_mutually_exclusive_group(required = True) # Mutually exclusive group for genus and user
    
    doms = ["bacteria", "archaea", "fungi", "viral", "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other", "invertebrate", "vertebrate_other"]
    parser.add_argument("-d", "--domain", dest = "domain",
                        help = f"What domain of life does the genus belong to (default is bacteria)? Choices are: {', '.join(doms)}",
                        choices = doms,
                        default = "bacteria",
                        metavar = '')
    group3.add_argument("-g", "--genus", dest = "genus", 
                        help="The genus you want to search within. E.g. 'Staphylococcus' OR 'Staphylococcus aureas' if wanting to use a species. Mutually exclusive with --user")
    
    group3.add_argument("-u", "--user", dest="user",
                        help = "Directory containing user defined set of genomes to run RibDif on. Give full path. Files canbe gzipped. Mutually exclusive with --genus")
    
    parser.add_argument("--ignore-sp", dest = "sp_ignore",
                        help = "Ignore genomes with unspecified species (i.e. their species is 'sp.')",
                        action = "store_true")
    
    parser.add_argument("-o", "--outdir", dest = "outdir",
                        help = "Output directory path. Default is current directory",
                        default = "False")
    
    group1.add_argument("-r", "--rerun", dest = "rerun",
                        help = "Re-run on same or different primers. Avoids having to redownload the same genomes or recopy to a new direcotry (if using --user). Mutually exclusive with --clobber",
                        action = "store_true") # Action store true will default to false when argument is not present
    
    group1.add_argument("-c", "--clobber", dest = "clobber",
                        help="Delete previous run if present in output directory. Mutually exclusive with --rerun", 
                        action = "store_true") 
    
    parser.add_argument("-p", "--primers", dest = "primers", 
                        help = "Path to custom primer-file, must be a tab separated file with name, forward and reverse primers. See default.primers",
                        default = "False")
    
    group2.add_argument("-a", "--ani", dest = "ANI", 
                        help = "ANI is off by default, turn on if you care about individual genomes. The summary.tsv files will only contain a list of genomes when off. Mutually exclusive with --whole-genome",
                        action = "store_true")
    
    parser.add_argument("-f", "--frag", dest = "frag", 
                        help = "Allow use of fragmented genomes. Full genomes are recommended/required for detecting all 16S-genes, use with caution. Off by default",
                        action = "store_true")
    
    parser.add_argument("-m", "--msa", dest = "msa",
                        help = "Make multiple sequence alignment and trees of amplicons (also whole 16S genes if not using --whole-genome)",
                        action = "store_true")
    
    parser.add_argument("-i", "--id", dest = "id",
                        help = "Identity to cluster amplicons at if not using the default 1.0. e.g. .99. Does not cluster at the genome level, so beware.",
                        default = 1)
    
    parser.add_argument("-t", "--threads", dest = "threads",
                        help = "Number of threads to use. Default is all available",
                        default = os.cpu_count(),
                        type = int)
    
    group2.add_argument("-w", "--whole-genome", dest = "whole", 
                        help = "Indicate the primers given are to be run on the whole genome so no barrnap or ani (Required if your primers are non 16S). Mutually exclusive with --ani",
                        action = "store_true")
    return parser.parse_args()

def arg_handling(args, workingDir, logger):
    
    logger.info("#= Parsing arguments =#\n\n")
    
    if not args.whole and args.domain != "bacteria":
        logger.error("If searching non bacteria life you must enable --whole-genome (-w) mode")
        return 1
    
    if args.genus:
        # Checking if user is running on a species
        if " " in args.genus: # checking is there is a space in genus/species name
            logger.info(f"Detected species {args.genus.split(' ')[1]}.\n\n")
            genus = args.genus.replace(" ", "_") # if so replacing it with an '_'
        else:
            genus = args.genus
    else: # If no genus if given (and thus args.user was) set genus to parent directory
        genus = Path(args.user).name
        
    # Checking user provided directory exists   
    if args.user and not Path(args.user).is_dir():
        logger.info(f"{args.user} is not a valid directory. Please check the path and try again")
        return 1
    
    # Checking if user provided own outdir and if not setting to called from directory
    if args.outdir == "False":
        outdir = Path(f"{Path.cwd()}/results/{genus}")
    else:
        outdir = Path(f"{args.outdir}/{genus}")

    # Resolving clobber argument
    try:
        if args.clobber:
            shutil.rmtree(Path(f"{outdir}")) # Remove genus and all subdirectories
            logger.info(f"Removing old run of {genus}")  
        elif Path(f"{outdir}").is_dir() and args.rerun == False: # catch if genus output already exists and rerun was not specified and clobber was not used
            raise FileExistsError()
    except FileNotFoundError: # catch if directory not found
        logger.info(f"{genus} folder does not exist, ignoring clobber request\n")
        pass
    except FileExistsError as err:
        logger.error(str(err), exc_info = True)
        logger.error(f"{outdir} folder already exists. Run again with -c/--clobber, -r/--rerun or set another output directory")
        return 1
    
    # Make the outdir
    Path.mkdir(Path(outdir), parents = True, exist_ok = True)
    
    # Replace logging file with one in outdir
    logging_config.replace_log_file(outdir, logger)
    
    #Resolving rerun argument
    if args.rerun:
        if Path(f"{outdir}").is_dir() == False:
            logger.info(f"No records of {genus} exists. Ignoring rerun request and downloading genomes\n")
            rerun = False
        else:
            logger.info(f"Reusing previously downloaded {genus} genomes\n")
            rerun = True
    else:
        rerun = False
    
    
    # Parsing the primers argument
    if args.primers == "False":
        try:
            if args.domain != "bacteria":
                raise IncompatiablityError()
            else:
                primer_file = Path(f"{workingDir}/default.primers")
        except IncompatiablityError as err:
            logger.exception(str(err)) 
            logger.error("You can't use the default 16S primers on non bacteria life.")
            return 1
    else:
        primer_file = args.primers
        
    # Check primer file exists and is not empty and is in correct format
    try:
        if Path(f"{primer_file}").is_file() == False: # Check if it exists
            raise FileNotFoundError() 
        elif os.stat(f"{primer_file}").st_size == 0: # Check if it is empty
            raise EmptyFileError()
        elif utils.primer_check(primer_file, logger): # Check it is in correct format
            raise IncorrectFormatError()
    # Catch the exceptions
    except EmptyFileError as err:
        logger.error(str(err), exc_info = True)
        logger.error(f"The provided primer file is empty.\n{primer_file}\nPlease provide a populated primer file")
        return 1
    except FileNotFoundError as err:
        logger.error(str(err), exc_info = True)
        logger.error(f"{primer_file} does not exist")
        return 1
    except IncorrectFormatError as err:
        logger.error(str(err), exc_info = True)
        return 1
    
    # Warning about primer usage
    if not args.whole and args.primers != "False":
        logger.info("You are using custom primers on only the 16S genes as you didn't enable whole-genome mode. This is not a problem (if they are 16S primers), but we are just letting you know\n")
    elif args.whole and args.primers == "False":
        logger.info("You are running in whole-genome mode but using the default primers. This is not a problem but will 'skip' barrnap and other potentially useful mectrics scrapped from the whole 16S genes\n" )

       
    logger.info("#= All arguments resolved =#\n\n")
    
    return genus, primer_file, outdir, rerun


def main():
    workingDir = Path(os.path.realpath(os.path.dirname(__file__))) # getting the path the script is running from
    
    args = parse_args()
    
    # Initialise the logging
    logger = logging_config.configure_logging(workingDir)
    logger.info(f"RibDif v{ribdif.__version__}\n")
    input_arg = sys.argv
    input_arg[0] = input_arg[0].split("/")[-1]
    logger.info("Input command:")
    logger.info(" ".join(input_arg))
    if args.genus:
        genus_line = f"#== RibDif2 is running on: {args.genus} in the {args.domain} domain ==#"
    elif args.user:
        genus_line = f"#== RibDif2 is running on user defined genomes at: {args.user} ==#"
    block_line = f"#{'=' * (len(genus_line)-2)}#"
    logger.info(f"\n{block_line}\n{genus_line}\n{block_line}\n\n")
    
    #Argument handeling
    try:
        genus, primer_file, outdir, rerun = arg_handling(args, workingDir, logger)
    except TypeError: # Catching the type error from "TypeError: cannot unpack non-iterable int object"
        sys.exit(1)
    
    log_dir = Path(outdir) / "ribdif_logs"
    Path(log_dir).mkdir(exist_ok = True, parents = True)
    

    # If rerun is false, download and handle genomes from NCBI
    if rerun == False:
        
        if args.genus:
            # Download genomes from NCBI
            status = ngd_download.genome_download(genus, outdir, args.threads, args.frag, args.sp_ignore, args.domain, logger)
            # Catching is any critical errors occured from downloading genomes
            if status == 1:
                sys.exit(status)
        
            # Un gziping fasta files
            with multiprocessing.Pool(args.threads) as pool: # Create a multiprocessing pool with #threads workers
                all_gz = [str(i) for i in list(Path(f"{outdir}/refseq/{args.domain}/").glob('**/*.gz'))]# Recursively search the directory for .gz files and convert path to string sotring in a list
                pool.map(utils.decompress, all_gz)
            
                
            # Remove unwanted characters from anywhere is file (should only be in fasta headers)
            logger.info("Modifying fasta headers.\n\n")
            with multiprocessing.Pool(args.threads) as pool:
                all_fna = [str(i) for i in list(Path(f"{outdir}/refseq/{args.domain}/").glob('*/*.fna'))]
                all_species = pool.map(utils.modify2, all_fna)
                genome_count = len(all_species)
                
        # Else if user defined  are given
        elif args.user:
            new_dir_path, genome_count = utils.own_genomes_copy(args.user, outdir, args.domain, logger) # copy to new location
            utils.own_genomes_gzip(new_dir_path) # decompress if needed
            utils.own_genomes_rename(new_dir_path, logger) # Rename fasta headers
            logger.info(f"{genome_count} user defined genomes were found\n\n")
    
    # Else get genome count of exising genomes      
    else:
        all_fna = [str(i) for i in list(Path(f"{outdir}/refseq/{args.domain}/").glob('*/*.fna'))]
        with multiprocessing.Pool(args.threads) as pool:
            all_species = pool.map(utils.sp_check, all_fna)
        genome_count = len(all_species)
        if args.genus:
            logger.info(f"{genome_count} previously downloaded genomes of {genus} were found\n\n")
        elif args.user:
            logger.info(f"{genome_count} previously user defined genomes were found\n\n")
    
    if args.genus:
        # getting unique species for later use in the overlap reports and removing "sp."
        unique_species = set(all_species)
        unique_species.discard("sp.")
    elif args.user:
        unique_species = set()
        all_species = ["sp."] * genome_count

            
    # If not using whole-genome mode assume the primers being used are 16S (which they are if default)
    if not args.whole:
        logger.info("#= Running barrnap on downloaded sequences =#\n\n")
        barrnap_run.barnap_call(outdir, threads = args.threads)
        
        # Processing barrnap output > fishing out 16S sequences
        with multiprocessing.Pool(args.threads) as pool:
            all_RNA = [str(i) for i in list(Path(f"{outdir}/refseq/{args.domain}/").glob('*/*.rRNA'))]
            #gene_num = [*range(len(all_RNA))] # adding gene num count here (in congruence with v1. Could also just use in silico pcr amp count)
            pool.map(barrnap_run.barrnap_process, all_RNA) # removed zip(gene_num) and was originally starmap
        
        # Concatinate all 16S to one file
        barrnap_run.barrnap_conc(genus, outdir)
        
        #If ANI is true calculate
        if args.ANI:
            # First need to split whole 16S sequences into seperate files
            with multiprocessing.Pool(args.threads) as pool:
                all_16S = [str(i) for i in list(Path(f"{outdir}/refseq/{args.domain}/").glob('*/*.16S'))]
                pool.map(barrnap_run.barrnap_split, all_16S)
               
            # Call pyani
            logger.info("Calculating intra-genomic mismatches and ANI for each genome.\n\n")
            pyani_run.pyani_call(outdir, args.threads, args.domain)
        else:
            logger.info("Skipping detailed intra-genomic analysis and ANI (if needed, use -a/--ANI).\n\n")

        # ALignment of full 16S genes recoverd from barrnap
        logger.info("Alligning full-length 16S genes within genomes with muscle.\n\n")
        msa_run.muscle_call_multi(outdir, args.threads, args.domain)
        
        summary_type = "16S"
        in_fna = f"{outdir}/full/{genus}.16S"
        summary_files.make_summary(in_fna, outdir, genus, args.whole, args.ANI, args.threads, summary_type, args.domain, args.user)
        
        # Running msa on concatinated 16S sequences
        if args.msa == True:
            logger.info(f"Alligning all {genus} 16S rRNA genes with muscle and building tree with fasttree.\n")
            infile , outAln, outTree = f"{outdir}/full/{genus}.16S", f"{outdir}/full/{genus}.16sAln", f"{outdir}/full/{genus}.16sTree" # Asigning in and out files
            msa_run.muscle_call_single(infile, outAln, outTree)
        else:
            logger.info("Skipping alignments and tree generation for 16S rRNA genes (if needed, use -m/--msa).\n\n")
            
        # PCR for default primers
        infile = f"{outdir}/full/{genus}.16S" # path to concatinated 16S barrnap output
        names =  pcr_run.pcr_call(infile, outdir, genus, primer_file, workingDir, logger)
        

        

    # PCR for custom primers   
    elif args.whole:

        names = pcr_run.pcr_parallel_call(outdir, genus, primer_file, workingDir, args.threads, logger, args.domain)
        
        # book keeping for parallel PCR
        for name in names:
            total_sum_dict = pcr_run.multi_cleaner(outdir, name)
            pcr_run.amplicon_cat(outdir, genus, name)
            pcr_run.sum_dict_write(outdir, genus, name, total_sum_dict)

        
    
    # Rename amplicon fasta headers to origin contig and removing any primers that did not amplify
    names = utils.amp_replace(outdir, genus, names, logger)
        
    # Catching if all amplification failed (empty lists evaluate to false)
    if not list(Path(f"{outdir}/amplicons/").rglob(f"{genus}-*.amplicons")):
        sys.exit("No amplification for any of the given primers was successfull. Try again with different primers")
    
    # Make summary file for whole genome mode (has to be after utils.amp_replace so cant have in main args.whole section)
    #if args.whole:
    for name in names:
        summary_type = f"{name}-amp"
        in_fna = f"{outdir}/amplicons/{name}/{genus}-{name}.amplicons"
        summary_files.make_summary(in_fna, outdir, genus, args.whole, args.ANI, args.threads, summary_type, args.domain, args.user)
        

    
    logger.info ("Making unique clusters with vsearch.\n\n")
    for name in names:
        vsearch_run.vsearch_call(outdir, genus, name, args.id, log_dir, args.threads, logger)
    
    
    
    # Generating the figures #
    Path.mkdir(Path(f"{outdir}/figures"), exist_ok = True)
    logger.info("Making reports and figures\n\n")
    if args.msa:
        logger.info("Aligning all amplicons for diversity calculation.\n")
        logger.info(f"Tree tip information can be found at: {outdir}/amplicons/{genus}-<name>-meta.tsv.\n\n")
    else:
        logger.info("Skipping total amplicon alignment and diversity calculation\n")
    
    with multiprocessing.Pool(args.threads) as pool: # Create a multiprocessing pool with #threads workersfor .gz files and convert path to string sotring in a list
        pool.starmap(utils.make_reports, zip(names, repeat(args.msa), repeat(outdir), repeat(genus), repeat(logger), repeat(args.user), repeat(unique_species), repeat(all_species), repeat(genome_count)))

    logger.info(f"You can find a saved version of the above at {outdir}/ribdif_log_file.log")
    
if __name__ == '__main__':
    main()