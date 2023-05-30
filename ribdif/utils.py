#!/usr/bin/env python3
import pandas as pd
import re
import gzip
import shutil
from pathlib import Path
import fileinput
import os
import logging
import chardet
from ribdif import overlaps, figures, msa_run, summary_files


# =============================================================================
# def modify(file_path):
#     # Open and read the contents of the file
#     with open(file_path, "r") as f_in:
#         contents = f_in.read()
#         
#     # Preform substitutions for unwanted characters
#     contents = re.sub(r"[:,/()=#\x27]", "", contents)
#     contents = re.sub(r"[: ]", "_", contents)
# 
#     # Write the modified contents back to file
#     with open(file_path, "w") as f_out:
#         f_out.write(contents)
#     return
# =============================================================================

# Remove unwanted characters from anywhere is file (should only be in fasta headers)
def modify2(file_path):
    GCF = str(Path(file_path).parent).split(r"/")[-1]
    # Open and read the contents of the file
    with fileinput.input(file_path, inplace = True) as f_in:
        for line in f_in:
            if line.startswith(">"):
                line = re.sub(r"[:,/()\[\]=#\x27]", "", line)
                line = re.sub(r"[: ]", "_", line)
                line = f"{line[0]}{GCF}_{line[1:]}" 
                print(line, end = '')
                current_sp = line.split("_")[5]
            else:
                print(line, end = '')
    return current_sp
# =============================================================================
# def modify3(file_path):
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

# open and check species of the genome
def sp_check(file):
    with open(file, "r") as f_in:
        line = f_in.readline()
        current_sp = line.split("_")[5]
    return current_sp

# Un gzip the downloaded NCBI genomes
def decompress(file_path):
    # Open the .gz file and decompress it
    with gzip.open(file_path, "rb") as f_in, open(file_path[:-3], "wb") as f_out:
        shutil.copyfileobj(f_in, f_out) # copy
    return

# Rename amplicon fasta headers to origin contig
def amp_replace(outdir, genus, names, logger):
    for name in names[:]:
        # Read in the summary dataframe
        df_sum = pd.read_csv(f"{outdir}/amplicons/{name}/{genus}-{name}.summary", sep = "\t", header = None, names = ["AmpId", "SequenceId", "PositionInSequence", "Length", "Misc"])
        dict_sum = dict(zip(df_sum.AmpId, df_sum.SequenceId)) # Make a dictionary of it
        # Open the temp amplicon file and the final amplicon file
        with open (f"{outdir}/amplicons/{name}/{genus}-{name}.temp.amplicons", "r") as f_in, open(f"{outdir}/amplicons/{name}/{genus}-{name}.amplicons", "w") as f_out:
           for line in f_in: # loop over lines in the file
               if ">amp" in line: # If it is a fasta header
                   # Replace with origin genome fasta header and then the amp count
                   amp = line.strip().strip(">")
                   line = line.replace(amp, dict_sum[amp] + f"_{amp.strip('amp_')}")
               f_out.write(line)
        Path.unlink(Path(f"{outdir}/amplicons/{name}/{genus}-{name}.temp.amplicons")) # remove temp file
        # Check if the primer resulted in any amplification and if not remove the file
        if os.stat(f"{outdir}/amplicons/{name}/{genus}-{name}.amplicons").st_size == 0:
            Path.unlink(Path(f"{outdir}/amplicons/{name}/{genus}-{name}.amplicons"))
            logger.info(f"{name} primer resulted in no amplification and will be excluded from further analysis. Are you sure the primer is correct?\n")
            names.remove(name)
    return names

def pairwise_to_csv(pairwise_match, gcf_species, outdir, genus, name):
    pairwise_save_df = pd.DataFrame(pairwise_match, index = gcf_species.values())
    pairwise_save_df.to_csv(f"{outdir}/amplicons/{name}/{genus}-{name}_confusion.csv", sep = ",", index = True)
    return

def detect_encode(file):
    with open(file, "rb") as f_in:
        return chardet.detect(f_in.readline())["encoding"]
    
def primer_check(primer_file, logger):
    encoding = detect_encode(primer_file)
    with open(primer_file, "r", encoding = encoding) as f_in:
        for line in f_in:
            if len(line.strip().split("\t")) != 4:
                logger.info(f"""Your primer file in in the incorrect format at the following line:
                            \t{line}
                            Please follow the 'name | forward | reverse | expected amplicon length' tab seperated format
                            Check the default.primers file or the README for more guidance""")
                return True
    return False

# Copy the user defined genomes to the output directory following ncbi-genome.download structure
def own_genomes_copy(dir_path, outdir, domain, logger):
    target_dir = f"{outdir}/refseq/{domain}" # Target directory
    Path.mkdir(Path(target_dir), exist_ok = True, parents = True) # Creating directory
    logger.info(f"Copying your genomes to {target_dir}\n\n")
    file_count = 0 # Initiate genome count
    # Looping over whole directory
    for file in Path(dir_path).iterdir():
        if file.is_file(): # if item is a file
            final_dir = f"{target_dir}/{file.stem.replace('_', '-')}"
            Path.mkdir(Path(final_dir))
            shutil.copy(file, f"{final_dir}/{file.stem.replace('_', '-')}.fna") # copy it replacing the file extension with '.fna'
            file_count += 1 # incriment file count
    return target_dir, file_count

# Decompress user defined files if they need it
def own_genomes_gzip(new_dir_path):
    for file in Path(new_dir_path).glob("*/*.fna"):
        if file.suffix == ".gz":
            decompress(file)
    return
    
# Rename the fasta headers of the file inplace
def own_genomes_rename(new_dir_path, logger):
    NZ_count = 1 # abritrary GCF
    logger.info("Changing the copied genomes headers to conform with NCBI format\n\n")
    # Loop throuh all user fna files in new directory
    for file in Path(new_dir_path).glob("*/*.fna"):
        with fileinput.input(file, inplace = True) as f_in: # open the file for inplace editing
            for line in f_in: # Loop through lines
                if line.startswith(">"): # if it is a fasta header
                    # Generate new fasta header    
                    genus = file.stem.replace("_", "-")
                    line = f">GCF_{genus}_NZ_CP{NZ_count}_{genus}_sp._placeholder\n" # generate random GCF
                    print(line, end = '')
                    NZ_count += 1
                else:
                    print(line, end = '') # else print the line (which would be actual sequence data)
    return


def make_reports(name, msa, outdir, genus, logger, user, unique_species, all_species, genome_count):

    if msa:
        # msa on all amplicons
        infile , outAln, outTree = f"{outdir}/amplicons/{name}/{genus}-{name}.amplicons", f"{outdir}/amplicons/{name}/{genus}-{name}.aln", f"{outdir}/amplicons/{name}/{genus}-{name}.tree" # Asigning in and out files
        msa_run.muscle_call_single(infile, outAln, outTree)
        msa_run.format_trees(outdir, genus, name)
        # Calculate shannon diversity across the primers
        shannon_div = summary_files.shannon_calc(outAln)
    else:
        shannon_div = "Was skipped"
    
    
    # Cleaning vsearch clustering data
    all_gcfs, uc_dict_clean, gcf_species, cluster_count = overlaps.uc_cleaner(outdir, genus, name)
    
    # Generate a dictionary (that will become a matrix) of GCF cluster membership  
    cluster_dict = overlaps.cluster_matrix(all_gcfs, uc_dict_clean, cluster_count)
    
    # Find all species overlap in the cluster_dict
    combinations = overlaps.species_overlap(cluster_dict, cluster_count, gcf_species)
    
    # If only one cluster is made then skip all the figure making. Could I do this with cluster count instead and avoid the two commands above?
    if len(cluster_dict) != 1:
        # Find all GCF overlaps in the cluster dictionary
        pairwise_match = overlaps.gcf_overlaps(all_gcfs, uc_dict_clean, gcf_species)
        pairwise_to_csv(pairwise_match, gcf_species, outdir, genus, name)
        
        # Generate metadata for heatmaps
        row_palette, species_series, species_palette = figures.heatmap_meta(gcf_species)
        
        # Plot the cluster matrix
        plot_clus, cluster_df = figures.cluster_heatmap(cluster_dict, row_palette, species_series)
        
        # Plot the GCF overlap matrix
        plot_dendo, pairwise_df = figures.pairwise_heatmap(pairwise_match, row_palette, species_series)
        
        plot_clus = figures.figure_fix(plot_clus)
        plot_dendo = figures.figure_fix(plot_dendo)
        
        # Save the heatmaps
        figures.pdf_save(plot_clus, plot_dendo, outdir, genus, name)
    
        # Generate graps from the pairwise dataframe
        adjacency_df = figures.create_adjacency(pairwise_df, cluster_df)
        graph_subs, n_subplots = figures.create_graph(adjacency_df)
        
        if n_subplots != 0:
            # Draw the generated graps into on plot
            figures.draw_graphs(graph_subs, n_subplots, species_palette, row_palette, outdir, genus, name)
        else:
            logger.info(f"Skipping graph making for {name} as no edges were found (even within a single species)\n")
    else:
        logger.info(f"Only one genome amplified for the {name} primer ({list(cluster_dict)[0]}) so we will skip making figures as they would be useless\n")
        cluster_df = overlaps.single_amp_df(cluster_dict)
    overlaps.overlap_report(combinations, gcf_species, cluster_df, genus, name, outdir, logger, shannon_div, unique_species, all_species, genome_count, user)