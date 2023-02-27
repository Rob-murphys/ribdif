#!/usr/bin/env python3

#from pathlib import Path
import multiprocessing
from glob import glob
from itertools import repeat

import numpy as np
import pandas as pd

from Bio import SeqIO
from collections import Counter

# =============================================================================
# from Bio import Phylo
# import matplotlib
# import matplotlib.pyplot as plt
# =============================================================================

def writer(outdir, genus, stats, summary_type):
    with open(f"{outdir}/{genus}_{summary_type}_summary.tsv", "a") as f_sum:
        f_sum.write("\t".join(stats) + "\n")

# =============================================================================
# def plot_tree(tree, pdf_out):
#     plt.figure(figsize=(25, 35), dpi=200)
#     Phylo.draw(tree, do_show = False)
#     plt.savefig(pdf_out)
#     return
# =============================================================================

def shannon_calc(alignment_path):
            # Read in alignment file and generate a dist of sequences
            seq_dict = {rec.id : rec.seq for rec in SeqIO.parse(alignment_path, "fasta")}
            n_seqs = len(seq_dict)
            seq_len = len(next(iter(seq_dict.values())))
            
            divs = [-1] * seq_len # empty list to record local shannon diversity on each nucleotide position
            # for each nucleotide position generate local shannon diversity
            for i in range(0,seq_len):
                pos = list(value[i] for value in seq_dict.values()) # get all values of a given position
                # Replicating the function of table() in R using a Counter
                c = Counter(pos)
                probs = np.array(list(c[x] for x in ["A", "G", "C", "T", "-"])) / n_seqs # getting count for each base in a numpy array
                divs[i] = -sum(probs * np.log(probs, where = probs != 0)) # log the probabilities where they are not 0 then sum them
            
            return sum(divs)

# =============================================================================
# def summary_16S_run(in_aln, outdir, genus):
#     
#     indir = Path(in_aln).parent
#     # Paths for all needed files
#     mismatch_path   = f"{indir}/ani/ANIm_similarity_errors.tab"
#     alignment_path = in_aln
#     #tree_path = in_aln.replace(".16sAln", ".16sTree")
#     #pdf_out = in_aln.replace(".16sAln", "_16S_div.pdf")
#     
#     # Getting sequences name details from first record of the alignment fasta
#     with open(alignment_path, "r") as f_in:
#         count_16S = 0
#         for value, line in enumerate(f_in):
#             if value == 0:
#                 splitname = line.strip().split("_")
#                 count_16S += 1
#             elif line.startswith(">"):
#                 count_16S += 1
#         
#         
#     
#     # Taking fasta header and getting important bits out
#     GCF = "_".join(splitname[:2]).strip(">")
#     #NZ = splitname[2] # Currently unused?
#     genera = splitname[4]
#     species = splitname[5]
#     
#     # Calculate total shanon diversity
#     total_div = str(shannon_calc(alignment_path))
#     
#     # If ANI was run then do following
#     if Path(mismatch_path).is_file():
#         
#         # Read in ANI simmilatiry errors (the number of unaligned or non-identical bases)
#         mismatch = pd.read_table(mismatch_path, sep="\t", index_col = 0)
#         
#         if len(mismatch) > 1:
#             # Get upper triange of this dataframe
#             indices  = np.triu_indices(mismatch.shape[0], k=1) # gets indices of upper triangle k=1 ofsets by 1 to avoid the diagonal center 
#             upper_mismatch = mismatch.values[indices] # gets the values of these indicies as a n array
#             
#             # Calculate stats on this array
#             mean_mis = str(round(np.mean(upper_mismatch), 2))
#             sd_mis = str(round(np.std(upper_mismatch, ddof = 1), 2))
#             max_mis = str(max(upper_mismatch))
#             min_mis = str(min(upper_mismatch))
#             
#             # roll_means_30 = np.convolve(divs, np.ones(30), "valid")/30   
#             
# # =============================================================================
# #             if fast_mode:
# #                 pass
# #             else:
# #                 tree = Phylo.read(tree_path, "newick")
# #                 plot_tree(tree, pdf_out)
# # =============================================================================
#             
#             stats = [GCF, genera, species, str(count_16S),  mean_mis, sd_mis, max_mis, min_mis, total_div]
#             writer(outdir, genus, stats)
# 
#         else:
#             mean_mis, sd_mis, max_mis, min_mis = (str(0), str(0), str(0), str(0))
#             stats = [GCF, genera, species, str(count_16S), mean_mis, sd_mis, max_mis, min_mis, total_div]
#             writer(outdir, genus, stats)
#     else:
#         mean_mis, sd_mis, max_mis, min_mis = ("-", "-", "-", "-")
#         stats = [GCF, genera, species, str(count_16S), mean_mis, sd_mis, max_mis, min_mis, total_div]
#         writer(outdir, genus, stats)
# =============================================================================
        
            
def summary_multiproc(outdir, genus, threads, file_list):
    with multiprocessing.Pool(threads) as pool:
        pool.starmap(make_summary, zip(file_list, repeat(outdir), repeat(genus)))
    return


def dict_parser(key, summary_dict, outdir, genus, whole_mode, ani_mode, summary_type, domain):
    # Looping through the generate dictionary
    value = summary_dict[key]
    if whole_mode == "Off": # is using barrnap (i.e. running ONLY on 16S genes)
        alignment_path = glob(f"{outdir}/refseq/bacteria/{key}/*.16sAln")[0] # generate path to alignment file
        summary_dict[key][7] = str(shannon_calc(alignment_path))# Calculate total shanon diversity
    if ani_mode == "On": # if using ani
        value = ani_stats(key, value, outdir, genus, domain) # get ani stats
    # Write it all to file
    value.insert(0, key) # append GCF to value list
    writer(outdir, genus, value, summary_type) # Write
    return
        
def make_summary(in_fna, outdir, genus, whole_mode, ani_mode, threads, summary_type, domain):
    
    # Initiate the summary file with headers
    with open(f"{outdir}/{genus}_{summary_type}_summary.tsv", "w") as f_out:
        f_out.write("GCF\tGenus\tSpecies\t#16S\tMean\tSD\tMin\tMax\tTotalDiv\n")
    # Paths for all needed files
    #tree_path = in_aln.replace(".16sAln", ".16sTree")
    #pdf_out = in_aln.replace(".16sAln", "_16S_div.pdf")
    
    
    summary_dict = {}
    # Loop through each fasta file looking for headers
    with open(in_fna, "r") as f_in:
        for line in f_in:
            if ">" in line:
                splitname = line.strip().split("_")
                # Taking fasta header and getting important bits out
                GCF = "_".join(splitname[:2]).strip(">")
                genera = splitname[4]
                species = splitname[5]
                # If the GCF already exists in the dict then just plus one to count
                if GCF in summary_dict:
                    summary_dict[GCF][2] = str(int(summary_dict[GCF][2]) + 1)
                # If it does not exist create it with space for all needed info
                else:
                    count = "1"
                    summary_dict[GCF] = [genera, species, count, "-", "-", "-", "-", "-"]

    # Multipocess the writing shannon and ani stuff out (probably dont need to do this)
    with multiprocessing.Pool(threads) as pool:
        pool.starmap(dict_parser, zip(summary_dict, repeat(summary_dict), repeat(outdir), repeat(genus), repeat(whole_mode), repeat(ani_mode), repeat(summary_type), repeat(domain)))
    return


def ani_stats(key, value, outdir, genus, domain):
    mismatch_path   = f"{outdir}/refseq/{domain}/{key}/ani/ANIm_similarity_errors.tab"
    
        
    # Read in ANI simmilatiry errors (the number of unaligned or non-identical bases)
    mismatch = pd.read_table(mismatch_path, sep="\t", index_col = 0)
    
    if len(mismatch) > 1:
        # Get upper triange of this dataframe
        indices  = np.triu_indices(mismatch.shape[0], k=1) # gets indices of upper triangle k=1 ofsets by 1 to avoid the diagonal center 
        upper_mismatch = mismatch.values[indices] # gets the values of these indicies as a n array
        
        # Calculate stats on this array
        mean_mis = str(round(np.mean(upper_mismatch), 2))
        sd_mis = str(round(np.std(upper_mismatch, ddof = 1), 2))
        max_mis = str(max(upper_mismatch))
        min_mis = str(min(upper_mismatch))
        
        # roll_means_30 = np.convolve(divs, np.ones(30), "valid")/30   
        
# =============================================================================
#             if fast_mode:
#                 pass
#             else:
#                 tree = Phylo.read(tree_path, "newick")
#                 plot_tree(tree, pdf_out)
# =============================================================================
        
        value[3:7] = mean_mis, sd_mis, max_mis, min_mis
        return value

    else:
        value[3:7] = (str(0), str(0), str(0), str(0))
        return value


