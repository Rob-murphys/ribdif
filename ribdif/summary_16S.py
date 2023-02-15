#!/usr/bin/env python3

from pathlib import Path
import multiprocessing
from glob import glob
from itertools import repeat

import numpy as np
import pandas as pd

from Bio import SeqIO
from collections import Counter

from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt

def writer(outdir, genus, stats):
    with open(f"{outdir}/{genus}_summary.tsv", "a") as f_app:
        f_app.write("\t".join(stats) + "\n")

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

def summary_16S_run(in_aln, outdir, genus):
    
    indir = Path(in_aln).parent
    # Paths for all needed files
    mismatch_path   = f"{indir}/ani/ANIm_similarity_errors.tab"
    alignment_path = in_aln
    tree_path = in_aln.replace(".16sAln", ".16sTree")
    pdf_out = in_aln.replace(".16sAln", "_16S_div.pdf")
    
    # Getting sequences name details from first record of the alignment fasta
    with open(alignment_path, "r") as f_in:
        count_16S = 0
        for value, line in enumerate(f_in):
            if value == 0:
                splitname = line.strip().split("_")
                count_16S += 1
            elif line.startswith(">"):
                count_16S += 1
        
        
    
    # Taking fasta header and getting important bits out
    GCF = "_".join(splitname[:2]).strip(">")
    #NZ = splitname[2] # Currently unused?
    genera = splitname[4]
    species = splitname[5]
    
    # If ANI was run then do following
    if Path(mismatch_path).is_file():
        
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
            
            total_div = str(shannon_calc(alignment_path))
            
            # roll_means_30 = np.convolve(divs, np.ones(30), "valid")/30   
            
# =============================================================================
#             if fast_mode:
#                 pass
#             else:
#                 tree = Phylo.read(tree_path, "newick")
#                 plot_tree(tree, pdf_out)
# =============================================================================
            
            stats = [GCF, genera, species, str(count_16S),  mean_mis, sd_mis, max_mis, min_mis, total_div]
            writer(outdir, genus, stats)

        else:
            mean_mis, sd_mis, max_mis, min_mis, total_div = (str(0), str(0), str(0), str(0), str(0))
            stats = [GCF, genera, species, str(count_16S), mean_mis, sd_mis, max_mis, min_mis, total_div]
            writer(outdir, genus, stats)
    else:
        mean_mis, sd_mis, max_mis, min_mis, total_div = ("-", "-", "-", "-", "-")
        stats = [GCF, genera, species, str(count_16S), mean_mis, sd_mis, max_mis, min_mis, total_div]
        writer(outdir, genus, stats)
        
            
def multiproc_sumamry(outdir, genus, threads):
    with multiprocessing.Pool(threads) as pool:
        all_aln = [str(i) for i in list(Path(f"{outdir}/refseq/bacteria/").glob('*/*.16sAln'))]
        pool.starmap(summary_16S_run, zip(all_aln, repeat(outdir), repeat(genus)))
    return

        



