#!/usr/bin/env python3

import multiprocessing
from pathlib import Path
import subprocess
import shlex
from Bio import Phylo
import pandas as pd
import re


# Spawning the shell call
def call_proc_muscle(infile):
    outfile = str(Path(infile).parent / Path(infile).stem)
    # Building the command
    #command1 = f"muscle -super5 {infile} -output {outfile}.16sAln -threads {threads} -nt" # Do I need to strip the alignement file of white space and commas?
    command1 = f"mafft --quiet {infile}"
    #command2 = f"fasttree -quiet -nopr -gtr -nt {outfile}.16sTree"
    # Passing the command to shell piping the stdout and stderr
    with open(f"{outfile}.16sAln", "w") as aln_out:
        subprocess.run(shlex.split(command1), stdout = aln_out, stderr = subprocess.PIPE)
# =============================================================================
#     with open(f"{outfile}.16sTree", "w") as f_std:
#         subprocess.run(shlex.split(command2), stdout = f_std, stderr = subprocess.PIPE)
# =============================================================================
    return


# Multithreading the msa calls
def muscle_call_multi(outdir, threads, domain):
    with multiprocessing.Pool(threads) as pool: # spawn the pool
        all_16S = [str(i) for i in list(Path(f"{outdir}/refseq/{domain}/").glob('*/*.16S'))]
        pool.map(call_proc_muscle, all_16S)
    return



# Calling Muscle for MSA of all 16S sequences
def muscle_call_single(infile, outAln, outTree):
    # Building the command
    #command1 = f"muscle -super5 {infile} -output {outAln} -threads {threads} -nt" # Do I need to strip the alignement file of white space and commas?
    command1 = f"mafft --quiet {infile}"
    command2 = f"fasttree -quiet -nopr -gtr -nt {outAln}"
    with open(f"{outAln}", "w") as aln_out:
        subprocess.run(shlex.split(command1), stdout = aln_out, stderr = subprocess.PIPE)
    with open(outTree, "w") as f_std:
        subprocess.run(shlex.split(command2), stdout = f_std, stderr = subprocess.PIPE)
    return

def format_trees(outdir, genus, name):
    uc_path = f"{outdir}/amplicons/{name}/{genus}-{name}.uc"
    tree_path = f"{outdir}/amplicons/{name}/{genus}-{name}.tree"
    out_path = f"{outdir}/amplicons/{name}/{genus}-{name}-meta.tsv"
    
    # read in cluster file
    uc_df = pd.read_csv(uc_path, sep = "\t", header = None)
    uc_df_clean = uc_df.loc[uc_df[0] != "C"]
    
    # read in tree and get lead labels
    tree = Phylo.read(tree_path, "newick")
    tree_terminals = tree.get_terminals() # getting all leaf info
    tips = [t.name for t in tree_terminals] # getting just leaf names
    tips_clean = [re.sub(r"_chromosome_.*|_complete_.*|_genome_.*", "", n) for n in tips] # cleaning leaf names

    
    # Getting individual name parts with list comprehension
    GCF = ["_".join(i.split("_")[:2]) for i in tips_clean]
    genus = [i.split("_")[4] for i in tips_clean]
    species = [i.split("_")[5] for i in tips_clean]
    strain = ["_".join(i.split("_")[6:])for i in tips_clean]
    n_gene = [i.split("_")[-1] for i in tips]
    
    # Getting cluster number for each item of tips in the order of tips # Why not just reorder tips to that of the dataframe?
    uc_df = uc_df.sort_values(by = 8, axis = 0)
    match_index =[uc_df_clean.index[uc_df_clean[8] == i] for i in tips]
    clusters = [uc_df.loc[m][1].values[0] for m in match_index]
    
    meta_df = pd.DataFrame({"name": tips, "GCF":GCF, "Genus":genus, "Species":species, "Strain":strain, "nGene":n_gene, "Cluster":clusters})
    meta_df.to_csv(out_path, sep = "\t", index = False)
    
  
