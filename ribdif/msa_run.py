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
    command1 = f"muscle -align {infile} -output {outfile}.16sAln" # Do I need to strip the alignement file of white space and commas?
    command2 = f"fasttree -quiet -nopr -gtr -nt {outfile}.16sAln"
    # Passing the command to shell piping the stdout and stderr
    subprocess.run(shlex.split(command1), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    with open(f"{outfile}.16sTree", "w") as f_std:
        subprocess.run(shlex.split(command2), stdout = f_std, stderr = subprocess.PIPE)
    
    return


# Multithreading the pyani calls
def muscle_call_multi(outdir, threads):
    with multiprocessing.Pool(threads) as pool: # spawn the pool
        all_16S = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").glob('*/*.16S'))]
        pool.map(call_proc_muscle, all_16S)
    return



# Calling Muscle for MSA of all 16S sequences
def muscle_call_single(infile, outAln, outTree):
    # Building the command
    command1 = f"muscle -align {infile} -output {outAln}" # Do I need to strip the alignement file of white space and commas? # Do we want to use -distance1 kbit20_3
    command2 = f"fasttree -quiet -nopr -gtr -nt {outAln}"
    subprocess.run(shlex.split(command1), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    with open(outTree, "w") as f_std:
        subprocess.run(shlex.split(command2), stdout = f_std, stderr = subprocess.PIPE)
    return

def format_trees(outdir, genus, name):
    uc_path = f"{outdir}/amplicons/{genus}-{name}.uc"
    tree_path = f"{outdir}/amplicons/{genus}-{name}.tree"
    out_path = f"{outdir}/amplicons/{genus}-{name}-meta.tsv"
    
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
    NZ = [i.split("_")[2] for i in tips_clean]
    species = [i.split("_")[4] for i in tips_clean]
    strain = ["_".join(i.split("_")[5:])for i in tips_clean]
    n_gene = [i.split("_")[-1] for i in tips]
    
    # Getting cluster number for each item of tips in the order of tips # Why not just reorder tips to that of the dataframe?
    uc_df = uc_df.sort_values(by = tips)
    match_index =[uc_df_clean.index[uc_df_clean[8] == i] for i in tips]
    clusters = [uc_df.loc[m][1].values[0] for m in match_index]
    
    meta_df = pd.DataFrame({"name": tips, "GCF":GCF, "NZ":NZ, "Species":species, "Strain":strain, "nGene":n_gene, "Cluster":clusters})
    meta_df.to_csv(out_path, sep = "\t", index = False)
    
  
