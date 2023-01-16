# -*- coding: utf-8 -*-

import multiprocessing
from pathlib import Path
import subprocess
import shlex

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
    command1 = f"muscle -align {infile} -output {outAln}" # Do I need to strip the alignement file of white space and commas?
    command2 = f"fasttree -quiet -nopr -gtr -nt {outAln}"
    subprocess.run(shlex.split(command1), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    with open(outTree, "w") as f_std:
        subprocess.run(shlex.split(command2), stdout = f_std, stderr = subprocess.PIPE)
    return

