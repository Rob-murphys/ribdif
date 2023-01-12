# -*- coding: utf-8 -*-

import multiprocessing
from pathlib import Path
import subprocess
import shlex

# Spawning the shell call
def call_proc_muscle(infile):
    outfile = str(Path(infile).stem)
    # Building the command
    command1 = f"'muscle -in {infile} -out {outfile}.16sAln -quiet" # Do I need to strip the alignement file of white space and commas?
    command2 = f"fasttree -quiet -nopr -gtr -nt {outfile}.16sAln > {outfile}.16sTree"
    # Passing the command to shell piping the stdout and stderr
    subprocess.run(shlex.split(command1), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with open(f"{outfile}.16sTree", "w") as f_std:
        subprocess.run(shlex.split(command2), stdout=f_std, stderr=subprocess.PIPE)
    
    return


# Multithreading the pyani calls
def muscle_call(outdir):
    Ncpu = multiprocessing.cpu_count()
    with multiprocessing.Pool(Ncpu) as pool: # spawn the pool
        all_16S = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.16S'))]
        pool.map(call_proc_muscle, all_16S)
    return
