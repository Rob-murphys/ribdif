# -*- coding: utf-8 -*-
from glob import glob
import multiprocessing
import subprocess
import shlex
from pathlib import Path

# Spawning the shell call
def call_proc_pyani(indir):
    aniDir = Path(indir).parent / "ani"
    # Building the command
    command = f"average_nucleotide_identity.py -i {indir} -o {aniDir}"
    # Passing the command to shell piping the stdout and stderr
    subprocess.run(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return


# Multithreading the pyani calls
def pyani_call(outdir):
    Ncpu = multiprocessing.cpu_count()
    with multiprocessing.Pool(Ncpu) as pool: # spawn the pool
        all_16S_dirs = [i for i in glob(f"{outdir}/genbank/bacteria/*/indiv_16S_dir")] # can we use a faster method than glob?
        pool.map(call_proc_pyani, all_16S_dirs)
    return

