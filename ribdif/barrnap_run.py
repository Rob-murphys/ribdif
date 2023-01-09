# -*- coding: utf-8 -*-

import multiprocessing
import subprocess
from pathlib import Path
import shlex


# Spawning the shell call
def call_proc(infile):
    # Building the command
    command = f"barrnap --kingdom bac --quiet --threads 1 --reject 0.90 -outseq {infile}.rRNA {infile}"
    # Passing the command to shell piping the stdout and stderr
    p = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out, err)

# Multithreading the barrnap calls
def barnap_call(outdir):
    print("\n\nRunning barrnap on downloaded sequences")
    Ncpu = multiprocessing.cpu_count()
    with multiprocessing.Pool(Ncpu) as pool: # spawn the pool
        all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.fna'))] # generate list of files ending in .rna
        pool.map(call_proc, all_fna) # 
    #pool.join() # The pool will close and wait for each running task to complete

def barrnap_split(RNA):
    count = 0
    dir_16S = (Path(RNA).parent / "indiv_16S_dir") # directory path for individual 16S files
    dir_16S.mkdir(exist_ok=True, parents=True) # create directory
    with open(RNA, "r") as f_in: 
        for line in f_in:
            if line.startswith(">16S_rRNA"):
                with open(f"{dir_16S}/{count}.fna", "w") as f_out:
                    f_out.write(">" + line.strip().strip(">16S_rRNA::") + f"_{count}\n") # write the fasta header to file removing >16S_rRNA:: and adding a counter
                    f_out.write(next(f_in)) # write next line also
                count += 1 # incriment count by one

                
    

