# -*- coding: utf-8 -*-

import multiprocessing
import subprocess


# Spawning the shell call
def call_proc(infile):
    # Building the command
    command = f'barrnap --kingdom "bac" --quiet --threads 1 --reject 0.90 -o "{infile}.rRNA" {infile}'
    # Passing the command to shell piping the stdout and stderr
    p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out, err)

# Multithreading the barrnap calls
def barnap_call(outdir):
    std = []
    Ncpu = multiprocessing.cpu_count()
    with multiprocessing.Pool(Ncpu) as pool: # spawn the pool
        all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.fna'))]
        std.append = pool.apply_async(call_proc, all_fna)
    pool.join() # The pool will close and wait for each running task to complete
    
    with open("barrnap.log", "w") as f_out:
        f_out.write(std)
    

