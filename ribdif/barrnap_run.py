#!/usr/bin/env python3
"""
To add:
    When a barrnap result file ends up empty(or just no 16S maybe?) take note and warn then user to check the barrnap logs and write to the logs which file had no 16S in it
"""
import multiprocessing
import subprocess
from pathlib import Path
import shlex


# Spawning the shell call
def call_proc_barrnap(infile):
    # Building the command
    command = f"barrnap --kingdom bac --quiet --threads 1 --reject 0.90 -outseq {infile}.rRNA {infile}"
    # Passing the command to shell piping the stdout and stderr
    subprocess.run(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return

# Multithreading the barrnap calls
def barnap_call(outdir, threads):
    with multiprocessing.Pool(threads) as pool: # spawn the pool
        all_fna = [str(i) for i in list(Path(f"{outdir}/refseq/bacteria/").glob('*/*.fna'))] # generate list of files ending in .rna
        pool.map(call_proc_barrnap, all_fna)
    return
        

# Fishes out 16S sequences and saves them to file
def barrnap_process(in_RNA):
    count = 0 # genome wide count: we are not going to use this
    GCF = str(Path(in_RNA).parent).split(r"/")[-1]
    with open(in_RNA, "r") as f_in, open(f"{in_RNA}.16S", "w") as f_out: 
        for line in f_in:
            if line.startswith(">16S_rRNA"):
                f_out.write(">" + GCF + "_" + line.strip().strip(">16S_rRNA::").split(":")[0] + f"_{count}\n") # write the fasta header to file removing >16S_rRNA:: and adding a counter
                f_out.write(next(f_in)) # write next line also
                count += 1 # incriment count by one
    return

# Splits barrnap output sequences from barrnap_process output
def barrnap_split(in_16S):
    dir_16S = (Path(in_16S).parent / "indiv_16S_dir") # directory path for individual 16S files
    dir_16S.mkdir(exist_ok=True, parents=True) # create directory
    with open(in_16S, "r") as f_in: 
        for line in f_in:
            if line.startswith(">"):
                count = line.strip()[-1]
                with open(f"{dir_16S}/{count}.fna", "w") as f_out:
                    f_out.write(line) # write the fasta header to file removing >16S_rRNA:: and adding a counter
                    f_out.write(next(f_in)) # write next line also
    return

# Concatinates fished out 16S barrnap output into a single file
def barrnap_conc(genus, outdir):
    all_16S = [str(i) for i in list(Path(f"{outdir}/refseq/bacteria/").glob('*/*.16S'))]
    full_path = Path(f"{outdir}/full")
    full_path.mkdir(parents = True, exist_ok = True)
    with open(f"{full_path}/{genus}.16S", "w") as f_out:
        for file in all_16S:
            with open(file, "r") as f_in:
                f_out.write(f_in.read())
    return
        
    