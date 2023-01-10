# -*- coding: utf-8 -*-

import multiprocessing
import subprocess
from pathlib import Path
import shlex


# Spawning the shell call
def call_proc_barrnap(infile):
    # Building the command
    command = f"barrnap --kingdom bac --quiet --threads 1 --reject 0.90 -outseq {infile}.rRNA {infile}"
    # Passing the command to shell piping the stdout and stderr
    p = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    print(out)
    print(err)
    return (out, err)

# Multithreading the barrnap calls
def barnap_call(outdir):
    print("\n\nRunning barrnap on downloaded sequences")
    Ncpu = multiprocessing.cpu_count()
    with multiprocessing.Pool(Ncpu) as pool: # spawn the pool
        all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.fna'))] # generate list of files ending in .rna
        pool.map(call_proc_barrnap, all_fna)
    return
        

# Fishes out 16S sequences and saves them to file
def barrnap_process(in_RNA):
    count = 0
    with open(in_RNA, "r") as f_in, open(f"{in_RNA}.16S", "w") as f_out: 
        for line in f_in:
            if line.startswith(">16S_rRNA"):
                f_out.write(">" + line.strip().strip(">16S_rRNA::") + f"_{count}\n") # write the fasta header to file removing >16S_rRNA:: and adding a counter
                f_out.write(next(f_in)) # write next line also
        count += 1 # incriment count by one
    return

# Splits barrnap output sequences from barrnap_process output
def barrnap_split(in_16S):
    dir_16S = (Path(in_16S).parent / "indiv_16S_dir") # directory path for individual 16S files
    dir_16S.mkdir(exist_ok=True, parents=True) # create directory
    with open(in_16S, "r") as f_in: 
        for line in f_in:
            if line.startswith(">16S_rRNA"):
                count = line[-1]
                with open(f"{dir_16S}/{count}.fna", "w") as f_out:
                    f_out.write(line) # write the fasta header to file removing >16S_rRNA:: and adding a counter
                    f_out.write(next(f_in)) # write next line also
    return

# Concatinates fished out 16S barrnap output into a single file
def barrnap_conc(genus, outdir):
    all_16S = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.16S'))]
    full_path = Path(f"{outdir}/full")
    full_path.mkdir(parents = True, exist_ok = True)
    with open(f"{full_path}/{genus}.16S", "w") as f_out:
        for file in all_16S:
            with open(file, "r") as f_in:
                f_out.write(f_in.read())
    return
        
    