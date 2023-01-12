# -*- coding: utf-8 -*-
#import multiprocessing
import subprocess
from pathlib import Path
#from itertools import repeat
import shlex
"""
Implement a producer and consumer setup for writing the pcr output whe multiprocessing: https://stackoverflow.com/questions/11196367/processing-single-file-from-multiple-processes

I removed the counter from all functions for labeling output files as it seems the cocatinated fasta ran pretty fast.
"""

# Spawning the shell call
def call_proc_pcr(infile, outdir, genus, name, fwd, rvs, workingDir): # removed , counter
    # Building the command
    command = f"perl {workingDir}/in_silico_PCR.pl -s {infile} -a {fwd} -b {rvs} -r -m -i > {outdir}/amplicons/{genus}-{name}.summary 2> {outdir}/amplicons/{genus}-{name}.temp.amplicons"
    print(command)
    # Passing the command to shell piping the stdout and stderr
    with open(f"{outdir}/amplicons/{genus}-{name}.summary", "w") as f_std, open(f"{outdir}/amplicons/{genus}-{name}.temp.amplicons", "w") as f_err:
        subprocess.run(shlex.split(command), stdout = f_std, stderr = f_err)
    return 

# Multithreading the in silico pcr calls
# =============================================================================
# def pcr_parallel_call(outdir, genus, primer_file, workingDir):
#     amplicon_dir = Path(f"{outdir}/amplicons")
#     amplicon_dir.mkdir(parents = True, exist_ok = True)
#     print("\n\nGenerating amplicon sequences")
#     Ncpu = multiprocessing.cpu_count()
#     with open(primer_file, "r") as f_primer:
#         for primer in f_primer:
#             name, fwd, rvs = primer.split("\t")
#             with multiprocessing.Pool(Ncpu) as pool: # spawn the pool
#                 all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.fna'))] # generate list of files ending in .fna
#                 counter = range(len(all_fna))
#                 pool.starmap(call_proc_pcr, zip(all_fna, repeat(outdir), repeat(genus), repeat(name), repeat(fwd), repeat(rvs), repeat(workingDir), counter))
#     return
# =============================================================================

def pcr_call(infile, outdir, genus, primer_file, workingDir):
    amplicon_dir = Path(f"{outdir}/amplicons")
    amplicon_dir.mkdir(parents = True, exist_ok = True)
    print("\n\nGenerating amplicon sequences")
    with open(primer_file, "r") as f_primer:
        for primer in f_primer:
            name, fwd, rvs = primer.strip().split("\t")
            call_proc_pcr(infile, outdir, genus, name, fwd, rvs, workingDir)
    return name

# =============================================================================
# def pcr_cleaner(outdir, primer_file, genus):
#     print("cleaning files")
#     amplicon_dir = Path(f"{outdir}/amplicons")
#     with open(primer_file, "r") as f_primer:
#         for primer in f_primer:
#             name = primer.split("\t")[0]
#             print(f"looking for {name}")
#             
#             # Summary file
#             summary_file = f"{amplicon_dir}/{genus}-{name}.summary"
#             print(f"master summary: {summary_file}")
#             all_sum = [str(i) for i in list(Path(f"{amplicon_dir}/").rglob(f"*{name}.summary*"))] # getting all summary files
#             for file in all_sum: # looping over them
#                 with open(file, "r") as f_in: # open each one
#                     if Path(summary_file).is_file(): # if the master summary file already exists
#                         with open(f"{summary_file}", "a") as f_sum: # open it in append mode
#                             lines = f_in.read().splitlines()
#                             f_sum.write("\n".join(lines[1:])) # and write everything from the second line
#                     else:
#                         with open(f"{summary_file}", "w") as f_sum:
#                             f_sum.write(f_in.read())
#             # Amplicon file
# 
# =============================================================================