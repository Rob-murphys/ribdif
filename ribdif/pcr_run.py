# -*- coding: utf-8 -*-
import multiprocessing
import subprocess
from pathlib import Path
import shlex


# Spawning the shell call
def call_proc_pcr(infile, outdir, genus, name, fwd, rvs, workingDir):
    # Building the command
    command = f"perl {workingDir}/in_silico_PCR.pl -s {infile} -a {fwd} -b {rvs} -r -m -i > {outdir}/amplicons/{genus}-{name}.summary 2> {outdir}/amplicons/{genus}-{name}.temp.amplicons"
    print(command)
    # Passing the command to shell piping the stdout and stderr
    with open(f"{outdir}/amplicons/{genus}-{name}.summary", "w") as f_std, open(f"{outdir}/amplicons/{genus}-{name}.temp.amplicons", "w") as f_err:
        subprocess.run(shlex.split(command), stdout = f_std, stderr = f_err)
    return

# Multithreading the in silico pcr calls
def pcr_parallel_call(outdir, genus, primer_file, workingDir):
    amplicon_dir = Path(f"{outdir}/amplicons")
    amplicon_dir.mkdir(parents = True, exist_ok = True)
    print("\n\nGenerating amplicon sequences")
    Ncpu = multiprocessing.cpu_count()
    with open(primer_file, "r") as f_primer:
        for primer in f_primer:
            name, fwd, rvs = primer.split("\t")
            with multiprocessing.Pool(Ncpu) as pool: # spawn the pool
                all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.fna'))] # generate list of files ending in .fna
                pool.apply(call_proc_pcr, args = (all_fna, outdir, genus, name, fwd, rvs, workingDir))
    return

def pcr_call(infile, outdir, genus, primer_file, workingDir):
    amplicon_dir = Path(f"{outdir}/amplicons")
    amplicon_dir.mkdir(parents = True, exist_ok = True)
    print("\n\nGenerating amplicon sequences")
    with open(primer_file, "r") as f_primer:
        for primer in f_primer:
            name, fwd, rvs = primer.split("\t")
            call_proc_pcr(infile, outdir, genus, name, fwd, rvs, workingDir)
    return

