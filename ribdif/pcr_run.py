# -*- coding: utf-8 -*-
import multiprocessing
import subprocess
from pathlib import Path
import shlex


# Spawning the shell call
def call_proc_pcr(infile, outdir, genus, name, fwd, rvs):
    # Building the command
    command = f"in_silico_PCR.pl -s -s {fwd} -b {rvs} -r -m -i > {outdir}/amplicons/{genus}-{name}.summary 2> {outdir}/amplicons/{genus}-{name}.temp.amplicons"
    #command = f"in_silico_PCR.pl -s $genus/full/$genus.16S -a $fwd    -b $rev -r -m -i > $genus/amplicons/$genus-$name.summary 2> $genus/amplicons/$genus-$name.temp.amplicons"
    # Passing the command to shell piping the stdout and stderr
    p = subprocess.Popen(shlex.split(command), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out, err)

# Multithreading the in silico pcr calls
def pcr_parallel_call(outdir, genus, primer_file):
    print("\n\nRunning barrnap on downloaded sequences")
    Ncpu = multiprocessing.cpu_count()
    with open(primer_file, "r") as f_primer:
        for primer in f_primer:
            name, fwd, rvs = primer.split("\t")
            with multiprocessing.Pool(Ncpu) as pool: # spawn the pool
                all_fna = [str(i) for i in list(Path(f"{outdir}/genbank/bacteria/").rglob('*.fna'))] # generate list of files ending in .fna
                pool.apply(call_proc_pcr, args = (all_fna, outdir, genus, name, fwd, rvs))
    return

def pcr_call(infile, outdir, genus, primer_file):
    print("\n\nRunning barrnap on downloaded sequences")
    with open(primer_file, "r") as f_primer:
        for primer in f_primer:
            name, fwd, rvs = primer.split("\t")
            call_proc_pcr(infile, outdir, genus, name, fwd, rvs)
    return

