#!/usr/bin/env python3

import subprocess
import shlex
from pathlib import Path
import os
import logging

def vsearch_call(outdir, genus, name, ident, log_dir, threads, logger):
    # Defining in and out files/dirs
    infile = f"{outdir}/amplicons/{name}/{genus}-{name}.amplicons"
    outfile = f"{outdir}/amplicons/{name}/{genus}-{name}.uc"
    clusdir = f"{outdir}/amplicons/{name}/{name}-clusters/{genus}-{name}-clus"
    
    Path.mkdir(Path(clusdir).parent, exist_ok = True, parents = True) # making outdir
    
    # Building the command
    command = f"vsearch --cluster_fast {infile} --id {ident} --strand both --uc {outfile} --clusters {clusdir} --quiet --threads {threads}"
    with open(f"{log_dir}/vsearch_{name}.out", "w") as out, open(f"{log_dir}/vsearch_{name}.err", "w") as err:
        # Passing the command to shell piping the stdout and stderr
        subprocess.run(shlex.split(command), stdout = out, stderr = err)
    
    if os.stat(f"{log_dir}/vsearch_{name}.err").st_size != 0:
        logger.warning(f"Vsearch did something non default. Hopefully this does not ruin down stream analysis but if you get an error check {log_dir}/vsearch_{name}.err")
    return
