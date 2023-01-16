# -*- coding: utf-8 -*-

import subprocess
import shlex
from pathlib import Path


def vsearch_call(outdir, genus, name):
    # Defining in and out files/dirs
    infile = f"{outdir}/amplicons/{genus}-{name}.amplicons"
    outfile = f"{outdir}/amplicons/{genus}-{name}.uc"
    outdir = f"{outdir}/amplicons/{name}-clusters/{genus}-{name}-clus"
    
    Path.mkdir(Path(outdir), exist_ok = True, parents = True) # making outdir
    
    # Building the command
    command = "vsearch -cluster_fast {infile} --id $id  -strand both --uc {outfile} --clusters {outdir} --quiet"
    
    # Passing the command to shell piping the stdout and stderr
    subprocess.run(shlex.split(command), stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    return
