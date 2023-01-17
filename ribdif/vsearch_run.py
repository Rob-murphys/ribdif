# -*- coding: utf-8 -*-

import subprocess
import shlex
from pathlib import Path
import os



def vsearch_call(outdir, genus, name, ident, log_dir, threads):
    # Defining in and out files/dirs
    infile = f"{outdir}/amplicons/{genus}-{name}.amplicons"
    outfile = f"{outdir}/amplicons/{genus}-{name}.uc"
    clusdir = f"{outdir}/amplicons/{name}-clusters/{genus}-{name}-clus"
    
    Path.mkdir(Path(clusdir).parent, exist_ok = True, parents = True) # making outdir
    
    # Building the command
    command = f"vsearch --cluster_fast {infile} --id {ident} --strand both --uc {outfile} --clusters {clusdir} --quiet --threads {threads}"
    with open(f"{log_dir}/vsearch_{name}.out", "w") as out, open(f"{log_dir}/vsearch_{name}.err", "w") as err:
        # Passing the command to shell piping the stdout and stderr
        subprocess.run(shlex.split(command), stdout = out, stderr = err)
    
    if os.stat(f"{log_dir}/vsearch_{name}.err").st_size != 0:
        raise Exception(f"Something went wrong with VSEARCH\n\nPlease check: {log_dir}/vsearch_{name}.err")
    return
