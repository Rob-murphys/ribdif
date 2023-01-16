# -*- coding: utf-8 -*-
import pandas as pd
import re
import gzip
import shutil
from pathlib import Path

# Remove unwanted characters from anywhere is file (should only be in fasta headers)
def modify(file_path):
    # Open and read the contents of the file
    with open(file_path, "r") as f_in:
        contents = f_in.read()
        
    # Preform substitutions for unwanted characters
    contents = re.sub(r"[:,/()=#\x27]", "", contents)
    contents = re.sub(r"[: ]", "_", contents)

    # Write the modified contents back to file
    with open(file_path, "w") as f_out:
        f_out.write(contents)
    return

# =============================================================================
# def modify2(file_path):
#     # Open and read the contents of the file
#     with open(file_path, "rb") as f_in:
#         contentsB = f_in.read()
#         
#     # Preform substitutions for unwanted characters
#     table = bytes().maketrans(b" ", b"_")
#     contentsS = contentsB.translate(None, b":,/()=#").translate(table).decode()
# 
#     # Write the modified contents back to file
#     with open(file_path, "w") as f_out:
#         f_out.write(contentsS)
# =============================================================================

# Un gzip the downloaded NCBI genomes
def decompress(file_path):
    # Open the .gz file and decompress it
    with gzip.open(file_path, "rb") as f_in, open(file_path[:-3], "wb") as f_out:
        shutil.copyfileobj(f_in, f_out) # copy
    return

# Rename amplicon fasta headers to origin contig
def amp_replace(outdir, genus, name):
    df_sum = pd.read_csv(f"{outdir}/amplicons/{genus}-{name}.summary", sep = "\t", header = None, names = ["AmpId", "SequenceId", "PositionInSequence", "Length", "Misc"])
    dict_sum = dict(zip(df_sum.AmpId, df_sum.SequenceId))
    with open (f"{outdir}/amplicons/{genus}-{name}.temp.amplicons", "r") as f_in, open(f"{outdir}/amplicons/{genus}-{name}.amplicons", "w") as f_out:
       for line in f_in:
           if ">amp" in line:
               amp = line.strip().strip(">")
               line = line.replace(amp, dict_sum[amp] + f"_{amp.strip('amp_')}")
           f_out.write(line)
    Path.unlink(f"{outdir}/amplicons/{genus}-{name}.temp.amplicons")
    return