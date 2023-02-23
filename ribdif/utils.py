#!/usr/bin/env python3
import pandas as pd
import re
import gzip
import shutil
from pathlib import Path
import fileinput
import os
import logging


# =============================================================================
# def modify(file_path):
#     # Open and read the contents of the file
#     with open(file_path, "r") as f_in:
#         contents = f_in.read()
#         
#     # Preform substitutions for unwanted characters
#     contents = re.sub(r"[:,/()=#\x27]", "", contents)
#     contents = re.sub(r"[: ]", "_", contents)
# 
#     # Write the modified contents back to file
#     with open(file_path, "w") as f_out:
#         f_out.write(contents)
#     return
# =============================================================================

# Remove unwanted characters from anywhere is file (should only be in fasta headers)
def modify2(file_path):
    GCF = str(Path(file_path).parent).split(r"/")[-1]
    # Open and read the contents of the file
    with fileinput.input(file_path, inplace = True) as f_in:
        for line in f_in:
            if line.startswith(">"):
                line = re.sub(r"[:,/()=#\x27]", "", line)
                line = re.sub(r"[: ]", "_", line)
                line = f"{line[0]}{GCF}_{line[1:]}" 
                print(line, end = '')
            else:
                print(line, end = '')
    return

# =============================================================================
# def modify3(file_path):
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
def amp_replace(outdir, genus, names, logger):
    for name in names:
        # Read in the summary dataframe
        df_sum = pd.read_csv(f"{outdir}/amplicons/{genus}-{name}.summary", sep = "\t", header = None, names = ["AmpId", "SequenceId", "PositionInSequence", "Length", "Misc"])
        dict_sum = dict(zip(df_sum.AmpId, df_sum.SequenceId)) # Make a dictionary of it
        # Open the temp amplicon file and the final amplicon file
        with open (f"{outdir}/amplicons/{genus}-{name}.temp.amplicons", "r") as f_in, open(f"{outdir}/amplicons/{genus}-{name}.amplicons", "w") as f_out:
           for line in f_in: # loop over lines in the file
               if ">amp" in line: # If it is a fasta header
                   # Replace with origin genome fasta header and then the amp count
                   amp = line.strip().strip(">")
                   line = line.replace(amp, dict_sum[amp] + f"_{amp.strip('amp_')}")
               f_out.write(line)
        Path.unlink(f"{outdir}/amplicons/{genus}-{name}.temp.amplicons") # remove temp file
        # Check if the primer resulted in any amplification and if not remove the file
        if os.stat(f"{outdir}/amplicons/{genus}-{name}.amplicons").st_size == 0:
            Path.unlink(f"{outdir}/amplicons/{genus}-{name}.amplicons")
            logger.info(f"{name} primer resulted in no amplification and will be excluded from further analysis. Are you sure the primer is correct?\n")
            names.remove(name)
    return names

def pairwise_to_csv(pairwise_match, gcf_species, outdir, genus, name):
    pairwise_save_df = pd.DataFrame(pairwise_match, index = gcf_species.values())
    pairwise_save_df.to_csv(f"{outdir}/amplicons/{genus}-{name}_confusion.csv", sep = ",", index = True)
    return

def cleaner(outdir):
   shutil.rmtree(f"{outdir}/refseq", ignore_errors = True)
   return
    
    