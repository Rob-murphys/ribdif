#!/usr/bin/env python3
import pandas as pd
import re
import gzip
import shutil
from pathlib import Path
import fileinput
import os
import logging
import chardet


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
                current_sp = line.split("_")[5]
            else:
                print(line, end = '')
    return current_sp
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

# open and check species of the genome
def sp_check(file):
    with open(file, "r") as f_in:
        line = f_in.readline()
        current_sp = line.split("_")[5]
    return current_sp

# Un gzip the downloaded NCBI genomes
def decompress(file_path):
    # Open the .gz file and decompress it
    with gzip.open(file_path, "rb") as f_in, open(file_path[:-3], "wb") as f_out:
        shutil.copyfileobj(f_in, f_out) # copy
    return

# Rename amplicon fasta headers to origin contig
def amp_replace(outdir, genus, names, logger):
    for name in names[:]:
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
        Path.unlink(Path(f"{outdir}/amplicons/{genus}-{name}.temp.amplicons")) # remove temp file
        # Check if the primer resulted in any amplification and if not remove the file
        if os.stat(f"{outdir}/amplicons/{genus}-{name}.amplicons").st_size == 0:
            Path.unlink(Path(f"{outdir}/amplicons/{genus}-{name}.amplicons"))
            logger.info(f"{name} primer resulted in no amplification and will be excluded from further analysis. Are you sure the primer is correct?\n")
            names.remove(name)
    return names

def pairwise_to_csv(pairwise_match, gcf_species, outdir, genus, name):
    pairwise_save_df = pd.DataFrame(pairwise_match, index = gcf_species.values())
    pairwise_save_df.to_csv(f"{outdir}/amplicons/{genus}-{name}_confusion.csv", sep = ",", index = True)
    return

def detect_encode(file):
    with open(file, "rb") as f_in:
        return chardet.detect(f_in.readline())["encoding"]
    
def primer_check(primer_file, logger):
    encoding = detect_encode(primer_file)
    with open(primer_file, "r", encoding = encoding) as f_in:
        for line in f_in:
            if len(line.strip().split("\t")) != 4:
                logger.info(f"""Your primer file in in the incorrect format at the following line:
                            \t{line}
                            Please follow the 'name | forward | reverse | expected amplicon length' tab seperated format
                            Check the default.primers file or the README for more guidance""")
                return True
    return False

# Copy the user defined genomes to the output directory following ncbi-genome.download structure
def own_genomes_copy(dir_path, outdir, domain, logger):
    target_dir = f"{outdir}/refseq/{domain}" # Target directory
    Path.mkdir(Path(target_dir), exist_ok = True, parents = True) # Creating directory
    logger.info(f"Copying your genomes to {target_dir}\n\n")
    file_count = 0 # Initiate genome count
    # Looping over whole directory
    for file in Path(dir_path).iterdir():
        if file.is_file(): # if item is a file
            Path.mkdir(Path(f"{target_dir}/GCF_{file.stem.replace('_', '-')}"))
            shutil.copy(file, f"{target_dir}/GCF_{file.stem}/{file.stem}.fna") # copy it replacing the file extension with '.fna'
            file_count += 1 # incriment file count
    return target_dir, file_count

# Decompress user defined files if they need it
def own_genomes_gzip(new_dir_path):
    for file in Path(new_dir_path).glob("*/*.fna"):
        if file.suffix == ".gz":
            decompress(file)
    return
    
# Rename the fasta headers of the file inplace
def own_genomes_rename(new_dir_path, logger):
    NZ_count = 1 # abritrary GCF
    logger.info("Changing the copied genomes headers to conform with NCBI format\n\n")
    # Loop throuh all user fna files in new directory
    for file in Path(new_dir_path).glob("*/*.fna"):
        with fileinput.input(file, inplace = True) as f_in: # open the file for inplace editing
            for line in f_in: # Loop through lines
                if line.startswith(">"): # if it is a fasta header
                    # Generate new fasta header    
                    genus = file.stem.replace("_", "-")
                    line = f">GCF_{genus}.1_NZ_CP{NZ_count}.1_{genus}_sp._placeholder\n" # generate random GCF
                    print(line, end = '')
                    NZ_count += 1
                else:
                    print(line, end = '') # else print the line (which would be actual sequence data)
    return


    