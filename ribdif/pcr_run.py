#!/usr/bin/env python3
import multiprocessing
import subprocess
from pathlib import Path
from itertools import repeat
import shlex
import csv
import fileinput
import logging
import shutil
from ribdif.utils import detect_encode
"""
Implement a producer and consumer setup for writing the pcr output whe multiprocessing: https://stackoverflow.com/questions/11196367/processing-single-file-from-multiple-processes

I removed the counter from all functions for labeling output files as it seems the cocatinated fasta ran pretty fast.
"""

# Spawning the shell call
def call_proc_pcr(infile, primer_path, genus, name, fwd, rvs, length, workingDir, multi): # removed , counter
    # Checking if running on multi processing mode
    if not multi:    
        # Building the command
        command = f"perl {workingDir}/in_silico_PCR.pl -s {infile} -a {fwd} -b {rvs} -l {length} -r -m -i > {primer_path}/{genus}-{name}.summary 2> {primer_path}/{genus}-{name}.temp.amplicons"
        # Passing the command to shell piping the stdout and stderr
        with open(f"{primer_path}/{genus}-{name}.summary", "w") as f_std, open(f"{primer_path}/{genus}-{name}.temp.amplicons", "w") as f_err:
            subprocess.run(shlex.split(command), stdout = f_std, stderr = f_err)
    # If running with multi then output is directed to seperate files for later processing
    else:
        outfile = Path(infile).stem
        command = f"perl {workingDir}/in_silico_PCR.pl -s {infile} -a {fwd} -b {rvs} -l {length} -r -m -i > {primer_path}/{outfile}.summary 2> {primer_path}/{outfile}.temp.amplicons"
        
        with open(f"{primer_path}/{outfile}_{name}.summary", "w") as f_std, open(f"{primer_path}/{outfile}_{name}.amplicons", "w") as f_err:
            subprocess.run(shlex.split(command), stdout = f_std, stderr = f_err)
    return 

# Multithreading the in silico pcr calls
def pcr_parallel_call(outdir, genus, primer_file, workingDir, threads, logger, domain):
    amplicon_dir = Path(f"{outdir}/amplicons") # path to amplicon directory
    amplicon_dir.mkdir(parents = True, exist_ok = True) # making the directory
    multi = True
    logger.info("Generating amplicon sequences\n\n")
    encoding = detect_encode(primer_file) # detecting the encoding of the primer file
    with open(primer_file, "r", encoding = encoding) as f_primer: # opening the primer file
        names = []
        for primer in f_primer: # looping through the lines (and thus primers)
            name, fwd, rvs, length = primer.strip().split("\t") # getting infor about each primer
            primer_path = Path(f"{amplicon_dir}/{name}") # directory for each primer output to be in
            if primer_path.is_dir():
                shutil.rmtree(primer_path, ignore_errors = False)
            primer_path.mkdir(parents = True, exist_ok = False)
            with multiprocessing.Pool(threads) as pool: # spawn the pool # opening the pool
                all_fna = [str(i) for i in list(Path(f"{outdir}/refseq/{domain}/").rglob('*.fna'))] # generate list of files ending in .fna
                #counter = range(len(all_fna))
                longer_length = int((float(length)+(float(length)*0.5)))
                pool.starmap(call_proc_pcr, zip(all_fna, repeat(primer_path), repeat(genus), repeat(name), repeat(fwd), repeat(rvs), repeat(longer_length), repeat(workingDir), repeat(multi))) # removed counter
            names.append(name)
            amplicon_filter(outdir, name, genus, length)
    return names

def pcr_call(infile, outdir, genus, primer_file, workingDir, logger):
    amplicon_dir = Path(f"{outdir}/amplicons")
    amplicon_dir.mkdir(parents = True, exist_ok = True)
    multi = False
    logger.info("#= Generating amplicon sequences =#\n\n")
    encoding = detect_encode(primer_file)
    with open(primer_file, "r", encoding = encoding) as f_primer:
        names = []
        for primer in f_primer:
            name, fwd, rvs, length = primer.strip().split("\t")
            primer_path = Path(f"{amplicon_dir}/{name}")
            if primer_path.is_dir():
                shutil.rmtree(primer_path, ignore_errors = False)
            primer_path.mkdir(parents = True, exist_ok = False)
            call_proc_pcr(infile, primer_path, genus, name, fwd, rvs, int((float(length)+(float(length)*0.5))), workingDir, multi)
            names.append(name)
            amplicon_filter(outdir, name, genus, length)
    return names


# House keeping of multithreaded pcr to ensure same naming convention as though they were PCRed in a from a single concatinated file
def multi_cleaner(outdir, name):
    amp_counter = 1
    total_sum_dict = {}
    # Loop through all amplicon files
    for file in Path(f"{outdir}/amplicons/{name}/").glob(f"*{name}.amplicons"):
        # Open the summary file and turn it into a dictionary
        with open(str(file).replace(".amplicons", ".summary")) as f_in:
             rows = ( line.strip().split('\t') for line in f_in )
             sum_dict = { row[0]:row[1:] for row in rows if "AmpId" not in row} # ignoring the header row
        with fileinput.input(file, inplace = True) as amp_in: # In place changing the amplicon file
            for line in amp_in:
                if ">amp_" in line: # if this is in the line
                    total_sum_dict[f"amp_{amp_counter}"] = sum_dict[line.strip().strip(">")]
                    #sum_dict[f"amp_{amp_counter}"] = sum_dict.pop(line.strip().strip(">")) # Update the respective row in the summary file dictionary
                    line = f">amp_{amp_counter}\n" # Change the line
                    print(line, end = '') # write it in place to the file
                    amp_counter += 1 # incrament by one
                else:
                    print(line, end = '') # Writing the nucleotide lines
        #total_sum_dict.update(sum_dict) # appending to the total dict
    return total_sum_dict
                      
# Just concatinating all corrected amplicon files  
def amplicon_cat(outdir, genus, name):
    with open(f"{outdir}/amplicons/{name}/{genus}-{name}.temp.amplicons", "w") as amp_out:
        for file in Path(f"{outdir}/amplicons/{name}/").glob(f"*{name}.amplicons"):
            with open(file, "r") as f_in:
                f_read = f_in.read()
                amp_out.write(f_read)
    return

# Writing the total summary dictionary to a tsv
def sum_dict_write(outdir, genus, name, total_sum_dict):
    with open(f"{outdir}/amplicons/{name}/{genus}-{name}.summary", "w", newline = '') as sum_out:
        writer = csv.writer(sum_out, delimiter = "\t") # generate a  csv writer
        writer.writerow(["AmpId", "SequenceId", "PositionInSequence", "Length", "Misc"]) # write the headers
        for key in total_sum_dict.keys():
            total_sum_dict[key].insert(0, key) # append the value with the key
            writer.writerow(total_sum_dict[key]) # write each value to file
    return

# Removes any amplicon that is not at least 90% of the expected length
def amplicon_filter(outdir, name, genus, length):
    for file in Path(f"{outdir}/amplicons/{name}/").glob(f"{genus}-{name}.amplicons"):
        with open(file, "r") as f_in:
            lines = f_in.readlines()
        modified_lines = [line for i, line in enumerate(lines) 
                          if (line.startswith(">") and len(lines[i+1]) >= (length * 0.9)) 
                          or (len(line) >= (length * 0.9) and lines[i-1].startswith(">"))]
        with open(file, "w") as f_in:
            f_in.writelines(modified_lines)