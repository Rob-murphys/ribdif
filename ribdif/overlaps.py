#!/usr/bin/env python3

import pandas as pd
import numpy as np



# Import and clean the cluster file
def uc_cleaner(outdir, genus, name):
    
    uc_path = f"{outdir}/amplicons/{name}/{genus}-{name}.uc"
    
    # read in cluster file
    uc_df = pd.read_csv(uc_path, sep = "\t", header = None)
    # remove cluster summaries
    uc_df_clean = uc_df.loc[uc_df[0] != "C"].copy() # Explicitly setting as a copy to avoid SettingWithCopyWarning warnings
    
    # Create new columns
    uc_df_clean.loc[:, "GCF"] = ["_".join(i.split("_")[:2]) for i in uc_df_clean[8]] # GCF column
    uc_df_clean.loc[:, "Species"] = [i.split("_")[5] for i in uc_df_clean[8]]  # Species column
    
    # Sorting dataframe
    #uc_df_clean = uc_df_clean.sort_values(by = ["GCF"]) # check if this matters later on
    uc_df_clean = uc_df_clean.sort_values(by = ["Species", "GCF"]) # Mikael orderd by both? why?
    
    # Moving unclassified species to the bottom
    if 'sp.' in uc_df_clean.Species.unique():
        uc_df_clean = pd.concat([uc_df_clean[uc_df_clean.Species != "sp."], uc_df_clean[uc_df_clean.Species == "sp."]])
    
    # Storing information needed later
    uc_dict_clean = uc_df_clean.T.to_dict("list") # Chaging the dataframe to a dictionary
    gcf_species = dict(zip(uc_df_clean.GCF, uc_df_clean.Species)) # dicttionary of GCF to species
    all_gcfs = uc_df_clean.GCF.unique() # getting an array of all GCFs
    #cluster_count = len(uc_df_clean[1].unique()) # getting a count of how many cluster there are
    cluster_count = max(uc_df_clean[1]) +1
    return all_gcfs, uc_dict_clean, gcf_species, cluster_count
    

    
# Generate a dictionary (that will become a matrix) of GCF cluster membership  
def cluster_matrix(all_gcfs, uc_dict_clean, cluster_count):
    # Creating a dictionary of lists where keys are GCFs and values are lists the length of cluster present
    cluster_dict = {key: [0]*cluster_count for key in all_gcfs}
    
    # Looping over unique GCFs
    for gcf in all_gcfs:
        # Fishing out current GCFs and converting to a dictionary of lists  where rows indicies in the df are keys
        current_GCF = dict([(k, v) for k, v in uc_dict_clean.items() if v[10] == gcf])
        #current_GCF = uc_df_clean[uc_df_clean.GCF == gcf].T.to_dict("list") # if using dataframe
        # Loop through fished out GCFs clusters
        for r in current_GCF.values():
            # For the current cluster for a given GCF add 1 to the count of time the given GCF is present in the given cluster
            cluster_dict[r[10]][r[1]] = cluster_dict[r[10]][r[1]] + 1
    return cluster_dict



# Find all species overlap in the cluster_dict
def species_overlap(cluster_dict, cluster_count, gcf_species):
    combinations = []
    # Loop through all clusters
    for i in range(cluster_count):
        # Get all GCF values that have more that one amplicon beloning to the current cluster
        keys = [k for k, v in cluster_dict.items() if v[i] > 1]
        # Get unique species names asigned for the GCFs recoverd above
        cluster_species = set([gcf_species[gcf] for gcf in keys])
        # If more than one species belongs to this cluster then save that combination
        if len(cluster_species) > 1: 
            combinations.append("/".join(cluster_species))
    return combinations


# Find all GCF overlaps in the cluster dictionary
def gcf_overlaps(all_gcfs, uc_dict_clean, gcf_species):
    # Dictionary of lists of length equal to all gcfs populated with list length of all gcfs
    pairwise_match = {key: [0]*len(all_gcfs) for key in all_gcfs}
    for gcf in all_gcfs:
        # Fish out current GCF
        current_GCF = dict([(k, v) for k, v in uc_dict_clean.items() if v[10] == gcf])
        # current_GCF = uc_df_clean[uc_df_clean.GCF == gcf].T.to_dict("list") # if using dataframe
        # get species of current GCF
        #species = gcf_species[gcf]
        # get unique cluster this GCF belongs to
        clusters = set([v[1] for v in current_GCF.values()])
        # Getting a list of other GCFs that are members of the clusters the current GCF belongs to
        clusMatchGCF = set([v[10] for v in uc_dict_clean.values() if v[1] in clusters])
        
        pairwise_match[gcf] = [1 if m in clusMatchGCF else 0 for m in all_gcfs]
        #for m in clusMatchGCF:
        #    gcf_index = np.where(all_gcfs==m)[0][0]
        #    pairwise_match[gcf][gcf_index] = 1
    return pairwise_match


def overlap_report(combinations, gcf_species, cluster_df, genus, name, outdir, logger, shannon_div, unique_species, all_species, genome_count, user):
    
    total_named_genomes = len([s for s in all_species if s != "sp."])
    total_unamed_genomes = genome_count - total_named_genomes
    total_unique_species = len(unique_species)
    
    amplify_count = len(cluster_df.index)
    amplify_named_genomes = len([s for s in gcf_species.values() if s != "sp."]) # Recovering species names that are not equal so "sp."
    amplify_unamed_genomes = len(gcf_species) - amplify_named_genomes # Subtracting names species from total specpes
    amplify_unique_species = len([i for i in set(gcf_species.values()) if i != "sp."]) # unique names speices
    # turning the cluster_df binary with np.where and summing rows then counting the times rowsum is greater that 1
    multi_allele = sum(np.where(cluster_df > 0, 1, 0).sum(axis = 1) > 1)
    if combinations: # are the unique entries into combinations greater than 0. Not sure we need the set.
        unq_combs = sorted(set(combinations), key = len, reverse = True) # getting unique combinations and sorting by reverse string length
        has_overlap = set([e for c in combinations for e in c.split("/")]) # getting a unique set of species that experienced overlap
        if "sp." in has_overlap:
            has_overlap.remove("sp.")
    else:
        unq_combs = []
        has_overlap = []
        
    count_overlap = len(has_overlap)
    with open(f"{outdir}/{genus}_{name}_overlap_report.txt", "w+") as f_out:
        f_out.write(f"""Summary of {genus} differentiation by {name} amplicons:\n
                    Genomes downloaded: {genome_count}
                    \tWith species name: {total_named_genomes}
                    \tWithout species name: {total_unamed_genomes}
                    \tUnique species names: {total_unique_species}\n
                    Genomes amplified by {name}: {amplify_count}
                    \tWith species name: {amplify_named_genomes}
                    \tWithout species name: {amplify_unamed_genomes}
                    \tUnique species names: {amplify_unique_species}\n
                    {multi_allele} of {amplify_count} ({round(100*multi_allele/amplify_count, 2)}%) genomes that amplified have multiple alleles.
                    {count_overlap} of {amplify_unique_species} ({0 if user else round(100*count_overlap/amplify_unique_species, 2)}%) species that experienced amplification have at least one overlap.\n
                    Total shannon diversity for {name} is: {shannon_div}\n\n""")
        if unq_combs:
            for i in range(len(unq_combs)):
                members = unq_combs[i].split("/")
                f_out.write(f"Group {i}:\t" + ' / '.join(members) + "\n")
        f_out.seek(0)
        logger.info(f_out.read())
        
#[i.replace("sp.", f"sp._{x}") for x, i in enumerate(combinations) if "sp." in i]

# Changing the cluster dictionary into a dataframe for when only a single genome amplified
def single_amp_df(cluster_dict):
    return pd.DataFrame.from_dict(cluster_dict).transpose()