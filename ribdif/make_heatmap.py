# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from matplotlib.backends.backend_pdf import PdfPages


# Import and clean the cluster file
def uc_cleaner(outdir, genus, name):
    
    uc_path = f"{outdir}/amplicons/{genus}-{name}.uc"
    
    # read in cluster file
    uc_df = pd.read_csv(uc_path, sep = "\t", header = None)
    # remove cluster summaries
    uc_df_clean = uc_df.loc[uc_df[0] != "C"].copy() # Explicitly setting as a copy to avoid SettingWithCopyWarning warnings
    
    # Create new columns
    uc_df_clean.loc[:, "GCF"] = ["_".join(i.split("_")[:2]) for i in uc_df_clean[8]] # GCF column
    uc_df_clean.loc[:, "Species"] = [i.split("_")[4] for i in uc_df_clean[8]]  # Species column
    
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
    cluster_count = len(uc_df_clean[1].unique()) # getting a count of how many cluster there are
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
        
        for m in clusMatchGCF:
            gcf_index = np.where(all_gcfs==m)[0][0]
            pairwise_match[gcf][gcf_index] = 1
    return pairwise_match

def heatmap_meta(gcf_species):
    # Turn gcf species cross dictionary into series
    species_series = pd.Series(gcf_species)
    
    # Generate colour palet from this
    palette  = sns.color_palette("husl", len(species_series.unique()))
    species_palette = dict(zip(species_series.unique(), palette))
    row_palette = species_series.map(species_palette)
    return row_palette, species_series


def cluster_heatmap(cluster_dict, row_palette, species_series):
    # Turn cluster dictionary into a dataframe and transpose
    cluster_df = pd.DataFrame.from_dict(cluster_dict).transpose()
    
    # Generate clustering on binary data
    row_clus = scipy.cluster.hierarchy.linkage(np.where(cluster_df > 0, 1, 0), method = "ward")
    col_clus = scipy.cluster.hierarchy.linkage(np.where(cluster_df.transpose() > 0, 1, 0), method = "ward")
    
    # Clustering heatmap
    plot_clus = sns.clustermap(cluster_df, standard_scale = None, 
                   row_linkage = row_clus, 
                   col_linkage = col_clus,
                   yticklabels = species_series,
                   row_colors = row_palette)
    
    # Change plot configeration #
    box_heatmap = plot_clus.ax_heatmap.get_position() # Get heatmaps position
    
    # Move row colours to left side
    ax_row_colors = plot_clus.ax_row_colors # get the axis
    box_cols = ax_row_colors.get_position() # recover its position
    ax_row_colors.set_position([box_heatmap.max[0], box_cols.y0, box_cols.width*2, box_cols.height]) # plot new position
    
    # Move dendogram to the left a bit
    ax_row_dendogram = plot_clus.ax_row_dendrogram
    box_dendo = ax_row_dendogram.get_position()
    ax_row_dendogram.set_position([box_dendo.x0+0.026, box_dendo.y0, box_dendo.width, box_dendo.height])
    
    return plot_clus

def pairwise_heatmap(pairwise_match, row_palette, species_series):
    # Turn pairwise dictionary into dataframe
    pairwise_df = pd.DataFrame(pairwise_match, index = pairwise_match.keys())
    
    # Clustering heatmap
    plot_dendo = sns.clustermap(pairwise_df, standard_scale = None,
                   row_colors = row_palette,
                   yticklabels = species_series,
                   method = "ward")
    
    # Change plot configeration #
    box_heatmap = plot_dendo.ax_heatmap.get_position() # Get heatmaps position
    
    # Move row colours to left side
    ax_row_colors = plot_dendo.ax_row_colors # get the axis
    box_cols = ax_row_colors.get_position() # recover its position
    ax_row_colors.set_position([box_heatmap.max[0], box_cols.y0, box_cols.width*2, box_cols.height]) # plot new position
    
    # Move dendogram to the left a bit
    ax_row_dendogram = plot_dendo.ax_row_dendrogram
    box_dendo = ax_row_dendogram.get_position()
    ax_row_dendogram.set_position([box_dendo.x0+0.026, box_dendo.y0, box_dendo.width, box_dendo.height])
    
    return plot_dendo, pairwise_df

def pdf_save(plots, outdir, genus, name):
    with PdfPages(f"{outdir}/amplicons/{genus}-{name}.pdf") as pdf_pages:
        for plot in plots:
            pdf_pages.savefig(plot.fig)
    return
        

    
            
