#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from matplotlib.backends.backend_pdf import PdfPages
import networkx as nx

def heatmap_meta(gcf_species):
    # Turn gcf species cross dictionary into series
    species_series = pd.Series(gcf_species)
    
    # Generate colour palet from this
    palette  = sns.color_palette("husl", len(species_series.unique()))
    species_palette = dict(zip(species_series.unique(), palette))
    row_palette = species_series.map(species_palette)
    return row_palette, species_series, species_palette

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
                   xticklabels = 1,
                   row_colors = row_palette)
    

    plot_clus.ax_heatmap.tick_params(axis='x', labelsize=8)
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
    
    return plot_clus, cluster_df

def pairwise_heatmap(pairwise_match, row_palette, species_series):
    # Turn pairwise dictionary into dataframe
    pairwise_df = pd.DataFrame(pairwise_match, index = pairwise_match.keys())
    
    # Clustering heatmap
    plot_dendo = sns.clustermap(pairwise_df, standard_scale = None,
                   row_colors = row_palette,
                   yticklabels = species_series,
                   xticklabels = 1,
                   method = "ward")
    
    plot_dendo.ax_heatmap.tick_params(axis='x', labelsize=8)
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

def pdf_save(plot_clus, plot_dendo, outdir, genus, name):
    plots = [plot_clus, plot_dendo]
    with PdfPages(f"{outdir}/amplicons/{genus}-{name}_heatmaps.pdf") as pdf_pages:
        for plot in plots:
            pdf_pages.savefig(plot.fig)
    return

def create_graph(pairwise_df):
    # copy the pairwise df so we can make it an adjacency df
    adjacency_df = pairwise_df.copy()
    np.fill_diagonal(adjacency_df.values, 0)# Ensure the diagonal is 0s otherwise we get looped networks
    
    # Create the graph
    graph = nx.from_pandas_adjacency(adjacency_df)
    graph.remove_nodes_from(list(nx.isolates(graph)))# Remove singletons
    #pos = nx.spring_layout(graph, k = 0.15, iterations = 15)
    #nx.draw(graph, pos, node_size = 100, font_size = 8)
    #network_plot = nx.draw_networkx(graph, pos, node_size = 100, font_size = 8)
    
    # Create a list of sub graphs
    graph_subs = [graph.subgraph(c) for c in nx.connected_components(graph)]
    n_subplots = len(graph_subs) # determine the number of sub graps
    return graph_subs, n_subplots

def draw_graphs(graph_subs, n_subplots, species_palette, row_palette, outdir, genus, name):

    
    # Determine the number of rows and columns needed to keep the figure a square
    n_rows = int(np.ceil(np.sqrt(n_subplots))) # rounding the sqrt of subplot number
    n_cols = int(np.ceil(n_subplots / n_rows))
    
    # Create a figure and subplots using a loop
    fig, axs = plt.subplots(n_rows, n_cols, figsize = (n_cols * 5, n_rows * 5))
    axs = axs.flatten()
    for i, graph in enumerate(graph_subs):
        node_cols = [row_palette[s] for s in graph.nodes] # get colours for current nodes
        nx.draw(graph, ax=axs[i], node_color = node_cols) # draw specific graph
        
    # Make any unused subplots blank
    for i in range(n_subplots, n_rows * n_cols):
        axs[i].axis('off')
        axs[i].set_xlim([0, 1])
        axs[i].set_ylim([0, 1])
    
    # Create a dumy legend
    for k, v in species_palette.items():
        plt.scatter([],[], color = v, label = k)
    plt.legend(ncol = n_cols)
    
    # Save the figure
    plt.savefig(f"{outdir}/amplicons/{genus}-{name}_graphs.pdf", bbox_inches = "tight")
    return

