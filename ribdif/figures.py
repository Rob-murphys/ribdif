#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
#import scipy
import fastcluster
from matplotlib.backends.backend_pdf import PdfPages
import networkx as nx
from itertools import combinations

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
    #row_count = len(cluster_df.index)
    # Generate clustering on binary data
    #row_clus = scipy.cluster.hierarchy.linkage(np.where(cluster_df > 0, 1, 0), method = "ward")
    #col_clus = scipy.cluster.hierarchy.linkage(np.where(cluster_df.transpose() > 0, 1, 0), method = "ward")
    row_clus = fastcluster.ward(np.where(cluster_df > 0, 1, 0))
    col_clus = fastcluster.ward(np.where(cluster_df.transpose() > 0, 1, 0))
    #plot_size = (16, 16) if row_count < 50 else ((row_count*0.2), (row_count*0.2))
  
    # Clustering heatmap
    plot_clus = sns.clustermap(cluster_df, standard_scale = None, 
                   row_linkage = row_clus, 
                   col_linkage = col_clus,
                   yticklabels = species_series,
                   xticklabels = 1,
                   row_colors = row_palette,
                   linecolor = "#bcc2cc",
                   linewidths = 0.1,
                   cmap = sns.cm.rocket_r)   

    return plot_clus, cluster_df

def pairwise_heatmap(pairwise_match, row_palette, species_series):
    # Turn pairwise dictionary into dataframe
    pairwise_df = pd.DataFrame(pairwise_match, index = pairwise_match.keys())
    
    # Clustering heatmap
    plot_dendo = sns.clustermap(pairwise_df, standard_scale = None,
                   row_colors = row_palette,
                   yticklabels = species_series,
                   xticklabels = 1,
                   method = "ward",
                   linecolor = "#bcc2cc",
                   linewidths = 0.1,
                   cbar_pos = None,
                   cmap = sns.cm.rocket_r)

    return plot_dendo, pairwise_df

def figure_fix(plot):
    
    plot.ax_heatmap.tick_params(axis='x', labelsize=8)
    plot.ax_row_colors.tick_params(bottom = False) # remove tickmark under row colours
    
    # Get info to change row_cols width
    yticklabels = plot.ax_heatmap.get_yticklabels() # getting the yticks
    max_width = max([label.get_window_extent().width for label in yticklabels]) # getting the max width of them
    fig_width = plot.fig.get_figwidth() * plot.fig.dpi # Getting figure width
    fraction_of_fig_width = (max_width / fig_width) + 0.01 # get fraction of figure of max ytick plus some margin
    
    # Change plot configeration #
    box_heatmap = plot.ax_heatmap.get_position() # Get heatmaps position
    
    # Move row colours to left side
    ax_row_colors = plot.ax_row_colors # get the axis
    box_cols = ax_row_colors.get_position() # recover its position
    ax_row_colors.set_position([box_heatmap.max[0], box_cols.y0, fraction_of_fig_width, box_cols.height]) # plot new position
    
    # Move dendogram to the left a bit
    ax_row_dendogram = plot.ax_row_dendrogram
    box_dendo = ax_row_dendogram.get_position()
    ax_row_dendogram.set_position([box_dendo.x0+0.026, box_dendo.y0, box_dendo.width, box_dendo.height])
    
    return plot

def pdf_save(plot_clus, plot_dendo, outdir, genus, name):
    plots = [plot_clus, plot_dendo]
    with PdfPages(f"{outdir}/figures/{genus}-{name}_heatmaps.pdf") as pdf_pages:
        for plot in plots:
            pdf_pages.savefig(plot.fig)
    return

def create_adjacency(pairwise_df, cluster_df):
    # Create new dataframe filled with the zame index and column as pairwise_df but filled with 0s
    adjacency_df = pd.DataFrame(np.zeros((len(pairwise_df), len(pairwise_df.columns))), index=pairwise_df.index, columns=pairwise_df.columns)
    # Get all possible pairs of indexs
    pairs = list(combinations(cluster_df.index,2))
    # For each pair mask where they both belong to the same cluster
    for i,j in pairs:
        mask = (cluster_df.loc[i] > 0) & (cluster_df.loc[j] > 0)
        adjacency_df[j][i] = mask.sum()
        adjacency_df[i][j] = mask.sum()
        
    return adjacency_df

def create_graph(adjacency_df):
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

    
    # Determine the number of rows and columns needed to keep the figure a rectangle
    n_rows = int(np.ceil(np.sqrt(n_subplots))) # rounding the sqrt of subplot number
    n_cols = int(np.ceil(n_subplots / n_rows))
    
    # Create a figure and subplots using a loop
    fig, axs = plt.subplots(n_rows, n_cols, figsize = (n_cols * 5, n_rows * 5))
    if n_rows > 1 and n_cols > 1:
        axs = axs.flatten()
    for i, graph in enumerate(graph_subs):
        #pos = nx.spring_layout(graph)
        #edge_labels = {(n1, n2): graph[n1][n2]['weight'] for n1, n2 in graph.edges()}
        weights = nx.get_edge_attributes(graph, "weight")
        node_cols = [row_palette[s] for s in graph.nodes] # get colours for current nodes
        try: 
            nx.draw(graph, ax=axs[i], node_color = node_cols, width = list(weights.values())) # draw specific graph
        except TypeError:
            nx.draw(graph, ax=axs, node_color = node_cols, width = list(weights.values()))
        #nx.draw_networkx_edge_labels(graph, pos, edge_labels = edge_labels)
        
    # Make any unused subplots blank
    for i in range(n_subplots, n_rows * n_cols):
        axs[i].axis('off')
        axs[i].set_xlim([0, 1])
        axs[i].set_ylim([0, 1])

    
    # Create a dumy legend
    for k, v in species_palette.items():
        plt.scatter([],[], color = v, label = k)
    plt.legend(ncol = n_cols)
    #plt.savefig("test_weights.pdf", bbox_inches="tight")
    # Save the figure
    plt.savefig(f"{outdir}/figures/{genus}-{name}_graphs.pdf", bbox_inches = "tight")
    return

