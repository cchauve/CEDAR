import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Input CSV file
# Format: tree,SPR_dist,order,dist_type,dist_value
INPUT_FILE = sys.argv[1]
OUT_DIR = sys.argv[2]
SIZE_TREES = int(sys.argv[3])
NB_TREES = int(sys.argv[4])
NB_ORDERS = int(sys.argv[5])

DISTANCES_DF = pd.read_csv(INPUT_FILE,header=0)

# Figure 1: SPR versus normalized RF
def SPR_vs_RF(distances, outfile):
    """
    Input:
    - distances: dataframe with columns tree,SPR_dist,order,dist_type,dist_value
    - outfile: PNG file with a plot of SPR (x-axis) versus RF distribution (y-axis)
    """
    RF = distances.loc[distances["dist_type"] == "RF"]
    SPR_dist = sorted(RF["SPR_dist"].unique())
    RF_values = {SPR: [] for SPR in SPR_dist}
    for _,row in RF.iterrows():
        RF_values[row["SPR_dist"]].append(row["dist_value"] / (2*SIZE_TREES))

    fig, ax = plt.subplots()
    ax.boxplot(RF_values.values())
    ax.set_xticklabels(RF_values.keys())
    plt.xlabel("SPR dist.")
    plt.ylabel("Normalizd RF dist.")
    plt.title(f"SPR vs RF over {NB_TREES} trees")
    plt.savefig(outfile)


# Figure 2: SPR versus mean normalized HOP
def SPR_vs_mean_HOP(distances, outfile):
    """
    Input:
    - distances: dataframe with columns tree,SPR_dist,order,dist_type,dist_value
    - outfile: PNG file with a plot of SPR (x-axis) versus mean (over orders) HOP distribution (y-axis)
    """
    HOP = distances.loc[distances["dist_type"] == "HOP"]
    SPR_dist = sorted(HOP["SPR_dist"].unique())
    trees = sorted(HOP["tree"].unique())
    HOP_values = {(SPR,tree): [] for SPR in SPR_dist for tree in trees}
    for _,row in HOP.iterrows():
        HOP_values[(row["SPR_dist"],row["tree"])].append(row["dist_value"] / SIZE_TREES)
    HOP_mean = {
        SPR: [np.mean(HOP_values[(SPR,tree)]) for tree in trees]
        for SPR in SPR_dist
    }
            
    fig, ax = plt.subplots()
    ax.boxplot(HOP_mean.values())
    ax.set_xticklabels(HOP_mean.keys())
    plt.xlabel("SPR dist.")
    plt.ylabel("Mean normalized HOP dist.")
    plt.title(f"SPR vs mean HOP over {NB_TREES} trees and {NB_ORDERS} random orders")
    plt.savefig(outfile)


# Figure 3: SPR versus standard deviation of normalized HOP
def SPR_vs_std_HOP(distances, outfile):
    """
    Input:
    - distances: dataframe with columns tree,SPR_dist,order,dist_type,dist_value
    - outfile: PNG file with a plot of SPR (x-axis) versus stdv (over orders) HOP distribution (y-axis)
    """
    HOP = distances.loc[distances["dist_type"] == "HOP"]
    SPR_dist = sorted(HOP["SPR_dist"].unique())
    trees = sorted(HOP["tree"].unique())
    HOP_values = {(SPR,tree): [] for SPR in SPR_dist for tree in trees}
    for _,row in HOP.iterrows():
        HOP_values[(row["SPR_dist"],row["tree"])].append(row["dist_value"] / SIZE_TREES)
    HOP_std = {
        SPR: [np.std(HOP_values[(SPR,tree)]) for tree in trees]
        for SPR in SPR_dist
    }
            
    fig, ax = plt.subplots()
    ax.boxplot(HOP_std.values())
    ax.set_xticklabels(HOP_std.keys())
    plt.xlabel("SPR dist.")
    plt.ylabel("Normalized HOP dist. standard deviation")
    plt.title(f"SPR vs HOP standard deviation over {NB_TREES} trees and {NB_ORDERS} random orders")
    plt.savefig(outfile)


# Figure 4: SPR versus normalized HOP
def SPR_vs_HOP(distances, outfile):
    """
    Input:
    - distances: dataframe with columns tree,SPR_dist,order,dist_type,dist_value
    - outfile: PNG file with a plot of SPR (x-axis) versus HOP distribution (y-axis)
    """
    HOP = distances.loc[distances["dist_type"] == "HOP"]
    SPR_dist = sorted(HOP["SPR_dist"].unique())
    HOP_values = {SPR: [] for SPR in SPR_dist}
    for _,row in HOP.iterrows():
        HOP_values[row["SPR_dist"]].append(row["dist_value"] / SIZE_TREES)
            
    fig, ax = plt.subplots()
    ax.boxplot(HOP_values.values())
    ax.set_xticklabels(HOP_values.keys())
    plt.xlabel("SPR dist.")
    plt.ylabel("Normalized HOP dist. ")
    plt.title(f"SPR vs HOP over {NB_TREES} trees and {NB_ORDERS} random orders")
    plt.savefig(outfile)
    
SPR_vs_RF(DISTANCES_DF,os.path.join(OUT_DIR,"exp4_SPR_vs_RF.png"))
SPR_vs_mean_HOP(DISTANCES_DF,os.path.join(OUT_DIR,"exp4_SPR_vs_mean_HOP.png"))
SPR_vs_std_HOP(DISTANCES_DF,os.path.join(OUT_DIR,"exp4_SPR_vs_std_HOP.png"))
SPR_vs_HOP(DISTANCES_DF,os.path.join(OUT_DIR,"exp4_SPR_vs_HOP.png"))
