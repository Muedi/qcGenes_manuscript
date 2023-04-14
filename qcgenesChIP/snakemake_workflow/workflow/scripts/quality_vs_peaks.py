#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 19:57:40 2020

@author: jannik
"""
# %%

import warnings
warnings.filterwarnings("ignore")
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%

XSMALL_SIZE = 14
SMALL_SIZE = 20
MEDIUM_SIZE = 30
BIGGER_SIZE = 36

plt.rc('font', size=XSMALL_SIZE)         # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)    # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

# %% 
if len(sys.argv) < 2:
    raise ValueError("Enter file paths of binned peak .bed files")

# retrieve input
filepaths = sys.argv[1:-1]
geo_id = sys.argv[-1]
dataframes = [pd.read_csv(file) for file in filepaths]

# %%

def save_bar_plot(df, agg, geo_id, bin_labels, agg_type = "",):
    if agg == "sum":
        ylab = "total peak count"
        row_id = "row_" + agg
        df[row_id] = df.loc[:,df.columns.str.startswith("SRR")].sum(axis=1).astype(int)
        grouped = df.groupby("bin").sum().reindex(bin_labels).reset_index()
    elif agg == "mean":
        ylab = "mean enrichment of bin " + agg_type
        row_id = "row_" + agg
        df[row_id] = df.loc[:,df.columns.str.startswith("SRR")].mean(axis=1)
        grouped = df.groupby("bin").mean().reindex(bin_labels).reset_index()
        grouped[row_id] = grouped[row_id].round(2)
    else:
        raise("invalid aggregation")
    p_values = df.groupby("bin").mean().reindex(bin_labels).reset_index().loc[:,"p"].round(2)
    fig, ax = plt.subplots(figsize=(18,12))
    plt.xticks(rotation=30)
    ax.set_xlabel("correlation bin")
    ax.set_ylabel(ylab)
    ax.set_title(ylab + " vs correlation bin for " + geo_id, pad = 20)
    ax.bar(x=grouped["bin"], height=grouped[row_id], width=0.8, align="center")
    max_y = max(grouped[row_id])
    ax.set_ylim(0, max_y + max_y / 10)
    for index, row in grouped.iterrows():
        plt.text(row["bin"], row[row_id] + np.ceil(max_y / 80), "n= " + str(row[row_id]), color='black', ha="center")
        plt.text(row["bin"], row[row_id] + np.ceil(max_y / 20), "p = " + str(p_values[index]), color="black", ha="center")
    plt.savefig(os.path.join("./output/plots", geo_id, "quality_vs_" + ylab).replace(" ", "_"))


# %%
bin_labels = ["[-1|-0.8)","[-0.8|-0.6)", "[-0.6|-0.4)", "[-0.4|-0.2)", "[-0.2|0)",
            "[0|0.2)","[0.2|0.4)", "[0.4|0.6)","[0.6|0.8)", "[0.8|1)"]

for df, file in zip(dataframes, filepaths):
    df.name = os.path.split(file)[1].split(".")[0]
    df.dropna(thresh=5, inplace=True)

# %%
for df in dataframes:
    if df.name == "peak_count":
        save_bar_plot(df, "sum", geo_id, bin_labels)
    else:
        save_bar_plot(df, "mean", geo_id, bin_labels, df.name.split("_")[1])
