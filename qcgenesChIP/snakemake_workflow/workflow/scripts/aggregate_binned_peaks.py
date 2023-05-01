#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 19:57:40 2020

@author: jannik
"""
import os
import sys
import pandas as pd
import numpy as np
from scipy import stats

# %%
if len(sys.argv) < 2:
    raise ValueError("Enter file paths of binned peak .bed files")
# set NA threshold (only rows with min number of non-NA values are allowed)
na_threshold = 3

# retrieve input
file_paths = sorted(sys.argv[1:-5])
scores_path = sys.argv[-5]
scores = pd.read_csv(scores_path, sep="\t", header=None).iloc[:,1].to_numpy()
sra_ids = [os.path.split(f)[1].split("_")[0] for f in file_paths]

# initialize dataframes
df_min = pd.DataFrame()
df_max = pd.DataFrame()
df_mean = pd.DataFrame()
df_count = pd.DataFrame()

dataframes = [df_min, df_max, df_mean, df_count]
indices = [3, 4, 5, 9]


# %%    
for col, df in zip(indices, dataframes):

    # find file that is not empty:
    i = 0
    while open(file_paths[i], "r").readline().isspace(): 
        print(i)
        i += 1 
    # then use this file to extract a multiindex (chromosome + genomic range):
    df["chromosome"] = pd.read_csv(file_paths[i], sep="\t", header=None).iloc[:,0]
    df["genomic_range"] = pd.read_csv(file_paths[i], sep="\t", header=None).iloc[:,1]
    # iterate over the same column in all files and aggregate into a dataframe:
    for sra_id, file in zip(sra_ids, file_paths):
        with open(file, "r") as f:
            # empty files are included as columns of NA-values: 
            if f.readline().isspace():
                df[sra_id] = np.full(df.shape[0], np.nan)
            # non-empty files are included as regular columns:
            else:   
                df[sra_id] = pd.read_csv(file, sep="\t", header=None, na_values=[".", 0]).iloc[:,col]
    # set the multiindex and remove rows that have less than a threshold of non-NA values.
    df.set_index(["chromosome", "genomic_range"], inplace=True)
    df.dropna(thresh=na_threshold, inplace=True)
    if(df.empty):
        continue
    # calculate pearson-rank correlation between peaks and scores for each row
    # (each row correspond to a bin of size 1000bp of ChIP-seq peaks)
    cor = df.apply(lambda row: stats.pearsonr(row[row.notna()], scores[row.notna()]), axis=1)
    df["cor"], df["p"] = zip(*cor)
    del(cor)
    # calculate bin for each row observation based on the correlation coefficient
    bin_labels = ["[-1|-0.8)","[-0.8|-0.6)", "[-0.6|-0.4)", "[-0.4|-0.2)", "[-0.2|0)", 
            "[0|0.2)","[0.2|0.4)", "[0.4|0.6)","[0.6|0.8)", "[0.8|1)"]
    df["bin"] = pd.cut(df["cor"], 10, labels=bin_labels)

# %% 
# write all dataframes to .csv files
for i, df in enumerate(dataframes):
    df.to_csv(sys.argv[i-4])