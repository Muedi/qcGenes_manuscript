"""
Script to load metadata and search for fitting subsets 
"""
# %%
import pandas as pd
import numpy as np
import glob
# import GEOparse
import re
import os

# %%
csv_folder = "../../config/metadata/"
Datasets_table = pd.read_csv(csv_folder + "Datasets_rev0.csv")
Datasets_table.set_index("GEO_Series", inplace=True)
old_mega = pd.read_csv(csv_folder + "Mega_SraRunTable_rev0.csv", index_col=False)
# old_mega["TISSUE"] = old_mega["tissue"]
mega_cols = old_mega.columns

MEGA_SRA_table = old_mega
MEGA_SRA_table = MEGA_SRA_table.loc[~MEGA_SRA_table.GEO_Series.str.contains("sub")]
# %%
data_ids = Datasets_table.index.tolist()
data_ids = [x for x in data_ids if "sub" not in x ]
# get stats for each dataset.
print("Get stats.")
dataset_stats = pd.DataFrame()
for DATAID in data_ids:
    file = os.path.join("../../output", "main", "qualityVsExp", DATAID, DATAID + ".cor.tsv")
    print(file)
    print("Exists: " + str(os.path.isfile(file)))
    if os.path.isfile(file):
        dataset_stats = pd.concat([dataset_stats, pd.read_csv(file, sep="\t")])
dataset_stats.set_index("dataset", inplace=True)

p_low = pd.DataFrame()
for DATAID in data_ids:
    file = os.path.join("../../output", "main", "scores", DATAID + ".scores.txt")
    print(file)
    print("Exists: " + str(os.path.isfile(file)))
    if os.path.isfile(file):
        p_low = pd.concat([p_low, pd.read_csv(file, sep="\t", header=None)])
p_low.columns = ["Run", "P_low", "anno"]
p_low = p_low.drop("anno", axis=1)
p_low.set_index("Run", inplace=True)

deg_stats = {}
for DATAID in data_ids:
    ctrl = Datasets_table.loc[DATAID, "Control"]
    treat = Datasets_table.loc[DATAID, "Treat"]
    file = os.path.join("../../output", "main", "RASflowResults", DATAID, "trans/dea", "DEA", "gene-level",
                        "deg_" + ctrl + "_" + treat + ".tsv") 
    print(file)
    print("Exists: " + str(os.path.isfile(file)))
    if os.path.isfile(file):
        with open(file) as f:
            deg_stats[DATAID] = sum(1 for line in f) -1 # (counts always one more than JFs script, most likeley endline)
    else:
        deg_stats[DATAID] = -1
# convert to df and join with dataset_stats
deg_stats = pd.DataFrame.from_dict(deg_stats, orient='index', columns=['deg'])
dataset_stats = dataset_stats.join(deg_stats)

#%%

## TODO: 
## Histograms or boxplots of all bigger datasets,
# to see where we have big quality differences between the samples
old_mega = old_mega.set_index("Run", drop=False).join(p_low)
# %%
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(25,25))
old_mega["P_low"].hist(by=old_mega["GEO_Series"], ax=ax)
fig.savefig("../../output/main/dist_plow.png")

#%%
# do the same but filter only big datasets with DEGs
# get datasets with more than 50 samples and more than 50 degs
big_datasets = old_mega.loc[
    old_mega.loc[
        old_mega.Selected.eq(1)].index].groupby("GEO_Series").count().index[
            old_mega.loc[
                old_mega.Selected.eq(1)].groupby("GEO_Series").count()["Selected"] > 50].tolist()
sets_gt_50_deg = deg_stats.loc[big_datasets].loc[(deg_stats.loc[big_datasets].deg.gt(50))]
sets_gt_50_deg = sets_gt_50_deg.index.tolist()

filtered_mega = old_mega.loc[old_mega.GEO_Series.isin(sets_gt_50_deg)]
# %%
fig, ax = plt.subplots(figsize=(25,25))
filtered_mega["P_low"].hist(by=filtered_mega["GEO_Series"], ax=ax)
fig.savefig("../../output/main/dist_plow_relevant.png")

# %%
filtered_mega = old_mega.loc[
    old_mega.GEO_Series.isin(["GSE105130", "GSE144269", "GSE74697"])]

# %%
### get subsets that are as heterogeneous as possible
GSE105130 = filtered_mega[filtered_mega.GEO_Series.eq("GSE105130")]
GSE144269 = filtered_mega[filtered_mega.GEO_Series.eq("GSE144269")]
GSE74697 = old_mega[old_mega.GEO_Series.eq("GSE74697")]

#GSE100925_homogeneous_plow_eq = 
# %%
series_list = ["GSE105130", "GSE144269", "GSE74697"]
for GSE in series_list:
    DATA = filtered_mega[filtered_mega.GEO_Series.eq(GSE)]
    treat_str = Datasets_table.loc[GSE, "Treat"]
    ctrl_str = Datasets_table.loc[GSE, "Control"]
    treat = DATA.loc[DATA.group.eq(treat_str)]
    ctrl = DATA.loc[DATA.group.eq(ctrl_str)]

    # no info given for this dataset, so we just orient on PLOW
    treat_lowestP_10 = treat.nsmallest(10, "P_low")
    treat_highestP_10 = treat.nlargest(10, "P_low")
    ctrl_lowestP_10 = ctrl.nsmallest(10, "P_low")
    ctrl_highestP_10 = ctrl.nlargest(10, "P_low")

    for i in range(3):
        treat_lowestP = treat_lowestP_10.sample(5, random_state=i)
        treat_highestP = treat_highestP_10.sample(5, random_state=i)
        ctrl_lowestP = ctrl_lowestP_10.sample(5, random_state=i)
        ctrl_highestP = ctrl_highestP_10.sample(5, random_state=i)
        simi_low = pd.concat([treat_lowestP,ctrl_lowestP])
        simi_low.GEO_Series = simi_low.GEO_Series + "sub_simi_low" + "_{}".format(i)
        simi_high = pd.concat([treat_highestP,ctrl_highestP])
        simi_high.GEO_Series = simi_high.GEO_Series + "sub_simi_high" + "_{}".format(i)
        diff_tu_low = pd.concat([treat_lowestP,ctrl_highestP])
        diff_tu_low.GEO_Series = diff_tu_low.GEO_Series + "sub_diff_tu_low" + "_{}".format(i)
        diff_tu_high = pd.concat([treat_highestP,ctrl_lowestP])
        diff_tu_high.GEO_Series = diff_tu_high.GEO_Series + "sub_diff_tu_high" + "_{}".format(i)
        subsets_output = pd.concat([simi_low,
                                    simi_high,
                                    diff_tu_low,
                                    diff_tu_high])
        subsets_output.to_csv(
                os.path.join("../../output", "main",  "{}subsets_{}.csv".format(GSE, i)),
                index=False
        )
        MEGA_SRA_table = pd.concat([MEGA_SRA_table, subsets_output])

MEGA_SRA_table.to_csv(csv_folder + "Mega_SraRunTable.csv", index=False)