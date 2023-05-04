# %%
import os
from numpy.core.numeric import NaN
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib_venn as venn
import sys
from seaborn.categorical import countplot
import yaml
from pathlib import Path
from itertools import chain
from collections import Counter
from gtfparse import read_gtf
import itertools

# %%
######## input variables
# mode = sys.argv[1]

# output
output_path = "output/main/overlap/"
if not os.path.exists(output_path):
    print("Output path did not exist yet, created folders.")
    Path(output_path).mkdir(parents=True, exist_ok=True)
########## Global variables
# config
config_path = "config/config.yaml"
# load yaml file to dict
with open(config_path, 'r') as config_file:
    config = yaml.load(config_file, Loader=yaml.FullLoader)
# gs2d table
gs2d_table = pd.read_csv(config["GS2D_PATH"], sep='\t')
# datasets infomation
datasets = pd.read_csv("config/metadata/Datasets.csv")
datasets.rename(columns={"mesh_terms (pipe-separated list)": "mesh_terms"}, inplace=True)
# filter out deselected sets.
datasets = datasets.loc[datasets.select.eq(1)]
datasets = datasets.loc[~datasets.GEO_Series.str.contains("sub")]
datasets.set_index("GEO_Series", inplace=True)
datasets.drop("select", axis=1, inplace=True) # drop select its not needed anymore
# cut offs for deg number and corr of quali vs sample
dataset_p_group_cor_cutoff = config["MAX_DATASET_P_VS_GROUP_CORRELATION"]
dataset_n_deg_cutoff = config["MIN_DATASET_N_DEG"]
# gene annot
# %%
# get GEO_Series for meshs with more or equal to 2 sets ande compare degs
# the for loop iterates over all keys with more or equal to 2 sets.
gtf_file = config["ANNOT_PATH"]
gtf = read_gtf(gtf_file)

geneid_to_name = gtf[["gene_id", "gene_name"]]
geneid_to_name = geneid_to_name.rename(columns={"gene_id":"geneid"})
geneid_to_name.set_index("geneid", inplace=True)
geneid_to_name = geneid_to_name.drop_duplicates()

# %%
########## statistics about quality and correlations for datasets
# get GEOids
data_ids = datasets.index.tolist()
# get stats for each dataset.
print("Get stats.")
dataset_stats = pd.DataFrame()
for DATAID in data_ids:
    file = os.path.join("output", "main", "qualityVsExp", DATAID, DATAID + ".cor.tsv")
    print(file)
    print("Exists: " + str(os.path.isfile(file)))
    if os.path.isfile(file):
        dataset_stats = pd.concat([dataset_stats, pd.read_csv(file, sep="\t")])
dataset_stats.set_index("dataset", inplace=True)

########## number of degs per dataset
print("Get number of differential expressed genes.")
deg_stats = {}
for DATAID in data_ids:
    ctrl = datasets.loc[DATAID, "Control"]
    treat = datasets.loc[DATAID, "Treat"]
    file = os.path.join("output", "main", "RASflowResults", DATAID, "trans/dea/DEA/gene-level",
                        "dea_" + ctrl + "_" + treat + ".tsv") 
    print(file)
    print("Exists: " + str(os.path.isfile(file)))
    if os.path.isfile(file):
        deg_tmp = pd.read_csv(file, index_col=0, sep='\t')
        deg_tmp = deg_tmp.loc[deg_tmp.log2FoldChange.ge(1) | deg_tmp.log2FoldChange.le(-1)]
        deg_tmp = deg_tmp.loc[deg_tmp.padj.le(0.05)]
        deg_stats[DATAID] = deg_tmp.shape[0] 
    else:
        deg_stats[DATAID] = -1
# convert to df and join with datasets_stats
deg_stats = pd.DataFrame.from_dict(deg_stats, orient='index', columns=['deg'])
dataset_stats = dataset_stats.join(deg_stats)
# possibility to filter out cells not meeting the threshold of deg or quali correlation
dataset_stats["selected"]= ~((dataset_stats.deg<dataset_n_deg_cutoff).to_numpy() | 
                             (dataset_stats.p_group_cor>dataset_p_group_cor_cutoff).to_numpy()
                            )
datasets = datasets.join(dataset_stats.selected)
# remove NA in mesh term
datasets = datasets.dropna(subset=["mesh_terms"])
########## summarize mesh terms and iterate over datsets with equal terms.

# get all meshterms and...
mesh_terms = datasets.loc[datasets.selected.eq(True), "mesh_terms"].dropna()
mesh_terms = mesh_terms.str.split("|") # list of lists
mesh_terms =  list(chain.from_iterable(mesh_terms))   # flatten
# ...count mesh occurences
counts_mesh = Counter(mesh_terms)

 # %%
# build a df with all datasets and all gene quali corralations
#  degs by GEO_series
GEO_Series = datasets.index.tolist()
degs = {} # all info
lfc_genes_dict = {} # only lfc
pval_genes_dict = {} # only pval adj

for DATAID in GEO_Series:
    ctrl = datasets.loc[DATAID, "Control"]
    treat = datasets.loc[DATAID, "Treat"]
    deg_file = os.path.join("output", "main", "RASflowResults", DATAID, "trans/dea/DEA/gene-level",
                            "deg_" + ctrl + "_" + treat + ".tsv")
    dea_file = os.path.join("output", "main", "RASflowResults", DATAID, "trans/dea/DEA/gene-level",
                            "dea_" + ctrl + "_" + treat + ".tsv")
    if os.path.isfile(deg_file):
        degs[DATAID] = pd.read_csv(dea_file, index_col=0, sep='\t')
        degs[DATAID] = degs[DATAID].loc[degs[DATAID].log2FoldChange.ge(1) | degs[DATAID].log2FoldChange.le(-1)]
        degs[DATAID] = degs[DATAID].loc[degs[DATAID].padj.le(0.05)]
        tmp_df = pd.read_csv(dea_file, sep='\t')
        lfc_genes_dict[DATAID] = tmp_df["log2FoldChange"]
        pval_genes_dict[DATAID] = tmp_df["padj"]


# genes and their correlation
"""########## get counted genes correlating with quality """
cor_file = os.path.join("output", "main", "qualityCorGenes", "pos.cor.genes.noBias.tsv")
neg_cor_file = os.path.join("output", "main", "qualityCorGenes", "neg.cor.genes.noBias.tsv")
if os.path.isfile(cor_file):
    # we need to use the complete corr table, since the others used cutoffs for corr coef
    cor_genes = pd.read_csv(cor_file, sep='\t')
    cor_genes = cor_genes[["geneid", "n"]].set_index("geneid")
    neg_cor_genes = pd.read_csv(neg_cor_file, sep='\t').set_index("geneid")

input_df = pd.concat(degs)
input_df = input_df.reset_index().rename(columns={"level_0":"Dataset", "level_1":"geneid"})
print(len(set(input_df.geneid.tolist())))


# add mesh terms
input_df['mesh'] = [datasets.loc[DATAID, 'mesh_terms'].split("|") for DATAID in input_df.Dataset]
input_df['mesh_string'] = [datasets.loc[DATAID, 'mesh_terms'] for DATAID in input_df.Dataset]

# literature genes vs non related genes with GS2D
# input_df = input_df.set_index("geneid")
#add symbol
input_df = input_df.join(geneid_to_name, on="geneid")    
input_df = input_df.rename(columns={"gene_name":"symbol"})
# add gs2d information in column "related_to_mesh":
# add related_to column if gene is related to the mesh, that the set belongs to
input_df['related_to'] = np.empty((len(input_df), 0)).tolist()
# use counts_mesh again to get all mesh terms
for mesh in counts_mesh.keys():
    sub_gs2d = gs2d_table.loc[gs2d_table.name.eq(mesh), "symbol" ]
    bool_array =  list(np.where(input_df.symbol.isin(sub_gs2d), True, False)) # checks where the related genes to given mesh are in df
    tmp_series = input_df.loc[bool_array, "related_to"].apply(lambda x: x+[mesh])
    input_df.loc[bool_array, "related_to"] = tmp_series

# add column, that says the gene is related to the mesh term given to the sample in the col mesh 
# input_df["is_related"] = [False] * len(input_df.index)

A = input_df["mesh"].tolist()
B = input_df["related_to"].tolist()
is_related = [False] * len(A)

for i in range(len(A)): 
    if not A or not B: 
        pass
    if any(b in A[i] for b in B[i]):
        is_related[i] = True

input_df["is_related"] = is_related
# %%
# investigate HCC
hcc_sets = input_df.loc[input_df.mesh_string.str.contains("Hepatocellular"),]
deg_hcc = hcc_sets.loc[~hcc_sets.Dataset.eq("GSE82177"),] # remove set with only 2 deg




# upsetplot
from upsetplot import UpSet, from_contents

grp = deg_hcc.groupby("Dataset")
sets = grp.apply(lambda x : set(x["geneid"]))
hcc_input_upset = from_contents(sets)



fig, ax = plt.subplots()
plt.axis('off')
plt.grid(b=None)
upset = UpSet(hcc_input_upset, subset_size='count', intersection_plot_elements=3, orientation="vertical", show_counts='%d')
# upset.add_catplot(value='median_value', kind='strip', color='blue')
# upset.add_catplot(value='AGE', kind='strip', color='black')
upset.plot(fig)

ax.set_title("Upsetplot for all DEGS in HCC Datsaets")
fig.set_size_inches(10, 22)
fig.savefig(output_path + "HCC-sets-deg-investigation.png", bbox_inches='tight')

# %%
from scipy.stats import fisher_exact
# do this for all possible combinations:
l = [False, True]
accessions = list(hcc_input_upset.index.names)
combinations = list(itertools.product(l, repeat=len(accessions)))

interesting_genes = gs2d_table.loc[(gs2d_table.name.eq("Carcinoma, Hepatocellular") ) #| gs2d_table.name.eq("Liver Neoplasms") )
                    & gs2d_table.fdr.le(0.05) 
                    & gs2d_table.count_name_in_gene_set.ge(5)
                    ] 

select_sets = 13
perc = np.floor(select_sets * 0.2)
print(perc)
quality_markers = cor_genes.loc[cor_genes.n.ge(perc)]
# quality_markers = quality_markers[~quality_markers.index.duplicated(keep='first')]
# print(quality_markers)

dict_all_subsets = {}
for combi in combinations:
    try:
        subset_overlap = hcc_input_upset.loc[combi, "id"].to_list()
        subset_overlap_symbol = [geneid_to_name.loc[geneid, "gene_name"] for geneid in subset_overlap if geneid in geneid_to_name.index ]
        gs2d_subset_overlap = interesting_genes.loc[interesting_genes.symbol.isin(subset_overlap_symbol)]
        quali_subset_overlap = quality_markers.reindex(subset_overlap).dropna()   
        dict_all_subsets[" + ".join(list(itertools.compress(accessions, combi)))] = {"num_degs":len(subset_overlap_symbol),
                                                                                 "num_gs2d":gs2d_subset_overlap.shape[0],
                                                                                 "ratio_gs2d":gs2d_subset_overlap.shape[0]/len(subset_overlap_symbol),
                                                                                 "num_quali":quali_subset_overlap.shape[0],
                                                                                 "ratio_quali":quali_subset_overlap.shape[0]/len(subset_overlap_symbol)} 
    except KeyError:
        pass
hcc_all_subsets = pd.DataFrame.from_dict(dict_all_subsets).T
hcc_all_subsets.sort_values("num_degs").to_csv(output_path + "hcc-all-subset-combinations.csv")


# number of genes: 21306
# ref: https://www.nature.com/articles/d41586-018-05462-w

### meta analysis of combination counts and ratios
subset_exclusive_genes = hcc_all_subsets.loc[list(set(deg_hcc.Dataset))]
bad_exclusive_genes = hcc_all_subsets.loc[["GSE105130", "GSE77314", "GSE105130 + GSE77314"]]
# comis that contain overlaps of more than 2 sets
combis_more_than_two = [combi for combi in combinations if sum(combi) > 2]
combis_more_than_two = [" + ".join(list(itertools.compress(accessions, combi))) for combi in combis_more_than_two]
combis_more_than_two = [combi for combi in combis_more_than_two if combi in hcc_all_subsets.index.tolist()]
more_than_2_subsets_genes = hcc_all_subsets.loc[combis_more_than_two]
# comis that contain overlaps of more than 3 sets
combis_more_than_three = [combi for combi in combinations if sum(combi) > 3]
combis_more_than_three = [" + ".join(list(itertools.compress(accessions, combi))) for combi in combis_more_than_three]
combis_more_than_three = [combi for combi in combis_more_than_three if combi in hcc_all_subsets.index.tolist()]
more_than_3_subsets_genes = hcc_all_subsets.loc[combis_more_than_three]
# comis that contain overlaps of more than 4 sets
combis_more_than_four = [combi for combi in combinations if sum(combi) > 4]
combis_more_than_four = [" + ".join(list(itertools.compress(accessions, combi))) for combi in combis_more_than_four]
combis_more_than_four = [combi for combi in combis_more_than_four if combi in hcc_all_subsets.index.tolist()]
more_than_4_subsets_genes = hcc_all_subsets.loc[combis_more_than_four]
# comis that contain overlaps of more than 5 sets
combis_more_than_five = [combi for combi in combinations if sum(combi) > 5]
combis_more_than_five = [" + ".join(list(itertools.compress(accessions, combi))) for combi in combis_more_than_five]
combis_more_than_five = [combi for combi in combis_more_than_five if combi in hcc_all_subsets.index.tolist()]
more_than_5_subsets_genes = hcc_all_subsets.loc[combis_more_than_five]


# divide into with and without bad influence
combis_more_than_two_bad = [combi for combi in combis_more_than_two if "GSE105130" in combi or "GSE77314" in combi]
more_than_2_subsets_genes_bad = hcc_all_subsets.loc[combis_more_than_two_bad]
combis_more_than_two_not_bad = [combi for combi in combis_more_than_two if combi not in combis_more_than_two_bad]
more_than_2_subsets_genes_not_bad = hcc_all_subsets.loc[combis_more_than_two_not_bad]

gene_sets = [subset_exclusive_genes,
            bad_exclusive_genes,
            more_than_2_subsets_genes,
            more_than_2_subsets_genes_bad,
            more_than_2_subsets_genes_not_bad,
            more_than_3_subsets_genes,
            more_than_4_subsets_genes,
            more_than_5_subsets_genes]

num_degs = []
num_gs2d = []
ratios = []
pvals_combis_subs = []
num_quali = []
ratios_quali = []
pvals_combis_subs_quali = []
for item in gene_sets:
    num_degs.append(item.sum().num_degs)
    num_gs2d.append(item.sum().num_gs2d)
    num_quali.append(item.sum().num_quali)

    ratio_item = item.sum().num_gs2d/item.sum().num_degs
    ratios.append(ratio_item)
    ratio_item = item.sum().num_quali/item.sum().num_degs
    ratios_quali.append(ratio_item)
    # fishers exact
    table = np.array([  [item.sum().num_degs, 21306 - item.sum().num_degs],
                        [item.sum().num_gs2d, interesting_genes.shape[0] - item.sum().num_gs2d]])
    odds_item, p_item = fisher_exact(table)
    pvals_combis_subs.append(p_item)

    table = np.array([  [item.sum().num_degs, 21306 - item.sum().num_degs],
                        [item.sum().num_quali, quality_markers.shape[0] - item.sum().num_quali]])
    odds_item, p_item = fisher_exact(table)
    pvals_combis_subs_quali.append(p_item)

combi_ratios = pd.DataFrame({
    "num_degs":num_degs,
    "num_gs2d":num_gs2d,
    "ratio_gs2d":ratios,
    "pval_gs2d":pvals_combis_subs,
    "num_quali":num_quali,
    "ratio_quali":ratios_quali,
    "pval_quali":pvals_combis_subs_quali
    }, 
    index=["set_exclusive_genes",
    "bad_set_genes",
    "genes_of_overlapped(>2)",
    "genes_of_overlapped(>2) bad influence",
    "genes_of_overlapped(>2) no bad infl",
    "genes_of_overlapped(>3)",
    "genes_of_overlapped(>4)",
    "genes_of_overlapped(>5)"] 
)

# pvals for each subset
pvals_excl_subs = []
pvals_excl_subs_quali = []
for DATAID in subset_exclusive_genes.index:
    table = np.array([  [subset_exclusive_genes.loc[DATAID].num_degs, 21306 - subset_exclusive_genes.loc[DATAID].num_degs],
                        [subset_exclusive_genes.loc[DATAID].num_gs2d, interesting_genes.shape[0] - subset_exclusive_genes.loc[DATAID].num_gs2d]])
    pvals_excl_subs.append(fisher_exact(table)[1])
    
    table = np.array([  [subset_exclusive_genes.loc[DATAID].num_degs, 21306 - subset_exclusive_genes.loc[DATAID].num_degs],
                        [subset_exclusive_genes.loc[DATAID].num_quali, quality_markers.shape[0] - subset_exclusive_genes.loc[DATAID].num_quali]])

    pvals_excl_subs_quali.append(fisher_exact(table)[1])
subset_exclusive_genes["pval_gs2d"] = pvals_excl_subs
subset_exclusive_genes["pval_quali"] = pvals_excl_subs_quali

output_ratios = pd.concat([subset_exclusive_genes, combi_ratios])
# add bonferroni correction, multiply by number of meshterms that are able to be detected in the gs2d table and number of quality markers
output_ratios["FDR_Bonferroni_gs2d"] = output_ratios["pval_gs2d"] * len(set(gs2d_table.name))
output_ratios["fold_change_gs2d"] = output_ratios["ratio_gs2d"] / (interesting_genes.shape[0]/21306)
output_ratios["FDR_Bonferroni_quali"] = output_ratios["pval_quali"] * quality_markers.shape[0]
output_ratios["fold_change_quali"] = output_ratios["ratio_quali"] / (quality_markers.shape[0]/21306)

output_ratios = output_ratios[['num_degs',  
                            'num_quali', 'ratio_quali', 'pval_quali', 'FDR_Bonferroni_quali', 'fold_change_quali',
                            'num_gs2d', 'ratio_gs2d', 'pval_gs2d', 'FDR_Bonferroni_gs2d', 'fold_change_gs2d']
                            ]

output_ratios.to_csv(output_path + "hcc-ratios-and-pvals-hcconly-default.csv")


#%%
# first plots count of all set's degs and add an venn diagram
fig, ax = plt.subplots(1, 1, figsize=(7, 15))
ax.set_title("HCC: Ratio of relevant genes according to gs2d and degs in subsets", fontsize=14)
sns.barplot(x=hcc_all_subsets.sort_values("num_degs")["ratio_gs2d"], 
            y=hcc_all_subsets.sort_values("num_degs").index , 
            color="darkorange",
            orient="h",
            ax=ax)
            
# annotate subsets which contain high p_group_cor sets (GSE105130, GSE77314)
list_of_labels = [ax.get_yticklabels()[i].get_text() for i in range(len(ax.get_yticklabels()))]
index_for_high_p_group = [i for i in range(len(list_of_labels)) if "GSE105130" in list_of_labels[i] and "GSE77314" in list_of_labels[i]]

for y in index_for_high_p_group:
    ax.annotate('p_group > 0.3', 
                xy=(hcc_all_subsets.sort_values("num_degs")["ratio_gs2d"][y] + 0.01, y), 
                xytext=(hcc_all_subsets.sort_values("num_degs")["ratio_gs2d"][y] + 0.05, y),
                arrowprops=dict(facecolor='red', shrink=0.05)
                )
ax.text(0.25, 5, "GSE105130, GSE77314 have\na p_group_cor of > 0.3", fontsize=12)
ax.set_xlim(0, 0.55)

fig.savefig(output_path + "HCC: Ratios of relevant genes vs deg.png", dpi=300, bbox_inches='tight')

#%%
# 
accessions = list(hcc_input_upset.index.names)
combinations = list(itertools.product(l, repeat=len(accessions)))
# query gs2d
base_url = "http://cbdm-01.zdv.uni-mainz.de/~jfontain/cgi-bin/genes2diseases.pl?"
base_url_bg = "http://cbdm-01.zdv.uni-mainz.de/~jfontain/cgi-bin/genes2diseases-test.pl?"
# mandatory query items:
# analysis type: in our case always geneset
analy_type = "analysis_type=geneset"
# list of genes
gene_list = "APP|BACE1|PSEN1|MAPT|APOE" # exemplary, will be overwriten with every iteration
gene_list_query = "items={}".format(gene_list) 
# type of output (always text)
output_type = "output_type=text"
max_fdr = "fdr_cutoff=0.001"
# background
bg = "bg_type=selection&bg_items={}".format("Carcinoma,%20Hepatocellular")
# unique_genes_combi = [  (True, True, True, True, True, True),
#                         (False, True, False, False, False, False),
#                         (False, False, True, False, False, False),
#                         (False, False, False, True, False, False),
#                         (False, False, False, False, True, False),
#                         (False, False, False, False, False, True)  ]

# query unique genes

dict_all_subsets = {}
dict_min_fdr = {}
for combi in combinations:
    try:
        subset_overlap = hcc_input_upset.loc[combi, "id"].to_list()
        subset_overlap_symbol = [geneid_to_name.loc[geneid, "gene_name"] for geneid in subset_overlap if geneid in geneid_to_name.index ]
        gene_list = "|".join(subset_overlap_symbol)
        gene_list_query1 = "items={}".format(gene_list[:round(len(gene_list)/2)]) 
        gene_list_query2 = "items={}".format(gene_list[round(len(gene_list)/2):]) 
        df1 = pd.read_table(base_url + "&".join([analy_type, gene_list_query1, output_type]), sep="\t")
        df2 = pd.read_table(base_url + "&".join([analy_type, gene_list_query2, output_type]), sep="\t")
        df = pd.concat([df1, df2])
        dict_all_subsets[" + ".join(list(itertools.compress(accessions, combi)))] = df[["Disease", "Genes count", "Fold change", "P-value", "FDR", "Gene symbols"]]
        dict_min_fdr[" + ".join(list(itertools.compress(accessions, combi)))] = df[["Disease", "Genes count", "Fold change", "P-value", "FDR", "Gene symbols"]].sort_values("FDR").iloc[:2]
    except KeyError:
        pass

# min_fdrs_subsets = pd.DataFrame(dict_all_subsets, index=["minFDR"]).T
fdrs_subsets = pd.concat(dict_all_subsets, ignore_index=False)
fdrs_subsets.index = fdrs_subsets.index.droplevel(1)
fdrs_subsets.to_csv(output_path + "hcc-all-fdr-gs2d.csv")
#relevant mesh
rel_mesh = ["Liver Neoplasms", "Carcinoma, Hepatocellular", "Neoplasms",
            "Liver Disease", "Liver Failure", "Liver Failure, Acute"]
fdrs_subsets.loc[fdrs_subsets.Disease.isin(rel_mesh),].to_csv(output_path + "hcc-rel-mesh-fdr-gs2d.csv")
# sort the top x
min_fdrs_subsets = pd.concat(dict_min_fdr, ignore_index=False)
min_fdrs_subsets.index = min_fdrs_subsets.index.droplevel(1)
min_fdrs_subsets.to_csv(output_path + "hcc-top-2-fdr-gs2d.csv")