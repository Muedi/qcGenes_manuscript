# %%
import os
from collections import Counter
import pandas as pd
# positively correlating genes
qcgenes_gene_positive = pd.read_csv("qcGenes/output/main_P1_default_DESEQ/qualityCorGenes/pos.cor.genes.tsv", sep='\t')
qcgenes_gene_positive_noBias = pd.read_csv("qcGenes/output/main_P1_default_DESEQ/qualityCorGenes/pos.cor.genes.noBias.tsv", sep='\t')
chip_gene_positive = pd.read_csv("qcgenesChIP/snakemake_workflow/output/positively_correlated_genes.csv")
chip_gene_positive = chip_gene_positive.sort_values("dataset_corr_count", ascending=0)
# negatively correlating genes 
qcgenes_gene_negative = pd.read_csv("qcGenes/output/main_P1_default_DESEQ/qualityCorGenes/neg.cor.genes.tsv", sep='\t')
qcgenes_gene_negative_noBias = pd.read_csv("qcGenes/output/main_P1_default_DESEQ/qualityCorGenes/neg.cor.genes.noBias.tsv", sep='\t')
chip_gene_negative = pd.read_csv("qcgenesChIP/snakemake_workflow/output/negatively_correlated_genes.csv")
chip_gene_negative = chip_gene_negative.sort_values("dataset_corr_count", ascending=0)

# %%
# compare gene lists pos
qcgenes = set(
    qcgenes_gene_positive_noBias.loc[qcgenes_gene_positive_noBias.n.gt(2), "genes"])
chip = set(
    chip_gene_positive.loc[chip_gene_positive.dataset_corr_count.gt(1), "SYMBOL"]
)
overlap_pos = set(qcgenes).intersection(set(chip))
print("length ol: ", len(overlap_pos))
print(overlap_pos)
print(qcgenes_gene_positive_noBias.loc[qcgenes_gene_positive_noBias.genes.isin(overlap_pos)].head(10))
print(chip_gene_positive.loc[chip_gene_positive.SYMBOL.isin(overlap_pos)].head(10))
qcgenes_gene_positive_noBias.loc[qcgenes_gene_positive_noBias.genes.isin(overlap_pos)].to_csv("revision_1_overlapping_low_quality_markers.tsv", sep="\t")

# %%
qcgenes = set(
    qcgenes_gene_negative.loc[qcgenes_gene_negative.n.gt(4), "genes"])
chip = set(
    chip_gene_negative.loc[chip_gene_negative.dataset_corr_count.gt(1), "SYMBOL"]
)
overlap_neg = set(qcgenes).intersection(set(chip))
print("length ol: ", len(overlap_neg))
print(overlap_neg)
print(qcgenes_gene_negative.loc[qcgenes_gene_negative.genes.isin(overlap_neg)].head(10))
print(chip_gene_negative.loc[chip_gene_negative.SYMBOL.isin(overlap_neg)].head(10))
qcgenes_gene_negative.loc[qcgenes_gene_negative.genes.isin(overlap_neg)].to_csv("revision_1_overlapping_high_quality_markers.tsv", sep="\t")
#%%
# psoitively correlating pathways

qcgenes_path_positive_noBias = pd.read_csv("qcGenes/output/main/qualityCorGenes/fgsea/fgsea.pos.cure.noBias.gsea.sets.tsv", sep='\t')
chip_path_positive = pd.read_csv("qcgenesChIP/snakemake_workflow/output/positively_correlated_pathways.csv")
chip_path_positive = chip_path_positive.sort_values("padj").head(100)
qcgenes_path_negative_noBias = pd.read_csv("qcGenes/output/main/qualityCorGenes/fgsea/fgsea.neg.cure.noBias.gsea.sets.tsv", sep='\t')
chip_path_negative = pd.read_csv("qcgenesChIP/snakemake_workflow/output/negatively_correlated_pathways.csv")
chip_path_negative = chip_path_negative.sort_values("padj").head(100)
#%%
qcgenes = set(qcgenes_path_positive_noBias["pathway"])
chip = set(chip_path_positive["pathway"])
overlap_pos = set(qcgenes).intersection(set(chip))
print("length ol: ", len(overlap_pos))
print(overlap_pos)
print(qcgenes_path_positive_noBias.loc[qcgenes_path_positive_noBias.pathway.isin(overlap_pos)].head(10))
print(chip_path_positive.loc[chip_path_positive.pathway.isin(overlap_pos)].head(10))
qcgenes_path_positive_noBias.loc[qcgenes_path_positive_noBias.pathway.isin(overlap_pos)].to_csv("revision_1_overlapping_low_quality_marker_pathways.tsv", sep="\t")
#%%
qcgenes = set(qcgenes_path_negative_noBias["pathway"])
chip = set(chip_path_negative["pathway"])
overlap_neg = set(qcgenes).intersection(set(chip))
print("length ol: ", len(overlap_neg))
print(overlap_neg)
print(qcgenes_path_negative_noBias.loc[qcgenes_path_negative_noBias.pathway.isin(overlap_neg)].head(10))
print(chip_path_negative.loc[chip_path_negative.pathway.isin(overlap_neg)].head(10))
qcgenes_path_negative_noBias.loc[qcgenes_path_negative_noBias.pathway.isin(overlap_neg)].to_csv("revision_1_overlapping_high_quality_marker_pathways.tsv", sep="\t")

