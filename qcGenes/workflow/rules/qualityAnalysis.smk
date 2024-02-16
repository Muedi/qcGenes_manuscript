# =============================================================================
# Snakemake workflow to analyze quality genes
# =============================================================================


# ------------------------------------------------------------------------------
rule batchedQualityAnalysis:
    input:
        expand("output/batched/qualityVsExp/{dataid}/{dataid}.cor.pdf", dataid=DATAIDS_BATCH)


# ------------------------------------------------------------------------------
rule qualityAnalysis:
    input:
        expand("output/main/qualityVsExp/{dataid}/{dataid}.cor.pdf", dataid=DATAIDS_MAIN),
	expand("output/batched/qualityVsExp/{dataid}/{dataid}.cor.pdf", dataid=DATAIDS_BATCH),
#        get_all_plot_files,
        config["MSIGDB_HALLMARK_PATH"],
        config["MSIGDB_POSITIONAL_PATH"],
        config["MSIGDB_CURATED_PATH"],
        config["MSIGDB_CURATED_CP_PATH"],
        config["MSIGDB_REGULATORY_PATH"],
        config["MSIGDB_CELLS_PATH"],
        config["GS2D_MSIGDB_PATH"]
    output:
        "output/main/qualityCorGenes/fgsea/fgsea.pos.cure.noBias.gsea.plot.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.pos.cure.gsea.plot.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.pos.cure.highBias.gsea.plot.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.neg.cure.noBias.gsea.plot.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.neg.cure.gsea.plot.pdf",
        "output/main/qualityCorGenes/pos.cor.genes.noBias.pdf",
        "output/main/qualityCorGenes/neg.cor.genes.noBias.pdf"
    group:
        "qualityAnalysis"
    log:
        "output/logs/main/qualityCorGenes.log"
    conda:
        "../envs/qualityCorGenes.yaml"
    shell:
        "Rscript workflow/scripts/quality_correlated_genes.r > {log} 2>&1"


# ------------------------------------------------------------------------------
rule quality_vs_exp:
    input:
        config["ENSID_PATH"],
        config["GS2D_PATH"],
        rules.aggregate_scores.output,
        rules.create_metadata_file.output,
        rules.run_pipeline.output,
        rules.post_pipeline.output
    output:
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.cor.pdf",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.pca.pdf",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.pca_diag.pdf",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.scores.pdf",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.cor_neg_genes.tsv",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.cor_pos_genes.tsv"
    group:
        "qualityAnalysis" 
    log:
        "output/logs/{analysis}/quality_vs_exp/{dataid}.quality_vs_exp.log"
    conda:
        "../envs/qualityAnalysis.yaml"
    shell:
        "Rscript workflow/scripts/quality_vs_expression.r {wildcards.dataid} {wildcards.analysis} > {log} 2>&1"

# ------------------------------------------------------------------------------
rule OL_analysis:
    input:
        rules.qualityAnalysis.output
        # "output/main/qualityCorGenes/pos.cor.genes.pdf"
    output:
        "output/main/overlap/hcc-top-2-fdr-gs2d.csv",
        "output/main/overlap/figure-6.pdf"
    group:
        "qualityAnalysis" 
    log:
        "output/logs/main/HCC_OL.log"
    conda:
        "../envs/OL_py.yaml"
    shell:
        # python3 workflow/overlap/OL_HCC_analysis.py main {output} > {log} 2>&1
        """
        Rscript workflow/overlap/process.gs2d.R  > {log} 2>&1
        Rscript workflow/overlap/overlap.analysis.R  > {log} 2>&1
        """
# ------------------------------------------------------------------------------
rule summary:
    input:
        "output/main/qualityCorGenes/pos.cor.genes.noBias.pdf",
        "output/main/qualityCorGenes/neg.cor.genes.noBias.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.pos.cure.gsea.plot.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.neg.cure.gsea.plot.pdf"
#	,
#	"output/batched/qualityVsExp/GSE163214/GSE163214.batch.pdf",
#	"output/batched/qualityVsExp/GSE173078/GSE173078.batch.pdf",
#	"output/batched/qualityVsExp/GSE81871/GSE81871.batch.pdf",
#	"output/batched/qualityVsExp/GSE82177/GSE82177.batch.pdf",
#	"output/batched/qualityVsExp/GSE120099/GSE120099.batch.pdf"
    output:
        "output/summary.html"
    log:
        "output/logs/main/summary.log"
    conda:
        "../envs/qualityAnalysis.yaml"
    shell:
        "Rscript -e \"rmarkdown::render('workflow/scripts/summary.Rmd', output_format='html_document')\" > {log} 2>&1 && "
        "mv workflow/scripts/summary.html output/"
#        "Rscript -e \"rmarkdown::render('summary.Rmd', output_format='html_document')\" > {log} 2>&1 && "
#        "Rscript -e \"rmarkdown::render('summary.Rmd', output_format='pdf_document')\" >> {log} 2>&1"

# ------------------------------------------------------------------------------
rule draft:
    input:
        "output/main/qualityCorGenes/fgsea/fgsea.pos.cure.noBias.gsea.plot.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.pos.cure.gsea.plot.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.pos.cure.highBias.gsea.plot.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.neg.cure.noBias.gsea.plot.pdf",
        "output/main/qualityCorGenes/fgsea/fgsea.neg.cure.gsea.plot.pdf",
        "output/main/qualityCorGenes/pos.cor.genes.noBias.pdf",
        "output/main/qualityCorGenes/neg.cor.genes.noBias.pdf",
        "output/main/overlap/hcc-top-2-fdr-gs2d.csv",
        "output/main/overlap/figure-6.pdf"
#	,
#	"output/batched/qualityVsExp/GSE163214/GSE163214.batch.pdf",
#	"output/batched/qualityVsExp/GSE173078/GSE173078.batch.pdf",
#	"output/batched/qualityVsExp/GSE81871/GSE81871.batch.pdf",
#	"output/batched/qualityVsExp/GSE82177/GSE82177.batch.pdf",
#	"output/batched/qualityVsExp/GSE120099/GSE120099.batch.pdf"
    output:
        "output/draft.html"
    log:
        "output/logs/main/draft.log"
    conda:
        "../envs/qualityAnalysis.yaml"
    shell:
        "Rscript -e \"rmarkdown::render('workflow/scripts/draft.Rmd', output_format='html_document')\" > {log} 2>&1 && "
        "mv workflow/scripts/draft.html output/"
#        "Rscript -e \"rmarkdown::render('summary.Rmd', output_format='html_document')\" > {log} 2>&1 && "
#        "Rscript -e \"rmarkdown::render('summary.Rmd', output_format='pdf_document')\" >> {log} 2>&1"


# ==============================================================================
# BATCH ANALYSIS
# ==============================================================================

# ------------------------------------------------------------------------------
#rule batchedQualityAnalysis:
#    input:
#        expand("output/batched/qualityVsExp/{dataid}/{dataid}.cor.pdf", dataid=DATAIDS_BATCH)
#        get_all_batched_plot_files

# ------------------------------------------------------------------------------
#rule quality_vs_exp:
#    input:
#        config["ENSID_PATH"],
#        config["GS2D_PATH"],
#        rules.aggregate_scores.output,
#        rules.create_metadata_file.output,
#        rules.run_pipeline.output,
#        rules.post_pipeline.output
#    output:
#        report("output/main/qualityVsExp/{dataid}/{dataid}.cor.pdf", caption="../../report/dataid.cor.rst", category="Quality and Genes"),
#        report("output/main/qualityVsExp/{dataid}/{dataid}.pca.pdf", caption="../../report/dataid.pca.rst", category="Outliers Removal PCA"),
#        "output/main/qualityVsExp/{dataid}/{dataid}.pca_diag.pdf",
#        report("output/main/qualityVsExp/{dataid}/{dataid}.scores.pdf", caption="../../report/dataid.scores.rst", category="Quality vs Size and Diagnosis"),
#        "output/main/qualityVsExp/{dataid}/{dataid}.cor_neg_genes.tsv",
#        "output/main/qualityVsExp/{dataid}/{dataid}.cor_pos_genes.tsv"
#    group:
#        "qualityAnalysis" 
#    log:
#        "output/logs/main/quality_vs_exp/{dataid}.quality_vs_exp.log"
#    conda:
#        "../envs/qualityAnalysis.yaml"
#    shell:
#        "Rscript workflow/scripts/quality_vs_expression.r {wildcards.dataid} main > {log} 2>&1"



# ------------------------------------------------------------------------------
#rule batched_quality_vs_exp:
#    input:
#        config["GS2D_PATH"],
#        rules.aggregate_batched_scores.output,
#        rules.create_batched_metadata_file.output,
#        rules.run_batched_pipeline.output,
#        rules.post_batched_pipeline.output
#    output:
#        report("output/batched/qualityVsExp/{dataid}/{dataid}.cor.pdf", caption="../../report/dataid.cor.rst", category="Batched Quality and Genes"),
#        report("output/batched/qualityVsExp/{dataid}/{dataid}.pca.pdf", caption="../../report/dataid.pca.rst", category="Batched Outliers Removal PCA"),
#        "output/batched/qualityVsExp/{dataid}/{dataid}.pca_diag.pdf",
#        report("output/batched/qualityVsExp/{dataid}/{dataid}.scores.pdf", caption="../../report/dataid.scores.rst", category="Batched Quality vs Size and Diagnosis"),
#        "output/batched/qualityVsExp/{dataid}/{dataid}.cor_neg_genes.tsv",
#        "output/batched/qualityVsExp/{dataid}/{dataid}.cor_pos_genes.tsv"
#    group:
#        "qualityAnalysis" 
#    log:
#        "output/logs/batched/quality_vs_exp/{dataid}.quality_vs_exp.log"
#    conda:
#        "../envs/qualityAnalysis.yaml"
#    shell:
#        "Rscript workflow/scripts/quality_vs_expression.r {wildcards.dataid} batched > {log} 2>&1"
