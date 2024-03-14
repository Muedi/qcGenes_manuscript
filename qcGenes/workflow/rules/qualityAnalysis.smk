# =============================================================================
# Snakemake workflow to analyze quality genes
# =============================================================================

# ------------------------------------------------------------------------------
rule vsExp:
    input:
        expand("output/main/qualityVsExp/{dataid}/{dataid}.cor_pos_genes.tsv", dataid=DATAIDS_MAIN),
        "output/main/overlap/figure-6.png"


# ------------------------------------------------------------------------------
rule quality_vs_exp:
    input:
        ensid=config["ENSID_PATH"],
        gs2d=config["GS2D_PATH"],
        scores=rules.aggregate_scores.output,
        rlog=rules.deseq2.output.rlog,
        diff=rules.deseq2.output.diff
    output:
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.cor001.png",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.pca001.png",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.pca_diag001.png",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.scores001.png",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.cor_neg_genes.tsv",
        "output/{analysis}/qualityVsExp/{dataid}/{dataid}.cor_pos_genes.tsv"
    group:
        "qualityAnalysis" 
    log: 
        "output/logs/{analysis}/quality_vs_exp/{dataid}.quality_vs_exp.log"
    params:
        outdir="output/{analysis}/qualityVsExp/{dataid}"
    shell:
        "Rscript workflow/scripts/quality_vs_expression.r {wildcards.dataid} {wildcards.analysis} {input.scores} {input.rlog} {input.diff} {params.outdir} > {log} 2>&1"


# ------------------------------------------------------------------------------
rule qualityAnalysis:
    input:
        get_q_vs_expr_output_files,
        config["MSIGDB_HALLMARK_PATH"],
        config["MSIGDB_POSITIONAL_PATH"],
        config["MSIGDB_CURATED_PATH"],
        config["MSIGDB_CURATED_CP_PATH"],
        config["MSIGDB_REGULATORY_PATH"],
        config["MSIGDB_CELLS_PATH"],
        config["GS2D_MSIGDB_PATH"]
    output:
        fgsea_pos_nobias="output/{analysis}/qualityCorGenes/fgsea/fgsea.pos.cure.noBias.gsea.plot.png",
        fgsea_pos="output/{analysis}/qualityCorGenes/fgsea/fgsea.pos.cure.gsea.plot.png",
        fgsea_pos_highbias="output/{analysis}/qualityCorGenes/fgsea/fgsea.pos.cure.highBias.gsea.plot.png",
        fgsea_neg_nobias="output/{analysis}/qualityCorGenes/fgsea/fgsea.neg.cure.noBias.gsea.plot.png",
        fgsea_neg="output/{analysis}/qualityCorGenes/fgsea/fgsea.neg.cure.gsea.plot.png",
        pos_genes_nobias="output/{analysis}/qualityCorGenes/pos.cor.genes.noBias.png",
        neg_genes_nobias="output/{analysis}/qualityCorGenes/neg.cor.genes.noBias.png",
        degs_vs_samples="output/{analysis}/qualityCorGenes/figure.deg.vs.samples.png",
        dataset_stats="output/{analysis}/qualityCorGenes/datasets.stats.tsv",
        pos_genes_nobias_table="output/{analysis}/qualityCorGenes/pos.cor.genes.noBias.tsv",
        neg_genes_nobias_table="output/{analysis}/qualityCorGenes/neg.cor.genes.noBias.tsv"
    group:
        "qualityAnalysis"
    log:
        "output/logs/{analysis}/qualityCorGenes.log"
    params:
        outdir="output/{analysis}/qualityCorGenes"
    shell:
        "Rscript workflow/scripts/quality_correlated_genes.r {params.outdir} > {log} 2>&1"


# ------------------------------------------------------------------------------
rule overlap:
    input:
        pos_genes=rules.qualityAnalysis.output.pos_genes_nobias_table,
        neg_genes=rules.qualityAnalysis.output.neg_genes_nobias_table,
        dataset_stats=rules.qualityAnalysis.output.dataset_stats,
        gs2d=config["GS2D_ENSEMBL_PATH"],
        files=get_deseq2_output_files
    output:
        "output/{analysis}/overlap/figure.6.png",
        "output/{analysis}/overlap/masterplot.degsVmarkers.small.datasets.png",
        "output/{analysis}/overlap/masterplot.degsVgs2d.small.datasets.png"
    log:
        "output/logs/{analysis}/overlap.log"
    params:
        outdir="output/{analysis}/overlap"
    shell:
        "Rscript workflow/scripts/overlap_analysis.R {input.pos_genes} {input.neg_genes} {input.dataset_stats} {input.gs2d} {params.outdir} {input.files} > {log} 2>&1"


# ------------------------------------------------------------------------------
rule qi_method_eval:
    input:
        get_q_vs_expr_output_files
    output:
        "output/{analysis}/metrics/outliers.png",
        "output/{analysis}/metrics/pvalues.png",
        "output/{analysis}/metrics/values.png"
    log:
        "output/logs/{analysis}/qi_method_eval.log"
    params:
        scores_dir="output/{analysis}/scores",
        outdir="output/{analysis}/metrics"
    shell:
        "Rscript workflow/scripts/qi_method_evaluation.R {params.scores_dir} {params.outdir} > {log} 2>&1"


# # ------------------------------------------------------------------------------
# rule summary:
#     input:
#         "output/main/qualityCorGenes/pos.cor.genes.noBias.png",
#         "output/main/qualityCorGenes/neg.cor.genes.noBias.png",
#         "output/main/qualityCorGenes/fgsea/fgsea.pos.cure.gsea.plot.png",
#         "output/main/qualityCorGenes/fgsea/fgsea.neg.cure.gsea.plot.png"
#     output:
#         "output/summary.html"
#     log:
#         "output/logs/main/summary.log"
#     shell:
#         "Rscript -e \"rmarkdown::render('workflow/scripts/summary.Rmd', output_format='html_document')\" > {log} 2>&1 && "
#         "mv workflow/scripts/summary.html output/"
# #        "Rscript -e \"rmarkdown::render('summary.Rmd', output_format='html_document')\" > {log} 2>&1 && "
# #        "Rscript -e \"rmarkdown::render('summary.Rmd', output_format='pdf_document')\" >> {log} 2>&1"

# ------------------------------------------------------------------------------
rule draft:
    input:
        "output/{analysis}/qualityCorGenes/fgsea/fgsea.pos.cure.noBias.gsea.plot.png",
        "output/{analysis}/qualityCorGenes/fgsea/fgsea.pos.cure.gsea.plot.png",
        "output/{analysis}/qualityCorGenes/fgsea/fgsea.pos.cure.highBias.gsea.plot.png",
        "output/{analysis}/qualityCorGenes/fgsea/fgsea.neg.cure.noBias.gsea.plot.png",
        "output/{analysis}/qualityCorGenes/fgsea/fgsea.neg.cure.gsea.plot.png",
        "output/{analysis}/qualityCorGenes/pos.cor.genes.noBias.png",
        "output/{analysis}/qualityCorGenes/neg.cor.genes.noBias.png",
        "output/{analysis}/qualityCorGenes/figure.deg.vs.samples.png",
        "output/{analysis}/overlap/masterplot.degsVmarkers.small.datasets.png",
        "output/{analysis}/overlap/masterplot.degsVgs2d.small.datasets.png",
        "output/{analysis}/overlap/figure.6.png",
        "output/{analysis}/metrics/outliers.png",
        "output/{analysis}/metrics/pvalues.png",
        "output/{analysis}/metrics/values.png",
        rules.config_backup.output.cfg_out_path
    output:
        "output/{analysis}/draft.html"
    log:
        "output/logs/{analysis}/draft.log"
    shell:
        "Rscript -e \"rmarkdown::render('workflow/scripts/draft.Rmd', output_format='html_document')\" > {log} 2>&1 && "
        "mv workflow/scripts/draft.html output/{wildcards.analysis}/"