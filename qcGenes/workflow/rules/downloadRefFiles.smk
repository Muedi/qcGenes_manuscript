# =============================================================================
# Snakemake workflow to download all reference files (data dependencies)
# =============================================================================
import os.path

# ------------------------------------------------------------------------------
rule downloadRefFiles:
    input:
        config["GS2D_PATH"],
        config["GS2D_MSIGDB_PATH"],
        config["GS2D_ENSEMBL_PATH"],
        config["ANNOT_PATH"],
        config["ENSID_PATH"],
        config["ENSIDFULLVER_PATH"],
        config["TRANS_PATH"],
        config["GENOM_PATH"],
        config["BOWTIE_IDX_FILE1"],
        config["MSIGDB_HALLMARK_PATH"],
        config["MSIGDB_POSITIONAL_PATH"],
        config["MSIGDB_CURATED_PATH"],
        config["MSIGDB_CURATED_CP_PATH"],
        config["MSIGDB_REGULATORY_PATH"],
        config["MSIGDB_CELLS_PATH"],
        config["PICARD_REFFLAT_FILE"]

# GS2D
# ------------------------------------------------------------------------------
rule get_gs2d_ref:
    output:
        config["GS2D_PATH"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/get_gs2d_ref.log"
    params:
        url=config["GS2D_URL"]
    shell:
        "wget -nv -m -nH -O {output} {params.url} > {log} 2>&1"

# GS2D_2_MSIGDB
# ------------------------------------------------------------------------------
rule gs2d_2_msigdb:
    input:
        config["GS2D_PATH"]
    output:
        config["GS2D_MSIGDB_PATH"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/get_gs2d_msigdb_ref.log"
    shell:
        "Rscript workflow/scripts/gs2d_2_msigdb.r {input} {output} > {log} 2>&1"

# GS2D_2_ENSEMBL
# ------------------------------------------------------------------------------
rule gs2d_2_ensembl:
    input:
        config["GS2D_PATH"]
    output:
        config["GS2D_ENSEMBL_PATH"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/get_gs2d_ensembl_ref.log"
    shell:
        "Rscript workflow/scripts/gs2d_2_ensembl.r {input} {output} > {log} 2>&1"

# GENOME
# ------------------------------------------------------------------------------
rule get_annot_ref:
    output:
        config["ANNOT_PATH"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/get_annot_ref.log"
    params:
        url=config["ANNOT_URL"]
    shell:
        "wget -nv -m -nH -O {output} {params.url} > {log} 2>&1"

rule gene_ids_file:
    input:
        config["ANNOT_PATH"]
    output:
        gene2tx=config["ENSID_PATH"],
        gene2tx2ver=config["ENSIDVER_PATH"],
        gene2tx2full=config["ENSIDFULLVER_PATH"]
    group:
        "downloadRefFiles" 
    log:
        "output/logs/ref/gene_ids_file.log"
    shell:
        "Rscript workflow/scripts/gene_ids_file.r > {log} 2>&1 && "
        "Rscript workflow/scripts/gene_ids_versions_file.r >> {log} 2>&1 && "
        "Rscript workflow/scripts/gene_ids_full_versions_file.r {input} {output.gene2tx2full} >> {log} 2>&1 "

rule get_trans_ref:
    output:
        config["TRANS_PATH"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/get_trans_ref.log"
    params:
        url=config["TRANS_URL"]
    shell:
        "wget -nv -m -nH -O {output} {params.url} > {log} 2>&1"

rule get_genom_ref:
    output:
        config["GENOM_PATH"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/get_genom_ref.log"
    params:
        url=config["GENOM_URL"]
    shell:
        "wget -nv -m -nH -O {output} {params.url} > {log} 2>&1"

# BOWTIE2
# ------------------------------------------------------------------------------
rule get_bowtie_idx:
    output:
        config["BOWTIE_IDX_PATH"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/get_bowtie_idx.log"
    params:
        url=config["BOWTIE_IDX_URL"]
    shell:
        "wget -nv -m -nH -O {output} {params.url} > {log} 2>&1"

rule extract_bowtie_idx:
    input:
        # config["BOWTIE_IDX_PATH"]
        rules.get_bowtie_idx.output
    output:
        config["BOWTIE_IDX_FILE1"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/extract_bowtie_idx.log"
    params:
        dir=config["BOWTIE_IDX_DIR"],
        ext=os.path.splitext(config["BOWTIE_IDX_PATH"])[1]
    shell:
        """
        if [ "{params.ext}" == ".zip"  ]; then
          unzip -j -n {input} -d {params.dir} > {log} 2>&1
        else
          tar xzf {input} -C {params.dir} > {log} 2>&1
        fi
        """

# PICARD
# ------------------------------------------------------------------------------
rule get_picard_ref:
    output:
        config["PICARD_REFFLAT_PATH"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/get_picard_ref.log"
    params:
        url=config["PICARD_REFFLAT_URL"]
    shell:
        "wget -nv -m -nH -O {output} {params.url} > {log} 2>&1"

rule extract_picard_ref:
    input:
        config["PICARD_REFFLAT_PATH"]
    output:
        config["PICARD_REFFLAT_FILE"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/extract_picard_ref.log"
    params:
        dir=config["PICARD_DIR"],
        ext=os.path.splitext(config["PICARD_REFFLAT_PATH"])[1]
    shell:
        "gzip -d {input} > {log} 2>&1"

# MSIGDB
# ------------------------------------------------------------------------------
rule get_msigdb_ref:
    output:
        path_hall=config["MSIGDB_HALLMARK_PATH"],
        path_posi=config["MSIGDB_POSITIONAL_PATH"],
        path_cure=config["MSIGDB_CURATED_PATH"],
        path_cucp=config["MSIGDB_CURATED_CP_PATH"],
        path_regu=config["MSIGDB_REGULATORY_PATH"],
        path_cell=config["MSIGDB_CELLS_PATH"]
    group:
        "downloadRefFiles"
    log:
        "output/logs/ref/get_msigdb_ref.log"
    params:
        url_hall=config["MSIGDB_HALLMARK_URL"],
        url_posi=config["MSIGDB_POSITIONAL_URL"],
        url_cure=config["MSIGDB_CURATED_URL"],
        url_cucp=config["MSIGDB_CURATED_CP_URL"],
        url_regu=config["MSIGDB_REGULATORY_URL"],
        url_cell=config["MSIGDB_CELLS_URL"]
    shell:
        "wget -nv -m -nH -nv -O {output.path_hall} {params.url_hall} > {log} 2>&1 && "
        "wget -nv -m -nH -nv -O {output.path_posi} {params.url_posi} > {log} 2>&1 && "
        "wget -nv -m -nH -nv -O {output.path_cure} {params.url_cure} > {log} 2>&1 && "
        "wget -nv -m -nH -nv -O {output.path_cucp} {params.url_cucp} > {log} 2>&1 && "
        "wget -nv -m -nH -nv -O {output.path_regu} {params.url_regu} > {log} 2>&1 && "
        "wget -nv -m -nH -nv -O {output.path_cell} {params.url_cell} > {log} 2>&1 "