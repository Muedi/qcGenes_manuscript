# Overlooked poor-quality patient samples in sequencing data impair reproducibility of published clinically relevant datasets
Workflows for correlation of Expression and ChIP-seq peaks with file quality.

This repository represents the snakemake workflows used in the paper "Overlooked poor-quality patient samples in sequencing data impair reproducibility of published clinically relevant datasets". 
The manuscript can be found here: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03331-6

The workflow has been run on Ubuntu 18 and 20, with different snakemake versions. The version of the papers output can be seeen in the qcGenes subfolder readme. 
The ChIPseq workflow has been originally implemented by Jannik Möllmann (@jmoellmann, https://github.com/jmoellmann), but adapted to newer software versions and datasets for the paper. 

#################### This workflow is a beta version! #####################

We plan to release a quality-informed RNAseq workflow, that incorporates our quality software seqQscorer, as well as the analyses shown in this paper and the related batch effect work. 
See: 
SeqQscorer, Quality Software: https://github.com/salbrec/seqQscorer
Batch effect paper: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-04775-y
