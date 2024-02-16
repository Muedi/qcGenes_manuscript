suppressMessages(library(fgsea, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(tidyverse, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(yaml, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggpmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
library(gplots)
library(stringr)
library(ComplexHeatmap)
library(grDevices)
library(circlize)

options(readr.num_columns = 0)
options(ggplot2.discrete.colour= c("#295D8A", "#A41720", "#4E4459"))
options(ggplot2.discrete.fill= c("#295D8A", "#A41720", "#4E4459"))


sra.path    <- file.path(".", "config", "metadata", "Mega_SraRunTable.csv")
datasets.path <- file.path(".", "config", "metadata", "Datasets.csv")

# ==============================================================================
sra.table <- read_csv(sra.path, show_col_types=F) 
# dataids.list <- unique(sra.table$GEO_Series)
datasets.table <- read_csv(datasets.path, show_col_types=F) %>% 
  rename(mesh_terms=`mesh_terms (pipe-separated list)`)
dataids.list <- (datasets.table %>% filter(select==1 & batches==0))$GEO_Series
#dataids.list <- (datasets.table %>% filter(select==1 & batches==0 & curators=="JF"))$GEO_Series
# remove subsets
dataids.list <- dataids.list[!grepl("sub", dataids.list)]

df <- data.frame(matrix(NA, nrow=length(dataids.list), ncol=3))
   
colnames(df) <- c("QI", "p.krusk", "size")
rownames(df) <- dataids.list

for (i in 1:length(dataids.list)) {
   p.low.path <- file.path(".", "output", "main", "scores", paste0(dataids.list[[i]], ".scores.txt"))
   p.low <- read_delim(p.low.path, delim="\t", col_names=F)
   p.low <- p.low %>% dplyr::select(-X3)
   colnames(p.low) <- c("Run", "p_low")

    sub.sra <- sra.table %>% 
        filter(GEO_Series==dataids.list[[i]]) %>%
        filter(Selected==1)

    sub.sra <- left_join(sub.sra, p.low, by="Run")
    p <- sub.sra$p_low
    group <- sub.sra$group %>% as.factor()
    group_num <- as.numeric(group)
    krusk.pval.p.group <- kruskal.test(p ~ group)$p.value
    p.group.cor <- abs(cor(p, group_num))


    df[dataids.list[[i]], "QI"] = p.group.cor
    df[dataids.list[[i]], "p.krusk"] = krusk.pval.p.group
    df[dataids.list[[i]], "size"] = length(p)
}

df <- round(df, 3)

metrics <- as.matrix(df[c("QI", "p.krusk")])
sizes <- as.matrix(df["size"])

metrics_hm <- Heatmap(metrics, 
                    name = "Metrics",
                    col = colorRamp2(c(0, 1), c("white", "darkorange")), 
                    column_names_gp=gpar(fontsize=10),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.3f", metrics[i, j]), x, y, gp = gpar(fontsize = 10))}
        )
sizes_hm <- Heatmap(sizes,
                    name = "Size",
                    col = colorRamp2(c(0, max(sizes)), c("white", "steelblue")), 
                    column_names_gp=gpar(fontsize=10),
                    cell_fun = function(j, i, x, y, width, height, fill) {
                            grid.text(sprintf("%.0f", sizes[i, j]), x, y, gp = gpar(fontsize = 10))}
                    )

pdf("review/heatmap-pvals.pdf", width = 7, height = 7) #, units = "in", res=300)
metrics_hm + sizes_hm
dev.off()

# png("review/heatmap.png")
# heatmap.2(as.matrix(df),
#           dendrogram = "none",  # To suppress row and column dendrograms
#           main = "Simple Heatmap",
#           xlab = "Columns",
#           ylab = "Rows",
#           key = TRUE,  # Show color key
#           key.title = NA,
#           cellnote = df,
#           notecol = "black",  # Set text color
# )
# dev.off()


# seqqscorer features

# df <- data.frame()
# for (i in 1:length(dataids.list)) {
#     sub.sra <- sra.table %>% 
#         filter(GEO_Series==dataids.list[[i]]) %>%
#         filter(Selected==1)
#     run_ids <- sub.sra$Run
#     for (j in 1:length(run_ids)) {
#           p.low.path <- file.path(".", "output", "qc", dataids.list[[i]], "scores",  paste0(run_ids[[j]], ".score.txt"))
#           p.low <- read_delim(p.low.path, delim="\t", col_names=F)
#         scorer.features.path <- file.path(".", "output", "qc", dataids.list[[i]], "scores",  paste0(run_ids[[j]], ".features.txt"))
#         scorer.features <- read_delim(scorer.features.path, delim="\t")
#         scorer.features["dataset"] <- dataids.list[[i]]
#         scorer.features["p_low"] <- p.low[[2]]
#         if (nrow(df)==0) {
#            df <- scorer.features
#         } else {
#             df <- rbind(df, scorer.features)
#         }
#     }
# }
# write_csv(df, "review/qc_metrics.csv")

df <- read_csv("review/qc_metrics.csv") %>% select(-1)

png("qcfeatures_heatmapfastqc.png")
heatmap.2(as.matrix(fastqc),
          dendrogram = "none",  # To suppress row and column dendrograms
          main = "Simple Heatmap",
          xlab = "Columns",
          ylab = "Rows",
          key = TRUE,  # Show color key
          key.title = NA,
          cellnote = df,
          notecol = "black",  # Set text color
)
dev.off()



# Sample data creation (replace with your actual data)
set.seed(123)
accessions <- c("A", "B", "C", "D", "E")
fastqc_data <- df %>%
    select("accession", colnames(df)[str_detect(colnames(df), "FastQC")]) %>%
    mutate_at(vars(colnames(df)[str_detect(colnames(df), "FastQC")]), as.factor)
mapping_rates <- df %>% select("accession", colnames(df)[str_detect(colnames(df), "Bowtie")])
quality_scores <- df %>% select("accession", "p_low")

# Merge the dataframes based on the Accession column
merged_data <- merge(merge(fastqc_data, mapping_rates, by = "accession"), quality_scores, by = "accession")

# Create a matrix with the desired data for each heatmap
fastqc_matrix <- as.matrix(merged_data[, grep("FastQC", names(merged_data))])
mapping_matrix <- as.matrix(merged_data[, grep("Bowtie", names(merged_data))])
quality_matrix <- as.matrix(merged_data$p_low)

# Set up column annotation
column_annotation <- data.frame(Mapping_Rate = merged_data$accession, Quality_Score = merged_data$p_low)

# Create separate heatmaps
heatmap_fastqc <- Heatmap(fastqc_matrix, name = "FastQC", col = c("0"="red", "1"="yellow", "2"="green"), column_names_gp=gpar(fontsize=10))
heatmap_mapping <- Heatmap(mapping_matrix, name = "Mapping Rate", column_names_gp=gpar(fontsize=10))
heatmap_quality <- Heatmap(quality_matrix, name = "P_low", col = colorRamp2(c(0, 1), c("white", "darkorange")), column_names_gp=gpar(fontsize=10))

# Combine heatmaps into a list

# Draw the aligned heatmaps
png("review/heatmaps-qc.png", width = 10, height = 10, units = "in", res=300)
heatmap_quality + heatmap_fastqc + heatmap_mapping
dev.off()

pdf("review/heatmaps-qc.pdf", width = 10, height = 10) #, units = "in", res=300)
heatmap_quality + heatmap_fastqc + heatmap_mapping
dev.off()

hm_oall <- grid.grabExpr(draw(heatmap_quality + heatmap_fastqc + heatmap_mapping))

merged_data["p_low_binned"] = cut(merged_data$p_low, breaks=4, lables=F)

dense_plot_a = ggplot(merged_data, aes(x = BowtieMI_uniquely, fill = factor(p_low_binned))) +
  geom_histogram(binwidth = 10, color = "black", alpha = 0.7, position = "stack") +
  # geom_density(alpha = 0.5) +
  labs(title = "Distribution of uniquely mapped reads",
       x = "BowtieMI_uniquely",
       y = "Count",
       fill = "p_low")

ggsave("review/density_plow_uniquely.png", plot=dense_plot_a, width=5, height=5, dpi=300)

dense_plot_b <- ggplot(merged_data, aes(x = FastQC_Sequence_Duplication_Levels, fill = factor(p_low_binned))) +
  geom_histogram(binwidth = 0.1, color = "black", alpha = 0.7, stat='count', position = "stack") +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution FastQC's Sequence Duplication Levels",
       x = "FastQC_Sequence_Duplication_Levels",
       y = "Count",
       fill = "p_low_binned")

ggsave("review/density_plow_seq_duplic.png", plot=dense_plot_b, width=5, height=5, dpi=300)


figure <- ggarrange(
  dense_plot_a,
  dense_plot_b,
  labels = c("A", "B"),
  # title = "Quality features against binned P_low", 
  font.label = list(size = 18, color = "black", face = "bold", family = NULL),
  common.legend = TRUE, legend = "bottom",
  ncol = 2, nrow = 1)
path <- file.path("review", "quality_histos.pdf")
ggsave(filename=path, plot=figure, width=10, height=5)

figure <- ggarrange(
  hm_oall,
  ggarrange( 
    dense_plot_a,
    dense_plot_b,
    labels = c("B", "C"),
    # title = "Quality features against binned P_low", 
    font.label = list(size = 18, color = "black", face = "bold", family = NULL),
    common.legend = TRUE, legend = "bottom", 
    ncol = 2, nrow = 1),
  ncol = 1, nrow = 2,
  heights = c(2, 1),
  labels = c("A", ""),
  font.label = list(size = 18, color = "black", face = "bold", family = NULL)
)
path <- file.path("review", "figure-S4-quality_histos_combi.pdf")
ggsave(filename=path, plot=figure, width=10, height=15)

