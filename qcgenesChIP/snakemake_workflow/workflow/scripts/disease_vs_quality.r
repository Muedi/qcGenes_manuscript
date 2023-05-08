library(ggplot2)
suppressMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
options(ggplot2.discrete.colour= c("#295D8A", "#A41720", "#4E4459"))

args = commandArgs(trailingOnly=TRUE)

df_path <- args[1]
scores_path <-args[2]
geo_id <- args[3]
  
df <- read.csv(file = df_path) %>% filter(is_sample == 1)
df$scores <- read.csv(file = scores_path, sep = "\t", header = FALSE)[[2]]

res <- cor.test(df$is_healthy, df$scores)

plot <- ggplot(data = df)  + 
  aes(x = reorder(Run, scores), y = scores, fill = condition) + 
  geom_bar(stat = "identity") +
  labs(x = "SRA Run", y = "Probability of Low Quality",
       title = sprintf("Correlation of Quality and Disease for %s", geo_id),
       subtitle = sprintf("Pearsons R: %1.3f, P-Value: %1.3e", res$estimate, res$p.value)) + 
  coord_flip()

ggsave(sprintf("output/plots/disease_vs_quality_%s.svg", geo_id), plot = plot, width = 10, height = 7)
# ggsave(sprintf("output/plots/disease_vs_quality_%s.png", geo_id), plot = plot, width = 10, height = 7)
