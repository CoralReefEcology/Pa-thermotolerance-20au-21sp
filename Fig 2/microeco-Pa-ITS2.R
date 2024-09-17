library(microeco)
set.seed(123)
library(ggplot2)
library(ggalluvial)
theme_set(theme_bw())
## reading data
asvtable <- read.csv("asv-table.csv", row.names = 1, header = TRUE, sep = ",")
asvtaxon <- read.csv("asv-tax.csv", row.names = 1, header = TRUE, sep = ",")
asvsample <- read.csv("metadata.tsv", row.names = 1, header = TRUE, sep = "\t")
dataset <- microtable$new(sample_table = asvsample, otu_table = asvtable, tax_table = asvtaxon)
dataset

dataset$sample_sums() %>% range
dataset$rarefy_samples(sample.size = 20000)
dataset$sample_sums() %>% range
dataset

##community structure
dataset$cal_abund()
class(dataset$taxa_abund)
dir.create("taxa_abund")
dataset$save_abund(dirpath = "taxa_abund")

t1 <- trans_abund$new(dataset = dataset, taxrank = "Species", ntaxa = 8)
p2 = t1$plot_bar(bar_type = "notfull", facet = "Time", use_alluvium = TRUE, clustering = TRUE, xtext_keep = FALSE, color_values = RColorBrewer::brewer.pal(8, "Set2"))
ggsave("Fig. 2A.pdf", p2, width = 6, height = 4, units = "in")

##RDA analysis
Pa_env_info <- read.csv("Pa_env_info.csv", row.names = 1, header = TRUE)
t2 <- trans_env$new(dataset = dataset, add_data = Pa_env_info)
t2$cal_ordination(method = "RDA", taxa_level = "Species")
t2$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
p3 = t2$plot_ordination(plot_color = "Time")
ggsave("Fig. 2B.pdf", p3, width = 6, height = 4, units = "in")
t2$cal_mantel(use_measure = "bray")
head(t2$res_mantel)
# by_group Variables mantel type Correlation method Correlation coefficient p.value p.adjusted Significance
# 1      All       DHW mantel test            pearson              0.07307839   0.003      0.006           **
# 2      All     Depth mantel test            pearson             -0.05834571   0.767      0.767 