library(microeco)
library(file2meco)
library(ape)
library(magrittr)
set.seed(123)
library(ggplot2)
library(cowplot)
theme_set(theme_bw())

## Reading data
dataset <- qiime2meco("table-filtered-contam-tree.qza", sample_table = "metadata.tsv", taxonomy_table = "classification.qza", phylo_tree = "asvs-tree.qza", rep_fasta = "rep-seqs-final.qza", auto_tidy = TRUE)
dataset
dataset$sample_sums() %>% range
dataset$rarefy_samples(sample.size = 27000)
dataset$sample_sums() %>% range

## alpha diversity
dataset$cal_abund()
class(dataset$taxa_abund)
dir.create("taxa_abund")
dataset$save_abund(dirpath = "taxa_abund")

t1 <- trans_abund$new(dataset = dataset, taxrank = "Phylum", ntaxa = 10)
community = t1$plot_bar(others_color = "grey70", facet = "Time", xtext_keep = FALSE, legend_text_italic = FALSE)
ggsave("Fig. S3A.pdf", community, width = 12, height = 5, units = "in")

library(picante)
dataset$cal_alphadiv(PD = TRUE)
class(dataset$alpha_diversity)
dir.create("alpha_diversity")
dataset$save_alphadiv(dirpath = "alpha_diversity")

t2 <- trans_alpha$new(dataset = dataset, group = "Time")
t2$data_stat
t2$alpha_stat[1:5, ]
t2$cal_diff(method = "KW")
t2$res_diff
t2$res_alpha_diff[1:5, ]
F1 = t2$plot_alpha(add_letter = TRUE, measure = "Chao1")
F2 = t2$plot_alpha(add_letter = TRUE, measure = "InvSimpson")
F3 = t2$plot_alpha(add_letter = TRUE, measure = "Shannon")
F4 = t2$plot_alpha(add_letter = TRUE, measure = "Pielou")
alpha = plot_grid(F1, F2, F3, F4, labels = c('A', 'B', 'C', 'D'))
ggsave("Fig. S3B.pdf", alpha, width = 7, height = 7, units = "in")

Pa_env_info <- read.csv("env_data.csv", header = TRUE, row.names = 1)
t1 <- trans_env$new(dataset = dataset, add_data = Pa_env_info)
t1$cal_cor(add_abund_table = dataset$alpha_diversity)
AlphaEnvs = t1$plot_cor(cluster_ggplot = "both")
ggsave("Fig. 3A.pdf", AlphaEnvs, width = 6, height = 4.5, units = "in")

## beta diversity
library(GUniFrac)
dataset$cal_betadiv(unifrac = TRUE)
class(dataset$beta_diversity)
dir.create("beta_diversity")
dataset$save_betadiv(dirpath = "beta_diversity")

t1 <- trans_beta$new(dataset = dataset, group = "Time", measure = "bray")
t1$cal_ordination(method = "PCoA")
class(t1$res_ordination)
Beta = t1$plot_ordination(plot_color = "Time", plot_shape = "Time", plot_type = c("point", "ellipse"))
ggsave("Fig. S3C.pdf", Beta, width = 6, height = 4.5, units = "in")
t1$cal_manova(manova_all = TRUE)
t1$res_manova
# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999
# 
# adonis2(formula = use_formula, data = metadata)
# Df SumOfSqs      R2      F Pr(>F)    
# Time      1    5.230 0.15973 16.538  0.001 ***
#   Residual 87   27.516 0.84027                  
# Total    88   32.746 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

t2 <- trans_env$new(dataset = dataset, add_data = Pa_env_info)
t2$cal_ordination(method = "RDA", taxa_level = "Genus")
t2$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)
RDA = t2$plot_ordination(plot_color = "Time")
ggsave("Fig. 3B.pdf", RDA, width = 6, height = 4.5, units = "in")
t2$cal_mantel(use_measure = "bray")
head(t2$res_mantel)
#   by_group Variables mantel type Correlation method Correlation coefficient p.value p.adjusted Significance
# 1      All     Depth mantel test            pearson              0.02509833   0.202      0.202             
# 2      All       DHW mantel test            pearson              0.44389423   0.001      0.003           **
# 3      All  D1_ratio mantel test            pearson              0.11087814   0.002      0.003           **

## differential bacteria
Pa_env_info <- read.csv("env_data_energy.csv", header = TRUE, row.names = 1)
t2 <- trans_env$new(dataset = dataset, add_data = Pa_env_info)
t2$cal_cor(use_data = "Genus", p_adjust_method = "fdr", p_adjust_type = "Env")
t2$plot_cor(filter_feature = c("", "*", "**"))
t3 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Time", taxa_level = "Genus")
t2$cal_cor(use_data = "other", p_adjust_method = "fdr", other_taxa = t3$res_diff$Taxa[1:30])
t2$plot_cor()
Sig_genus_envs = t2$plot_cor(pheatmap = TRUE, filter_feature = "", color_palette = rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))
ggsave("Fig. 3C.pdf", Sig_genus_envs, width = 6, height = 6, units = "in")

## The relationship between metabolism trait and environmental factors
t2 <- trans_func$new(dataset)
t2$cal_spe_func(prok_database = "FAPROTAX")
t2$res_spe_func[1:5, 1:2]
t2$cal_spe_func_perc(abundance_weighted = FALSE)
t2$res_spe_func[1:5, 1:2]

t3 <- trans_env$new(dataset = dataset, add_data = Pa_env_info)
t3$cal_cor(add_abund_table = t2$res_spe_func_perc, cor_method = "spearman")
t3$plot_cor(pheatmap = TRUE)
MetabolismEnvs = t3$plot_cor(pheatmap = TRUE, filter_feature = "")
ggsave("Fig. 3D.pdf", MetabolismEnvs, width = 6, height = 6, units = "in")




