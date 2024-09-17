Pa_data = read.csv("Pa-Physical-XS-20au-21sp-total.csv", header = TRUE, row.names = 1)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(car)
library(patchwork)

shapiro.test(Pa_data$Protein)
# Shapiro-Wilk normality test
# data:  Pa_data$Protein
# W = 0.98353, p-value = 0.3215
leveneTest(Pa_data$Protein, as.factor(Pa_data$Time))
# Levene's Test for Homogeneity of Variance (center = median)
# Df F value Pr(>F)
# group 1 0.0304 0.862
# 87
t.test(Protein~Time, var.equal=T, data=Pa_data)
# Two Sample t-test
# 
# data: Protein by Time
# t = -10.282, df = 87, p-value < 2.2e-16
# alternative hypothesis: true difference in means between group 20au and group 21sp is not equal to 0
# 95 percent confidence interval:
# -0.3242498 -0.2191990
# sample estimates:
# mean in group 20au mean in group 21sp
# 0.3500851 0.6218095


shapiro.test(Pa_data$Carbohydrate)
# Shapiro-Wilk normality test
# data: Pa_data$Carbohydrate
# W = 0.81603, p-value = 3.666e-09
shapiro.test(Pa_data$Lipid)
# Shapiro-Wilk normality test
# data: Pa_data$Lipid
# W = 0.96369, p-value = 0.01368
kruskal.test(Carbohydrate~Time, data = Pa_data)
# Kruskal-Wallis rank sum test
# data: Carbohydrate by Time
# Kruskal-Wallis chi-squared = 26.135, df = 1, p-value = 3.183e-07
kruskal.test(Lipid~Time, data = Pa_data)
# Kruskal-Wallis rank sum test
# data: Lipid by Time
# Kruskal-Wallis chi-squared = 29.379, df = 1, p-value = 5.953e-08
nutrient_protein = ggplot(Pa_data, aes(x=Time, y=Protein, fill=Time))+geom_bar(stat="summary", fun=mean, color="black", position=position_dodge())+stat_summary(fun.data='mean_se', geom="errorbar", colour="black", width=0.25)+scale_fill_discrete(guide=FALSE)+ylab("Protein content")+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))
nutrient_lipid = ggplot(Pa_data, aes(x=Time, y=Lipid, fill=Time))+geom_bar(stat="summary", fun=mean, color="black", position=position_dodge())+stat_summary(fun.data='mean_se', geom="errorbar", colour="black", width=0.25, position=position_dodge(0.9))+scale_fill_discrete(guide=FALSE)+ylab("Lipid content")+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))
nutrient_carbohydrate = ggplot(Pa_data, aes(x=Time, y=Carbohydrate, fill=Time))+geom_bar(stat="summary", fun=mean, color="black", position=position_dodge())+stat_summary(fun.data='mean_se', geom="errorbar", colour="black", width=0.25, position=position_dodge(0.9))+scale_fill_discrete(guide=FALSE)+ylab("Carbohydrate content")+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))
nutrient <- nutrient_protein+nutrient_carbohydrate+nutrient_lipid
ggsave("Fig. S3.pdf", nutrient, width = 5, height = 3, units = "in")
