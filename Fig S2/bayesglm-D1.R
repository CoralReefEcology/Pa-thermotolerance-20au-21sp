setwd("D:/BaiduSyncdisk/Paper_submitting/Ms-Pa-xisha-2020au-2021sp/Physiological-parameter-data")
Pa_data = read.csv("Pa-Physical-XS-20au-21sp-total.csv", header = TRUE, row.names = 1)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(visreg)
library(arm)

d1_ratio_dhw=bayesglm(cbind(D_numbers, C_numbers)~DHW, data=Pa_data, family=quasibinomial)
summary(d1_ratio_dhw)
# Call:
# bayesglm(formula = cbind(D_numbers, C_numbers) ~ DHW, family = quasibinomial,
# data = Pa_data)
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.814306 0.209378 3.889 0.000196 ***
# DHW 0.009491 0.029928 0.317 0.751899
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasibinomial family taken to be 30705.34)
# 
# Null deviance: 2947674 on 88 degrees of freedom
# Residual deviance: 2944482 on 87 degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 7
D1DHW = visreg(d1_ratio_dhw, gg=TRUE, ylab="Proportion of Durusdinium D1", scale="response", line=list(col="darkgreen"), whitespace=0.8)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_y_continuous(limits=c(0.55,0.85))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 2C.pdf", D1DHW, width = 5, height = 3, units = "in")

d1_ratio_dhw=glm(cbind(D_numbers, C_numbers)~DHW, data=Pa_data, family=quasibinomial)
summary(d1_ratio_dhw)
# Call:
# glm(formula = cbind(D_numbers, C_numbers) ~ DHW, family = quasibinomial,
# data = Pa_data)
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.813301 0.210807 3.858 0.000219 ***
# DHW 0.009792 0.030393 0.322 0.748078
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasibinomial family taken to be 30704.82)
# 
# Null deviance: 2947674 on 88 degrees of freedom
# Residual deviance: 2944479 on 87 degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 4

d1_ratio_depth=bayesglm(cbind(D_numbers, C_numbers)~Depth, data=Pa_data, family=quasibinomial)
summary(d1_ratio_depth)
# Call:
# bayesglm(formula = cbind(D_numbers, C_numbers) ~ Depth, family = quasibinomial,
# data = Pa_data)
# 
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)
# (Intercept) 0.99675 0.25752 3.871 0.000209 ***
# Depth -0.02685 0.03946 -0.681 0.497984
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasibinomial family taken to be 30852.35)
# 
# Null deviance: 2947674 on 88 degrees of freedom
# Residual deviance: 2933121 on 87 degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 7
D1Depth = visreg(d1_ratio_depth, gg=TRUE, ylab="Proportion of Durusdinium D1", scale="response", line=list(col="darkgreen"), whitespace=0.8)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_y_continuous(limits=c(0.3,0.9))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. S2.pdf", D1Depth, width = 5, height = 3, units = "in")

d1_ratio_depth=glm(cbind(D_numbers, C_numbers)~Depth, data=Pa_data, family=quasibinomial)
summary(d1_ratio_depth)
# Call:
#   glm(formula = cbind(D_numbers, C_numbers) ~ Depth, family = quasibinomial, 
#       data = Pa_data)
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  1.00144    0.26017   3.849 0.000226 ***
#   Depth       -0.02768    0.04005  -0.691 0.491375    
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for quasibinomial family taken to be 30862.11)
# 
# Null deviance: 2947674  on 88  degrees of freedom
# Residual deviance: 2933108  on 87  degrees of freedom
# AIC: NA
# 
# Number of Fisher Scoring iterations: 4