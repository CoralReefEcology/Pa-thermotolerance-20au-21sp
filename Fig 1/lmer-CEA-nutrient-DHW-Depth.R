Pa_data = read.csv("Pa-Physical-XS-20au-21sp-total.csv", header = TRUE, row.names = 1)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(lmerTest)
library(performance)
library(ggpmisc)
library(visreg)
library(arm)
library(sjPlot)

## Effects of DHW and Depth on the celluar energy allocation of the coral host
CEA_DHW_Depth=lmer(log(CEA)~DHW+Depth+(1|Station), data=Pa_data)
check_model(CEA_DHW_Depth)
check_normality(CEA_DHW_Depth)
# OK: residuals appear as normally distributed (p = 0.810).
check_heteroscedasticity(CEA_DHW_Depth)
# OK: Error variance appears to be homoscedastic (p = 0.989).
check_outliers(CEA_DHW_Depth)
# OK: No outliers detected.
# - Based on the following method and threshold: cook (0.8).
# - For variable: (Whole model)
check_collinearity(CEA_DHW_Depth)
# Check for Multicollinearity
# Low Correlation
# Term VIF VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# DHW 1.03 [1.00, 18.90] 1.02 0.97 [0.05, 1.00]
# Depth 1.03 [1.00, 18.90] 1.02 0.97 [0.05, 1.00]
summary(CEA_DHW_Depth)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: log(CEA) ~ DHW + Depth + (1 | Station)
# Data: Pa_data
# 
# REML criterion at convergence: 102.7
# 
# Scaled residuals:
# Min 1Q Median 3Q Max
# -2.13423 -0.54820 -0.05468 0.56334 2.33913
# 
# Random effects:
# Groups Name Variance Std.Dev.
# Station (Intercept) 0.05668 0.2381
# Residual 0.12708 0.3565
# Number of obs: 89, groups: Station, 19
# 
# Fixed effects:
# Estimate Std. Error df t value Pr(>|t|)
# (Intercept) 4.90812 0.09249 43.98652 53.068 <2e-16 ***
# DHW -0.01335 0.01083 62.36589 -1.233 0.222
# Depth -0.02829 0.01073 83.79406 -2.636 0.010 **
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) DHW
# DHW -0.242
# Depth -0.537 -0.179
tab_model(CEA_DHW_Depth)
figCEADHW <- visreg(CEA_DHW_Depth, "DHW", ylab="log(CEA)", line=list(col="brown"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 1A.pdf", figCEADHW, width = 5, height = 3, units = "in")

Ea_DHW_Depth=lmer(log(Ea)~DHW+Depth+(1|Station), data=Pa_data)
check_model(Ea_DHW_Depth)
check_collinearity(Ea_DHW_Depth)
# # Check for Multicollinearity
# Low Correlation
# Term VIF VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# DHW 1.04 [1.00, 11.49] 1.02 0.96 [0.09, 1.00]
# Depth 1.04 [1.00, 11.49] 1.02 0.96 [0.09, 1.00]
check_normality(Ea_DHW_Depth)
# OK: residuals appear as normally distributed (p = 0.120).
check_outliers(Ea_DHW_Depth)
# OK: No outliers detected.- Based on the following method and threshold: cook (0.8).
# - For variable: (Whole model)
check_heteroscedasticity(Ea_DHW_Depth)
# OK: Error variance appears to be homoscedastic (p = 0.766).
summary(Ea_DHW_Depth)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: log(Ea) ~ DHW + Depth + (1 | Station)
# Data: Pa_data
# 
# REML criterion at convergence: 86.8
# 
# Scaled residuals:
# Min 1Q Median 3Q Max
# -2.45115 -0.42706 0.08145 0.41861 2.42292
# 
# Random effects:
# Groups Name Variance Std.Dev.
# Station (Intercept) 0.06175 0.2485
# Residual 0.10165 0.3188
# Number of obs: 89, groups: Station, 19
# 
# Fixed effects:
# Estimate Std. Error df t value Pr(>|t|)
# (Intercept) 4.231365 0.088161 34.109240 47.996 < 2e-16 ***
# DHW -0.051835 0.010105 68.039868 -5.130 2.60e-06 ***
# Depth -0.043373 0.009707 81.327912 -4.468 2.52e-05 ***
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) DHW
# DHW -0.215
# Depth -0.509 -0.189
tab_model(Ea_DHW_Depth)
figEaDHW <- visreg(Ea_DHW_Depth, "DHW", ylab="log(Ea)", line=list(col="brown"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 1B.pdf", figEaDHW, width = 5, height = 3, units = "in")

Ec_DHW_Depth=lmer(log(Ec)~DHW+Depth+(1|Station), data=Pa_data)
check_model(Ec_DHW_Depth)
check_collinearity(Ec_DHW_Depth)
# Check for Multicollinearity
# Low Correlation
# Term VIF VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# DHW 1.04 [1.00, 14.49] 1.02 0.97 [0.07, 1.00]
# Depth 1.04 [1.00, 14.49] 1.02 0.97 [0.07, 1.00]
check_normality(Ec_DHW_Depth)
# OK: residuals appear as normally distributed (p = 0.267).
check_outliers(Ec_DHW_Depth)
# OK: No outliers detected.
# - Based on the following method and threshold: cook (0.8).
# - For variable: (Whole model)
check_heteroscedasticity(Ec_DHW_Depth)
# OK: Error variance appears to be homoscedastic (p = 0.765).
summary(Ec_DHW_Depth)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: log(Ec) ~ DHW + Depth + (1 | Station)
# Data: Pa_data
# 
# REML criterion at convergence: 60.5
# 
# Scaled residuals:
# Min 1Q Median 3Q Max
# -2.90227 -0.43674 0.07306 0.60876 2.30921
# 
# Random effects:
# Groups Name Variance Std.Dev.
# Station (Intercept) 0.03983 0.1996
# Residual 0.07635 0.2763
# Number of obs: 89, groups: Station, 19
# 
# Fixed effects:
# Estimate Std. Error df t value Pr(>|t|)
# (Intercept) 6.227379 0.073951 44.529773 84.210 < 2e-16 ***
# DHW -0.040594 0.008580 68.875245 -4.731 1.15e-05 ***
# Depth -0.013740 0.008366 83.445705 -1.642 0.104
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) DHW
# DHW -0.228
# Depth -0.523 -0.184
tab_model(Ec_DHW_Depth)
figEcDHW <- visreg(Ec_DHW_Depth, "DHW", ylab="log(Ec)", line=list(col="brown"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 1C.pdf", figEcDHW, width = 5, height = 3, units = "in")

## Effects of DHW and Depth on the protein, carbohydrate, and lipid content in the coral host
Protein_DHW_Depth=lmer(Protein~DHW+Depth+(1|Station), data=Pa_data)
check_model(Protein_DHW_Depth)
check_normality(Protein_DHW_Depth)
# OK: residuals appear as normally distributed (p = 0.225).
check_heteroscedasticity(Protein_DHW_Depth)
# Warning: Heteroscedasticity (non-constant error variance) detected (p < .001).Warning message:
# In sqrt(insight::get_deviance(x)/(insight::n_obs(x) - sum(!is.na(estimates)))) :
# NaNs produced
check_outliers(Protein_DHW_Depth)
# OK: No outliers detected.
# - Based on the following method and threshold: cook (0.8).
# - For variable: (Whole model)
check_collinearity(Protein_DHW_Depth)
# Check for Multicollinearity
# Low Correlation
# Term VIF VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# DHW 1.04 [1.00, 6.31] 1.02 0.96 [0.16, 1.00]
# Depth 1.04 [1.00, 6.31] 1.02 0.96 [0.16, 1.00]
summary(Protein_DHW_Depth)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Protein ~ DHW + Depth + (1 | Station)
# Data: Pa_data
# 
# REML criterion at convergence: -112.9
# 
# Scaled residuals:
# Min 1Q Median 3Q Max
# -3.4914 -0.5993 0.0060 0.6348 2.6071
# 
# Random effects:
# Groups Name Variance Std.Dev.
# Station (Intercept) 0.009448 0.09720
# Residual 0.009258 0.09622
# Number of obs: 89, groups: Station, 19
# 
# Fixed effects:
# Estimate Std. Error df t value Pr(>|t|)
# (Intercept) 0.640926 0.030312 32.693994 21.145 < 2e-16 ***
# DHW -0.025880 0.003245 83.095540 -7.976 7.2e-12 ***
# Depth -0.006246 0.002981 79.323982 -2.095 0.0393 *
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) DHW
# DHW -0.172
# Depth -0.453 -0.206
tab_model(Protein_DHW_Depth)
figProDHW <- visreg(Protein_DHW_Depth, "DHW", ylab="Protein content", line=list(col="brown"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 1D.pdf", figProDHW, width = 5, height = 3, units = "in")

Carbohydrate_DHW_Depth=lmer(log(Carbohydrate)~DHW+Depth+(1|Station), data=Pa_data)
check_model(Carbohydrate_DHW_Depth)
check_normality(Carbohydrate_DHW_Depth)
# Warning: Non-normality of residuals detected (p = 0.047).
check_heteroscedasticity(Carbohydrate_DHW_Depth)
# OK: Error variance appears to be homoscedastic (p = 0.996).
check_outliers(Carbohydrate_DHW_Depth)
# OK: No outliers detected.
# - Based on the following method and threshold: cook (0.8).
# - For variable: (Whole model)
check_collinearity(Carbohydrate_DHW_Depth)
# Check for Multicollinearity
# Low Correlation
# Term VIF VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# DHW 1.02 [1.00, 101.18] 1.01 0.98 [0.01, 1.00]
# Depth 1.02 [1.00, 101.18] 1.01 0.98 [0.01, 1.00]
summary(Carbohydrate_DHW_Depth)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: log(Carbohydrate) ~ DHW + Depth + (1 | Station)
# Data: Pa_data
# 
# REML criterion at convergence: 161.1
# 
# Scaled residuals:
# Min 1Q Median 3Q Max
# -3.12624 -0.52693 -0.03515 0.54230 2.25947
# 
# Random effects:
# Groups Name Variance Std.Dev.
# Station (Intercept) 0.05339 0.2311
# Residual 0.27190 0.5214
# Number of obs: 89, groups: Station, 19
# 
# Fixed effects:
# Estimate Std. Error df t value Pr(>|t|)
# (Intercept) -1.472436 0.119431 41.902623 -12.329 1.58e-15 ***
# DHW -0.067730 0.014053 29.732817 -4.819 3.96e-05 ***
# Depth -0.007447 0.015222 85.749530 -0.489 0.626
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) DHW
# DHW -0.313
# Depth -0.591 -0.155
tab_model(Carbohydrate_DHW_Depth)
figCarDHW <- visreg(Carbohydrate_DHW_Depth, "DHW", ylab="log(Carbohydrate content)", line=list(col="brown"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 1E.pdf", figCarDHW, width = 5, height = 3, units = "in")

Lipid_DHW_Depth=lmer(sqrt(Lipid)~DHW+Depth+(1|Station), data=Pa_data)
check_model(Lipid_DHW_Depth)
check_normality(Lipid_DHW_Depth)
# OK: residuals appear as normally distributed (p = 0.299).
check_heteroscedasticity(Lipid_DHW_Depth)
# Warning: Heteroscedasticity (non-constant error variance) detected (p < .001).Warning message:
# In sqrt(insight::get_deviance(x)/(insight::n_obs(x) - sum(!is.na(estimates)))) :
# NaNs produced
check_outliers(Lipid_DHW_Depth)
# OK: No outliers detected.
# - Based on the following method and threshold: cook (0.8).
# - For variable: (Whole model)
check_collinearity(Lipid_DHW_Depth)
# Check for Multicollinearity
# Low Correlation
# Term VIF VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# DHW 1.04 [1.00, 6.53] 1.02 0.96 [0.15, 1.00]
# Depth 1.04 [1.00, 6.53] 1.02 0.96 [0.15, 1.00]
summary(Lipid_DHW_Depth)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: sqrt(Lipid) ~ DHW + Depth + (1 | Station)
# Data: Pa_data
# 
# REML criterion at convergence: -21.9
# 
# Scaled residuals:
# Min 1Q Median 3Q Max
# -2.10535 -0.42338 -0.05018 0.52034 2.32362
# 
# Random effects:
# Groups Name Variance Std.Dev.
# Station (Intercept) 0.02635 0.1623
# Residual 0.02683 0.1638
# Number of obs: 89, groups: Station, 19
# 
# Fixed effects:
# Estimate Std. Error df t value Pr(>|t|)
# (Intercept) 1.058018 0.051053 27.747721 20.724 < 2e-16 ***
# DHW -0.015802 0.005501 81.465836 -2.873 0.005188 **
# Depth -0.018358 0.005068 77.718031 -3.622 0.000519 ***
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) DHW
# DHW -0.175
# Depth -0.458 -0.205
tab_model(Lipid_DHW_Depth)
figLipDHW <- visreg(Lipid_DHW_Depth, "DHW", ylab="sqrt(Lipid content)", line=list(col="brown"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 1F.pdf", figLipDHW, width = 5, height = 3, units = "in")
