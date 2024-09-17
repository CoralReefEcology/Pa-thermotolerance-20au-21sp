Pa_data = read.csv("Pa-Physical-XS-20au-21sp-total.csv", header = TRUE, row.names = 1)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(lmerTest)
library(performance)
library(ggpmisc)
library(visreg)
library(sjPlot)

##Effects of algal symbiont density, D1_ratio and coral C13 on the nutrient contents
Protein_Density_D1_C13=lmer(Protein~Density+D1_ratio+C13+(1|Station), data=Pa_data)
check_model(Protein_Density_D1_C13)
check_normality(Protein_Density_D1_C13)
# OK: residuals appear as normally distributed (p = 0.204).
check_heteroskedasticity(Protein_Density_D1_C13)
# Warning: Heteroscedasticity (non-constant error variance) detected (p < .001).Warning message:
# In sqrt(insight::get_deviance(x)/(insight::n_obs(x) - sum(!is.na(estimates)))) :
# NaNs produced
check_outliers(Protein_Density_D1_C13)
# OK: No outliers detected.
# - Based on the following method and threshold: cook (0.846).
# - For variable: (Whole model)
check_collinearity(Protein_Density_D1_C13)
# Check for Multicollinearity
# Low Correlation
# Term VIF VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# Density 1.00 [1.00, Inf] 1.00 1.00 [0.00, 1.00]
# D1_ratio 1.00 [1.00, Inf] 1.00 1.00 [0.00, 1.00]
# C13 1.00 [1.00, Inf] 1.00 1.00 [0.00, 1.00]
summary(Protein_Density_D1_C13)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Protein ~ Density + D1_ratio + C13 + (1 | Station)
# Data: Pa_data
# 
# REML criterion at convergence: -70.4
# 
# Scaled residuals:
# Min 1Q Median 3Q Max
# -2.8292 -0.5802 -0.1173 0.5046 2.4380
# 
# Random effects:
# Groups Name Variance Std.Dev.
# Station (Intercept) 0.008954 0.09463
# Residual 0.016525 0.12855
# Number of obs: 89, groups: Station, 19
# 
# Fixed effects:
# Estimate Std. Error df t value Pr(>|t|)
# (Intercept) -0.181432 0.203850 84.999732 -0.890 0.37596
# Density 0.014281 0.006271 83.797164 2.277 0.02533 *
# D1_ratio 0.033367 0.049011 84.344432 0.681 0.49787
# C13 -0.032061 0.010146 84.098152 -3.160 0.00219 **
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) Densty D1_rat
# Density -0.097
# D1_ratio -0.184 -0.024
# C13 0.966 0.049 -0.019
tab_model(Protein_Density_D1_C13)
figProDensity <- visreg(Protein_Density_D1_C13, "Density", ylab="Protein content", line=list(col="brown"), fill=list(fill="green"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 2D.pdf", figProDensity, width = 5, height = 3, units = "in")
figProC13 <- visreg(Protein_Density_D1_C13, "C13", ylab="Protein content", line=list(col="yellow"), fill=list(fill="lightblue"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 4E.pdf", figProC13, width = 5, height = 3, units = "in")


Carbohydrate_Density_D1_C13=lmer(log(Carbohydrate)~Density+D1_ratio+C13+(1|Station), data=Pa_data)
check_model(Carbohydrate_Density_D1_C13)
check_normality(Carbohydrate_Density_D1_C13)
# Warning: Non-normality of residuals detected (p = 0.014).
check_heteroscedasticity(Carbohydrate_Density_D1_C13)
# OK: Error variance appears to be homoscedastic (p = 0.966).
check_outliers(Carbohydrate_Density_D1_C13)
# OK: No outliers detected.
# - Based on the following method and threshold: cook (0.846).
# - For variable: (Whole model)
check_collinearity(Carbohydrate_Density_D1_C13)
# Check for Multicollinearity
# Low Correlation
# Term VIF VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# Density 1.01 [1.00, 2.65e+14] 1.00 0.99 [0.00, 1.00]
# D1_ratio 1.00 [1.00, Inf] 1.00 1.00 [0.00, 1.00]
# C13 1.00 [1.00, Inf] 1.00 1.00 [0.00, 1.00]
summary(Carbohydrate_Density_D1_C13)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: log(Carbohydrate) ~ Density + D1_ratio + C13 + (1 | Station)
# Data: Pa_data
# 
# REML criterion at convergence: 172.4
# 
# Scaled residuals:
# Min 1Q Median 3Q Max
# -3.2969 -0.5263 -0.1248 0.6021 2.3301
# 
# Random effects:
# Groups Name Variance Std.Dev.
# Station (Intercept) 0.07598 0.2756
# Residual 0.31436 0.5607
# Number of obs: 89, groups: Station, 19
# 
# Fixed effects:
# Estimate Std. Error df t value Pr(>|t|)
# (Intercept) -3.673604 0.839935 76.189513 -4.374 3.82e-05 ***
# Density -0.009322 0.026271 84.792146 -0.355 0.7236
# D1_ratio 0.434082 0.205652 84.828438 2.111 0.0377 *
# C13 -0.086186 0.042287 82.584210 -2.038 0.0447 *
# ---
# Signif. codes: 0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
# (Intr) Densty D1_rat
# Density -0.088
# D1_ratio -0.171 -0.045
# C13 0.968 0.055 -0.004
tab_model(Carbohydrate_Density_D1_C13)
figCarD1 <- visreg(Carbohydrate_Density_D1_C13, "D1_ratio", ylab="log(Carbohydrate content)", line=list(col="brown"), fill=list(fill="green"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 2E.pdf", figCarD1, width = 5, height = 3, units = "in")
figCarC13 <- visreg(Carbohydrate_Density_D1_C13, "C13", ylab="log(Carbohydrate content)", line=list(col="yellow"), fill=list(fill="lightblue"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig. 4F.pdf", figCarC13, width = 5, height = 3, units = "in")


Lipid_Density_D1_C13=lmer(Lipid~Density+D1_ratio+C13+(1|Station), data=Pa_data)
check_model(Lipid_Density_D1_C13)
check_normality(Lipid_Density_D1_C13)
# OK: residuals appear as normally distributed (p = 0.083).
check_heteroscedasticity(Lipid_Density_D1_C13)
# OK: Error variance appears to be homoscedastic (p = 0.842).
check_outliers(Lipid_Density_D1_C13)
# OK: No outliers detected.
# - Based on the following method and threshold: cook (0.846).
# - For variable: (Whole model)
check_collinearity(Lipid_Density_D1_C13)
# Check for Multicollinearity
# Low Correlation
# Term VIF VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# Density 1.00 [1.00, Inf] 1.00 1.00 [0.00, 1.00]
# D1_ratio 1.00 [1.00, Inf] 1.00 1.00 [0.00, 1.00]
# C13 1.00 [1.00, Inf] 1.00 1.00 [0.00, 1.00]
summary(Lipid_Density_D1_C13)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: Lipid ~ Density + D1_ratio + C13 + (1 | Station)
# Data: Pa_data
# 
# REML criterion at convergence: 77.5
# 
# Scaled residuals:
# Min 1Q Median 3Q Max
# -2.51538 -0.43626 -0.04681 0.47427 2.57403
# 
# Random effects:
# Groups Name Variance Std.Dev.
# Station (Intercept) 0.09794 0.3129
# Residual 0.08425 0.2903
# Number of obs: 89, groups: Station, 19
# 
# Fixed effects:
# Estimate Std. Error df t value Pr(>|t|)
# (Intercept) 0.335473 0.481113 82.408006 0.697 0.488
# Density 0.009926 0.014635 79.440746 0.678 0.500
# D1_ratio -0.065197 0.115216 82.313128 -0.566 0.573
# C13 -0.029274 0.023664 78.694279 -1.237 0.220
# 
# Correlation of Fixed Effects:
# (Intr) Densty D1_rat
# Density -0.101
# D1_ratio -0.199 -0.007
# C13 0.960 0.047 -0.035
tab_model(Lipid_Density_D1_C13)


##Effects of D1 proportion on Ec in the coral host
Ec_Density_D1=lmer(sqrt(Ec)~Density + D1_ratio+(1|Station), data=Pa_data)
check_model(Ec_Density_D1)
check_normality(Ec_Density_D1)
#OK: residuals appear as normally distributed (p = 0.155).
check_heteroskedasticity(Ec_Density_D1)
# OK: Error variance appears to be homoscedastic (p = 0.523).
check_outliers(Ec_Density_D1)
# OK: No outliers detected.
# - Based on the following method and threshold: cook (0.8).
# - For variable: (Whole model)
check_collinearity(Ec_Density_D1)
# Low Correlation
# 
# Term  VIF  VIF 95% CI Increased SE Tolerance Tolerance 95% CI
# Density 1.00 [1.00, Inf]         1.00      1.00     [0.00, 1.00]
# D1_ratio 1.00 [1.00, Inf]         1.00      1.00     [0.00, 1.00]
summary(Ec_D1)
# Linear mixed model fit by REML. t-tests use Satterthwaite's method ['lmerModLmerTest']
# Formula: sqrt(Ec) ~ Density + D1_ratio + (1 | Station)
#    Data: Pa_data
# 
# REML criterion at convergence: 443.2
# 
# Scaled residuals: 
#      Min       1Q   Median       3Q      Max 
# -2.61289 -0.46026 -0.08076  0.44575  1.96366 
# 
# Random effects:
#  Groups   Name        Variance Std.Dev.
#  Station  (Intercept) 6.944    2.635   
#  Residual             6.760    2.600   
# Number of obs: 88, groups:  Station, 19
# 
# Fixed effects:
#             Estimate Std. Error      df t value Pr(>|t|)    
# (Intercept)  18.0920     1.2011 67.7075  15.063   <2e-16 ***
# Density       0.3166     0.1303 81.0661   2.429   0.0173 *  
# D1_ratio      1.8147     1.0559 83.8351   1.719   0.0894 .  
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Correlation of Fixed Effects:
#          (Intr) Densty
# Density  -0.521       
# D1_ratio -0.621 -0.005
tab_model(Ec_Density_D1)
figEcDensity <- visreg(Ec_Density_D1, "Density", ylab="sqrt(Ec)", line=list(col="brown"), fill=list(fill="yellow"), points=list(cex=12, pch=19, size=3), gg=TRUE)+theme(axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.text = element_text(size = 10))+scale_color_npg()+scale_fill_npg()
ggsave("Fig.2F.pdf", figEcDensity, width = 5, height = 3, units = "in")
