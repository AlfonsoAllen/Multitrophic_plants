setwd("~/Tesis/R_repositorios/Multitrophic_plants")

library(tidyverse)
library(cxr)
library(ggplot2)
#load data
comp <- read_csv2("total_comp_19_20.csv")
cov <- read_csv2("covariates_salt_h_fv_20_19.csv")

####################                    coexistance metrics             #########################
#first: list of competition and covariates
comp <- comp[,c("focal","fitness","BEMA", "CHFU", "CHMI", "CETE", "HOMA","LEMA","PAIN","PLCO","POMA","POMO",
        "PUPA","SASO","SCLA","SPRU","MESU","MEEL","MEPO","SOAS","FRPU","SUSP","COSQ","RAPE")]
comp.list <- split(comp, interaction(comp$focal))
cov <- cov[,c("Plant", "herb1","fv1","sal")]
cov.list <- split(cov, interaction(cov$Plant))
##
#select the species (again)
my.sp2 <- c("HOMA", "POMO","POMA", "LEMA","CHFU", "PUPA", "SOAS", "SCLA", "MESU", "SPRU","BEMA",
           "CETE", "PAIN", "PLCO")

obs_3sp <- comp.list[my.sp2]
any(is.na(obs_3sp))#check NAs

# discard focal colum
for(i in 1:length(obs_3sp)){
    obs_3sp[[i]] <- obs_3sp[[i]][,2:length(obs_3sp[[i]])]
}
# load covariates: salinity+herb+fv
cov1 <- cov.list[my.sp2]
# discard the plant column          
for(i in 1:length(cov1)){
    cov1[[i]] <- cov1[[i]][,2:4]
}
names(obs_3sp)
# observation data
head(obs_3sp[[1]])
# number of fitness observations
nrow(obs_3sp[[1]])
# salinity data
head(cov1[[1]])
# number of covariate observations
nrow(cov1[[1]])

#function parameters
model_family <- "BH"
covariates <- cov1
# bobyqa is generally more robust than other bounded methods
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "global"
alpha_cov_form <- "pairwise"
# note how lambda_cov and alpha_cov
# have different initial values for each covariate effect
# the commented assignations are also possible, 
# giving equal initial values to all parameters
initial_values = list(lambda = 1,
                      alpha_intra = 0.1,
                      alpha_inter = 0.1,
                      lambda_cov = c(0.1,0.2,0.1),
                      alpha_cov = c(0.1,0.2,0.1))


# same with boundaries
lower_bounds = list(lambda = 0,
                    alpha_intra = -1,
                    alpha_inter = -1,
                    lambda_cov = c(-1,0,-1),
                    alpha_cov = c(-1,0,-1))

upper_bounds = list(lambda = 100,
                    alpha_intra = 1,
                    alpha_inter = 1,
                    lambda_cov = c(1,2,1),
                    alpha_cov = c(1,2,1))


fixed_terms <- NULL
bootstrap_samples <- 3


fit_multi_cov <- cxr_pm_multifit(data = obs_3sp,
                                 focal_column = my.sp2,
                                 model_family = model_family,
                                 covariates = covariates,
                                 optimization_method = optimization_method,
                                 alpha_form = alpha_form,
                                 lambda_cov_form = lambda_cov_form,
                                 alpha_cov_form = alpha_cov_form,
                                 initial_values = initial_values,
                                 lower_bounds = lower_bounds,
                                 upper_bounds = upper_bounds,
                                 fixed_terms = fixed_terms,
                                 bootstrap_samples = bootstrap_samples) #error: function returned is infinite or NA


fit_3sp <- cxr_pm_multifit(data = obs_3sp,
                           focal_column = my.sp2,
                           model_family = "RK",
                           # here we use a bounded method for demonstration purposes
                           optimization_method = "bobyqa", 
                           covariates = cov1,
                           alpha_form = "pairwise",
                           lambda_cov_form = "global", # effect of covariates over lambda
                           alpha_cov_form = "pairwise", # effect of covariates over alpha
                           initial_values = list(lambda = 1,
                                                 alpha_intra = 0.1,
                                                 alpha_inter = 0.1,
                                                 lambda_cov = 0.1,
                                                 alpha_cov = 0.1),
                           lower_bounds = list(lambda = 0,
                                               alpha_intra = 0,
                                               alpha_inter = -1,
                                               lambda_cov = 0,
                                               alpha_cov = 0),
                           upper_bounds = list(lambda = 100,
                                               alpha_intra = 1,
                                               alpha_inter = 1,
                                               lambda_cov = 1,
                                               alpha_cov = 1),
                           # no standard errors
                           bootstrap_samples = 0) #dan errores que no entiendo. Pero corre

summary(fit_3sp)
fit_3sp$log_likelihood #son valores negativos lo qe devuelve aunque aparezca
#en positivo
