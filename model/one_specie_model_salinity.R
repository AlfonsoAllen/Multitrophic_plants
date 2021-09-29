#In this script i'm going to do the model to obtain lambda, alphas and the cov for the nexts models
#in this script i am going to use the SALINITY
#############
#######                     SALINITY
##########

#install.packages("remotes")
#intalar actualizacion de david
#remotes::install_github("RadicalCommEcol/cxr")

library(tidyverse)
library(cxr)
library(ggplot2)
#load data
comp <- read_csv2("total_comp_19_20.check.csv")
cov <- read_csv2("covariates_salt_h_fv_20_19.csv")

#
nombres <- c("LEMA", "BEMA","HOMA", "SOAS","CHFU","SCLA","SPRU","POMA","POMO",
             "PAIN","CETE","PLCO","PUPA","MESU")
output <- list()#creo una lista donde voy a guardar todos los output necesarios del loop

for (i in 1:length(nombres)){
    #ka_1 [i] <- which(colnames(comp)==nombres [i])
    matriz.subset <- subset(comp, focal == nombres [i])
    #crear nombres para que no me sobreescriba la matriz
    matriz.subset$intra <- rowSums(matriz.subset[,c(nombres [i])])
    #nombres unicos para cada especie. no sobreescribir
    matriz.subset$inter <- rowSums(select(matriz.subset, -c("focal", "fitness","intra", nombres [i])))
    output[[i]]<- matriz.subset  #doble corche para las hojas de la lista. Así se me guarda cada especie por hoja
    
}

final.m <- do.call("rbind",output)#esto me sirve para pegar cada uno de los output y guardarlos
#final.m[which(final.m$focal== "BEMA"),]#esto me sirve para comprobar en la lista que la especie focal que he guardado
#                                       en la lista corresponde con la que se ha guardado en la lista general
data.comp <- final.m 
comp.list <- split(data.comp, interaction(data.comp$focal))

### ahora de esta lista es de la que voy a hacer el subset de las species: species primeras a probar 
#LEMA
lema <- comp.list$LEMA #ya está hecho el subset de las especies
cov <- cov[,c("Plant", "herb1","fv1","sal")]
cov.list <- split(cov, interaction(cov$Plant))
# ahora selecciono mis species
sp <- c("LEMA")
any(is.na(lema))
#checks
# discard focal colum
lema<- lema[,c("fitness","LEMA", "BEMA","HOMA", "SOAS","CHFU","SCLA","SPRU","POMA","POMO",
                  "PAIN","CETE","PLCO","PUPA","MESU","inter", "intra")]
lema.junto <- lema[,c("fitness","LEMA", "BEMA","HOMA", "SOAS","CHFU","SCLA","SPRU","POMA","POMO",
        "PAIN","CETE","PLCO","PUPA","MESU")]
lema.junto.list <- list(lema.junto)
lema.intrainter <- lema[,c("fitness","intra", "inter")]
lema.intrainter.list <- list(lema.intrainter)

lema.list <- list(lema)
# load covariates: salinity+herb+fv
cov1.lema <- cov.list[sp]
# discard the plant column          
for(i in 1:length(cov1.lema)){
    cov1.lema[[i]] <- cov1.lema[[i]][,2:4]
}
names(lema)
# observation data
head(lema[[1]])
# number of fitness observations
nrow(lema)
# salinity data
head(cov1.lema)
# number of covariate observations
nrow(cov1.lema[[1]])
#1.salinity ----
# 1.1 lema
# load covariates: salinity  
salinity <- cov.list[sp]
# discard the plant column  and herb and fv  
sal.lema <- cov.list$LEMA["sal"]
sal.lema <- list(sal.lema)

# discard intra e inter colum
mis_species_junto <- mis_species #hago una copia de mis_species que es el original (no sobreescribir)
lema.junto <- data.frame(lema.junto.list)
names(lema.junto.list[1])
# observation data
head(lema.junto.list[[1]])
# number of fitness observations
nrow(lema.junto.list[[1]])
# salinity data
head(sal.lema[[1]])
# number of salinity observations
nrow(sal.lema[[1]])
any(is.na(lema.junto.list))
sal.lema <- data.frame(sal.lema)

#modelo
modelo.sal.lema.junto <- cxr_pm_fit(data = lema.junto,
                              focal_column = sp,
                              model_family = "BH",
                              # here we use a bounded method for demonstration purposes
                              optimization_method = "GenSA", 
                              covariates = sal.lema,
                              alpha_form = "pairwise",
                              lambda_cov_form = "global", # effect of covariates over lambda
                              alpha_cov_form = "global", # effect of covariates over alpha
                              initial_values = list(lambda = 1,
                                                    alpha_intra = 0.1,
                                                    alpha_inter = 0.1,
                                                    lambda_cov = 0.1,
                                                    alpha_cov = 0.1),
                              lower_bounds = list(lambda = 0,
                                                  alpha_intra = -10,
                                                  alpha_inter = -10,
                                                  lambda_cov = -10,
                                                  alpha_cov = -10),
                              upper_bounds = list(lambda = 1500,
                                                  alpha_intra = 10,
                                                  alpha_inter = 10,
                                                  lambda_cov = 10,
                                                  alpha_cov = 10),
                              # no standard errors
                              bootstrap_samples = 0) #salen problemas con los boundaries

sp.intra.lema <- c("intra")
modelo.sal.intrainter <- cxr_pm_fit(data = lema.intrainter,
                                    focal_column = sp.intra.lema,
                                    model_family = "BH",
                                    # here we use a bounded method for demonstration purposes
                                    optimization_method = "GenSA", 
                                    covariates = sal.lema,
                                    alpha_form = "pairwise",
                                    lambda_cov_form = "global", # effect of covariates over lambda
                                    alpha_cov_form = "global", # effect of covariates over alpha
                                    initial_values = list(lambda = 1,
                                                          alpha_intra = 0.1,
                                                          alpha_inter = 0.1,
                                                          lambda_cov = 0.1,
                                                          alpha_cov = 0.1),
                                    lower_bounds = list(lambda = 0,
                                                        alpha_intra = -15,
                                                        alpha_inter = -10,
                                                        lambda_cov = -10,
                                                        alpha_cov = -10),
                                    upper_bounds = list(lambda = 1500,
                                                        alpha_intra = 15,
                                                        alpha_inter = 10,
                                                        lambda_cov = 10,
                                                        alpha_cov = 10),
                                    # no standard errors
                                    bootstrap_samples = 0) #salen problemas con los boundaries
                #no fitea alpha intra, no da el valor

summary(modelo.sal.intrainter)

#chfu
chfu <- comp.list$CHFU #ya está hecho el subset de la especies
sp1 <- c("CHFU")
any(is.na(chfu))
#checks
# discard focal colum
chfu<- chfu[,c("fitness","LEMA", "BEMA","HOMA", "SOAS","CHFU","SCLA","SPRU","POMA","POMO",
               "PAIN","CETE","PLCO","PUPA","MESU","inter", "intra")]
chfu.junto <- chfu[,c("fitness","LEMA", "BEMA","HOMA", "SOAS","CHFU","SCLA","SPRU","POMA","POMO",
                      "PAIN","CETE","PLCO","PUPA","MESU")]
chfu.intrainter <- chfu[,c("fitness","intra", "inter")]

# load covariates: salinity+herb+fv
cov1.chfu <- cov.list[sp1]
# discard the plant column  and herb and fv  
sal.chfu <- cov1.chfu$CHFU["sal"]

names(chfu.junto)
# observation data
head(chfu.junto[[1]])
# number of fitness observations
nrow(chfu.junto[1])
# salinity data
head(sal.chfu[[1]])
# number of salinity observations
nrow(sal.chfu[1])
any(is.na(chfu.junto))

#modelo
modelo.sal.chfu.junto <- cxr_pm_fit(data = chfu.junto,
                         focal_column = sp1,
                         model_family = "BH",
                         # here we use a bounded method for demonstration purposes
                         optimization_method = "GenSA", 
                         covariates = sal.chfu,
                         alpha_form = "pairwise",
                         lambda_cov_form = "global", # effect of covariates over lambda
                         alpha_cov_form = "global", # effect of covariates over alpha
                         initial_values = list(lambda = 1,
                                               alpha_intra = 0.1,
                                               alpha_inter = 0.1,
                                               lambda_cov = 0.1,
                                               alpha_cov = 0.1),
                         lower_bounds = list(lambda = 0,
                                             alpha_intra = -10,
                                             alpha_inter = -10,
                                             lambda_cov = -50,
                                             alpha_cov = -10),
                         upper_bounds = list(lambda = 1500,
                                             alpha_intra = 10,
                                             alpha_inter = 15,
                                             lambda_cov = 50,
                                             alpha_cov = 10),
                         # no standard errors
                         bootstrap_samples = 0) #salen problemas con los boundaries

spi.intra <- c("intra")
modelo.sal.chfu.intrainter <- cxr_pm_fit(data = chfu.intrainter,
                                    focal_column = spi.intra,
                                    model_family = "BH",
                                    # here we use a bounded method for demonstration purposes
                                    optimization_method = "GenSA", 
                                    covariates = sal.chfu,
                                    alpha_form = "pairwise",
                                    lambda_cov_form = "global", # effect of covariates over lambda
                                    alpha_cov_form = "global", # effect of covariates over alpha
                                    initial_values = list(lambda = 1,
                                                          alpha_intra = 0.1,
                                                          alpha_inter = 0.1,
                                                          lambda_cov = 0.1,
                                                          alpha_cov = 0.1),
                                    lower_bounds = list(lambda = 0,
                                                        alpha_intra = -10,
                                                        alpha_inter = -10,
                                                        lambda_cov = -50,
                                                        alpha_cov = -10),
                                    upper_bounds = list(lambda = 1500,
                                                        alpha_intra = 10,
                                                        alpha_inter = 15,
                                                        lambda_cov = 50,
                                                        alpha_cov = 10),
                                    # no standard errors
                                    bootstrap_samples = 0) #salen problemas con los boundaries
summary(modelo.sal.chfu.intrainter)

#PUPA
pupa <- comp.list$PUPA #ya está hecho el subset de la especies
sp1.pupa <- c("PUPA")
any(is.na(pupa))
#checks
# discard focal colum
pupa<- pupa[,c("fitness","LEMA", "BEMA","HOMA", "SOAS","CHFU","SCLA","SPRU","POMA","POMO",
               "PAIN","CETE","PLCO","PUPA","MESU","inter", "intra")]
pupa.junto <- pupa[,c("fitness","LEMA", "BEMA","HOMA", "SOAS","CHFU","SCLA","SPRU","POMA","POMO",
                      "PAIN","CETE","PLCO","PUPA","MESU")]
pupa.intrainter <- pupa[,c("fitness","intra", "inter")]

# load covariates: salinity+herb+fv
cov1.pupa <- cov.list[sp1.pupa]
# discard the plant column  and herb and fv  
sal.pupa <- cov1.pupa$PUPA["sal"]

names(pupa.junto)
# observation data
head(pupa.junto)
# number of fitness observations
nrow(pupa.junto)
# salinity data
head(sal.pupa[[1]])
# number of salinity observations
nrow(sal.pupa[1])
any(is.na(chfu.junto))

modelo.sal.pupa.junto <- cxr_pm_fit(data = pupa.junto,
                                    focal_column = sp1.pupa,
                                    model_family = "BH",
                                    # here we use a bounded method for demonstration purposes
                                    optimization_method = "GenSA", 
                                    covariates = sal.pupa,
                                    alpha_form = "pairwise",
                                    lambda_cov_form = "global", # effect of covariates over lambda
                                    alpha_cov_form = "global", # effect of covariates over alpha
                                    initial_values = list(lambda = 1,
                                                          alpha_intra = 0.1,
                                                          alpha_inter = 0.1,
                                                          lambda_cov = 0.1,
                                                          alpha_cov = 0.1),
                                    lower_bounds = list(lambda = 0,
                                                        alpha_intra = -10,
                                                        alpha_inter = -10,
                                                        lambda_cov = -50,
                                                        alpha_cov = -10),
                                    upper_bounds = list(lambda = 1500,
                                                        alpha_intra = 10,
                                                        alpha_inter = 15,
                                                        lambda_cov = 50,
                                                        alpha_cov = 10),
                                    # no standard errors
                                    bootstrap_samples = 0) #salen problemas con los boundaries
sp1.pupa.intra <- c("intra")
modelo.sal.pupa.intrainter <- cxr_pm_fit(data = pupa.intrainter,
                                    focal_column = sp1.pupa.intra,
                                    model_family = "BH",
                                    # here we use a bounded method for demonstration purposes
                                    optimization_method = "GenSA", 
                                    covariates = sal.pupa,
                                    alpha_form = "pairwise",
                                    lambda_cov_form = "global", # effect of covariates over lambda
                                    alpha_cov_form = "global", # effect of covariates over alpha
                                    initial_values = list(lambda = 1,
                                                          alpha_intra = 0.1,
                                                          alpha_inter = 0.1,
                                                          lambda_cov = 0.1,
                                                          alpha_cov = 0.1),
                                    lower_bounds = list(lambda = 0,
                                                        alpha_intra = -10,
                                                        alpha_inter = -10,
                                                        lambda_cov = -50,
                                                        alpha_cov = -10),
                                    upper_bounds = list(lambda = 1500,
                                                        alpha_intra = 10,
                                                        alpha_inter = 15,
                                                        lambda_cov = 50,
                                                        alpha_cov = 10),
                                    # no standard errors
                                    bootstrap_samples = 0) #salen problemas con los boundaries

#3 species a la vez
my.sp <- c("LEMA","CHFU","PUPA")
obs_3sp <- comp.list[my.sp]
# discard ID column
for(i in 1:length(obs_3sp)){
    obs_3sp[[i]] <- obs_3sp[[i]][,2:length(obs_3sp[[i]])]
}
# load covariates: salinity
head(cov.list)
salinity <- cov.list[my.sp]
# keep only salinity column
for(i in 1:length(salinity)){
    salinity[[i]] <- as.matrix(salinity[[i]][,4])
    colnames(salinity[[i]]) <- "salinity"
}

modelo.fit_3sp <- cxr_pm_multifit(data = obs_3sp,
                           focal_column = my.sp,
                           model_family = "BH",
                           # here we use a bounded method for demonstration purposes
                           optimization_method = "GenSA", 
                           covariates = salinity,
                           alpha_form = "pairwise",
                           lambda_cov_form = "global", # effect of covariates over lambda
                           alpha_cov_form = "global", # effect of covariates over alpha
                           initial_values = list(lambda = 1,
                                                 alpha_intra = 0.1,
                                                 alpha_inter = 0.1,
                                                 lambda_cov = 0.1,
                                                 alpha_cov = 0.1),
                           lower_bounds = list(lambda = 0,
                                               alpha_intra = -10,
                                               alpha_inter = -10,
                                               lambda_cov = -10,
                                               alpha_cov = -10),
                           upper_bounds = list(lambda = 1500,
                                               alpha_intra = 10,
                                               alpha_inter = 10,
                                               lambda_cov = 10,
                                               alpha_cov = 10),
                           # no standard errors
                           bootstrap_samples = 0) #sale error en alpha matrix
