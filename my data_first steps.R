setwd("~/Tesis/R_repositorios/Multitrophic_plants")

library(tidyverse)
library(cxr)


#este script es para probar el codigo con los datos de 2020 de Caracoles. 
#Para este año en concreto debo recordar de que no hay abundancias en el plot 4
#para todas las especies

fv <- read_csv2("data/pollinator_2020.csv")
h <- read_csv2("data/foodweb_2020.csv")
comp <- read_csv2("data/competition_2020.csv")
salt <- read_csv2("data/Salinity_moisture_2020.csv")


#pasos a dar
#primero lo que tendré que hacer es coger los datos de competicion, asignarle un numero a cada fila que sea único obs_ID, 
#luego lo que tendré que hacer es separarlo por especies, y por último me tendré que quedar solo con seed y las diferentes
#species de plantas. 
comp$ID_obs = 1:dim(comp)[1] #creo el ID_obs, que debe ser 1 por coordenada espacial
comp.s <- comp[,c("ID_obs",  "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
        ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ","ANAR", "RAPE", "ACHI")]
comp.list <- split(comp.s, interaction(comp$Focal)) #lista de neighbors

#######DUDA: PARA QUE ES LA LISTA LLAMADA SPATIAL SAMPLING?; otra duda, en todas las listas el id_obs
#va relacionado con el lugar espacial que ocupa?respuesta: Sí, está relacionado con el lugar que ocupa. David lo
#mantiene, pero el paquete en si no lo necesita. 

#SALINIDAD: me tengo que quedar con Obs_ID y la salinity
salt$ID_obs = 1:dim(salt)[1] #creo el ID_obs
salt.s <- salt [,c("ID_obs", "Salinity")] # necesito luego tener una sublista por especie

#importante, los datos de salinidad los tengo que tener sumados, nos interesa saber
#el valor de todo el año
sal.simple <- salt %>% group_by(Plot, Subplot) %>% summarise(sal = sum(Salinity))


#tengo que cargar los datos de seeds survival y germination


#floral visitors. tengo que hacer el sumatorio de todas las visitas por plot y subplot
fv.simple <- fv %>% group_by(Plot, Subplot, Plant) %>% summarise (visits = sum(number_visits))%>%
    ungroup()
fv.simple$Focal <-fv.simple$Plant


####
a<-left_join(sal.simple, comp)
b<-a[,c("Plot","Subplot","Focal", "sal")]#sal con focal


c <- left_join(b,fv.simple)
c <-subset(c, select = -c(Plant))

c$visits[is.na(c$visits)]<- 0
data <- c
covdata.list <- split(data, interaction(data$Focal)) #asi son las covariables, pero
#ademas necesito tener mi lista de vecinos.

solosal <- data[,c("Focal", "sal")]
sal.list <- split(solosal, interaction(solosal$Focal))

###1: comencemos solo por 1 sp ----
#Para ello solo necesito los datos de competencia de una especie
#concreta, y no necesito el id_obs. A través de esto voy a sacar
#lambda y la matriz de competencia (alpha). Por eso solo
#necesito los vecinos y el fitness. Ahora lo voy a hacer
#solo con una sp, pero cuando quiera hacerlo con varias 
#tengo que usar la funcion cxr_pm_multifit.

#alla vamoos
my.sp <- "HOMA"
#ahora saco los datos de esta especie de la lista que cree
obs_homa <- comp.list[[my.sp]]
#como no necesito la id_obs la borro
obs_homa <- subset(obs_homa, select = -c(ID_obs))
head(obs_homa)

#corremos el modelo para sacar lambda y alpha

#?cxr_pm_fit #check the help file for a description of the arguments
fit_homa <- cxr_pm_fit(data = obs_homa,
                       focal_column = my.sp,
                       model_family = "RK",
                       covariates = NULL,
                       optimization_method = "Nelder-Mead",
                       alpha_form = "pairwise",
                       lambda_cov_form = "none",
                       alpha_cov_form = "none",
                       initial_values = list(lambda = 1,
                                             alpha_intra = .1,
                                             alpha_inter = .1),
                       #not aplicable to this optimazation method
                       # lower_bounds = list(lambda = 0, 
                       #                     alpha_intra = 0,
                       #                     alpha_inter = 0),
                       # upper_bounds = list(lambda = 10,
                       #                     alpha_intra = 1,
                       #                     alpha_inter = 1),
                       fixed_terms = NULL,
                       # a low number of bootstrap samples
                       # for demonstration purposes, 
                       # increase it for robust results.
                       bootstrap_samples = 3) 
summary(fit_homa)
names(fit_homa) #list of all available elements.
fit_homa$lambda
# intraspecific interaction
fit_homa$alpha_intra
# interspecific interactions
fit_homa$alpha_inter

#ahora voy a probar lo mismo pero con muchas especies, funcion cxr_pm_multifit. 
#primero lo que tengo que hacer es descomponer en sublistas la lista de salinidad. 
#cómo lo puedo hacer? junto las abundancias, species y plot y subplot

#Dato sobre la sal: voy a coger el sumatorio, recuerdalo. 

#mas de 1 sp ----

my.sp <- c("HOMA", "CETE", "LEMA")
obs_3sp <- comp.list[my.sp]
any(is.na(obs_3sp)) #compruebo que no hay ningun NA

# discard ID column
for(i in 1:length(obs_3sp)){
    obs_3sp[[i]] <- obs_3sp[[i]][,2:length(obs_3sp[[i]])]
}
# load covariates: salinity

salinity <- sal.list[my.sp]

# keep only salinity column           
for(i in 1:length(salinity)){
    salinity[[i]] <- as.matrix(salinity[[i]][,2:length(salinity[[i]])])
    colnames(salinity[[i]]) <- "salinity"
    
}
names(obs_3sp)
# observation data
head(obs_3sp[[1]])
# number of fitness observations
nrow(obs_3sp[[1]])
# salinity data
head(salinity[[1]])
# number of covariate observations
nrow(salinity[[1]])

fit_3sp <- cxr_pm_multifit(data = obs_3sp,
                           focal_column = my.sp,
                           model_family = "RK",
                           # here we use a bounded method for demonstration purposes
                           optimization_method = "bobyqa", 
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
                           bootstrap_samples = 0)
