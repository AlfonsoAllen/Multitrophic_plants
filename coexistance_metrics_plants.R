setwd("~/Tesis/R_repositorios/Multitrophic_plants")

library(tidyverse)
library(cxr)
library(ggplot2)
#load data
comp <- read_csv2("total_comp_19_20.check.csv")
cov <- read_csv2("covariates_salt_h_fv_20_19.csv")

#################data distribution 
hist(cov$herb1)# I have a lot of 0 in my data.
hist(log(cov$herb1)) 
hist(cov$fv1) #muchos ceros
hist(log(cov$fv1))
hist(cov$sal) #parece que mas o menos puede cumplir una distribucion normal sin hacer el log


####################                    coexistance metrics             #########################
#first: list of competition and covariates
comp<- comp[,c("focal","fitness","BEMA", "CHFU", "CHMI", "CETE", "HOMA","LEMA","PAIN","PLCO","POMA","POMO",
        "PUPA","SASO","SCLA","SPRU","MESU","MEEL","MEPO","SOAS","FRPU","SUSP","COSQ","RAPE")]

########## abundancias de datos ----
#adv <- subset(comp, focal =="LEMA")#en este subset he ido cambiando las especies para calcular cuantas focales tengo
#de cada
#adv1<- adv[,c("BEMA", "CHFU", "CHMI", "CETE", "HOMA","LEMA","PAIN","PLCO","POMA","POMO",
  #     "PUPA","SASO","SCLA","SPRU","MESU","MEEL","MEPO","SOAS","FRPU","SUSP","COSQ","RAPE")]
#VARIOS CHECKS
#vb <-rowSums(adv1)
#min(vb)
#nrow(adv)
###

prueba1 <- gather(comp,"BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
       ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE",
       key = "Plant", value ="abundance" )

prueba1$ID <- paste(prueba1$focal, prueba1$Plant, sep="_")
numero.combinacion.vecinos <-  prueba1 %>% group_by(ID) %>% summarise(total = sum(abundance))%>%
    ungroup()
numero.combinacion.vecinos$ID <- as.factor(numero.combinacion.vecinos$ID)
bb.mean <- numero.combinacion.vecinos %>% group_by(ID) %>% summarise(mean_neigh = mean(total))%>%
    ungroup()
bb.sd <- numero.combinacion.vecinos %>% group_by(ID) %>% summarise(sd_neigh = sd(total))%>%
    ungroup()
combi.neigh <- left_join(bb.mean, bb.sd) #base de datos con la media y desviacion estandar de las combinaciones
#       que ocurren en el campo. Quien es vecino de quien y la abundancia


#numero max y min de cada una de las covariables
max(cov$herb1)
min(cov$herb1)
max(cov$fv1)
min(cov$fv1)
max(cov$sal)
min(cov$sal)

######separar vecinos intra de los inter. Una columna con los intra y el resto los inter juntos y en otra los inter separados
head(comp)
#esto es para una especie para pensar el loop
#l.comp.mesu <- subset(comp, focal == "MESU")
#l.comp.mesu$n.intra <- l.comp.mesu[,c("MESU")]
#l.comp.mesu$n.inter <- rowSums(select(l.comp.mesu, -c("focal", "fitness", "n.intra", "MESU")))

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
###


################################### competition
comp.list <- split(data.comp, interaction(data.comp$focal))
cov <- cov[,c("Plant", "herb1","fv1","sal")]
cov.list <- split(cov, interaction(cov$Plant))
##
#select the species (again)
my.sp2 <- c("HOMA", "POMO","POMA", "LEMA","CHFU", "PUPA", "SOAS", "SCLA", "MESU", "SPRU","BEMA",
           "CETE", "PAIN", "PLCO")

mis_species<- comp.list[my.sp2]
any(is.na(mis_species))#check NAs

# discard focal colum
for(i in 1:length(my.sp2)){
    mis_species[[i]] <- mis_species[[i]][,2:length(mis_species[[i]])]
}
# load covariates: salinity+herb+fv
cov1 <- cov.list[my.sp2]
# discard the plant column          
for(i in 1:length(cov1)){
    cov1[[i]] <- cov1[[i]][,2:4]
}
names(mis_species)
# observation data
head(mis_species[[1]])
# number of fitness observations
nrow(mis_species[[4]])
# salinity data
head(cov1[[1]])
# number of covariate observations
nrow(cov1[[4]])

#       each separate covariable per model----
#1.salinity ----
# load covariates: salinity  
salinity <- cov.list[my.sp2]
# discard the plant column  and herb and fv         
for(i in 1:length(my.sp2)){
    salinity[[i]] <- salinity[[i]][,4] #he sobreescrito aqui, OJO!
}
# discard intra e inter colum
mis_species_junto <- mis_species #hago una copia de mis_species que es el original (no sobreescribir)
for(i in 1:length(nombres)){
  mis_species_junto[[i]]<- mis_species_junto[[i]][,1:23]
}
head(mis_species_junto[[1]])
# number of fitness observations
nrow(mis_species_junto[[2]])
# salinity data
head(salinity[[1]])
# number of salinity observations
nrow(salinity[[2]])
any(is.na(mis_species_junto))

modelo.sal <- cxr_pm_multifit(data = mis_species_junto,
                           focal_column = my.sp2,
                           model_family = "BH",
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
                                               alpha_intra = -1,
                                               alpha_inter = -1,
                                               lambda_cov = -1,
                                               alpha_cov = -1),
                           upper_bounds = list(lambda = 100,
                                               alpha_intra = 1,
                                               alpha_inter = 1,
                                               lambda_cov = 1,
                                               alpha_cov = 1),
                           # no standard errors
                           bootstrap_samples = 0) #salen problemas con los boundaries


summary(modelo.sal) #modelo con todas las especies seraradas (no inter+intra, ni intra+rest of species)
modelo.sal$log_likelihood 
#sal con intra e inter

# select intra e inter
mis_species_intra_inter <- mis_species #copia de mis_species. Asi no sobreescribo
for(i in 1:length(my.sp2)){
  mis_species_intra_inter[[i]] <- mis_species[[i]][,c(1,24,25 )]
}
head(mis_species_intra_inter[[1]])
# number of fitness observations
nrow(mis_species_intra_inter[[1]])
# salinity data
head(salinity[[1]])
# number of salinity observations
nrow(salinity[[1]])
any(is.na(mis_species_intra_inter))

modelo.sal.inter.intra <- cxr_pm_multifit(data = mis_species_intra_inter,
                              focal_column = my.sp2,
                              model_family = "BH",
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
                                                  alpha_intra = -1,
                                                  alpha_inter = -1,
                                                  lambda_cov = -1,
                                                  alpha_cov = -1),
                              upper_bounds = list(lambda = 100,
                                                  alpha_intra = 1,
                                                  alpha_inter = 1,
                                                  lambda_cov = 1,
                                                  alpha_cov = 1),
                              # no standard errors
                              bootstrap_samples = 0) 


summary(modelo.sal.inter.intra) #modelo con Inter +intra
#saveRDS(modelo.sal.inter.intra, file = "modelo.sal.inter.intra.rds") #esto es para guardar los outputs
                                          #de los modelos
modelo.sal.inter.intra$log_likelihood #estos datos voy a tener que meterlos en una tabla para ver cual es el mejor modelo

#2.herb ----
#load covariates: herb
#herb con todas las especies por separado
herb <- cov.list[my.sp2]
# discard the plant column          
for(i in 1:length(herb)){
    herb[[i]] <- herb[[i]][,2]
}
 
head(mis_species_junto[[1]]) #para los modelos los datos que van en el argument de data son los mismos si considero
#                           todas las especies juntas y luego cunaod consideros las de intra+inter
# number of fitness observations
nrow(mis_species_junto[[1]])
# salinity data
head(herb[[1]])
# number of salinity observations
nrow(herb[[1]])
modelo.herb.junto <- cxr_pm_multifit(data = mis_species_junto,
                              focal_column = my.sp2,
                              model_family = "BH",
                              # here we use a bounded method for demonstration purposes
                              optimization_method = "bobyqa", 
                              covariates = herb,
                              alpha_form = "pairwise",
                              lambda_cov_form = "global", # effect of covariates over lambda
                              alpha_cov_form = "global", # effect of covariates over alpha
                              initial_values = list(lambda = 1,
                                                    alpha_intra = 0.1,
                                                    alpha_inter = 0.1,
                                                    lambda_cov = 0.1,
                                                    alpha_cov = 0.1),
                              lower_bounds = list(lambda = 0,
                                                  alpha_intra = -1,
                                                  alpha_inter = -1,
                                                  lambda_cov = -1,
                                                  alpha_cov = -1),
                              upper_bounds = list(lambda = 100,
                                                  alpha_intra = 1,
                                                  alpha_inter = 1,
                                                  lambda_cov = 1,
                                                  alpha_cov = 1),
                              # no standard errors
                              bootstrap_samples = 0) #salen problemas con los boundaries

summary(modelo.herb.junto)
modelo.herb.junto$log_likelihood 
log_likelihood_herb <- enframe(modelo.herb.junto$log_likelihood)
log_likelihood_herb$log_likelihood_spjuntas <- log_likelihood_herb$value

# herb con vecinos intra + inter
modelo.herb.intra.inter <- cxr_pm_multifit(data = mis_species_intra_inter,
                                     focal_column = my.sp2,
                                     model_family = "BH",
                                     # here we use a bounded method for demonstration purposes
                                     optimization_method = "bobyqa", 
                                     covariates = herb,
                                     alpha_form = "pairwise",
                                     lambda_cov_form = "global", # effect of covariates over lambda
                                     alpha_cov_form = "global", # effect of covariates over alpha
                                     initial_values = list(lambda = 1,
                                                           alpha_intra = 0.1,
                                                           alpha_inter = 0.1,
                                                           lambda_cov = 0.1,
                                                           alpha_cov = 0.1),
                                     lower_bounds = list(lambda = 0,
                                                         alpha_intra = -1,
                                                         alpha_inter = -1,
                                                         lambda_cov = -1,
                                                         alpha_cov = -1),
                                     upper_bounds = list(lambda = 100,
                                                         alpha_intra = 1,
                                                         alpha_inter = 1,
                                                         lambda_cov = 1,
                                                         alpha_cov = 1),
                                     # no standard errors
                                     bootstrap_samples = 0) #salen problemas con los boundaries

summary(modelo.herb.intra.inter)
modelo.herb.intra.inter$log_likelihood 
###(aquí iría copia1)-----

#3.fv----
# covariate: floral visitors 
fv <- cov.list[my.sp2]
# discard the plant column          
for(i in 1:length(fv)){
    fv[[i]] <- fv[[i]][,2]
}
head(mis_species_junto[[1]])
# number of fitness observations
nrow(mis_species_junto[[1]])
# salinity data
head(fv[[1]])
# number of salinity observations
nrow(fv[[1]])
modelo.fv <- cxr_pm_multifit(data = mis_species_junto,
                               focal_column = my.sp2,
                               model_family = "BH",
                               # here we use a bounded method for demonstration purposes
                               optimization_method = "bobyqa", 
                               covariates = fv,
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
                                                   lambda_cov = -1,
                                                   alpha_cov = -1),
                               upper_bounds = list(lambda = 100,
                                                   alpha_intra = 1,
                                                   alpha_inter = 1,
                                                   lambda_cov = 1,
                                                   alpha_cov = 1),
                               # no standard errors
                               bootstrap_samples = 0) #salen problemas con los boundaries

#summary(modelo.fv)
modelo.fv$log_likelihood 

modelo.fv.intra.inter <- cxr_pm_multifit(data = mis_species_intra_inter,
                             focal_column = my.sp2,
                             model_family = "BH",
                             # here we use a bounded method for demonstration purposes
                             optimization_method = "bobyqa", 
                             covariates = fv,
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
                                                 lambda_cov = -1,
                                                 alpha_cov = -1),
                             upper_bounds = list(lambda = 100,
                                                 alpha_intra = 1,
                                                 alpha_inter = 1,
                                                 lambda_cov = 1,
                                                 alpha_cov = 1),
                             # no standard errors
                             bootstrap_samples = 0) #salen problemas con los boundaries

#summary(modelo.fv.intra.inter)
modelo.fv.intra.inter$log_likelihood 

#extraer outputs----
###extraer valores de los outputs de los modelos para ponerlos como valores iniciales en 
#el modelo global

### 1.1 extrar datos de sal----
#outputspara todas las especies juntas y separadas (inter,intra)
summary(modelo.sal)
#modelo.sal$alpha_matrix
matriz.sal.junto <- modelo.sal$alpha_matrix
matriz.sal.intrainter <- modelo.sal.inter.intra$alpha_matrix
n <- rownames(matriz.sal.junto)
v1 <- list()
v2 <- list()
v4 <- list()
for (i in 1: length(n)){ 
  v1 [[i]]<- matriz.sal.junto[which(rownames(matriz.sal.junto)== n [i]), which (colnames(matriz.sal.junto)== n [i])]
  v2 <- matriz.sal.junto[which(rownames(matriz.sal.junto) != n [i]), which (colnames(matriz.sal.junto) != n [i])]
  v4 [[i]] <- mean(v2 [i]) #medias inter por sp
}
sal_alpha_junto <- do.call("rbind",v4)
intra_sal_junto <- do.call("rbind",v1)
alpha_sal_junto <- do.call(rbind, Map(data.frame, alpha_inter_junto=sal_alpha_junto, alpha_intra_junto=intra_sal_junto))
alpha_sal_junto <- cbind(alpha_sal_junto, Species=n)

d <- cbind(alpha_sal_junto, matriz.sal.intrainter)
d$alpha_inter_intrainter <- d$inter
d$alpha_intra_intrainter <- d$intra
alpha_sal <- d[,c("Species", "alpha_inter_junto","alpha_intra_junto", "alpha_inter_intrainter", "alpha_intra_intrainter")]

lambda_sal_junto <- modelo.sal$lambda
lambda_sal_interintra <- modelo.sal.inter.intra$lambda
h <- modelo.sal$alpha_cov
h1 <- modelo.sal.inter.intra$alpha_cov

w <- do.call(rbind.data.frame, h)
alpha_cov_sal_juntas <- mean(as.matrix(w)) #media de alphas cov species juntas
w1 <- do.call(rbind.data.frame, h1)
alpha_cov_sal_interintra <- mean(as.matrix(w1)) #media de alphas cov interintra


########    2.1  extraer datos de los herb----

#con especies juntas y separadas (inter , intra)
summary(modelo.herb.junto)
#modelo.sal$alpha_matrix
matriz.herb.junto <- modelo.herb.junto$alpha_matrix
matriz.herb.interintra <- modelo.herb.intra.inter$alpha_matrix
n.herb <- rownames(matriz.herb.junto)
v1.herb <- list()
v2.herb <- list()
v4.herb <- list()
for (i in 1: length(n.herb)){ 
  v1.herb [[i]]<- matriz.herb.junto[which(rownames(matriz.herb.junto)== n.herb [i]), 
                                   which (colnames(matriz.herb.junto)== n [i])]
  v2.herb <- matriz.herb.junto[which(rownames(matriz.herb.junto) != n.herb [i]), 
                              which (colnames(matriz.herb.junto) != n.herb [i])]
  v4.herb [[i]] <- mean(v2.herb [i]) #medias inter por sp
}
herb_alpha_junto <- do.call("rbind",v4.herb)
intra_herb_junto <- do.call("rbind",v1.herb)
alpha_herb_junto <- do.call(rbind, Map(data.frame, alpha_inter_junto=herb_alpha_junto, alpha_intra_junto=intra_herb_junto))
alpha_herb_junto <- cbind(alpha_herb_junto, Species=n.herb)

d.herb <- cbind(alpha_herb_junto, matriz.herb.interintra)
d.herb$alpha_inter_intrainter <- d.herb$inter
d.herb$alpha_intra_intrainter <- d.herb$intra
alpha_herb <- d.herb[,c("Species", "alpha_inter_junto","alpha_intra_junto", "alpha_inter_intrainter", "alpha_intra_intrainter")]

lambda_herb_junto <- modelo.herb.junto$lambda
lambda_herb_intrainter<-  modelo.herb.intra.inter$lambda
h.herb <- modelo.herb.junto$alpha_cov
h1.herb <- modelo.herb.intra.inter$alpha_cov

w.herb <- do.call(rbind.data.frame, h.herb)
alpha_cov_herb_juntas <- mean(as.matrix(w.herb)) #media de alphas cov species juntas
w1.herb <- do.call(rbind.data.frame, h1.herb)
alpha_cov_herb_interintra <- mean(as.matrix(w1.herb)) #media de alphas cov interintra

######## 3.1 extraer datos de los fv----
#con especies juntas y separadas (inter , intra)
summary(modelo.fv)
#modelo.sal$alpha_matrix
matriz.fv.junto <- modelo.fv$alpha_matrix
matriz.fv.interintra <- modelo.fv.intra.inter$alpha_matrix
n.fv <- rownames(matriz.fv.junto)
medias.intra.sp.1 <- list()
v2.fv <- list()
medias.inter.sp.1 <- list()
for (i in 1: length(n.fv)){ 
  medias.intra.sp.1 [[i]]<- matriz.fv.junto[which(rownames(matriz.fv.junto)== n.fv [i]), 
                                    which (colnames(matriz.fv.junto)== n [i])]
  v2.fv <- matriz.fv.junto[which(rownames(matriz.fv.junto) != n.fv [i]), 
                               which (colnames(matriz.fv.junto) != n.fv [i])]
  medias.inter.sp.1 [[i]] <- mean(v2.fv [i]) #medias inter por sp
}
fv_alpha_junto <- do.call("rbind",medias.inter.sp.1)
intra_fv_junto <- do.call("rbind",medias.intra.sp.1)
alpha_fv_junto <- do.call(rbind, Map(data.frame, alpha_inter_junto=fv_alpha_junto, alpha_intra_junto=intra_fv_junto))
alpha_fv_junto <- cbind(alpha_fv_junto, Species=n.fv)

d.fv <- cbind(alpha_fv_junto, matriz.fv.interintra)
d.fv$alpha_inter_intrainter <- d.fv$inter
d.fv$alpha_intra_intrainter <- d.fv$intra
alpha_fv <- d.fv[,c("Species", "alpha_inter_junto","alpha_intra_junto", "alpha_inter_intrainter", "alpha_intra_intrainter")]

lambda_fv_junto <- modelo.fv$lambda
lambda_fv_intrainter<-  modelo.fv.intra.inter$lambda
h.fv <- modelo.fv$alpha_cov
h1.fv <- modelo.fv.intra.inter$alpha_cov

w.fv <- do.call(rbind.data.frame, h.fv)
alpha_cov_fv_juntas <- mean(as.matrix(w.fv)) #media de alphas cov species juntas
w1.fv <- do.call(rbind.data.frame, h1.fv)
alpha_cov_fv_interintra <- mean(as.matrix(w1.fv)) #media de alphas cov interintra

#modificaciones de datos para el modelo. necesito tener solo 1 lambda por especie, asi que hago la media
mean.lambda.species.junto <- (lambda_sal_junto + lambda_herb_junto+lambda_fv_junto)/3
mean.lambda.species.junto.vec <- as.vector(mean.lambda.species.junto)
#alpha inter e intra junto----
mean.alpha.intra.junto <- (alpha_sal$alpha_intra_junto+alpha_herb$alpha_intra_junto+ alpha_fv$alpha_intra_junto)/3
mean.alpha.inter.junto <- (alpha_sal$alpha_inter_junto+alpha_herb$alpha_inter_junto+ alpha_fv$alpha_inter_junto)/3
#lambda junto
mean.lambda.junto <- (lambda_sal_junto+lambda_herb_junto+ lambda_fv_junto)/3
mean.lambda.junto.vec <- as.vector(mean.lambda.junto)

#datos lambda cov junto ----
lambda_cov_sal_junto <-mean(as.vector(modelo.sal$lambda_cov))
lambda_cov_herb_junto <- mean(as.vector(modelo.herb.junto$lambda_cov))
lambda_cov_fv_junto <- mean(as.vector(modelo.fv$lambda_cov))
lambda_cov_junto <- c(lambda_cov_sal_junto,lambda_cov_herb_junto,lambda_cov_fv_junto)

#alpha inter e intra de interintra----
mean.alpha.intra.interintra <- (alpha_sal$alpha_intra_intrainter+alpha_herb$alpha_intra_intrainter+ alpha_fv$alpha_intra_intrainter)/3
mean.alpha.inter.interintra <- (alpha_sal$alpha_inter_intrainter+alpha_herb$alpha_inter_intrainter+ alpha_fv$alpha_inter_intrainter)/3


#lambda intrainter ----
lambda_intrainter <- as.vector((lambda_sal_interintra+lambda_herb_intrainter+lambda_fv_intrainter)/3)

#lambda cov intrainter----
lambda_cov_sal_intrainter <-mean(as.vector(modelo.sal.inter.intra$lambda_cov))
lambda_cov_herb_intrainter <- mean(as.vector(modelo.herb.intra.inter$lambda_cov))
lambda_cov_fv_intrainter <- mean(as.vector(modelo.fv.intra.inter$lambda_cov))
lambda_cov_intrainter <- c(lambda_cov_sal_intrainter,lambda_cov_herb_intrainter,lambda_cov_fv_intrainter)


##########MODELO GLOBAL----
#function parameters
model_family <- "BH"
covariates <- cov1
# bobyqa is generally more robust than other bounded methods
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "global"
alpha_cov_form <- "global"
# note how lambda_cov and alpha_cov
# have different initial values for each covariate effect
# the commented assignations are also possible, 
# giving equal initial values to all parameters
initial_values_juntas = list(lambda = mean.lambda.junto.vec,
                      alpha_intra = mean.alpha.intra.junto,
                      alpha_inter = mean.alpha.inter.junto,
                      lambda_cov = lambda_cov_junto,
                      alpha_cov = c(alpha_cov_sal_juntas,alpha_cov_herb_juntas, alpha_cov_fv_juntas))


# same with boundaries
lower_bounds = list(lambda = 0,
                    alpha_intra = -1,
                    alpha_inter = -1,
                    lambda_cov = -1,
                    alpha_cov = c(-1,-1,-1))

upper_bounds = list(lambda = 10000,
                    alpha_intra = 1,
                    alpha_inter = 1,
                    lambda_cov = 1,
                    alpha_cov = c(1,1,1))


fixed_terms <- NULL
bootstrap_samples <- 1000

any(is.na(mis_species))
any(is.na(covariates))

fit_multi_cov.juntas <- cxr_pm_multifit(data = mis_species_junto,
                                 focal_column = my.sp2,
                                 model_family = model_family,
                                 covariates = covariates,
                                 optimization_method = optimization_method,
                                 alpha_form = alpha_form,
                                 lambda_cov_form = lambda_cov_form,
                                 alpha_cov_form = alpha_cov_form,
                                 initial_values = initial_values_juntas,
                                 lower_bounds = lower_bounds,
                                 upper_bounds = upper_bounds,
                                 fixed_terms = fixed_terms,
                                 bootstrap_samples = bootstrap_samples) #error: check the data, initial values and bounds
#                                   parameter fitting failed for all focal species

summary(fit_multi_cov)

#modelo general intrainter----
initial_values_intrainter = list(lambda = lambda_intrainter,
                             alpha_intra = mean.alpha.intra.interintra,
                             alpha_inter = mean.alpha.inter.interintra,
                             lambda_cov = lambda_cov_intrainter,
                             alpha_cov = c(alpha_cov_sal_interintra,alpha_cov_herb_interintra, alpha_cov_fv_interintra))

fit_multi_cov.intrainter <- cxr_pm_multifit(data = mis_species_intra_inter,
                                        focal_column = my.sp2,
                                        model_family = model_family,
                                        covariates = covariates,
                                        optimization_method = optimization_method,
                                        alpha_form = alpha_form,
                                        lambda_cov_form = lambda_cov_form,
                                        alpha_cov_form = alpha_cov_form,
                                        initial_values = initial_values_intrainter,
                                        lower_bounds = lower_bounds,
                                        upper_bounds = upper_bounds,
                                        fixed_terms = fixed_terms,
                                        bootstrap_samples = bootstrap_samples)
