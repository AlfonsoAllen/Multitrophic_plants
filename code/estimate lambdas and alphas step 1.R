library(tidyverse)
library(cxr)
library(ggplot2)

comp<-read.csv("data/total_comp_19_20.check.csv", header=T, sep=";")
#env<-read.csv("data/covariates_salt_h_fv_20_19.csv", header=T, sep=";")
load("C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/cov_total.Rda")
env <- total
all.sp<-unique(comp$focal)
all.sp<-all.sp[-c(2,15,16)]#he quitado bema, achi, anar

comp$neigh_total<-rowSums(comp[, 3:24], na.rm=T)
min(comp$neigh_total)
#first is to create a list per species. 
neigh_total<-rowSums(comp[, 3:24], na.rm=T)
comp$neigh_total <- neigh_total
comp1 <- subset(comp, neigh_total>=1) #minimo 1 vecino por especie focal
comp1 <- subset(comp1, select = - neigh_total)

# add the neighbours column
obs_total<- list()
for (i in 1:length(all.sp)){
  x <- as_tibble(subset(comp1, focal==all.sp[i]))
  x <- subset( x, select = - focal) 
  x <- x[, which(colSums(x, na.rm=T) != 0)] #remove no neighbours so pairs of species that are not seen to interact
  if(all.sp[i] %in% names(x) =="TRUE"){
    x<-x %>% select(fitness, all.sp[i], everything())
  }
  obs_total[[i]] <- x
}
names(obs_total) <- unique(all.sp) 


#1. create an overall alpha model with no covariates----

d1<-as.data.frame(cbind(comp1[, 2], neigh_total))
names(d1)<-c("fitness", "neighbours")
d1<-subset(d1, fitness<8000)
d1<-subset(d1, neighbours<58)
d1<-subset(d1, neighbours>=1)



fit1 <- cxr_pm_fit(data = d1,
                   model_family = "BH",
                   # here we use a bounded method for demonstration purposes
                   optimization_method = "bobyqa",
                   alpha_form = "global",
                   initial_values = list(lambda = 1000, alpha_inter = 0),
                   lower_bounds = list(lambda = 10,alpha_inter = -5),
                   upper_bounds = list(lambda = 2000,alpha_inter = 5),
                   bootstrap_samples = 0)
summary(fit1)
fit1$log_likelihood
#check wether makes sense fitted values 
ggplot(d1, aes(neighbours , fitness)) + 
  geom_point() +
  stat_function(fun = function(x) fit1$lambda/(1+fit1$alpha_inter*x), lwd = 1.5, colour = "blue")

#2. create an overall intra and inter with no covariates----
obs2<-list()
for(i in 1:length(obs_total)){
  if(names(obs_total[i]) %in% names(obs_total[[i]][2])=="TRUE"){
    obs2[[i]]<-cbind(obs_total[[i]][1:2],rowSums(obs_total[[i]][-c(1,2)], na.rm=T)) #there is here intra and inter
  } else {
    obs2[[i]]<-NULL #there is no intra observed
  }
}  
names(obs2) <- unique(all.sp) 
all.sp2<-names(obs2) # this is later to have the focal column

for(i in 1:length(obs2)){
  names(obs2[[i]])<-c("fitness", names(obs2[i]), "neighbours")
}  


fit2 <- cxr_pm_multifit(data = obs2,
                        focal_column = all.sp2,
                        model_family = "BH",
                        optimization_method = "bobyqa",
                        alpha_form = "pairwise",
                        initial_values = list(lambda = fit1$lambda, alpha_intra =0.1, alpha_inter = 0.1),
                        lower_bounds = list(lambda = 10, alpha_intra = -1, alpha_inter = -2),
                        upper_bounds = list(lambda = 2000, alpha_intra = 1, alpha_inter = 2),
                        bootstrap_samples = 0)

summary(fit2)
#matrix<-fit2$alpha_matrix
#colnames(matrix)<-sub("^X+", "", colnames(matrix)) #remove the x associated

#4. create a pairwise with no covariates----
#We stick to the species already analyzed

#which species are common to the two dataset, this will be used for the next analysis 

obs4 <- obs_total
all.sp4<-names(obs4)

fit2_matrix <- fit2$alpha_matrix

fit4<-list()
for(i in 1:length(all.sp4)){
  obs<-obs4[i]
  sp <- all.sp4[i]
  x <- as.data.frame(obs4[[i]][,1])
  lambda <- mean(x[,1])
  #inter <- 0.1 ### rellenar por lo que sale del modelo 2
  #intra <- 0.1 #### rellenar por lo que sale del modelo 2
  inter <- fit2_matrix[sp,which(colnames(fit2_matrix)== c("neighbours"))]
  intra <- fit2_matrix[sp,sp]

fit <- cxr_pm_multifit(data = obs,
                        focal_column = sp,
                        model_family = "BH",
                        optimization_method = "bobyqa",
                        alpha_form = "pairwise",
                        initial_values = list(lambda = lambda, alpha_intra =intra, alpha_inter = inter),
                        lower_bounds = list(lambda = 10, alpha_intra = -2, alpha_inter = -2),
                        upper_bounds = list(lambda = 4000, alpha_intra = 2, alpha_inter = 2),
                        bootstrap_samples = 0)

fit4[[i]]<- fit

}

names(fit4)<-all.sp4

#5.fit 5 include each covariable separately ----
obs5<-obs_total
all.sp5<-names(obs5)


x<-list()
for(i in 1:length(fit4)){
  y <- as.data.frame(fit4[[i]]$alpha_matrix)
  y <- y %>% select(all.sp5)
  x[[i]]<-y
}

sp_matrix<-do.call("rbind", x)


salinity <-list()
for (i in 1:length(all.sp5)){
  x <- as_tibble(subset(env, Plant==all.sp5[i]))
  salinity[[i]] <- subset(x , select= sal) 
}
names(salinity) <- unique(all.sp5) 

herb <-list()
for (i in 1:length(all.sp5)){
  x <- as_tibble(subset(env, Plant==all.sp5[i]))
  herb[[i]] <- subset(x , select= herb1)
}
names(herb) <- unique(all.sp5) 

pol <-list()
for (i in 1:length(all.sp5)){
  x <- as_tibble(subset(env, Plant==all.sp5[i]))
  pol[[i]] <- subset(x , select= fv1) 
}
names(pol) <- unique(all.sp5) 


lower_bounds = list(lambda = 10,
                    alpha_intra = -2,
                    alpha_inter = -2,
                    lambda_cov = -4,
                    alpha_cov = -4)

upper_bounds = list(lambda = 4000,
                    alpha_intra = 5,
                    alpha_inter = 5,
                    lambda_cov = 5,
                    alpha_cov = 5)



#5.1. salinity----
#voy a comprobar que tanto en competencia como en salinidad haya las mismas filas de datos
head(obs5[[1]])
# number of fitness observations
nrow(obs5[[2]])
# salinity data
head(salinity[[1]])
# number of covariate observations
nrow(salinity[[2]])
any(is.na(salinity)) # ho hay Nas ni en salinidad ni en obs5

fit5_salinity<-list()
for(i in 1:length(all.sp5)){
  obs<-obs5[i]
  sp <- all.sp5[i]
  env_list <-salinity[i]
  lambda <- fit4[[i]]$lambda
  inter <- sp_matrix[sp,which(colnames(sp_matrix)!=sp)]
  inter <-mean(inter[,1])
  intra <- sp_matrix[sp,sp]
  
  fit <- cxr_pm_multifit(data = obs,
                         focal_column = sp,
                         model_family = "BH",
                         covariates = env_list,
                         optimization_method = "bobyqa",
                         alpha_form = "pairwise",
                         lambda_cov_form = "global", # effect of covariates over lambda
                         alpha_cov_form = "global", # effect of covariates over alpha
                         initial_values = list(lambda = lambda, alpha_intra = intra, alpha_inter = inter,
                                               lambda_cov = 0, alpha_cov = 0),
                         lower_bounds = lower_bounds,
                         upper_bounds = upper_bounds,
                         bootstrap_samples = 0)
  
  fit5_salinity[[i]]<- fit
  
}


names(fit5_salinity)<-all.sp5


x.sal<-list()
for(i in 1:length(fit5_salinity)){
  y <- as.data.frame(fit5_salinity[[i]]$alpha_matrix)
  y <- y %>% select(all.sp5)
  x.sal[[i]]<-y
}

sp_matrix_sal<-do.call("rbind", x.sal)





# 5.2 Herb----
  
fit5_herb<-list()
for(i in 1:length(all.sp5)){
  obs<-obs5[i]
  sp <- all.sp5[i]
  env_list <-herb[i]
  lambda <- fit4[[i]]$lambda
  inter <- sp_matrix[sp,which(colnames(sp_matrix)!=sp)]
  inter <-mean(inter[,1])
  intra <- sp_matrix[sp,sp]
  
  fit <- cxr_pm_multifit(data = obs,
                         focal_column = sp,
                         model_family = "BH",
                         covariates = env_list,
                         optimization_method = "bobyqa",
                         alpha_form = "pairwise",
                         lambda_cov_form = "global", # effect of covariates over lambda
                         alpha_cov_form = "global", # effect of covariates over alpha
                         initial_values = list(lambda = lambda, alpha_intra =intra, alpha_inter = inter,
                                               lambda_cov = 0, alpha_cov = 0),
                         lower_bounds = lower_bounds,
                         upper_bounds = upper_bounds,
                         bootstrap_samples = 0)
  
  fit5_herb[[i]]<- fit
  
}


names(fit5_herb)<-all.sp5


x.herb<-list()
for(i in 1:length(fit5_herb)){
  y <- as.data.frame(fit5_herb[[i]]$alpha_matrix)
  y <- y %>% select(all.sp5)
  x.herb[[i]]<-y
}

sp_matrix_herb<-do.call("rbind", x.herb)





# 5.3 floral visitors ----
fit5_pol<-list()
for(i in 1:length(all.sp5)){
  obs<-obs5[i]
  sp <- all.sp5[i]
  env_list <-pol[i]
  lambda <- fit4[[i]]$lambda
  inter <- sp_matrix[sp,which(colnames(sp_matrix)!=sp)]
  inter <-mean(inter[,1])
  intra <- sp_matrix[sp,sp]
  
  fit <- cxr_pm_multifit(data = obs,
                         focal_column = sp,
                         model_family = "BH",
                         covariates = env_list,
                         optimization_method = "bobyqa",
                         alpha_form = "pairwise",
                         lambda_cov_form = "global", # effect of covariates over lambda
                         alpha_cov_form = "global", # effect of covariates over alpha
                         initial_values = list(lambda = lambda, alpha_intra =intra, alpha_inter = inter,
                                               lambda_cov = 0, alpha_cov = 0),
                         lower_bounds = lower_bounds,
                         upper_bounds = upper_bounds,
                         bootstrap_samples = 0)
  
  fit5_pol[[i]]<- fit
  
}

names(fit5_pol)<-all.sp5



x.pol<-list()
for(i in 1:length(fit5_pol)){
  y <- as.data.frame(fit5_pol[[i]]$alpha_matrix)
  y <- y %>% select(all.sp5)
  x.pol[[i]]<-y
}#esto es para hacer las matrices cuadradas

sp_matrix_pol<-do.call("rbind", x.pol) 



env_total <- list()
# 6. model with the three covariates

for (i in 1:length(all.sp5)){
  x <- as_tibble(subset(env, Plant==all.sp5[i]))
  env_total[[i]] <- subset(x , select = c(sal, herb1, fv1))
}

names(env_total)<-all.sp5


obs6 <- obs_total





fit6<-list()
lambda_sal <- list()
lambda_herb <- list()
lambda_pol <- list()
alpha_cov_sal <- list()
mean_alpha_cov_sal <- list()
alpha_cov_herb<- list()
alpha_cov_pol <- list()
inter_sal <- list()
inter_herb <- list()
inter_pol <- list()
intra_sal <- list()
intra_herb <- list()
intra_pol <- list()
lambda_cov <- list()

for(i in 1:length(all.sp5)){
  obs<-obs6[i]
  sp <- all.sp5[i]
  env_list <-env_total[i]
  #lambda <- mean(c(fit5_salinity[[i]]$lambda, fit5_herb[[i]]$lambda, fit5_pol[[i]]$lambda))
  lambda_sal [[i]] <- fit5_salinity[[i]]$lambda
  lambda_herb [[i]] <- fit5_herb[[i]]$lambda
  lambda_pol [[i]] <- fit5_pol[[i]]$lambda
  inter_sal [[i]] <- mean(as.numeric(sp_matrix_sal[sp, which(colnames(sp_matrix_sal)!=sp)]))
  inter_herb [[i]] <-  mean(as.numeric(sp_matrix_herb[sp, which(colnames(sp_matrix_herb)!=sp)]))      
  inter_pol [[i]] <- mean(as.numeric(sp_matrix_pol[sp, which(colnames(sp_matrix_pol)!=sp)]))    
  intra_sal[[i]] <- sp_matrix_sal[sp,sp]
  intra_herb [[i]]<- sp_matrix_herb[sp,sp]
  intra_pol [[i]] <- sp_matrix_pol[sp,sp]
  lambda_cov [[i]]<- c((mean(fit5_salinity[[i]]$lambda_cov)), mean(fit5_herb[[i]]$lambda_cov), 
                  (mean(fit5_pol[[i]]$lambda_cov)))
 alpha_cov_sal [[i]] <- ((fit5_salinity[[i]]$alpha_cov)$sal)
 alpha_cov_sal [[i]]  <- mean(alpha_cov_sal [[i]])
 alpha_cov_herb [[i]]<- ((fit5_herb[[i]]$alpha_cov)$herb1)
 alpha_cov_herb [[i]]<- mean(alpha_cov_herb [[i]])
 alpha_cov_pol [[i]]<- ((fit5_pol[[i]]$alpha_cov)$fv1)
 alpha_cov_pol[[i]] <- mean(alpha_cov_pol [[i]])
 cova  <- data.frame(cbind(alpha_cov_sal,alpha_cov_herb, alpha_cov_pol))
 cova <- c(mean(cova$alpha_cov_sal), mean(cova$alpha_cov_herb), mean(cova$alpha_cov_pol))
 
  
}

#6.1 sacar parametros para el modelo global
#lambda
lambda_sal <- do.call("rbind",lambda_sal)
lambda_sal <- data.frame(lambda_sal)
names(lambda_sal)<- c("lambda_sal")

lambda_herb <- do.call("rbind",lambda_herb)
lambda_herb <- data.frame(lambda_herb)
names(lambda_herb)<- c("lambda_herb")

lambda_pol <- do.call("rbind",lambda_pol)
lambda_pol <- data.frame(lambda_pol)
names(lambda_pol)<- c("lambda_pol")

lambdas <- cbind(lambda_sal,lambda_herb,lambda_pol)
lambdas$lambda_sal <-mean(as.numeric(lambdas$lambda_sal))
lambdas$lambda_herb <-mean(as.numeric(lambdas$lambda_herb))
lambdas$lambda_pol <-mean(as.numeric(lambdas$lambda_pol))
lambdas <- unique(lambdas)
lambdas_mean <- mean(as.numeric(lambdas))

#alpha inter

inter_sal <- do.call("rbind",inter_sal)
inter_sal <- data.frame(inter_sal)
names(inter_sal)<- c("inter_sal")

inter_herb <- do.call("rbind",inter_herb)
inter_herb <- data.frame(inter_herb)
names(lambda_herb)<- c("inter_herb")

inter_pol <- do.call("rbind",inter_pol)
inter_pol <- data.frame(inter_pol)
names(inter_pol)<- c("inter_pol")

inter_alpha<- cbind(inter_sal,inter_herb,inter_pol)
inter_alpha$inter_sal <-mean(as.numeric(inter_alpha$inter_sal))
inter_alpha$inter_herb <-mean(as.numeric(inter_alpha$inter_herb))
inter_alpha$inter_pol <-mean(as.numeric(inter_alpha$inter_pol))
inter_alpha <- unique(inter_alpha)
inter_alpha_mean <- mean(as.numeric(inter_alpha))


#intra alpha
intra_sal <- do.call("rbind",intra_sal)
intra_sal <- data.frame(intra_sal)
names(intra_sal)<- c("intra_sal")

intra_herb <- do.call("rbind",intra_herb)
intra_herb <- data.frame(intra_herb)
names(intra_herb)<- c("intra_herb")

intra_pol <- do.call("rbind",intra_pol)
intra_pol <- data.frame(intra_pol)
names(intra_pol)<- c("intra_pol")

intra_alpha<- cbind(intra_sal,intra_herb,intra_pol)
intra_alpha$intra_sal <-mean(as.numeric(intra_alpha$intra_sal))
intra_alpha$intra_herb <-mean(as.numeric(intra_alpha$intra_herb))
intra_alpha$intra_pol <-mean(as.numeric(intra_alpha$intra_pol))
intra_alpha <- unique(intra_alpha)
intra_alpha_mean <- mean(as.numeric(intra_alpha))
#lambda cov
lambda_cov_1 <- do.call("rbind",lambda_cov)
lambda_cov_salt <- mean(as.numeric(lambda_cov_1[,1]))
lambda_cov_herb <- mean(as.numeric(lambda_cov_1[,2]))
lambda_cov_pol <- mean(as.numeric(lambda_cov_1[,3]))
lambda_cov_vector <- c(lambda_cov_salt,lambda_cov_herb,lambda_cov_pol )

#alpha_cov

alpha_cov_sal <- do.call("rbind",alpha_cov_sal)
alpha_cov_sal <- data.frame(alpha_cov_sal)
names(alpha_cov_sal)<- c("alpha_cov_sal")

alpha_cov_herb <- do.call("rbind",alpha_cov_herb)
alpha_cov_herb <- data.frame(alpha_cov_herb)
names(alpha_cov_herb)<- c("alpha_cov_herb")

alpha_cov_pol <- do.call("rbind",alpha_cov_pol)
alpha_cov_pol <- data.frame(alpha_cov_pol)
names(alpha_cov_pol)<- c("alpha_cov_pol")

alpha_cov_vector<- cbind(alpha_cov_sal,alpha_cov_herb,alpha_cov_pol)
mean_alpha_cov_sal <-mean(as.numeric(alpha_cov_vector$alpha_cov_sal))
mean_alpha_cov_herb <-mean(as.numeric(alpha_cov_vector$alpha_cov_herb))
mean_alpha_cov_pol <-mean(as.numeric(alpha_cov_vector$alpha_cov_pol))
alpha_cov_final<- c(mean_alpha_cov_sal,mean_alpha_cov_herb,mean_alpha_cov_pol)


names(alpha_cov_final) <- c("salt","herb", "pol")
names(lambda_cov_vector) <- c("salt","herb", "pol")

initial_values = list(lambda = lambdas_mean,
                                 alpha_intra = intra_alpha_mean,
                                 alpha_inter = inter_alpha_mean,
                                 lambda_cov = lambda_cov_vector,
                                 alpha_cov = alpha_cov_final)

# same with boundaries
lower_bounds = list(lambda = 10,
                    alpha_intra = -2,
                    alpha_inter = -2,
                    lambda_cov = -4,
                    alpha_cov = c(-4, -4,-4))

upper_bounds = list(lambda = 4000,
                    alpha_intra = 5,
                    alpha_inter = 5,
                    lambda_cov = 5,
                    alpha_cov = c(5,5,5))
#function parameters
model_family <- "BH"
covariates <- env_total
# bobyqa is generally more robust than other bounded methods
optimization_method <- "bobyqa"
alpha_form <- "pairwise"
lambda_cov_form <- "global"
alpha_cov_form <- "global"
bootstrap_samples <- 0
fixed_terms <- NULL

fit_multi_cov <- cxr_pm_multifit(data = obs6,
                                            focal_column = all.sp5,
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
                                            bootstrap_samples = bootstrap_samples)

summary(fit_multi_cov)

x.multi<-list()
for(i in 1:length(fit_multi_cov)){
  y <- as.data.frame(fit_multi_cov$alpha_matrix)
  y <- y %>% select(all.sp5)
  x.multi<-y
}#esto es para hacer las matrices cuadradas



x <- d1$neighbours
caso.lema.inter <- (x.multi["LEMA",which(colnames(x.multi)!="LEMA")])
caso.lema.inter <- rowSums(caso.lema.inter)
lambda_cov_lema <- data.frame(fit_multi_cov$lambda_cov[1,])
fv_lema_lambda_cov <- lambda_cov_lema[1,]

alpha_cov.plot <- fit_multi_cov$alpha_cov
alpha_cov.plotfv <- alpha_cov.plot$fv1["LEMA","LEMA"]



ggplot(d1, aes(neighbours , fitness)) + 
  geom_point() +
  stat_function(fun = function(x) fit_multi_cov$lambda["LEMA"]/(1+caso.lema.inter*x), 
                lwd = 1.5, colour = "blue") +
  stat_function(fun = function(x) fit_multi_cov$lambda["LEMA"]*(1+(fv_lema_lambda_cov)*8)/((1+caso.lema.inter+alpha_cov.plotfv*8)*x), 
                lwd = 1.5, colour = "red")

