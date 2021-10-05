library(tidyverse)
library(cxr)
library(ggplot2)


comp<-read.csv("data/total_comp_19_20.check.csv", header=T, sep=";")
#env<-read.csv("covariates_salt_h_fv_20_19.csv", header=T, sep=";")
all.sp<-unique(comp$focal)
all.sp<-all.sp[-c(15:16)]
load("data/cov_total.Rda")
env <- total

#first is to create a list per species. 

# add the neighbours column
obs_total<- list()
for (i in 1:length(all.sp)){
    x <- as_tibble(subset(comp, focal==all.sp[i]))
    x <- subset( x, select = - focal) 
    x <- x[, which(colSums(x, na.rm=T) != 0)] #remove no neighbours so pairs of species that are not seen to interact
    if(all.sp[i] %in% names(x) =="TRUE"){
        x<-x %>% select(fitness, all.sp[i], everything())
    }
    obs_total[[i]] <- x
}
names(obs_total) <- unique(all.sp) 


#1. create an overall alpha model with no covariates----

neigh_total<-rowSums(comp[, 3:24], na.rm=T)
d1<-as.data.frame(cbind(comp[, 2], neigh_total))
names(d1)<-c("fitness", "neighbours")
d1<-subset(d1, fitness<8000)
d1<-subset(d1, neighbours<58)
d1<-subset(d1, neighbours>=1)



fit1 <- cxr_pm_fit(data = d1,
                   model_family = "BH",
                   # here we use a bounded method for demonstration purposes
                   optimization_method = "bobyqa",
                   alpha_form = "global",
                   initial_values = list(lambda = 1000, alpha_inter = 0.001),
                   lower_bounds = list(lambda = 10, alpha_inter = 0.001),
                   upper_bounds = list(lambda = 2000, alpha_inter = 5),
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
    } 
}  
names(obs2) <- unique(all.sp) 
all.sp2<-names(obs2) # this is later to have the focal column

for(i in 1:length(obs2)){
    names(obs2[[i]])<-c("fitness", names(obs2[i]), "neighbours")
}  


fit2_intra<-list()
for(i in 1:length(obs2)){
    obs<-obs2[[i]][c(1,2)]
    sp<-all.sp2[i]
    lambda<-mean(obs$fitness)
    
    fit <- cxr_pm_fit(data = obs,
                      #focal_column = sp,
                      model_family = "BH",
                      optimization_method = "bobyqa",
                      alpha_form = "global",
                      initial_values = list(lambda = lambda, alpha_inter = 0.5),
                      lower_bounds = list(lambda = 10, alpha_inter = 0.001),
                      upper_bounds = list(lambda = 3000, alpha_inter = 10),
                      bootstrap_samples = 0)
    fit2_intra[[i]]<-fit
}

names(fit2_intra)<-all.sp2

fit2_inter<-list()
for(i in 1:length(obs2)){
    obs<-obs2[[i]][c(1,3)]
    sp<-all.sp2[i]
    lambda<-mean(obs$fitness)
    
    fit <- cxr_pm_fit(data = obs,
                      #focal_column = sp,
                      model_family = "BH",
                      optimization_method = "bobyqa",
                      alpha_form = "global",
                      initial_values = list(lambda = lambda, alpha_inter = 0.5),
                      lower_bounds = list(lambda = 10, alpha_inter = 0.001),
                      upper_bounds = list(lambda = 3000, alpha_inter = 10),
                      bootstrap_samples = 0)
    fit2_inter[[i]]<-fit
}

names(fit2_inter)<-all.sp2


#matrix<-fit2$alpha_matrix
#colnames(matrix)<-sub("^X+", "", colnames(matrix)) #remove the x associated

#4. create a pairwise with no covariates----
#We stick to the species already analyzed

#which species are common to the two dataset, this will be used for the next analysis 

obs4 <- obs_total
all.sp4<-names(obs4)

fit4<-list()
for(i in 1:length(all.sp4)){
    obs<-obs4[i]
    sp <- all.sp4[i]
    x <- as.data.frame(obs4[[i]][,1])
    lambda <- mean(fit2_inter[[1]]$lambda, fit2_intra[[1]]$lambda)
    inter <- fit2_inter[[1]]$alpha_inter
    intra <- fit2_intra[[1]]$alpha_inter
    
    fit <- cxr_pm_multifit(data = obs,
                           focal_column = sp,
                           model_family = "BH",
                           optimization_method = "bobyqa",
                           alpha_form = "pairwise",
                           initial_values = list(lambda = lambda, alpha_intra =intra, alpha_inter = inter),
                           lower_bounds = list(lambda = 10, alpha_intra = 0.001, alpha_inter = 0.001),
                           upper_bounds = list(lambda = 4000, alpha_intra = 10, alpha_inter = 10),
                           bootstrap_samples = 0)
    
    fit4[[i]]<- fit
    
}


names(fit4)<-all.sp4

#fit 5 include each covariable separately ----
obs5<-obs_total
all.sp5<-all.sp

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
    salinity[[i]] <- subset(x , select= sal) #divided by 1000 to allow models converge later on
}
names(salinity) <- unique(all.sp5) 

herb <-list()
for (i in 1:length(all.sp5)){
    x <- as_tibble(subset(env, Plant==all.sp5[i]))
    herb[[i]] <- subset(x , select= herb1) #divided by 1000 to allow models converge later on
}
names(herb) <- unique(all.sp5) 

pol <-list()
for (i in 1:length(all.sp5)){
    x <- as_tibble(subset(env, Plant==all.sp5[i]))
    pol[[i]] <- subset(x , select= fv1) #divided by 1000 to allow models converge later on
}
names(pol) <- unique(all.sp5) 


lower_bounds = list(lambda = 10,
                    alpha_intra = 0.001,
                    alpha_inter = 0.001,
                    lambda_cov = -4,
                    alpha_cov = -4)

upper_bounds = list(lambda = 4000,
                    alpha_intra = 10,
                    alpha_inter = 10,
                    lambda_cov = 4,
                    alpha_cov = 4)


nrow(obs5[[1]])
nrow(salinity[[1]])

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
                           initial_values = list(lambda = lambda, alpha_intra =intra, alpha_inter = inter,
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
                                                 lambda_cov = 1, alpha_cov = 1),
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
    lower_bounds <- if(unique(env_list[[1]]$fv1==0)==TRUE){list(lambda = 10,
                                                                alpha_intra = 0.001,
                                                                alpha_inter = 0.001,
                                                                lambda_cov = 0.000001,
                                                                alpha_cov = 0.000001)}else{
      list(lambda = 10, alpha_intra = 0.001,alpha_inter = 0.001, lambda_cov = -4, alpha_cov = -4)}
    
    upper_bounds <- if(unique(env_list[[1]]$fv1==0)==TRUE){list(lambda = 4000,
                                                                alpha_intra = 10,
                                                                alpha_inter = 10,
                                                                lambda_cov = 0.000002,
                                                                alpha_cov = 0.000002)}else{
      list(lambda = 4000, alpha_intra = 10,alpha_inter = 10, lambda_cov = 4, alpha_cov = 4)}
    
    initial_values <- if(unique(env_list[[1]]$fv1==0)==TRUE){list(lambda = lambda,
                                                                alpha_intra = intra,
                                                                alpha_inter = inter,
                                                                lambda_cov = 0.000001,
                                                                alpha_cov = 0.000001)}else{
                                                                  list(lambda = lambda, alpha_intra = intra,alpha_inter = inter, lambda_cov = 1, alpha_cov = 1)}
    
    
    fit <- cxr_pm_multifit(data = obs,
                           focal_column = sp,
                           model_family = "BH",
                           covariates = env_list,
                           optimization_method = "bobyqa",
                           alpha_form = "pairwise",
                           lambda_cov_form = "global", # effect of covariates over lambda
                           alpha_cov_form = "global", # effect of covariates over alpha
                           initial_values = initial_values,
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

#6. modelo global----

obs6 <- obs_total

fit6<-list()
for (i in 1:length(all.sp5)){
  obs<-obs6[i]
  sp <- all.sp5[i]
  env_list <- env_total[i]
  lambda <- mean(fit5_salinity[[i]]$lambda, fit5_herb[[i]]$lambda, fit5_pol[[i]]$lambda)
  inter <-mean(mean(fit5_salinity[[i]]$alpha_matrix), mean(fit5_herb[[i]]$alpha_matrix), mean(fit5_pol[[i]]$alpha_matrix))
  intra <- mean(fit5_salinity[[i]]$alpha_matrix[,sp], fit5_herb[[i]]$alpha_matrix[,sp], fit5_pol[[i]]$alpha_matrix[,sp])
  
  lambda_cov<-c(if(is.null(fit5_salinity[[i]]$lambda_cov)==TRUE){1}else{fit5_salinity[[i]]$lambda_cov},
                if(is.null(fit5_herb[[i]]$lambda_cov)==TRUE){1}else{fit5_herb[[i]]$lambda_cov},
                if(is.null(fit5_pol[[i]]$lambda_cov)==TRUE){1}else{fit5_pol[[i]]$lambda_cov})
  alpha_cov <-c(if(is.na(mean(fit5_salinity[[i]]$alpha_cov$sal))==TRUE){1}else{mean(fit5_salinity[[i]]$alpha_cov$sal)},
                if(is.na(mean(fit5_herb[[i]]$alpha_cov$herb1))==TRUE){1}else{mean(fit5_herb[[i]]$alpha_cov$herb1)},
                if(is.na(mean(fit5_pol[[i]]$alpha_cov$fv1))==TRUE){1}else{mean(fit5_pol[[i]]$alpha_cov$fv1)})
  lambda_cov<-if(any(lambda_cov<0)==TRUE){c(1,1,1)}else{lambda_cov}
  alpha_cov <-if(any(alpha_cov<0)==TRUE){c(1,1,1)}else{alpha_cov}
  
  lower_bounds <- if(unique(env_list[[1]]$fv1==0)==TRUE){list(lambda = 10,
                                                              alpha_intra = 0.001,
                                                              alpha_inter = 0.001,
                                                              lambda_cov = c(-4, -4, 0.000001), alpha_cov = c(-4, -4, 0.000001))
  }else {list(lambda = 10, alpha_intra = 0.001,alpha_inter = 0.001, lambda_cov =c(-4, -4, -4), alpha_cov = c(-4, -4, -4))}

  upper_bounds <- if(unique(env_list[[1]]$fv1==0)==TRUE){list(lambda = 4000,
                                                              alpha_intra = 10,
                                                              alpha_inter = 10,
                                                              lambda_cov = c(4, 4, 0.000002), alpha_cov = c(4, 4, 0.000002))
  }else {list(lambda = 4000, alpha_intra = 10,alpha_inter = 10, lambda_cov =c(4, 4, 4), alpha_cov = c(4, 4, 4))}
  
  initial_values <- if(unique(env_list[[1]]$fv1==0)==TRUE){list(lambda = lambda,
                                                                alpha_intra = intra,
                                                                alpha_inter = inter,
                                                                lambda_cov = 0.000001,
                                                                alpha_cov = 0.000001)}else{
                                                                  list(lambda = lambda, alpha_intra = intra,alpha_inter = inter, lambda_cov = lambda_cov, alpha_cov = alpha_cov)}
  
  fit <- cxr_pm_multifit(data = obs,
                         focal_column = sp,
                         model_family = "BH",
                         covariates = env_list,
                         optimization_method = "bobyqa",
                         alpha_form = "pairwise",
                         lambda_cov_form = "global", # effect of covariates over lambda
                         alpha_cov_form = "global", # effect of covariates over alpha
                         initial_values = initial_values,
                         lower_bounds = lower_bounds,
                         upper_bounds = upper_bounds,
                         bootstrap_samples = 0)
  
     fit6[[i]]<- fit
  
}

names(fit6)<-all.sp5
  




#selecionar parametros del modelo a guardar
#6.1 alpha matrix----
alpha <- list()
for (i in 1:length(all.sp5)) {
  alpha [[i]] <- (as_tibble(fit6[[i]]$alpha_matrix)) %>% select(all.sp5)
  
}

alpha.matrix <- data.frame(do.call("rbind",alpha))
colnames(alpha.matrix) <- all.sp5
rownames(alpha.matrix) <- all.sp5
#save(alpha.matrix, file="data/alpha_matrix_final.Rda")


#6.2 lambda ----
lambda <- list()
for (i in 1:length(all.sp5)) {
  lambda [[i]] <- (as_tibble(fit6[[i]]$lambda)) 
  
}
lambdas <- data.frame(do.call("rbind", lambda))
colnames(lambdas) <- "lambda"
rownames(lambdas) <- all.sp5
#save(lambdas, file="data/lambdas_final.Rda")

#6.3 lambda cov ----
lambda_cov <- list()
for (i in 1:length(all.sp5)) {
  lambda_cov [[i]] <- (as_tibble(fit6[[i]]$lambda_cov)) 
  
}
lambdas_cov <- data.frame(do.call("rbind", lambda_cov))
colnames(lambdas_cov) <- c("pol", "herb","salt")
rownames(lambdas_cov) <- all.sp5
#save(lambdas_cov, file="data/lambdas_cov_final.Rda")


#6.4 alpha cov ----
alpha_cov_pol <- list()
alpha_cov_herb <- list()
alpha_cov_salt <- list()
for (i in 1:length(all.sp5)) {
  alpha_cov_pol [[i]] <- mean(as.numeric(fit6[[i]]$alpha_cov$fv1)) 
  alpha_cov_herb [[i]] <- mean(as.numeric(fit6[[i]]$alpha_cov$herb1))
  alpha_cov_salt [[i]] <- mean(as.numeric(fit6[[i]]$alpha_cov$sal)) 
  
}
alpha_cov_pol <- data.frame(do.call("rbind", alpha_cov_pol))
alpha_cov_herb <- data.frame(do.call("rbind", alpha_cov_herb))
alpha_cov_salt <- data.frame(do.call("rbind", alpha_cov_salt))
alpha_cov <- cbind(alpha_cov_pol,alpha_cov_herb,alpha_cov_salt)
colnames(alpha_cov) <- c("pol", "herb","salt")
rownames(alpha_cov) <- all.sp5
#save(alpha_cov, file="data/alphas_cov_final.Rda")


# 7. Plot de fitness y neighbors con pol, salt y herb ----
x <- d1$neighbours
env_list.herb<-list()
env_list.pol <- list()
caso.inter <- list()
caso.intra <- list()
lambda_cov <- list()
lambda_cov_sal <- list()
lambda_cov_herb <- list()
lambda_cov_fv <- list()
alpha_cov.d <- list()
alpha_cov.sal <- list()
alpha_cov.herb <- list()
alpha_cov.pol <- list()
lambda_sp1 <- list()
g.pol <- list()
g.sal <- list()
g.herb <- list()

for(i in 1:length(all.sp5)){
  obs <- obs6[i]
  sp <- all.sp5[i]
  env_list <- env_total[i]
  env_list.herb [[i]] <- env_total[[i]][,"herb1"]
  env_list.pol [[i]] <- env_total[[i]][,"fv1"]
  caso.inter [[i]] <- mean(as.numeric(alpha.matrix[sp,which(colnames(alpha.matrix)!= sp)]))
  lambda_cov[[i]] <- (lambdas_cov[sp,])
  lambda_cov_sal[[i]] <- (lambdas_cov[sp,"salt"])
  lambda_cov_herb[[i]] <- (lambdas_cov[sp,"herb"])
  lambda_cov_fv [[i]] <- (lambdas_cov[sp,"pol"])
  alpha_cov.d [[i]] <- (alpha_cov[sp,])
  alpha_cov.sal [[i]]<- (alpha_cov[sp,"salt"])
  alpha_cov.herb[[i]] <- (alpha_cov[sp,"herb"])
  alpha_cov.pol[[i]] <- (alpha_cov[sp,"pol"])
  
  lambda_sp1 [[i]]<- (lambdas[sp,])
  
  #for (j in  1: length(env_list.pol[[i]])) {
  g.pol[[i]]<-  ggplot(d1, aes(neighbours , fitness)) + 
    geom_point() +
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])/(1+(as.numeric(caso.inter[i]))*x), 
                  lwd = 1.5, colour = "blue") +
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_fv[i]))*2)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.pol[i]))*2))*x), 
                  lwd = 1.5, colour = "red")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_fv[i]))*6)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.pol[i]))*6))*x), 
                  lwd = 1.5, colour = "green")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_fv[i]))*10)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.pol[i]))*10))*x), 
                  lwd = 1.5, colour = "pink")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_fv[i]))*15)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.pol[i]))*15))*x), 
                  lwd = 1.5, colour = "yellow")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_fv[i]))*21)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.pol[i]))*21))*x), 
                  lwd = 1.5, colour = "purple")+
    ylim(0,200)
  
  g.herb[[i]]<-  ggplot(d1, aes(neighbours , fitness)) + 
    geom_point() +
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])/(1+(as.numeric(caso.inter[i]))*x), 
                  lwd = 1.5, colour = "blue") +
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_herb[i]))*3)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.herb[i]))*3))*x), 
                  lwd = 1.5, colour = "red")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_herb[i]))*8)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.herb[i]))*8))*x), 
                  lwd = 1.5, colour = "green")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_herb[i]))*15)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.herb[i]))*15))*x), 
                  lwd = 1.5, colour = "pink")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_herb[i]))*23)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.herb[i]))*23))*x), 
                  lwd = 1.5, colour = "yellow")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_herb[i]))*31)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.herb[i]))*31))*x), 
                  lwd = 1.5, colour = "purple")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_herb[i]))*48)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.herb[i]))*48))*x), 
                  lwd = 1.5, colour = "orange")+
    ylim(0,200)
  
  g.sal[[i]]<-  ggplot(d1, aes(neighbours , fitness)) + 
    geom_point() +
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])/(1+(as.numeric(caso.inter[i]))*x), 
                  lwd = 1.5, colour = "blue") +
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_sal[i]))*0.0894)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.sal[i]))*0.0894))*x), 
                  lwd = 1.5, colour = "red")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_sal[i]))*0.6966)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.sal[i]))*0.6966))*x), 
                  lwd = 1.5, colour = "green")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_sal[i]))*0.87)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.sal[i]))*0.87))*x), 
                  lwd = 1.5, colour = "pink")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_sal[i]))*1.1670)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.sal[i]))*1.1670))*x), 
                  lwd = 1.5, colour = "yellow")+
    stat_function(fun = function(x) as.numeric(lambda_sp1[i])*(1+(as.numeric(lambda_cov_sal[i]))*1.977)/
                    ((1+(as.numeric(caso.inter[i]))+((as.numeric(alpha_cov.sal[i]))*1.977))*x), 
                  lwd = 1.5, colour = "purple")+
    ylim(0,200)
  
  
}

