library(tidyverse)
library(cxr)
library(ggplot2)


comp<-read.csv("total_comp_19_20.check.csv", header=T, sep=";")
env<-read.csv("covariates_salt_h_fv_20_19.csv", header=T, sep=";")
all.sp<-unique(comp$focal)
all.sp<-all.sp[-c(15:16)]

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
                   lower_bounds = list(lambda = 10, alpha_inter = -5),
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


  



