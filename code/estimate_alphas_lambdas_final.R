library(tidyverse)
library(cxr)
library(ggplot2)


comp<-read.csv("data/total_comp_19_20.check.csv", header=T, sep=";")
#env<-read.csv("covariates_salt_h_fv_20_19.csv", header=T, sep=";")
all.sp<-unique(comp$focal)
all.sp<-all.sp[-c(15:16)]
load("C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/cov_total.Rda")
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
alpha_cov_sal1 <- list()
alpha_cov_herb1 <- list()
alpha_cov_pol1 <- list()

for(i in 1:length(all.sp5)){
    obs<-obs6[i]
    sp <- all.sp5[i]
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
    alpha_cov_sal1 [[i]]  <- mean(alpha_cov_sal [[i]])
    alpha_cov_herb [[i]]<- ((fit5_herb[[i]]$alpha_cov)$herb1)
    alpha_cov_herb1 [[i]]<- mean(alpha_cov_herb [[i]])
    alpha_cov_pol [[i]]<- ((fit5_pol[[i]]$alpha_cov)$fv1)
    alpha_cov_pol1[[i]] <- mean(alpha_cov_pol [[i]])
    cova  <- data.frame(cbind(alpha_cov_sal,alpha_cov_herb, alpha_cov_pol))
    cova <- c(mean(cova$alpha_cov_sal), mean(cova$alpha_cov_herb), mean(cova$alpha_cov_pol))
    
    
}

#6.1 sacar parametros para el modelo global----
#estos datos son para meter todas las especies a la vez 
#lambda
lambda_sal <- do.call("rbind",lambda_sal)
lambda_sal <- data.frame(lambda_sal)
names(lambda_sal)<- c("lambda_sal")
rownames(lambda_sal) <- all.sp5

lambda_herb <- do.call("rbind",lambda_herb)
lambda_herb <- data.frame(lambda_herb)
names(lambda_herb)<- c("lambda_herb")
rownames(lambda_herb) <- all.sp5

lambda_pol <- do.call("rbind",lambda_pol)
lambda_pol <- data.frame(lambda_pol)
names(lambda_pol)<- c("lambda_pol")
rownames(lambda_pol) <- all.sp5

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

alpha_cov_sal <- do.call("rbind",alpha_cov_sal1)
alpha_cov_sal <- data.frame(alpha_cov_sal)
names(alpha_cov_sal)<- c("alpha_cov_sal")

alpha_cov_herb <- do.call("rbind",alpha_cov_herb1)
alpha_cov_herb <- data.frame(alpha_cov_herb)
names(alpha_cov_herb)<- c("alpha_cov_herb")

alpha_cov_pol <- do.call("rbind",alpha_cov_pol1)
alpha_cov_pol <- data.frame(alpha_cov_pol)
names(alpha_cov_pol)<- c("alpha_cov_pol")

alpha_cov_vector<- cbind(alpha_cov_sal,alpha_cov_herb,alpha_cov_pol)
mean_alpha_cov_sal <-mean(as.numeric(alpha_cov_vector$alpha_cov_sal))
mean_alpha_cov_herb <-mean(as.numeric(alpha_cov_vector$alpha_cov_herb))
mean_alpha_cov_pol <-mean(as.numeric(alpha_cov_vector$alpha_cov_pol))
alpha_cov_final<- c(mean_alpha_cov_sal,mean_alpha_cov_herb,mean_alpha_cov_pol)


names(alpha_cov_final) <- c("salt","herb1", "fv1")
names(lambda_cov_vector) <- c("salt","herb1", "fv1")

initial_values = list(lambda = lambdas_mean,
                      alpha_intra = intra_alpha_mean,
                      alpha_inter = inter_alpha_mean,
                      lambda_cov = lambda_cov_vector,
                      alpha_cov = alpha_cov_final)

# same with boundaries
lower_bounds = list(lambda = 10,
                    alpha_intra = 0.001,
                    alpha_inter = 0.001,
                    lambda_cov = -2,
                    alpha_cov = -2)

upper_bounds = list(lambda = 5000,
                    alpha_intra = 10,
                    alpha_inter = 10,
                    lambda_cov = 5,
                    alpha_cov = 5)

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




# 6.2 modelo global pero por especies ----
lower_bounds = list(lambda = 10,
                    alpha_intra = 0.001,
                    alpha_inter = 0.001,
                    lambda_cov = -2,
                    alpha_cov = -2)

upper_bounds = list(lambda = 5000,
                    alpha_intra = 10,
                    alpha_inter = 10,
                    lambda_cov = 5,
                    alpha_cov = 5)

bootstrap_samples <- 0
fixed_terms <- NULL

fit_global <- list()

for(i in 1:length(all.sp5)){
    obs<-obs6[i]
    sp <- all.sp5[i]
    env_list <-env_total[i]
    lambda_sp <- mean((mean(as.numeric(lambda_sal[i,]))),(mean(as.numeric(lambda_herb[i,]))), (mean(as.numeric(lambda_pol[i,]))))
    inter_sp <- mean((mean(as.numeric(inter_sal[i,]))),(mean(as.numeric(inter_herb[i,]))), (mean(as.numeric(inter_pol[i,]))))
    intra_sp <- mean((mean(as.numeric(intra_sal[i,]))),(mean(as.numeric(intra_herb[i,]))), (mean(as.numeric(intra_pol[i,]))))
    lambda_cov_sp <- as.numeric(lambda_cov[[i]])
    alpha_cov_sal_sp <- alpha_cov_sal1[[i]]
    alpha_cov_herb_sp <- alpha_cov_herb1 [[i]]
    alpha_cov_pol_sp <- alpha_cov_pol1 [[i]]
    alpha_cov_total <- c(alpha_cov_sal_sp,alpha_cov_herb_sp,alpha_cov_pol_sp )
    
    
    fit <- cxr_pm_multifit(data = obs,
                           focal_column = sp,
                           model_family = "BH",
                           covariates = env_list,
                           optimization_method = "bobyqa",
                           alpha_form = "pairwise",
                           lambda_cov_form = "global", # effect of covariates over lambda
                           alpha_cov_form = "global", # effect of covariates over alpha
                           initial_values = list(lambda = lambda_sp, alpha_intra =intra_sp, alpha_inter = inter_sp,
                                                 lambda_cov = lambda_cov_sp, alpha_cov = alpha_cov_total),
                           lower_bounds = lower_bounds,
                           upper_bounds = upper_bounds,
                           bootstrap_samples = 0)
    
    fit_global[[i]]<- fit
}

names(fit_global)<-all.sp5

#voy a seleccionar las especies que han convergido en el modelo: LEMA, HOMA, SOAS, CHFU, SCLA, SPRU, POMA, POMO, MESU

subset.species.global <- c("LEMA", "HOMA", "SOAS","CHFU", "SCLA","SPRU","POMA","POMO","MESU")

lema<- (as_tibble(fit_global$LEMA$alpha_matrix)) %>% select(subset.species.global)
homa <- (as_tibble(fit_global$HOMA$alpha_matrix))  %>% select(subset.species.global)
soas <- (as_tibble(fit_global$SOAS$alpha_matrix))  %>% select(subset.species.global)
chfu <- (as_tibble(fit_global$CHFU$alpha_matrix))  %>% select(subset.species.global)
scla <- (as_tibble(fit_global$SCLA$alpha_matrix))  %>% select(subset.species.global)
spru <- (as_tibble(fit_global$SPRU$alpha_matrix))  %>% select(subset.species.global)
poma <- (as_tibble(fit_global$POMA$alpha_matrix))  %>% select(subset.species.global)
pomo <- (as_tibble(fit_global$POMO$alpha_matrix))  %>% select(subset.species.global)
mesu <- (as_tibble(fit_global$MESU$alpha_matrix))  %>% select(subset.species.global)

sp_alpha_matrix.subset <- data.frame (rbind(lema, homa, soas, chfu, scla, spru, poma, pomo, mesu))
rownames(sp_alpha_matrix.subset) <- subset.species.global
#save(sp_alpha_matrix.subset, file="data/pairwise_alpha_matrix.Rda")

lema.lambda.cov<- (as_tibble(fit_global$LEMA$lambda_cov)) 
homa.lambda.cov <- (as_tibble(fit_global$HOMA$lambda_cov))  
soas.lambda.cov <- (as_tibble(fit_global$SOAS$lambda_cov))  
chfu.lambda.cov <- (as_tibble(fit_global$CHFU$lambda_cov))  
scla.lambda.cov <- (as_tibble(fit_global$SCLA$lambda_cov))  
spru.lambda.cov <- (as_tibble(fit_global$SPRU$lambda_cov))  
poma.lambda.cov <- (as_tibble(fit_global$POMA$lambda_cov))  
pomo.lambda.cov <- (as_tibble(fit_global$POMO$lambda_cov))  
mesu.lambda.cov <- (as_tibble(fit_global$MESU$lambda_cov))  
sp_lambda_cov.subset <- data.frame (rbind(lema.lambda.cov, homa.lambda.cov, soas.lambda.cov, chfu.lambda.cov, 
                                          scla.lambda.cov, spru.lambda.cov, poma.lambda.cov, pomo.lambda.cov, mesu.lambda.cov))
rownames(sp_lambda_cov.subset) <- subset.species.global

lema.lambda<- (as_tibble(fit_global$LEMA$lambda)) 
homa.lambda <- (as_tibble(fit_global$HOMA$lambda))  
soas.lambda <- (as_tibble(fit_global$SOAS$lambda))  
chfu.lambda <- (as_tibble(fit_global$CHFU$lambda))  
scla.lambda <- (as_tibble(fit_global$SCLA$lambda))  
spru.lambda <- (as_tibble(fit_global$SPRU$lambda))  
poma.lambda <- (as_tibble(fit_global$POMA$lambda))  
pomo.lambda <- (as_tibble(fit_global$POMO$lambda))  
mesu.lambda <- (as_tibble(fit_global$MESU$lambda))  
sp_lambda_.subset <- data.frame (rbind(lema.lambda, homa.lambda, soas.lambda, chfu.lambda, scla.lambda, spru.lambda,
                                       poma.lambda, pomo.lambda, mesu.lambda))
rownames(sp_lambda_.subset) <- subset.species.global
colnames(sp_lambda_.subset)<-"lambda"


lema.alpha_cov<- c((mean(as.numeric(fit_global$LEMA$alpha_cov$sal))), (mean(as.numeric(fit_global$LEMA$alpha_cov$herb1))),
                   (mean(as.numeric(fit_global$LEMA$alpha_cov$fv1))))
homa.alpha_cov <- c((mean(as.numeric(fit_global$HOMA$alpha_cov$sal))), (mean(as.numeric(fit_global$HOMA$alpha_cov$herb1))),
                    (mean(as.numeric(fit_global$HOMA$alpha_cov$fv1))))
soas.alpha_cov <- c((mean(as.numeric(fit_global$SOAS$alpha_cov$sal))), (mean(as.numeric(fit_global$SOAS$alpha_cov$herb1))),
                    (mean(as.numeric(fit_global$SOAS$alpha_cov$fv1)))) 
chfu.alpha_cov <- c((mean(as.numeric(fit_global$CHFU$alpha_cov$sal))), (mean(as.numeric(fit_global$CHFU$alpha_cov$herb1))),
                    (mean(as.numeric(fit_global$CHFU$alpha_cov$fv1))))
scla.alpha_cov <- c((mean(as.numeric(fit_global$SCLA$alpha_cov$sal))), (mean(as.numeric(fit_global$SCLA$alpha_cov$herb1))),
                    (mean(as.numeric(fit_global$SCLA$alpha_cov$fv1))))
spru.alpha_cov <- c((mean(as.numeric(fit_global$SPRU$alpha_cov$sal))), (mean(as.numeric(fit_global$SPRU$alpha_cov$herb1))),
                    (mean(as.numeric(fit_global$SPRU$alpha_cov$fv1)))) 
poma.alpha_cov <- c((mean(as.numeric(fit_global$POMA$alpha_cov$sal))), (mean(as.numeric(fit_global$POMA$alpha_cov$herb1))),
                    (mean(as.numeric(fit_global$POMA$alpha_cov$fv1))))
pomo.alpha_cov<- c((mean(as.numeric(fit_global$POMO$alpha_cov$sal))), (mean(as.numeric(fit_global$POMO$alpha_cov$herb1))),
                   (mean(as.numeric(fit_global$POMO$alpha_cov$fv1))))
mesu.alpha_cov <- c((mean(as.numeric(fit_global$MESU$alpha_cov$sal))), (mean(as.numeric(fit_global$MESU$alpha_cov$herb1))),
                    (mean(as.numeric(fit_global$MESU$alpha_cov$fv1))))  

sp_alpha_cov_subset <- data.frame (rbind(lema.alpha_cov, homa.alpha_cov, soas.alpha_cov, chfu.alpha_cov, scla.alpha_cov, spru.alpha_cov,
                                         poma.alpha_cov, pomo.alpha_cov, mesu.alpha_cov))
rownames(sp_alpha_cov_subset) <- subset.species.global
colnames(sp_alpha_cov_subset)<- c("salt", "herb", "pol")

x <- d1$neighbours
sp1 <- c("LEMA", "HOMA", "SOAS","CHFU", "SCLA","SPRU","POMA","POMO","MESU")
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

for(i in 1:length(sp1)){
    obs <- obs6[i]
    sp <- sp1[i]
    env_list <- env_total[i]
    env_list.herb [[i]] <- env_total[[i]][,"herb1"]
    env_list.pol [[i]] <- env_total[[i]][,"fv1"]
    caso.inter [[i]] <- mean(as.numeric(sp_alpha_matrix.subset[sp,which(colnames(sp_alpha_matrix.subset)!= sp)]))
    lambda_cov[[i]] <- (sp_lambda_cov.subset[sp,])
    lambda_cov_sal[[i]] <- (sp_lambda_cov.subset[sp,"sal"])
    lambda_cov_herb[[i]] <- (sp_lambda_cov.subset[sp,"herb1"])
    lambda_cov_fv [[i]] <- (sp_lambda_cov.subset[sp,"fv1"])
    alpha_cov.d [[i]] <- (sp_alpha_cov_subset[sp,])
    alpha_cov.sal [[i]]<- (sp_alpha_cov_subset[sp,"salt"])
    alpha_cov.herb[[i]] <- (sp_alpha_cov_subset[sp,"herb"])
    alpha_cov.pol[[i]] <- (sp_alpha_cov_subset[sp,"pol"])
    
    lambda_sp1 [[i]]<- (sp_lambda_.subset[sp,])
    
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

