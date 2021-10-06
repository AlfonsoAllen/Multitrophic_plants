library(tidyverse)

load("data/alpha_matrix_final.Rda")
load("data/lambdas_final.Rda")
load("data/lambdas_cov_final.Rda")
load("data/alphas_cov_final.Rda")
alpha_matrix<- as.matrix(alpha.matrix)
alpha_cov <- alpha_cov
lambda <- lambdas
lambda_cov <- lambdas_cov
load("data/cov_total.Rda")
env <- total


comp<-read.csv("data/total_comp_19_20.check.csv", header=T, sep=";")

pairwise_alphas <- alpha_matrix[1:14,1:14]

species<- colnames(pairwise_alphas)

niche_diff<- matrix(NA, nrow= 14, ncol=14)

rownames(niche_diff)<-species

colnames(niche_diff)<-species



for( i in 1:14){
    
    for(j in 1:14){
        
        
        
        niche_diff[i,j] <- 1-(sqrt((pairwise_alphas[i,j]*pairwise_alphas[j,i])/(pairwise_alphas[i,i]*pairwise_alphas[j,j])))
        
    }
    
}

#fitness differences are a multiplicative interaction between demographic differences and competitive response differences

#First compute demographic differences


effects_all <- comp %>% group_by(focal)%>% summarise(fit = mean(fitness)) #effects_all lo he entendido como el fitness total


effects_all <- effects_all%>% filter(focal != "ACHI", focal != "ANAR")
    
fitness <- as.vector(effects_all[1:14,2])



fitness2<-replicate(14, fitness)
fitness2 <- do.call("rbind", fitness2)

dem_diff<- matrix(NA, nrow=14, ncol=14)

rownames(dem_diff)<-species

colnames(dem_diff)<-species



for( i in 1:14){
    
    for(j in 1:14){
        
        
        
        dem_diff[i,j] <- (fitness2[i,j]/fitness2[j,i])
        
    }
    
}


#Second compute competitive response differences

comp_res_diff<- matrix(NA, nrow= 14, ncol=14)

rownames(comp_res_diff)<-species

colnames(comp_res_diff)<-species



for( i in 1:14){
    
    for(j in 1:14){
        
        
        
        comp_res_diff[i,j] <- sqrt((pairwise_alphas[j,i]*pairwise_alphas[j,j])/(pairwise_alphas[i,i]*pairwise_alphas[i,j]))
        
    }
    
}


# Multiply demographic differences by competitive response differences to obtain fitness differences


fitness_diff<- comp_res_diff * dem_diff



# Ok, once calculated niche and fitness differences, plot them, 

#Remove diagonal

diag(niche_diff)=NA

diag(fitness_diff)=NA



plot(niche_diff, log(fitness_diff), pch=1, lwd=2, xlab="Niche differences", ylab="Fitness differences (Log. Transformed)", main= "a) Floral Visitors", ylim=c(-5,10), xlim=c(-60,1))

curve(log(1/(1-x)), add=T, col="red", lwd=3)

arrows(x0=0.9899, x1=1, y0=0, y1=10, length=0, lty=1, lwd=3, col="red")

text(x=-8.0,y=9.5, "Exclusion", font=3)

text(x=-4,y=-3, "Coexistence", font=3)



##Add the effect of pollinators on niche and fitness differences
###pollinators----

## build a matrix of interactions 

alphas_pol <- matrix(0, nrow=14,ncol=14)
row.names(alphas_pol)<- species
colnames(alphas_pol)<- species

for (i in 1:length(alpha_cov$pol)){ #this loop is to create a matrix with the effect of pol on changing per capita interactions
    alphas_pol[i,1:14]<-alpha_cov$pol[i]
}

alpha_matrix <-as.matrix(alpha_matrix)

alpha_matrix_pol <- list()

for (i in 1: max(total$fv1)){
    alpha_matrix_pol[[i]]<- alpha_matrix + alphas_pol * i
} #this is to create the 21 matrices 
     
#now calculate niche differences for each matrix 
niche_diff_pol <- list()
for (k in 1:length(alpha_matrix_pol)){
    
    xx<-alpha_matrix_pol[[k]]
    
    for( i in 1:14){
        
        for(j in 1:14){
            niche_diff[i,j] <- 1-(sqrt((xx[i,j]*xx[j,i])/(xx[i,i]*xx[j,j])))
        }
    }
    niche_diff_pol[[k]]<-niche_diff
}

#now fitness differences
#first the demographic differences
lambda_pol <- list()

for (i in 1: max(total$fv1)){
    lambda_pol[[i]]<- lambda*(1 + lambda_cov$pol * i)
}

demo_diff_pol <- list()
for (k in 1:length(lambda_pol)){
    xx<- as.data.frame(lambda_pol[[k]])
    xx2 <- t(replicate(length(species), xx$lambda))
    for( i in 1:14){
        for(j in 1:14){
            dem_diff[i,j] <- (xx2[i,j]/xx2[j,i])
        }
    }
    
    demo_diff_pol[[k]]<-dem_diff
}



#Second compute competitive response differences
comp_res_diff_pol <- list()
for (k in 1:length(alpha_matrix_pol)){
    
    xx<-alpha_matrix_pol[[k]]
    
    for( i in 1:14){
        
        for(j in 1:14){
            comp_res_diff[i,j] <- sqrt((xx[j,i]*xx[j,j])/(xx[i,i]*xx[i,j]))
        }
    }
    comp_res_diff_pol[[k]]<-comp_res_diff
}


fitness_diff_pol <- list()
for (k in 1:length(alpha_matrix_pol)){

    fitness_diff_pol[[k]]<- comp_res_diff_pol[[k]] * demo_diff_pol[[k]]
}

#now that we have the effect of polinators from 1 to 21 abundances (min and max), let see how the landscape changes
#select a species pair for an example for instance LEMA SPRU

example_niche<- list()
for (i in 1:length(niche_diff_pol)){
    example_niche[[i]] <- niche_diff_pol[[i]]["LEMA", "SPRU"] 
}
example_niche <- unlist(example_niche)

example_fitness<- list()
for (i in 1:length(fitness_diff_pol)){
    example_fitness[[i]] <- fitness_diff_pol[[i]]["LEMA", "SPRU"] 
}
example_fitness <- unlist(example_fitness)


#according to the data we see that pollinator in this case change competitive outcomes from one species to the another
#it can be better seen in this graph
boun_df<-data.frame(niche_overlap=c(seq(0,5, 0.02))) # creating a vector with niche overlap
boun_df$niche_diff<-(1-boun_df$niche_overlap) # calculating stabilizating differences from niche overlap 1-rho
boun_df$fitness_differences_sp_1<-(1/boun_df$niche_overlap) # coexistence line
boun_df$fitness_differences_sp_temp<- 1-boun_df$fitness_differences_sp_1 #this is an intermediate step to see the differences above one 
#which is later incorporated into the 2 species
boun_df$fitness_differences_sp_2<- 1+ boun_df$fitness_differences_sp_temp
boun_df<-boun_df[, -4]
#remove the intermediate step 
plot(example_niche, log(1/example_fitness),  pch=1, lwd=2, xlim=c(0, 1), ylim=c(0,max(log(1/example_fitness))),
     xlab="Niche differences", ylab="Fitness differences (Log. Transformed)",  main= "A) Floral Visitors")


points(boun_df$niche_diff, log(boun_df$fitness_differences_sp_1))
lines(boun_df$niche_diff, log(boun_df$fitness_differences_sp_1), type = "l", lty = 1, col="red")
lines(boun_df$niche_diff, log(boun_df$fitness_differences_sp_2), type = "l", lty = 1, col="blue")
abline(h=0)
text(x=-0.3, y=1, "Priority effect", cex=.8)
text(x=0.3, y=1, "Coexistence", cex=.8)
text(x=0.3, y=2.5, "SPRU excluded", cex=.8)
text(x=-0.2, y=0.1, "LEMA excluded", cex=.8)


#This can be done for other pair of species, the idea now is to look for a generalized measure that allows to capture these effects



##Until this point is ok 

##herbivores ----
## build a matrix of interactions 

alphas_herb <- matrix(0, nrow=14,ncol=14)
row.names(alphas_herb)<- species
colnames(alphas_herb)<- species

for (i in 1:length(alpha_cov$herb)){ #this loop is to create a matrix with the effect of herb on changing per capita interactions
    alphas_herb[i,1:14]<-alpha_cov$herb[i]
}

alpha_matrix <-as.matrix(alpha_matrix)

alpha_matrix_herb <- list()

for (i in 1: max(total$herb1)){
    alpha_matrix_herb[[i]]<- alpha_matrix + alphas_herb * i
} #this is to create the 21 matrices 

#now calculate niche differences for each matrix 
niche_diff_herb <- list()
for (k in 1:length(alpha_matrix_herb)){
    
    xx<-alpha_matrix_herb[[k]]
    
    for( i in 1:14){
        
        for(j in 1:14){
            niche_diff[i,j] <- 1-(sqrt((xx[i,j]*xx[j,i])/(xx[i,i]*xx[j,j])))
        }
    }
    niche_diff_herb[[k]]<-niche_diff
}

#now fitness differences
#first the demographic differences
lambda_herb <- list()

for (i in 1: max(total$herb1)){
    lambda_herb[[i]]<- lambda + lambda_cov$herb * i
}

demo_diff_herb <- list()
for (k in 1:length(lambda_herb)){
    xx<- as.data.frame(lambda_herb[[k]])
    xx2 <- t(replicate(length(species), xx$lambda))
    for( i in 1:14){
        for(j in 1:14){
            dem_diff[i,j] <- (xx2[i,j]/xx2[j,i])
        }
    }
    
    demo_diff_herb[[k]]<-dem_diff
}



#Second compute competitive response differences
comp_res_diff_herb <- list()
for (k in 1:length(alpha_matrix_herb)){
    
    xx<-alpha_matrix_herb[[k]]
    
    for( i in 1:14){
        
        for(j in 1:14){
            comp_res_diff[i,j] <- sqrt((xx[j,i]*xx[j,j])/(xx[i,i]*xx[i,j]))
        }
    }
    comp_res_diff_herb[[k]]<-comp_res_diff
}


fitness_diff_herb <- list()
for (k in 1:length(alpha_matrix_herb)){
    
    fitness_diff_herb[[k]]<- comp_res_diff_herb[[k]] * demo_diff_herb[[k]]
}

#now that we have the effect of herbivores from 1 to 54 abundances (min and max), let see how the landscape changes
#select a species pair for an example for instance LEMA HOMA

example_niche1<- list()
for (i in 1:length(niche_diff_herb)){
    example_niche1[[i]] <- niche_diff_herb[[i]]["LEMA", "HOMA"] 
}
example_niche1 <- unlist(example_niche1)

example_fitness1<- list()
for (i in 1:length(fitness_diff_herb)){
    example_fitness1[[i]] <- fitness_diff_herb[[i]]["LEMA", "HOMA"] 
}
example_fitness1 <- unlist(example_fitness1)


#according to the data we see that herbivores in this case change competitive outcomes from one species to the another
#it can be better seen in this graph
boun_df1<-data.frame(niche_overlap=c(seq(0,2, 0.05))) # creating a vector with niche overlap
boun_df1$niche_diff<-(1-boun_df1$niche_overlap) # calculating stabilizating differences from niche overlap 1-rho
boun_df1$fitness_differences_sp_1<-(1/boun_df1$niche_overlap) # solid line in your graph this is ok
boun_df1$fitness_differences_sp_temp<- 1-boun_df1$fitness_differences_sp_1 #this is an intermediate step to see the differences above one 
#which is later incorporated into the 2 species
boun_df1$fitness_differences_sp_2<- 1+ boun_df1$fitness_differences_sp_temp
boun_df1<-boun_df1[, -4]
#remove the intermediate step 
plot(example_niche1, log(example_fitness1), xlim=c(-0.8, 0.5), pch=1, lwd=2, xlab="Niche differences", ylab="Fitness differences (Log. Transformed)", main= "A) Floral Visitors")

points(boun_df1$niche_diff, boun_df1$fitness_differences_sp_1)
lines(boun_df1$niche_diff, boun_df1$fitness_differences_sp_1, type = "l", lty = 1, col="red")
lines(boun_df1$niche_diff, boun_df1$fitness_differences_sp_2, type = "l", lty = 1, col="blue")
text(x=-0.3, y=1, "Priority effect", cex=.8)
text(x=0.3, y=1, "Coexistence", cex=.8)
text(x=0.3, y=2.5, "HOMA excluded", cex=.8)
text(x=-0.2, y=0.1, "LEMA excluded", cex=.8)


#3. salinidad ----
## build a matrix of interactions 

alphas_salt <- matrix(0, nrow=14,ncol=14)
row.names(alphas_salt)<- species
colnames(alphas_salt)<- species

for (i in 1:length(alpha_cov$salt)){ #this loop is to create a matrix with the effect of salt on changing per capita interactions
    alphas_salt[i,1:14]<-alpha_cov$salt[i]
}

alpha_matrix <-as.matrix(alpha_matrix)

alpha_matrix_salt <- list()

for (i in 1: max(total$sal)){
    alpha_matrix_salt[[i]]<- alpha_matrix + alphas_salt * i
} #this is to create the 21 matrices 

#now calculate niche differences for each matrix 
niche_diff_salt <- list()
for (k in 1:length(alpha_matrix_salt)){
    
    xx<-alpha_matrix_salt[[k]]
    
    for( i in 1:14){
        
        for(j in 1:14){
            niche_diff[i,j] <- 1-(sqrt((xx[i,j]*xx[j,i])/(xx[i,i]*xx[j,j])))
        }
    }
    niche_diff_salt[[k]]<-niche_diff
}

#now fitness differences
#first the demographic differences
lambda_salt <- list()

for (i in 1: max(total$sal)){
    lambda_salt[[i]]<- lambda + lambda_cov$salt * i
}

demo_diff_salt <- list()
for (k in 1:length(lambda_salt)){
    xx<- as.data.frame(lambda_salt[[k]])
    xx2 <- t(replicate(length(species), xx$lambda))
    for( i in 1:14){
        for(j in 1:14){
            dem_diff[i,j] <- (xx2[i,j]/xx2[j,i])
        }
    }
    
    demo_diff_salt[[k]]<-dem_diff
}



#Second compute competitive response differences
comp_res_diff_salt <- list()
for (k in 1:length(alpha_matrix_salt)){
    
    xx<-alpha_matrix_salt[[k]]
    
    for( i in 1:14){
        
        for(j in 1:14){
            comp_res_diff[i,j] <- sqrt((xx[j,i]*xx[j,j])/(xx[i,i]*xx[i,j]))
        }
    }
    comp_res_diff_salt[[k]]<-comp_res_diff
}


fitness_diff_salt <- list()
for (k in 1:length(alpha_matrix_salt)){
    
    fitness_diff_salt[[k]]<- comp_res_diff_salt[[k]] * demo_diff_salt[[k]]
}

#now that we have the effect of salt from 0.0894 to 1.977 abundances (min and max), let see how the landscape changes
#select a species pair for an example for instance LEMA SPRU

example_niche2<- list()
for (i in 1:length(niche_diff_salt)){
    example_niche2[[i]] <- niche_diff_salt[[i]]["LEMA", "HOMA"] 
}
example_niche2 <- unlist(example_niche2)

example_fitness2<- list()
for (i in 1:length(fitness_diff_salt)){
    example_fitness2[[i]] <- fitness_diff_salt[[i]]["LEMA", "HOMA"] 
}
example_fitness2 <- unlist(example_fitness2)


#according to the data we see that pollinator in this case change competitive outcomes from one species to the another
#it can be better seen in this graph
boun_df2<-data.frame(niche_overlap=c(seq(0,2, 0.05))) # creating a vector with niche overlap
boun_df2$niche_diff<-(1-boun_df2$niche_overlap) # calculating stabilizating differences from niche overlap 1-rho
boun_df2$fitness_differences_sp_1<-(1/boun_df2$niche_overlap) # solid line in your graph this is ok
boun_df2$fitness_differences_sp_temp<- 1-boun_df2$fitness_differences_sp_1 #this is an intermediate step to see the differences above one 
#which is later incorporated into the 2 species
boun_df2$fitness_differences_sp_2<- 1+ boun_df2$fitness_differences_sp_temp
boun_df2<-boun_df2[, -4]
#remove the intermediate step 
plot(example_niche2, log(example_fitness2), xlim=c(-0.5, 0.5), pch=1, lwd=2, xlab="Niche differences", ylab="Fitness differences (Log. Transformed)", main= "A) Floral Visitors")

points(boun_df2$niche_diff, boun_df2$fitness_differences_sp_1)
lines(boun_df2$niche_diff, boun_df2$fitness_differences_sp_1, type = "l", lty = 1, col="red")
lines(boun_df2$niche_diff, boun_df2$fitness_differences_sp_2, type = "l", lty = 1, col="blue")
text(x=-0.3, y=1, "Priority effect", cex=.8)
text(x=0.3, y=1, "Coexistence", cex=.8)
text(x=0.3, y=2.5, "HOMA excluded", cex=.8)
text(x=-0.2, y=0.1, "LEMA excluded", cex=.8)

#coexistence metrics

#this example is with pollinators 
#changes in competitive dominance
comp_domin <-list()
change_comp_domin<-list()
xx<-rep("NA", times=21)

for(i in 1:length(fitness_diff_pol)){
    comp_domin[[i]] <-t(fitness_diff_pol[[i]])[lower.tri(t(fitness_diff_pol[[i]]))]
}

for(i in 1:91){ #this is because there are 91 sp. combinations 
    for(k in 1:length(fitness_diff_pol)){
        xx[k]<- comp_domin[[k]][i]
    }
    change_comp_domin[[i]]<-as.numeric(xx)
}        

for(i in 1:length(change_comp_domin)){
    if(change_comp_domin[[i]]>0 | change_comp_domin[[i]]<0) == TRUE{}
}

    
    


