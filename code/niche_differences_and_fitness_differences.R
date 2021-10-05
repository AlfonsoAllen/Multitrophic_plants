load("C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/data/pairwise_alpha_matrix.Rda")
load("C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/data/lambda_cov_matrix.Rda")
load("C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/data/lambdas.Rda")
load("C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/data/alpha_cov_matrix.Rda")
alpha_matrix<- sp_alpha_matrix.subset
alpha_cov <- sp_alpha_cov_subset
lambda <- sp_lambda_.subset
lambda_cov <- sp_lambda_cov.subset
load("C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/data/cov_total.Rda")
env <- total


comp<-read.csv("data/total_comp_19_20.check.csv", header=T, sep=";")

pairwise_alphas <- alpha_matrix[1:9,1:9]

species<- c("LEMA","HOMA","SOAS","CHFU","SCLA","SPRU","POMA","POMO","MESU")

niche_diff<- matrix(NA, nrow= 9, ncol=9)

rownames(niche_diff)<-species

colnames(niche_diff)<-species



for( i in 1:9){
    
    for(j in 1:9){
        
        
        
        niche_diff[i,j] <- 1-(sqrt((pairwise_alphas[i,j]*pairwise_alphas[j,i])/(pairwise_alphas[i,i]*pairwise_alphas[j,j])))
        
    }
    
}



#fitness differences are a multiplicative interaction between demographic differences and competitive response differences



#First compute demographic differences


effects_all <- comp %>% group_by(focal)%>% summarise(fit = mean(fitness)) #effects_all lo he entendido como el fitness total


effects_all <- effects_all%>% filter(focal != "ACHI", focal != "ANAR", focal != "BEMA", focal != "CETE", focal!= "PAIN",focal != "PLCO",
                     focal != "PUPA" )
    


fitness <- as.vector(effects_all[1:9,2])



fitness2<-replicate(9, fitness)
fitness2 <- do.call("rbind", fitness2)

dem_diff<- matrix(NA, nrow=9, ncol=9)

rownames(dem_diff)<-species

colnames(dem_diff)<-species



for( i in 1:9){
    
    for(j in 1:9){
        
        
        
        dem_diff[i,j] <- (fitness2[i,j]/fitness2[j,i])
        
    }
    
}





#Second compute competitive response differences



comp_res_diff<- matrix(NA, nrow= 9, ncol=9)

rownames(comp_res_diff)<-species

colnames(comp_res_diff)<-species



for( i in 1:9){
    
    for(j in 1:9){
        
        
        
        comp_res_diff[i,j] <- sqrt((pairwise_alphas[j,i]*pairwise_alphas[j,j])/(pairwise_alphas[i,i]*pairwise_alphas[i,j]))
        
    }
    
}





# Multiply demographic differences by competitive response differences to obtain fitness differences



fitness_diff<- comp_res_diff * dem_diff



## Values lower than 1 for fitness differences need to be removed



fitness_diff[fitness_diff<1]<-NA

# Ok, once calculated niche and fitness differences, plot them, 

#Remove diagonal

diag(niche_diff)=NA

diag(fitness_diff)=NA



plot(niche_diff, log(fitness_diff), pch=1, lwd=2, xlab="Niche differences", ylab="Fitness differences (Log. Transformed)", main= "a) Floral Visitors", ylim=c(0,10), xlim=c(0,1))

curve(log(1/(1-x)), add=T, col="red", lwd=3)

arrows(x0=0.9899, x1=1, y0=4.6, y1=10, length=0, lty=1, lwd=3, col="red")

text(x=0.8,y=9.5, "Exclusion", font=3)

text(x=0.8,y=0.2, "Coexistence", font=3)



##Add the effect of pollinators on niche and fitness differences----



## build a matrix of interactions 

pairwise_alphas_pol <- matrix(0, nrow=9,ncol=9)

row.names(pairwise_alphas_pol)<- species

colnames(pairwise_alphas_pol)<- species

alpha_matrix <-as.matrix(alpha_matrix)

#effects_all <-as.matrix(effects_all[1:14,1:14]) #aqui me pierdo

pairwise_alphas_pol[1,]<-alpha_matrix[1,1:9] + alpha_cov[1,3]

pairwise_alphas_pol[2,]<-alpha_matrix[2,1:9] + alpha_cov[2,3]

pairwise_alphas_pol[3,]<-alpha_matrix[3,1:9] + alpha_cov[3,3]

pairwise_alphas_pol[4,]<-alpha_matrix[4,1:9] + alpha_cov[4,3]

pairwise_alphas_pol[5,]<-alpha_matrix[5,1:9] + alpha_cov[5,3]

pairwise_alphas_pol[6,]<-alpha_matrix[6,1:9] + alpha_cov[6,3]

pairwise_alphas_pol[7,]<-alpha_matrix[7,1:9] + alpha_cov[7,3]

pairwise_alphas_pol[8,]<-alpha_matrix[8,1:9] + alpha_cov[8,3]

pairwise_alphas_pol[9,]<-alpha_matrix[9,1:9] + alpha_cov[9,3]





niche_diff_pol<- matrix(NA, nrow= 9, ncol=9)

rownames(niche_diff_pol)<-species

colnames(niche_diff_pol)<-species



for( i in 1:9){
    
    for(j in 1:9){
        
        
        
        niche_diff_pol[i,j] <- 1-(sqrt((pairwise_alphas_pol[i,j]*pairwise_alphas_pol[j,i])/(pairwise_alphas_pol[i,i]*pairwise_alphas_pol[j,j])))
        
    }
    
}



#fitness differences are a multiplicative interaction between demographic differences and competitive response differences



#First compute demographic differences


fitness_pol <- matrix(0, nrow=9, ncol=1)

row.names(fitness_pol)<- species

colnames(fitness_pol)<- "lambda"

fitness_pol[1,]<-lambda[1,1] + lambda_cov[1,1]

fitness_pol[2,]<-lambda[2,1] + lambda_cov[2,1]

fitness_pol[3,]<-lambda[3,1] + lambda_cov[3,1]

fitness_pol[4,]<-lambda[4,1] + lambda_cov[4,1]

fitness_pol[5,]<-lambda[5,1] + lambda_cov[5,1]

fitness_pol[6,]<-lambda[6,1] + lambda_cov[6,1]

fitness_pol[7,]<-lambda[7,1] + lambda_cov[7,1]

fitness_pol[8,]<-lambda[8,1] + lambda_cov[8,1]

fitness_pol[9,]<-lambda[9,1] + lambda_cov[9,1]

pol <- seq(0,21, by=1)
fitness_pol <- as.vector(fitness_pol)

fitness_pol.visits <- list()
fitness2_pol1 <- list()


for (i in 1: length(pol)){
num <- pol[i]
fitness_pol.visits [[i]] <- (fitness_pol * num)
fitness2_pol1 [[i]] <-replicate(9, fitness_pol.visits [[i]])

    for( k in 1:9){
    
         for(j in 1:9){
    
            
            dem_diff_pol[[i]][k,j]<- as.matrix((fitness2_pol1[[i]][k,j])/(fitness2_pol1[[i]][j,k]))
          
    }
    
}


}

dem_diff_pol<- matrix(NA, nrow= 9, ncol=9)

rownames(dem_diff_pol)<-species

colnames(dem_diff_pol)<-species




#Second compute competitive response differences

comp_res_diff_pol<- matrix(NA, nrow= 9, ncol=9)

rownames(comp_res_diff_pol)<-species

colnames(comp_res_diff_pol)<-species

comp_res_diff_pol <- list()

for( i in 1:9){
    
    for(j in 1:9){
        for( k in 1:length(pol)){
        
        
        comp_res_diff_pol[i,j] <- as.matrix(sqrt((pairwise_alphas_pol[j,i]*pairwise_alphas_pol[j,j]))*k/
                                                ((pairwise_alphas_pol[i,i]*pairwise_alphas_pol[i,j])*k))
        }    
    }
    
}



# Multiply demographic differences by competitive response differences to obtain fitness differences

fitness_diff_pol<- comp_res_diff_pol * dem_diff_pol



fitness_diff_pol[fitness_diff_pol<1]<-NA

# Ok, once calculated niche and fitness differences, plot them, 

diag(niche_diff_pol)=NA

diag(fitness_diff_pol)=NA



points(niche_diff_pol, log(fitness_diff_pol), pch=16, lwd=1, col="black")

legend("topleft", c("No floral visitors", "Floral visitors"), col = c(1, 1), lwd=c(1,1), lty=c(0,0), 
       
       pch=c(1,16), bty="n",
       
       merge = TRUE)



#↑esta ultima parte mirar detenidamente, ya que segun los tratamientos que ponga y las species seran una cosa u otra

#Use arrows to conect dots between treatments. 

fitness_diff[6,1]<-fitness_diff[1,6]

fitness_diff[1,6]<-NA

fitness_diff[3,2]<-fitness_diff[2,3]

fitness_diff[2,3]<-NA

fitness_diff[5,3]<-fitness_diff[3,5]

fitness_diff[3,5]<-NA

fitness_diff[6,4]<-fitness_diff[4,6]

fitness_diff[4,6]<-NA





fitness_diff_pol[5,3]<-fitness_diff_pol[3,5]

fitness_diff_pol[3,5]<-NA

fitness_diff_pol[2,6]<-fitness_diff_pol[6,2]

fitness_diff_pol[6,2]<-NA



arrows(x0=niche_diff, x1=niche_diff_pol, y0=log(fitness_diff), y1=log(fitness_diff_pol), length=0.10, lty=2)
