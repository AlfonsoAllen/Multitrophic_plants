pairwise_alphas <- alpha_matrix[1:14,1:14]

species<- c("LEMA","BEMA","HOMA","SOAS","CHFU","SCLA","SPRU","POMA","POMO","PAIN","CETE","PLCO","PUPA","MESU")

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



fitness <- as.vector(effects_all[1:14,1])



fitness2<-replicate(14, fitness)

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

pairwise_alphas_pol <- matrix(0, nrow=14,ncol=14)

row.names(pairwise_alphas_pol)<- species

colnames(pairwise_alphas_pol)<- species

alpha_matrix <-as.matrix(alpha_matrix)

effects_all <-as.matrix(effects_all[1:14,1:14])

pairwise_alphas_pol[1,]<-alpha_matrix[1,1:14] + effects_all[1,5]

pairwise_alphas_pol[2,]<-alpha_matrix[2,1:14] + effects_all[2,5]

pairwise_alphas_pol[3,]<-alpha_matrix[3,1:14] + effects_all[3,5]

pairwise_alphas_pol[4,]<-alpha_matrix[4,1:14] + effects_all[4,5]

pairwise_alphas_pol[5,]<-alpha_matrix[5,1:14] + effects_all[5,5]

pairwise_alphas_pol[6,]<-alpha_matrix[6,1:14] + effects_all[6,5]

pairwise_alphas_pol[7,]<-alpha_matrix[7,1:14] + effects_all[7,5]

pairwise_alphas_pol[8,]<-alpha_matrix[8,1:14] + effects_all[8,5]

pairwise_alphas_pol[9,]<-alpha_matrix[9,1:14] + effects_all[9,5]

pairwise_alphas_pol[10,]<-alpha_matrix[10,1:14] + effects_all[10,5]

pairwise_alphas_pol[11,]<-alpha_matrix[11,1:14] + effects_all[11,5]

pairwise_alphas_pol[12,]<-alpha_matrix[12,1:14] + effects_all[12,5]

pairwise_alphas_pol[13,]<-alpha_matrix[13,1:14] + effects_all[13,5]

pairwise_alphas_pol[14,]<-alpha_matrix[14,1:14] + effects_all[14,5]



niche_diff_pol<- matrix(NA, nrow= 14, ncol=14)

rownames(niche_diff_pol)<-species

colnames(niche_diff_pol)<-species



for( i in 1:14){
    
    for(j in 1:14){
        
        
        
        niche_diff_pol[i,j] <- 1-(sqrt((pairwise_alphas_pol[i,j]*pairwise_alphas_pol[j,i])/(pairwise_alphas_pol[i,i]*pairwise_alphas_pol[j,j])))
        
    }
    
}



#fitness differences are a multiplicative interaction between demographic differences and competitive response differences



#First compute demographic differences



fitness_pol <- matrix(0, nrow=14, ncol=1)

row.names(fitness_pol)<- species

colnames(fitness_pol)<- "lambda"

fitness_pol[1,]<-effects_all[1,1] + effects_all[1,2]

fitness_pol[2,]<-effects_all[2,1] + effects_all[2,2]

fitness_pol[3,]<-effects_all[3,1] + effects_all[3,2]

fitness_pol[4,]<-effects_all[4,1] + effects_all[4,2]

fitness_pol[5,]<-effects_all[5,1] + effects_all[5,2]

fitness_pol[6,]<-effects_all[6,1] + effects_all[6,2]

fitness_pol[7,]<-effects_all[7,1] + effects_all[7,2]

fitness_pol[8,]<-effects_all[8,1] + effects_all[8,2]

fitness_pol[9,]<-effects_all[9,1] + effects_all[9,2]

fitness_pol[10,]<-effects_all[10,1] + effects_all[10,2]

fitness_pol[11,]<-effects_all[11,1] + effects_all[11,2]

fitness_pol[12,]<-effects_all[12,1] + effects_all[12,2]

fitness_pol[13,]<-effects_all[13,1] + effects_all[13,2]

fitness_pol[14,]<-effects_all[14,1] + effects_all[14,2]




fitness_pol <- as.vector(fitness_pol)



fitness2_pol<-replicate(6, fitness_pol)

dem_diff_pol<- matrix(NA, nrow= 14, ncol=14)

rownames(dem_diff_pol)<-species

colnames(dem_diff_pol)<-species



for( i in 1:14){
    
    for(j in 1:14){
        
        
        
        dem_diff_pol[i,j] <- (fitness2_pol[i,j]/fitness2_pol[j,i])
        
    }
    
}





#Second compute competitive response differences



comp_res_diff_pol<- matrix(NA, nrow= 14, ncol=14)

rownames(comp_res_diff_pol)<-species

colnames(comp_res_diff_pol)<-species



for( i in 1:14){
    
    for(j in 1:14){
        
        
        
        comp_res_diff_pol[i,j] <- sqrt((pairwise_alphas_pol[j,i]*pairwise_alphas_pol[j,j])/(pairwise_alphas_pol[i,i]*pairwise_alphas_pol[i,j]))
        
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