setwd("~/Tesis/R_repositorios/Multitrophic_plants")
#libraries
library(tidyverse)

#load data
comp <- read_csv2("total_comp_19_20.csv") #son 2 años de campo

prueba1 <- gather(comp,"BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
                  ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE",
                  key = "Plant", value ="abundance" )

prueba1$ID <- paste(prueba1$focal, prueba1$Plant, sep="_")
numero.combinacion.vecinos <-  prueba1 %>% group_by(ID) %>% summarise(total = sum(abundance))%>%
    ungroup() #base de datos con las abundancias por pares de especies
numero.combinacion.vecinos$ID <- as.factor(numero.combinacion.vecinos$ID)
bb.mean <- numero.combinacion.vecinos %>% group_by(ID) %>% summarise(mean_neigh = mean(abundance))%>%
    ungroup()
bb.sd <- numero.combinacion.vecinos %>% group_by(ID) %>% summarise(sd_neigh = sd(abundance))%>%
    ungroup()
combi.neigh <- left_join(bb.mean, bb.sd) #base de datos con la media y desviacion estandar de las combinaciones
#       que ocurren en el campo 2019-2020