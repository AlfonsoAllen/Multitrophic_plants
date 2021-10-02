setwd("~/Tesis/R_repositorios/Multitrophic_plants")
#this script if for cleaning the data that i'm going to use for the multitrophic plants
library(tidyverse)
library(cxr)
library(ggplot2)

#load the data
#2020
fv.20 <- read_csv2("data/raw_data/pollinator_2020.csv")
h.20 <-read.table(file = "data/raw_data/foodweb_2020.csv", header = TRUE, sep = ";")
comp.20 <- read_csv2("data/raw_data/competition_2020.csv")
salt.20 <- read_csv2("data/raw_data/Salinity_moisture_2020.csv")
#2019
fv<-read.table(file = "data/raw_data/pollinator_2016_2019.csv", header = TRUE, sep = ";")# i need to keep only 2019
h.19<-read.table(file = "data/raw_data/foodweb_2019.csv", header = TRUE, sep = ";")
comp.19<-read.table(file = "data/raw_data/competition_2019.csv", header = TRUE, sep = ";")
salt.19 <- read.table(file = "data/raw_data/Salinity_moisture_2019.csv", header = TRUE, sep = ";")
#position plots
position <- read.table(file="data/raw_data/caracolesplotposition.csv", header= TRUE, sep=";")
position <- na.omit(position)
position$Plot <- position$plot
position$Subplot <- position$position
position <- position[,c("Plot","Subplot")]

#1. clean the data ----
# floral visitors and herbivores. ----
#Keep only the 2019 and 2020 years. 
#fv 2020
ggplot(data=fv.20, aes(x= Plant, y=number_visits, fill= Group)) + 
    geom_boxplot(aes(x = Plant, y = number_visits)) +
    ggtitle("visitor abundance per Plant")+
    # facet_grid(Group~.)+
    NULL
fv.20 <-fv.20 %>% filter(Group!= "Wasp", Group!= "0", Plant!= "NA", Plant!="OUT",ID_Simple!= "Cassida_sp.",
                         ID_Simple!="Bee_nest", ID_Simple!= "0",Subplot!= "OUT",
                         ID_Simple!= "Elateridae",Plant!= "MEPO", Plant!= "ANAR", 
                         Plant!= "CHMI", Plant!="RAPE", Plant!= "SASO")
fv.20.simple <- fv.20 %>% group_by(Plot, Subplot, Plant) %>% summarise (fv = sum(number_visits))%>%
    ungroup()
fv.20.group <- fv.20 %>% group_by(Plot, Subplot, Plant, Group) %>% summarise (fv = sum(number_visits))%>%
    ungroup()

fv.20.simple$Focal <-fv.20.simple$Plant
#fv 2019
fv.19 <- subset(fv, Year == "2019")
fv.19 <- fv.19 %>% filter(Group!= "NA", Group!= "Mantodea", Group!= "Larva", Group!= "Bug",
                          Group!= "Grasshopper",Plant!="0", Plant!= "HOMA",
                          Plant!= "NA", Plant!="OUT", Plant!= "Lysimachia_arvensis",
                          ID_Simple!= "Cassida_sp.", ID_Simple!="Bee_nest", ID_Simple!= "0",
                          ID_Simple!= "Elateridae", ID_Simple!= "Coccinella_septempunctata", 
                          ID_Simple!= "Larva", Subplot!="OUT", Plant!= "MEPO", Plant!= "ANAR", 
                          Plant!= "CHMI", Plant!="RAPE", Plant!= "SASO")
fv.19$number_visits <-as.numeric(fv.19$number_visits)
fv.19 <- fv.19 %>% group_by(Plot, Subplot, Plant, Group) %>% summarise (fv = sum(number_visits))%>%
    ungroup()
ggplot(data=fv.19, aes(x= Plant, y=fv, fill= Group)) + 
    geom_boxplot(aes(x = Plant, y = fv)) +
    ggtitle("visitor abundance per Plant")+
    # facet_grid(Group~.)+
    NULL #clean: grasshopper, mantodea, bug, larva, NA



#herbivores 2020
ggplot(data=h.20, aes(x= Plant_simple, y=number_animal, fill= Group)) + 
    geom_boxplot(aes(x = Plant_simple, y = number_animal)) +
    ggtitle("herb abundance per Plant")+
    # facet_grid(Group~.)+
    NULL
h1.20 <- h.20 %>% group_by(Plot, Subplot, Plant_simple) %>% summarise (herb = sum(number_animal))%>%
    ungroup()
h1.20$Focal <- h1.20$Plant_simple
h1.20 <- h1.20[,c("Plot","Subplot","Focal", "herb")]


h.20.filter <- h.20 %>% filter(!is.na(Group),!is.na(Plot), !is.na(Subplot), Group!="Dragonfly",Group!="Mantis",Group!="Butterfly", Group!="Neuroptera",
                         Group!="Spider", Group!="NA", Plant_simple!= "AIR", Plant_simple!="Air", Plant_simple!= "Ground",
                         Plant_simple!= "Stick", Plant_simple!= "Salicornia", Plant_simple!= "NA", Group!="Fly", 
                         Plant_simple!= "New_grass", Plant_simple!="ground", Plant_simple!="Junco", Plant_simple!="Artemisia",
                         Plant_simple!="ACHI", Plant_simple!= "salicornia", Plant_simple!= "New_Grass", 
                         raw_species_id!= "Brassicogethes", raw_species_id!= "Psilothrix_viridicoerulea",
                         raw_species_id!= "Cerambycidae", raw_species_id!="Cryptocephalus_sp", 
                         raw_species_id!="Curculionidae", raw_species_id!= "Malachius", 
                         raw_species_id!= "Malachius_bipustulatus", raw_species_id!= "Melyridae", 
                         raw_species_id!="Mordellidae", raw_species_id!= "Oedemeridae",
                         raw_species_id!= "Lagorina_sericea", raw_species_id!= "Lepidoptera_eggs",
                         raw_species_id!= "lepidoptera_coccon",Plant!= "MEPO", Plant!= "ANAR", 
                         Plant!= "CHMI", Plant!="RAPE", Plant!= "SASO"
                         )
h.20.group <- h.20.filter%>% group_by(Plot, Subplot, Plant_simple, Group) %>% summarise (herb = sum(number_animal))%>%
    ungroup()
ggplot(data=h.20.group, aes(x= Plant_simple, y=herb, fill= Group)) + 
    geom_boxplot(aes(x = Plant_simple, y = herb)) +
    ggtitle("herb abundance per Plant 2020")+
    # facet_grid(Group~.)+
    NULL

#herb 2019
h1.19 <- h.19 %>% group_by(Plot, Subplot, Plant_Simple) %>% summarise (herb = sum(number_animal))%>%
    ungroup()
h.19.group <- h.19%>% group_by(Plot, Subplot, Plant_Simple, Group) %>% summarise (herb = sum(number_animal))%>%
    ungroup()
ggplot(data=h.19.group, aes(x= Plant_Simple, y=herb, fill= Group)) + 
    geom_boxplot(aes(x = Plant_Simple, y = herb)) +
    ggtitle("herb abundance per Plant")+
    # facet_grid(Group~.)+
    NULL
h.19.filter <- h.19 %>% filter(!is.na(Group),!is.na(Plot), !is.na(Subplot), Group!="NA", Plant_Simple!= "NA",
                      Plant_Simple!= "LYTR", Plant_Simple!= "Rush", Plant_Simple!= "Salicornia",
                      Plant_Simple!= "Stick",raw_species_id!= "Brassicogethes", raw_species_id!= "Psilothrix_viridicoerulea",
                      raw_species_id!= "Cerambycidae", raw_species_id!="Cryptocephalus_sp", 
                      raw_species_id!="Curculionidae", raw_species_id!= "Malachius", 
                      raw_species_id!= "Malachius_bipustulatus", raw_species_id!= "Melyridae", 
                      raw_species_id!="Mordellidae", raw_species_id!= "Oedemeridae",
                      raw_species_id!= "Lagorina_sericea",Plant!= "MEPO", Plant!= "ANAR", 
                      Plant!= "CHMI", Plant!="RAPE", Plant!= "SASO")
h.19.group.1 <- h.19.filter%>% group_by(Plot, Subplot, Plant_Simple, Group) %>% summarise (herb = sum(number_animal))%>%
    ungroup()

ggplot(data=h.19.group.1, aes(x= Plant_Simple, y=herb, fill= Group)) + 
    geom_boxplot(aes(x = Plant_Simple, y = herb)) +
    ggtitle("herb abundance per Plant 2019")+
    # facet_grid(Group~.)+
    NULL
#salt ----
#2020
salt.20$Salinity <- as.numeric(salt.20$Salinity)
sal.20.simple <- salt.20 %>% group_by(Plot, Subplot) %>% summarise(sal = sum(Salinity))%>%
    ungroup()
any(is.na(sal.20.simple))
which(is.na(sal.20.simple$sal))
str(sal.20.simple)
sal.prueba2 <- sal.20.simple
sal.prueba2[78,3]<-0.1231
sal.prueba2[141,3]<-0.0894
sal.prueba2[187,3]<-0.1069
sal.prueba2[322,3]<- 0.122
any(is.na(sal.prueba2))
max(sal.prueba2$sal)
#2019
str(salt.19) #There are a lot of Nas in the salinity, so I'm going to calculate the average per Plot
nas <- salt.19
any(is.na(nas$Salinity))
which(is.na(nas$Salinity))
nas <- nas %>% filter(!is.na(Salinity))

s <- nas %>% group_by(Plot) %>% summarise(sal = mean(Salinity ))%>%
    ungroup()
#

sal.19.simple <- salt.19 %>% group_by(Plot, Subplot) %>% summarise(sal = sum(Salinity_mean))%>%
    ungroup()
which(is.na(sal.19.simple$sal))#ya tengo la base de datos sin NAs
max(sal.19.simple$sal)
#competition ----
#el fitness tiene que ser numero de semillas en 1 fruto o numero de semillas por indv?----
#creo que lo correcto es numero de frutos por semillas que hay dentro de 1 fruto.
#comp.20.completa <- left_join(position, comp.20)#asi me aseguro de que esté todos los plots
comp.20 <- comp.20[,c("Plot", "Subplot", "Focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
                  ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]
comp.20.1 <- comp.20 %>% filter(Focal!= "RAPE", Focal!="MEPO", Focal!="COSQ" , Focal!= "SUSP", 
                          Focal!="FRPU",Focal!="MEEL", Focal!= "CHMI", Focal!= "SASO")
comp.20.1$focal <- comp.20.1$Focal
comp.20.1.join <- comp.20.1[,c("Plot", "Subplot", "focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
                             ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]
comp.20.1<- comp.20.1[,c("focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
           ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]

comp.20.1$check1 <- rowSums(comp.20.1[,3:24])
min(comp.20.1$check1) #hay filas enteras con 0 vecinos. Hya que eliminarlos. ya que el modelo no funciona sino
#comp.20.final <- comp.20.1 %>% filter(check1!= "0")
comp.20.final <- subset(comp.20.1, check1>=1)
comp.20.final <- comp.20.final[,c("focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
             ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]

comp.20.1.join$check1 <- rowSums(comp.20.1.join[,5:26])
min(comp.20.1.join$check1) #hay filas enteras con 0 vecinos. Hya que eliminarlos. 
#comp.20.1.join <- comp.20.1.join %>% filter(check1!= "0")
comp.20.final.junto <- subset(comp.20.1.join, check1>=1)
comp.20.final.junto <- comp.20.final.junto[,c("Plot","Subplot","focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
                         ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]


#position$plot <- position$Plot
#position$subplot <- position$Subplot
#comp.19.completa <- left_join(position, comp.19)#asi me aseguro de que esté todos los plots
comp.19 <- comp.19[,c("plot", "subplot" ,"focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
                      ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]
comp.19.1 <-comp.19 %>% filter(focal!= "RAPE", focal!="MEPO", focal!="COSQ" , 
                                 focal!= "SUSP", focal!="FRPU",focal!="MEEL", focal!= "CHMI", focal!= "SASO")
comp.19.1.join <- comp.19.1[,c("plot", "subplot" ,"focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
             ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]
comp.19.1 <- comp.19.1[,c("focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
           ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]
comp.19.1$check1 <- rowSums(comp.19.1[,3:24]) #hay filas con 0 vecinos, hay que eliminarlos
min(comp.19.1$check1)
#comp.19.final <- comp.19.1 %>% filter(check1!= "0")
comp.19.final <- subset(comp.19.1, check1>=1)
comp.19.final <- comp.19.final[,c("focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
                          ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]

comp.19.1.join$check1 <- rowSums(comp.19.1.join[,5:26]) #hay filas con 0 vecinos, hay que eliminarlos
min(comp.19.1.join$check1)
#comp.19.1.join <- comp.19.1.join %>% filter(check1!= "0")
comp.19.final.junto <- subset(comp.19.1.join, check1>=1)
comp.19.final.junto <- comp.19.final.junto[,c("plot", "subplot","focal", "fitness", "BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
                     ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE")]


#join the competition data
c.20 <- colnames(comp.20.final)
c.19 <- colnames(comp.19.final)
common_names <- intersect(c.20, c.19)
total_comp.final <- rbind(comp.20.final[common_names], comp.19.final[common_names])
ghj<- rowSums(total_comp.final [,3:24]) #comprobacion de que no haya 0 vecinos para una focal
min(ghj)
c.20.1 <- colnames(comp.20.final.junto)
comp.20.final.junto$plot <- comp.20.final.junto$Plot
comp.20.final.junto$subplot <- comp.20.final.junto$Subplot
c.20.1 <- colnames(comp.20.final.junto)
c.19.1 <- colnames(comp.19.final.junto)
common_names1 <- intersect(c.20.1, c.19.1)


total_comp_withplots.final <- rbind(comp.20.final.junto[common_names1], comp.19.final.junto[common_names1])

ghj1<- rowSums(total_comp_withplots.final [,3:24])
min(ghj1)
#write.csv2(total_comp.final, file ="C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/total_comp_19_20.final.csv", row.names = FALSE)
#write.csv2(total_comp_withplots.final, file ="C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/total_comp_19_20_withplot_subplot.final.csv", row.names = FALSE)


#exploratory analyses
#competencia: filtrar por especies. una vez lo tenga hecho por especies tengo que sumar los vecinos intra+ inter
#y plotear numero de vecinos segun el fitness
comp.20.junto <- gather(comp.20.final.junto,"BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
                        ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE",
                        key = "Plant", value ="abundance" )
comp.20.junto <- comp.20.junto %>% filter(focal!= "RAPE", focal!="MEPO", focal!="COSQ" , focal!= "SUSP", 
                                          focal!= "SUSP", focal!="FRPU")

comp.19.junto <- gather(comp.19.final.junto,"BEMA","CHFU","CHMI","CETE","HOMA","LEMA","PAIN", "PLCO","POMA","POMO", "PUPA","SASO","SCLA"
                        ,"SPRU","MESU","MEEL","MEPO","SOAS", "FRPU", "SUSP","COSQ", "RAPE",
                        key = "Plant", value ="abundance" )
comp.19.junto <- comp.19.junto %>% filter(focal!= "RAPE", focal!="MEPO", focal!="COSQ" , focal!= "SUSP", 
                                           focal!="FRPU", focal!="MEEL", focal!= "CHMI", focal!= "SASO")

# the species that I'm going to select for the study are: HOMA, POMO, POMA, LEMA, CHFU, PUPA, SOAS, SCLA, MESU, SPRU, BEMA, CETE
    #PAIN,  and PLCO. 14 species in total. 


#aqui podría crear una columna donde pusiera vecinos intra y vecinos inter. esto tiene que
#ser a partir de un condicionante if, y de alguna funcion mutate
neigh.20 <- comp.20.junto %>% group_by(Plot, Subplot, focal, fitness) %>% summarise (neigh = sum(abundance))%>%
    ungroup()
neigh.20 <- subset(neigh.20, focal %in% c("SOAS","CHFU","LEMA", "SCLA","BEMA","CETE","MESU","PUPA", "SPRU", "HOMA","POMA","POMO",
                                         "PAIN", "PLCO" ))

neigh.19 <- comp.19.junto %>% group_by(plot, subplot, focal, fitness) %>% summarise (neigh = sum(abundance))%>%
    ungroup()
neigh.19<- subset(neigh.19, focal %in% c("SOAS","CHFU","LEMA", "SCLA","BEMA","CETE","MESU","PUPA", "SPRU", "HOMA", "POMA","POMO",
                                         "PAIN", "PLCO"))

ggplot(neigh.19, aes(x = neigh, y = fitness, fill= focal))+
    geom_point()+
    geom_smooth(method = "lm")+
    xlab("Abundancia de plantas")+
    ylab("Numero de visitas de polinizadores")+
    ggtitle ("relacion entre el numero de abundancias de plantas y el fitness 2019")+
    NULL


ggplot(neigh.19, aes(x = neigh, y = fitness, group = focal))+
    geom_point(aes(color = focal))+
    geom_smooth(method = "lm", aes(color = focal))+
    ylab(" viable seed")+
    xlab("Number of neighbors")+
    ylim(1, 2000)+
    xlim(1, 70)+
    ggtitle("Fitness in relation with the abundance of neighbors 2019")+
    NULL
#nota, podría distinguir vecinos intra e inter


ggplot(neigh.20, aes(x = neigh, y = log(fitness), group = focal))+
    geom_point(aes(color = focal))+
    geom_smooth(method = "lm", aes(color = focal))+
    ylab("log seed")+
    xlab("Number of neighbors")+
    ylim(0, 10)+
    ggtitle("Fitness in relation with the abundance of neighbors 2020")+
    NULL
#join all the data ----
#2020
comp.20.final.junto$Plant_simple <- comp.20.final.junto$focal
agroup.h.20 <- h.20.group %>% group_by(Plot, Subplot, Plant_simple) %>% summarise (herb1 = sum(herb))%>%
    ungroup()
h.20.comp <- left_join(comp.20.final.junto, agroup.h.20)

s.herb <- h.20.comp[,c("Plot", "Subplot", "herb1", "Plant_simple")]
s.herb$Plant <- s.herb$Plant_simple
s.herb$herb1 <- as.numeric(s.herb$herb1)
s.herb$herb1[is.na(s.herb$herb1)]<- 0 #cambio los NAS por 0
agrup.fv.20 <- fv.20.group %>% group_by(Plot, Subplot, Plant) %>% summarise (fv1 = sum(fv))%>%
    ungroup()
s.herb.fv <- left_join(s.herb,agrup.fv.20 )
s.herb.fv$fv1[is.na(s.herb.fv$fv1)]<- 0 #cambio los NAS por 0 #datos de 2020 covariates
s.herb.fv.sal <- left_join(s.herb.fv,sal.prueba2 )
any(is.na(s.herb.fv.sal))
cov.clean.20 <-s.herb.fv.sal [,c("Plot", "Subplot", "herb1", "Plant","fv1", "sal")]

#2019 
agrup.h.19 <- h.19.group.1 %>% group_by(Plot, Subplot, Plant_Simple) %>% summarise (herb1 = sum(herb))%>%
    ungroup()
agrup.fv.19 <- fv.19 %>% group_by(Plot, Subplot, Plant) %>% summarise (fv1 = sum(fv))%>%
    ungroup()
comp.19.final.junto$Plot <- comp.19.final.junto$plot
comp.19.final.junto$Subplot <- comp.19.final.junto$subplot
comp.19.final.junto$Plant_Simple <- comp.19.final.junto$focal
h.19.comp <- left_join(comp.19.final.junto, agrup.h.19)
s.herb.19 <- h.19.comp[,c("Plot", "Subplot", "herb1", "Plant_Simple")]
s.herb.19$Plant <- s.herb.19$Plant_Simple
s.herb.19$herb1 <- as.numeric(s.herb.19$herb1)
s.herb.19$herb1[is.na(s.herb.19$herb1)]<- 0 #cambio los NAS por 0
s.herb.19$Plot <- as.numeric(s.herb.19$Plot)
agrup.fv.19$Plot <- as.numeric(agrup.fv.19$Plot)
s.herb.fv.19 <- left_join(s.herb.19,agrup.fv.19 )
s.herb.fv.19$fv1[is.na(s.herb.fv.19$fv1)]<- 0 #cambio los NAS por 0 #datos de 2020 covariates
s.herb.fv.19.sal <- left_join(s.herb.fv.19, sal.19.simple)#no hay NAs

cov.clean.19 <-s.herb.fv.19.sal [,c("Plot", "Subplot", "herb1", "Plant","fv1","sal")]


#join both dataset: 2019 and 2020
total <- rbind(cov.clean.20,cov.clean.19)


#final data of competition and covariates. Save them and reload them in the next script.


write.csv2(total, file ="C:/Users/maria/Documents/Tesis/R_repositorios/Multitrophic_plants/covariates_salt_h_fv_20_19.f.csv",
           row.names = FALSE)

save(total, file = "cov_total.Rda")
