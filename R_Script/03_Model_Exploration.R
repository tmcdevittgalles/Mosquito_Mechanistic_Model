###### Powell Center: Phenological patterns of mosquitoes #######
# Travis McDevitt-Galles
# 06/30/2021
# title: 03_Model_Exploration

# Mechanistic model for mosquito abundance patterns. 

library(here)
library(dplyr)
library(ggplot2)
library(patchwork)

# loading model functions
source( here("R_Script", "01_Mosquito_Abundance_Functions.R" ))

# loading temperature data 
temp.df <- readRDS( here("Data", "temp.rdata") )



################################################################################
################################################################################
###################                                #############################
###################    Mosquito 2 stage model      #############################
###################                                #############################
###################  Thermal sensitive Dev rate    #############################
###################                                #############################
################################################################################
################################################################################


####### Exploring the model 


parameter <- c(
  
  #### Demographic rates 
  
  3, ## Fecundity
  0.55, ## adult mortality
  0.4, ## larval mortality
  0.15, ## Development rate
  
  #### Thermal performance curve traits
  
  20, ## development optimal temp
  25, ## development maximum temp
  1.5 ## Sigma
)

inits <- c(
  1000, # initial juvenile population size
  0 # initial adult population size
)


### exploring the initial model 

## Selecting Harvard forest site and 2017 year

dum.df <- filter(temp.df , fYear == 2017 & Site == "HARV" )

at1 <- Mosq_TwoStage(time=dum.df$DOY, init=inits, pars = parameter,
                       temp =  dum.df$MaxTemp)



pop_dynamics <- at1 %>% filter(DOY >124 ) %>% 
  ggplot() + geom_line(aes(x=DOY, y= Adult), 
                       color='Navy Blue',size=2,alpha=.5) +
  geom_line(aes(x=DOY, y= Juv),
            color ="Brown",size=2,alpha=.5)

pop_dynamics



dev_dynamics <- at1 %>% filter(DOY >124 ) %>% 
  ggplot() + geom_line(aes(x=DOY, y= DevRate), 
                       color='Navy Blue',size=2,alpha=.5) +
  geom_line(aes(x=DOY, y= (Temp/max(at1$Temp))*.2 ),
            color ="Brown",size=2,alpha=.5)

dev_dynamics



################################################################################
################################################################################
###################                                #############################
###################    Mosquito 2 stage model      #############################
###################                                #############################
###################  Thermal sensitive all traits  #############################
###################                                #############################
################################################################################
################################################################################


Base_Par <- c(
 2, ## Fecundity
  0.5, ## adult mortality
  0.4, ## larval mortality
  0.15 ## Development rate
)

TPC_Par <- c(
    # Fecundity
  20, ## development optimal temp
  25, ## development maximum temp
  1.5, ## Sigma
    # Adult mortality
  15, ## development optimal temp
  25, ## development maximum temp
  5.5, ## Sigma
    # Larval morality
  15, ## development optimal temp
  25, ## development maximum temp
  5.5, ## Sigma
    # Larval development rate
  20, ## development optimal temp
  25, ## development maximum temp
  1.5 ## Sigma
)

inits <- c(
  10, # initial juvenile population size
  0 # initial adult population size
)


dum.df <- filter(temp.df , fYear == 2017 & Site == "HARV" )

at1 <- Mosq_AllTPC(time=dum.df$DOY, init=inits, ParBase =  Base_Par,
                   TPC = TPC_Par, temp =  dum.df$MaxTemp)


pop_dynamics <- at1 %>% filter(DOY >124 ) %>%
  ggplot() + geom_line(aes(x=DOY, y= Adult),
                       color='Navy Blue',size=2,alpha=.5) +
  geom_line(aes(x=DOY, y= Juv),
            color ="Brown",size=2,alpha=.5)

pop_dynamics


dev_dynamics <- at1 %>% filter(DOY >124 & DOY < 300) %>%
  ggplot() + geom_line(aes(x=DOY, y= 1/DevRate),
                       color='Navy Blue',size=2,alpha=.5) +
  ylab("Days to adult") + theme_classic()


fec_dynamics <- at1 %>% filter(DOY >124 & DOY < 300) %>%
  ggplot() + geom_line(aes(x=DOY, y= FecRate),
                       color='dark red',size=2,alpha=.5)+
  ylab("Number of eggs per day per mosquito") + theme_classic()

lmort_dynamics <- at1 %>% filter(DOY >124 & DOY < 300) %>%
  ggplot() + geom_line(aes(x=DOY, y=lMortRate),
                       color='Dark green',size=2,alpha=.5)+
  ylab("Larvae mortality rate")+ theme_classic()


amort_dynamics <- at1 %>% filter(DOY >124 & DOY < 300) %>%
  ggplot() + geom_line(aes(x=DOY, y= aMortRate),
                       color='dark orange',size=2,alpha=.5) +
  ylab("Adult mortality rate") + theme_classic()


(fec_dynamics | amort_dynamics)/(dev_dynamics|lmort_dynamics)
