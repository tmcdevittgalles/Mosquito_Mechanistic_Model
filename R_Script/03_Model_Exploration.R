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


##### Exploring the model 


parameter <- c(
  
  #### Demographic rates 
  
  2, ## Fecundity
  0.4, ## adult mortality
  0.6, ## larval mortality
  0.2, ## Development rate
  
  #### Thermal performance curve traits
  
  25, ## development optimal temp
  40, ## development maximum temp
  2 ## Sigma
)

inits <- c(
  1000, # initial juvenile population size
  0 # initial adult population size
)


### exploring the initial model 

## Selecting Harvard forest site and 2017 year

dum.df <- filter(temp.df , fYear == 2017 & Site == "HARV" )

predicted <- Mosq_TwoStage(time=dum.df$DOY, init=inits, pars = parameter,
                       temp =  dum.df$MaxTemp)



pop_dynamics <- predicted %>% filter(DOY >124 ) %>% 
  ggplot() + geom_line(aes(x=DOY, y= Adult,color='1'), 
                       size=2,alpha=.5)+
  geom_line(aes(x=DOY, y= Juv ,color='2'),
            size=2,alpha=.5)+ ylab("Abundance") +
  scale_color_manual(values = c( "1" = "Navy Blue","2"= "Brown"), 
                     labels=c("Adult", "Juv"), name= "Life stage")+
  theme_classic()

pop_dynamics


##
dev_dynamics <- predicted %>% 
  ggplot() + geom_line(aes(x=DOY, y= DevRate,color='1'), 
                       size=2,alpha=.5) +
  geom_line(aes(x=DOY, y= (Temp/max(Temp))*.2 ,color ="2"),
            size=2,alpha=.5)+ theme_classic()+
  scale_color_manual(values = c( "1" = "Navy Blue","2"= "Brown"), 
                     labels=c("Development rate", "Temperature"), name= "") +
  scale_y_continuous(
    "Development rate", 
    sec.axis = sec_axis(~ (. /.20) *max(predicted$Temp), name = "Temperature (C)")
  )

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
  0.3, ## adult mortality
  0.6, ## larval mortality
  0.2 ## Development rate
)

TPC_Par <- c(
    # Fecundity
  28, ## development optimal temp
  40, ## development maximum temp
  1.5, ## Sigma
    # Adult mortality
  20, ## development optimal temp
  40, ## development maximum temp
  2.5, ## Sigma
    # Larval morality
  26, ## development optimal temp
  30, ## development maximum temp
  2.5, ## Sigma
    # Larval development rate
  25, ## development optimal temp
  40, ## development maximum temp
  1.5 ## Sigma
)

inits <- c(
  10, # initial juvenile population size
  0 # initial adult population size
)


dum.df <- filter(temp.df , fYear == 2017 & Site == "HARV" )

predicted <- Mosq_AllTPC(time=dum.df$DOY, init=inits, ParBase =  Base_Par,
                   TPC = TPC_Par, temp =  dum.df$MaxTemp)



pop_dynamics <- predicted %>% filter(DOY >124 ) %>% 
  ggplot() + geom_line(aes(x=DOY, y= Adult,color='1'), 
                       size=2,alpha=.5)+
  geom_line(aes(x=DOY, y= Juv ,color='2'),
            size=2,alpha=.5)+ ylab("Abundance") +
  scale_color_manual(values = c( "1" = "Navy Blue","2"= "Brown"), 
                     labels=c("Adult", "Juv"), name= "Life stage")+
  theme_classic()

pop_dynamics


dev_dynamics <-predicted %>% filter(DOY >124 & DOY < 300) %>%
  ggplot() + geom_line(aes(x=DOY, y= 1/DevRate),
                       color='Navy Blue',size=2,alpha=.5) +
  ylab("Days to adult") + theme_classic()


fec_dynamics <-predicted %>% filter(DOY >124 & DOY < 300) %>%
  ggplot() + geom_line(aes(x=DOY, y= FecRate),
                       color='dark red',size=2,alpha=.5)+
  ylab("Number of eggs per day per mosquito") + theme_classic()

lmort_dynamics <- predicted %>% filter(DOY >124 & DOY < 300) %>%
  ggplot() + geom_line(aes(x=DOY, y=lMortRate),
                       color='Dark green',size=2,alpha=.5)+
  ylab("Larvae mortality rate")+ theme_classic()


amort_dynamics <- predicted %>% filter(DOY >124 & DOY < 300) %>%
  ggplot() + geom_line(aes(x=DOY, y= aMortRate),
                       color='dark orange',size=2,alpha=.5) +
  ylab("Adult mortality rate") + theme_classic()


(fec_dynamics | amort_dynamics)/(dev_dynamics|lmort_dynamics)
