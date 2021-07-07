###### Powell Center: Phenological patterns of mosquitoes #######
# Travis McDevitt-Galles
# 06/30/2021
# title: 01_Mosquito_Abundance_Functions.R

# Mechanistic model for mosquito abundance patterns. 

## The TPC function comes from Childress and Letcher 2017 - Ecology
## Estimated parameters (Toptim and CTmax) estimated based on reported values
## from Mordecai et al. 2019


#### Building a thermal performance curve for the different thermal sensitive
#### Traits

TPC <- function( Temp, Toptim, CTmax, sigma){
  
  perform <- rep(NA, length(Temp)) # Creating a vector of trait value based on
  # given temperature
  
  for( i in 1:length(Temp)){ # Looping through each temperature value
    if( Temp[i] <= Toptim){ # If temperature is less than or equal to Toptim
      perform[i] <- exp(1)^( (Temp[i] - Toptim)/ (2 *sigma)^2  ) 
    }else{ # if temperature is greater than the optim temp
      perform[i] <- 1 - ( (Temp[i] - Toptim) / (Toptim - CTmax) )^2
    }
  } 
  
  perform[perform < 0] <- 0 
  return(perform)
  
}
# 
# ## testing the TPC function
# temp <- seq(-5, 45, length.out = 100) # initial temperature range
# 
# at1 <- TPC( Temp = temp, Toptim = 22, CTmax= 35, sigma = 1)
# 
# # plotting the thermal performance curve
# 
# plot(x= temp, y= at1)
# doesn't look that good but we can carry on with it


################################################################################
################################################################################
###################                                #############################
###################    Mosquito 2 stage model      #############################
###################                                #############################
###################  Thermal sensitive Dev rate    #############################
###################                                #############################
################################################################################
################################################################################


Mosq_TwoStage <- function( time, init, pars, temp ){
  
  # warning messages
  
 if( length(time) != length(temp) ){
  
     stop("Error: Temperature and time steps are not equal") 
 } 
 if( length(init) != 2){

     stop("Error: not enought initial conditions, need 2")
 } 
 if( length(pars) != 7){
   
     stop("Error: not enought parameters, need 7") 
 }

   
   ## setting data frame to store predicted values 
   
   m.df <- data.frame( Juv = as.numeric( rep(NA, length= length(time))),
                       Adult = as.numeric( rep(NA, length = length(time))),
                       DevRate = as.numeric( rep(NA, length = length(time))))
   
   m.df$DOY <- time
   
   m.df$Temp <- temp
   ## Setting inital conditions
   
   m.df$Juv[1] <- init[1]
   m.df$Adult[1] <- init[2]
   
   ## Defining demographic parameters 
  
   fec <- pars[1]      ## number of new 'juvenile' mosquitoes added per day per 
                       ## female adult mosquito
   aMortal <- pars[2]  ## adult mortality rate: proportion of adults that die 
                       ## per day (24 hr)
   lMortal <- pars[3]  ## juvenile mortality rate: proportion of larvae that die 
                       ## per day (24 hr)
   lDevelop <- pars[4] ## max development rate: proportion of juveniles that 
                       ## develop into the next stage (Adult)
   
   ## Defining developmental thermal performance curve parameter
   
   Toptim <- pars[5]  ## temperature where development rate is maximized
   
   CTmax <-  pars[6]  ## Upper threshold temperature where above it parameter
                      ## drops to 0
   sigma <-  pars[7]  ## rate of increase to Toptim 
   
   
   ## Discrete step model ##
   
   for( i in 1:(nrow(m.df)-1)){
     ## Larval dynamics
     m.df$Juv[i+1] <-  m.df$Adult[i] * fec + ## new Juv from adults
                        m.df$Juv[i] * ( 1 - lMortal) - # Proportion surviving
       m.df$Juv[i] *  TPC( 20,  Toptim = Toptim , # proportion devel.
                                            CTmax =  CTmax , # into adults
                                            sigma = sigma) * lDevelop 

     ## Adult dynamics
     m.df$Adult[i+1] <-  m.df$Adult[i] * ( 1-aMortal) + # adult survival
       m.df$Juv[i] *  TPC( temp[i],  Toptim = Toptim , # proportion devel.
                            CTmax =  CTmax , # into adults
                            sigma = sigma) * lDevelop
     
     # tracking the development rate across time
     m.df$DevRate[i] <-  TPC( m.df$Temp[i],  Toptim = Toptim , # proportion devel.
                                             CTmax =  CTmax , # into adults
                                             sigma = sigma) * lDevelop 
   
   }

  
  return(m.df)
  
}

####### Exploring the model 


# parameter <- c(
#   
#   #### Demographic rates 
#   
#   3, ## Fecundity
#   0.55, ## adult mortality
#   0.4, ## larval mortality
#   0.15, ## Development rate
#   
#   #### Thermal performance curve traits
#   
#   20, ## development optimal temp
#   25, ## development maximum temp
#   1.5 ## Sigma
# )
# 
# inits <- c(
#   1000, # initial juvenile population size
#   0 # initial adult population size
# )
# 
# 
# ### exploring the initial model 
# 
# 
# dum.df <- select(pred.df,c("DOY", "Temp"))
# 
# at1 <- Mosq_TwoStage(time=dum.df$DOY, init=inits, pars = parameter,
#                      temp =  dum.df$Temp)
# 
# 
# 
# pop_dynamics <- at1 %>% filter(DOY >124 ) %>% 
#   ggplot() + geom_line(aes(x=DOY, y= Adult), 
#                        color='Navy Blue',size=2,alpha=.5) +
#   geom_line(aes(x=DOY, y= Juv),
#             color ="Brown",size=2,alpha=.5)
# 
# pop_dynamics
# 
# 
# 
# dev_dynamics <- at1 %>% filter(DOY >124 ) %>% 
#   ggplot() + geom_line(aes(x=DOY, y= DevRate), 
#                        color='Navy Blue',size=2,alpha=.5) +
#   geom_line(aes(x=DOY, y= (Temp/max(at1$Temp))*.2 ),
#             color ="Brown",size=2,alpha=.5)
# 
# dev_dynamics

################################################################################
################################################################################
###################                                #############################
###################    Mosquito 2 stage model      #############################
###################                                #############################
###################  Thermal sensitive all traits  #############################
###################                                #############################
################################################################################
################################################################################


Mosq_AllTPC <- function( time, init, ParBase, TPC ,temp ){
  
  # warning messages
  
  if( length(time) != length(temp) ){
    
    
    stop("Error: Temperature and time steps are not equal") 
  }
    if( length(init) != 2){
      
      stop("Error: not enought initial conditions, need 2")
    }
      if( length(ParBase) != 4 ){
        
        stop("Error: not enought base parameters, need 4") 
      }
      if( length(TPC) != 12){
        
        stop("Error: not enought TPC parameters, need 7") 
      }
    
    ## setting data frame to store predicted values 
    
    m.df <- data.frame( Juv = as.numeric( rep(NA, length= length(time))),
                        Adult = as.numeric( rep(NA, length = length(time))),
                        FecRate = as.numeric( rep(NA, length = length(time))),
                        aMortRate = as.numeric( rep(NA, length = length(time))),
                        lMortRate = as.numeric( rep(NA, length = length(time))),
                        DevRate = as.numeric( rep(NA, length = length(time))))
    
    m.df$DOY <- time
    
    m.df$Temp <- temp
    ## Setting inital conditions
    
    m.df$Juv[1] <- init[1]
    m.df$Adult[1] <- init[2]
    
    ## Defining demographic parameters 
    
    fec <- ParBase[1]      ##  max number of new 'juvenile' mosquitoes added 
                           ## per day per female adult mosquito
    aMortal <- ParBase[2]  ## adult mortality rate: proportion of adults that die 
                           ## per day (24 hr)
    lMortal <- ParBase[3]  ## juvenile mortality rate: proportion of larvae that die 
                           ## per day (24 hr)
    lDevelop <- ParBase[4] ## max development rate: proportion of juveniles that 
                           ## develop into the next stage (Adult)
    
    ## Defining developmental thermal performance curve parameter

    ########################## Fecundity #####################################
    
    FecToptim <- TPC[1]   ## temperature where fecundity is maximized
    
    FecCTmax <-  TPC[2]   ## Upper threshold temperature where above it parameter
                           ## drops to 0
    
    Fecsigma <-  TPC[3]   ## rate of increase to Toptim 
    
    ########################## Adult mortality #################################
    
    aMortToptim <- TPC[4]  ## temperature where development rate is maximized
    
    aMortCTmax <-  TPC[5]  ## Upper threshold temperature where above it parameter
                            ## drops to 0
    
    aMortsigma <-  TPC[6]  ## rate of increase to Toptim 
    
    ########################## Juv mortality #################################
    
    lMortToptim <- TPC[7]  ## temperature where development rate is maximized
    
    lMortCTmax <-  TPC[8]  ## Upper threshold temperature where above it parameter
                            ## drops to 0
    
    lMortsigma <-  TPC[9]  ## rate of increase to Toptim 
    
    ########################## development rate ################################
    
    DevToptim <- TPC[10]  ## temperature where development rate is maximized
    
    DevCTmax <-  TPC[11]  ## Upper threshold temperature where above it parameter
                           ## drops to 0
    
    Devsigma <-  TPC[12]  ## rate of increase to Toptim 
    
    ## Discrete step model ##
    
    for( i in 1:(nrow(m.df)-1)){
      
      m.df$DevRate[i] <- TPC( m.df$Temp[i],
                                 Toptim =  DevToptim ,
                                 CTmax= DevCTmax,
                                 sigma = Devsigma) * lDevelop
      
      
      m.df$FecRate[i] <- TPC( m.df$Temp[i],
                                 Toptim = FecToptim,
                                 CTmax= FecCTmax,
                                 sigma = Fecsigma ) * fec
      
      
      m.df$aMortRate[i] <- ( TPC( m.df$Temp[i],
                                       Toptim = aMortToptim,
                                       CTmax= aMortCTmax,
                                       sigma = aMortsigma ) ) * aMortal
      
      
      m.df$lMortRate[i] <- ( TPC( m.df$Temp[i],
                                       Toptim =  lMortToptim,
                                       CTmax=  lMortCTmax,
                                       sigma = lMortsigma ) ) * lMortal
      
      # tracking the development rate across time
      m.df$DevRate[i] <-  TPC( m.df$Temp[i],  Toptim = DevToptim , # proportion devel.
                               CTmax =  DevCTmax , # into adults
                               sigma = Devsigma) * lDevelop 
      
      
      ######### Discrete time step model ##########
      
      
      ## Larval dynamics
      m.df$Juv[i+1] <-  m.df$Adult[i] * m.df$FecRate[i] + # new Juv from adults
       m.df$Juv[i] - m.df$Juv[i] * (  m.df$lMortRate[i] ) - # Propotion surviving
        m.df$Juv[i] *  m.df$DevRate[i] # Proportion emerging
      
      ## Adult dynamics
      m.df$Adult[i+1] <-m.df$Adult[i]  - m.df$Adult[i] * ( m.df$aMortRate[i]) + # adult survival
        m.df$Juv[i] *  m.df$DevRate[i] # Proportion emerging
      
  }
  
  return(m.df)
  
}


# 
# Base_Par <- c(
#  2, ## Fecundity
#   0.5, ## adult mortality
#   0.4, ## larval mortality
#   0.15 ## Development rate
# )
# 
# TPC_Par <- c(
#     # Fecundity 
#   20, ## development optimal temp
#   25, ## development maximum temp
#   1.5, ## Sigma
#     # Adult mortality
#   15, ## development optimal temp
#   25, ## development maximum temp
#   5.5, ## Sigma
#     # Larval morality
#   15, ## development optimal temp
#   25, ## development maximum temp
#   5.5, ## Sigma
#     # Larval development rate
#   20, ## development optimal temp
#   25, ## development maximum temp
#   1.5 ## Sigma
# )
# 
# inits <- c(
#   10, # initial juvenile population size
#   0 # initial adult population size
# )
# 
# 
# dum.df <- select(pred.df,c("DOY", "Temp"))
# 
# at1 <- Mosq_AllTPC(time=dum.df$DOY, init=inits, ParBase =  Base_Par,
#                    TPC = TPC_Par, temp =  dum.df$Temp)
# 
# 
# pop_dynamics <- at1 %>% filter(DOY >124 ) %>% 
#   ggplot() + geom_line(aes(x=DOY, y= Adult), 
#                        color='Navy Blue',size=2,alpha=.5) +
#   geom_line(aes(x=DOY, y= Juv),
#             color ="Brown",size=2,alpha=.5)
# 
# pop_dynamics
# 
# 
# dev_dynamics <- at1 %>% filter(DOY >124 & DOY < 300) %>% 
#   ggplot() + geom_line(aes(x=DOY, y= 1/DevRate), 
#                        color='Navy Blue',size=2,alpha=.5) 
# 
# 
# fec_dynamics <- at1 %>% filter(DOY >124 & DOY < 300) %>% 
#   ggplot() + geom_line(aes(x=DOY, y= FecRate), 
#                        color='dark red',size=2,alpha=.5) 
# 
# 
# lmort_dynamics <- at1 %>% filter(DOY >124 & DOY < 300) %>% 
#   ggplot() + geom_line(aes(x=DOY, y=lMortRate), 
#                        color='Dark green',size=2,alpha=.5) 
# 
# 
# amort_dynamics <- at1 %>% filter(DOY >124 & DOY < 300) %>% 
#   ggplot() + geom_line(aes(x=DOY, y= aMortRate), 
#                        color='dark orange',size=2,alpha=.5) 
# 
# 
# library(patchwork)
# 
# (fec_dynamics | amort_dynamics)/(dev_dynamics|lmort_dynamics)
