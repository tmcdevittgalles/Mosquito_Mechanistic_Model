###### Powell Center: Phenological patterns of mosquitoes #######
# Travis McDevitt-Galles
# 06/30/2021
# title: 02_WM_Temperature_Profile

# Extracting predicted temperature patterns for the different domains in NEON

##  Using GAM approachs to model seasonal max daily temperature for all neon Domains
## across the SIX years of our study. The predicted temperautre values will be 
## used to simulate seasonal mosquito population dynamics 

## Loading required libraries


library( dplyr )
library( tidyr )
library( ggplot2 )
library( gamm4 )

## setting wd to upload data
setwd("C:/Users/tmcdevitt-galles/powell-mosquito-phenology")

# input combined prism data downloaded and merged in "03_MP_PRISM_Data.R" 

# Daily Prism data
load("./Data/DailyPrismMod.Rda")

# Count data to get the domains and lat and long information
load("./Data/Mosquito_Data_Clean.Rda")

########### modeling climate patterns across domain ###########

site.df <- unique(select(ungroup(complete.df), c("Site", "Plot", "Domain")))

# Need to add Site info to prism data 
contig.df <- left_join(contigus.df, site.df , by="Plot")


### Modeling four differnt sites , DSNY == Southeast,
###                                HARV == Northeast
###                                TEAK == Southwest
###                                UNDE == Great lakes


######################### DSNY ###################################

focal.df <- filter( contig.df, Site == "DSNY" )  # Southeast
focal.df$fYear <- as.factor(focal.df$Year)

### Temperature
gam.temp <- gam( Tmax7 ~ s(DOY,by=interaction(fYear), bs="cc"),
                   data=focal.df, family="gaussian")


summary(gam.temp)
## lets try and plot this with new data

pred.df <- select( focal.df, c("DOY", "fYear", "Site"))

dum.df <-  unique(pred.df)

dum.df$Pred <- predict(gam.temp, newdata = dum.df)

# Looking at predicted number
dum.df %>% filter( fYear != "2013") %>% 
  ggplot(  aes(x=DOY, y=(Pred), color=fYear))+ 
  geom_point(data=focal.df, aes(x= DOY,y=Tmax7, color=fYear),size=1)+
  #geom_point( data=at1, aes(x= DOY, y=PPT14))+
  geom_line(size=2,alpha=.75)+ theme_classic()+ 
  facet_wrap(~Site)

# putting predicted temp data to new data set
temp.df <- dum.df



######################### HARV ###################################

focal.df <- filter( contig.df, Site == "HARV")  # Southeast
focal.df$fYear <- as.factor(focal.df$Year)

### Temperature
gam.temp <- gam( Tmax7 ~ s(DOY,by=interaction(fYear), bs="cc"),
                 data=focal.df, family="gaussian")


summary(gam.temp)
## lets try and plot this with new data

pred.df <- select( focal.df, c("DOY", "fYear", "Site"))

dum.df <-  unique(pred.df)

dum.df$Pred <- predict(gam.temp, newdata = dum.df)

# Looking at predicted number
dum.df %>% filter( fYear != "2013") %>% 
  ggplot(  aes(x=DOY, y=(Pred), color=fYear))+ 
  geom_point(data=focal.df, aes(x= DOY,y=Tmax7, color=fYear),size=1)+
  geom_line(size=2,alpha=.75)+ theme_classic()+ 
  facet_wrap(~Site)

# putting predicted temp data to new data set
temp.df <- rbind.data.frame(temp.df, dum.df)

######################### TEAK ###################################

focal.df <- filter( contig.df, Site == "TEAK")  # Southwest
focal.df$fYear <- as.factor(focal.df$Year)

### Temperature
gam.temp <- gam( Tmax7 ~ s(DOY,by=interaction(fYear), bs="cc"),
                 data=focal.df, family="gaussian")


summary(gam.temp)
## lets try and plot this with new data

pred.df <- select( focal.df, c("DOY", "fYear", "Site"))

dum.df <-  unique(pred.df)

dum.df$Pred <- predict(gam.temp, newdata = dum.df)

# Looking at predicted number
dum.df %>% filter( fYear != "2013") %>% 
  ggplot(  aes(x=DOY, y=(Pred), color=fYear))+ 
  geom_point(data=focal.df, aes(x= DOY,y=Tmax7, color=fYear),size=1)+
  geom_line(size=2,alpha=.75)+ theme_classic()+ 
  facet_wrap(~Site)

# putting predicted temp data to new data set
temp.df <- rbind.data.frame(temp.df, dum.df)


######################### UNDE ###################################

focal.df <- filter( contig.df, Site == "UNDE")  # Southwest
focal.df$fYear <- as.factor(focal.df$Year)

### Temperature
gam.temp <- gam( Tmax7 ~ s(DOY,by=interaction(fYear), bs="cc"),
                 data=focal.df, family="gaussian")


summary(gam.temp)
## lets try and plot this with new data

pred.df <- select( focal.df, c("DOY", "fYear", "Site"))

dum.df <-  unique(pred.df)

dum.df$Pred <- predict(gam.temp, newdata = dum.df)

# Looking at predicted number
dum.df %>% filter( fYear != "2013") %>% 
  ggplot(  aes(x=DOY, y=(Pred), color=fYear))+ 
  geom_point(data=focal.df, aes(x= DOY,y=Tmax7, color=fYear),size=1)+

  geom_line(size=2,alpha=.75)+ theme_classic()+ 
  facet_wrap(~Site)

# putting predicted temp data to new data set
temp.df <- rbind.data.frame(temp.df, dum.df)



## Plotting all data

temp.df %>% filter( fYear != "2013") %>% 
  ggplot(  aes(x=DOY, y=(Pred), color=Site))+ 
  geom_line(size=2,alpha=.75)+ theme_classic()+ 
  facet_wrap(~fYear)


colnames(temp.df)[4] <- "MaxTemp"

## Saving temperature data as an rData file

saveRDS(temp.df, "temp.rds")
