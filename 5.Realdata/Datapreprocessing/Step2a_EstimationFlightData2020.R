library(tidyverse)
library(dplyr)
library(readxl)

setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/Realdata")

load("flightdataJanFeb2020.Rdata")

data = flightdataJanFeb2020
##################Estimate March 2020 based on Jan 2020 and flights reduction
flightdataMarch2020 = list()
for(d in 1:31){
flightmatrix = data[[d]]
countries = colnames(flightmatrix)
numbercountries = length(countries)
flightmatrix = as.matrix(flightmatrix)

###########Estimate number of travelers in 2020-03

file_name2 = paste("flight_reduction_scaling_factors_march20vsjanuary20.csv",sep="")
factor030120 <- read.csv(file_name2)
#rename for use 
factor030120 <- factor030120%>%dplyr::rename( depart = origin_country_iso_code, 
                                              arrive = destination_country_iso_code, scale = scaling_factor)%>%
  dplyr::select(-X)

commonfactor = factor030120%>%dplyr::ungroup()%>%dplyr::summarise(ave = mean(scale))
commonfactor = as.numeric(commonfactor)


# #######
flightmatrix1 = matrix(0, numbercountries, numbercountries)

for(i in 1:numbercountries){
  for(j in 1:numbercountries){
    departcountry = countries[i]
    arrivecountry = countries[j]
    factortmp <- factor030120%>%dplyr::filter(depart==departcountry, arrive==arrivecountry)
    factortmp = as.data.frame(factortmp)
    tmp = as.numeric(factortmp[,"scale"])
    #print(length(tmp))
    if(length(tmp>0)){flightmatrix1[i,j]= tmp*flightmatrix[i,j]
    } else{
      flightmatrix1[i,j]= commonfactor*flightmatrix[i,j]
    }
  }
}



colnames(flightmatrix1) = countries
rownames(flightmatrix1) = countries

flightdataMarch2020[[d]] = flightmatrix1}


save(flightdataMarch2020, file = "flightdataMarch2020.Rdata")

#################################################

##################Estimate April 2020 based on Jan 2020 and flights reduction
flightdataApril2020 = list()
for(d1 in 1:30){
  flightmatrix = data[[d1]]
  countries = colnames(flightmatrix)
  numbercountries = length(countries)
  flightmatrix = as.matrix(flightmatrix)
  
  ###########Estimate number of travelers in 2020-04
  
  file_name2 = paste("flight_reduction_scaling_factors_april20vsjanuary20.csv",sep="")
  factor030120 <- read.csv(file_name2)
  #rename for use 
  factor030120 <- factor030120%>%dplyr::rename( depart = origin_country_iso_code, 
                                                arrive = destination_country_iso_code, scale = scaling_factor)%>%
    dplyr::select(-X)
  
  commonfactor = factor030120%>%dplyr::ungroup()%>%dplyr::summarise(ave = mean(scale))
  commonfactor = as.numeric(commonfactor)
  
  
  # #######
  flightmatrix1 = matrix(0, numbercountries, numbercountries)
  
  for(i in 1:numbercountries){
    for(j in 1:numbercountries){
      departcountry = countries[i]
      arrivecountry = countries[j]
      factortmp <- factor030120%>%dplyr::filter(depart==departcountry, arrive==arrivecountry)
      factortmp = as.data.frame(factortmp)
      tmp = as.numeric(factortmp[,"scale"])
      #print(length(tmp))
      if(length(tmp>0)){flightmatrix1[i,j]= tmp*flightmatrix[i,j]
      } else{
        flightmatrix1[i,j]= commonfactor*flightmatrix[i,j]
      }
    }
  }
  
  
  
  colnames(flightmatrix1) = countries
  rownames(flightmatrix1) = countries
  
  flightdataApril2020[[d1]] = flightmatrix1}


save(flightdataApril2020, file = "flightdataApril2020.Rdata")

#################################################








