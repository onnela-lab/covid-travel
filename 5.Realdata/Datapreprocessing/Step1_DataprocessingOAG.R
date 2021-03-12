library(tidyverse)
library(dplyr)
library(readxl)
#This code is created to process the OAG data into matrix form
#setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/Realdata")

file_name = paste("OAGdata.csv",sep="")
flightdat <- read.csv(file_name)
OAG = flightdat

file_name1 = paste("PairFlightcountriesandCovidAndAlpha3.xlsx",sep="")
alpha3 = read_excel(file_name1)
alpha3flight = alpha3[,-2]
alpha3flight = alpha3flight%>%dplyr::rename(country = 'Countries Flightdata')

#Rename and remove internal flight at step 1
OAG <- OAG %>% dplyr::rename( depart = Dep.IATA.Country.Name, arrive = Arr.IATA.Country.Name, 
                           freq = Frequency, seats = Seats..Total. , day = Time.series)%>%
               dplyr::rowwise() %>%
               dplyr::filter(depart != arrive)

#Standardizing ALPHA 3 by merging territories
#standardize depart

OAG_depart <- OAG %>% 
  dplyr::select(day, seats, country = depart) %>%
  dplyr::left_join(alpha3flight)
#standardize arrive
OAG_arrive <- OAG %>% 
  dplyr::select(day, country = arrive) %>%
  dplyr::left_join(alpha3flight)
####Pull back depart and arrive
OAG_depart$arrivecountry = OAG_arrive$country
OAG_depart$arrivecountryALPHA3 = OAG_arrive$ALPHA3
#############select only Alpha 3 depart, arrive, day, total number of seats
OAGdata <- OAG_depart%>%dplyr::select(day,seats,ALPHA3, arrivecountryALPHA3)%>%
  dplyr::rename(depart= ALPHA3, arrive = arrivecountryALPHA3)

##Replace others by ZZZ code, our own notation to make OTHERS at the end
OAGdata <- OAGdata%>%dplyr::mutate(depart = gsub("OTHERS","ZZZ", depart), arrive = gsub("OTHERS","ZZZ", arrive) )
##Make sure all arrive and depart are different after standardized
OAGdata <- OAGdata%>%
  dplyr::rowwise() %>%
  dplyr::filter(depart != arrive)

##Summation number of passengers, group by day, depart, arrive

dates <- dplyr::arrange(unique(OAGdata[, "day"]))

flightdataJanFeb2020 = list()

for(d in 1:nrow(dates)){

day1 = dates[d,]



OAGdata1 = OAGdata%>%dplyr::filter(day==day1)
#########

OAGdata1 = OAGdata1%>%
  dplyr::group_by(day, depart, arrive) %>%
  dplyr::summarise(seats = sum(seats))%>%
  dplyr::arrange(depart) #make sure things in Alphabet order at least for departure
##Fill values of OAGdata1 in a matrix with alphabet order of ALPHA3
countriesdepart = unique(as.data.frame(OAGdata1[,"depart"]))
countriesarrive = unique(as.data.frame(OAGdata1[,"arrive"]))

countriesdepart = countriesdepart[,"depart"] # extract names only
countriesarrive = countriesarrive[,"arrive"]
countries = sort(unique(union(countriesdepart,countriesarrive )))

numbercountries = length(countries)

flightmatrix = matrix(0, numbercountries, numbercountries)


for(i in 1:numbercountries){
  for(j in 1:numbercountries){
    departcountry = countries[i]
    arrivecountry = countries[j]
    OAGdata2 <- OAGdata1%>%dplyr::filter(depart==departcountry, arrive==arrivecountry)
    OAGdata2 = as.data.frame(OAGdata2)
    tmp = OAGdata2[,"seats"]
    if(length(tmp>0)){flightmatrix[i,j]= tmp}
     }
}


colnames(flightmatrix) = countries
rownames(flightmatrix) = countries

flightdataJanFeb2020[[d]] = flightmatrix
print(d/60) }

save(flightdataJanFeb2020, file="flightdataJanFeb2020.RData")

















