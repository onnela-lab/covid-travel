library(tidyverse)
library(dplyr)
library(readxl)
setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/Realdata/Datapreprocessing")
load("flightdataJanFeb2020.Rdata")
load("flightdataMarch2020.Rdata")
load("flightdataApril2020.Rdata")
load("flightdataMay2020.Rdata")
load("flightdataJune2020.Rdata")

file_name_active <- paste("time_series_covid19_confirmed_global.csv",sep="")
active_data <- read.csv(file_name_active)

file_name_recovered <- paste("time_series_covid19_recovered_global.csv",sep="")
recovered_data <- read.csv(file_name_recovered)

file_name_death <- paste("time_series_covid19_deaths_global.csv",sep="")
death_data <- read.csv(file_name_death)

file_name_ISO <- paste("PairFlightcountriesandCovidAndAlpha3.xlsx",sep="")
alpha3 = read_excel(file_name_ISO)
alpha3flight = alpha3[,-2]
alpha3flight = alpha3flight%>%dplyr::rename(country = 'Countries Flightdata')
alpha3covid = alpha3[,-1]
alpha3covid = alpha3covid%>%dplyr::rename(country = 'CountriesCovidData')%>%
              dplyr::distinct()


file_name_population  <- paste("countriespopulation.csv",sep="")
population_data <- read.csv(file_name_population)

population_data <- population_data%>%dplyr::rename(ALPHA3 = alpha3)%>%
  dplyr::select(ALPHA3,population)

#Define a function to convert the date format
date_converter <- function(data){
  # This function takes the active, recovered, death data and converts the format of the date
  # names in the columns to be consistent with the travel data ie. 2020-01-23. Original format
  # must be of the form "X1.23.20" 
  # data
  # Either the recoverd, active, or death data whose column names you want to change
  
  act_dates <- colnames(data)[3:dim(data)[2]] # the dates start at the 3rd column for this data
  act_dates <- substr(act_dates, 2, nchar(act_dates))
  act_dates_split <- strsplit(act_dates, ".", fixed = TRUE)
  act_dates_changed <- NULL
  for(i in 1:length(act_dates_split)){
    dte <- act_dates_split[[i]]
    n1 <- which(nchar(dte) == 1)
    if(length(n1) > 0){ # turning 1.23.20 --> 2020-01-23
      dte[n1] <- paste(0, dte[n1], sep = "")
    }
    act_dates_changed[i] <- paste("", dte[3], "-", dte[1], "-", dte[2], sep = "")
  }
  colnames(data)[3:dim(data)[2]] <- act_dates_changed
  return(data)
}
# Changing dates in the colnames from "X2.29.20" to a "2020-02-29" format
active_data <- date_converter(active_data)
recovered_data <- date_converter(recovered_data)
death_data <- date_converter(death_data)

# Keep only data up to the interested point
dates = seq(as.Date("2020/01/22"), as.Date("2020/06/30"), "days") #04/30/2020 as the cut off point
totaldays = length(dates)
cutoff = totaldays +2 #Plus 2 for including the first two column

active_data = active_data[,2:cutoff] # remove province and after cutoff
recovered_data = recovered_data[,2:cutoff]
death_data = death_data[,2:cutoff]

########Sum of cases for all countries with data reported at the province level
active_data <- active_data %>% dplyr::rename( country = Country.Region)%>%
              dplyr::group_by(country)%>% summarise_all(funs(sum))

recovered_data <- recovered_data %>% dplyr::rename( country = Country.Region)%>%
  dplyr::group_by(country)%>% summarise_all(funs(sum))

death_data <- death_data %>% dplyr::rename( country = Country.Region)%>%
  dplyr::group_by(country)%>% summarise_all(funs(sum))


#####Make Alpha 3 code for all countries
active_data = active_data%>%dplyr::left_join(alpha3covid)%>%dplyr::relocate(ALPHA3)
recovered_data = recovered_data%>%dplyr::left_join(alpha3covid)%>%dplyr::relocate(ALPHA3)
death_data = death_data%>%dplyr::left_join(alpha3covid)%>%dplyr::relocate(ALPHA3)


##Change OTHERS to ZZZ code
active_data <- active_data%>%dplyr::mutate(ALPHA3 = gsub("OTHERS","ZZZ", ALPHA3))%>%
  dplyr::select(-country)%>%dplyr::arrange(ALPHA3)

recovered_data <- recovered_data%>%dplyr::mutate(ALPHA3 = gsub("OTHERS","ZZZ", ALPHA3))%>%
  dplyr::select(-country)%>%dplyr::arrange(ALPHA3)



death_data <- death_data%>%dplyr::mutate(ALPHA3 = gsub("OTHERS","ZZZ", ALPHA3))%>%
             dplyr::select(-country)%>%dplyr::arrange(ALPHA3)



##Extract countries name with more than 100 confirmed cases at March 15 or column 56 
tt = dplyr::filter(active_data, active_data[,56]>100)
tt = tt[,"ALPHA3"]
countriesvector = rep(0,nrow(tt))
for(i in 1:nrow(tt)){
  countriesvector[i] = as.character(tt[i,])
}
countriesvector = sort(unique(countriesvector))
countriesvector = countriesvector[-length(countriesvector)]

#####################Extract data for countries with more than 100 confirmed cases at March 15 or column 56 

active_names = active_data[,"ALPHA3"]

index = c()
for(i in 1:nrow(active_names)){
  tmp = as.character(active_names[i,])
  if( tmp %in% countriesvector){
    index=c(index,i)
  }
  
}

active_data1 = active_data[index,]
recovered_data1 = recovered_data[index,]
death_data1 = death_data[index,] 

active_data1 = active_data1%>%dplyr::left_join(population_data)%>%dplyr::relocate(ALPHA3,population)

recovered_data1 = recovered_data1%>%dplyr::left_join(population_data)%>%dplyr::relocate(ALPHA3,population)

death_data1 = death_data1%>%dplyr::left_join(population_data)%>%dplyr::relocate(ALPHA3,population)

##########define a function collecting A,R,D
ncolumn = dim(active_data1)[2]
duration = ncolumn - 2
datalist = list()
numbercountrieshighcases =  length(countriesvector)

datafunction = function(country){
  matrix = matrix(0,nrow = duration,3)
  activeconfirmed = active_data1[country,3:ncolumn] - recovered_data1[country,3:ncolumn] - death_data1[country,3:ncolumn]
  tmp1 = unlist(activeconfirmed)
  tmp1 = as.vector(tmp1)
  matrix[,1] = tmp1
  
  recovered = recovered_data1[country,3:ncolumn]
  tmp2 = unlist(recovered)
  tmp2 = as.vector(tmp2)
  matrix[,2] = tmp2
  
  death = death_data1[country,3:ncolumn]
  tmp3 = unlist(death)
  tmp3 = as.vector(tmp3)
  matrix[,3] = tmp3
  
  return(matrix)
}

for (i in 1:numbercountrieshighcases){
  datalist[[countriesvector[i]]] =  datafunction(i)
}

#####Adding ZZZ as all 0 for consistent
datalist[["ZZZ"]] = matrix(0,nrow = duration,3)

###Adding name to keep track
datalist[["country"]] = c(countriesvector,"ZZZ")

###Adding corresponding population, ZZZ corresponding to 0
population = active_data1[,"population"]
population =  unlist(population)
population = as.vector(population)
datalist[["population"]] = c(population,0)
coviddataJanJune = datalist
#Save data for later use
save(coviddataJanJune, file ="coviddataJanJune.Rdata")
############################MAKE FLIGHT MATRIX BECOME COUNTRIES BELONG TO THE LIST TO OTHERS####

#CHANGE ALL NAMES NOT BELONG TO countries vector

flightdata1 = c(flightdataJanFeb2020, flightdataMarch2020, flightdataApril2020, flightdataMay2020, flightdataJune2020)

flightdata_rearrangeJanJune20  = list()

for(d in 1:length(flightdata1)){



dat = flightdata1[[d]]
flight_names = colnames(dat)


for(i in 1:length(flight_names)){
  tmp = flight_names[i]
  if( tmp %in% countriesvector){
    flight_names[i] = tmp
  } else{flight_names[i] = "ZZZ"}
  
}




colnames(dat) = flight_names
rownames(dat) = flight_names

###Travelers between high cases countries

dat1 = dat[countriesvector,countriesvector]

###Travelers from each countries to ZZZ #########


index=c()
for(i in 1:length(flight_names)){
  
  tmp = flight_names[i]
  if(tmp %in% countriesvector){
    index = c(index,i)
  }
}

dat12 = dat[index,]# take only departure from countries with high covid
dat12 = dat12[,-index] # remove the arrival of high cases
dat12 = rowSums(dat12)

########Travelers from ZZZ to each countries #########
dat21 = dat[-index,]# take only departure from countries with high covid
dat21 = dat21[,index] # remove the arrival of high cases
dat21 = colSums(dat21)

########Create the new matrix
n1 = length(countriesvector)  #
n2 = n1 + 1
flightrearrange = matrix(0,n2,n2)
flightrearrange[1:n1,1:n1] = dat1
flightrearrange[1:n1,n2] = dat12
flightrearrange[n2,1:n1] = dat21
rearrangenames = c(countriesvector, "ZZZ")
rownames(flightrearrange) = rearrangenames
colnames(flightrearrange) = rearrangenames

flightdata_rearrangeJanJune20[[d]] = round(flightrearrange, digits=0)

}

save(flightdata_rearrangeJanJune20, file="flightdata_rearrangeJanJune20.Rdata")
