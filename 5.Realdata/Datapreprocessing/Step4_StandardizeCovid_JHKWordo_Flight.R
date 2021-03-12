##This section is created to bring JHK data, Worldodata, and JHK fixed data together

library(tidyverse)
library(dplyr)
library(readxl)

setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/Realdata/Datapreprocessing")
load("flightdataJanFeb2020.Rdata")
load("flightdataMarch2020.Rdata")
load("flightdataApril2020.Rdata")
load("flightdataMay2020.Rdata")
load("flightdataJune2020.Rdata")
load("flightdataMarch2019.Rdata")
load("flightdataApril2019.Rdata")
load("flightdataMay2019.Rdata")
load("flightdataJune2019.Rdata")

######################

cutoffday = 88#decide the day where all countries extract has more than 100 cases, 56= March 15, 72 = March31,88 = April 15
# Keep only data up to the interested point
dates = seq(as.Date("2020/01/22"), as.Date("2020/06/30"), "days") # cutoff 04/30 can be change later


############
file_name_active <- paste("time_series_covid19_confirmed_global.csv",sep="")
active_data <- read.csv(file_name_active)

file_name_recovered <- paste("time_series_covid19_recovered_global.csv",sep="")
recovered_data <- read.csv(file_name_recovered)

file_name_death <- paste("time_series_covid19_deaths_global.csv",sep="")
death_data <- read.csv(file_name_death)

file_name_active_worldo <- paste("cumcases_worldometers.csv",sep="")
active_data_worldo <- read.csv(file_name_active_worldo)


file_name_death_worldo <- paste("cumdeath_worldometers.csv",sep="")
death_data_worldo <- read.csv(file_name_death_worldo)
countriesworldo = read.csv("jhkvsworldocountrynames.csv", header =T)


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
active_data_worldo <- date_converter(active_data_worldo)
death_data_worldo <- date_converter(death_data_worldo)

active_data_worldo[is.na(active_data_worldo)] <- 0
death_data_worldo[is.na(death_data_worldo)] <- 0


########################

totaldays = length(dates)
cutoff = totaldays +2 #Plus 2 for including the first two column

active_data = active_data[,2:cutoff] # remove province and after cutoff
recovered_data = recovered_data[,2:cutoff]
death_data = death_data[,2:cutoff]
active_data_worldo = active_data_worldo[,2:cutoff] # remove province and after cutoff
death_data_worldo = death_data_worldo[,2:cutoff]


########Sum of cases for all countries with data reported at the province level
active_data <- active_data %>% dplyr::rename( country = Country.Region)%>%
              dplyr::group_by(country)%>% summarise_all(funs(sum))

recovered_data <- recovered_data %>% dplyr::rename( country = Country.Region)%>%
  dplyr::group_by(country)%>% summarise_all(funs(sum))

death_data <- death_data %>% dplyr::rename( country = Country.Region)%>%
  dplyr::group_by(country)%>% summarise_all(funs(sum))
##########################




active_data_worldo <- active_data_worldo %>% dplyr::rename( country1 = Country.Region)%>%
  dplyr::left_join(countriesworldo)%>%dplyr::relocate(country)%>%
  dplyr::select(- one_of(c("country1","X")))%>%
  dplyr::group_by(country)%>% summarise_all(funs(sum))


death_data_worldo <- death_data_worldo %>% dplyr::rename( country1 = Country.Region)%>%
  dplyr::left_join(countriesworldo)%>%dplyr::relocate(country)%>%
  dplyr::select(- one_of(c("country1","X")))%>%
  dplyr::group_by(country)%>% summarise_all(funs(sum))
#####Make Alpha 3 code for all countries

#########################
active_data = active_data%>%dplyr::left_join(alpha3covid)%>%dplyr::relocate(ALPHA3)
recovered_data = recovered_data%>%dplyr::left_join(alpha3covid)%>%dplyr::relocate(ALPHA3)
death_data = death_data%>%dplyr::left_join(alpha3covid)%>%dplyr::relocate(ALPHA3)



active_data_worldo = active_data_worldo%>%dplyr::left_join(alpha3covid)%>%dplyr::relocate(ALPHA3)
death_data_worldo = death_data_worldo%>%dplyr::left_join(alpha3covid)%>%dplyr::relocate(ALPHA3)

##Change OTHERS to ZZZ code
active_data <- active_data%>%dplyr::mutate(ALPHA3 = gsub("OTHERS","ZZZ", ALPHA3))%>%
  dplyr::select(-country)%>%dplyr::arrange(ALPHA3)

recovered_data <- recovered_data%>%dplyr::mutate(ALPHA3 = gsub("OTHERS","ZZZ", ALPHA3))%>%
  dplyr::select(-country)%>%dplyr::arrange(ALPHA3)

death_data <- death_data%>%dplyr::mutate(ALPHA3 = gsub("OTHERS","ZZZ", ALPHA3))%>%
             dplyr::select(-country)%>%dplyr::arrange(ALPHA3)

active_data_worldo <- active_data_worldo%>%dplyr::mutate(ALPHA3 = gsub("OTHERS","ZZZ", ALPHA3))%>%
  dplyr::select(-country)%>%dplyr::arrange(ALPHA3)



death_data_worldo <- death_data_worldo%>%dplyr::mutate(ALPHA3 = gsub("OTHERS","ZZZ", ALPHA3))%>%
  dplyr::select(-country)%>%dplyr::arrange(ALPHA3)

########Make all NA become ZZZ
TMP = active_data_worldo
names1 = TMP$ALPHA3
for(i in 1:length(names1)){
  tmp =  names1[i]
  if (is.na(tmp)==1){
    names1[i] = "ZZZ"
  }
}

active_data_worldo$ALPHA3 = names1 


##########
TMP1 = death_data_worldo
names2 = TMP1$ALPHA3
for(i in 1:length(names2)){
  tmp =  names2[i]
  if (is.na(tmp)==1){
    names2[i] = "ZZZ"
  }
}

death_data_worldo$ALPHA3 = names2

####################################

##Extract countries name with more than 100 confirmed cases at the cutoffday
tt = dplyr::filter(active_data, active_data[,cutoffday]>500)  
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
###########
active_names1 = active_data_worldo[,"ALPHA3"]

index1 = c()
for(i in 1:nrow(active_names1)){
  tmp = as.character(active_names1[i,])
  if( tmp %in% countriesvector){
    index1=c(index1,i)
  }
  
}

active_data1_worldo = active_data_worldo[index1,]
death_data1_worldo = death_data_worldo[index1,]


################################
active_data1 = active_data1%>%dplyr::left_join(population_data)%>%dplyr::relocate(ALPHA3,population)

recovered_data1 = recovered_data1%>%dplyr::left_join(population_data)%>%dplyr::relocate(ALPHA3,population)

death_data1 = death_data1%>%dplyr::left_join(population_data)%>%dplyr::relocate(ALPHA3,population)

active_data1_worldo = active_data1_worldo%>%dplyr::left_join(population_data)%>%dplyr::relocate(ALPHA3,population)

death_data1_worldo = death_data1_worldo%>%dplyr::left_join(population_data)%>%dplyr::relocate(ALPHA3,population)



##########define a function collecting A,R,D
ncolumn = dim(active_data1)[2]
duration = ncolumn - 2
datalist = list()
numbercountrieshighcases =  length(countriesvector)

datafunction = function(country){
  matrix = matrix(0,nrow = duration,8)
  confirmed = active_data1[country,3:ncolumn] 
  tmp1 = unlist(confirmed)
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
datalist[["ZZZ"]] = matrix(0,nrow = duration,8)

###Adding name to keep track
datalist[["country"]] = c(countriesvector,"ZZZ")

###Adding corresponding population, ZZZ corresponding to 0
population = active_data1[,"population"]
population =  unlist(population)
population = as.vector(population)
datalist[["population"]] = c(population,0)
#########################################3
countrynames = datalist$country
numbercountries1 = length(countrynames) -1

for(i in 1:numbercountries1){
  

  tmp = countrynames[i]
  index1  = which(active_data1_worldo$ALPHA3==tmp)
  index2  = which(death_data1_worldo$ALPHA3==tmp)
  
  if(min(length(index1), length(index2))>0){
  
  
  activeworldo = active_data1_worldo[index1,3:ncolumn]
  tmp1 = unlist(activeworldo)
  tmp1 = as.vector(tmp1)
  datalist[[i]][,4] = tmp1
  
  deathworldo = death_data1_worldo[index2,3:ncolumn]
  tmp2 = unlist(deathworldo)
  tmp2 = as.vector(tmp2)
  datalist[[i]][,5] = tmp2}
  
  
}

coviddata1 = datalist

##Round1. Confirmed cases need to be fixed , JHK
indexc = c()
for(i in 1:numbercountries1){
  tt = datalist[[i]]
  t1 = which(diff(tt[,1],1)<0)
  if(length(t1)>0){
    indexc = c(indexc,i)
  }
  }
####Round1. Death cases need to be fixed , JHK
indexd = c()
for(i in 1:numbercountries1){
  tt = datalist[[i]]
  t1 = which(diff(tt[,3],1)<0)
  if(length(t1)>0){
    indexd = c(indexd,i)
  }
}  
#######Round1. Recover cases need to be fixed , JHK
indexr = c()
for(i in 1:numbercountries1){
  tt = datalist[[i]]
  t1 = which(diff(tt[,2],1)<0)
  if(length(t1)>0){
    indexr = c(indexr,i)
  }
}  

##########Fixing 1st round confirmed cases based on worldometers
for ( ic in 1:length(indexc)){
k = indexc[ic]
tmp = datalist[[k]]
ts = tmp[,1]
ts1 = tmp[,4]
lts1 = length(ts)-1

for(i in 1:lts1){
  
  i1 = i+1
  d = ts[i1] - ts[i]
  d1 = ts1[i1] - ts1[i]
  
  if(d<0&& d1>=0){
    ts[i1] = ts[i] + d1
  }
}

datalist[[k]][,1] = ts
}
##########Fixing 1st round death cases based on worldometers
for ( id in 1:length(indexd)){
  k = indexd[id]
  tmp = datalist[[k]]
  ts = tmp[,3]
  ts1 = tmp[,5]
  lts1 = length(ts)-1
  for(i in 1:lts1){
    
    i1 = i+1
    d = ts[i1] - ts[i]
    d1 = ts1[i1] - ts1[i]
    
    if(d<0 && d1>=0){
      ts[i1] = ts[i] + d1
    }
  }
  
  datalist[[k]][,3] = ts
}

#####################
##Round2. Confirmed cases need to be fixed , JHK
indexc1 = c()
for(i in 1:numbercountries1){
  tt = datalist[[i]]
  t1 = which(diff(tt[,1],1)<0)
  if(length(t1)>0){
    indexc1 = c(indexc1,i)
  }
}
####Round2. Death cases need to be fixed , JHK
indexd1 = c()
for(i in 1:numbercountries1){
  
  tt = datalist[[i]]
  t1 = which(diff(tt[,3],1)<0)
  if(length(t1)>0){
    indexd1 = c(indexd1,i)
  }
} 
###############

##########Fixing 2nd round confirmed cases based on the data itself, average middle for big small big gaps 
if(length(indexc1)>0){
for ( ic in 1:length(indexc1)){
 
  k = indexc1[ic]
  tmp = datalist[[k]]
  ts = tmp[,1]
  
  lts1 = length(ts)-1
  for(i in 1:lts1){
    

    i1 = i+1
    d = ts[i1] - ts[i]
    i2 = i-1
    d1 = round((ts[i1] + ts[i2])/2,digits=0)
   
    if(d<0&&d1>=0){
      ts[i] =  d1
    }
  }
  
  datalist[[k]][,1] = ts
}

indexc2 = c()
for(i in 1:numbercountries1){
  tt = datalist[[i]]
  t1 = which(diff(tt[,1],1)<0)
  if(length(t1)>0){
    indexc2 = c(indexc2,i)
  }
}
} else{indexc2=c()}
############Fixing 3rd round confirmed cases based on the data itself, average sliding for the tough guys 
if(length(indexc2)>0){

for ( ic in 1:length(indexc2)){
  
  k = indexc2[ic]
  tmp = datalist[[k]]
  ts = tmp[,1]
  lts1 = length(ts)-1
  d1 = diff(ts,1)
  k2 = min(which(d1<0))
  k2a=k2+1
  p = ts - ts[k2]
  k3 = min(which(p>0))
  k3a =k3-1
  d = round((ts[k3]-ts[k2])/(k3-k2+1), digits=0)
  for(i in k2a:k3a){
    ts[i] = ts[k2]+d
  }
 
 
  datalist[[k]][,1] = ts
}

indexc3 = c()
for(i in 1:numbercountries1){
  tt = datalist[[i]]
  t1 = which(diff(tt[,1],1)<0)
  if(length(t1)>0){
    indexc3 = c(indexc3,i)
  }
}

} else{indexc3=c()}

indexc3

##########Fixing 2nd round confirmed cases based on the data itself, average middle for big small big gaps 

if(length(indexd1)>0){
for ( id in 1:length(indexd1)){
  
  k = indexd1[id]
  tmp = datalist[[k]]
  ts = tmp[,3]
  lts1 = length(ts)-1
  for(i in 1:lts1){
    
    
    i1 = i+1
    d = ts[i1] - ts[i]
    i2 = i-1
    d1 = round((ts[i1] + ts[i2])/2,digits=0)
    
    if(d<0&&d1>=0){
      ts[i] =  d1
    }
  }
  
  datalist[[k]][,3] = ts
}
######
indexd2 = c()
for(i in 1:numbercountries1){
  tt = datalist[[i]]
  t1 = which(diff(tt[,3],1)<0)
  if(length(t1)>0){
    indexd2 = c(indexd2,i)
  }
}
#####
} else{indexd2=c()}
############Fixing 3rd round confirmed cases based on the data itself, average sliding for the tough guys 
if(length(indexd2)>0){
  for ( id in 1:length(indexd2)){
  
  k = indexd2[id]
  tmp = datalist[[k]]
  ts = tmp[,3]
  lts1 = length(ts)-1
  d1 = diff(ts,1)
  k2 = min(which(d1<0))
  k2a=k2+1
  p = ts - ts[k2]
  k3 = min(which(p>0))
  k3a =k3-1
  d = round((ts[k3]-ts[k2])/(k3-k2+1), digits=0)
  for(i in k2a:k3a){
    ts[i] = ts[k2]+d
  }
  
  
  datalist[[k]][,3] = ts
}

indexd3 = c()
for(i in 1:numbercountries1){
  tt = datalist[[i]]
  t1 = which(diff(tt[,3],1)<0)
  if(length(t1)>0){
    indexd3 = c(indexd3,i)
  }
}
} else{indexd3=c()}


####################################3
coviddata = datalist

tmpn = length(coviddata) - 2
for (i in 1:tmpn){
  coviddata1[[i]][,6:8] = coviddata[[i]][,1:3]
}

coviddataJanJune20 = coviddata1

#save(coviddataJanJune20 , file ="coviddataJanJune20.Rdata")
############################MAKE FLIGHT MATRIX BECOME COUNTRIES BELONG TO THE LIST TO OTHERS####
####using flight data 2020
#CHANGE ALL NAMES NOT BELONG TO countries vector

flightdata1 = c(flightdataJanFeb2020, flightdataMarch2020, flightdataApril2020,  flightdataMay2020, flightdataJune2020)

flightdataJanJune20  = list()

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


####Idea fixing is first move things not common to the corner then process and match back

###Need to think is this worth to fix this
# countriesvector1 = intersect(countriesvector, flight_names)
# countriesvectordiff = setdiff(countriesvector, flight_names)
# 
# ###Travelers between high cases countries
# dat2 = dat[countriesvector1,countriesvector1]
# tmp = matrix(0,length(countriesvector),length(countriesvector))

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

flightdataJanJune20[[d]] = round(flightrearrange, digits=0)


}  #end loop for d

#save(flightdataJanJune20, file="flightdataJanJune20.Rdata")

###########################MAKE FLIGHT MATRIX BECOME COUNTRIES BELONG TO THE LIST TO OTHERS####
####using flight data 2020
#CHANGE ALL NAMES NOT BELONG TO countries vector

flightdata1 = c(flightdataJanFeb2020, flightdataMarch2019, flightdataApril2019,  flightdataMay2019, flightdataJune2019)

flightdataJanJune20mixed19  = list()

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
  
  
  ####Idea fixing is first move things not common to the corner then process and match back
  
  ###Need to think is this worth to fix this
  # countriesvector1 = intersect(countriesvector, flight_names)
  # countriesvectordiff = setdiff(countriesvector, flight_names)
  # 
  # ###Travelers between high cases countries
  # dat2 = dat[countriesvector1,countriesvector1]
  # tmp = matrix(0,length(countriesvector),length(countriesvector))
  
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
  
  flightdataJanJune20mixed19[[d]] = round(flightrearrange, digits=0)
  
  
}  #end loop for d

save(flightdataJanJune20mixed19, file="flightdataJanJune20mixed19.Rdata")


