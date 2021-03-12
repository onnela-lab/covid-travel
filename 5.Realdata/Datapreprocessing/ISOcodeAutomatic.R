##########################################################################################
# This script has two purposes

# 1. Takes Travel data from OAGdata.csv and constructs a list of matricies. Each list is
# for a day. Each matrix has outgoing country as the row and incoming country as the column.
# The entries are the total number of seats on all flights that day for 2 pairs of countries.

# 2. Takes Covid data from time_series_covid19_confirmed_global.csv, 
# time_series_covid19_deaths_global.csv, and time_series_covid19_recovered_global.csv and 
# constructs a list. Each entry of the list is a matrix for a country. 

# First column is the active confired for the country
# Second column is the cumulative recovered
# Third column is the cumulative deaths
# Each row is a day.

# Regional data is ignored 


# Date started:  01/06/2020
# Last edited:  01/08/2020
##########################################################################################



##################
# 1
##################
#file_name <- "~/Desktop/Research/JP/Covid/Data/Data_Processing /CONETDATA/OAGdata.csv"
setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/Realdata")


file_name = paste("OAGdata.csv",sep="")


daily_flight_data <- function(file_name){
  ##################################################################
  # This function generates the travel data and outputs our list.
  # file_name:
  # location of the file with the travel data
  ##################################################################
  
  OAG <- read.csv(file_name)
  
  countries <- unique(c(as.character(OAG[, "Dep.IATA.Country.Name"]),
                        as.character(OAG[, "Arr.IATA.Country.Name"])))
  dates <- sort(as.character(unique(OAG[, "Time.series"])))
  
  # There are 60 unique dates
  d_n <- 1
  
  # There are 230 unique contries in the dataset
  c_n <- length(unique(countries))
  
  # list that maintains the matrices
  OAG_list <- list()
  
  for(i in 1:d_n){
    # intializing matrix
    OAG_list[[dates[i]]] <- matrix(0, c_n, c_n)
    colnames(OAG_list[[dates[i]]]) <- rownames(OAG_list[[dates[i]]]) <- countries
    
    # This is all the data for one day
    OAG_day_data <- OAG[which(OAG[, "Time.series"] == dates[i]), ]
    
    # Filling in the country data using the day data
    for(j in 1:c_n){
      for(k in 1:c_n){
        depart <- rownames(OAG_list[[dates[i]]])[j]
        arrive <- colnames(OAG_list[[dates[i]]])[k]
        
        # Temporary data to fill in each entry. Constantly modified in loop.
        t_n <- which(OAG_day_data[, "Dep.IATA.Country.Name"] == depart & 
                       OAG_day_data[, "Arr.IATA.Country.Name"] == arrive)
        if(length(t_n) == 0){
          # if there is no flights from country A-->B leave as 0.
          next
        }
        
        t_data <- OAG_day_data[t_n, ]
        
        # summing up the data from all the flights
        flight_total <- 0
        for(x in 1:dim(t_data)[1]){
          flight_total <- flight_total + t_data[x, "Frequency"]*t_data[x, "Seats..Total."]
        }
        
        # Placing entry in matrix
        OAG_list[[dates[i]]][j,k] <- flight_total
      }
    }
    
    print(i/d_n)
  }
  
  return(OAG_list)
}



flightdata <- daily_flight_data(file_name)

aa = flightdata[[1]]
View(aa)

save(flightdata, file="flightdata.RData")

#saveRDS(results, file = "inital_travel_results.RDS")


###########################################################
# 2
##################
file_name_active <- paste("time_series_covid19_confirmed_global.csv",sep="")

file_name_recovered <- paste("time_series_covid19_recovered_global.csv",sep="")

file_name_death <- paste("time_series_covid19_deaths_global.csv",sep="")



date_converter <- function(data){
  # This function takes the active, recovered, death data and converts the format of the date
  # names in the columns to be consistent with the travel data ie. 2020-01-23. Original format
  # must be of the form "X1.23.20" 
  # data
  # Either the recoverd, active, or death data whose column names you want to change
  
  act_dates <- colnames(data)[5:dim(data)[2]] # the dates start at the 5th column for this data
  act_dates <- substr(act_dates, 2, nchar(act_dates))
  act_dates_split <- strsplit(act_dates, ".", fixed = TRUE)
  act_dates_changed <- NULL
  for(i in 1:length(act_dates_split)){
    dte <- act_dates_split[[i]]
    n1 <- which(nchar(dte) == 1)
    if(length(n1) > 0){ # turning 1.23.20 --> 2020-01-23
      dte[n1] <- paste(0, dte[n1], sep = "")
    }
    act_dates_changed[i] <- paste("20", dte[3], "-", dte[1], "-", dte[2], sep = "")
  }
  colnames(data)[5:dim(data)[2]] <- act_dates_changed
  return(data)
}

covid_infec_data <- function(file_name_active, file_name_recovered, file_name_death){
  # This function processes our Covid data and outputs a list
  
  # file_name_active:
  # file location that contains the active Covid cases
  # file_name_recovered:
  # file location that contains the covid recovery data
  # file_name_death:
  # file location that contains the covid death data
  
  # Reading in Data
  active_data <- read.csv(file_name_active)
  recovered_data <- read.csv(file_name_recovered)
  death_data <- read.csv(file_name_death)
  
  # Changing dates in the colnames from "X2.29.20" to a "2020-02-29" format
  active_data <- date_converter(active_data)
  recovered_data <- date_converter(recovered_data)
  death_data <- date_converter(death_data)
  
  ### Processing the Covid data
  
  # There are 191 unique countries
  countries <- unique(c(as.character(active_data[, "Country.Region"]),
                        as.character(recovered_data[, "Country.Region"]),
                        as.character(death_data[, "Country.Region"])))
  c_n <- length(countries)
  
  # There are 349 unique dates
  dates <- sort(unique(c(grep("-", colnames(active_data), value = TRUE),
                         grep("-", colnames(recovered_data), value = TRUE),
                         grep("-", colnames(death_data), value = TRUE))))
  #d_n <- length(dates)
  d_n <- 68  #only 39 days based on flight data available
  # List that maintains the matrices
  Covid_list <- list()
  
  for(i in 1:c_n){
    # isolating country
    a_n <- which(active_data[, "Country.Region"] == countries[i])
    t_active_data <- active_data[a_n, ]
    
    r_n <- which(recovered_data[, "Country.Region"] == countries[i])
    t_recovered_data <- recovered_data[r_n, ]
    
    dea_n <- which(death_data[, "Country.Region"] == countries[i])
    t_death_data <- death_data[dea_n, ]
    
    # initalizing matrix
    Covid_list[[countries[i]]] <- matrix(0, d_n, 3)
    colnames(Covid_list[[countries[i]]]) <- c("active_confirmed", 
                                              "cumulative_recovered",
                                              "cumulative_deaths")
    rownames(Covid_list[[countries[i]]]) <- dates[1:d_n]
    
    # Filling in Covid data 
    for(j in 1:d_n){
      day <- dates[j]
      
      a_n_day <- which(colnames(t_active_data) == day) 
      r_n_day <- which(colnames(t_recovered_data) == day) 
      dea_n_day <- which(colnames(t_death_data) == day) 
      
      Covid_list[[countries[i]]][j, "active_confirmed"] <- sum(t_active_data[, a_n_day]) - sum(t_recovered_data[, r_n_day]) -sum(t_death_data[, dea_n_day])  
      # dates start at 5th column
      Covid_list[[countries[i]]][j, "cumulative_recovered"] <- sum(t_recovered_data[, r_n_day]) 
      Covid_list[[countries[i]]][j, "cumulative_deaths"] <- sum(t_death_data[, dea_n_day])
    }
    
    print(i/c_n)
  }
  
  return(Covid_list)
}

coviddata <- covid_infec_data(file_name_active,
                              file_name_recovered,
                              file_name_death)


coviddata[1]


index = c()
for(i in 1: length(coviddata)){
  tmp = coviddata[[i]]
  if(tmp[32,3]>0){
    print(i)
    index=c(index,i)
    
  }
}
########
countries[index]
k=10
t = index[k]

coviddata[t]


length(coviddata)
coviddata[[36]]


coviddata[[178]]

test = diff(coviddata$Australia[,3],1)

summary(test)
data = read.csv(file_name_death,header=T)
data[,2]
length(coviddata)
dim(flightdata[[1]])


tail(coviddata$US,9)
cc = tail(coviddata$US,11)
cc[,1]
diff(cc[,1],1)

diff(tt,1)

save(coviddata, file="coviddata.RData")



tt1 = names(coviddata)

tt2 = flightdata[[1]]
tt2 = rownames(tt2)

TT = intersect(tt1,tt2)
variable = TT
write.csv(TT, "commoncountries.csv")
data1 = read.csv("commoncountries.csv")
iso_data = read.csv("countriesISOcode.csv")
iso_data[,3]

for(i in 1:length(variable)){
  var1 = variable[i]
  n0 <- which(iso_data[, "Country"] == var1)
  if(length(n0) == 0){
    next
  }else{
    variable[i] <- as.character(iso_data[n0, "Alpha.3.code"])
  }
}


variable[1]
variable[2]

variable1 <- gsub("\"", "",variable)
countries = data.frame(TT,TT,variable1)
write.csv(countries, "commoncountries.csv")


