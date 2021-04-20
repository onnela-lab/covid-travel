#--- script for cleaning the Opensky data, to rp
setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/Realdata/Datapreprocessing")
library(dplyr)
library(here)
library(countrycode)
library(tidyverse)
library(lubridate)



open_sky1   <- read.csv("opensky_january2019.csv") #need to change
open_sky2   <- read.csv("opensky_january2020.csv")  #need to change




#--- reading in airport lookup table
airport_lookup <- read.delim(here( "airports.dat"), sep = ",") %>%
  dplyr::select(country = Country, three_letter_code, four_letter_code)

#--- cleaning flight data, keeping only the total flights between countries

clean_flight_data <- function(open_sky_data)
{
  #########Clean date 
  
  
  open_sky_data[,"day"] = gsub("00:", "", open_sky_data[, "day"])
  open_sky_data[,"day"] = gsub("00", "", open_sky_data[, "day"])
  open_sky_data[,"day"] = str_remove_all(open_sky_data[, "day"], "[+]")
  
  
  
  ############
  flights_by_origin <- open_sky_data %>% 
    dplyr::mutate(date = lubridate::ymd(day))%>%
    dplyr::select(date, four_letter_code = origin) %>%
    dplyr::left_join(airport_lookup) %>%
    #tidyr::drop_na() %>%
    dplyr::rename(origin_country = country,
                  origin_code = four_letter_code) %>%
    dplyr::select(-three_letter_code)
  
  flights_by_destination <- open_sky_data %>% 
    dplyr::mutate(date = lubridate::ymd(day))%>%
    dplyr::select(date, four_letter_code = destination) %>%
    dplyr::left_join(airport_lookup) %>%
    #tidyr::drop_na() %>%
    dplyr::rename(destination_country = country,
                  destination_code = four_letter_code) %>%
    dplyr::select(-three_letter_code)
  
  
  flights_by_origin$destination_code <- flights_by_destination$destination_code
  flights_by_origin$destination_country <- flights_by_destination$destination_country
  
  flight_data_out <- flights_by_origin %>%
    tidyr::drop_na() %>%
    dplyr::select(date, origin_code, destination_code, origin_country, destination_country) %>%
    dplyr::rowwise() %>%
    dplyr::filter(origin_country != destination_country)
  
  return(flight_data_out)
  
}  



################Define the reduction function

reductionflights  = function(data1, data2){
  
flight_data1  <- clean_flight_data(data1)

flight_data2  <- clean_flight_data(data2)




#--- summarising data1
total_flights1 <- flight_data1 %>%
  dplyr::filter(origin_country != destination_country) %>%
  dplyr::group_by(origin_country, destination_country) %>%
  dplyr::summarise(total_flights1 = dplyr::n()) %>%
  dplyr::mutate(origin_country_iso_code = countrycode::countrycode(origin_country, "country.name", "iso3c"),
                destination_country_iso_code = countrycode::countrycode(destination_country, "country.name", "iso3c")) %>%
  dplyr::ungroup() %>%
  dplyr::select(origin_country_iso_code, destination_country_iso_code, total_flights1)

#--- summarising data2
total_flights2 <- flight_data2 %>%
  dplyr::filter(origin_country != destination_country) %>%
  dplyr::group_by(origin_country, destination_country) %>%
  dplyr::summarise(total_flights2 = dplyr::n()) %>%
  dplyr::mutate(origin_country_iso_code = countrycode::countrycode(origin_country, "country.name", "iso3c"),
                destination_country_iso_code = countrycode::countrycode(destination_country, "country.name", "iso3c")) %>%
  dplyr::ungroup() %>%
  dplyr::select(origin_country_iso_code, destination_country_iso_code, total_flights2)

#--- putting data1 and 2 together and calculating all scaling factors, capped at 1 where 
#--- for some reason more flights are occurring
total_flights12 <- total_flights1 %>% 
  dplyr::ungroup() %>%
  dplyr::left_join(total_flights2) %>%
  dplyr::mutate(scaling_factor = total_flights2/total_flights1) %>%
  dplyr::mutate(scaling_factor = dplyr::case_when(scaling_factor > 1  ~ 1,
                                                  scaling_factor <= 1 ~ scaling_factor)) %>%
  dplyr::select(origin_country_iso_code, destination_country_iso_code, scaling_factor) %>%
  tidyr::drop_na()

##
return(total_flights12)
}



#--- writing the results to a .csv
reduction = reductionflights(open_sky1, open_sky2)

#write.csv(reduction , "flight_reduction_scaling_factors_march20vsmarch19.csv") # need to change



#--- calculating the mean reduction over all flight paths so that we can use the mean
#--- in the absence of an estimate. Comes out to 0.372

reduction %>%
  dplyr::ungroup() %>%
  dplyr::summarise(average = mean(scaling_factor))