% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/totaltravelout_regulated.R
\name{totaltravelout_regulated}
\alias{totaltravelout_regulated}
\title{This function gives a list of total travel data allowed from 1 country to
another for a given travel policy.}
\usage{
totaltravelout_regulated(traveloutdat, policymatrix, P)
}
\arguments{
\item{traveloutdat}{is the data of travel out from each country without pandemic}

\item{policymatrix}{is a where each row is policy travel allow,
row range from 11,12,..1n, 21,22,...,2n, ..., n1,n2,...,nn
In each row, we have travel rate allow for each day during the duration of travel
in theory rows 11, 22, ..,nn should be 0, but can also be any values as convenience due to
the code itself can take care this fact}

\item{P}{vector population of countries in order as travel data}
}
\value{
A list of number travelers can go from 1 country to another each day during the duration
}
\description{
This function gives a list of total travel data allowed from 1 country to
another for a given travel policy.
}
\examples{
\dontrun{
set.seed(1)
P1 = 10^7
P2 = 3*10^6
P3 = 2*10^6
P = c(P1, P2, P3) #population of 3 countries
traveloutdat = travelout_3dat
numberCountries = ncol(traveloutdat)
#create a policy matrix where each row is policy travel allow,
#row range from 11,12,..1n, 21,22,...,2n, ..., n1,n2,...,nn
#In each row, we have travel rate allow for each day during the duration of travel
Here we ignored the fact that at row 11, 22,...,nn the rate should be 0, since the travel out
divide in the code already take care the 0 part
total = nrow(traveloutdat)*numbercountries^2
policy = matrix(runif(total, 0, 1), ncol = nrow(traveloutdat), nrow = numbercountries^2)
totaltravelout_regulated(traveloutdat, policy, P)
}
}
