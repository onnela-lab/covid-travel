#' This function gives a list of total travel data allowed from 1 country to
#' another for a given travel policy.
#' @param traveloutdat is the data of travel out from each country without pandemic
#' @param P vector population of countries in order as travel data
#' @param policymatrix is a where each row is policy travel allow,
#' row range from 11,12,..1n, 21,22,...,2n, ..., n1,n2,...,nn
#' In each row, we have travel rate allow for each day during the duration of travel
#' in theory rows 11, 22, ..,nn should be 0, but can also be any values as convenience due to
#' the code itself can take care this fact
#' @return  A list of number travelers can go from 1 country to another each day during the duration
#' @examples
#' \dontrun{
#' set.seed(1)
#' P1 = 10^7
#' P2 = 3*10^6
#' P3 = 2*10^6
#' P = c(P1, P2, P3) #population of 3 countries
#' traveloutdat = travelout_3dat
#' numberCountries = ncol(traveloutdat)
#' #create a policy matrix where each row is policy travel allow,
#' #row range from 11,12,..1n, 21,22,...,2n, ..., n1,n2,...,nn
#' #In each row, we have travel rate allow for each day during the duration of travel
#' Here we ignored the fact that at row 11, 22,...,nn the rate should be 0, since the travel out
#' divide in the code already take care the 0 part
#' total = nrow(traveloutdat)*numbercountries^2
#' policy = matrix(runif(total, 0, 1), ncol = nrow(traveloutdat), nrow = numbercountries^2)
#' totaltravelout_regulated(traveloutdat, policy, P)
#' }


#' @import matrixcalc
#' @export


totaltravelout_regulated =  function(traveloutdat,policymatrix, P){
##Construct a divide list travel out
numberCountries = ncol(traveloutdat)

traveloutDivide = list()

for (j in 1:nrow(traveloutdat)){
  tmp = traveloutdat[j,]
  ##CREATE MATRIX TRAVELOUT EACH ROW
  travelmat = matrix(0,nrow = numberCountries, ncol = numberCountries)
  for(i in 1:numberCountries){
    ones = rep(1,numberCountries)
    travelmat[i,] = kronecker(ones,as.numeric(tmp[i]) )
  }
  ##CREATE THE PROPORTION OF POPULATION MATRIX
  ones = t(t(rep(1,numberCountries)))
  popmat = kronecker(ones,t(P))
  ##
  denommat = matrix(0,nrow = numberCountries, ncol = numberCountries)
  for(i in 1:numberCountries){
    tmp = sum(P[-i])
    ones = rep(1,numberCountries)
    denommat[i,] = kronecker(ones,as.numeric(tmp) )
  }
  ##Take product for proportion
  promat = popmat*denommat^(-1)
  diag(promat) = 0
  # Take Hadamard product for actual number
  traveloutDivide[[j]] = hadamard.prod(promat, travelmat)
}


##Matrix list of traveling
durationtravel = nrow(traveloutdat)

policymatlist = list()

for(i in 1:durationtravel){
  tmp =  matrix(0,nrow = numberCountries, ncol = numberCountries)
  tmp = matrix(policymatrix[,i], nrow = numberCountries, ncol = numberCountries)
  policymatlist[[i]] = tmp

}

############Matrix travel when percentage of travel in regulated
traveloutDivideRegulated = list()
for (i in 1:nrow(traveloutdat)){

  tmp1 = as.matrix(traveloutDivide[[i]])
  tmp2 = as.matrix(policymatlist[[i]])

  traveloutDivideRegulated[[i]] = round(hadamard.prod(tmp1, tmp2), digits=0)
}

return(traveloutDivideRegulated)
}
######
