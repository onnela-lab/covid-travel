#' This function gives the estimation for the dynamic population.
#' @param country list include total travel in: total_in, total travel out each day,
#' total_out, durationtravel, population P
#' @return  the dynamic population of the coutry with travel take in to account
#' @examples
#' \dontrun{
#' library(CONETTravel)
#' P= 10^6
#' durationtravel = 84
#' total_in = round(runif(84,300,400),digits=0)
#' total_out = round(runif(84,300,400),digits=0)
#' country = list(P =P, total_in = total_in, total_out= total_out,durationtravel=durationtravel )
#' populationdynamicfunc(country)
#' }
#' @export

populationdynamicfunc = function(country){
  fluc1 = country$total_in - country$total_out
  pop = rep(0,country$durationtravel)
  pop[1] = country$P

  for(i in 2:country$durationtravel){
    x = pop[i-1]
    pop[i] = x+ fluc1[i] #update total population
  }
  return(pop)

}
