#' This function gives the average realization of 1 given country for a given
#' A (active confirmed), R(Recovered confirmed), D(Confirmed Deceased) sequences.
#' @param country is a list include a data matrix with A (active confirmed), R(Recovered confirmed),
#'  D(Confirmed Deceased), populationdynamic,and infect
#' @param theta parameter
#' @examples
#' \dontrun{
#' library(CONETTravel)
#' library(rlist)
#' data = datas_3travel[[1]][,3:5]
#' betas = betafunc(data)
#' x = c(9999620, 250,130,0,0,0)
#' P= sum(x)
#' durationtravel = 84
#' w21= travelout_3dat[1,1]/(travelout_3dat[1,1] + travelout_3dat[1,3])
#' w31= travelout_3dat[1,1]/(travelout_3dat[1,1] + travelout_3dat[1,2])
#' out21 = travelout_3dat[,2]
#' out31 = travelout_3dat[,3]
#' total_in = out21*w21 + out31*w31
#' total_out = travelout_3dat[,1]
#' country = list(P =P, total_in = total_in, total_out= total_out,durationtravel=durationtravel)
#' populationdynamic = populationdynamicfunc(country)
#' infect = c(infectfunc(data, x),0)
#' country =list.append(country, data=data, populationdynamic=populationdynamic,
#'                     infect= infect)
#'                     theta = thetas_3travel[[1]][1:6]
#'                     recoverfunc(country,theta)
#' }
#' @export
recoverfunc = function(country, theta){

  data = country$data
  pop = country$populationdynamic

  beta = theta[3]

  Ut = rowSums(data)
  #Difference Ut
  dUt = c(diff(Ut, differences = 1),0)
  ########

  recovermat = matrix(0,nrow = nrow(data),ncol=6)

  infected = country$infect

  dRus = rep(0,nrow(data))

  for (i in 1:nrow(data)){

    dRus[i] = dUt[i]*beta

  }

  Ruseq1 = cumsum(dRus)

  Ruseq = c(0,Ruseq1[-length(Ruseq1)])


  Sseq = pop - Ut - Ruseq - infected

  recovermat = data.frame(Sseq, infected, data[,1], data[,2], data[,3], Ruseq)



  return(recovermat)

}
