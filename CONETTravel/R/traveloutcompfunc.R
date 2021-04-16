#' This function gives travel out compartment in 6 categories from 1 country to another.
#' @param total total travelers out
#' @param country a list include data of A(active confirmed), R(Recovered confirmed),D (confirmed deceased) sequences, the dynamic population pop, infect.
#' @return  The average realization of the country during the period
#' @examples
#' \dontrun{
#'library(CONETTravel)
#'library(rlist)
#'data = datas_3travel[[1]][,3:5]
#'x = c(9999620, 250,130,0,0,0)
#'P= sum(x)
#'durationtravel = 84
#'w21= travelout_3dat[1,1]/(travelout_3dat[1,1] + travelout_3dat[1,3])
#'w31= travelout_3dat[1,1]/(travelout_3dat[1,1] + travelout_3dat[1,2])
#'out21 = travelout_3dat[,2]
#'out31 = travelout_3dat[,3]
#'total_in = out21*w21 + out31*w31
#'total_out = travelout_3dat[,1]
#'country = list(P=P, total_in =  total_in, total_out= total_out,durationtravel=durationtravel)
#'populationdynamic = populationdynamicfunc(country)
#'infect = c(infectfunc(data, x),0)
#'travelin_compartments = matrix(0,nrow=country$durationtravel,6)
#'S_in = round(total_in*.97,digits=0)
#'I_in = round(total_in*.03,digits=0)
#'travelin_compartments[,1] = S_in
#'travelin_compartments[,2] = I_in
#'country =list.append(country,data=data, travelin_compartments = travelin_compartments,
#'pop=populationdynamic,infect= infect, x_ini = x)
#'total = total_out
#'traveloutcompfunc(total, country)}
#' @export

traveloutcompfunc = function(total,country)
{
  mydat = country$data
  pop = country$pop
  infect = country$infect
#######

  drecover = diff(mydat[,2], differences = 1)

  active = mydat[,1][-length(mydat[,1])]

  index0 = which(active<=.1)

  if (length(index0)!=0){
    drecover = drecover[-index0]
    active = active[-index0]
  }

  betas = drecover/active


  Ut = rowSums(mydat)
  #Difference Ut
  dUt = c(diff(Ut, differences = 1),0)
  ########
  dRus = rep(0,nrow(mydat))

  for (i in 1:(nrow(mydat)-1)){

    dRus[i] = infect[i]*median(betas)

  }


  Ruseq1 = cumsum(dRus)

  Ruseq = c(0,Ruseq1[-length(Ruseq1)])


  Sseq = pop - Ut - Ruseq - infect



  compartments = cbind(Sseq, infect, Ruseq)


  update = matrix(0,ncol=6,nrow=nrow(mydat))

  for (i in 2:nrow(mydat)){
    i1= i-1
    x = compartments[i1,]

    ####
    if(min(x)>0){
      update[i,] = c(total[i]*x[1]/(x[1]+x[2]), total[i]*x[2]/(x[1]+x[2]), 0, 0, 0, 0)
    } else{
      update[i,] = c(total[i], 0, 0, 0, 0, 0)
    }
  }
  return(update)
}

