#' This function gives a rough estimation of of the transmission rate alpha
#'  and used to improve the estimation of alpha.
#' @param data data of A (active confirmed), R(Recovered confirmed), D(Confirmed Deceased)
#' @param country a list include betas, infect,
#' populationdynamic, travelin_compartments, travelout_compartments
#'
#' @return  a sequence of estimating values of alpha
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
#' travelin_compartments = matrix(0,nrow=country$durationtravel,6)
#' travelout_compartments = matrix(0,nrow=country$durationtravel,6)
#' S_in = round(total_in*.97,digits=0)
#' I_in = round(total_in*.03,digits=0)
#' travelin_compartments[,1] = S_in
#' travelin_compartments[,2] = I_in
#' S_out = round(total_out*.99,digits=0)
#' I_out = round(total_out*.01,digits=0)
#' travelout_compartments[,1] = S_out
#' travelout_compartments[,2] = I_out
#' country =list.append(country, betas = betas, populationdynamic=populationdynamic,
#'                      infect= infect, travelout_compartments = travelout_compartments,
#'                                           travelin_compartments = travelin_compartments)
#'alphafunc(data,country)
#' }
#' @export
alphafunc = function(data,country){


  beta = median(country$betas)
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


  Sseq = country$populationdynamic - Ut - Ruseq - infected


  SIseq = Sseq*infected

  #########Construct dSseq
  Sseq_p = rep(0, nrow(data))

  for (i in 1:nrow(data)){
    Sseq_p[i] = Sseq[i] - country$travelin_compartments[i,1] + country$travelout_compartments[i,1]
  }
  ######S+ and Spre sequence for alpha
  Splus_al = Sseq[-c(nrow(data)-1, nrow(data))]

  Spre_al = Sseq_p[-c(1,nrow(data))]

  dSseq_al = Splus_al - Spre_al

  SIseq_al = SIseq[-c(nrow(data)-1, nrow(data))]

  pop_al = country$populationdynamic[-c(nrow(data)-1, nrow(data))]


  #####Remove 0s in the denom
  index1 = which(SIseq_al<=.5)
  ####Remove negative difference
  index2 = which(dSseq_al<=0.5)
  ###########
  index = c(index1,index2)
  index = unique(index)
  index = sort(index)

  ######
  if (length(index!=0)){
    dSseq_al = dSseq_al[-index]
    SIseq_al = SIseq_al[-index]
    pop_al = pop_al[-index]
  }


  alphas_sim = dSseq_al/SIseq_al*pop_al

  ########This help to avoid empty array
  if (length(alphas_sim) == 0){alphas_sim = 10^8}
  ###########


  return(alphas_sim)

}


