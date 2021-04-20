#' This function gives a rough estimation of of the recovered rate beta
#'  and used to improve the estimation of beta.
#' @param data data of A(active confirmed), R (recovered confirmed),D(confirmed deceased)
#' @return  a sequence of estimating values of beta
#' @examples
#' \dontrun{
#' library(CONETTravel)
#' data = datas_3travel[[1]][,3:5]
#' betafunc(data)
#' }


#' @export
betafunc = function(data){
  drecover_sim = data[,2] - c(0,data[,2][-length(data[,2])])
  active_sim = c(1,data[,1][-length(data[,1])])

  index5 = which(active_sim==0)
  index6 = c(1,index5)
  drecover_sim = drecover_sim[-index6]
  active_sim = active_sim[-index6]
  betas_sim = drecover_sim/active_sim
  ########This help to avoid empty array
  if (length(betas_sim) == 0){betas_sim = 10^8}
  ###########


  return(betas_sim)

}
