#' This function gives a rough estimation of of the death rate delta
#'  and used to improve the estimation of delta.
#' @param data data of A (active confirmed), R(Recovered confirmed), D(Confirmed Deceased)
#' @return  a sequence of estimating values of delta
#' @examples
#' \dontrun{
#' library(CONETTravel)
#' data = datas_3travel[[1]][,3:5]
#' deltafunc(data)
#' }

#' @export
deltafunc = function(data){
  ddeath_sim = data[,3] - c(0,data[,3][-length(data[,3])])
  active_sim = c(1,data[,1][-length(data[,1])])

  index7 = which(active_sim==0)
  index8 = c(1,index7)
  ddeath_sim =  ddeath_sim[-index8]
  active_sim = active_sim[-index8]

  deltas_sim = ddeath_sim/active_sim
  ########This help to avoid empty array
  if (length(deltas_sim) == 0){deltas_sim = 10^8}
  ###########
  return(deltas_sim)

}
