#' This function estimate number infected in each country for a given data set
#'  under our model assumption.
#' @param data A (active confirmed), R(Recovered confirmed), D(Confirmed Deceased)
#' @param x six initial states
#' @return  a sequence of estimating values of delta
#' @examples
#' \dontrun{
#' library(CONETTravel)
#' data = datas_3travel[[1]][,3:5]
#' x = c(9999620, 250,130,0,0,0)
#' infectfunc(data, x)
#' }
#' @export
infectfunc = function(data,x){
  U = t(t(rowSums(data)))
  diff1 = diff(U, differences = 1)
  diffa = diff1[-1]
  diffb = diff1[-length(diff1)]
  infect = rep(0,length(diff1))

  ###############
  infect[1] = x[2]

  for (i in 2:length(diff1)){
    i1 = i-1
    tmp = infect[i1]


    if (diffb[i1]>.4){
      infect[i] = tmp*diffa[i1]/diffb[i1]

    } else{
      infect[i] = diff1[i]


    }

  }

  return (infect)

}

