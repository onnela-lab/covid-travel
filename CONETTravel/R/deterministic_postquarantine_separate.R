#' This function gives deterministic compartments of the evolving process during
#'  the quarantine period. Outputs include travelers status done
#'  quarantine and new active confirmed each day during the quarantine period.
#' @param theta parameter
#' @param inp is a list include durationquarantine : number of days and ini: initial compartments of arrival
#' @return  The average status of arrivals compartments after complete quarantine
#' @examples
#' \dontrun{## Initial Condition
#' S1 = 900
#' I1 = 100
#' A1 = 0
#' x1 = c(S1,I1,A1,0,0,0)#
#' days= 14
#' inp = list(durationquarantine = days, ini = x1)
#' theta0 = c(0,0,1/14,3/100,1,.5)
#' deterministic_postquarantine_separate(theta0,inp)}
#' @export


deterministic_postquarantine_separate = function (theta, inp) {
  n1 = 2 + inp$durationquarantine
  status_matrix = matrix(0, nrow = n1, ncol = 6)
  activeeachday_matrix = matrix(0, nrow = n1, ncol = 6)
  status_matrix[1, ] = inp$ini
  harzard2 = function(x, theta) {
    h2 = theta[6] * x[2]
    names(h2) = c("hazard2")
    return(h2)
  }


  harzard5 = function(x, theta) {
    h5 = theta[5] * theta[3] * x[2]
    names(h5) = c("hazard5")
    return(h5)
  }

  for (i in 2:n1) {
    x = status_matrix[(i - 1), ]
    y2 = harzard2(x, theta)

    y5 = harzard5(x, theta)
    if (x[2] - y2 - y5 >= 0) {
      x[2] = x[2] - y2 - y5
    }
    else {
      y2 = x[2] - y5
      x[2] = 0
    }
    if (y2 < 0) {
      y2 = 0
      y5 = x[2]
    }

    x[3] =  y2

    x[4] =  0
    x[5] = 0
    x[6] = x[6] + y5
    status_matrix[i, ] = x
  }
  n2 = n1 - 1
  tmp =  status_matrix[n2, ]# lastday categorizes
  tmp[3] = 0 #replace A by 0, since A can be added active confirm each day
  activeeachday_matrix[,3] = status_matrix[,3] #only keep the active confirm category
  return(list(donequarantine = tmp, activeconfirm_eachday =activeeachday_matrix[1:n2,] ))
}
