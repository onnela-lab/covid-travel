#' This function gives stochastic compartments evolved process during the
#' quarantine period.
#' Outputs include travelers status done
#'  quarantine and new active confirmed during the quarantine period.
#' @param theta parameter
#' @param inp is a list include durationquarantine : number of days and ini: initial compartments of arrivals
#' @return  The stochastic status of arrivals compartments after complete quarantine
#' @importFrom stats rpois
#' @examples
#' \dontrun{## Initial Condition
#' S1 = 900
#' I1 = 100
#' A1 = 0
#' x1 = c(S1,I1,A1,0,0,0)#
#' days= 14
#' inp = list(durationquarantine = days, ini = x1)
#' theta0 = c(0,0,1/14,3/100,1,.5)
#' stochastic_postquarantine_separate(theta0,inp)}
#' @export



stochastic_postquarantine_separate = function (theta, inp)
{

  n1 = 2 + inp$durationquarantine
  status_matrix = matrix(0, nrow = n1, ncol = 6)
  status_matrix1 = matrix(0, nrow = n1, ncol = 6)
  activeeachday_matrix = matrix(0, nrow = n1, ncol = 6)

  status_matrix[1, ] = inp$ini
  status_matrix1[1, ] = inp$ini


  harzard2 = function(x, theta) {
    h2 = theta[6] * x[2]
    names(h2) = c("hazard2")
    return(h2)
  }
  harzard3 = function(x, theta) {
    h3 = theta[3] * x[3]
    names(h3) = c("hazard3")
    return(h3)
  }
  harzard4 = function(x, theta) {
    h4 = theta[4] * x[3]
    names(h4) = c("hazard4")
    return(h4)
  }
  harzard5 = function(x, theta) {
    h5 = theta[5] * theta[3] * x[2]
    names(h5) = c("hazard5")
    return(h5)
  }
  for (i in 2:n1) {


    x = status_matrix[(i - 1), ]
    x1 = status_matrix1[(i - 1), ]
    y3 = rpois(1, harzard3(x1, theta)) # active to recover
    y4 = rpois(1, harzard4(x1, theta)) # active to death

    y2 = rpois(1, harzard2(x, theta))
    y5 = rpois(1, harzard5(x, theta))

    #####update active confirmed each day
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
    x[3] = y2
    x[4] = 0
    x[5] = 0
    x[6] = x[6] + y5
    status_matrix[i, ] = x
    ###################
    #####update all status each day statusmatrix1
    if (x1[2] - y2 - y5 >= 0) {
      x1[2] = x1[2] - y2 - y5
    }
    else {
      y2 = x1[2] - y5
      x1[2] = 0
    }
    if (y2 < 0) {
      y2 = 0
      y5 = x1[2]
    }
    if (y2 - y3 - y4 + x1[3] >= 0) {
      x1[3] = x1[3] + y2 - y3 - y4
    }
    else {
      y3 = y2 - y4 + x1[3]
      x1[3] = 0
    }
    if (y3 < 0) {
      y3 = 0
      y4 = x1[3] + y2
    }
    x1[4] = x1[4] + y3
    x1[5] = x1[5] + y4

    x1[6] = x1[6] + y5
    status_matrix1[i, ] = x1


  }
  n2 = n1 - 1
  tmp = status_matrix1[n2, ]

  activeeachday_matrix[, 3] = status_matrix[, 3]

  activeeachday_matrix[n2, 3] = 0  ###hide the status of active confirmed in last day due to jump in to donequarantine
  return(list(donequarantine = tmp, activeconfirm_eachday = activeeachday_matrix[1:n2,
  ]))
}
