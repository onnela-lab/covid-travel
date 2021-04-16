#' This function gives stochastic compartments after completed quarantine.
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
#' stochastic_postquarantine(theta0,inp)}
#' @export




stochastic_postquarantine = function(theta,inp){
  n1 = 2 + inp$durationquarantine # Adjust for 0 and 1
  status_matrix = matrix(0,nrow = n1, ncol=6)
  status_matrix[1,] = inp$ini


  #New confirmed rate, gamma, c(alpha0,alpha, beta, delta, eta, gamma)
  harzard2 = function(x,theta){
    h2 = theta[6]*x[2]
    names(h2)=c("hazard2")
    return(h2)
  }
  #New confirmed recover, beta, c(alpha0,alpha, beta, delta, eta, gamma)
  harzard3 = function(x,theta){
    h3 = theta[3]*x[3]
    names(h3)=c("hazard3")
    return(h3)
  }
  #New confirmed death, delta, c(alpha0,alpha, beta, delta, eta, gamma)
  harzard4 = function(x,theta){
    h4 = theta[4]*x[3]
    names(h4)=c("hazard4")
    return(h4)
  }
  #New unconfirmed recover, eta*beta, c(alpha0,alpha, beta, delta, eta, gamma)
  harzard5 = function(x,theta){
    h5 = theta[5]*theta[3]*x[2]
    names(h5)=c("hazard5")
    return(h5)
  }

  for (i in 2:n1){

    x = status_matrix[(i-1),]


    #
    y2 =  rpois(1, harzard2(x,theta))
    #
    y3 =  rpois(1, harzard3(x,theta))
    #
    y4 =  rpois(1, harzard4(x,theta))
    #
    y5 =  rpois(1, harzard5(x,theta))

    #######Infect
    if(x[2] - y2 - y5 >= 0){
      x[2] = x[2] - y2 - y5} else{
        y2 =  x[2]  - y5
        x[2] =0
      }

    if(y2 <0){
      y2 = 0
      y5 = x[2]
    }

    #######Active
    if(y2-y3-y4+x[3] >= 0){
      x[3] = x[3] + y2 - y3 - y4} else{
        y3 = y2-y4+x[3]
        x[3] =0
      }

    if(y3 <0){
      y3 = 0
      y4 = x[3] + y2
    }

    #Recover
    x[4] =x[4] + y3
    #Death
    x[5]= x[5] + y4
    #Recover Unconfirmed
    x[6]= x[6] + y5





    status_matrix[i,] = x
  }
  n2 = n1-1
  lastdayquarantine = status_matrix[n2, ]
  lastdayquarantine = round(lastdayquarantine,digits=0)
  return(lastdayquarantine)
}
