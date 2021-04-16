#' This function gives a stochastic realization for a given parameter and
#' initial condition.
#' @param theta parameter
#' @param inp is a list include duration : number of days and ini: initial compartments of the country
#' @return  The average realization of the country during the period
#' @importFrom stats rpois
#' @examples
#' \dontrun{## Initial Condition
#' P1 = 10^7
#' I1 = 250
#' A1 = 130
#' S1 = P1 - I1 - A1
#' x1 = c(S1,I1,A1,0,0,0)#
#' days= 84
#' inp = list(duration =days, ini = x1)
#' k= 1
#' theta0 = as.numeric(thetas_3travel[[k]][1:6])
#' stochasticmodel_1country(theta0,inp)}
#' @export

stochasticmodel_1country = function(theta,inp){
  status_matrix = matrix(0,nrow = inp$duration,ncol=6)
  status_matrix[1,] = inp$ini

  #Transmission rate, alpha
  harzard1 = function(x,theta){
    h1 = (theta[1] + theta[2])*x[1]*x[2]/sum(x)
    names(h1)=c("hazard1")
    return(h1)
  }

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

  for (i in 2:inp$duration){

    x = status_matrix[(i-1),]

    y1 = rpois(1, harzard1(x,theta))
    #
    y2 =  rpois(1, harzard2(x,theta))
    #
    y3 = rpois(1, harzard3(x,theta))
    #
    y4 =  rpois(1, harzard4(x,theta))
    #
    y5 = rpois(1, harzard5(x,theta))



    ##Susceptible
    if(y1<= x[1]){
      x[1] = x[1] - y1} else{
        y1 = x[1]
        x[1] =0
      }



    #######Infect
    if(y1-y2-y5+x[2]>= 0){
      x[2] = x[2] + y1 - y2 - y5} else{
        y2 =  x[2] + y1 - y5
        x[2] =0
      }

    if(y2 <0){
      y2 = 0
      y5 = x[2] + y1
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

  return(round(status_matrix, digits=0))
}

