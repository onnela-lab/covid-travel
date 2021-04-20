#' This function gives stochastic realization for n countries with a given regulated
#' strategy and quarantine duration that each destination country required for
#'  all other countries entered it's authority.
#' @param thetamatrix is a matrix of parameters, parameters of each country is on 1 row
#' @param inp is a list include durationtravel : durationtravel (days),
#'  durationquarantine_adjustedin : number of days people travel in have to quarantine based on each country policy,
#' travelregulated: a list of travel allowed from 1 country to another during the duration,
#' initialmatrix is a matrix of initial compartments of countries, each country is on 1 row, and
#' quarantinerate is the rate people follow quarantine
#' @importFrom stats rpois
#' @importFrom stats rmultinom
#' @return  The stochastic realization of n countries with travel data regulated and quarantine in

#' @examples
#' \dontrun{
#' library(CONETTravel)
#' P1 = 10^7
#' I1 = 250
#' A1 = 130
#' S1 = P1 - I1 - A1
#' x1 = c(S1,I1,A1,0,0,0) # State corresponding S,I,A,R,D,Ru country 1
#' P2 = 3*10^6
#' I2 = 20
#' A2 = 10
#'  S2 = P2 - I2 - A2
#'  x2 = c(S2,I2,A2,0,0,0) # State corresponding S,I,A,R,D,Ru country 2
#'  P3 = 2*10^6
#' I3 = 15
#' A3 = 15
#' S3 = P3 - I3 - A3
#'  x3 = c(S3,I3,A3,0,0,0) # State corresponding S,I,A,R,D,Ru country 3
#'  travelout_data = travelout_3dat
#'  initial_corona = as.matrix(rbind(x1,x2,x3) )#initial conditions of 3 countries
#'  P = c(P1, P2, P3) #population 3 countries
#'  k = 13
#'  theta0 = rbind(thetas_3travel[[k]][1:6],thetas_3travel[[k]][7:12],thetas_3travel[[k]][13:18])
#'  ratein = 1 # policy that allows full rate of travel in
#'  traveloutDivideRegulated = totaltravelout_samerate_regulated(travelout_data, ratein, P)
#'  inp = list(durationtravel = nrow(travelout_data), travelregulated = traveloutDivideRegulated,
#'            initialmatrix = initial_corona, quarantinerate = 1, durationquarantine_adjustedin = c(14,14,14))
#'  stochasticmodel_inadjusted_trafficregulated_quarantine(theta0, inp)

#' }

#' @export


stochasticmodel_inadjusted_trafficregulated_quarantine =  function(thetamatrix, inp){

  ##################Defining harzard functions
  #New infected rate, alphas,c(alpha0,alpha, beta, delta, eta, gamma)
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
  ##############
  numbercountries = nrow(inp$initialmatrix)
  compartments = ncol(inp$initialmatrix)
  status_matrix = matrix(0, nrow = inp$durationtravel, ncol = numbercountries*compartments)

  f_in = matrix(0, nrow =  inp$durationtravel, ncol = numbercountries*compartments)
  f_out = matrix(0, nrow =  inp$durationtravel, ncol = numbercountries*compartments)

  totalduration = inp$durationtravel + max(inp$durationquarantine_adjustedin)
  f_in_donequarantine =  matrix(0, nrow = totalduration ,ncol = numbercountries*compartments) # Create a matrix to contain quarantine once done


  initial = rep(0, numbercountries*compartments)

  for(val3 in 1:numbercountries){
    h1 = 1 + (val3 - 1)*6
    h2 = 6 + (val3 - 1)*6
    initial[h1:h2] = inp$initialmatrix[val3,]

  }

  status_matrix[1,] = initial



  for (i in 2: inp$durationtravel){
    traveloutregulated = as.matrix(inp$travelregulated[[i]])
    totaltravelout = rowSums(traveloutregulated)

    for (j in 1:numbercountries){

      c1 = (j-1)*6 + 1
      c2 = j*6
      x = status_matrix[(i-1),c1:c2]
      ##Updated traveling flow
      #Number out from country j
      out = totaltravelout[j]
      if (x[1]+x[2]  > 0){
        outj = c(round(out*x[1]/(x[1]+x[2] ),digits=0), round(out*x[2]/(x[1]+x[2]),digits=0), 0,0,0,0)
        ###############
        f_out[i,c1:c2] = outj
      }else{

        f_out[i,c1:c2] = c(0, 0,0,0,0,0)
      }

      #travel out from country i with sick and susceptible

      theta = thetamatrix[j,]
      ##Generating Poisson values based on hazard functions
      y1 =  rpois(1, harzard1(x,theta))
      #
      y2 =  rpois(1, harzard2(x,theta))
      #
      y3 = rpois(1, harzard3(x,theta))
      #
      y4 =  rpois(1, harzard4(x,theta))
      #
      y5 = rpois(1, harzard5(x,theta))
      #


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



      status_matrix[i,c1:c2] = x

    }

    ##Construct matrix out of compartments for 3 countries

    f_outmat = matrix(0, nrow = numbercountries, ncol = numbercountries*compartments )

    for (val in 1:numbercountries){
      d1 = (val-1)*6 + 1
      d2 = val*6
      f_outtotal = f_out[i,][d1:d2]
      d3 = d1+1
      #Distribute number of infectious from country val to other countries

      infect_outtotal = f_out[i,][d3]# total infect go out from country i
      probdistribute = rep(0,nrow(traveloutregulated))

      for (val6 in 1:nrow(traveloutregulated)){
        if(sum(traveloutregulated[val,])>0){
          probdistribute[val6] = traveloutregulated[val,val6]/sum(traveloutregulated[val,])
        }else{
          probdistribute[val6] = 0
        }
      }

      #Random assign number infectious from the val-country to other countries
      if(sum(traveloutregulated[val,])>0){
        infect_outdistribute = rmultinom(1, size = infect_outtotal, prob = probdistribute)

      } else {
        infect_outdistribute = rep(0,numbercountries)
      }
      ##########

      for (val1 in 1:numbercountries){
        e1 = (val1 -1)*6 + 1
        e2 = val1*6
        if(sum(traveloutregulated[val,])>0){
          #Average out from val to val1
          tmp = round(f_outtotal*traveloutregulated[val,val1]/sum(traveloutregulated[val,]), digits = 0)

          #Adjust susceptible of average out by using infect_outdistribute
          suseptible = tmp[1] + tmp[2] - infect_outdistribute[val1]
          f_outmat[val,e1:e2] = c(suseptible, infect_outdistribute[val1], tmp[3], tmp[4], tmp[5], tmp[6])

        }else{
          f_outmat[val,e1:e2] = rep(0,6)
        }

      }
    }


    ##Construct matrix in of compartments for n countries
    for (val2 in 1:numbercountries){
      f1 = (val2 -1)*6 + 1
      f2 = val2*6
      f_in[i, f1:f2] = colSums(f_outmat[,f1:f2])

    }

    f_in[i,] = round(f_in[i,], digits=0)
    ##Quarantine f_in
    for( qua in 1:numbercountries){
      a1 = 1 + (qua-1)*6
      a2 = 6 + (qua-1)*6
      quarantineinp = inp$quarantinerate*f_in[i,a1:a2]
      inp1 = list(durationquarantine = inp$durationquarantine_adjustedin[qua], ini = quarantineinp )
      theta1 = thetamatrix[qua,]
      i1 = i + inp$durationquarantine_adjustedin[qua]
      f_in_donequarantine[i1,a1:a2] = stochastic_postquarantine(theta1,inp1)

    }

    #######
    update = status_matrix[i,]  +  f_in_donequarantine[i,] + (1 -inp$quarantinerate)*f_in[i,] - f_out[i,]
    update[update<0.5]=0
    status_matrix[i,] = update

  }
  return(round(status_matrix, digits=0))
}
