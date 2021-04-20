#' This function gives a stochastic realization for n countries with a given regulated
#' strategy and quarantine duration required by the destination countries
#' for all travel out countries.
#' It also keeps track the number of new imported active confirmed cases and traveler status before and after done quarantine.
#' @param thetamatrix is a matrix of parameters, parameters of each country is on 1 row
#' @param inp is a list include durationtravel : durationtravel (days),
#'  durationquarantine_adjustedin : number of days people travel in have to quarantine based on each country policy,
#' travelregulated: a list of travel allowed from 1 country to another during the duration,
#' initialmatrix is a matrix of initial compartments of countries, each country is on 1 row, and
#' quarantinerate is the rate people follow quarantine
#'
#' @return  The a stochastic realization of n countries with travel data regulated, also returns
#' travelers status before and after quarantine, and active confirm update during the quarantine time

#' @examples
#' \dontrun{
#' library(CONETTravel)

#' ######function generating parameters with R0 in a given range
#' thetagenerating = function(lowerbound, upperbound){
#' tmp2 = 1 # need for kick off
#' while(tmp2 >0){
#'    theta = c( alpha0 = 0,alpha = runif(1,0,1),beta = runif(1,0,.25), delta=runif(1,0,.25),
#'            eta=1, gamma=runif(1,0,1) )
#' tmp1 = theta[2]/(theta[3] + theta[6])
#'  tmp2 = (tmp1 - lowerbound)*(tmp1 - upperbound)
#' }
#' return(theta)
#' }
#' ############ function generate theta matrix
#' thetafunction <- function(numbercountries){

#' thetamat = matrix(0, nrow=numbercountries, ncol=6)
#' for(i in 1:numbercountries){
#'   thetamat[i,] = thetagenerating(0.47,6.47) # R0 belongs to .47, 6.47
#' }
#' return(thetamat)
#' }

#' ###########initial matrix function
#' initialmatrix_func =  function(numbercountries){
#'  initialmatrix = matrix(0, numbercountries, 6)
#' for (country in 1:numbercountries){
#'  P = round(runif(1, 50000, 20000000000), digits=0)
#' I = round(runif(1, 0, 2000), digits=0)
#'  S = P - I
#'   initialmatrix[country,] = c(S, I, 0, 0,0,0)
#' }
#'  return(initialmatrix)
#'}
#'############function generate travel data
#'traveldata_func = function(P, numbercountries, travelrate, durationtravel){
#'  traveldata = matrix(0, nrow = durationtravel, ncol = numbercountries)
#' for( day in 1:durationtravel){
#'  for (country in 1:numbercountries){
#'    Totaltravel = P[country]*travelrate
#'      SdTotaltravel = Totaltravel*.05
#'      traveldata[day,country] = round(rnorm(1, Totaltravel, SdTotaltravel), digits=0)
#'    }
#'  }
#'  return(traveldata)
#' }
#' #############
#' numbercountries = 3 # choose the number of countries
#' initial_corona =  initialmatrix_func(numbercountries)
#' P = rowSums(initial_corona)
#' travelrate = 40/(365*328) #a given travel rate each day
#' durationtravel = 84 # number of days travel
#generate total travel out data for each country

#' travelout_data = traveldata_func(P, numbercountries, travelrate, durationtravel)
#' #generate theta matrix for each countries

#' thetamatrix = thetafunction(numbercountries)

#' ratein = 1 # policy that allows full rate of travel in
#' traveloutDivideRegulated = totaltravelout_samerate_regulated(travelout_data, ratein, P)
#' inp = list(durationtravel = durationtravel, travelregulated = traveloutDivideRegulated,
#'           initialmatrix = initial_corona, quarantinerate = 1, durationquarantine_adjustedin = rep(14,numbercountries))

#' stochasticmodel_inadjust_pandemictravel(thetamatrix, inp)


#' }

#' @export




stochasticmodel_inadjust_pandemictravel = function (thetamatrix, inp)
{

  harzard1 = function(x, theta) {
    h1 = (theta[1] + theta[2]) * x[1] * x[2]/sum(x)
    names(h1) = c("hazard1")
    return(h1)
  }
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

  numbercountries = nrow(inp$initialmatrix)
  compartments = ncol(inp$initialmatrix)

  status_matrix = matrix(0, nrow = inp$durationtravel, ncol = numbercountries *
                           compartments)
  f_in = matrix(0, nrow = inp$durationtravel, ncol = numbercountries *
                  compartments)
  f_out = matrix(0, nrow = inp$durationtravel, ncol = numbercountries *
                   compartments)


  totalduration = inp$durationtravel + max(inp$durationquarantine_adjustedin) +
    inp$durationtravel - 1




  f_in_donequarantine = matrix(0, nrow = totalduration, ncol = numbercountries *
                                 compartments)  ##Note 2. this is a new addition
  f_in_prequarantine = matrix(0, nrow = totalduration, ncol = numbercountries *
                                compartments)  ##Note 2. this is a new addition

  f_in_activequarantine = matrix(0, nrow = totalduration, ncol = numbercountries *
                                   compartments)  ##Note 2. this is a new addition

  f_in_activenoquarantine = matrix(0, nrow = totalduration, ncol = numbercountries *
                                     compartments)  ##Note 2. this is a new addition

  ##This need since if quarantine required, then during this time active confirmed also updated till the end.
  f_in_activequarantine_add= matrix(0,  nrow = totalduration, ncol = numbercountries * compartments)

  ##This also needed since if no quarantine require then eventually things also become active confirmed
  f_in_activenoquarantine_add= matrix(0, nrow = totalduration, ncol = numbercountries * compartments)


  initial = rep(0, numbercountries * compartments)


  for (val3 in 1:numbercountries) {
    h1 = 1 + (val3 - 1) * 6
    h2 = 6 + (val3 - 1) * 6
    initial[h1:h2] = inp$initialmatrix[val3, ]
  }


  status_matrix[1, ] = initial



  ###########################################
  ##start loop for i at each time step
  for (i in 2:inp$durationtravel) {






    traveloutregulated = as.matrix(inp$travelregulated[[i]])
    totaltravelout = rowSums(traveloutregulated)
    ##internal pandemic at each country at the i-th time step
    for (j in 1:numbercountries) {
      c1 = (j - 1) * 6 + 1
      c2 = j * 6
      x = status_matrix[(i - 1), c1:c2]
      out = totaltravelout[j]
      if (x[1] + x[2] > 0) {
        outj = c(round(out * x[1]/(x[1] + x[2]), digits = 0),
                 round(out * x[2]/(x[1] + x[2]), digits = 0),
                 0, 0, 0, 0)
        f_out[i, c1:c2] = outj
      }
      else {
        f_out[i, c1:c2] = c(0, 0, 0, 0, 0, 0)
      }
      theta = thetamatrix[j, ]
      y1 = rpois(1, harzard1(x, theta))
      y2 = rpois(1, harzard2(x, theta))
      y3 = rpois(1, harzard3(x, theta))
      y4 = rpois(1, harzard4(x, theta))
      y5 = rpois(1, harzard5(x, theta))
      if (y1 <= x[1]) {
        x[1] = x[1] - y1
      }

      else {
        y1 = x[1]
        x[1] = 0
      }
      if (y1 - y2 - y5 + x[2] >= 0) {
        x[2] = x[2] + y1 - y2 - y5
      }
      else {
        y2 = x[2] + y1 - y5
        x[2] = 0
      }
      if (y2 < 0) {
        y2 = 0
        y5 = x[2] + y1
      }
      if (y2 - y3 - y4 + x[3] >= 0) {
        x[3] = x[3] + y2 - y3 - y4
      }
      else {
        y3 = y2 - y4 + x[3]
        x[3] = 0
      }
      if (y3 < 0) {
        y3 = 0
        y4 = x[3] + y2
      }
      x[4] = x[4] + y3
      x[5] = x[5] + y4
      x[6] = x[6] + y5
      status_matrix[i, c1:c2] = x

    }  #end loop for j, internal pandemic done
    ##start update external pandemic by travelling

    ##distribute compartments of travel from one country to another
    f_outmat = matrix(0, nrow = numbercountries, ncol = numbercountries *
                        compartments)
    for (val in 1:numbercountries) {



      d1 = (val - 1) * 6 + 1
      d2 = val * 6
      f_outtotal = f_out[i, ][d1:d2]
      d3 = d1 + 1
      infect_outtotal = f_out[i, ][d3]
      probdistribute = rep(0, numbercountries)
      ####define the weight for multinomial based on traffic data
      for (val6 in 1:numbercountries) {
        if (sum(traveloutregulated[val, ]) > 0) {
          probdistribute[val6] = traveloutregulated[val,
                                                    val6]/sum(traveloutregulated[val, ])
        }else{
          probdistribute[val6] = 0
        }
      }  #end loop for val6

      if (sum(traveloutregulated[val, ]) > 0) {
        infect_outdistribute = rmultinom(1, size = infect_outtotal,
                                         prob = probdistribute)
      }else{
        infect_outdistribute = rep(0, numbercountries)
      }
      ###end multinomial distribute

      ##################
      for (val1 in 1:numbercountries) {
        e1 = (val1 - 1) * 6 + 1
        e2 = val1 * 6
        if (sum(traveloutregulated[val, ]) > 0) {
          tmp = round(f_outtotal * traveloutregulated[val,
                                                      val1]/sum(traveloutregulated[val, ]), digits = 0)
          suseptible = tmp[1] + tmp[2] - infect_outdistribute[val1]
          f_outmat[val, e1:e2] = c(suseptible, infect_outdistribute[val1],
                                   tmp[3], tmp[4], tmp[5], tmp[6])
        }
        else {
          f_outmat[val, e1:e2] = rep(0, 6)
        }
      } # end loop for val1
      #############

    }  #end loop for val, done with distribute compartments of travel from one country to another
    ########################################
    #####################3. Note3: Add up number enter each country

    for (val2 in 1:numbercountries) {
      f1 = (val2 - 1) * 6 + 1
      f2 = val2 * 6
      f_in[i, f1:f2] = colSums(f_outmat[, f1:f2])
    }  ##end val2 group, count total enter the country at step i-th

    f_in[i, ] = round(f_in[i, ], digits = 0)
    #####################################################
    #################Note 4: Add up things after arrivals done with quarantine both active confirmed and added compartments
    for (qua in 1:numbercountries) {
      a1 = 1 + (qua - 1) * 6
      a2 = 6 + (qua - 1) * 6
      a3 = 3 + (qua - 1) * 6
      quarantineinp = inp$quarantinerate * f_in[i, a1:a2]
      ##input with quarantine duration required
      inp1 = list(durationquarantine = inp$durationquarantine_adjustedin[qua],
                  ini = quarantineinp)

      ##inp with no quarantine require or the prequarantine list
      inp2 = list(durationquarantine = 0, ini = quarantineinp)

      ##parameter of the destination country
      theta1 = thetamatrix[qua, ]
      ############
      i1 = i + inp$durationquarantine_adjustedin[qua]
      i2 = i + 0
      ########inp3 to keep track imported active confirmed eventually no matter quarantine or not
      iadd = (inp$durationtravel - 1) + inp$durationquarantine_adjustedin[qua]
      i3 = i + iadd  # this helps to keep track the quarantine imported which eventually become active confirmed even quarantine required
      inp3 = list(durationquarantine = iadd, ini = quarantineinp)
      ###
      tmp = stochastic_postquarantine_separate(theta1, inp1)  #quarantine require
      #############################

      tmp2 = stochastic_postquarantine_separate(theta1,  inp2) #no quarantine require

      tmp3 = stochastic_postquarantine_separate(theta1, inp3) # use to check A evolving eventually
      #####################
      f_in_prequarantine[i2, a1:a2] = tmp2$donequarantine
      f_in_activenoquarantine[i:i3, a3] = tmp3$activeconfirm_eachday[, 3]

      ##Batch back in case no quarantine require
      if (inp$durationquarantine_adjustedin[qua] > 0) {
        f_in_donequarantine[i1, a1:a2] = tmp$donequarantine

        f_in_activequarantine[i:i1, a3] = tmp$activeconfirm_eachday[, 3]
      }
      else {
        f_in_donequarantine[i, a1:a2] = tmp$donequarantine
        f_in_activequarantine[i, a3] = tmp$activeconfirm_eachday[3]
      }

    }



    ############################
    ##########################


    f_in_activequarantine_add = f_in_activequarantine_add + f_in_activequarantine

    f_in_activenoquarantine_add = f_in_activenoquarantine_add + f_in_activenoquarantine

    ############################
    ###########################update with quarantine in after done with quarantine and rate of quarantine
    update = status_matrix[i, ] + f_in_donequarantine[i,
    ] + (1 - inp$quarantinerate) * f_in[i, ] - f_out[i,
    ] + f_in_activequarantine_add[i, ]

    ##############################
    update[update < 0.5] = 0
    status_matrix[i, ] = update













  } # end loop for time step i, i =2 ..84


  mylist = list(model_output = round(status_matrix, digits = 0), activeconfirm_imported = round(f_in_activenoquarantine_add[1:inp$durationtravel,
  ], digits = 0), travelarrival_postquarantine = round(f_in_donequarantine[1:inp$durationtravel,
  ], digits = 0), travelarrival_prequarantine = round(f_in_prequarantine[1:inp$durationtravel,
  ], digits = 0))

  return (mylist)}







