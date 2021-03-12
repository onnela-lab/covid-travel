args <- commandArgs(trailingOnly = TRUE)
iteration <- as.integer(args[1])
library(CONETTravel)

#######generating travel data
traveldata_func = function(P, numbercountries, travelrate, durationtravel){
  traveldata = matrix(0, nrow = durationtravel, ncol = numbercountries)
  for( day in 1:durationtravel){
    for (country in 1:numbercountries){
      Totaltravel = P[country]*travelrate
      SdTotaltravel = Totaltravel*.05
      traveldata[day,country] = round(rnorm(1, Totaltravel, SdTotaltravel), digits=0)
    }
  }
  return(traveldata)
}

############### Generate parameters for n countries, R0 from .47 to 6.47

#R0 = (alpha*S/P)/(gamma + beta)

###########generating parameters with R0 in  a range
thetagenerating = function(lowerbound, upperbound){
  tmp2 = 1 # need for kick off
  while(tmp2 >0){
    theta = c( alpha0 = 0,alpha = runif(1,0,1),beta = runif(1,0,.25), delta=runif(1,0,.25),
               eta=1, gamma=runif(1,0,1) )
    tmp1 = theta[2]/(theta[3] + theta[6])
    tmp2 = (tmp1 - lowerbound)*(tmp1 - upperbound)
  }
  return(theta)
}

######### generate theta matrix for all countries

#########countrygroup1
thetafunction1 <- function(numbercountries){
  
  thetamat = matrix(0, nrow=numbercountries, ncol=6)
  for(i in 1:numbercountries){
    thetamat[i,] = thetagenerating(1.1,6.47) # R0 belongs to .47, 6.47
  }
  return(thetamat)
}
#########countrygroup2
thetafunction2 <- function(numbercountries){
  
  thetamat = matrix(0, nrow=numbercountries, ncol=6)
  for(i in 1:numbercountries){
    thetamat[i,] = thetagenerating(1,1.1) # R0 belongs to .47, 6.47
  }
  return(thetamat)
}
#########countrygroup3
thetafunction3 <- function(numbercountries){
  
  thetamat = matrix(0, nrow=numbercountries, ncol=6)
  for(i in 1:numbercountries){
    thetamat[i,] = thetagenerating(0.9,1) # R0 belongs to .47, 6.47
  }
  return(thetamat)
}
#########countrygroup4
thetafunction4 <- function(numbercountries){
  
  thetamat = matrix(0, nrow=numbercountries, ncol=6)
  for(i in 1:numbercountries){
    thetamat[i,] = thetagenerating(0.47,0.9) # R0 belongs to .47, 6.47
  }
  return(thetamat)
}
############
#generating initial outbreak
initialmatrix_func1 =  function(numbercountries){
  initialmatrix = matrix(0, numbercountries, 6)
  for (country in 1:numbercountries){
    P = round(runif(1, 50*10^4, 100*10^6), digits=0)
    I = round(runif(1, 0, 200), digits=0)
    A = round(runif(1, 0, 10), digits=0)
    S = P - I - A
    initialmatrix[country,] = c(S, I, A, 0,0,0)
  }
  
  return(initialmatrix)
}

initialmatrix_func2 =  function(numbercountries){
  initialmatrix = matrix(0, numbercountries, 6)
  for (country in 1:numbercountries){
    P = 5*10^6
    I = round(runif(1, 0, 200), digits=0)
    A = round(runif(1, 0, 10), digits=0)
    S = P - I - A
    initialmatrix[country,] = c(S, I, A, 0,0,0)
  }

  return(initialmatrix)
}





######################
datalist = list()
thetalist = list()
traveldatalist = list()
initialmatrixlist = list()

starttime = proc.time()
######################
#####
numbercountries = 4 # change to any number as we want
theta0 = matrix(0,numbercountries,6)


while(length(thetalist)==0){


  
  theta0[1,] = thetagenerating(1.1,6.47)
  theta0[2,] = thetagenerating(1,1.1)
  theta0[3,] = thetagenerating(0.9,1)
  theta0[4,] = thetagenerating(0.47,0.9)
  
  
  initial_corona = matrix(0,numbercountries,6)

  initial_corona = initialmatrix_func1(numbercountries) 
   
  ###initial_corona[5:8,] = initialmatrix_func2(4) #The last 4 countries with samller Populations




  P = rowSums(initial_corona) #population 
  travelrate = 40/(365*328)
  durationtravel = 84 # 12 weeks
  travelout =  traveldata_func(P, numbercountries, travelrate, durationtravel)
  #######
  ratein = 1 # policy that allows full rate of travel in
  traveloutDivideRegulated = totaltravelout_samerate_regulated(travelout, ratein, P)
  durationquarantine = 0 # quarantine time
  
  #create input to feed in the model
  inp = list(durationtravel = durationtravel, travelregulated = traveloutDivideRegulated,
             initialmatrix = initial_corona, quarantinerate = 1, 
             durationquarantine_adjustedin = rep(durationquarantine,numbercountries))
  ########
  mydata_stochastic =  stochasticmodel_inadjusted_trafficregulated_quarantine(theta0, inp)
  
  #################Stochastic check
  maxinfected_stochastic = rep(0,numbercountries)
  checkcountries_stochastic = rep(0, numbercountries)
  Rt = rep(0,numbercountries)
  #######
  for(i in 1:numbercountries){
    
    #######for country i
    a1 = 1 + (i-1)*6 # susceptible
    a2 = 2 + (i-1)*6 # unobserved infected
    a3 = 3 + (i-1)*6 # active confirmed
    a4 = 4 + (i-1)*6 # recover confirmed cumulative
    a5 = 5 + (i-1)*6 # death cumulative    a6 = 6 + (i-1)*6 #recover unconfirmed cumulative
    d = inp$durationtravel # duration travel
    d1 = d -42  # before the end of travel duration 28 days
    
    cumulativeconfirmed_stochastic = mydata_stochastic[d, a3] + 
      mydata_stochastic[d, a4] + mydata_stochastic[d, a5]
   
    ##Make sure Rt still in range at the end of the period
    Rt[i] = theta0[i,2]/(theta0[i,3]+theta0[i,6])*mydata_stochastic[d1,a1]/P[i]

 
    
    maxinfected_stochastic[i] = max(mydata_stochastic[,a2])
   
    if (mydata_stochastic[d, a1] > 1 &&  max(mydata_stochastic[d1:d, a2]) > 1 && 
        mydata_stochastic[d, a5]< .4*cumulativeconfirmed_stochastic && 
        mydata_stochastic[d, a4] > .6*cumulativeconfirmed_stochastic
      ) {
      
      checkcountries_stochastic[i] = 1
      
    } #end loop for each country, if the country satisfy all conditions then return 1 
    
    ############
  } # end loop for all countries
  
  ##RATIO MAX CASE OF COUNTRY 1 OVER THE OTHER COUNTRIES TO MAKE SURE EXIST ONE COUNTRY HEAVILY INFECTED
  
  r1 = (maxinfected_stochastic[1:4]/maxinfected_stochastic[1])^(-1)
  r1 = min(r1[-1])
 
  ###################Make sure Rt for all groups still in the range
  R1 = (Rt[1] -1.1)*(Rt[1] - 6.47)
  R2 = (Rt[2] -1.1)*(Rt[2] - 1)
  R3 = (Rt[3] -1)*(Rt[3] - 0.9)
  R4 = (Rt[4] -0.9)*(Rt[4] - 0.47)
  R = max(R1,R2,R3,R4)

  ##############
  if (r1 > 5 && sum(checkcountries_stochastic) == numbercountries&& R<0) {
    thetalist = theta0
    datalist =  mydata_stochastic
    traveldatalist = inp$travelregulated
    initialmatrixlist = inp$initialmatrix
  }
  }


mylist = list(theta=thetalist, data = datalist, travel=traveldatalist, initial = initialmatrixlist)

fname = paste('mylist4countriesiteration_',iteration,".Rdata",sep="")


save(mylist, file=fname)

# 
checktime = proc.time() - starttime
print(checktime)

