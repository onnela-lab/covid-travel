
##This code is created on March - 11, 2021
##Last edit: March 11 10h11 p.m.
##Purpose: Prediction for average control ######

args <- commandArgs(trailingOnly = TRUE)
numbertestday  <- as.integer(args[1])
trials <- as.integer(args[2])

##setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/Realdata/EstimationandPrediction")
library(protoABC)
library(ggplot2)
library(dplyr)
library(matrixcalc)
library(tidyr)
library(rlist)
library(xtable)
library(matrixStats)
library(CONETTravel)
library(cowplot)
tmptmp = proc.time()

load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/coviddataJanJune20.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/flightdataJanJune20.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/flightdataJanJune20mixed19.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/Allinitialsforpredictionchangepoint.Rdata")


#######STEP1 SETTING TRAVEL DATA AND STANDARDIZE DATA FROM THE DAY WITH MORE THAN 500 CASES
travelout_datadivided = flightdataJanJune20mixed19 
sce = 4 #travel as average control
inflation = .1 #10% different compare to shutdown
##GENERAL SETTING FOR EXTRACT OUTPUT LATER
percentile_lower = .025
percentile_upper = .975

###I. Standradize covid data and flight data with covid data in terms of time
mycompletedata = coviddataJanJune20 
temp = length(travelout_datadivided) - nrow(mycompletedata[[1]])
temp1 = 1:temp
travelout_datadivided = travelout_datadivided[-temp1] #REMOVED FLIGHT DATA BEFORE 22JAN, THE 1ST GLOBAL COVID DAY REPORTED 

numbercountries = ncol(travelout_datadivided[[1]] )
numberactualcountries = numbercountries-1 ##number actual remove the fake 1 ZZZ


####Replace all 8 columns by the 3 of the corresponding data used, and put back in A,R,D structure

for( i in 1:numbercountries)
{
  tmp = mycompletedata[[i]][,1:3]##Use JHK data, NEED TO CHANGE 1:3 FOR OTHER if want other data
  tmp[,1] = tmp[,1] - tmp[,2] -tmp[,3] #put back in A,R,D structure
  mycompletedata[[i]] = tmp
}



P = mycompletedata$population #population all countries
worldpopulation = 7.8*10^9  #WORLD POPULATION
P[numbercountries] = worldpopulation - sum(P)

########Make travel data and real data during the testing period 

mydatatest = list() # create a list of data for each country from 1 to n = 93 this case
startday = nrow(mycompletedata[[1]]) - 30 #startday of the testing period
endday = startday + numbertestday #endday of the testing period


for(i in 1:numbercountries){
  mydatatest[[i]] =  mycompletedata[[i]][startday:endday,]
}  

travelout_datadivided = travelout_datadivided[startday:endday]


### ACHIEVED. ALL DATA AND TRAVEL DATA NOW ARE TRIM TO THE CONSIDERING COUNTRY




#################STEP2. setting initial for starting the predition period

##make the 93 country have something with very small R0 to make the model run
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
thetazzz = thetagenerating(0.2,0.4)


############Get best thetahat for each country
thetahat = matrix(0,numbercountries,ncol=6)
for (i in 1:numbercountries){
  if(i == numbercountries){
    thetahat[i,] = thetazzz
  }else{
    thetahat[i,] = as.matrix(initials[[i]][2,])
  }
}

###

initialprediction = matrix(0,nrow=numbercountries, ncol = 6)

for (i in 1:numbercountries){
  if(i == numbercountries){
    initialprediction[i,] = c(P[numbercountries],0,0,0,0,0)
  }else{
    initialprediction[i,] = as.matrix(initials[[i]][1,])
  }
}



durationprediction = endday - startday + 1

##########PART 1. Define functions for the outputs: ################

#####################1. Mean realization#########


infectedprediction_mean =  function( theta, inp){
  
  
  #Predictive deterministic realization
  activeconfirmed_imported = rep(0, numbercountries)
  unobservedinfected_imported_prequarantine = rep(0, numbercountries)
  unobservedinfected_imported_postquarantine = rep(0, numbercountries)
  
  unobservedinfected_max = rep(0, numbercountries) 
  activeconfirmed_max = rep(0, numbercountries) 
  
  activeconfirmed_change = matrix(0, numbercountries, 4)
  newcases_change = matrix(0, numbercountries, 4)
  
  Rtchange = matrix(0, numbercountries, 4)
  
  u1 = deterministicmodel_outadjust_pandemictravel(thetahat, inp)
  u = u1$model_output
  
  
  
  
  for (country in 1:numbercountries){
    
    
    
    a1 = 1 + (country -1)*6
    a2 = 2 + (country -1)*6
    a3 = 3 + (country -1)*6
    a4 = 4 + (country -1)*6
    a5 = 5 + (country -1)*6
    a6 = 6 + (country -1)*6
    #Max case
    unobservedinfected_max[country] = max(u[,a2])
    activeconfirmed_max[country] = max(u[,a3])
    
    ##Total unobserved infected
    unobservedinfected_imported_prequarantine[country] = sum(u1$travelarrival_prequarantine[,a2])
    
    unobservedinfected_imported_postquarantine[country] = sum (u1$travelarrival_postquarantine[,a2])
    
    ####Total new imported active confirmed cases
    activeconfirmed_imported[country] =  sum(u1$activeconfirm_imported[,a3])
    
    
    ###Total new confirmed case
    tmp = rowSums(u[,c(a3,a4,a5)])
    A1 = tmp[1]
    A2 = tmp[length(tmp)] 
    Aslope = (A2 - A1)/A1
    activeconfirmed_change[country,]  = c(A1, A2,A2-A1, Aslope)
    ###Total new case
    tmp1 = rowSums(u[,c(a2,a3,a4,a5,a6)])
    N1 = tmp1[1]
    N2 = tmp1[length(tmp1)] 
    Nslope = (N2 - N1)/N1
    newcases_change[country,]  = c(N1, N2,N2-N1, Nslope)
    ###Changed in Rt
    R0 = theta[country,][2]/(theta[country,][3] + theta[country,][6])
    R1 = R0*u[1, a1]/sum(u[1,a1:a6]) 
    R2 = R0*u[nrow(u), a1]/sum(u[nrow(u),a1:a6])
    R0slope = (R2 -R1)/R1 
    Rtchange[country,]  = c(R1, R2, R2-R1, R0slope)
    
  }
  
  
  
  
  
  output = list(unobservedinfected_imported_prequarantine = unobservedinfected_imported_prequarantine,
                unobservedinfected_imported_postquarantine = unobservedinfected_imported_postquarantine,
                max_unobservedinfected =  unobservedinfected_max,
                max_activeconfirmed = activeconfirmed_max,
                activeconfirmed_imported = activeconfirmed_imported,
                activeconfirmed_change = activeconfirmed_change, newcases_change = newcases_change,
                Rtchange = Rtchange)
  
  return(output)
  
} #end infectedprediction_mean definition

####################



####2. Out quarantine, Predict percentile infected of k stochastic realizations, out adjust, a fair
#####criteria especially for high cases of COVID


infectedprediction_percentile =  function(k, percentile_lower, percentile_upper, theta, inp){
  
  ############
  #########
  
  
  
  #Predictive stochastic realization
  
  CI_activeconfirmed_imported = matrix(0, numbercountries, 2)
  CI_unobservedinfected_imported_prequarantine = matrix(0, numbercountries, 2)
  CI_unobservedinfected_imported_postquarantine = matrix(0, numbercountries, 2)
  CI_unobservedinfected_max = matrix(0, numbercountries, 2)
  CI_activeconfirmed_max = matrix(0, numbercountries, 2)
  
  
  CI_activeconfirmed_change = matrix(0, numbercountries, 8)
  CI_newcases_change = matrix(0, numbercountries, 8)
  CI_Rtchange = matrix(0, numbercountries, 8)
  #########
  allquantiles = numbercountries*3
  CI_totalconfirmedcases = matrix(0, allquantiles, inp$durationtravel)
  CI_totaldeathcases = matrix(0, allquantiles, inp$durationtravel)
  #############
  Aimported_list = list()
  Uimported_prequarantine_list = list()
  Uimported_postquarantine_list = list()
  Imax_list = list()
  Amax_list = list()
  A_list = list()
  N_list = list()
  R_list = list()
  Totalconfirmed_list = list()
  Totaldeath_list = list()
  for( country in 1:numbercountries){
    
    Aimported_list[[country]] =  rep(0, k)
    Uimported_prequarantine_list[[country]] =   rep(0, k)
    Uimported_postquarantine_list[[country]] =   rep(0, k)
    Imax_list[[country]] = rep(0, k)
    Amax_list[[country]] = rep(0, k)
    
    R_list[[country]] = matrix(0, k, 4)
    N_list[[country]] = matrix(0, k, 4)
    A_list[[country]] = matrix(0, k, 4)
    
    Totalconfirmed_list[[country]] = matrix(0, k, inp$durationtravel)
    Totaldeath_list[[country]] = matrix(0, k, inp$durationtravel)
  }
  
  
  for(i in 1:k){
    
    
    u1 = stochasticmodel_outadjust_pandemictravel(thetahat, inp)
    u = u1$model_output
    
    for (country in 1:numbercountries){
      
      
      a1 = 1 + (country -1)*6
      a2 = 2 + (country -1)*6
      a3 = 3 + (country -1)*6
      a4 = 4 + (country -1)*6
      a5 = 5 + (country -1)*6
      a6 = 6 + (country -1)*6
      #Max case
      Imax_list[[country]][i] = max(u[,a2])
      Amax_list[[country]][i] = max(u[,a3])
      
      
      ##Total unobserved infected
      Aimported_list[[country]][i] = sum(u1$activeconfirm_imported[,a3])
      
      Uimported_prequarantine_list[[country]][i] = sum(u1$travelarrival_prequarantine[,a2])
      
      Uimported_postquarantine_list[[country]][i] = sum (u1$travelarrival_postquarantine[,a2])
      
      ###Total new confirmed case
      tmp = rowSums(u[,c(a3,a4,a5)])
      ########Extract total confirmed and death for each country, each row is for 1 realization
      Totalconfirmed_list[[country]][i,] = tmp
      Totaldeath_list[[country]][i,] = u[,a5]
      ##########
      A1 = tmp[1]
      A2 = tmp[length(tmp)] 
      Aslope = (A2 - A1)/A1
      A_list[[country]][i,]  = c(A1, A2,A2-A1, Aslope)
      
      ###Total new case
      tmp1 = rowSums(u[,c(a2,a3,a4,a5,a6)])
      N1 = tmp1[1]
      N2 = tmp1[length(tmp1)] 
      Nslope = (N2 - N1)/N1
      N_list[[country]][i,]  = c(N1, N2,N2-N1, Nslope)
      
      ###Changed in Rt
      R0 = theta[country,][2]/(theta[country,][3] + theta[country,][6])
      R1 = R0*u[1, a1]/sum(u[1,a1:a6]) 
      R2 = R0*u[nrow(u), a1]/sum(u[nrow(u),a1:a6])
      R0slope = (R2 -R1)/R1 
      R_list[[country]][i,]  = c(R1, R2,R2-R1, R0slope)
      
    }##end loop for country
  }##end loop for i, realizations from 1 to k
  
  
  
  for(country in 1:numbercountries){
    
    
    f1 = (country-1)*3+1
    f2 = (country-1)*3+2
    f3 = (country-1)*3+3
    CI_totalconfirmedcases[f1,] = quantile(Totalconfirmed_list[[country]], probs = percentile_lower)
    CI_totalconfirmedcases[f2,] = quantile(Totalconfirmed_list[[country]], probs = .5)
    CI_totalconfirmedcases[f3,] = quantile(Totalconfirmed_list[[country]], probs = percentile_upper)
    
    CI_totaldeathcases[f1,] = quantile(Totaldeath_list[[country]], probs = percentile_lower)
    CI_totaldeathcases[f2,] = quantile(Totaldeath_list[[country]], probs = .5)
    CI_totaldeathcases[f3,] = quantile(Totaldeath_list[[country]], probs = percentile_upper)
    
    
    
    
    
    CI_unobservedinfected_max[country,] = c(quantile(Imax_list[[country]], probs = percentile_lower),
                                            quantile(Imax_list[[country]], probs = percentile_upper) )
    
    
    
    
    
    
    
    CI_activeconfirmed_max[country,] = c(quantile(Amax_list[[country]], probs = percentile_lower),
                                         quantile(Amax_list[[country]], probs = percentile_upper) )
    
    
    CI_activeconfirmed_imported[country,] = c(quantile(Aimported_list[[country]], probs = percentile_lower),
                                              quantile(Aimported_list[[country]], probs = percentile_upper) )
    
    
    CI_unobservedinfected_imported_prequarantine[country,] = c(quantile(Uimported_prequarantine_list[[country]], probs = percentile_lower),
                                                               quantile(Uimported_prequarantine_list[[country]], probs = percentile_upper) )
    
    CI_unobservedinfected_imported_postquarantine[country,] = c(quantile(Uimported_postquarantine_list[[country]], probs = percentile_lower),
                                                                quantile(Uimported_postquarantine_list[[country]], probs = percentile_upper) )
    
    
    CI_activeconfirmed_change[country,] = c(colQuantiles(A_list[[country]], probs = percentile_lower),colQuantiles(A_list[[country]], probs = percentile_upper) )
    
    
    CI_newcases_change[country,] = c(colQuantiles(N_list[[country]], probs = percentile_lower),colQuantiles(N_list[[country]], probs = percentile_upper) )
    
    
    CI_Rtchange[country,] = c(colQuantiles(R_list[[country]], probs = percentile_lower),colQuantiles(R_list[[country]], probs = percentile_upper) )
    
    
    
  }
  
  output = list(CI_activeconfirmed_imported = CI_activeconfirmed_imported,
                CI_unobservedinfected_imported_prequarantine= CI_unobservedinfected_imported_prequarantine,
                CI_unobservedinfected_imported_postquarantine = CI_unobservedinfected_imported_postquarantine,
                CI_max_unobservedinfected = CI_unobservedinfected_max,
                CI_max_activeconfirmed = CI_activeconfirmed_max,
                CI_activeconfirmed_change = CI_activeconfirmed_change, CI_newcases_change = CI_newcases_change,
                CI_Rtchange = CI_Rtchange, CI_totaldeathcases = CI_totaldeathcases,
                CI_totalconfirmedcases = CI_totalconfirmedcases )
  
  
  return(output)
  
} #end infectedprediction_percentile function




###inp1 for fully open

inp1 = list(durationtravel = durationprediction, travelregulated =  travelout_datadivided,
            initialmatrix = initialprediction, quarantinerate = 1, 
            durationquarantine_adjustedout =  rep(0, numbercountries))


#############Shutdown
travelout_datadividedshut = inp1$travelregulated


for(i in 1:length(travelout_datadivided)){
  travelout_datadividedshut[[i]] =  as.matrix(travelout_datadivided[[i]])*0
}
inp2 = inp1 #borrow inp1 information
inp2$travelregulated = travelout_datadividedshut #replace by 0 traffic data
infectmean_shut = infectedprediction_mean(thetahat, inp2) 
#####Average control 

##########define an indicator function making p belongs to (0,1)


indicatorfunc = function(x,lower,upper){
  hlow = x -lower
  hup = x - upper
  if(hup>0){
    x=1
  } 
  
  if(hlow<0){
    x=0
  }
  
  return(x)
}
#borrow inp for the initial conditions and some   
inp = inp1
##Design a function for average control based on data divided as full traffic
myfunction1 = function(traveldata,policymatrix){
  numbercountries =  nrow(traveldata[[1]])
  durationtravel = length(traveldata)
  policymatlist = list()
  for (i in 1:durationtravel){
    tmp = matrix(0, nrow = numbercountries, ncol = numbercountries)
    tmp = matrix(policymatrix[, i], nrow = numbercountries, 
                 ncol = numbercountries)
    policymatlist[[i]] = tmp
  }
  
  traveloutDivideRegulated = list()
  for (i in 1:durationtravel) {
    tmp1 = as.matrix(traveldata[[i]])
    tmp2 = as.matrix(policymatlist[[i]])
    traveloutDivideRegulated[[i]] = round(hadamard.prod(tmp1, 
                                                        tmp2), digits = 0)
  }
  return(traveloutDivideRegulated)
}


averagecontrol_travelout_simplifyregulated = function(countryconsider, threshold){
  
  
  
  thetaconsider = thetahat[countryconsider,]
  initialconsider = inp$initialmatrix[countryconsider,]
  psi = 1 + thetaconsider[2]*initialconsider[1]/P[countryconsider] - thetaconsider[3] -thetaconsider[6]
  initialinfected_consider = initialconsider[2]
  #create a sequence options for p from the inequalities
  pseq = rep(0,inp$durationtravel)
  for (i in 1:inp$durationtravel){
    pseq[i] =((threshold/initialinfected_consider)^(1/i))/psi - 1
  }
  
  
  p = min(pseq)
  
  ################
  p = indicatorfunc(p,0,1) # make sure things fall in 0,1
  
  
  infectconsider  = rep(0,inp$durationtravel)
  
  for(i in 1:inp$durationtravel){
    i1 = i -1
    infectconsider[i] = ((1+p)*psi)^i1 * initialinfected_consider
    
  }
  
  infecttravelin_capacity = infectconsider*p/(p+1)
  infecttravelin_capacity_divide = infecttravelin_capacity/(numbercountries - 1)
  
  #######SET UP THE RATE MATRIX FOR AVERAGE CONTROL
  #Be conservative by using .975 percentile Prediction of other countries, 
  # and Full  travel for prediction of infected
  ratein = 1 # policy that allows full rate of travel in
  durationquarantine = rep(0, numbercountries) #no quarantine
  travelout_regulated = inp$travelregulated
  inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
             initialmatrix = initialprediction, quarantinerate = 1, 
             durationquarantine_adjustedout = durationquarantine)
  
  average_fullTravel = deterministicmodel_outadjust_pandemictravel(thetahat, inp)$model_output
  
  travelrate_matrix = matrix(1,nrow = numbercountries^2, ncol = inp$durationtravel)
  ## Decides how much travel the considering country let other countries travel in
  ##Construct rate matrix where row 1 is 1 let 1 in, row2 is 1 let 2 in, row3 1 let 3 in, so on
  for ( country in 1:numbercountries){
    
    a1 = 1+(country-1)*6
    a2 = 2+ (country-1)*6
    infectrate_Predictioncountry = average_fullTravel[,a2]/rowSums(average_fullTravel[,c(a1,a2)])
    infecttravel_country = rep(0, inp$durationtravel)
    
    for(i in 1:inp$durationtravel){
      travelmat = inp$travelregulated[[i]]
      infecttravel_country[i] = travelmat[country,countryconsider]*infectrate_Predictioncountry[i]
    }
    
    ratein = rep(0,inp$durationtravel)
    
    
    for (i in 1:inp$durationtravel){
      
      if(infecttravel_country[i] < infecttravelin_capacity_divide[i]){
        ratein[i] = 1 
      } else if(infecttravel_country[i] == 0){
        ratein[i] = 0} else {
          ratein[i] = infecttravelin_capacity_divide[i]/ infecttravel_country[i]}
      
    }
    
    
    ############ Design cut-off to simplify the approach
    
    if(min(ratein)< 1/3){
      ratein = rep(0, inp$durationtravel) 
      
    } else{
      if((min(ratein) - 1/3)*(min(ratein) - .5)<0) {
        ratein = rep(1/3, inp$durationtravel)
        
      } else {
        if((min(ratein) -.5)*(min(ratein) - 1)<0) {
          ratein = rep(.5, inp$durationtravel)
          
        } else (ratein = rep(1, inp$durationtravel))
        
      }
    }
    
    #plug in the rate matrix
    a3 = country + (countryconsider-1)*numbercountries
    
    
    travelrate_matrix[a3,] = ratein
    
    
    
  } # end loop of letting rate from other countries enter the consider country
  
  
  travelout_regulated = myfunction1(inp$travelregulated, travelrate_matrix) #average control for travelers
 
  return(travelout_regulated)
}


######Make a list for coordinate
travelin_averagecontrol_list = list()
for(countryconsider in 1:numbercountries){
  
  threshold = max(infectmean_shut$max_unobservedinfected[countryconsider]*(1+inflation),1) 
  ### inlation alllow compare to shutdown, use max for case 0 in max
  travelout_regulated1 = averagecontrol_travelout_simplifyregulated(countryconsider, threshold)
  travelin_averagecontrol_list[[countryconsider]] = travelout_regulated1
  print(countryconsider)
}
#########
mylist = list()

for( t1 in 1:inp$durationtravel){
  mymat = matrix(0,numbercountries, numbercountries)
  for(country in 1:numbercountries){
    
    
    tmpmat = travelin_averagecontrol_list[[country]][[t1]]
    mymat[,country] = tmpmat[,country]
    
  }
  mylist[[t1]] = mymat
  
}



inp4 = inp1
inp4$travelregulated = mylist#coordination travel list by all countries, average  control



infectpercentile = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp4)



totaltravelers = rep(0,numbercountries)

for(inbound in 1:numbercountries){
  for(timestep in 1:length(inp4$travelregulated)){
    tmp = inp4$travelregulated[[timestep]]
    totaltravelers[inbound]  = totaltravelers[inbound] + sum(tmp[,inbound])
  }
}




output = list(travel = totaltravelers, infectpercentile = infectpercentile, data = mydatatest)
fname = paste('Differentinboundchangepointtravelsce_',sce,'_testduration_',durationprediction,'_trials_',trials,".Rdata",sep="")
save(output, file=fname)


proc.time() - tmptmp



