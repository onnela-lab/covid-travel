#The code is created on Feb 7 to run on ABC for Coordinate Effectiveness studies
# args <- commandArgs(trailingOnly = TRUE)
# k <- as.integer(args[1])
# prob1 <- as.integer(args[2])




setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/4.CoordinateEffectiveness")
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





tmp = proc.time()
iteration = 234
fname = paste("mylistiteration_",iteration,".Rdata",sep="")

# load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravelSmall/thetas_3travel.RData")
# load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravelSmall/datas_3travel.RData")
# load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravelSmall/travelout_3dat.RData")
# 

load(fname)
names(mylist)
theta0 = mylist$theta
data = mylist$data
travelout_datadivided = mylist$travel
numbercountries = ncol(travelout_datadivided[[1]] )
initialmat = mylist$initial
theta0



policytime = 14 #use 2 weeks for evaluation
inflation = .1 # inflation allow for the proposed method
##################


trials = 10 # number trials of stochastic relization, need to change to 10000

countryconsider = 2

 
prob1 = 10

prob =prob1/1000
particles = 10 # Need to change to 1000







P = rowSums(initialmat) #population 3 countries
########
thetatruth_list = list()

for (i in 1:numbercountries){
  a1 = 1 + (i-1)*6
  b1 = 6 + (i-1)*6
  thetatruth_list[[i]] = theta0[i,]
}


#######Seperate 3 countries choose A, R, D
mydat = list()

for (i in 1:numbercountries){
  a1 = 3 + (i-1)*6
  b1 = 5 + (i-1)*6
  mydat[[i]] = data[,a1:b1]
}

############## Use the first d1 days as training
d1 = 42 #choose 42 days here

mydat_train = list() # create a list of data for each country from 1 to n

for(i in 1:numbercountries){
  mydat_train[[i]] = mydat[[i]][1:d1,]
}

############# Recover total travel out from each country
travelout_data = matrix(0,nrow=84,ncol=numbercountries)
for ( i in 1:84){
  travelout_data[i,] = rowSums(travelout_datadivided[[i]])
}


travelout_traindata =  travelout_data[1:d1,]
###########Prediction period b1:b2
b1 = nrow(travelout_traindata)-1
b2 = b1 + policytime
#####Proposed method default


probcontrol = .9

percentile_lower = .25
percentile_upper = .975
############

##Number travelers from 1 country to another during the traveling period
# ratein = 1
# traveloutDivideRegulated = totaltravelout_samerate_regulated(travelout_traindata, ratein, P)

traveloutDivideRegulated = travelout_datadivided[1:d1]
###########################################################
###I. ESTIMATION PARAMETERS BASED ON TRAIN DATA: traveloutdat
#############################################################
############# Estimate unobserved infected each country

infect =list() # list of number infected estimated from each country from 1 to n
for (i in 1:numbercountries){
  mydat = mydat_train[[i]]
  x = initialmat[i,]
  tmp = infectfunc(mydat,x)
  tmp1 = c(tmp,0)
  infect[[i]] = tmp1
}


###############Reconstruct Total travel in and out at each time point, each column for one country
inmat = matrix(0, nrow = nrow(travelout_traindata), ncol = nrow(travelout_traindata))
outmat = matrix(0, nrow = nrow(travelout_traindata), ncol = nrow(travelout_traindata))

for (i in 1:length(traveloutDivideRegulated)){
  for(j in 1:numbercountries){
    inmat[i,j] = sum(traveloutDivideRegulated[[i]][,j])
    outmat[i,j] = sum(traveloutDivideRegulated[[i]][j,])
  }
  
}
##############

countries = list() # create a list of inputs for country 1 to n
for (i in 1:numbercountries){
  tmp = list(P = P[i], total_in = inmat[,i], total_out = outmat[,i], infect = infect[[i]],
             durationtravel = length(infect[[i]]), data = mydat_train[[i]] , 
             betas = betafunc(mydat_train[[i]]), deltas = deltafunc(mydat_train[[i]]))
  countries[[i]] = tmp
}




##########append population dynamic for each country
for (i in 1:length(countries)){
  tmp = countries[[i]]
  tmp1 = populationdynamicfunc(tmp)
  countries[[i]] = list.append(countries[[i]], populationdynamic = tmp1)
}


#########################################################

travelin_compartments =  list()  # create a travel in compartment list from other countries to country 1, ...,n

for (i in 1:numbercountries){
  tmp2 = list()
  for (j in 1:numbercountries){
    
    tmp = rep(0, length(traveloutDivideRegulated) )#get number travel from j to i sequence
    for (k in 1:length(traveloutDivideRegulated)){
      tmp[k] = traveloutDivideRegulated[[k]][j,i]
    }
    
    tmp1 = traveloutcompfunc(tmp, countries[[j]]) # compartments travel  from j to i during the travelduration
    
    
    tmp2[[j]] = tmp1
    
  }
  travelin_compartments[[i]] = Reduce(`+`, tmp2)
}



#########################################################

travelout_compartments =  list()  # create a travel in compartment list from other countries to country 1, ...,n

for (i in 1:numbercountries){
  
  tmp = rep(0, length(traveloutDivideRegulated) )#get number travel from j to i sequence
  for (k in 1:length(traveloutDivideRegulated)){
    tmp[k] = sum(traveloutDivideRegulated[[k]][i,])
  }
  
  travelout_compartments[[i]] = traveloutcompfunc(tmp, countries[[i]]) # compartments travel  from country i during the travelduration
}




for (i in 1:length(countries)){
  countries[[i]] = list.append(countries[[i]], x_ini = initialmat[i,],
                               travelin_compartments = travelin_compartments[[i]], 
                               travelout_compartments = travelout_compartments[[i]],
                               thetatruth = thetatruth_list[[i]] )
}

############



########Define an ABC_function to estimate each country
ABC_function =  function(country){
  
  ##Uniform prior for parameters
  prior <- function(n){
    data.frame(
      alpha0 = 0,
      alpha = runif(n,0,2),
      beta = runif(n,0,1),
      delta=runif(n,0,1),
      eta=1,
      gamma=runif(n,0,1)
      
    )
  }
  #########Then prior density
  prior_eval <- function(theta){
    prior_value <-  dunif(theta["alpha"], 0, 2)*dunif(theta["beta"], 0, 1)*dunif(theta["delta"], 0, 1)*dunif(theta["gamma"], 0, 1)
    
    return(prior_value)
  }
  ##########
  ############3. Learn alpha based on available data, beta, and gamma##
  
  alpha_function = function(data,theta,country){
    
    beta = theta[3]
    gamma = theta[6]
    Ut = rowSums(data)
    #Difference Ut
    dUt = c(diff(Ut, differences = 1),0)
    ########
    recovermat = matrix(0,nrow = nrow(data),ncol=6)
    
    infected = rep(0,nrow(data))
    
    dRus = rep(0,nrow(data))
    
    for (i in 1:nrow(data)){
      
      infected[i] = dUt[i]/gamma
      dRus[i] = dUt[i]/gamma*beta
      
    }
    
    
    Ruseq1 = cumsum(dRus)
    
    Ruseq = c(0,Ruseq1[-length(Ruseq1)])
    Sseq = country$populationdynamic - Ut - Ruseq - infected
    SIseq = Sseq*infected
    ##Construct dSseq
    Sseq_p = rep(0, nrow(data))
    
    for (i in 1:nrow(data)){
      Sseq_p[i] = Sseq[i] - country$travelin_compartments[i,1] + country$travelout_compartments[i,1]
    }
    ######S+ and Spre sequence for alpha
    Splus_al = Sseq[-c(nrow(data)-1, nrow(data))]
    
    Spre_al = Sseq_p[-c(1,nrow(data))]
    
    dSseq_al = Splus_al - Spre_al
    
    SIseq_al = SIseq[-c(nrow(data)-1, nrow(data))]
    
    pop_al = country$populationdynamic[-c(nrow(data)-1, nrow(data))]
    
    
    #####Remove 0s in the denom
    index1 = which(SIseq_al<=.5)
    ####Remove negative difference
    index2 = which(dSseq_al<=0.5)
    ###########
    index = c(index1,index2)
    index = unique(index)
    index = sort(index)
    
    ######
    if (length(index!=0)){
      dSseq_al = dSseq_al[-index]
      SIseq_al = SIseq_al[-index]
      pop_al = pop_al[-index]
    }
    
    
    alphas_sim = dSseq_al/SIseq_al*pop_al
    
    ########This help to avoid empty array
    if (length(alphas_sim) == 0){alphas_sim = 10^8}
    ###########
    
    
    return(alphas_sim)
    
  }
 
  ###
###########
  standardize_seq = function(k,country){
    prior <- function(n){
      data.frame(
        alpha0 = 0,
        alpha = runif(n,0,2),
        beta = runif(n,0,1),
        delta=runif(n,0,1),
        eta=1,
        gamma=runif(n,0,1)
        #gamma=1.25/15
      )
    }
    para = prior(k)
    sd = rep(0, country$durationtravel) 
    mylist = list()
    Smat = matrix(0,k,country$durationtravel)
    Imat = matrix(0,k,country$durationtravel)
    Amat = matrix(0,k,country$durationtravel)
    Rmat = matrix(0,k,country$durationtravel)
    Dmat = matrix(0,k,country$durationtravel)
    Rumat = matrix(0,k,country$durationtravel)
    for(i in 1:k){
      theta = as.numeric(para[i,])
      u = stochastic_marginalestimate(theta,country)
      mylist[[i]] = u
      Smat[i,] = t(u[,1])
      Imat[i,] = t(u[,2])
      Amat[i,] = t(u[,3])
      Rmat[i,] = t(u[,4])
      Dmat[i,] = t(u[,5])
      Rumat[i,] = t(u[,6])
      
    }
    S_sd = apply(Smat,2,sd)
    I_sd = apply(Imat,2,sd)
    A_sd = apply(Amat,2,sd)
    R_sd = apply(Rmat,2,sd)
    D_sd = apply(Dmat,2,sd)
    Ru_sd = apply(Rumat,2,sd)
    return(rbind(S_sd, I_sd, A_sd, R_sd, D_sd, Ru_sd))
    
  }
  
  sd = standardize_seq(3000,country)[,-1]
  
  
  country = list.append(country, sd=sd)
  
  
  
  sd_function1 = function(vec1, vec2, country){
    sd = country$sd
    output1 <- sum(((vec1[-1]-vec2[-1])/t(sd[3,]))^2)
    return(output1)
    
  }
  sd_function2 = function(vec1, vec2, country){
    sd = country$sd
    output1 <- sum(((vec1[-1]-vec2[-1])/t(sd[4,]))^2)
    return(output1)
    
  }
  sd_function3 = function(vec1, vec2, country){
    sd = country$sd
    output1 <- sum(((vec1[-1]-vec2[-1])/t(sd[5,]))^2)
    return(output1)
    
  }
  
  #####################
  ##1. L4 DISTANCE
  
  
  distance1 = function(theta, country){
    
    sim <- stochastic_marginalestimate(theta, country)
    ########## Learn beta, delta
    sim1 = sim[,3:5]
    
    betas_sim = betafunc(sim1)
    deltas_sim = deltafunc(sim1)
    ###Learn alpha#######
    alphas = alpha_function(country$data, theta, country)
    
    alphas_sim = alpha_function(sim1,theta,country)
    
    #############
    output1 <- abs(median(betas_sim) - median(country$betas))
    output2 <-  abs(median(deltas_sim) - median(country$deltas))
    output3 <- abs(median(alphas_sim) - median(alphas))
    ############
    output4 <- sd_function1(country$data[,1],sim[,3],country)
    output5 <- sd_function2(country$data[,2],sim[,4],country)
    output6 <- sd_function3(country$data[,3],sim[,5],country)
    #############
    
    output <- (output1 + output2 + output3)^.5 + (1/nrow(country$data)*(output4 + output5+ output6))^.5
    
    return(output)
  }
  
  ##
  
  abc_post_1 <- abc_start(
    prior,
    distance1, 
    distance_args = country,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = TRUE)
  )
  
  thetahat = apply(abc_post_1, 2, median)
  
  return(thetahat)
}#end the estimation function

#######Estimate parameters all countries
# thetahat = matrix(0,nrow=numbercountries, ncol = 6)
# for(i in 1:numbercountries){
#   thetahat[i,] = as.numeric(ABC_function(countries[[i]]))
# }

# Replace by truth
thetahat = theta0


###########################################################
###II. DIFFERENT COORDINATE SCENARIOS DEPLOYED AND PREDICTION OF UNOBSERVED CASES
#############################################################



############Recover the average path each country
initialprediction = matrix(0,nrow=numbercountries, ncol = 6)
for(i in 1:numbercountries){
avereal = recoverfunc(countries[[i]],thetahat[i,])
initialprediction[i,] =  as.numeric(avereal[nrow(avereal) -1,])
}

###########

travelprediction = travelout_data[b1:b2,]
durationprediction = nrow(travelprediction)

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
  
}

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
  
  
  Aimported_list = list()
  Uimported_prequarantine_list = list()
  Uimported_postquarantine_list = list()
  Imax_list = list()
  Amax_list = list()
  A_list = list()
  N_list = list()
  R_list = list()
  
  for( country in 1:numbercountries){
    
    Aimported_list[[country]] =  rep(0, k)
    Uimported_prequarantine_list[[country]] =   rep(0, k)
    Uimported_postquarantine_list[[country]] =   rep(0, k)
    Imax_list[[country]] = rep(0, k)
    Amax_list[[country]] = rep(0, k)
    
    R_list[[country]] = matrix(0, k, 4)
    N_list[[country]] = matrix(0, k, 4)
    A_list[[country]] = matrix(0, k, 4)
    
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
      
    }
  }
  
  
  
  for(country in 1:numbercountries){
    
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
                CI_Rtchange = CI_Rtchange)
  
  
  return(output)
  
}



############################
####3 & 4. In quarantine, Predict mean and percentile infected of k stochastic realizations, a one way
####criteria and not fair
###############################
infectedprediction_mean_inquarantine =  function( theta, inp){
  
  
  #Predictive deterministic realization
  activeconfirmed_imported = rep(0, numbercountries)
  unobservedinfected_imported_prequarantine = rep(0, numbercountries)
  unobservedinfected_imported_postquarantine = rep(0, numbercountries)
  
  unobservedinfected_max = rep(0, numbercountries) 
  activeconfirmed_max = rep(0, numbercountries) 
  
  activeconfirmed_change = matrix(0, numbercountries, 4)
  newcases_change = matrix(0, numbercountries, 4)
  
  Rtchange = matrix(0, numbercountries, 4)
  
  u1 = deterministicmodel_inadjust_pandemictravel(thetahat, inp)
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
  
}





###########################

infectedprediction_percentile_inquarantine =  function(k, percentile_lower, percentile_upper, theta, inp){
  
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
  
  
  Aimported_list = list()
  Uimported_prequarantine_list = list()
  Uimported_postquarantine_list = list()
  Imax_list = list()
  Amax_list = list()
  A_list = list()
  N_list = list()
  R_list = list()
  
  for( country in 1:numbercountries){
    
    Aimported_list[[country]] =  rep(0, k)
    Uimported_prequarantine_list[[country]] =   rep(0, k)
    Uimported_postquarantine_list[[country]] =   rep(0, k)
    Imax_list[[country]] = rep(0, k)
    Amax_list[[country]] = rep(0, k)
    
    R_list[[country]] = matrix(0, k, 4)
    N_list[[country]] = matrix(0, k, 4)
    A_list[[country]] = matrix(0, k, 4)
    
  }
  
  
  for(i in 1:k){
    
    
    u1 = stochasticmodel_inadjust_pandemictravel(thetahat, inp)
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
      
    }
  }
  
  
  
  for(country in 1:numbercountries){
    
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
                CI_Rtchange = CI_Rtchange)
  
  
  return(output)
  
}




###########Scene 1. Full traffic###########
ratein = 1 # policy that allows full rate of travel in
durationquarantine = rep(0, numbercountries) #no quarantine

#travelout_regulated = totaltravelout_samerate_regulated(travelprediction, ratein, P)
#use data from available already for full travel
travelout_regulated =  travelout_datadivided[b1:b2]
inp1 = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
           initialmatrix = initialprediction, quarantinerate = 1, 
           durationquarantine_adjustedout = durationquarantine)

#########PREDICTION OF SCENE 1############################

infectmean_scene1 = infectedprediction_mean(thetahat, inp1) 

infectpercentile_scene1 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp1)

totaltravelers_scene1 = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = inp1$travelregulated[[timestep]]
  totaltravelers_scene1 = totaltravelers_scene1 + sum(tmp)
}


###########Scene 2. Shutdown###########

ratein = 0 # shutdown
durationquarantine = rep(0, numbercountries) #no quarantine
travelout_regulated = totaltravelout_samerate_regulated(travelprediction, ratein, P)

inp2 = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
           initialmatrix = initialprediction, quarantinerate = 1, 
           durationquarantine_adjustedout = durationquarantine)


infectmean_scene2 = infectedprediction_mean(thetahat, inp2) 


infectpercentile_scene2 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp2)


totaltravelers_scene2 = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = inp2$travelregulated[[timestep]]
  totaltravelers_scene2 = totaltravelers_scene2 + sum(tmp)
}

################

###########Scene 3. Open but require 14 days quarantine for all enter###########

inp3 = inp1 ##reuse most of the part of in 1
durationquarantine = rep(14, numbercountries) #14 days quarantine

inp3$durationquarantine_adjustedout = durationquarantine

infectmean_scene3 = infectedprediction_mean(thetahat, inp3) 


infectpercentile_scene3 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp3)


totaltravelers_scene3 = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = inp3$travelregulated[[timestep]]
  totaltravelers_scene3 = totaltravelers_scene3 + sum(tmp)
}

totaltravelers_scene3_adjusted = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = inp3$travelregulated[[timestep]]
  totaltravelers_scene3_adjusted = totaltravelers_scene3_adjusted + round(sum(tmp)*.05, digits=0)
}



# ########## Policy 4. Quarantine for countries with number new confirmed cases over than a certain threshold
# ratein = 1 
# dailyaverage_ratethreshold = 20/100000
# 
# for (country in 1:numbercountries){
#   acumulate_confirmcases = rowSums(countries[[country]]$data)
#   daily_confirmcases = diff(acumulate_confirmcases,  differences = 1)
#   dailyrate_confirmcases = daily_confirmcases/P[country]
#   a1 = length(dailyrate_confirmcases) - 14
#   a2 = length(dailyrate_confirmcases)
#   if(mean(dailyrate_confirmcases[a1:a2])>dailyaverage_ratethreshold){
#     durationquarantine[country] = 14 #quarantile for these countries
#   }else{
#     durationquarantine[country] = 0 #no quarantine
#   }
# }
# 
# #travelout_regulated = totaltravelout_samerate_regulated(travelprediction, ratein, P)
# travelout_regulated =  travelout_datadivided[b1:b2]
# inp4 = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
#            initialmatrix = initialprediction, quarantinerate = 1, 
#            durationquarantine_adjustedout = durationquarantine)
# 
# infectmean_scene4 = infectedprediction_mean(thetahat, inp4) 
# 
# infectpercentile_scene4 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp4)
# 
# 
# totaltravelers_scene4 = 0
# for(timestep in 1:length(travelout_regulated)){
#   tmp = travelout_regulated[[timestep]]
#   totaltravelers_scene4 = totaltravelers_scene4 + sum(tmp[,countryconsider])
# }
# ##
# totaltravelers_scene4_adjusted = 0
# 
# for(timestep in 1:length(travelout_regulated)){
#   tmp = travelout_regulated[[timestep]]
#   for (country in 1:numbercountries){
#     if(durationquarantine[country]>0){
#       totaltravelers_scene4_adjusted = totaltravelers_scene4_adjusted + round(tmp[country ,countryconsider]*.05, digits=0)
#       
#     } else {
#       totaltravelers_scene4_adjusted = totaltravelers_scene4_adjusted + tmp[country ,countryconsider]
#     }
#     }
#   }
# 


#################Scene 5. Simplify of average control policy##########
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
    travelout_regulated = totaltravelout_samerate_regulated(travelprediction, ratein, P)
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
    
 
    
    
    
    travelout_regulated = totaltravelout_regulated(travelprediction, travelrate_matrix, P) #total travelers regulated by policy 3
    return(travelout_regulated)
}


######Make a list for coordinate
travelin_averagecontrol_list = list()
for(countryconsider in 1:numbercountries){
  
  threshold = infectmean_scene2$max_unobservedinfected[countryconsider]*(1+inflation) # inlation alllow compare to shutdown
  travelout_regulated1 = averagecontrol_travelout_simplifyregulated(countryconsider, threshold)
  travelin_averagecontrol_list[[countryconsider]] = travelout_regulated1
  
}
#########
mylist = list()
###inp$durationtravel =15, time prediction

for( t1 in 1:inp$durationtravel){
  mymat = matrix(0,numbercountries, numbercountries)
  for(country in 1:numbercountries){
  

  tmpmat = travelin_averagecontrol_list[[country]][[t1]]
  mymat[,country] = tmpmat[,country]
  
  }
  mylist[[t1]] = mymat
  
}

travelout_regulated = mylist#coordination travel list by all countries, average  control

######################
durationquarantine = rep(0, numbercountries) #No quarantine needed

inp4 = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
             initialmatrix = initialprediction, quarantinerate = 1, 
             durationquarantine_adjustedout = durationquarantine)
  
  
  infectpercentile_scene4 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp4)

  infectmean_scene4 = infectedprediction_mean(thetahat, inp4) 
  
  totaltravelers_scene4 = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = inp4$travelregulated[[timestep]]
    totaltravelers_scene4  = totaltravelers_scene4 + sum(tmp)
  }
  

####################Scene 5, countries 1,2,3,4  require all travel in must be 14 days quarantine
###countries 5,6,7,8 let things fully open #############

inp5 = inp3 #reused most of the part inp3, all countries ask for quarantine, travel data as full
  
inp5$durationquarantine_adjustedin = c(rep(14,4),0,0,0,0) 
  
infectmean_scene5 = infectedprediction_mean_inquarantine(thetahat, inp5) 
infectpercentile_scene5 = infectedprediction_percentile_inquarantine(trials,percentile_lower ,percentile_upper, thetahat, inp5)
  
  
  totaltravelers_scene5 = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = inp5$travelregulated[[timestep]]
    totaltravelers_scene5 = totaltravelers_scene5 + sum(tmp)
  }
  ###########5,6,7,8 with full
  ###1,2,3,4 with 5%
  totaltravelers_scene5_adjusted = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = inp5$travelregulated[[timestep]]
    totaltravelers_scene5_adjusted = totaltravelers_scene5_adjusted + round(sum(tmp[,1:4])*.05, digits=0)+round(sum(tmp[,5:8]), digits=0)
  }

############################################################  
  ####################Scene 6, countries 1,2,3,4  require all travel in must be 14 days quarantine
  ###countries 5,6,7,8 shutdown #############
  
  inp6 = inp3 #reused most of the part inp3, all countries ask for quarantine, travel data as full
  
  
  ###Replace travel in of countries 5 to 8 by 0
  traveldata6 = inp6$travelregulated
  for (i in 1:length(traveldata6)){
    tmp= traveldata6[[i]]
    tmp[,5:8] = tmp[,5:8]*0
    traveldata6[[i]] = tmp
  }
  
  inp6$travelregulated = traveldata6
  #######################################

  inp6$durationquarantine_adjustedin = c(rep(14,4),0,0,0,0) 
  
  infectmean_scene6 = infectedprediction_mean_inquarantine(thetahat, inp6) 
  infectpercentile_scene6 = infectedprediction_percentile_inquarantine(trials,percentile_lower ,percentile_upper, thetahat, inp6)
  
  
  totaltravelers_scene6 = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = inp6$travelregulated[[timestep]]
    totaltravelers_scene6 = totaltravelers_scene6 + sum(tmp)
  }
  ###########5,6,7,8 with full
  ###1,2,3,4 with 5%
  totaltravelers_scene6_adjusted = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = inp6$travelregulated[[timestep]]
    totaltravelers_scene6_adjusted = totaltravelers_scene6_adjusted + round(sum(tmp[,1:4])*.05, digits=0)+round(sum(tmp[,5:8]), digits=0)
  }
  
  ############################################################  
  #############Scene 7, Countries 1,2,3,4 average control, countries 5,6,7,8 fully open
  inp7 = inp4 #reused most of the part inp3, all countries ask for quarantine, travel data as full
  
  
  infectmean_scene7 = infectedprediction_mean(thetahat, inp7) 
  infectpercentile_scene7 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp7)
  
  
  totaltravelers_scene7 = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = inp7$travelregulated[[timestep]]
    totaltravelers_scene7 = totaltravelers_scene7 + sum(tmp)
  }
  
  ####################Scene 8, countries 1,2,3,4  require average control
  ###countries 5,6,7,8 shutdown #############
  
  inp8 = inp4 #reused most of the part inp4, all countries ask for quarantine, travel data as full
  
  
  ###Replace travel in of countries 5 to 8 by 0
  traveldata8 = inp8$travelregulated
  for (i in 1:length(traveldata6)){
    tmp= traveldata8[[i]]
    tmp[,5:8] = tmp[,5:8]*0
    traveldata8[[i]] = tmp
  }
  
  inp8$travelregulated = traveldata8
  #######################################
  

  infectmean_scene8 = infectedprediction_mean(thetahat, inp8) 
  infectpercentile_scene8 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp8)
  
  
  totaltravelers_scene8 = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = inp8$travelregulated[[timestep]]
    totaltravelers_scene8 = totaltravelers_scene8 + sum(tmp)
  }
#######################################################################################
#################################3. REPORT########################
######1. TOTAL TRAVELERS
#We create a matrix of 3 rows, and 8 columns, each column is corresponding to one policy
#1st row, number people can travel
#2nd row, percentage compare to full travel
#3rd row, expected number travelers due to the restriction

  
##SCENE 1: OPEN FULL, SCENE 2 SHUTDOWN, SCENE 3: ALL 14 DAYS QUARANTINE, SCENE 4: ALL AVERAGE CONTROL
##SCENE 5: 1-4 14 days quarantine, 5-8 open, SCENE 6: 1-4 14 days quarantine, 5-8 shutdown
##SCENE 7: 1-4 average control, 5-8 open, SCENE 8: 1-4 average control, 5-8 shutdown
  
travelmatrix = matrix(0,3,8)

###First row travelers total
for(i in 1:8){
  
  
  tmp = get( paste("totaltravelers_scene", i, sep="") ) 
  tmp
  
  
  
  if (is.na(tmp) == 0){
    travelmatrix[1,i] =  tmp
    
  } else{
    travelmatrix[1,i] =  NA
  }
}
  
###2nd row, percentage compare to full travel
travelmatrix[2,] = round(travelmatrix[1,]/travelmatrix[1,1], digits = 2)
##### #3rd row, expected number travelers due to the restriction, SCENE 3, 5,6 NEED ADJUSTED DUE TO QUARANTINE
travelmatrix[3,] = travelmatrix[1,]
travelmatrix[3,c(3,5,6)] =  c(totaltravelers_scene3_adjusted, totaltravelers_scene5_adjusted, totaltravelers_scene6_adjusted) #replaced by the expected values
travelmatrix[3,] = round(travelmatrix[3,]/travelmatrix[1,1], digits =2)
totaltravelers = travelmatrix[1,1] # total travel in the country considered under the full load


###################################################
################ CI  COMPARE###################
#We return a matrix of 8 columns, each column corresponding to one policy
#1th and 2nd row, lower and upper percentile percent unobserved infected enter the country prequarantine
#3rd and 4th row, lower and upper percentile percent unobserved infected enter the country post quarantine
#5th and 6tth row, lower and upper percentile percent total active confirmed imported by travel
#7th and 8th row, lower and upper percentile relative change in new active confirmed
#9th and 10th row, lower and upper percentile relative change in new cases change
#11th and 12th row, lower and upper percentile reduction in Rt, magnitude
#13th and 14th row, lower and upper percentile reduction in Rt, relative

for (countryconsider in 1:numbercountries){

CI_policieseffect = matrix(0,14,8)

for(i in 1:8){
  
 
  tmp = get( paste("infectpercentile_scene", i, sep="") ) 
  names(tmp)
  
  
  if(length(tmp) > 1){
    
    CI_policieseffect[1:2, i] =  tmp$CI_unobservedinfected_imported_prequarantine[countryconsider,]/totaltravelers 
   
    CI_policieseffect[3:4, i] =  tmp$CI_unobservedinfected_imported_postquarantine[countryconsider,]/totaltravelers 
    
    
    CI_policieseffect[5:6, i] =  tmp$CI_activeconfirmed_imported[countryconsider,]/totaltravelers 
    
    CI_policieseffect[7:8, i] =  tmp$CI_activeconfirmed_change[countryconsider, c(4,8)]
    
    CI_policieseffect[9:10, i] =  tmp$CI_newcases_change[countryconsider, c(4,8)]
    
    CI_policieseffect[11:12, i] =  tmp$CI_Rtchange[countryconsider, c(3,6)]
    
    CI_policieseffect[13:14, i] =  tmp$CI_Rtchange[countryconsider, c(4,8)]
    
    
    
    
    
  } else{
    CI_policieseffect[,i] =  NA
    
  }
}



############################################

# ## 2. CI  COMPARE, DIFFERENCE
#We return a matrix of 8 columns, each column corresponding to one policy
#1th and 2nd row, lower and upper percentile  unobserved infected enter the country prequarantine
#3rd and 4th row, lower and upper percentile unobserved infected enter the country post quarantine
#5th and 6tth row, lower and upper percentile total active confirmed imported by travel
#7th and 8th row, lower and upper percentile  change in new active confirmed
#9th and 10th row, lower and upper percentile  change in new cases change
#11th and 12th row, lower and upper percentile  in Rt at the beginning
#13th and 14th row, lower and upper percentile reduction in Rt at the end

CI_policieseffect_difference = matrix(0,14,8)

for(i in 1:8){
  
  
  tmp = get( paste("infectpercentile_scene", i, sep="") ) 
  
  
  
  if(length(tmp) > 1){
    
    CI_policieseffect_difference[1:2, i] =  tmp$CI_unobservedinfected_imported_prequarantine[countryconsider,]
    
    CI_policieseffect_difference[3:4, i] =  tmp$CI_unobservedinfected_imported_postquarantine[countryconsider,]
    
    
    CI_policieseffect_difference[5:6, i] =  tmp$CI_activeconfirmed_imported[countryconsider,]
    
    CI_policieseffect_difference[7:8, i] =  tmp$CI_activeconfirmed_change[countryconsider, c(3,7)]
    
    CI_policieseffect_difference[9:10, i] =  tmp$CI_newcases_change[countryconsider, c(3,7)]
    
    CI_policieseffect_difference[11:12, i] =  tmp$CI_Rtchange[countryconsider, c(1,5)]
    
    CI_policieseffect_difference[13:14, i] =  tmp$CI_Rtchange[countryconsider, c(2,6)]
    
  } else{
    CI_policieseffect_difference[,i] =  NA
    
  }
}


#######################################
finaloutput = rbind(travelmatrix,CI_policieseffect,CI_policieseffect_difference)

fname = paste('Coordinateeffectcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")

write.table(finaloutput,file=fname, row.names = F)
}
#########################################################
#########collect back countries by groups of Rt
group1 = c(1,5) # country highly infected
group2 = c(2,6) # country with Rt belong 1 and 1.1
group3 = c(3,7) #country with Rt belong .9 and 1
group4 = c(4,8) #country with Rt less than .9
datagroup1 =  matrix(0,31,numbercountries) 
datagroup2 =  matrix(0,31,numbercountries) 
datagroup3 =  matrix(0,31,numbercountries) 
datagroup4 =  matrix(0,31,numbercountries) 


for (countryconsider in 1:numbercountries){

if(countryconsider %in% c(1,5)){
  fname1 = paste('Coordinateeffectcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
  tmp = read.table(fname1, header =T)
  datagroup1 = datagroup1 + tmp
  datagroup1 = datagroup1/2
  fname2 = paste('Coordinateeffectcountrygroup1_durationpolicy_',policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
  write.table(datagroup1,file=fname2, row.names = F)
}
###########
  if(countryconsider %in% c(2,6)){
    fname1 = paste('Coordinateeffectcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
    tmp = read.table(fname1, header =T)
    datagroup2 = datagroup2 + tmp
    datagroup2 = datagroup2/2
    fname2 = paste('Coordinateeffectcountrygroup2_durationpolicy_',policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
    write.table(datagroup2,file=fname2, row.names = F)
  }
  ###########
  if(countryconsider %in% c(3,7)){
    fname1 = paste('Coordinateeffectcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
    tmp = read.table(fname1, header =T)
    datagroup3 = datagroup3 + tmp
    datagroup3 = datagroup3/2
    fname2 = paste('Coordinateeffectcountrygroup3_durationpolicy_',policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
    write.table(datagroup3,file=fname2, row.names = F)
  }
 ################ 
  ###########
  if(countryconsider %in% c(4,8)){
    fname1 = paste('Coordinateeffectcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
    tmp = read.table(fname1, header =T)
    datagroup4 = datagroup4 + tmp
    datagroup4 = datagroup4/2
    fname2 = paste('Coordinateeffectcountrygroup4_durationpolicy_',policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
    write.table(datagroup4,file=fname2, row.names = F)
  }
  
}



