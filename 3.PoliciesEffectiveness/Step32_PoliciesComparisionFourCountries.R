args <- commandArgs(trailingOnly = TRUE)
iteration <- as.integer(args[1])
policytime <- as.integer(args[2])
inflation <- as.integer(args[3])
 



inflation = inflation/100 # inflation allow for the proposed method






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




fname = paste("mylist4countriesiteration_",iteration,".Rdata",sep="")
load(fname)
theta0 = mylist$theta
data = mylist$data
travelout_datadivided = mylist$travel
numbercountries = ncol(travelout_datadivided[[1]] )
initialmat = mylist$initial

trials = 1000 # number trials of stochastic relization, need to change to 1000

for(countryconsider in 1:numbercountries )

{
  
theta0 = mylist$theta
data = mylist$data
travelout_datadivided = mylist$travel  
 
prob1 = 10

prob =prob1/1000
particles = 1000 # Need to change to 1000







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
d = nrow(mydat[[1]])
d1 = d - 42 #choose back 42 days here

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

percentile_lower = .025
percentile_upper = .975
############

##Number travelers from 1 country to another during the traveling period
ratein = 1


traveloutDivideRegulated = totaltravelout_samerate_regulated(travelout_traindata, ratein, P)


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
}

#Estimate parameters all countries
 thetahat = matrix(0,nrow=numbercountries, ncol = 6)
 for(i in 1:numbercountries){
   thetahat[i,] = as.numeric(ABC_function(countries[[i]]))
 }



###########################################################
###II. DIFFERENT POLICIES DEPLOYED AND PREDICTION OF UNOBSERVED CASES
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

##########PART 1. POLICIES: Set traffic to m % (0%: Shutdown, 100%: Fully open) and 14 days quarantine
##Return: n percentile of k stochastic realizations#

#####################Mean of infected#########


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



#Predict percentile infected of k stochastic realizations 


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


infectedprediction_percentile_inquarantine =  function(k, percentile_lower, percentile_upper, theta, inp){

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

###########################
############################


infectedprediction_percentile_outtheninquarantine =  function(k, percentile_lower, percentile_upper, theta, inp){

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


    u1 = stochasticmodel_outtheninadjust_pandemictravel(thetahat, inp)
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




###########Policy 1. Full traffic###########
ratein = 1 # policy that allows full rate of travel in
durationquarantine = rep(0, numbercountries) #no quarantine
travelout_regulated = totaltravelout_samerate_regulated(travelprediction, ratein, P)
inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
           initialmatrix = initialprediction, quarantinerate = 1, 
           durationquarantine_adjustedout = durationquarantine)

#########PREDICTION OF POLICY 1############################

infectmean_policy1 = infectedprediction_mean(thetahat, inp) 

infectpercentile_policy1 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp)

totaltravelers_policy1 = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = travelout_regulated[[timestep]]
  totaltravelers_policy1 = totaltravelers_policy1 + sum(tmp[,countryconsider])
}


###########Policy 2. Shutdown###########

ratein = 0 # shutdown
durationquarantine = rep(0, numbercountries) #no quarantine
travelout_regulated = totaltravelout_samerate_regulated(travelprediction, ratein, P)
inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
           initialmatrix = initialprediction, quarantinerate = 1, 
           durationquarantine_adjustedout = durationquarantine)


infectmean_policy2 = infectedprediction_mean(thetahat, inp) 


infectpercentile_policy2 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp)


totaltravelers_policy2 = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = travelout_regulated[[timestep]]
  totaltravelers_policy2 = totaltravelers_policy2 + sum(tmp[,countryconsider])
}

################

###########Policy 3. Open but require 14 days quarantine for all enter###########

ratein = 1 # Fully open
durationquarantine = rep(0, numbercountries) #0 days quarantine for all
durationquarantine[countryconsider] = 14 #replace the position at the considering country by quarantine 14 days 
travelout_regulated = totaltravelout_samerate_regulated(travelprediction, ratein, P)
inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
           initialmatrix = initialprediction, quarantinerate = 1, 
           durationquarantine_adjustedout = durationquarantine, durationquarantine_adjustedin = durationquarantine)


infectpercentile_policy3 = infectedprediction_percentile_inquarantine(trials,percentile_lower ,percentile_upper, thetahat, inp)


totaltravelers_policy3 = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = travelout_regulated[[timestep]]
  totaltravelers_policy3 = totaltravelers_policy3 + sum(tmp[,countryconsider])
}

totaltravelers_policy3_adjusted = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = travelout_regulated[[timestep]]
  totaltravelers_policy3_adjusted = totaltravelers_policy3_adjusted + round(sum(tmp[,countryconsider])*.05, digits=0)
}



########## Policy 4. Quarantine for countries with number new confirmed cases over than a certain threshold
ratein = 1 
dailyaverage_ratethreshold = 20/100000

for (country in 1:numbercountries){
  acumulate_confirmcases = rowSums(countries[[country]]$data)
  daily_confirmcases = diff(acumulate_confirmcases,  differences = 1)
  dailyrate_confirmcases = daily_confirmcases/P[country]
  a1 = length(dailyrate_confirmcases) - 14
  a2 = length(dailyrate_confirmcases)
  if(mean(dailyrate_confirmcases[a1:a2])>dailyaverage_ratethreshold){
    durationquarantine[country] = 14 #quarantile for these countries
  }else{
    durationquarantine[country] = 0 #no quarantine
  }
}

travelout_regulated = totaltravelout_samerate_regulated(travelprediction, ratein, P)
countryinrequire = rep(0,numbercountries)
countryinrequire[countryconsider]=1

inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
           initialmatrix = initialprediction, quarantinerate = 1, 
           durationquarantine_adjustedout = durationquarantine,
           durationquarantine_adjustedout1 = rep(0,numbercountries),
           countryinrequire=countryinrequire)

 

infectpercentile_policy4 = infectedprediction_percentile_outtheninquarantine(trials,percentile_lower ,percentile_upper, thetahat, inp)


totaltravelers_policy4 = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = travelout_regulated[[timestep]]
  totaltravelers_policy4 = totaltravelers_policy4 + sum(tmp[,countryconsider])
}
##
totaltravelers_policy4_adjusted = 0

for(timestep in 1:length(travelout_regulated)){
  tmp = travelout_regulated[[timestep]]
  for (country in 1:numbercountries){
    if(durationquarantine[country]>0){
      totaltravelers_policy4_adjusted = totaltravelers_policy4_adjusted + round(tmp[country ,countryconsider]*.05, digits=0)
      
    } else {
      totaltravelers_policy4_adjusted = totaltravelers_policy4_adjusted + tmp[country ,countryconsider]
    }
    }
  }


#########################################
#######################################

#######Define this function to make things between 0 and 1
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


##Same threshold for Pol5 to Pol8
threshold = max(infectmean_policy2$max_unobservedinfected[countryconsider]*(1+inflation),1) # inlation alllow compare to shutdown
 #make sure threshold at least 1

#################4. Average control policy##########
#1. Exact control
##Define the traffic regulate function

averagecontrol_travelout_exactregulated = function(countryconsider, threshold){
  
  
  
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
##############
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
##Construct rate matrix where row 1 is 1 let 1 in, row2 is 1 let 2 in, row3 1 let 3 in, so on
for ( country in 1:numbercountries){
  
  a1 = 1+(country-1)*6
  a2 = 2+ (country-1)*6
  infectedPrediction_country = infectpercentile_policy1$unobservedinfected[country,]
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
  
  
  #plug in the rate matrix
  a3 = country + (countryconsider-1)*numbercountries
  
    
  travelrate_matrix[a3,] = ratein
  
} # end loop for country

travelout_regulated = totaltravelout_regulated(travelprediction, travelrate_matrix, P) #total travelers regulated by policy 3

return(travelout_regulated)
} 
######



travelout_regulated = averagecontrol_travelout_exactregulated(countryconsider, threshold)



  durationquarantine = rep(0, numbercountries) #No quarantine needed
  
  inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
             initialmatrix = initialprediction, quarantinerate = 1, 
             durationquarantine_adjustedout = durationquarantine)
  
  
  infectpercentile_policy5 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp)

  infectmean_policy5 = infectedprediction_mean(thetahat, inp) 
  
  totaltravelers_policy5 = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = travelout_regulated[[timestep]]
    totaltravelers_policy5 = totaltravelers_policy5 + sum(tmp[,countryconsider])
  }
  



############################

#b. Simplify control
##Define the traffic regulate function

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
    ##Construct rate matrix where row 1 is 1 let 1 in, row2 is 1 let 2 in, row3 1 let 3 in, so on
    for ( country in 1:numbercountries){
      a1 = 1+(country-1)*6
      a2 = 2+ (country-1)*6
      infectedPrediction_country = infectpercentile_policy1$unobservedinfected[country,]
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
      
    }
    
    travelout_regulated = totaltravelout_regulated(travelprediction, travelrate_matrix, P) #total travelers regulated by policy 3
    return(travelout_regulated)
  } 
######



travelout_regulated = averagecontrol_travelout_simplifyregulated(countryconsider, threshold)



  durationquarantine = rep(0, numbercountries) #No quarantine needed
  
  inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
             initialmatrix = initialprediction, quarantinerate = 1, 
             durationquarantine_adjustedout = durationquarantine)
  
  
  infectpercentile_policy6 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp)

  infectmean_policy6 = infectedprediction_mean(thetahat, inp) 
  
  totaltravelers_policy6 = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = travelout_regulated[[timestep]]
    totaltravelers_policy6 = totaltravelers_policy6 + sum(tmp[,countryconsider])}
  



############################ Probability control



probabilitycontrol_travelout_exactregulated = function(countryconsider, threshold, probcontrol){
  
  thetaconsider = thetahat[countryconsider,]
  initialconsider = inp$initialmatrix[countryconsider,]
  psi = 1 + thetaconsider[2]*initialconsider[1]/P[countryconsider] - thetaconsider[3] -thetaconsider[6]
  initialinfected_consider = initialconsider[2]
  ############
  m1 =  initialinfected_consider/(1-probcontrol) #this number control the I* sequence
  #create a sequence options for p from the inequalities
  pseq = rep(0,inp$durationtravel)
  for (i in 1:inp$durationtravel){
    pseq[i] =((threshold/m1)^(1/i))/psi - 1
  }
  #######
  
  
  p = min(pseq)
  ###############
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
    ##Construct rate matrix where row 1 is 1 let 1 in, row2 is 1 let 2 in, row3 1 let 3 in, so on
    for ( country in 1:numbercountries){
      a1 = 1+(country-1)*6
      a2 = 2+ (country-1)*6
      infectedPrediction_country = infectpercentile_policy1$unobservedinfected[country,]
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
   
     
      

      #plug in the rate matrix
      a3 = country + (countryconsider-1)*numbercountries
      
      
      travelrate_matrix[a3,] = ratein
      
    }
    
    travelout_regulated = totaltravelout_regulated(travelprediction, travelrate_matrix, P) #total travelers regulated by policy 3
    return(travelout_regulated)
  }
######


travelout_regulated = probabilitycontrol_travelout_exactregulated(countryconsider, threshold, probcontrol)


durationquarantine = rep(0, numbercountries) #No quarantine needed

inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
           initialmatrix = initialprediction, quarantinerate = 1, 
           durationquarantine_adjustedout = durationquarantine)


infectpercentile_policy7 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp)

infectmean_policy7 = infectedprediction_mean(thetahat, inp) 

totaltravelers_policy7 = 0
for(timestep in 1:length(travelout_regulated)){
  tmp = travelout_regulated[[timestep]]
  totaltravelers_policy7  = totaltravelers_policy7  + sum(tmp[,countryconsider])
}





####################Simplify version of probability

probabilitycontrol_travelout_simplifyregulated = function(countryconsider, threshold, probcontrol){
  thetaconsider = thetahat[countryconsider,]
  initialconsider = inp$initialmatrix[countryconsider,]
  psi = 1 + thetaconsider[2]*initialconsider[1]/P[countryconsider] - thetaconsider[3] -thetaconsider[6]
  initialinfected_consider = initialconsider[2]
  ############
  m1 =  initialinfected_consider/(1-probcontrol) #this number control the I* sequence
  #create a sequence options for p from the inequalities
  pseq = rep(0,inp$durationtravel)
  for (i in 1:inp$durationtravel){
    pseq[i] =((threshold/m1)^(1/i))/psi - 1
  }
  #######
  
  
  p = min(pseq)
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
    ##Construct rate matrix where row 1 is 1 let 1 in, row2 is 1 let 2 in, row3 1 let 3 in, so on
    for ( country in 1:numbercountries){
      a1 = 1+(country-1)*6
      a2 = 2+ (country-1)*6
      infectedPrediction_country = infectpercentile_policy1$unobservedinfected[country,]
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
      
    }
    
    travelout_regulated = totaltravelout_regulated(travelprediction, travelrate_matrix, P) #total travelers regulated by policy 3
    return(travelout_regulated)
  }
######



travelout_regulated = probabilitycontrol_travelout_simplifyregulated(countryconsider, threshold, probcontrol)


  durationquarantine = rep(0, numbercountries) #No quarantine needed
  
  inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
             initialmatrix = initialprediction, quarantinerate = 1, 
             durationquarantine_adjustedout = durationquarantine)
  
  
  infectpercentile_policy8 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp)

  infectmean_policy8 = infectedprediction_mean(thetahat, inp) 
  
  totaltravelers_policy8 = 0
  for(timestep in 1:length(travelout_regulated)){
    tmp = travelout_regulated[[timestep]]
    totaltravelers_policy8  = totaltravelers_policy8  + sum(tmp[,countryconsider])
  }
  
  



#################################3. REPORT########################



policies = 8 #

travelinbound = rep(0,policies)

for(i in 1:policies){
  tmp = get( paste("totaltravelers_policy", i, sep="") )
  travelinbound[i] =  tmp
     }
################ The adjusted inbound travel in case quarantine needed

adjustedindex = c(3,4)

travelinbound1 = travelinbound

for(i in 1:length(adjustedindex)){
  i1 = adjustedindex[i]
  tmp = get( paste("totaltravelers_policy", i1,"_adjusted", sep="") )
  travelinbound1[i1] = tmp
}




 
#############POLICIES  COMPARE###########

#1st and 2nd row, lower and upper percentile number unobserved infected enter the country prequarantine
#3rd and 4th row, lower and upper percentile percent unobserved infected enter the country prequarantine

#Rows 5-6: lower and upper number unobserved infected enter the country post quarantine
#7 th and 8th row, lower and upper percentile percent unobserved infected enter the country post quarantine

#9th and 10th row, lower and upper percentile total active confirmed imported
#Row 11 and 12 percent of active confirmed imported people


#Row 13 and 14, lower and upper percentile percent change in new cases
#Row 15th and 16th, lower and upper percentile relative change in new cases

#Row 17  and 18, lower and upper percentile percent change in active confirmed
#Row 19th and 20th, lower and upper percentile relative change in active confirmed


#21st and 22nd row, lower and upper percentile reduction in Rt
#23th and 24th row, lower and upper percentile relative reduction in Rt




CI_policieseffect = matrix(0,24,policies)

for(i in 1:policies){
  tmp = get( paste("infectpercentile_policy", i, sep="") )
  inboundtravel = travelinbound[i]

    CI_policieseffect[1:2,i] = tmp$CI_unobservedinfected_imported_prequarantine[countryconsider,]
 if(inboundtravel >0){
    CI_policieseffect[3:4,i] = tmp$CI_unobservedinfected_imported_prequarantine[countryconsider,]/inboundtravel
    } else {
    CI_policieseffect[3:4,i] = tmp$CI_unobservedinfected_imported_prequarantine[countryconsider,]*0
    }


   CI_policieseffect[5:6,i] = tmp$CI_unobservedinfected_imported_postquarantine[countryconsider,]

    if(inboundtravel >0){
    CI_policieseffect[7:8,i] = tmp$CI_unobservedinfected_imported_postquarantine[countryconsider,]/inboundtravel
    } else {
    CI_policieseffect[7:8,i] = tmp$CI_unobservedinfected_imported_postquarantine[countryconsider,]*0
    }



   CI_policieseffect[9:10,i] = tmp$CI_activeconfirmed_imported[countryconsider,]

    if(inboundtravel >0){
      CI_policieseffect[11:12,i] = tmp$CI_activeconfirmed_imported[countryconsider,]/inboundtravel
    } else {
      CI_policieseffect[11:12,i] = tmp$CI_activeconfirmed_imported[countryconsider,]*0
    }

    CI_policieseffect[13:14,i] = tmp$CI_newcases_change[countryconsider,c(3,7)]/P[countryconsider]
    CI_policieseffect[15:16,i] = tmp$CI_newcases_change[countryconsider,c(4,8)]

    CI_policieseffect[17:18,i] = tmp$CI_activeconfirmed_change[countryconsider,c(3,7)]/P[countryconsider]
    CI_policieseffect[19:20,i] = tmp$CI_activeconfirmed_change[countryconsider,c(4,8)]
    CI_policieseffect[21:22,i] = tmp$CI_Rtchange[countryconsider,c(3,7)]
    CI_policieseffect[23:24,i] = tmp$CI_Rtchange[countryconsider,c(4,8)]


 }

  
#######################################
finaloutput = rbind(travelinbound,travelinbound1, CI_policieseffect)

fname = paste('Policieseffectestimationcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
write.table(finaloutput,file=fname, row.names = F)



}


