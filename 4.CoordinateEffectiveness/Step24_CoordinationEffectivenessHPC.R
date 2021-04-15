#The code is created on Jan 27 to estimate PoliciesEfectiveness to run on ABC for Traffic mode, 4 countries

args <- commandArgs(trailingOnly = TRUE)
iteration <- as.integer(args[1])
policytime <- as.integer(args[2])
inflation <- as.integer(args[3])



inflation = inflation/100 # inflation allow for the proposed method






#setwd("C:/Users/thl902/Desktop/ABC_3countriesTrafficRegulating/3countries")

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


# load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravelSmall/thetas_3travel.RData")
# load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravelSmall/datas_3travel.RData")
# load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravelSmall/travelout_3dat.RData")
# 

# load("thetaarray.RData")
# load("dataarray.RData")
# load("traveldataarray.RData")
# load("initialmatrixarray.RData")
# # 


load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravel/thetaarray.RData")
load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravel/dataarray.RData")
load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravel/traveldataarray.RData")
load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravel/initialmatrixarray.RData")



####
theta0 = thetaarray[[iteration]]
data = dataarray[[iteration]]
travelout_datadivided = traveldataarray[[iteration]]
numbercountries = ncol(travelout_datadivided[[1]] )
initialmat = initialmatrixarray[[iteration]]
trials = 10  # Need to change to 1000

for(countryconsider in 1:numbercountries )
  
{ 
  
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
  travelout_data = matrix(0,nrow=84,ncol=3)
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
  
  # Replace by truth
  #thetahat = theta0
  
  
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
      activeconfirmed_imported[country] =  sum(u1$activeconfirm_importednoquarantine[,a3])
      
      
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
        Aimported_list[[country]][i] = sum(u1$activeconfirm_importednoquarantine[,a3])
        
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
  durationquarantine = rep(14, numbercountries) #14 days quarantine
  travelout_regulated = totaltravelout_samerate_regulated(travelprediction, ratein, P)
  inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
             initialmatrix = initialprediction, quarantinerate = 1, 
             durationquarantine_adjustedout = durationquarantine)
  
  infectmean_policy3 = infectedprediction_mean(thetahat, inp) 
  
  
  infectpercentile_policy3 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp)
  
  
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
  inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
             initialmatrix = initialprediction, quarantinerate = 1, 
             durationquarantine_adjustedout = durationquarantine)
  
  infectmean_policy4 = infectedprediction_mean(thetahat, inp) 
  
  infectpercentile_policy4 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp)
  
  
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
    if(p>=0){
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
          } else{
            ratein[i] = infecttravelin_capacity_divide[i]/ infecttravel_country[i]
          }
        }
        
        #plug in the rate matrix
        a3 = country + (countryconsider-1)*numbercountries
        
        
        travelrate_matrix[a3,] = ratein
        
      }
      
      travelout_regulated = totaltravelout_regulated(travelprediction, travelrate_matrix, P) #total travelers regulated by policy 3
      return(travelout_regulated)
    } else{
      return(list(travel = NA))
    }
  }
  ######
  
  
  threshold = infectmean_policy2$max_unobservedinfected[countryconsider]*(1+inflation) # inlation alllow compare to shutdown
  
  travelout_regulated = averagecontrol_travelout_exactregulated(countryconsider, threshold)
  
  
  
  if(length(travelout_regulated)>1){
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
    
  } else{
    infectpercentile_policy5 = list(travel = NA)
    infectmean_policy5 = list(travel = NA)
    totaltravelers_policy5 = NA
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
    if(p>=0){
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
          } else{
            ratein[i] = infecttravelin_capacity_divide[i]/ infecttravel_country[i]
          }
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
    } else{
      return(list(travel = NA))
    }
  }
  ######
  
  
  threshold = infectmean_policy2$max_unobservedinfected[countryconsider]*(1+inflation) # inlation alllow compare to shutdown
  
  travelout_regulated = averagecontrol_travelout_simplifyregulated(countryconsider, threshold)
  
  
  
  if(length(travelout_regulated)>1){
    durationquarantine = rep(0, numbercountries) #No quarantine needed
    
    inp = list(durationtravel = durationprediction, travelregulated = travelout_regulated,
               initialmatrix = initialprediction, quarantinerate = 1, 
               durationquarantine_adjustedout = durationquarantine)
    
    
    infectpercentile_policy6 = infectedprediction_percentile(trials,percentile_lower ,percentile_upper, thetahat, inp)
    
    infectmean_policy6 = infectedprediction_mean(thetahat, inp) 
    
    totaltravelers_policy6 = 0
    for(timestep in 1:length(travelout_regulated)){
      tmp = travelout_regulated[[timestep]]
      totaltravelers_policy6 = totaltravelers_policy6 + sum(tmp[,countryconsider])
    }
    
  } else{
    infectpercentile_policy6 = list(travel = NA)
    infectmean_policy6 = list(travel = NA)
    totaltravelers_policy6 = NA
  }
  
  
  
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
    
    #Only process when the case is meaningful
    if(p>=0){
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
          } else{
            ratein[i] = infecttravelin_capacity_divide[i]/ infecttravel_country[i]
          }
        }
        
        
        
        
        #plug in the rate matrix
        a3 = country + (countryconsider-1)*numbercountries
        
        
        travelrate_matrix[a3,] = ratein
        
      }
      
      travelout_regulated = totaltravelout_regulated(travelprediction, travelrate_matrix, P) #total travelers regulated by policy 3
      return(travelout_regulated)
    } else{
      return(list(travel = NA))
    }
  }
  ######
  
  
  threshold = infectmean_policy2$max_unobservedinfected[countryconsider]*(1+inflation) # inlation allow compare to shutdown
  
  travelout_regulated = probabilitycontrol_travelout_exactregulated(countryconsider, threshold, probcontrol)
  
  if(length(travelout_regulated)>1){
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
    
    
  } else{
    infectpercentile_policy7 = list(travel = NA)
    infectmean_policy7 = list(travel = NA)
    totaltravelers_policy7 = NA
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
    
    #Only process when the case is meaningful
    if(p>=0){
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
          } else{
            ratein[i] = infecttravelin_capacity_divide[i]/ infecttravel_country[i]
          }
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
    } else{
      return(list(travel = NA))
    }
  }
  ######
  
  
  
  threshold = infectmean_policy2$max_unobservedinfected[countryconsider]*(1+inflation) # inlation allow compare to shutdown
  
  
  travelout_regulated = probabilitycontrol_travelout_simplifyregulated(countryconsider, threshold, probcontrol)
  
  if(length(travelout_regulated) > 1){
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
    
    
  } else{
    infectpercentile_policy8 = list(travel = NA)
    infectmean_policy8 = list(travel = NA)
    totaltravelers_policy8 = NA
  }
  
  
  
  
  
  #################################3. REPORT########################
  ######1. TOTAL TRAVELERS
  #We create a matrix of 3 rows, and 8 columns, each column is corresponding to one policy
  #1st row, number people can travel
  #2nd row, percentage compare to full travel
  #3rd row, expected number travelers due to the restriction
  
  travelmatrix = matrix(0,3,8)
  
  
  for(i in 1:8){
    
    
    tmp = get( paste("totaltravelers_policy", i, sep="") ) 
    
    if (is.na(tmp) == 0){
      travelmatrix[1,i] =  tmp
      
    } else{
      travelmatrix[1,i] =  NA
    }
  }
  travelmatrix[2,] = round(travelmatrix[1,]/travelmatrix[1,1], digits = 2)
  travelmatrix[3,] = travelmatrix[1,]
  travelmatrix[3,3:4] =  c(totaltravelers_policy3_adjusted, totaltravelers_policy4_adjusted) #replaced by the expected values
  travelmatrix[3,] = round(travelmatrix[3,]/travelmatrix[1,1], digits =2)
  totaltravelers = travelmatrix[1,1] # total travel in the country considered under the full load
  
  
  ####################################
  # ## 2. MEAN COMPARE, STANDARDIZE
  #We return a matrix of 8 columns, each column corresponding to one policy
  #1st row, percent unobserved infected enter the country prequarantine
  #2nd row, percent unobserved infected enter the country post quarantine
  #3rd row, percent total active confirmed imported by travel
  #4th row, relative change in new active confirmed
  #5th row, relative change in new cases change
  #6th row, reduction in Rt, magnitude
  #7th row, reduction in Rt, relative
  
  mean_policieseffect = matrix(0,7,8)
  
  for(i in 1:8){
    
    
    tmp = get( paste("infectmean_policy", i, sep="") ) 
    
    if(length(tmp) > 1){
      
      mean_policieseffect[1,i] =  tmp$unobservedinfected_imported_prequarantine[countryconsider]/totaltravelers 
      mean_policieseffect[2,i] =  tmp$unobservedinfected_imported_postquarantine[countryconsider]/totaltravelers 
      
      mean_policieseffect[3,i] =  tmp$activeconfirmed_imported[countryconsider]/totaltravelers 
      
      mean_policieseffect[4,i] =  tmp$activeconfirmed_change[countryconsider,4]
      mean_policieseffect[5,i] =  tmp$newcases_change[countryconsider,4]
      mean_policieseffect[6,i] =  tmp$Rtchange[countryconsider,3]
      mean_policieseffect[7,i] =  tmp$Rtchange[countryconsider,4]
      
      
      
      
      
      
    } else{
      mean_policieseffect[,i] =  NA
      
    }
  }
  #####################
  # ## 2. MEAN COMPARE, DIFFERENCE
  #We return a matrix of 8 columns, each column corresponding to one policy
  #1st row, unobserved infected enter the country prequarantine
  #2nd row, unobserved infected enter the country post quarantine
  #3rd row, total active confirmed imported by travel
  #4th row, difference change in new active confirmed
  #5th row, difference change in new cases change
  #6th row, Rt beginning
  #7th row, Rt end of policies
  
  
  
  
  mean_policieseffect_difference = matrix(0,7,8)
  
  for(i in 1:8){
    
    
    tmp = get( paste("infectmean_policy", i, sep="") ) 
    
    
    
    if(length(tmp) > 1){
      mean_policieseffect_difference[1,i] =  tmp$unobservedinfected_imported_prequarantine[countryconsider]
      mean_policieseffect_difference[2,i] =  tmp$unobservedinfected_imported_postquarantine[countryconsider]
      
      mean_policieseffect_difference[3,i] =  tmp$activeconfirmed_imported[countryconsider] 
      
      mean_policieseffect_difference[4,i] =  tmp$activeconfirmed_change[countryconsider,3]
      mean_policieseffect_difference[5,i] =  tmp$newcases_change[countryconsider,3]
      mean_policieseffect_difference[6,i] =  tmp$Rtchange[countryconsider,1]
      mean_policieseffect_difference[7,i] =  tmp$Rtchange[countryconsider,2]
      
      
      
      
    } else{
      mean_policieseffect_difference[,i] =  NA
      
    }
  }
  
  
  
  #############
  
  
  
  
  ####################################
  # ## 2. CI  COMPARE
  #We return a matrix of 8 columns, each column corresponding to one policy
  #1th and 2nd row, lower and upper percentile percent unobserved infected enter the country prequarantine
  #3rd and 4th row, lower and upper percentile percent unobserved infected enter the country post quarantine
  #5th and 6tth row, lower and upper percentile percent total active confirmed imported by travel
  #7th and 8th row, lower and upper percentile relative change in new active confirmed
  #9th and 10th row, lower and upper percentile relative change in new cases change
  #11th and 12th row, lower and upper percentile reduction in Rt, magnitude
  #13th and 14th row, lower and upper percentile reduction in Rt, relative
  
  CI_policieseffect = matrix(0,14,8)
  
  for(i in 1:8){
    
    
    tmp = get( paste("infectpercentile_policy", i, sep="") ) 
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
    
    
    tmp = get( paste("infectpercentile_policy", i, sep="") ) 
    
    
    
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
  finaloutput = rbind(travelmatrix,mean_policieseffect,mean_policieseffect_difference,CI_policieseffect,CI_policieseffect_difference)
  
  fname = paste('Policieseffectestimationcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")
  write.table(finaloutput,file=fname, row.names = F)
  
}
