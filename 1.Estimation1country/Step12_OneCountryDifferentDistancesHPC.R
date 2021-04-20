
args <- commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1])
prob1 <- as.integer(args[2])

prob =prob1/1000
load("/n/home08/tle/ABCConet/thetas_list.RData")
load("/n/home08/tle/ABCConet/datas_list.RData")
###########
library(protoABC)
library(ggplot2)
library(dplyr)
library(tidyr)

########
tmp = proc.time()
d = 84
n1 = 5
particles = 1000 ##Need to change to 1000###
ac1s = ac2s = ac3s = matrix(0,n1,8)
ac4s = ac5s = ac6s = matrix(0,n1,8)
ac7s = ac8s = ac9s= matrix(0,n1,8)


theta_cover1 = theta_cover2 = theta_cover3 = theta_cover4 = theta_cover5 = matrix(0,n1,4)
theta_cover6 = theta_cover7 = theta_cover8 = theta_cover9 = theta_truth = matrix(0,n1,4)


IQR1s = IQR2s = IQR3s = IQR4s = IQR5s = IQR6s = IQR7s = IQR8s = IQR9s = matrix(0,n1,4)





P =10^7
I1 = 15
A1 = 13
S1 = P-I1-A1
x1 = c(S1,I1,A1,0,0,0)# State corresponding S,I,A,R,D,Ru


for(j in 1:n1){
  j1 = j + k*5
  
  theta0 = as.numeric(thetas_list[[j1]])
  data = datas_list[[j1]]
  mydat1 = data[,3:5]
  
  nrep = nrow(mydat1)
  drecover = mydat1[,2] - c(0,mydat1[,2][-length(mydat1[,2])])
  ddeath = mydat1[,3] - c(0,mydat1[,3][-length(mydat1[,3])])
  active = c(1,mydat1[,1][-length(mydat1[,1])])
  
  index0 = which(active==0) 
  index = c(1,index0)
  drecover = drecover[-index]
  ddeath =  ddeath[-index]
  active = active[-index]
  
  
  betas = drecover/active
  deltas = ddeath/active
  
  inp <- list(x1 = x1, data = mydat1, nrep = nrep, beta = betas, delta = deltas, P=P)
  
  status_matrix1 = matrix(0,nrow = d  ,ncol=6)
  status_matrix1[1,] = x1
  
  #########################################################################
  ##########STOCHASTIC MODEL######################################
  #We define a generating function stochastic driven#############
  
  sto_generating_func = function(theta,inp){
    status_matrix1 = matrix(0,nrow = inp$nrep,ncol=6)
    status_matrix1[1,] = inp$x1
    
    ##################Deefining harzard functions
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
    
    for (i in 2:inp$nrep){
      
      x = status_matrix1[(i-1),]
      
      
      
      
      ##Generating Poisson values based on hazard functions
      y1 = rpois(1, harzard1(x,theta))
      
      y2 =  rpois(1, harzard2(x,theta))
      
      
      y3 = rpois(1, harzard3(x,theta))
      
      
      y4 =  rpois(1, harzard4(x,theta))
      
      
      
      y5 = rpois(1, harzard5(x,theta))
      
      
      
      
      if(y1>x[1]){y1=x[1]}
      if(y1-y2-y5+x[2]<0){y2=x[2]+y1-y5}
      if(y2-y3-y4+x[3]<0){y3=x[3]+y2-y4}
      
      y=c(y1, y2, y3, y4, y5)
      ##Transition matrix
      tran_mat = matrix(c(-1,1,0,0,0,0,0,-1,1,0,0,0,0,0,-1,1,0,0,0,0,-1,0,1,0,0,-1,0,0,0,1),nrow=6,ncol=5)
      ##Updating values
      val = tran_mat%*%y
      x = x+ t(val)
      status_matrix1[i,] = x
    }
    
    return(status_matrix1)
  }
  
  
  ################################Distance 1, uniform prior, with calibration for beta and delta########
  ########1a.##############
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
  #########Then prior density
  prior_eval <- function(theta){
    prior_value <-  dunif(theta["alpha"], 0, 2)*dunif(theta["beta"], 0, 1)*dunif(theta["delta"], 0, 1)*dunif(theta["gamma"], 0, 1)
    
    return(prior_value)
  }
  
  ##########
  
  
  
  ##############################################################################
  ############3. Learn alpha based on available data, beta, and gamma##
  
  
  alpha_function = function(data,theta,P){
    beta = theta[3]
    gamma = theta[6]
    Ut = data[,1] + data[,2] + data[,3]
    dUt = c(Ut[-1],0)- Ut
    # #########################
    # 
    infected= rep(0,nrow(data))
    dRus = rep(0,nrow(data))
    for (i in 1:nrow(data)){
      infected[i] = dUt[i]/gamma
      dRus[i] = dUt[i]/gamma*beta
    }
    
    Ruseq1 = cumsum(dRus)
    Ruseq = c(0,Ruseq1[-length(Ruseq1)])
    Sseq = rep(P,nrow(data)) - Ut - Ruseq - infected
    
    dSseq = (Sseq - c(0,Sseq[-length(Sseq)]))*(-1)
    SIseq = Sseq*infected
    SIseq = c(1, SIseq[-length(SIseq)])
    index1 = which(SIseq==0)
    index2 = c(1,index1,length(SIseq))
    dSseq = dSseq[-index2]
    SIseq = SIseq[-index2]
    
    alphas_sim = dSseq/SIseq*P
    ########This help to avoid empty array
    if (length(alphas_sim) == 0){alphas_sim = 10^8}
    ###########
    
    
    return(alphas_sim)
    
  }
  
  
  
  
  #########Learn beta functions
  beta_function = function(data){
    drecover_sim = data[,2] - c(0,data[,2][-length(data[,2])])
    active_sim = c(1,data[,1][-length(data[,1])])
    
    index5 = which(active_sim==0) 
    index6 = c(1,index5)
    drecover_sim = drecover_sim[-index6]
    active_sim = active_sim[-index6]
    
    
    betas_sim = drecover_sim/active_sim
    return(betas_sim)
    
  }
  ########
  #########Learn delta functions
  delta_function = function(data){
    ddeath_sim = data[,3] - c(0,data[,3][-length(data[,3])])
    active_sim = c(1,data[,1][-length(data[,1])])
    
    index7 = which(active_sim==0) 
    index8 = c(1,index7)
    ddeath_sim =  ddeath_sim[-index8]
    active_sim = active_sim[-index8]
    
    deltas_sim = ddeath_sim/active_sim
    return(deltas_sim)
    
  }
  ########
  
  
  ########################################
  ############################ Traditional Distances
  #########################################
  #######Parallel Claim Environment#####
  
  
  relative_rate = function(vec1,vec2){
    dvec = vec1 - vec2
    index = which(vec2 == 0)
    if (length(index) != 0){
      vec2 = vec2[-index]
      dvec = dvec[-index]
    }
    rrate = dvec/vec2
    return(rrate)
  }
  
  relativemax_rate = function(vec1,vec2){
    dvec = vec1 - vec2
    rrate = dvec/max(vec2)
    return(rrate)
  }
  ##################Recover the average path based on the observed data
  recover_function = function(data,theta,P){
    beta = theta[3]
    gamma = theta[6]
    Ut = data[,1] + data[,2] + data[,3]
    dUt = c(Ut[-1],0)- Ut
    # #########################
    # 
    recovermat = matrix(0,nrow = nrow(data),ncol=6)
    infected= rep(0,nrow(data))
    dRus = rep(0,nrow(data))
    for (i in 1:nrow(data)){
      infected[i] = dUt[i]/gamma
      dRus[i] = dUt[i]/gamma*beta
    }
    
    Ruseq1 = cumsum(dRus)
    
    Ruseq = c(0,Ruseq1[-length(Ruseq1)])
    Sseq = rep(P,nrow(data)) - Ut - Ruseq - infected
    
    recovermat = data.frame(Sseq, infected, data[,1], data[,2], data[,3], Ruseq)
    
    dSseq = (Sseq - c(0,Sseq[-length(Sseq)]))*(-1)
    SIseq = Sseq*infected
    SIseq = c(1, SIseq[-length(SIseq)])
    index1 = which(SIseq==0)
    index2 = c(1,index1,length(SIseq))
    dSseq = dSseq[-index2]
    SIseq = SIseq[-index2]
    
    alphas_sim = dSseq/SIseq*P
    ########This help to avoid empty array
    if (length(alphas_sim) == 0){alphas_sim = 10^8}
    ###########
    
    
    return(recovermat)
    
  }
  
  #################the deterministic model, corresponding to the average realization##############
  
  generating_func = function(theta,inp){
    status_matrix = matrix(0,nrow = inp$nrep,ncol=6)
    status_matrix[1,] = inp$x1
    
    ##################Deefining harzard functions
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
    
    for (i in 2:inp$nrep){
      
      x = status_matrix[(i-1),]
      
      
      y1 = harzard1(x,theta)
      
      y2 =   harzard2(x,theta)
      
      
      y3 = harzard3(x,theta)
      
      
      y4 =  harzard4(x,theta)
      
      
      
      y5 = harzard5(x,theta)
      
      if(y1>x[1]){y1=x[1]}
      if(y1-y2-y5+x[2]<0){y2=x[2]+y1-y5}
      if(y2-y3-y4+x[3]<0){y3=x[3]+y2-y4}
      
      y=c(y1, y2, y3, y4, y5)
      ##Transition matrix
      tran_mat = matrix(c(-1,1,0,0,0,0,0,-1,1,0,0,0,0,0,-1,1,0,0,0,0,-1,0,1,0,0,-1,0,0,0,1),nrow=6,ncol=5)
      ##Updating values
      val = tran_mat%*%y
      x = x+ t(val)
      status_matrix[i,] = x
    }
    
    return(status_matrix)
  }
  
  ##########################Standardized##########
  
  standardize_seq = function(k,inp){
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
    sd = rep(0, inp$nrep) 
    mylist = list()
    Smat = matrix(0,k,inp$nrep)
    Imat = matrix(0,k,inp$nrep)
    Amat = matrix(0,k,inp$nrep)
    Rmat = matrix(0,k,inp$nrep)
    Dmat = matrix(0,k,inp$nrep)
    Rumat = matrix(0,k,inp$nrep)
    for(i in 1:k){
      theta = as.numeric(para[i,])
      u = sto_generating_func(theta,inp)
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
  sd = standardize_seq(3000,inp)[,-1]
  
  inp <- list(x1 = x1, data = mydat1, nrep = nrep, beta = betas, delta = deltas, P=P, sd = sd)
  
  
  sd_function1 = function(vec1, vec2,inp){
    sd = inp$sd
    output1 <- sum(((vec1[-1]-vec2[-1])/t(sd[3,]))^2)
    return(output1)
    
  }
  sd_function2 = function(vec1, vec2,inp){
    sd = inp$sd
    output1 <- sum(((vec1[-1]-vec2[-1])/t(sd[4,]))^2)
    return(output1)
    
  }
  sd_function3 = function(vec1, vec2,inp){
    sd = inp$sd
    output1 <- sum(((vec1[-1]-vec2[-1])/t(sd[5,]))^2)
    return(output1)
    
  }
  
  #####################
  library(parallel)
  cl <- makeCluster(4)
  clusterExport(cl,c("sto_generating_func","generating_func","sd_function1", "sd_function2", "sd_function3",
                     "recover_function","relative_rate","relativemax_rate","alpha_function","beta_function","delta_function"))
  
  ########BASED ON PARAMETERS##########
  distance1 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    sim_u = rowSums(sim[,3:5])
    ########## Learn beta, delta
    mydat2 = sim[,3:5]
    
    betas_sim = beta_function(mydat2)
    deltas_sim = delta_function(mydat2)
    
    ##########Learn alpha#######
    alphas = alpha_function(inp$data,theta,inp$P)
    alphas_sim = alpha_function(mydat2,theta,inp$P)
    #############
    output1 <- abs(median(betas_sim) - median(inp$beta))
    output2 <-  abs(median(deltas_sim) - median(inp$delta))
    output3 <- abs(median(alphas_sim) - median(alphas))
    output <- (output1 +output2+output3)^.5
    return(output)
  }
  
  
  
  abc_post_1 <- abc_start(
    prior,
    distance1,cl=cl,
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = FALSE)
  )
  
  
  
  # ##### Euclidean A(t)
  
  distance2 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    #####
    r_active = sum((sim[,3] - inp$data[,1])^2)
    r_recover = sum((sim[,4] - inp$data[,2])^2)
    r_death = sum((sim[,5]-inp$data[,3])^2)
    output <- sqrt(1/nrow(inp$data)*(r_active + r_recover + r_death ))
    
    return(output)
  }
  
  
  abc_post_2 <- abc_start(
    prior,
    distance2, cl=cl,
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = FALSE)
  )
  
  # ##########################   STANDARDDIZE DISTANCES, A(t)   #################
  
  
  distance3 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    
    
    
    output1 <- sd_function1(inp$data[,1],sim[,3],inp)
    output2 <- sd_function2(inp$data[,2],sim[,4],inp)
    output3 <- sd_function3(inp$data[,3],sim[,5],inp)
    output <- sqrt(1/nrow(inp$data)*(output1 + output2+ output3))
    
    
    
    return(output)
  }
  
  abc_post_3 <- abc_start(
    prior,
    distance3,cl=cl,
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = FALSE)
  )
  
  
  ## RELATIVE DISTANCE#####################
  
  distance4 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    sim_u = rowSums(sim[,3:5])
    #####
    r_active = relative_rate(sim[,3], inp$data[,1])
    r_recover = relative_rate(sim[,4], inp$data[,2])
    r_death = relative_rate(sim[,5], inp$data[,3])
    output <- sqrt(1/nrow(inp$data)*(sum(r_active^2) + sum(r_recover^2) + sum(r_death^2 )))
    
    return(output)
  }
  
  abc_post_4 <- abc_start(
    prior,
    distance4, cl=cl,
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = FALSE)
  )
  
  ##   Relative max, A(t)   #################
  
  
  
  distance5 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    sim_u = rowSums(sim[,3:5])
    #####
    r_active = relativemax_rate(sim[,3], inp$data[,1])
    r_recover = relativemax_rate(sim[,4], inp$data[,2])
    r_death = relativemax_rate(sim[,5], inp$data[,3])
    output <- sqrt(1/nrow(inp$data)*(sum(r_active^2) + sum(r_recover^2) + sum(r_death^2 )))
    
    return(output)
  }
  
  
  
  abc_post_5 <- abc_start(
    prior,
    distance5, cl=cl,
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = FALSE)
  )
  
  #######sqrt(Relative Rate) x sqrt(Parameters)######################## 
  
  distance6 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    sim_u = rowSums(sim[,3:5])
    ########## Learn beta, delta
    mydat2 = sim[,3:5]
    
    betas_sim = beta_function(mydat2)
    deltas_sim = delta_function(mydat2)
    
    ##########Learn alpha#######
    
    
    alphas = alpha_function(inp$data,theta,inp$P)
    alphas_sim = alpha_function(mydat2,theta,inp$P)
    
    #############
    output1 <- abs(median(betas_sim) - median(inp$beta))
    output2 <-  abs(median(deltas_sim) - median(inp$delta))
    output3 <- abs(median(alphas_sim) - median(alphas))
    ############
    r_active = relative_rate(sim[,3], inp$data[,1])
    r_recover = relative_rate(sim[,4], inp$data[,2])
    r_death = relative_rate(sim[,5], inp$data[,3])
    output4 <- (1/nrow(inp$data)*(sum(r_active^2) + sum(r_recover^2) + sum(r_death^2 )))
    #############
    
    output <- (output1 +output2+output3)^.5*output4^.5
    return(output)
  }
  
  
  ##
  
  abc_post_6 <- abc_start(
    prior,
    distance6,cl=cl,
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = FALSE)
  )
  
  ###############
  ##\sqrt(RelativeRate) + \sqrt(Parameter)
  
  
  distance7 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    sim_u = rowSums(sim[,3:5])
    ########## Learn beta, delta
    mydat2 = sim[,3:5]
    
    betas_sim = beta_function(mydat2)
    deltas_sim = delta_function(mydat2)
    
    ##########Learn alpha#######
    
    
    alphas = alpha_function(inp$data,theta,inp$P)
    alphas_sim = alpha_function(mydat2,theta,inp$P)
    
    #############
    output1 <- abs(median(betas_sim) - median(inp$beta))
    output2 <-  abs(median(deltas_sim) - median(inp$delta))
    output3 <- abs(median(alphas_sim) - median(alphas))
    ############
    r_active = relative_rate(sim[,3], inp$data[,1])
    r_recover = relative_rate(sim[,4], inp$data[,2])
    r_death = relative_rate(sim[,5], inp$data[,3])
    output4 <- (1/nrow(inp$data)*(sum(r_active^2) + sum(r_recover^2) + sum(r_death^2 )))
    #############
    
    output <- (output1 +output2+output3)^.5 + output4^.5
    return(output)
  }
  #########
  
  abc_post_7 <- abc_start(
    prior,
    distance7,cl=cl,
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = FALSE)
  )
  
  
  ################
  ##\sqrt(Sd) x \sqrt(Parameter)
  
  
  distance8 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    sim_u = rowSums(sim[,3:5])
    ########## Learn beta, delta
    mydat2 = sim[,3:5]
    
    betas_sim = beta_function(mydat2)
    deltas_sim = delta_function(mydat2)
    
    ##########Learn alpha#######
    
    
    alphas = alpha_function(inp$data,theta,inp$P)
    alphas_sim = alpha_function(mydat2,theta,inp$P)
    
    #############
    output1 <- abs(median(betas_sim) - median(inp$beta))
    output2 <-  abs(median(deltas_sim) - median(inp$delta))
    output3 <- abs(median(alphas_sim) - median(alphas))
    ############
    output4 <- sd_function1(inp$data[,1],sim[,3],inp)
    output5 <- sd_function2(inp$data[,2],sim[,4],inp)
    output6 <- sd_function3(inp$data[,3],sim[,5],inp)
    #############
    
    output <- (output1 +output2+output3)^.5*(1/nrow(inp$data)*(output4 + output5+ output6))^.5
    return(output)
  }
  
  ##
  abc_post_8 <- abc_start(
    prior,
    distance8,cl=cl,
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = FALSE)
  )
  
  ###############
  
  
  ##\sqrt(Sd) + \sqrt(Parameter)############
  
  
  distance9 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    sim_u = rowSums(sim[,3:5])
    ########## Learn beta, delta
    mydat2 = sim[,3:5]
    
    betas_sim = beta_function(mydat2)
    deltas_sim = delta_function(mydat2)
    
    ##########Learn alpha#######
    
    
    alphas = alpha_function(inp$data,theta,inp$P)
    alphas_sim = alpha_function(mydat2,theta,inp$P)
    
    #############
    output1 <- abs(median(betas_sim) - median(inp$beta))
    output2 <-  abs(median(deltas_sim) - median(inp$delta))
    output3 <- abs(median(alphas_sim) - median(alphas))
    ############
    output4 <- sd_function1(inp$data[,1],sim[,3],inp)
    output5 <- sd_function2(inp$data[,2],sim[,4],inp)
    output6 <- sd_function3(inp$data[,3],sim[,5],inp)
    #############
    
    output <- (output1 +output2+output3)^.5 + (1/nrow(inp$data)*(output4 + output5+ output6))^.5
    return(output)
  }
  
  ##
  
  abc_post_9 <- abc_start(
    prior,
    distance9, cl=cl,
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = FALSE)
  )
  
  #####################Accuracy median
  ac_med_func = function(data,theta){
    val11 = (median(data$alpha) - theta[2])/theta[2]
    val12 =  (median(data$beta) - theta[3])/theta[3]
    val13 = (median(data$delta) - theta[4])/theta[4]
    val14 = (median(data$gamma) - theta[6])/theta[6]
    return(c(abs(val11), abs(val12), abs(val13), abs(val14)))
  }
  
  
  ######Accuracy mean
  ac_mean_func = function(data,theta){
    val11 = (mean(data$alpha) - theta[2])/theta[2]
    val12 =  (mean(data$beta) - theta[3])/theta[3]
    val13 = (mean(data$delta) - theta[4])/theta[4]
    val14 = (mean(data$gamma) - theta[6])/theta[6]
    return(c(abs(val11), abs(val12), abs(val13), abs(val14)))
  }
  
  ac1s[j,] = c(ac_med_func(abc_post_1,theta0), ac_mean_func(abc_post_1,theta0))
  ac2s[j,] = c(ac_med_func(abc_post_2,theta0), ac_mean_func(abc_post_2,theta0))
  ac3s[j,] = c(ac_med_func(abc_post_3,theta0), ac_mean_func(abc_post_3,theta0))
  ac4s[j,] = c(ac_med_func(abc_post_4,theta0), ac_mean_func(abc_post_4,theta0))
  ac5s[j,] = c(ac_med_func(abc_post_5,theta0), ac_mean_func(abc_post_5,theta0))
  ac6s[j,] = c(ac_med_func(abc_post_6,theta0), ac_mean_func(abc_post_6,theta0))
  ac7s[j,] = c(ac_med_func(abc_post_7,theta0), ac_mean_func(abc_post_7,theta0))
  ac8s[j,] = c(ac_med_func(abc_post_8,theta0), ac_mean_func(abc_post_8,theta0))
  ac9s[j,] = c(ac_med_func(abc_post_9,theta0), ac_mean_func(abc_post_9,theta0))
  
  ########Coverage######
  coverage_func = function(data,theta){
    val11 = (quantile(data$alpha,.25) - theta[2])*(quantile(data$alpha,.75) - theta[2])
    val12 =  (quantile(data$beta,.25) - theta[3])*(quantile(data$beta,.75) - theta[3])
    val13 = (quantile(data$delta,.25) - theta[4])*(quantile(data$delta,.75) - theta[4])
    val14 = (quantile(data$gamma,.25) - theta[6])*(quantile(data$gamma,.75) - theta[6])
    #########
    if(val11 <0){ a1 = 1
    } else {a1 = 0 }
    
    if(val12 <0){ a2 = 1
    } else {a2 = 0 }
    
    if(val13 <0){ a3 = 1
    } else {a3 = 0 }
    
    if(val14 <0){ a4 = 1
    } else {a4 = 0 }
    ######
    
    return(c(a1, a2,a3,a4))
  }
  
  theta_cover1[j,] = coverage_func(abc_post_1,theta0)
  
  theta_cover2[j,] = coverage_func(abc_post_2,theta0)
  
  theta_cover3[j,] = coverage_func(abc_post_3,theta0)
  
  theta_cover4[j,] = coverage_func(abc_post_4,theta0)
  
  theta_cover5[j,] = coverage_func(abc_post_5,theta0)
  
  theta_cover6[j,] = coverage_func(abc_post_6,theta0)
  
  theta_cover7[j,] = coverage_func(abc_post_7,theta0)
  
  theta_cover8[j,] = coverage_func(abc_post_8,theta0)
  
  theta_cover9[j,] = coverage_func(abc_post_9,theta0)
  
  #########IQR
  IQR_func = function(data,theta){
    val11 = quantile(data$alpha,.75) - quantile(data$alpha,.25)
    val12 =  quantile(data$beta,.75) - quantile(data$beta,.25)
    val13 = quantile(data$delta,.75) - quantile(data$delta,.25)
    val14 = quantile(data$gamma,.75) - quantile(data$gamma,.25)
    
    return(c(val11,val12,val13,val14))
  }
  
  IQR1s[j,] = IQR_func(abc_post_1,theta0)
  IQR2s[j,] = IQR_func(abc_post_2,theta0)
  IQR3s[j,] = IQR_func(abc_post_3,theta0)
  
  IQR4s[j,] = IQR_func(abc_post_4,theta0)
  IQR5s[j,] = IQR_func(abc_post_5,theta0)
  IQR6s[j,] = IQR_func(abc_post_6,theta0)
  
  IQR7s[j,] = IQR_func(abc_post_7,theta0)
  IQR8s[j,] = IQR_func(abc_post_8,theta0)
  IQR9s[j,] = IQR_func(abc_post_9,theta0)
  
  
  theta_truth[j,] = theta0[-c(1,5)]
  
  print(j1)
  
  
}

output1 = cbind(ac1s, IQR1s, theta_cover1, theta_truth)

output2 = cbind(ac2s, IQR2s, theta_cover2, theta_truth)

output3 = cbind(ac3s, IQR3s, theta_cover3, theta_truth)

output4 = cbind(ac4s, IQR4s, theta_cover4, theta_truth)

output5 = cbind(ac5s, IQR5s, theta_cover5, theta_truth)

output6 = cbind(ac6s, IQR6s, theta_cover6, theta_truth)

output7 = cbind(ac7s, IQR7s, theta_cover7, theta_truth)

output8 = cbind(ac8s, IQR8s, theta_cover8, theta_truth)

output9 = cbind(ac9s, IQR9s, theta_cover9, theta_truth)



output = rbind(output1,output2,output3,output4,output5,output6,output7,output8,output9)

fname =paste("Distanceoutput_particles",particles,"_pacc",prob,"_loop",k,".txt",sep="")
write.table(output,file=fname, row.names = F)

proc.time() -tmp
