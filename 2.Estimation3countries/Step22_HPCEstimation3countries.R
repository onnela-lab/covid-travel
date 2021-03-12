#The code is created on Oct 23 to run on ABC for Traffic mode, 3 countries
args <- commandArgs(trailingOnly = TRUE)
k <- as.integer(args[1])
prob1 <- as.integer(args[2])

prob =prob1/1000


library(protoABC)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rlist)
tmp = proc.time()



load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravelSmall/thetas_3travel.RData")
load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravelSmall/datas_3travel.RData")
load("/n/holyscratch01/onnela_lab/ThienLe/ABCTravelSmall/travelout_3dat.RData")



particles = 1000 # Need to change to 1000



##################
theta0 = thetas_3travel[[k]]
data = datas_3travel[[k]]
traveloutdat = travelout_3dat

#######Seperate 3 countries
mydat1 = data[,3:5]
mydat2 = data[,9:11]
mydat3 = data[,15:17]
theta01 = theta0[1:6]
theta02 = theta0[7:12]
theta03 = theta0[13:18]
##############
##########
d=84


## COUNTRY 1, US
P1 = 10^7
I1 = 250
A1 = 130
S1 = P1 - I1 - A1
x1 = c(S1,I1,A1,0,0,0)# State corresponding S,I,A,R,D,Ru

#############
## COUNTRY 2, ITA
P2 = 3*10^6
I2 = 20
A2 = 10
S2 = P2 - I2 - A2
x2 = c(S2,I2,A2,0,0,0)
###################
## COUNTRY 3, CAN
P3 = 2*10^6
I3 = 15
A3 = 15
S3 = P3 - I3 - A3
x3 = c(S3,I3,A3,0,0,0)
###################################
initial_corona = c(x1,x2,x3)
#############
##Reconstruct number infected from each country
##Define function to recover infected##
## First way based on each time step
infect_func = function(U,x){
  diff1 = diff(U, differences = 1)
  diffa = diff1[-1]
  diffb = diff1[-length(diff1)]
  infect = rep(0,length(diff1))
  
  ###############
  infect[1] = x[2]
  
  for (i in 2:length(diff1)){
    i1 = i-1
    tmp = infect[i1]
    
    
    if (diffb[i1]>.4){
      infect[i] = tmp*diffa[i1]/diffb[i1]
      
    } else{
      infect[i] = diff1[i]
      
      
    }
    
  }
  
  return (infect)
  
}


##########
U1 = t(t(rowSums(mydat1)))

U2 = t(t(rowSums(mydat2)))

U3 = t(t(rowSums(mydat3)))

infect1 = infect_func(U1,x1)
infect1 = c(infect1,0)
infect2 = infect_func(U2,x2)
infect2 = c(infect2,0)
infect3 = infect_func(U3,x3)
infect3 = c(infect3,0)
###############Reconstruct Total each time point
in1 = traveloutdat[,2]*P1/(P1+P3) + traveloutdat[,3]*P1/(P1+P2)
out1 = traveloutdat[,1]

in2 = traveloutdat[,1]*P2/(P2+P3) + traveloutdat[,3]*P2/(P1+P2)
out2 = traveloutdat[,2]

in3 = traveloutdat[,1]*P3/(P2+P3) + traveloutdat[,2]*P3/(P1+P3)
out3 = traveloutdat[,3]


############

country1 = list(P= P1, data = mydat1, nrep= nrow(mydat1), total_in = in1, total_out = out1, infect =infect1)
country2 = list(P= P2, data = mydat2,  nrep= nrow(mydat2), total_in = in2, total_out = out2, infect =infect2)
country3 = list(P= P3, data = mydat3,  nrep= nrow(mydat3), total_in = in3, total_out = out3, infect =infect3)



pop_func = function(country){
  fluc1 = country$total_in - country$total_out
  pop = rep(0,d)
  pop[1] = country$P
  
  for(i in 2:d){
    x = pop[i-1]
    pop[i] = x+ fluc1[i] #update total population
  }
  return(pop)
  
}

pop1 = pop_func(country1)


pop2 = pop_func(country2)


pop3 = pop_func(country3)


#########Construct Travel Compartments
country1 = list.append(country1, pop = pop1)
country2 = list.append(country2, pop = pop2)
country3 = list.append(country3,  pop = pop3)
######
##########

####### Define a function to reconstruct the travel out mechanism from each country 




travelout_func = function(total,country)
{
  
  # country = country2
  # total= travel21
  mydat = country$data
  pop = country$pop
  infect = country$infect
  
  
  #######
  
  drecover = diff(mydat[,2], differences = 1)
  
  active = mydat[,1][-length(mydat[,1])]
  
  index0 = which(active<=.1) 
  
  if (length(index0)!=0){
    drecover = drecover[-index0]
    active = active[-index0]
  }
  
  betas = drecover/active
  
  
  Ut = rowSums(mydat)
  #Difference Ut
  dUt = c(diff(Ut, differences = 1),0)
  ########
  dRus = rep(0,nrow(mydat))
  
  for (i in 1:(nrow(mydat)-1)){
    
    dRus[i] = infect[i]*median(betas)
    
  }
  
  
  Ruseq1 = cumsum(dRus)
  
  Ruseq = c(0,Ruseq1[-length(Ruseq1)])
  
  
  h1 =  Ruseq - data[,12]
  max(abs(h1))
  
  Sseq = pop - Ut - Ruseq - infect
  
  
  
  compartments = cbind(Sseq, infect, Ruseq)
  
  
  update = matrix(0,ncol=6,nrow=nrow(mydat))
  
  for (i in 2:nrow(mydat)){
    i1= i-1
    x = compartments[i1,]
    
    ####
    if(min(x)>0){
      update[i,] = c(total[i]*x[1]/(x[1]+x[2]), total[i]*x[2]/(x[1]+x[2]), 0, 0, 0, 0)
    } else{
      update[i,] = c(total[i], 0, 0, 0, 0, 0)
    }
    
  }
  
  
  return(update)
}



######################


travel12 = traveloutdat[,1]*P2/(P2+P3)

travelout12 = travelout_func(travel12,country1)


travel13 = traveloutdat[,1]*P3/(P2+P3)

travelout13 = travelout_func(travel13,country1)

travelout1 = travelout12 + travelout13
##########


travel21 = traveloutdat[,2]*P1/(P1+P3)

travelout21 = travelout_func(travel21, country2)


travel23 = traveloutdat[,2]*P3/(P1+P3)

travelout23 = travelout_func(travel23,country2)

travelout2 = travelout21 + travelout23



###########
travel31 = traveloutdat[,3]*P1/(P1+P2)

travelout31 = travelout_func(travel31,country3)

travel32= traveloutdat[,3]*P2/(P1+P2)

travelout32 = travelout_func(travel32,country3)

travelout3 = travelout31 + travelout32 

#################Travel in sequences
travelin1 = travelout21 + travelout31

travelin2 = travelout12 + travelout32

travelin3 = travelout13 + travelout23

###########
########Define an ABC_function to estimate each country
country1 = list.append(country1, x_ini = x1, travelin = travelin1, travelout = travelout1, thetatruth = theta01)

country2 = list.append(country2, x_ini = x2, travelin = travelin2, travelout = travelout2, thetatruth = theta02)

country3 = list.append(country3, x_ini = x3, travelin = travelin3, travelout = travelout3, thetatruth = theta03)


######################

ABC_function = function(country){
  mydat = country$data
  drecover = mydat[,2] - c(0,mydat[,2][-length(mydat[,2])])
  ddeath = mydat[,3] - c(0,mydat[,3][-length(mydat[,3])])
  active = c(1,mydat[,1][-length(mydat[,1])])
  index0 = which(active==0) 
  index = c(1,index0)
  drecover = drecover[-index]
  ddeath = ddeath[-index]
  active = active[-index]
  betas = drecover/active
  deltas = ddeath/active
  
  inp <- list(x_ini = country$x_ini, P = country$P, data = country$data, nrep = country$nrep,
              beta = betas, delta = deltas, travelin = country$travelin,
              travelout = country$travelout, total_out = country$total_out, pop = country$pop)
  
  
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
      
    )
  }
  #########Then prior density
  prior_eval <- function(theta){
    prior_value <-  dunif(theta["alpha"], 0, 2)*dunif(theta["beta"], 0, 1)*dunif(theta["delta"], 0, 1)*dunif(theta["gamma"], 0, 1)
    
    return(prior_value)
  }
  
  
  ##############################################################################
  ############3. Learn alpha based on available data, beta, and gamma##
  
  
  alpha_function = function(data,theta,inp){
    
    
    nr = nrow(data)
    beta = theta[3]
    gamma = theta[6]
    Ut = rowSums(data)
    #Difference Ut
    dUt = c(diff(Ut, differences = 1),0)
    ########
    recovermat = matrix(0,nrow = nr,ncol=6)
    
    infected = rep(0,nr)
    
    dRus = rep(0,nr)
    
    for (i in 1:nr){
      
      infected[i] = dUt[i]/gamma
      dRus[i] = dUt[i]/gamma*beta
      
    }
    
    
    Ruseq1 = cumsum(dRus)
    
    Ruseq = c(0,Ruseq1[-length(Ruseq1)])
    
    
    Sseq = inp$pop - Ut - Ruseq - infected
    
    
    SIseq = Sseq*infected
    
    #########Construct dSseq
    Sseq_p = rep(0, nr)
    
    for (i in 1:nr){
      Sseq_p[i] = Sseq[i] - inp$travelin[i,1] + inp$travelout[i,1]
    }
    ######S+ and Spre sequence for alpha
    Splus_al = Sseq[-c(nr-1, nr)]
    
    Spre_al = Sseq_p[-c(1,nr)]
    
    dSseq_al = Splus_al - Spre_al
    
    SIseq_al = SIseq[-c(nr-1, nr)]
    
    pop_al = inp$pop[-c(nr-1, nr)]
    
    
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
  
  
  #########Learn beta functions
  beta_function = function(data){
    drecover_sim = data[,2] - c(0,data[,2][-length(data[,2])])
    active_sim = c(1,data[,1][-length(data[,1])])
    
    index5 = which(active_sim==0) 
    index6 = c(1,index5)
    drecover_sim = drecover_sim[-index6]
    active_sim = active_sim[-index6]
    betas_sim = drecover_sim/active_sim
    ########This help to avoid empty array
    if (length(betas_sim) == 0){betas_sim = 10^8}
    ###########
    
    
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
    ########This help to avoid empty array
    if (length(deltas_sim) == 0){deltas_sim = 10^8}
    ###########
    return(deltas_sim)
    
  }
  
  
  
  ##################Recover the average path based on the observed data
  recover_function = function(inp,theta){
    
    data = inp$data
    pop = inp$pop
    
    beta = theta[3]
    gamma = theta[6]
    Ut = data[,1] + data[,2] + data[,3]
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
    
    
    Sseq = pop - Ut - Ruseq - infected
    
    recovermat = data.frame(Sseq, infected, data[,1], data[,2], data[,3], Ruseq)
    
    
    
    return(recovermat)
    
  }
  
  
  ##############
  stofunc_travel =  function(adjtheta,inp){
    
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
    status_matrix = matrix(0, nrow = d,ncol = 6)
    status_matrix[1,] = inp$x_ini
    
    f_in = matrix(0, nrow = d,ncol = 6)
    
    f_out = matrix(0, nrow = d, ncol = 6)
    
    for (i in 2:d){
      x = status_matrix[(i-1),]
      ##Updated traveling flow
      #Number out from country 2
      out = inp$total_out[i]
      if (x[1]+x[2] > 0){
        out_compartment = c(round(out*x[1]/(x[1]+x[2]),digits=0), round(out*x[2]/(x[1]+x[2]),digits=0), 0,0,0,0)
      } else{
        
        out_compartment = c(out, 0,0,0,0,0)
      }
      
      ##Generating Poisson values based on hazard functions
      y1 = rpois(1, harzard1(x,adjtheta))
      # 
      y2 =  rpois(1, harzard2(x,adjtheta))
      # 
      y3 = rpois(1, harzard3(x,adjtheta))
      # 
      y4 =  rpois(1, harzard4(x,adjtheta))
      # 
      y5 = rpois(1, harzard5(x,adjtheta))
      
      
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
      
      #update the flow traffic
      
      update = x  + inp$travelin[i,] - out_compartment
      
      update[update<0.1]=0
      
      status_matrix[i,] = update
      
      
    }
    return(status_matrix)
  }  
  
  
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
      u = stofunc_travel(theta,inp)
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
  
  
  inp = list.append(inp, sd=sd)
  
  
  
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
  ##1. L4 DISTANCE
  
  
  distance1 = function(theta,inp){
    
    sim <- stofunc_travel(theta,inp)
    ########## Learn beta, delta
    sim1 = sim[,3:5]
    
    betas_sim = beta_function(sim1)
    deltas_sim = delta_function(sim1)
    ###Learn alpha#######
    
    
    alphas = alpha_function(inp$data,theta,inp)
    alphas_sim = alpha_function(sim1,theta,inp)
    
    #############
    output1 <- abs(median(betas_sim) - median(inp$beta))
    output2 <-  abs(median(deltas_sim) - median(inp$delta))
    
    output3 <- abs(median(alphas_sim) - median(alphas))
    ############
    output4 <- sd_function1(inp$data[,1],sim[,3],inp)
    output5 <- sd_function2(inp$data[,2],sim[,4],inp)
    output6 <- sd_function3(inp$data[,3],sim[,5],inp)
    #############
    
    output <- (output1 + output2 + output3)^.5 + (1/nrow(inp$data)*(output4 + output5+ output6))^.5
    
    return(output)
  }
  
  ##
  
  abc_post_1 <- abc_start(
    prior,
    distance1, 
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = TRUE)
  )
  
  theta0 = country$thetatruth
  #####Collect output with descriptive stat
  ###############################
  coverage_func = function(data,theta){
    val11 = (quantile(data$alpha,.25) - theta[2])*(quantile(data$alpha,.75) - theta[2])
    val12 =  (quantile(data$beta,.25) - theta[3])*(quantile(data$beta,.75) - theta[3])
    val13 = (quantile(data$delta,.25) - theta[4])*(quantile(data$delta,.75) - theta[4])
    val14 = (quantile(data$gamma,.25) - theta[6])*(quantile(data$gamma,.75) - theta[6])
    #########
    val21 = (quantile(data$alpha,.025) - theta[2])*(quantile(data$alpha,.975) - theta[2])
    val22 =  (quantile(data$beta,.025) - theta[3])*(quantile(data$beta,.975) - theta[3])
    val23 = (quantile(data$delta,.025) - theta[4])*(quantile(data$delta,.975) - theta[4])
    val24 = (quantile(data$gamma,.025) - theta[6])*(quantile(data$gamma,.975) - theta[6])
    
    ####################
    if(val11 <0){ a1 = 1
    } else {a1 = 0 }
    
    if(val12 <0){ a2 = 1
    } else {a2 = 0 }
    
    if(val13 <0){ a3 = 1
    } else {a3 = 0 }
    
    if(val14 <0){ a4 = 1
    } else {a4 = 0 }
    ######
    if(val21 <0){ a5 = 1
    } else {a5 = 0 }
    
    if(val22 <0){ a6 = 1
    } else {a6 = 0 }
    
    if(val23 <0){ a7 = 1
    } else {a7 = 0 }
    
    if(val24 <0){ a8 = 1
    } else {a8 = 0 }
    
    return(c(a1, a2,a3,a4,a5,a6,a7,a8))
  }
  
  #####################Accuracy median
  ac_med_func = function(data,theta){
    val11 = (median(data$alpha) - theta[2])/theta[2]
    val12 =  (median(data$beta) - theta[3])/theta[3]
    val13 = (median(data$delta) - theta[4])/theta[4]
    val14 = (median(data$gamma) - theta[6])/theta[6]
    return(c(val11, val12, val13, val14))
  }
  
  
  ######Accuracy mean
  ac_mean_func = function(data,theta){
    val11 = (mean(data$alpha) - theta[2])/theta[2]
    val12 =  (mean(data$beta) - theta[3])/theta[3]
    val13 = (mean(data$delta) - theta[4])/theta[4]
    val14 = (mean(data$gamma) - theta[6])/theta[6]
    return(c(val11, val12, val13, val14))
  }
  
  
  
  IQR_func = function(data,theta){
    val11 = quantile(data$alpha,.75) - quantile(data$alpha,.25)
    val12 =  quantile(data$beta,.75) - quantile(data$beta,.25)
    val13 = quantile(data$delta,.75) - quantile(data$delta,.25)
    val14 = quantile(data$gamma,.75) - quantile(data$gamma,.25)
    
    #############
    val21 = quantile(data$alpha,.975) - quantile(data$alpha,.025)
    val22 =  quantile(data$beta,.975) - quantile(data$beta,.025)
    val23 = quantile(data$delta,.975) - quantile(data$delta,.025)
    val24 = quantile(data$gamma,.975) - quantile(data$gamma,.025)
    ##################
    
    
    return(c(val11,val12,val13,val14,val21,val22,val23,val24))
  }
  
  IQR1 = IQR_func(abc_post_1,theta0)
  theta_cover1 = coverage_func(abc_post_1, theta0)
  ac1as = c(ac_med_func(abc_post_1,theta0), ac_mean_func(abc_post_1,theta0))
  
  
  out1 = c( ac1as,theta_cover1,IQR1, theta0)
  
  
  distance2 = function(theta,inp){
    
    sim <- stofunc_travel(theta,inp)
    ########## Learn beta, delta
    sim1 = sim[,3:5]
    
    betas_sim = beta_function(sim1)
    deltas_sim = delta_function(sim1)
    
    #############
    output1 <- abs(median(betas_sim) - median(inp$beta))
    output2 <-  abs(median(deltas_sim) - median(inp$delta))
    
    
    ############
    output4 <- sd_function1(inp$data[,1],sim[,3],inp)
    output5 <- sd_function2(inp$data[,2],sim[,4],inp)
    output6 <- sd_function3(inp$data[,3],sim[,5],inp)
    #############
    
    output <- (output1 + output2 )^.5 + (1/nrow(inp$data)*(output4 + output5+ output6))^.5
    
    return(output)
  }
  
  ##
  
  abc_post_2 <- abc_start(
    prior,
    distance2, 
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = TRUE)
  )
  
  IQR2 = IQR_func(abc_post_2,theta0)
  theta_cover2 = coverage_func(abc_post_2, theta0)
  ac2as = c(ac_med_func(abc_post_2,theta0), ac_mean_func(abc_post_2,theta0))
  
  
  out2 = c( ac2as,theta_cover2,IQR2, theta0)
  
  #Euclidean
  ###3. EUCLIDEAN DISTANCE
  
  
  distance3 = function(theta,inp){
    
    sim <- stofunc_travel(theta,inp)
    
    r_active = sum((sim[,3] - inp$data[,1])^2)
    r_recover = sum((sim[,4] - inp$data[,2])^2)
    r_death = sum((sim[,5]-inp$data[,3])^2)
    
    
    output <-  (1/nrow(inp$data)*(r_active + r_recover+ r_death))^.5
    
    return(output)
  }
  
  ##
  
  abc_post_3 <- abc_start(
    prior,
    distance3, 
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = TRUE)
  )
  
  IQR3 = IQR_func(abc_post_3,theta0)
  theta_cover3 = coverage_func(abc_post_3, theta0)
  ac3as = c(ac_med_func(abc_post_3,theta0), ac_mean_func(abc_post_3,theta0))
  out3 = c(ac3as,theta_cover3,IQR3, theta0)
  
  #################
  #############################ADD NAIVE FOR A BACKUP ADVANTAGE#####
  
  #############naive visit############
  sto_generating_func = function(theta,inp){
    status_matrix1 = matrix(0,nrow = inp$nrep,ncol=6)
    status_matrix1[1,] = inp$x_ini
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
      
      status_matrix1[i,] = x
    }
    
    return(status_matrix1)
  }
  
  
  
  ##############################################################################
  ############3. Learn alpha based on available data, naive model
  
  
  alpha_function1 = function(data,theta,P){
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
  
  
  
  
  ##########################Standardized##########
  
  standardize_seq1 = function(k,inp){
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
  
  sd1 = standardize_seq1(3000,inp)[,-1]
  
  
  
  sd_function1 = function(vec1, vec2){
    
    output1 <- sum(((vec1[-1]-vec2[-1])/t(sd1[3,]))^2)
    return(output1)
    
  }
  sd_function2 = function(vec1, vec2){
    output1 <- sum(((vec1[-1]-vec2[-1])/t(sd1[4,]))^2)
    return(output1)
    
  }
  sd_function3 = function(vec1, vec2){
    output1 <- sum(((vec1[-1]-vec2[-1])/t(sd1[5,]))^2)
    return(output1)
    
  }
  
  
  
  #########
  ###4. EUCLIDEAN NAIVE
  
  
  distance4 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    
    r_active = sum((sim[,3] - inp$data[,1])^2)
    r_recover = sum((sim[,4] - inp$data[,2])^2)
    r_death = sum((sim[,5]-inp$data[,3])^2)
    
    
    output <-  (1/nrow(inp$data)*(r_active + r_recover+ r_death))^.5
    
    return(output)
  }
  
  ##
  
  abc_post_4 <- abc_start(
    prior,
    distance4, 
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = TRUE)
  )
  
  IQR4 = IQR_func(abc_post_4,theta0)
  theta_cover4 = coverage_func(abc_post_4, theta0)
  ac4as = c(ac_med_func(abc_post_4,theta0), ac_mean_func(abc_post_4,theta0))
  out4 = c(ac4as,theta_cover4,IQR4, theta0)
  ##########
  #######5. L4 Naive
  distance5 = function(theta,inp){
    
    sim <- sto_generating_func(theta,inp)
    output1 <- sd_function1(inp$data[,1],sim[,3])
    output2 <- sd_function2(inp$data[,2],sim[,4])
    output3 <- sd_function3(inp$data[,3],sim[,5])
    output4 <- 1/nrow(inp$data)*(output1 + output2 + output3)
    ################################
    mydat2 = sim[,3:5]
    betas_sim = beta_function(mydat2)
    deltas_sim = delta_function(mydat2)
    
    ##########Learn alpha#######
    
    
    alphas = alpha_function1(inp$data,theta,inp$P)
    alphas_sim = alpha_function1(mydat2,theta,inp$P)
    
    #############
    output5 <- abs(median(betas_sim) - median(inp$beta))
    output6 <-  abs(median(deltas_sim) - median(inp$delta))
    output7 <- abs(median(alphas_sim) - median(alphas))
    ############
    output = output4^.5 + (output5 + output6 + output7)^.5
    
    return(output)
  }
  #
  abc_post_5 <- abc_start(
    prior,
    distance5, 
    distance_args = inp,
    method = "RABC",
    control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
    output_control = list(print_output = TRUE)
  )
  
  IQR5 = IQR_func(abc_post_5,theta0)
  theta_cover5 = coverage_func(abc_post_5, theta0)
  ac5as = c(ac_med_func(abc_post_5,theta0), ac_mean_func(abc_post_5,theta0))
  out5 = c( ac5as,theta_cover5,IQR5, theta0)
  
  
  
  
  out = c(out1,out2,out3,out4,out5)
  return(out)
}
##########


#########Run ABC for 3 countries
out_1 = ABC_function(country1)
out_2 = ABC_function(country2)
out_3 = ABC_function(country3)
# #######Merge output
out = c(out_1,out_2,out_3)

fname =paste("Travelestimateoutput_particles",particles,"_pacc",prob,"_loop",k,".txt",sep="")
write.table(out,file=fname, row.names = F)

proc.time() - tmp

