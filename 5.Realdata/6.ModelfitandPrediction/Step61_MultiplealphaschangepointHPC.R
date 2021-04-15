#This code is used to extract the initial condition and posterior each country independently with alpha varied change point detection
#Date started: March 28
#Last edited: March 28, 12h52 p.m.

args <- commandArgs(trailingOnly = TRUE)
countryconsider <- as.integer(args[1])
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/coviddataJanJune20.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/flightdataJanJune20.Rdata")

tmptmp = proc.time()

theconsidercountry = coviddataJanJune20$country[countryconsider]
#1. set up Acceptance rate and number of particles to run ABC
prob1 = 10 #need to change to 10
prob =prob1/1000
particles = 1000 # Need to change to 1000

trials = 1000 # change to 1000


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
library(parallel)

##########################################################
############STEP1:LEARN I0####################
mycompletedata = coviddataJanJune20

P = mycompletedata$population[countryconsider] #population considering country
########Choose the start day to process data for the country when number of cases large enough
countryconsiderdata = mycompletedata[[countryconsider ]]


#CHOOSE CUTOFF AS 500 CONFIRMED CASES FOR THE START DAY WITH JHK Original
temp = countryconsiderdata[,1] ## CAN CHOOSE OTHER COLUMN IF WANT WORLDOMETERS AS THE CUTOFF
casecutoff = 500
startday = min(which(temp>casecutoff))
countryconsiderdata = countryconsiderdata[startday:nrow(countryconsiderdata),1:3] #USING JHK
dataperiod = nrow(countryconsiderdata)


############## Use the first d1 days as training

trainingtime = dataperiod - 30 ## Use the first d1 days as training


countryconsiderdata_train = countryconsiderdata[1:trainingtime,]
########


ave7func = function(a1){
  a11 = length(a1)-6
  for(i in 1:a11){
    i1 = i+1
    i2 = i+2 
    i3 = i+3
    i4 = i+4
    i5 = i+5 
    i6 = i+6
    a1[i] = (a1[i]+a1[i1]+ a1[i2] +a1[i3]+a1[i4]+a1[i5]+ a1[i6])/7
  }
  a12 = a11+1
  a13 =a11+2
  a14 =a11+3
  a15 = a11+4
  a16 =a11+5
  a17 =a11+6
  a1[a12] = (a1[a12]+ a1[a13]+a1[a14]+ a1[a15]+ a1[a16]+a1[a17])/6
  a1[a13] = ( a1[a13]+a1[a14]+ a1[a15]+ a1[a16]+a1[a17])/5
  a1[a14] = ( a1[a14]+ a1[a15]+ a1[a16]+a1[a17])/4
  a1[a15] = (  a1[a15]+ a1[a16]+a1[a17])/3
  a1[a16] = ( a1[a16]+a1[a17])/2
  
  return(a1)
}

##Extract the change point
data = countryconsiderdata_train
confirmedreport = ave7func(data[,1])
lag1 = c(0,diff(confirmedreport,1))
lag2 = c(0,diff(lag1,1))
index1 = which(lag2<0)

if(length(index1)>1){
  indexchange = max(index1[2],14) #use 2 to guard against random drop if happened and constraint of 14 to avoid too short period
}else{
  indexchange=trainingtime
}


#####Get the initial value for I0 from capture recapture, same start day with the US 

initialmat = rep(0,6)
I0 = countryconsiderdata_train[1,1] #use for a temporary value
A0 = countryconsiderdata_train[1,1] - countryconsiderdata_train[1,2] - countryconsiderdata_train[1,3]
R0 = countryconsiderdata_train[1,2] 
D0 = countryconsiderdata_train[1,3] 


initialmat = c((P - I0 - A0 - R0 - D0), I0, A0, R0, D0, 0 )

inp = list(duration= nrow(countryconsiderdata_train), ini = initialmat,
           data = countryconsiderdata_train, P = P, changepoint = indexchange)



######################
prior1 <- function(n){
  data.frame(
    alpha0 = runif(n,0,1),
    alpha = runif(n,0,1),
    beta = runif(n,0,.25),
    delta=runif(n,0,.25),
    eta=1,
    gamma=runif(n,0,1),
    kappa = runif(n,0,50)

  )
}
#########3. Then prior density
prior_eval1 <- function(theta){
  prior_value <-  dunif(theta["alpha0"], 0, 1)*dunif(theta["alpha"], 0, 1)*dunif(theta["beta"], 0, .25)*dunif(theta["delta"], 0, .25)*dunif(theta["gamma"], 0, 1)*dunif(theta["kappa"], 0, 50)

  return(prior_value)
}

#################Define function with kappa

stochasticmodel_1countrykappa_alphachange = function (theta, inp)
{
  status_matrix = matrix(0, nrow = inp$duration, ncol = 6)
  inp$ini[2] = sum(inp$ini[3:5])*theta[7]
  status_matrix[1, ] = inp$ini
  ##This is alpha at the 1st period
  harzard1a = function(x, theta) {
    h1a = theta[1] * x[1] * x[2]/sum(x)
    names(h1a) = c("hazard1a")
    return(h1a)
  }
  ##This is alpha at the 2nd period
  harzard1b = function(x, theta) {
    h1b = theta[2] * x[1] * x[2]/sum(x)
    names(h1b) = c("hazard1b")
    return(h1b)
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
  for (i in 2:inp$duration) {
    x = status_matrix[(i - 1), ]
    if(i<inp$changepoint){
      y1 = rpois(1, harzard1a(x, theta))
    } else{
      y1 = rpois(1, harzard1b(x, theta))
    }

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
    status_matrix[i, ] = x
  }
  return(round(status_matrix, digits = 0))
}



############use the indicator function to remove negative jumps
indicator = function(a){
  a[a<0]=0
  return(a)
}




########BASED ON PARAMETERS##########



######PARALLEL TO ACCELERATE
distance1 = function(theta, inp){
  sim <-  stochasticmodel_1countrykappa_alphachange(theta,inp)
  
  # #############REPLACE BY EUCLIDEAN
  U = rowSums(sim[,3:5])
  U1 = diff(U,1)
  Ur = diff(inp$data[,1],1)
  Ur = indicator(Ur)
  Ur = ave7func(Ur)
  D1 = diff(sim[,5],1)
  Dr = diff(inp$data[,3],1)
  Dr = indicator(Dr)
  Dr =  ave7func(Dr)  
  
  
  
  output1 <- sum((U1 - Ur)^2) # Total Confirmed daily
  output2 <- sum((D1 - Dr)^2) #Total Death daily
  
  output <- 1/nrow(inp$data)*(output1^.5 + output2^.5)
  
  return(output)
}

###########
##################
abcposterior1 <- abc_start(
  prior1,
  distance1,
  distance_args = inp,
  method = "RABC",
  control = list(prior_eval = prior_eval1,  n = particles, pacc_final = prob),
  output_control = list(print_output = TRUE)
)


########7. Estimate the parameters###################

thetahat1 = apply(abcposterior1, 2, median)
inp$ini[2] = sum(inp$ini[3:5])*thetahat1["kappa"]

##############################
##DONE STEP1, EXTRACT I0#############
####STEP2. ExTRACT POSTERIOR FOR STD###############2. Define uniform prior for parameters
prior2 <- function(n){
  data.frame(
    alpha0 = runif(n,0,1),
    alpha = runif(n,0,1),
    beta = runif(n,0,.25),
    delta=runif(n,0,.25),
    eta=1,
    gamma=runif(n,0,1)

  )
}
#########3. Then prior density
prior_eval2 <- function(theta){
  prior_value <-  dunif(theta["alpha0"], 0, 1)*dunif(theta["alpha"], 0, 1)*dunif(theta["beta"], 0, .25)*dunif(theta["delta"], 0, .25)*dunif(theta["gamma"], 0, 1)

  return(prior_value)
}
###############################################

stochasticmodel_1country_alphachange = function (theta, inp)
{
  status_matrix = matrix(0, nrow = inp$duration, ncol = 6)
  status_matrix[1, ] = inp$ini
  ##This is alpha at the 1st period
  harzard1a = function(x, theta) {
    h1a = theta[1] * x[1] * x[2]/sum(x)
    names(h1a) = c("hazard1a")
    return(h1a)
  }
  ##This is alpha at the 2nd period
  harzard1b = function(x, theta) {
    h1b = theta[2] * x[1] * x[2]/sum(x)
    names(h1b) = c("hazard1b")
    return(h1b)
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
  for (i in 2:inp$duration) {
    x = status_matrix[(i - 1), ]
    if(i<inp$changepoint){
      y1 = rpois(1, harzard1a(x, theta))
    } else{
      y1 = rpois(1, harzard1b(x, theta))
    }

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
    status_matrix[i, ] = x
  }
  return(round(status_matrix, digits = 0))
}







#############
distance2 = function(theta, inp){
  
  sim <- stochasticmodel_1country_alphachange(theta,inp)
  # #############REPLACE BY EUCLIDEAN
  U = rowSums(sim[,3:5])
  U1 = diff(U,1)
  Ur = diff(inp$data[,1],1)
  Ur = indicator(Ur)
  Ur = ave7func(Ur) 
  D1 = diff(sim[,5],1)
  Dr = diff(inp$data[,3],1)
  Dr = indicator(Dr)
  Dr = ave7func(Dr)
  output1 <- sum((U1 - Ur)^2) # Total Confirmed daily
  output2 <- sum((D1 - Dr)^2) #Total Death daily
  
  output <- 1/nrow(inp$data)*(output1^.5 + output2^.5)
  
  return(output)
}



#############################


abcposterior2 <- abc_start(
  prior2,
  distance2,
  distance_args = inp,
  method = "RABC",
  control = list(prior_eval = prior_eval2,  n = particles, pacc_final = prob),
  output_control = list(print_output = TRUE)
)

#######################################
abcposterior2 = as.matrix(abcposterior2)


#######DONE WITH STEP 2 GET POSTERIOR TO TO EXTRACT STD

output = list(initial =  inp$ini, posteriorstep1 = abcposterior2)

fname = paste('Posteriorandinitialfromstep1withmultiplealphachangepointcountry',countryconsider,".Rdata",sep="")
save(output, file=fname)

proc.time() - tmptmp



