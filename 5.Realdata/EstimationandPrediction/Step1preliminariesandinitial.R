##This code is used to get the initial conditions and initial ABC
#Date started: Feb 24
#Last edited: March 3, 15h27p.m.

args <- commandArgs(trailingOnly = TRUE)
countryconsider <- as.integer(args[1])



#1. set up Acceptance rate and number of particles to run ABC
prob1 = 5
prob =prob1/1000
particles = 1000 # Need to change to 1000




library(protoABC)
library(CONETTravel)
library(matrixStats)
library(matrixcalc)
library(cowplot)
library(ggplot2)
tmp = proc.time()

load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/coviddataJanJune20.Rdata")
##########################################################
###############ROUND 1: JHK data################
mycompletedata = coviddataJanJune20

P = mycompletedata$population[countryconsider] #population considering country
########Choose the start day to process data for the country when number of cases large enough
countryconsiderdata = mycompletedata[[countryconsider ]]


#CHOOSE CUTOFF AS 500 CONFIRMED CASES FOR THE START DAY WITH JHK Original
temp = countryconsiderdata[,1] ## CAN CHOOSE OTHER COLUMN IF WANT WORLDOMETERS AS THE CUTOFF
casecutoff = 1000
startday = min(which(temp>casecutoff))
countryconsiderdata = countryconsiderdata[startday:nrow(countryconsiderdata),1:3] #USING JHK
dataperiod = nrow(countryconsiderdata)


############## Use the first d1 days as training

trainingtime = dataperiod - 30 ## Use the first d1 days as training


countryconsiderdata_train = countryconsiderdata[1:trainingtime,]



#####Get the initial value for I0 from capture recapture, same start day with the US 

initialmat = rep(0,6)
I0 = countryconsiderdata_train[1,1] #use for a temporary value
A0 = countryconsiderdata_train[1,1] - countryconsiderdata_train[1,2] - countryconsiderdata_train[1,3]
R0 = countryconsiderdata_train[1,2] 
D0 = countryconsiderdata_train[1,3] 


initialmat = c((P - I0 - A0 - R0 - D0), I0, A0, R0, D0, 0 )

inp = list(duration= nrow(countryconsiderdata_train), ini = initialmat, data = countryconsiderdata_train, P = P)



######################
prior1 <- function(n){
  data.frame(
    alpha0 = 0,
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
  prior_value <-  dunif(theta["alpha"], 0, 1)*dunif(theta["beta"], 0, .25)*dunif(theta["delta"], 0, .25)*dunif(theta["gamma"], 0, 1)*dunif(theta["kappa"], 0, 50)
  
  return(prior_value)
}


#################Define function with kappa
stochasticmodel_1countrykappa = function (theta, inp) 
{
  status_matrix = matrix(0, nrow = inp$duration, ncol = 6)
  inp$ini[2] = sum(inp$ini[3:5])*theta[7]
  status_matrix[1, ] = inp$ini
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
  for (i in 2:inp$duration) {
    x = status_matrix[(i - 1), ]
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
    status_matrix[i, ] = x
  }
  return(round(status_matrix, digits = 0))
}

##########Define a function to smooth data due to delay report of lag 4

ave4func = function(a1)
{
  
a11 = length(a1)-3
for(i in 1:a11){
  i1 = i+1
  i2 = i+2 
  i3 = i+3
  a1[i] = (a1[i]+a1[i1]+ a1[i2] +a1[i3])/4
}
a12 = a11+1
a13 =a11+2
a14 =a11+3
a1[a12] = (a1[a12]+ a1[a13]+a1[a14])/3
a1[a13] = ( a1[a13]+a1[a14])/2
return(a1)
}


############use the indicator function to remove negative jumps
indicator = function(a){
  a[a<0]=0
  return(a)
}






######PARALLEL TO ACCELERATE
distance1 = function(theta, inp){
  sim <-  stochasticmodel_1countrykappa(theta,inp)
  
  # #############REPLACE BY EUCLIDEAN
  U = rowSums(sim[,3:5])
  U1 = diff(U,1)
  Ur = diff(inp$data[,1],1)
  Ur = indicator(Ur)
  Ur = ave4func(Ur)
  D1 = diff(sim[,5],1)
  Dr = diff(inp$data[,3],1)
  Dr = indicator(Dr)
  Dr =  ave4func(Dr)  



  output1 <- sum((U1 - Ur)^2) # Total Confirmed daily
  output2 <- sum((D1 - Dr)^2) #Total Death daily
  
  output <- 1/nrow(inp$data)*(output1^.5 + output2^.5)
  
  return(output)
}

###########
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
inp$ini[2] = sum(inp$ini[3:5])*thetahat1[7]

##############################

####STEP2. GO WITH OUR MODEL###############2. Define uniform prior for parameters
prior <- function(n){
  data.frame(
    alpha0 = 0,
    alpha = runif(n,0,1),
    beta = runif(n,0,.25),
    delta=runif(n,0,.25),
    eta=1,
    gamma=runif(n,0,1)
    
  )
}
#########3. Then prior density
prior_eval <- function(theta){
  prior_value <-  dunif(theta["alpha"], 0, 1)*dunif(theta["beta"], 0, .25)*dunif(theta["delta"], 0, .25)*dunif(theta["gamma"], 0, 1)
  
  return(prior_value)
}


#############
distance = function(theta, inp){
  
  sim <- CONETTravel::stochasticmodel_1country(theta,inp)
  # #############REPLACE BY EUCLIDEAN
  U = rowSums(sim[,3:5])
  U1 = diff(U,1)
  Ur = diff(inp$data[,1],1)
  Ur = indicator(Ur)
  Ur = ave4func(Ur) 
  D1 = diff(sim[,5],1)
  Dr = diff(inp$data[,3],1)
  Dr = indicator(Dr)
  Dr = ave4func(Dr)
  output1 <- sum((U1 - Ur)^2) # Total Confirmed daily
  output2 <- sum((D1 - Dr)^2) #Total Death daily
  
  output <- 1/nrow(inp$data)*(output1^.5 + output2^.5)
  
  return(output)
}


abcposterior <- abc_start(
  prior,
  distance,
  distance_args = inp,
  method = "RABC",
  control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
  output_control = list(print_output = TRUE)
)

#######################################
abcposterior = as.matrix(abcposterior)
posteriorandinitialmat = rbind(inp$ini, abcposterior)
fname = paste('C1000InitialandPosteriorprelim_JHKdata_countryconsider_',countryconsider,".txt",sep="")
write.table(posteriorandinitialmat,file=fname, row.names = F)

