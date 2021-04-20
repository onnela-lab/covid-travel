#This code is used to plot and save posterior with alpha varied every 4 weeks
#Date started: March 21
#Last edited: March 21, 4h50 p.m.

args <- commandArgs(trailingOnly = TRUE)
countryconsider <- as.integer(args[1])
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/coviddataJanJune20.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/flightdataJanJune20.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/AllPreliminariesposteriofromstep1changepoint.Rdata")
fname1 = paste('Allinitialfromstep1changepoint',".txt",sep="")
initialmat  = read.table(fname1,header=T)
initialmat = as.matrix(initialmat)


tmptmp = proc.time()

theconsidercountry = coviddataJanJune20$country[countryconsider]
#1. set up Acceptance rate and number of particles to run ABC
prob1 = 5 #change to 5
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

####Extract the change point
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






indicator = function(a){
  a[a<0]=0
  return(a)
}




###############2. Define uniform prior for parameters
##Step 2. ExTRACT POSTERIOR FOR STD###############2. Define uniform prior for parameters
prior <- function(n){
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
prior_eval <- function(theta){
  prior_value <-  dunif(theta["alpha0"], 0, 1)*dunif(theta["alpha"], 0, 1)*dunif(theta["beta"], 0, .25)*dunif(theta["delta"], 0, .25)*dunif(theta["gamma"], 0, 1)

  return(prior_value)
}





## Extract STD
####################################################
travelout_datadivided = flightdataJanJune20

#Standradize covid data and flight data with covid data in terms of time
mycompletedata = coviddataJanJune20 
temp = length(travelout_datadivided) - nrow(mycompletedata[[1]])
temp1 = 1:temp
travelout_datadivided = travelout_datadivided[-temp1]
###########Use the covid data and flight data up to a cutoff time we want
numbercountries = ncol(travelout_datadivided[[1]] )
numberactualcountries = numbercountries-1 ##number actual remove the fake 1 ZZZ
cutofftime = nrow(mycompletedata[[1]])   #
temp2 = 1:cutofftime
travelout_datadivided = travelout_datadivided[temp2]

for(i in 1:numbercountries){
  mycompletedata[[i]] = mycompletedata[[i]][temp2,]
}


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


##All start when it own first pass 500 cases
mydatacutoffadapt = list()

for ( i in 1:numberactualcountries){
  tmp = mycompletedata[[i]]
  startday1 = min(which( rowSums(tmp)>casecutoff))
  mydatacutoffadapt[[i]] = tmp[startday1:nrow(mycompletedata[[1]]),]
}

########Make all other country start at the same day at the considering country

countryconsiderdata1 = mycompletedata[[countryconsider ]]
startday = min(which(rowSums(countryconsiderdata1)>casecutoff))
countryconsiderdata = countryconsiderdata1[startday:nrow(mycompletedata[[1]]),]


mydatacutoff = list() # create a list of data for each country from 1 to n = 47 this case

for(i in 1:numbercountries){
  mydatacutoff[[i]] = mycompletedata[[i]][startday:nrow(mycompletedata[[1]]),]
}  ##all start as the same day as the considering country has more than 100 cases.

#2nd round extracting the corresponding flight data with covid data, cutoff based on the considering data first pass 100 cases
travelout_datadivided = travelout_datadivided[startday:nrow(mycompletedata[[1]])]


### ACHIEVED. ALL DATA AND TRAVEL DATA NOW ARE TRIM TO THE CONSIDERING COUNTRY

dataperiod = nrow(mydatacutoff[[countryconsider]])
testduration = 30
trainingtime = dataperiod - testduration # Training data up to May 31 
starttestingtime = trainingtime +1
mydatacutoff_train = list() # create a list of data for each country from 1 to n = 93 this case

for(i in 1:numbercountries){
  mydatacutoff_train[[i]] = mydatacutoff[[i]][1:trainingtime,]
  
}


mydatacutoffadapt_train = list()
for(i in 1:numberactualcountries){
  tmp = mydatacutoffadapt[[i]]
  h1 = nrow(tmp) - testduration
  mydatacutoffadapt_train[[i]] =  tmp[1:h1,]
}


mydatacutoff_test = list() # create a list of data for each country from 1 to n = 47 this case

for(i in 1:numbercountries){
  mydatacutoff_test[[i]] = mydatacutoff[[i]][starttestingtime:dataperiod,]
}


travelout_datadivided_train = travelout_datadivided[1:trainingtime]  #Number travelers from 1 country to another during the training

travelout_datadivided_test = travelout_datadivided[starttestingtime:dataperiod]
##Get initials from the independent model with kappa


##########Notice initial mat here corresponding to first time pass 500 cases of each country
############# Recover total travel out from each country
travelout_data = matrix(0,nrow=dataperiod,ncol=numbercountries)
for ( i in 1:dataperiod){
  travelout_data[i,] = rowSums(travelout_datadivided[[i]])
}

travelout_traindata =  travelout_data[1:trainingtime,]


#############################################################
############# TASK 2. Estimate unobserved infected each country ???

infectadapt =list() # list of number infected estimated from each country from 1 to n
for (i in 1:numberactualcountries){
  
  mydat = mydatacutoffadapt[[i]]
  mydat[,1]= mydat[,1]-mydat[,2] -mydat[,3] # restructure to get the shape required by infectfunc
  x = initialmat[i,]
  tmp = infectfunc(mydat,x)
  tmp1 = c(tmp,0)
  tmp1[tmp1<0] = 1 #make one for weird thing less than 0
  infectadapt[[i]] = round(tmp1,digits=0)
} #infect estimate for 46 countries with the adapted version


##Borrow infect adapt for infect regulated by the considering country
mydatamerge = list()
for(i in 1:numberactualcountries){
  
  datatmp1 = matrix(0,nrow=nrow(mydatacutoffadapt[[i]]),ncol=4)
  datatmp1[,2:4] = mydatacutoffadapt[[i]]
  datatmp1[,1] = infectadapt[[i]]
  
  datatmp2 = matrix(0,nrow=nrow(mydatacutoff[[i]]),ncol=4)
  datatmp2[,2:4]  = mydatacutoff[[i]]
  
  if(nrow(datatmp2)>nrow(datatmp1)){
    indextmp =c()
    for (j in 1:nrow(datatmp2)){
      tmp = datatmp2[j,2:4] - datatmp1[1,2:4]
      tmp1 = sum(tmp^2)
      if(tmp1==0){indextmp = c(indextmp,j)}
    }
    i1 = min(indextmp)
    datatmp2[i1:nrow(datatmp2),1] = datatmp1[,1]
  }else{
    i1 = nrow(datatmp1) - nrow(datatmp2)+1
    datatmp2[,1] = datatmp1[i1:nrow(datatmp1),1]
  }
  
  mydatamerge[[i]] = datatmp2
  
}



###############TASK 3. Reconstruct Total travel in and out at each time point, each column for one country
inmat = matrix(0, nrow = nrow(travelout_traindata), ncol = numbercountries)
outmat = matrix(0, nrow = nrow(travelout_traindata), ncol = numbercountries)

for (i in 1:length(travelout_datadivided_train)){
  for(j in 1:numbercountries){
    inmat[i,j] = sum(travelout_datadivided_train[[i]][,j]) #get in country j at day i
    outmat[i,j] = sum(travelout_datadivided_train[[i]][j,]) #get out country j at day i
  }
  
}
##############

countries = list() # create a list of inputs for country 1 to n
for (countryth in 1:numbercountries){
  if (countryth == numbercountries){
    countries[[countryth]] = list(P = P[countryth], total_in = inmat[,countryth], total_out = outmat[,countryth],
                                  infect = rep(0,nrow(mydatacutoff_train[[1]])),durationtravel = nrow(mydatacutoff_train[[1]]),
                                  data = mydatacutoff_train[[countryth]] , 
                                  betas1 = 0,
                                  deltas1 = 0)
  }else{
    tmp =  list(P = P[countryth], total_in = inmat[,countryth], total_out = outmat[,countryth],
                infect = mydatamerge[[countryth]][1:trainingtime,1],durationtravel = nrow(mydatacutoff_train[[1]]),
                data = mydatacutoff_train[[countryth]] , 
                betas1 = betafunc(mydatacutoffadapt_train[[countryth]]),
                deltas1 = deltafunc(mydatacutoffadapt_train[[countryth]]))
    countries[[countryth]] = tmp
  }
  
}

#################Remove weird betas and deltas candidates by late in report

for (countryth in 1:numbercountries){
  if(countryth == numbercountries){
    countries[[countryth]] = list.append(countries[[countryth]], betas = 0, deltas = 0)
  } else{
    betas = countries[[countryth]]$betas1
    deltas = countries[[countryth]]$deltas1
    betas = betas[betas>-.001]
    deltas = deltas[deltas>-.001] #remove negative numbers
    countries[[countryth]] = list.append(countries[[countryth]], betas = betas, deltas = deltas)
  }
  
}
##################Implement betas and deltas if they are less than 0 by data entry miistakes

##########TASK 4. append population dynamic for each country
for (countryth in 1:numbercountries){
  tmp = countries[[countryth]]
  tmp1 = populationdynamicfunc(tmp)
  countries[[countryth]] = list.append(countries[[countryth]], populationdynamic = tmp1)
}


########################TASK 5. travel in and out compartments

travelout_compartments =  list()  # create a travel in compartment list from other countries to country 1, ...,n

for (i in 1:numbercountries){
  
  
  tmp = rep(0, length(travelout_datadivided_train) )#get number travel from j to i sequence
  for (k in 1:length(travelout_datadivided_train)){
    tmp[k] = sum(travelout_datadivided_train[[k]][i,])
  }
  
  
  ##########
  if(i == numbercountries){
    ZZZ = matrix(0,trainingtime,6)
    ZZZ[,1] = tmp
    tmp1 = ZZZ
  }else{
    tmp1 = traveloutcompfunc(tmp, countries[[i]]) 
  }
  
  
  travelout_compartments[[i]] = round(tmp1,digits=0) # compartments travel  from country i during the travelduration
  
}

#########################################
travelineachday = matrix(0,length(travelout_datadivided_train), numbercountries*6) #storage matrix of travel in
###each row is 1 day, col (i-1)*6+1, i*6: for each coutry
#get number travel from j to i sequence
for (k in 1:length(travelout_datadivided_train)){
  
  
  ##Distribute undetected infected to different countries at day k
  f_outmat = matrix(0,numbercountries,numbercountries*6 )
  for (countryout in 1:numbercountries) {
    
    
    
    d1 = (countryout  - 1) * 6 + 1
    d2 = countryout  * 6
    f_outtotal = travelout_compartments[[countryout]][k,] ##compartments go out day k
    infect_outtotal = f_outtotal[2]
    
    travelouttmp = travelout_datadivided_train[[k]]
    
    ##get weights based on travel data
    if (sum(travelouttmp[countryout, ]) > 0){
      probdistribute = travelouttmp[countryout,]/sum(travelouttmp[countryout, ])
    }else{
      probdistribute = rep(0, numbercountries)
    }
    
    ##distribute undetected infected cases
    
    if (sum(travelouttmp[countryout, ]) > 0) {
      infect_outdistribute = rmultinom(1, size = infect_outtotal, 
                                       prob = probdistribute)
    }else{
      infect_outdistribute = rep(0, numbercountries)
    }
    
    ##Distribute undocumented infected cases for countries
    for (countryin in 1:numbercountries) {
      e1 = (countryin - 1) * 6 + 1
      e2 = countryin * 6
      if (sum(travelouttmp[countryout, ]) > 0) {
        tmp = round(f_outtotal * travelouttmp[countryout, 
                                              countryin]/sum(travelouttmp[countryout, ]), digits = 0)
        suseptible = tmp[1] + tmp[2] - infect_outdistribute[countryin]
        f_outmat[countryout, e1:e2] = c(suseptible, infect_outdistribute[countryin], 
                                        tmp[3], tmp[4], tmp[5], tmp[6])
      }else{
        f_outmat[countryout, e1:e2] = rep(0, 6)
      }
    }  ###end countryin loop
  }  ###end countryout loop
  
  travelineachday[k,] = colSums(f_outmat) 
  
  
}  ###end k for each day loop

travelin_compartments =  list()  # create a travel in compartment list from other countries to country 1, ...,n

for (i in 1:numbercountries){
  a1 = (i-1)*6 +1
  a2 =i*6
  
  travelin_compartments[[i]] = travelineachday[,a1:a2]
}

###############################
for (i in 1:numbercountries){
  countries[[i]] = list.append(countries[[i]], x_ini = initialmat[i,],
                               travelin_compartments = travelin_compartments[[i]], 
                               travelout_compartments = travelout_compartments[[i]] )
}



#######Part II. #Estimate the parameters of the considering country.

country = countries[[countryconsider]]
country$x_ini = initialmat[countryconsider,] # get the initial condition as changepoint 
country$changepoint = indexchange

prelimposterior =  as.matrix(JHKposterior[[countryconsider]]) #Borrow preliminaries particles to learn the std
########Define the marginal function

stochastic_marginalestimate_alphachange =function (theta, country){
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
  ##
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
  status_matrix = matrix(0, nrow = country$durationtravel,
                         ncol = 6)
  status_matrix[1, ] = country$x_ini
  f_in = matrix(0, nrow = country$durationtravel, ncol = 6)
  f_out = matrix(0, nrow = country$durationtravel, ncol = 6)
  for (i in 2:country$durationtravel) {
    x = status_matrix[(i - 1), ]
    out = country$total_out[i]
    if (x[1] + x[2] > 0) {
      out_compartments = c(round(out * x[1]/(x[1] + x[2]),
                                 digits = 0), round(out * x[2]/(x[1] + x[2]),
                                                    digits = 0), 0, 0, 0, 0)
    }
    else {
      out_compartments = c(out, 0, 0, 0, 0, 0)
    }
    if(i<country$changepoint){
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
   update = x + country$travelin_compartments[i, ] - out_compartments
    update[update < 0.1] = 0
    status_matrix[i, ] = update
  }
  return(status_matrix)
}








###########5. Define the standardize function to get sd of A, R, D
thresholddeath = .05
standardize_seq = function(k,country){
  
  sd = rep(0, country$durationtravel) 
  
  
  Umat = matrix(0,k,country$durationtravel)
  Dmat = matrix(0,k,country$durationtravel)
  
  
  i <- 1
  while (i <= k ) {
    
    index= sample(1:nrow(prelimposterior) ,1)
    theta = prelimposterior[index,]
    h1 = as.numeric(theta["alpha0"]/2)
    theta["alpha0"] = theta["alpha0"]+runif(1,-h1,h1)
    h2 = as.numeric(theta["alpha"]/2)
    theta["alpha"] = theta["alpha"]+runif(1,-h2,h2)
        
    
    
    
    
    u = stochastic_marginalestimate_alphachange(theta,country)
    #########making sure data generated not too far away from the truth
    checkval = abs(u[nrow(u),5] - country$data[nrow(u),3])/country$data[nrow(u),3]    
    checkval1 = abs(sum(u[nrow(u),3:5]) - sum(country$data[nrow(u),]))/sum(country$data[nrow(u),])
    
    CHECK = max(checkval, checkval1)
    
    if (CHECK <thresholddeath ){  
      
      Umat[i,] = c(0,diff(rowSums(u[,3:5]),1))  #use daily to avoid accumulative error
      Dmat[i,] = c(0,diff(u[,5],1)) #use daily to avoid accumulative error
      
      i = i+1
    }
    
  }
  
  
  U_sd = apply(Umat,2,sd)
  
  D_sd = apply(Dmat,2,sd)
  
  return(rbind(U_sd, D_sd))
  
}





sd = standardize_seq(trials,country)[,-1]


country = list.append(country, sd=sd)



sdU_function = function(vec1, vec2, country){
  sdU = country$sd[1,]
  output <- sum(((vec1-vec2)/t(sdU))^2)
  return(output)
}

sdD_function = function(vec1, vec2, country){
  sdD = country$sd[2,]
  output <- sum(((vec1-vec2)/t(sdD))^2)
  return(output)
  
}
####################Define the distance


distance = function(theta, country){
  
  
  sim <- stochastic_marginalestimate_alphachange(theta, country)
  ################Learn beta, delta########################
  U = rowSums(sim[,3:5])
  U1 = diff(U,1) #total cases
  Ur = diff(rowSums(country$data),1)
  Ur = indicator(Ur)
  Ur = ave7func(Ur) 
  D1 = diff(sim[,5],1) #deaths
  Dr = diff(country$data[,3],1)
  Dr = indicator(Dr)
  Dr = ave7func(Dr)
  ##Standardized version
  output1 <- sdU_function(Ur,U1,country)
  output2 <- sdD_function( Dr, D1,country)
  
  ###EU distances
  # output1 <- sum((U1 - Ur)^2) # Total Confirmed daily
  # output2 <- sum((D1 - Dr)^2) #Total Death daily
  ####################
  output <- 1/nrow(country$data)*(output1^.5 + output2^.5)
  
  return(output)
}
#############

#####6. Run ABC
abcposterior <- abc_start(
  prior,
  distance,
  distance_args = country,
  method = "RABC",
  control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
  output_control = list(print_output = TRUE)
)
print("Tenten")
abcposterior = as.matrix(abcposterior)
fname = paste('Posteriorproposedwithmultiplealphachangepointcountryconsider_',countryconsider,".txt",sep="")
write.table(abcposterior, file=fname, row.names = F)



############DONE PART II, THE ESTIMATION STEP##########
#########PART III. MATCHING#############################
######################################################################


mylistall = list()
realizations = 100 #number realizations
averageall = rep(0,nrow(abcposterior))
for(particle in 1:nrow(abcposterior)){
  
  mylistall[[particle]] = list()
  average = 0
  theta = abcposterior[particle,]
  for (i in 1:realizations){
    
    data = stochastic_marginalestimate_alphachange(theta, country)
    tmp = sum((rowSums(data[,3:5]) - rowSums(country$data))^2) + sum((data[,5] - country$data[,3])^2)
    
    average = average + tmp
    mylistall[[particle]][[i]] = data
  }
  averageall[particle] = average/realizations
  
} #end loop for paricle

indexbest = min(which.min(averageall)) # choose any smallest in case there are 2 or more with minimum
########################
mylist = mylistall[[indexbest]]

############Matching for daily active confirmed
trainingtime = nrow(mylist[[1]])
activemat =  matrix(0,nrow=trainingtime, ncol=realizations)

for (i in 1:realizations){
  activemat[,i] = c(0,diff(rowSums(mylist[[i]][,3:5]),1))
  #activemat[,i] = rowSums(mylist[[i]][,3:5])
  
}

###############
scalefactor = 1


activemat = activemat/scalefactor # scale to to smaller number
quantileactive = rowQuantiles( activemat, probs =c(.025,.5,.975))
############# real data
confirmed =  c(0,diff(rowSums(country$data)/scalefactor,1))
#confirmed =  inp$data[,1]/scalefactor
confirmedupper = max(confirmed)*1.5
#########create data set
confirmeddata = cbind(quantileactive, confirmed)


colnames(confirmeddata) = c("lower","median", "upper", "datareport")
confirmeddata1 = log(confirmeddata)/log(10)

#############add time

t = 1:trainingtime
wholeperiod = rev(seq(as.Date("2020-06-30"), length =  dataperiod, by = "-1 day"))
trainingperiod = wholeperiod[1:trainingtime]

######Create the data
dataconfirmed = data.frame( trainingperiod,confirmeddata)
dataconfirmed = dataconfirmed[-1,]

dataconfirmed1 = data.frame(trainingperiod,confirmeddata1)
dataconfirmed1 = dataconfirmed1[-1,]

##Plot real data vs estimation ##
plt1a = ggplot(data=dataconfirmed, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+ theme_minimal()+
  geom_line(aes(y=median), colour="blue") +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs(title = theconsidercountry, x = "", y = "Daily confirmed")+ylim(0,confirmedupper)


#################################
############Matching for daily death confirmed

deathmat =  matrix(0,nrow=trainingtime,ncol=realizations)

for (i in 1:realizations){
  deathmat[,i] = c(0, diff(mylist[[i]][,5],1))
  #deathmat[,i] =  mylist[[i]][,5]
}


quantiledeath = rowQuantiles( deathmat, probs =c(.025,.5,.975))
############# real data
death =  c(0, diff(country$data[,3]),1)
#death =  inp$data[,3]
deathupper = max(death)*1.5
#########create data set
deathdata = cbind(quantiledeath, death)

colnames(deathdata) = c("lower","median", "upper", "datareport")
deathdata1 = log(deathdata)/log(10)
######Create the data
datadeath = data.frame( trainingperiod,deathdata)
datadeath  = datadeath[-1,]
datadeath1 = data.frame( trainingperiod,deathdata1)
datadeath1 = datadeath1[-1,]
##Plot real data vs estimation ##
plt1b = ggplot(data=datadeath, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+ theme_minimal()+
  geom_line(aes(y=median), colour="blue") + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs( x = "", y = "Daily death")+ylim(0,deathupper)




fname1a = paste('Dailyconfirmedcasesfitmultiplealphachangepointcountryconsider_',countryconsider,".png",sep="")
ggsave(fname1a,plot=plt1a)

fname1b = paste('Dailyconfirmeddeathfitmultiplealphachangepointcountryconsider_',countryconsider,".png",sep="")
ggsave(fname1b,plot=plt1b)


###################MATCHING WITH ACCUMULATED
############Matching for daily active confirmed
trainingtime = nrow(mylist[[1]])
activemat =  matrix(0,nrow=trainingtime, ncol=realizations)

for (i in 1:realizations){
  
  activemat[,i] = rowSums(mylist[[i]][,3:5])
  
}

###############
scalefactor = 1


activemat = activemat/scalefactor # scale to to smaller number
quantileactive = rowQuantiles( activemat, probs =c(.025,.5,.975))
############# real data

confirmed =  rowSums(country$data)/scalefactor
confirmedupper = max(confirmed)*1.5
#########create data set
confirmeddata = cbind(quantileactive, confirmed)


colnames(confirmeddata) = c("lower","median", "upper", "datareport")
confirmeddata1 = log(confirmeddata)/log(10)

#############add time

t = 1:trainingtime
wholeperiod = rev(seq(as.Date("2020-06-30"), length =  dataperiod, by = "-1 day"))
trainingperiod = wholeperiod[1:trainingtime]

######Create the data
dataconfirmed = data.frame( trainingperiod,confirmeddata)

dataconfirmed1 = data.frame(trainingperiod,confirmeddata1)

##Plot real data vs estimation ##
plt2a = ggplot(data=dataconfirmed, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+ theme_minimal()+
  geom_line(aes(y=median), colour="blue") +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs(title = theconsidercountry, x = "", y = "Accumulated confirmed")+ylim(0,confirmedupper)


#################################
############Matching for daily death confirmed

deathmat =  matrix(0,nrow=trainingtime,ncol=realizations)

for (i in 1:realizations){
 
  deathmat[,i] =  mylist[[i]][,5]
}


quantiledeath = rowQuantiles( deathmat, probs =c(.025,.5,.975))
############# real data

death =  country$data[,3]
deathupper = max(death)*1.5
#########create data set
deathdata = cbind(quantiledeath, death)

colnames(deathdata) = c("lower","median", "upper", "datareport")
deathdata1 = log(deathdata)/log(10)
######Create the data
datadeath = data.frame( trainingperiod,deathdata)
datadeath1 = data.frame( trainingperiod,deathdata1)

##Plot real data vs estimation ##
plt2b = ggplot(data=datadeath, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+ theme_minimal()+
  geom_line(aes(y=median), colour="blue") + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs( x = "", y = "Accumulated death")+ylim(0,deathupper)


plt2 = plot_grid( plt2a, plt2b,ncol=2)
fname2a = paste('Accumulatedconfirmedcasesfitmultiplealphachangepointcountryconsider_',countryconsider,".png",sep="")
ggsave(fname2a,plot=plt2a)


fname2b = paste('Accumulatedconfirmeddeathfitmultiplealphachangepointcountryconsider_',countryconsider,".png",sep="")
ggsave(fname2b,plot=plt2b)
##############################
###Save one best initial and parameters for the prediction step
alldistances = rep(0, realizations)

for (i in 1:realizations){

  data = as.matrix(mylist[[i]])
  tmp = sum((rowSums(data[,3:5]) - rowSums(country$data))^2) + sum((data[,5] - country$data[,3])^2)
  alldistances[i] = tmp
}

indexbest1 = min(which.min(alldistances)) # choose any smallest in case there are 2 or more minimum values

newinitial = mylist[[indexbest1]][nrow(mylist[[1]]),]
newinitial[3:5] = country$data[nrow(country$data),] # replace A,R,D by real data for a better fit
thetahat1 = abcposterior[indexbest,]
indexchange1 = indexchange +14
alphaindex =c(1,2)
if(indexchange1<trainingtime){
thetahat = c(0,thetahat1[-1])
}else{
thetahat = c(0, thetahat1[1],thetahat1[-alphaindex])
}



tmp = rbind(newinitial, thetahat)
fname3 = paste('Bestinitialandparametersforpredictionchangepointcountryconsider_',countryconsider,".txt",sep="")
write.table(tmp,file=fname3, row.names = F)







proc.time() - tmptmp



