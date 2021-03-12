###This code is used to get the initial conditions and initial for prediction with the standardized distance
#Date started: Feb 26
#Last edited: March 11, 16h05 a.m.

args <- commandArgs(trailingOnly = TRUE)
countryconsider <- as.integer(args[1])



prob =10/1000
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

tmp = proc.time()




load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/coviddataJanJune20.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/flightdataJanJune20.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/AllPreliminariesposteriorandinitialStep1.Rdata")



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

#################
casecutoff = 500 #CHOOSE CUTOFF AS 500 CONFIRMED CASES FOR THE START DAY

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

initialmat = matrix(0, numbercountries, 6)



##2 for JHK fixed data, 3 for Wordo data
for(i in 1:numberactualcountries){
  initialmat[i,] = round(as.matrix(JHKinitial[[i]])[1,],digits=0)
}
##########Notice initial mat here corresponding to first time pass 100 cases of each country
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



############DONE PART I, THE PREPARATION STEP ###########



#######Part II. #Estimate the parameters of the considering country.


country = countries[[countryconsider]]


###############2. Define uniform prior for parameters
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


#################Borrowing preliminaries parameters from the naive model
prelimposterior = as.matrix(JHKinitial[[countryconsider]][-1,])


###########5. Define the standardize function to get sd of A, R, D
thresholddeath = .05
thresholdconfirmed = .05
standardize_seq = function(k,country){
  
  sd = rep(0, country$durationtravel) 
  
  
  Umat = matrix(0,k,country$durationtravel)
  Dmat = matrix(0,k,country$durationtravel)
  
  
  i <- 1
  while (i <= k ) {
    
    index= sample(1:nrow(prelimposterior) ,1)
    theta = prelimposterior[index,]
    h1 = theta[2]/2
    theta[2] = theta[2]+runif(1,-h1,h1)
    u = stochastic_marginalestimate(theta,country)
    #########making sure data generated not too far away from the truth
    checkval = abs(u[nrow(u),5] - country$data[nrow(u),3])/country$data[nrow(u),3]    
    checkval1 = abs(sum(u[nrow(u),3:5]) - sum(country$data[nrow(u),]))/sum(country$data[nrow(u),])
    K1 = round(nrow(u)/2,digits=0)
    checkval2 = abs(u[K1,5] - country$data[K1,3])/country$data[K1,3]
    checkval3 = abs(sum(u[K1,3:5]) - sum(country$data[K1,]))/sum(country$data[K1,])
    CHECK = max(checkval, checkval1,checkval2,checkval3)
 
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

###########Define the indicator function to remove weird jumps
indicator = function(a){
  a[a<0]=0
  return(a)
}

##########Define the ave4 function to smooth daily cases
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
####################Define the distance


distance1 = function(theta, country){
  
  
  sim <- CONETTravel::stochastic_marginalestimate(theta, country)
  ################Learn beta, delta########################
  U = rowSums(sim[,3:5])
  U1 = diff(U,1) #total cases
  Ur = diff(rowSums(country$data),1)
  Ur = indicator(Ur)
  Ur = ave4func(Ur) 
  D1 = diff(sim[,5],1) #deaths
  Dr = diff(country$data[,3],1)
  Dr = indicator(Dr)
  Dr = ave4func(Dr)
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

#####6. Run ABC

abcposterior1 <- abc_start(
  prior,
  distance1, 
  distance_args = country,
  method = "RABC",
  control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
  output_control = list(print_output = TRUE)
)

####################Use the Euclidean distance
distance2 = function(theta, country){


  sim <- CONETTravel::stochastic_marginalestimate(theta, country)
  ################Learn beta, delta########################
  U = rowSums(sim[,3:5])
  U1 = diff(U,1) #total cases
  Ur = diff(rowSums(country$data),1)
  Ur = indicator(Ur)
  Ur = ave4func(Ur)
  D1 = diff(sim[,5],1) #deaths
  Dr = diff(country$data[,3],1)
  Dr = indicator(Dr)
  Dr = ave4func(Dr)
  ###EU distances
   output1 <- sum((U1 - Ur)^2) # Total Confirmed daily
   output2 <- sum((D1 - Dr)^2) #Total Death daily
  ####################
  output <- 1/nrow(country$data)*(output1^.5 + output2^.5)

  return(output)
}

abcposterior2 <- abc_start(
  prior,
  distance2,
  distance_args = country,
  method = "RABC",
  control = list(prior_eval = prior_eval,  n = particles, pacc_final = prob),
  output_control = list(print_output = TRUE)
)





#########CHOOSING THE BEST PARTICLE FROM POSTERIOR##############
abcposterior1 = as.matrix(abcposterior1)
abcposterior2 = as.matrix(abcposterior2)

#################
posteriorandinitialmat = rbind(country$x_ini, abcposterior1, abcposterior2)
fname = paste('GlobalPosteriorandInitial_JHKdata_countryconsider_',countryconsider,'_thresholddeath1',thresholddeath,".txt",sep="")
write.table(posteriorandinitialmat,file=fname, row.names = F)


