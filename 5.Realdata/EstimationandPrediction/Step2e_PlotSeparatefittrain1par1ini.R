#Created:March 5
#Last edit: March 7 4h35 p.m.
##Purpose: Plotting the fit based on the outputs at Step 2 and save the best particle and initial for step 3######

args <- commandArgs(trailingOnly = TRUE)
DUMMY  <- as.integer(args[1])




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



 




load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/AllposteriorandinitialforPredictionStep2.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/coviddataJanJune20.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/flightdataJanJune20.Rdata")
load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/AllPreliminariesposteriorandinitialStep1.Rdata")
# ##Get initials from the independent model with kappa

countryindex = c(1:92)

for(tt in 1:length(countryindex)){  

countryconsider = countryindex[tt]
#obtain the posterior of the considering country
###fname = paste('MarginalGlobalInitialforPredictionandPosterior_JHKdata_countryconsider_',countryconsider,".txt",sep="")

JHKposteriorandini = JHKPosIniStep2[[countryconsider]]
###

travelout_datadivided = flightdataJanJune20

#Standradize covid data and flight data with covid data in terms of time
mycompletedata = coviddataJanJune20 
temp = length(travelout_datadivided) - nrow(mycompletedata[[1]])
temp1 = 1:temp
travelout_datadivided = travelout_datadivided[-temp1]

###########Standardize the covid data and flight data up to a cutoff time we want
numbercountries = ncol(travelout_datadivided[[1]] )
numberactualcountries = numbercountries-1 ##number actual remove the fake 1 ZZZ
cutofftime = nrow(mycompletedata[[1]])   
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

###Get initial conditions from Preliminaries can also do JHKIniStep2
initialmat = matrix(0, numbercountries, 6)

for(i in 1:numberactualcountries){
  initialmat[i,] = round(as.matrix(JHKinitial[[i]])[1,],digits=0) 
}


##########
casecutoff = 500 #CHOOSE CUTOFF AS 500 CONFIRMED CASES FOR THE START DAY

##All start when it own first pass 100 cases
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


mydatacutoff = list() # create a list of data for each country from 1 to n = 93 this case

for(i in 1:numbercountries){
  mydatacutoff[[i]] = mycompletedata[[i]][startday:nrow(mycompletedata[[1]]),]
}  ##all start as the same day as the considering country has more than 100 cases.

#2nd round extracting the corresponding flight data with covid data, cutoff based on the considering data first pass 100 cases
travelout_datadivided = travelout_datadivided[startday:nrow(mycompletedata[[1]])]


### ACHIEVED. ALL DATA AND TRAVEL DATA NOW ARE TRIM TO THE CONSIDERING COUNTRY

dataperiod = nrow(mydatacutoff[[countryconsider]]) -16+1
trainingtime = dataperiod # use the whole data as t get all the information at one
mydatacutoff_train = list() # create a list of data for each country from 1 to n = 93 this case

for(i in 1:numbercountries){
  mydatacutoff_train[[i]] = mydatacutoff[[i]][1:trainingtime,]
  
}



mydatacutoffadapt_train = list()
for(i in 1:numberactualcountries){
  tmp = mydatacutoffadapt[[i]]
  h1 = nrow(tmp) # use the whole data first
  mydatacutoffadapt_train[[i]] =  tmp[1:h1,]
}



travelout_datadivided_train = travelout_datadivided[1:trainingtime]  #Number travelers from 1 country to another during the training


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
#########PART II. MATCHING#############################

###Matching by simple best quantile
##1a,b. Daily matching, Extract the best particles based on alpha
country = countries[[countryconsider]]
abcposterior = as.matrix(JHKposteriorandini[2:1001,])
theconsidercountry = coviddataJanJune20$country[countryconsider]
###############replaced here by a loop to get the best particle with the smallest average distancereal
d11 = nrow(country$data) - 14
country_train =  list(data = country$data[1:d11,], x_ini = country$x_ini, 
                      travelin_compartments = country$travelin_compartments[1:d11,],
                      total_out = country$total_out[1:d11], durationtravel = d11)


mylistall = list()
realizations = 100 #number realizations
averageall = rep(0,nrow(abcposterior))
for(particle in 1:nrow(abcposterior)){
  
  mylistall[[particle]] = list()
  average = 0
  theta = abcposterior[particle,]
  for (i in 1:realizations){
    # inp = list(ini = country$x_ini, duration = country$durationtravel)
    # data = stochasticmodel_1country(theta, inp)
    data = stochastic_marginalestimate(theta, country_train)
    tmp = sum((rowSums(data[,3:5]) - rowSums(country_train$data))^2) + sum((data[,5] - country_train$data[,3])^2)
    
    average = average + tmp
    mylistall[[particle]][[i]] = data
  }
  averageall[particle] = average/realizations
  
} #end loop for paricle

indexbest = min(which.min(averageall)) # choose any smallest in case there are 2 or more with minimum
########################
mylist = mylistall[[indexbest]]

########add time


wholeperiod = rev(seq(as.Date("2020-06-14"), length =  dataperiod, by = "-1 day"))
trainingtime = length(wholeperiod)-14
ttrain = 1:trainingtime
trainingperiod = wholeperiod[ttrain]



##1a,b. Accumulated matching, Extract the best particles based on alpha

activemat =  matrix(0,nrow=nrow(mylist[[1]]),ncol=realizations)

for (i in 1:realizations){
  
  activemat[,i] = rowSums(mylist[[i]][,3:5]) # Use accumulated to fixed
  
}

###############

scalefactor = 1

activemat = activemat/scalefactor # scale to to samller number
quantileactive = rowQuantiles( activemat, probs =c(.025,.5,.975))
############# real data

active =  rowSums(country$data)/scalefactor  # Use accumulated to fixed
#active =  c(0,diff(rowSums(country$data)/scalefactor,1))



#########create data set
activedata = cbind(quantileactive, active)
activedata1 = log(activedata)/log(10)


colnames(activedata) = c("lower","median", "upper", "datareport")
colnames(activedata1) = c("lower","median", "upper", "datareport")

dataactive = data.frame( trainingperiod,activedata)

dataactive1 = data.frame( trainingperiod,activedata1)
##Now extract only the training part



##Plot real data vs estimation ##



plt1a = ggplot(data=dataactive, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+ theme_minimal()+
  geom_line(aes(y=median), colour="blue") +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs(title = theconsidercountry, x = "", y = "Accumulated confirmed")

plt1b = ggplot(data=dataactive1, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+ theme_minimal()+
  geom_line(aes(y=median), colour="blue") +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs( x = "Time period", y = "Log scale ccumulated confirmed")



#################################
############Matching for daily death confirmed

deathmat =  matrix(0,nrow=nrow(mylist[[1]]),ncol=realizations)

for (i in 1:realizations){
  deathmat[,i] = mylist[[i]][,5] #use accumulated to fix
}


quantiledeath = rowQuantiles( deathmat, probs =c(.025,.5,.975))

death =  country$data[,3] # real death data

deathdata = cbind(quantiledeath, death)
deathdata1 = log(deathdata)/log(10)
colnames(deathdata) = c("lower","median", "upper", "datareport")
colnames(deathdata1) = c("lower","median", "upper", "datareport")

######Create the data
datadeath = data.frame( trainingperiod,deathdata)
datadeath1 = data.frame( trainingperiod,deathdata1)

datadeath = datadeath[1:d11,] #extract only the training data
datadeath1 = datadeath1[1:d11,] #extract only the training data

##Plot real data vs estimation ##


plt1c = ggplot(data=datadeath, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+ theme_minimal()+
  geom_line(aes(y=median), colour="blue") + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs( x = "", y = "Accumulated death")

plt1d = ggplot(data=datadeath1, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+ theme_minimal()+
  geom_line(aes(y=median), colour="blue") + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs( x = "Time period", y = "Log scale accumulated death")

plt1 = plot_grid( plt1a,plt1c,plt1b,plt1d, nrow=2)


##############################################################
fname = paste('Fitoneparoneinicountry_',countryconsider,".png",sep="")
ggsave(fname,plot=plt1)




###2. Save the best initial and best parameters for testing based on training


country1 = countries[[countryconsider]]
durationtest = 14 +1 # add 1 to make things go back to May 31 as the initial day
starttest= nrow(country1$data) - 14 
trainingperiod = nrow(country1$data) - 14
##Extract the best realization in mylist, the corresponding smallest average 
alldistances = rep(0, realizations)

for (i in 1:realizations){
  
  data = as.matrix(mylist[[i]])
  tmp = sum((rowSums(data[,3:5]) - rowSums(country_train$data))^2) + sum((data[,5] - country_train$data[,3])^2)
  alldistances[i] = tmp
}

indexbest1 = min(which.min(alldistances)) # choose any smallest in case there are 2 or more with minimum

newinitial = mylist[[indexbest1]][nrow(mylist[[1]]),]
newinitial[3:5] = country1$data[starttest,] # replace A,R,D by real data for a better fit
thetahat = abcposterior[indexbest,]
tmp = rbind(newinitial, thetahat)
fname1 = paste('Bestinitialandparametersforstep3_',countryconsider,".txt",sep="")
write.table(tmp,file=fname1, row.names = F)

}







