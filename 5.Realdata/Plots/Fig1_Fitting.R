#This code is used to plot and save posterior with alpha varied every 4 weeks
#Date started: March 21
#Last edited: March 21, 4h50 p.m.
library(ggplot2)
library(gridExtra)

myfunctionUS = function(countryconsider, n1,confirmupper, deathupper, size1,sizetext, scalefactor,scalefactor1){

load("coviddataJanJune20.Rdata")
load("flightdataJanJune20.Rdata")
load("AllPreliminariesposteriofromstep1changepoint.Rdata")
fname1 = paste('Allinitialfromstep1changepoint',".txt",sep="")

initialmat  = read.table(fname1,header=T)
initialmat = as.matrix(initialmat)



theconsidercountry = coviddataJanJune20$country[countryconsider]
#1. set up Acceptance rate and number of particles to run ABC



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
library(scales)
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



fname = paste('Posteriorproposedwithmultiplealphachangepointcountryconsider_',countryconsider,".txt",sep="")
abcposterior = read.table(fname,header=T)
abcposterior = as.matrix(abcposterior)


############DONE PART II, THE ESTIMATION STEP##########
#########PART III. MATCHING#############################
######################################################################


mylistall = list()
realizations = 100 #number realizations
averageall = rep(0,n1)
for(particle in 1:n1){
  
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


###################MATCHING WITH ACCUMULATED
############Matching for daily active confirmed
trainingtime = nrow(mylist[[1]])
activemat =  matrix(0,nrow=trainingtime, ncol=realizations)

for (i in 1:realizations){
  
  activemat[,i] = rowSums(mylist[[i]][,3:5])
  
}




activemat = activemat/scalefactor1 # scale to to smaller number
quantileactive = rowQuantiles( activemat, probs =c(.025,.5,.975))
############# real data

confirmed =  rowSums(country$data)/scalefactor1
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

##############
scientific <- function(x){
  ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
}



##Plot real data vs estimation ##
plt1a = ggplot(data=dataconfirmed, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+
  geom_line(aes(y=median), colour="blue") +
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs(title = theconsidercountry, x = "", y = "Confirmed")+
  scale_y_continuous(breaks = c(0,1000000, 2000000),
                      label = c("0","1e+6","2e+6"), limits=c(0,confirmupper))+
  theme_classic()+theme(plot.title = element_text(face = "bold",size=sizetext),
    axis.text=element_text(size=size1),
    axis.title.x = element_text( size=sizetext, face ="bold"),
    axis.title.y = element_text(size=sizetext, face="bold"),
    legend.title = element_text( size=sizetext, 
                                 face="bold"), axis.text.x=element_blank())+  
  scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                          as.Date('2020-05-31')),
               label = c("Mar 1","Apr 15", "May 31"),
               limits = as.Date(c('2020-03-01','2020-06-01')))



###################################################
#################################
############Matching for daily death confirmed

deathmat =  matrix(0,nrow=trainingtime,ncol=realizations)

for (i in 1:realizations){
  
  deathmat[,i] =  mylist[[i]][,5]
}


quantiledeath = rowQuantiles( deathmat, probs =c(.025,.5,.975))
############# real data

death =  country$data[,3]

#########create data set
deathdata = cbind(quantiledeath, death)
deathdata = deathdata/scalefactor
colnames(deathdata) = c("lower","median", "upper", "datareport")

######Create the data
datadeath = data.frame( trainingperiod,deathdata)



##Plot real data vs estimation ##
plt1b = ggplot(data=datadeath, aes(x=trainingperiod, y=datareport)) +
  geom_line(aes(y=datareport), colour="red")+ 
  geom_line(aes(y=median), colour="blue") + 
  geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
  labs( x = "", y = "Deaths ")+
  scale_y_continuous(breaks = c(0,50000, 100000),
                      label = c("0","5e+4","1e+5"), limits=c(0,deathupper))+
  theme_classic()+theme(
    axis.text=element_text(size=size1),
    axis.title.x = element_text( size=sizetext, face ="bold"),
    axis.title.y = element_text(size=sizetext, face="bold"),
    legend.title = element_text( size=sizetext, 
                                 face="bold"), axis.text.x=element_blank())+  
  scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                          as.Date('2020-05-31')),
               label = c("Mar 1","Apr 15", "May 31"),
               limits = as.Date(c('2020-03-01','2020-06-01')))

return(list(plt1a,plt1b))}




myplotfunc = function(countryconsider, n1,confirmupper, deathupper,size1,sizetext, scalefactor, scalefactor1){
  setwd("C:/Users/thl902/Desktop/Investigating/changepoint/PlottingMaintext")
  load("coviddataJanJune20.Rdata")
  load("flightdataJanJune20.Rdata")
  load("AllPreliminariesposteriofromstep1changepoint.Rdata")
  fname1 = paste('Allinitialfromstep1changepoint',".txt",sep="")
  
  initialmat  = read.table(fname1,header=T)
  initialmat = as.matrix(initialmat)
  
  
  
  theconsidercountry = coviddataJanJune20$country[countryconsider]
  #1. set up Acceptance rate and number of particles to run ABC
  
  
  
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
  library(scales)
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
  
  
  
  fname = paste('Posteriorproposedwithmultiplealphachangepointcountryconsider_',countryconsider,".txt",sep="")
  abcposterior = read.table(fname,header=T)
  abcposterior = as.matrix(abcposterior)
  
  
  ############DONE PART II, THE ESTIMATION STEP##########
  #########PART III. MATCHING#############################
  ######################################################################
  
  
  mylistall = list()
  realizations = 100 #number realizations
  averageall = rep(0,n1)
  for(particle in 1:n1){
    
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
  

  ###################MATCHING WITH ACCUMULATED
  ############Matching for daily active confirmed
  trainingtime = nrow(mylist[[1]])
  activemat =  matrix(0,nrow=trainingtime, ncol=realizations)
  
  for (i in 1:realizations){
    
    activemat[,i] = rowSums(mylist[[i]][,3:5])
    
  }
  
  ###############
 
  
  activemat = activemat/scalefactor1 # scale to to smaller number
  quantileactive = rowQuantiles( activemat, probs =c(.025,.5,.975))
  ############# real data
  
  confirmed =  rowSums(country$data)/scalefactor1
  #########create data set
  confirmeddata = cbind(quantileactive, confirmed)
  
  colnames(confirmeddata) = c("lower","median", "upper", "datareport")

  
  #############add time
  
  t = 1:trainingtime
  wholeperiod = rev(seq(as.Date("2020-06-30"), length =  dataperiod, by = "-1 day"))
  trainingperiod = wholeperiod[1:trainingtime]
  
  ######Create the data
  dataconfirmed = data.frame( trainingperiod,confirmeddata)
  

  
  ##############
  scientific <- function(x){
    ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", " %*% 10^", scientific_format()(x)))))
  }
  
  
  
  ##Plot real data vs estimation ##
  plt1a = ggplot(data=dataconfirmed, aes(x=trainingperiod, y=datareport)) +
    geom_line(aes(y=datareport), colour="red")+
    geom_line(aes(y=median), colour="blue") +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
    labs(title = theconsidercountry, x = "", y = "Confirmed")+
    scale_y_continuous(breaks = c(0,300000, 600000),
                        label = c("0","3e+5","6e+5"), limits=c(0,confirmupper))+
    theme_classic()+theme(
      axis.text=element_text(size=size1),plot.title = element_text(face = "bold",size=sizetext),
      axis.title.x = element_text( size=sizetext, face ="bold"),
      axis.title.y = element_text(size=sizetext, face="bold"),
      legend.title = element_text( size=sizetext, 
                                   face="bold"), axis.text.x=element_blank())+  
    scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                            as.Date('2020-05-31')),
                 label = c("Mar 1","Apr 15", "May 31"),
                 limits = as.Date(c('2020-03-01','2020-06-01')))
  
  
  plt1b = ggplot(data=dataconfirmed, aes(x=trainingperiod, y=datareport)) +
    geom_line(aes(y=datareport), colour="red")+
    geom_line(aes(y=median), colour="blue") +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
    labs(title = theconsidercountry, x = "", y = "")+
    scale_y_continuous(breaks = c(0,300000, 600000),
                        label = c("0","3e+5","6e+5"), limits=c(0,confirmupper))+
    theme_classic()+theme(
      axis.text=element_text(size=size1),plot.title = element_text(face = "bold",size=sizetext),
      axis.title.x = element_text( size=sizetext, face ="bold"),
      axis.title.y = element_text(size=sizetext, face="bold"),
      legend.title = element_text( size=sizetext, 
                                   face="bold"), axis.text.x=element_blank())+  
    scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                            as.Date('2020-05-31')),
                 label = c("Mar 1","Apr 15", "May 31"),
                 limits = as.Date(c('2020-03-01','2020-06-01')))
  
  plt1b1 = ggplot(data=dataconfirmed, aes(x=trainingperiod, y=datareport)) +
    geom_line(aes(y=datareport), colour="red")+
    geom_line(aes(y=median), colour="blue") +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
    labs(title = theconsidercountry, x = "", y = "")+ylim(0,confirmupper)+
    theme_classic()+theme(plot.title = element_text(face = "bold",size=sizetext),
      axis.text=element_text(size=size1),
      axis.title.x = element_text( size=sizetext, face ="bold"),
      axis.title.y = element_text(size=sizetext, face="bold"),
      legend.title = element_text( size=sizetext, 
                                   face="bold"), axis.text.y=element_blank(), axis.text.x=element_blank())+  
    scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                            as.Date('2020-05-31')),
                 label = c("Mar 1","Apr 15", "May 31"),
                 limits = as.Date(c('2020-03-01','2020-06-01')))
  ###################################################
  #################################
  ############Matching for daily death confirmed
  
  deathmat =  matrix(0,nrow=trainingtime,ncol=realizations)
  
  for (i in 1:realizations){
    
    deathmat[,i] =  mylist[[i]][,5]
  }
  
  
  quantiledeath = rowQuantiles( deathmat, probs =c(.025,.5,.975))
  ############# real data
  
  death =  country$data[,3]
  #########create data set
  deathdata = cbind(quantiledeath, death)
  deathdata = deathdata/scalefactor
  colnames(deathdata) = c("lower","median", "upper", "datareport")
  
  ######Create the data
  datadeath = data.frame( trainingperiod,deathdata)
  
  
  
  ##Plot real data vs estimation ##
  plt1c = ggplot(data=datadeath, aes(x=trainingperiod, y=datareport)) +
    geom_line(aes(y=datareport), colour="red")+ 
    geom_line(aes(y=median), colour="blue") + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
    labs( x = "", y = "")+scale_y_continuous(breaks = c(0,20000, 40000),
                                             label = c("0","2e+4","4e+4"), limits=c(0,deathupper))+
    theme_classic()+theme(
      axis.text=element_text(size=size1),
      axis.title.x = element_text( size=sizetext, face ="bold"),
      axis.title.y = element_text(size=sizetext, face="bold"),
      legend.title = element_text( size=sizetext, 
                                   face="bold"), axis.text.x=element_blank())+  
    scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                            as.Date('2020-05-31')),
                 label = c("Mar 1","Apr 15", "May 31"),
                 limits = as.Date(c('2020-03-01','2020-06-01')))
  ##Plot real data vs estimation ##
  plt1c1 = ggplot(data=datadeath, aes(x=trainingperiod, y=datareport)) +
    geom_line(aes(y=datareport), colour="red")+ 
    geom_line(aes(y=median), colour="blue") + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
    labs( x = "", y = "Deaths")+scale_y_continuous(breaks = c(0,20000, 40000),
                                                   label = c("0","2e+4","4e+4"), limits=c(0,deathupper))+
    theme_classic()+theme(
      axis.text=element_text(size=size1),
      axis.title.x = element_text( size=sizetext, face ="bold"),
      axis.title.y = element_text(size=sizetext, face="bold"),
      legend.title = element_text( size=sizetext, 
                                   face="bold"), axis.text.x=element_blank())+  
    scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                            as.Date('2020-05-31')),
                 label = c("Mar 1","Apr 15", "May 31"),
                 limits = as.Date(c('2020-03-01','2020-06-01')))
  ##Plot real data vs estimation ##
  plt1d = ggplot(data=datadeath, aes(x=trainingperiod, y=datareport)) +
    geom_line(aes(y=datareport), colour="red")+ 
    geom_line(aes(y=median), colour="blue") + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
    labs( x = "", y = "")+ylim(0,deathupper)+ 
    theme_classic()+theme(
      axis.text=element_text(size=size1),
      axis.title.x = element_text( size=sizetext, face ="bold"),
      axis.title.y = element_text(size=sizetext, face="bold"),
      legend.title = element_text( size=sizetext,  
                                   face="bold"), axis.text.y=element_blank(),axis.text.x=element_blank())+  
    scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                            as.Date('2020-05-31')),
                 label = c("Mar 1","Apr 15", "May 31"),
                 limits = as.Date(c('2020-03-01','2020-06-01')))
  
  #############################################################
  plt1e = ggplot(data=datadeath, aes(x=trainingperiod, y=datareport)) +
    geom_line(aes(y=datareport), colour="red")+ 
    geom_line(aes(y=median), colour="blue") + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
    labs( x = "", y = "Deaths")+scale_y_continuous(breaks = c(0,20000, 40000),
                                                   label = c("0","2e+4","4e+4"), limits=c(0,deathupper))+
    theme_classic()+theme(
      axis.text=element_text(size=size1),
      axis.title.x = element_text( size=sizetext, face ="bold"),
      axis.title.y = element_text(size=sizetext, face="bold"),
      legend.title = element_text( size=sizetext,  
                                   face="bold"))+  
    scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                            as.Date('2020-05-31')),
                 label = c("Mar 1","Apr 15", "May 31"),
                 limits = as.Date(c('2020-03-01','2020-06-01')))
  
  plt1f = ggplot(data=datadeath, aes(x=trainingperiod, y=datareport)) +
    geom_line(aes(y=datareport), colour="red")+ 
    geom_line(aes(y=median), colour="blue") + 
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.4)+
    labs( x = "", y = "")+ylim(0,deathupper)+
     theme_classic()+theme(
      axis.text=element_text(size=size1),
      axis.title.x = element_text( size=sizetext, face ="bold"),
      axis.title.y = element_text(size=sizetext, face="bold"),
      legend.title = element_text( size=sizetext, 
                                   face="bold"), axis.text.y=element_blank())+  
    scale_x_date(breaks = c(as.Date('2020-03-01'),as.Date('2020-04-15'),
                            as.Date('2020-05-31')),
                 label = c("Mar 1","Apr 15", "May 31"),
                 limits = as.Date(c('2020-03-01','2020-06-01')))

  ##############################
  ###Save one best initial and parameters for the prediction step
  return(list(plt1a,plt1b,plt1b1,plt1c,plt1c1,plt1d,plt1e,plt1f))}


trials = 10
scalefactor =1
scalefactor1 =1
confirmupper = 2000000/scalefactor1

deathupper =  100000/scalefactor

confirmupper1 = 600000/scalefactor1

deathupper1 = 40000/scalefactor
#######################
tmp1 = proc.time()
size1 =  14 
size2 = 22
##1.USA:90
countryconsider = 90 
plt1 = myfunctionUS(countryconsider, trials, confirmupper, deathupper, size1,size2,scalefactor,scalefactor1)
plt1[[2]]
##2.BRA:16
countryconsider = 16 
plt2 = myplotfunc(countryconsider, trials, confirmupper1, deathupper1, size1,size2,scalefactor,scalefactor1)


##3.RUS:78
countryconsider = 78 
plt3 = myplotfunc(countryconsider, trials, confirmupper1, deathupper1, size1,size2,scalefactor,scalefactor1)
##4.GBR:39
countryconsider = 39 
plt4 = myplotfunc(countryconsider, trials, confirmupper1, deathupper1, size1,size2,scalefactor,scalefactor1)
##5.ESP:35
countryconsider = 35 
plt5 = myplotfunc(countryconsider, trials, confirmupper1, deathupper1, size1,size2,scalefactor,scalefactor1)
##6.ITA:51
countryconsider = 51 
plt6 = myplotfunc(countryconsider, trials, confirmupper1, deathupper1, size1,size2,scalefactor,scalefactor1)
##7.FRA:38
countryconsider = 38 
plt7 = myplotfunc(countryconsider, trials, confirmupper1, deathupper1,size1,size2,scalefactor,scalefactor1)
##8.IND:45
countryconsider = 45 
plt8 = myplotfunc(countryconsider, trials, confirmupper1, deathupper1, size1,size2,scalefactor,scalefactor1)
# ##9.DEU:28
# countryconsider = 28 
# plt9 = myplotfunc(countryconsider, trials, confirmupper1, deathupper1,size1,scalefactor,scalefactor1)
# ##10.PER:72
# countryconsider = 72
# plt10= myplotfunc(countryconsider, trials, confirmupper1, deathupper1,size1,scalefactor,scalefactor1)

#####1.USA:90
plt1a = plt1[[1]]
plt1c1 = plt1[[2]]

##2.BRA:16
plt2a = plt2[[1]]
plt2b = plt2[[2]]
plt2b1 = plt2[[3]]
plt2c = plt2[[4]]
plt2c1 = plt2[[5]]
plt2d = plt2[[6]]
plt2e = plt2[[7]]
plt2f = plt2[[8]]
##############
##3.RUS:78
plt3a = plt3[[1]]
plt3b = plt3[[2]]
plt3b1 = plt3[[3]]
plt3c = plt3[[4]]
plt3c1 = plt3[[5]]
plt3d = plt3[[6]]
plt3e = plt3[[7]]
plt3f = plt3[[8]]
##4.GBR:39
plt4a = plt4[[1]]
plt4b = plt4[[2]]
plt4b1 = plt4[[3]]
plt4c = plt4[[4]]
plt4c1 = plt4[[5]]
plt4d = plt4[[6]]
plt4e = plt4[[7]]
plt4f = plt4[[8]]

##5.ESP:35
plt5a = plt5[[1]]
plt5b = plt5[[2]]
plt5b1 = plt5[[3]]
plt5c = plt5[[4]]
plt5c1 = plt5[[5]]
plt5d = plt5[[6]]
plt5e = plt5[[7]]
plt5f = plt5[[8]]
##6.ITA:51
plt6a = plt6[[1]]
plt6b = plt6[[2]]
plt6b1 = plt6[[3]]
plt6c = plt6[[4]]
plt6c1 = plt6[[5]]
plt6d = plt6[[6]]
plt6e = plt6[[7]]
plt6f = plt6[[8]]
##7.FRA:38
plt7a = plt7[[1]]
plt7b = plt7[[2]]
plt7b1 = plt7[[3]]
plt7c = plt7[[4]]
plt7c1 = plt7[[5]]
plt7d = plt7[[6]]
plt7e = plt7[[7]]
plt7f = plt7[[8]]
##8.IND:45
plt8a = plt8[[1]]
plt8b = plt8[[2]]
plt8b1 = plt8[[3]]
plt8c = plt8[[4]]
plt8c1 = plt8[[5]]
plt8d = plt8[[6]]
plt8e = plt8[[7]]
plt8f = plt8[[8]]
##9.DEU:28
plt9a = plt9[[1]]
plt9b = plt9[[2]]
plt9b1 = plt9[[3]]
plt9c = plt9[[4]]
plt9c1 = plt9[[5]]
plt9d = plt9[[6]]
plt9e = plt9[[7]]
plt9f = plt9[[8]]
# ##10.PER:72
# plt10a = plt10[[1]]
# plt10b = plt10[[2]]
# plt10b1 = plt10[[3]]
# plt10c = plt10[[4]]
# plt10c1 = plt10[[5]]
# plt10d = plt10[[6]]
# plt10e = plt10[[7]]
# plt10f = plt10[[8]]
#############Confirmed panel1
#plt1aaUS = plt1a + theme(plot.margin = unit(c(1, 1, 1, 1), "lines"))
plt1aa = plt1a + theme(plot.margin = unit(c(0.5, 1, 1, 1), "lines"))
plt2aa <- plt2b + theme(plot.margin = unit(c(0.5, 1, 1, -2), "lines"))
plt3aa <- plt3b1 + theme(plot.margin = unit(c(0.5, 1, 1, -.5), "lines"))
plt4aa <- plt4b1 + theme(plot.margin = unit(c(0.5, 1, 1, -.5), "lines"))

#############Deaths panel1
#plt1caUS = plt1c1 + theme(plot.margin = unit(c(-1.75, 1, 1, 1), "lines"))
plt1ca = plt1c1 + theme(plot.margin = unit(c(-1.75, 1, 1, 1), "lines"))
plt2ca <- plt2c + theme(plot.margin = unit(c(-1.75, 1, 1, -2), "lines"))
plt3ca <- plt3d + theme(plot.margin = unit(c(-1.75, 1, 1, -.5), "lines"))
plt4ca <- plt4d + theme(plot.margin = unit(c(-1.75, 1, 1, -.5), "lines"))

#####

############confirmed panel 2

plt5aa <- plt5a + theme(plot.margin = unit(c(-1.75, 1, 1, 1), "lines"))
plt6aa <- plt6b1 + theme(plot.margin = unit(c(-1.75, 1, 1, 0.7), "lines"))
plt7aa <- plt7b1 + theme(plot.margin = unit(c(-1.75, 1, 1, -.5), "lines"))
plt8aa <- plt8b1 + theme(plot.margin = unit(c(-1.75, 1, 1, -.5), "lines"))
plta2 = grid.arrange(plt5aa, plt6aa, plt7aa, plt8aa,ncol=4)
#############Deaths panel2
plt5ca <- plt5e + theme(plot.margin = unit(c(-1.75, 1, 1, 1), "lines"))
plt6ca <- plt6f + theme(plot.margin = unit(c(-1.75, 1, 1, 0.7), "lines"))
plt7ca <- plt7f+ theme(plot.margin = unit(c(-1.75, 1, 1, -.5), "lines"))
plt8ca <- plt8f+ theme(plot.margin = unit(c(-1.75, 1, 1, -.5), "lines"))
#####

grid.arrange(plt1aa,plt2aa, plt3aa, plt4aa,plt1ca, plt2ca, plt3ca, plt4ca, plt5aa,
             plt6aa, plt7aa, plt8aa, plt5ca,plt6ca, plt7ca, plt8ca, ncol=4)

# grid.arrange(plt1aaUS,plt1caUS, plt1aa,plt2aa, plt3aa, plt4aa,plt1ca, plt2ca, 
#                                             plt3ca, plt4ca, plt5aa, plt6aa, plt7aa, plt8aa,
#                                             plt5ca,plt6ca, plt7ca, plt8ca, 
#              layout_matrix = rbind(c(1,1,2,2), c(3,4,5,6),c(7,8,9,10),c(11,12,13,14),
#                                    c(15,16,17,18)))
# fname1 = paste('Zfitting',".pdf",sep="")
# ggsave(fname1,plot=plt,scale=1)
# 
