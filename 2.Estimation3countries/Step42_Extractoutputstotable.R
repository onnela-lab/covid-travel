#The code is used to check HPC OUTPUT FOR 3 COUNTRIES, CREATED OCT 19.
setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/Estimation3countries")
particles = 1000
prob=.01

#library(tidyr)

fname1 = paste("ZC3countriesDistanceoutput_particles",particles,"_pacc",prob,".txt",sep="")
data = as.matrix(read.table(fname1,header=T))
data = abs(data)

dim(data)
data[1:7,1:7]
country1 = data[,1:150]
country2 = data[,151:300]
country3 = data[,301:450]

##Column 1:4  : Median Relative Bias; Column 5:8: Mean Relative Bias
##Column 9:16: theta cover; 
#column 17:24: Coverage;
#column 25:30: Truth; 


####################Print Text output
library(xtable)


biasfun = function(mydat, k){
  
  mydat1 = mydat[,1:30]
  mydat2 = mydat[,31:60]
  mydat3 = mydat[,61:90]
  mydat4 = mydat[,91:120]
  mydat5 = mydat[,121:150]
 
  a1 =1+k*4
  b1 = (k+1)*4

  mydat1 = mydat1[,a1:b1]*mydat1[,c(26,27,28,30)]
  mydat3 = mydat3[,a1:b1]*mydat3[,c(26,27,28,30)]
  mydat4 = mydat4[,a1:b1]*mydat4[,c(26,27,28,30)]
  
  
  
  
  H1a = rbind(apply(mydat1,2,mean))
  h1 = apply(mydat1,1,mean)
  H1b  = t(t(c(mean(h1))))
  H1 = cbind(H1b, H1a)
  ####
 
  
  H3a = rbind(apply(mydat3,2,mean))
  h3 = apply(mydat3,1,mean)
  H3b  = t(t(c(mean(h3))))
  H3 = cbind(H3b, H3a)
  ########
  H4a = rbind(apply(mydat4,2,mean))
  h4 = apply(mydat4,1,mean)
  H4b  = t(t(c(mean(h4))))
  H4 = cbind(H4b, H4a)
  ##########
 
  
  
  H = rbind(H4,H3,H1)
  
  
  colnames(H) = c("Ave", "alpha", "beta", "delta", "gamma")
  # Mean and Sd each parameter
  
  
  H = as.data.frame(H)
  tab <-xtable(H)
  
  digits(tab) <- rep(3,6)
  
  ##############
  
  
  
  
  print(tab,include.rownames=FALSE)
  
}





##############Bias#########
rebiasfun = function(mydat,k){
  
  
  
  
  mydat1 = mydat[,1:30]
  mydat2 = mydat[,31:60]
  mydat3 = mydat[,61:90]
  mydat4 = mydat[,91:120]
  mydat5 = mydat[,121:150]
  a1 =1+k*4
  b1 = (k+1)*4
  #########RELATIVE BIAS#########
  H1a = rbind(apply(mydat1[,a1:b1],2,mean))
  # Mean and Sd of avarage
  h1 = apply(mydat1[,a1:b1],1,mean)
  H1b  = t(t(c(mean(h1))))
  H1 = cbind(H1b, H1a)
  ####
  
  
  H3a = rbind(apply(mydat3[,a1:b1],2,mean))
  h3 = apply(mydat3[,a1:b1],1,mean)
  H3b  = t(t(c(mean(h3))))
  H3 = cbind(H3b, H3a)
  ########
  H4a = rbind(apply(mydat4[,a1:b1],2,mean))
  h4 = apply(mydat4[,a1:b1],1,mean)
  H4b  = t(t(c(mean(h4))))
  H4 = cbind(H4b, H4a)
  ######

  
  
  
  H = rbind(H4,H3,H1)
  
  
  colnames(H) = c("Ave", "alpha", "beta", "delta", "gamma")
  # Mean and Sd each parameter
  
  
  H = as.data.frame(H)
  tab <-xtable(H)
  
  digits(tab) <- rep(3,6)
  
  ##############
  
  
  
  
  print(tab,include.rownames=FALSE)
  
}
#########Prcocess each country
#Accurate

country = country2
#Median Accuracy
biasfun(country,0)

rebiasfun(country,0)

#####Mean Accuracy
biasfun(country,1)

rebiasfun(country,1)


##Cover Range
rebiasfun(country,4)

rebiasfun(country,5)

##Cover rate
rebiasfun(country,2)

rebiasfun(country,3)
#########Combine country, equal weight

country = 1/3*(country1 + country2 + country3)

#Median Accuracy
biasfun(country,0)

rebiasfun(country,0)

#####Mean Accuracy
biasfun(country,1)

rebiasfun(country,1)


##Cover Range
rebiasfun(country,4)

rebiasfun(country,5)

##Cover rate
rebiasfun(country,2)

rebiasfun(country,3)

##################

#########Combine country, equal weight
P1 = 10^7
P2 = 3*10^6
P3 =2*10^6

country = P1/(P1+P2+P3)*country1 +  P2/(P1+P2+P3)*country2 +  P3/(P1+P2+P3)*country3 

#Median Accuracy
biasfun(country,0)

rebiasfun(country,0)

#####Mean Accuracy


##Cover Range
rebiasfun(country,4)

rebiasfun(country,5)

##Cover rate
rebiasfun(country,2)

rebiasfun(country,3)




