args <- commandArgs(trailingOnly = TRUE)
dummy  <- as.integer(args[1])
n1  <- as.integer(args[2])

set.seed(101)
prior1 <- function(n){
  data.frame(
    alpha0 = 0,
    alpha = runif(n,0,2),
    beta = runif(n,0,1),
    delta=runif(n,0,1),
    eta=1,
    gamma=runif(n,0,1),alpha0 = 0,
    alpha = runif(n,0,2),
    beta = runif(n,0,1),
    delta=runif(n,0,1),
    eta=1,
    gamma=runif(n,0,1),alpha0 = 0,
    alpha = runif(n,0,2),
    beta = runif(n,0,1),
    delta=runif(n,0,1),
    eta=1,
    gamma=runif(n,0,1)
    #gamma=1.25/15
  )
}

thetas = prior1(n1)
d=84
#############Travelers array, this is a matrix of number of travel out each day from 1 country to another
## Average travel rate: 40/(365*328) = .0003341129 ~= 3*10^(-4)
P1 = 10^7
P2 = 3*10^6
P3 = 2*10^6


T1 =  P1*40/(365*328) # Total travel out
T2 =  P2*40/(365*328)
T3 =  P3*40/(365*328)
SdT1 = T1*.05 # Vary the Sd by 5%
SdT2 = T2*.05
SdT3 = T3*.05

##Create travelers array for each day
prior2 <- function(n){
  data.frame(
    Travel1 = rnorm(n, T1, SdT1),
    Travel2 = rnorm(n, T2, SdT2),
    Travel3 = rnorm(n, T3, SdT3)
  )
}
traveloutdat = prior2(d)
traveloutdat = round(traveloutdat,digits=0)

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




tmp = proc.time()

stofunc_travel =  function(adjtheta){

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




status_matrix = matrix(0, nrow = d,ncol = length(initial_corona))
status_matrix1 =  status_matrix

f_in = matrix(0, nrow = d,ncol = length(initial_corona))
f_in
f_out = matrix(0, nrow = d, ncol = length(initial_corona))
status_matrix[1,] = initial_corona



for (i in 2:d){
  for (j in 1:3){
    
    c1 = (j-1)*6 + 1
    c2 = j*6
    x = status_matrix[(i-1),c1:c2]
    ##Updated traveling flow
    #Number out from country j
    out = traveloutdat[i,j]
    if (x[1]+x[2] > 0){
      outj = c(round(out*x[1]/(x[1]+x[2]),digits=0), round(out*x[2]/(x[1]+x[2]),digits=0), 0,0,0,0)
      #########
      f_out[i,c1:c2] = outj 
    }
    #travel out from country i with sick and susceptible
    
    theta = adjtheta[c1:c2]
    ##Generating Poisson values based on hazard functions
    y1 = rpois(1, harzard1(x,theta))
    # 
    y2 =  rpois(1, harzard2(x,theta))
    # 
     y3 = rpois(1, harzard3(x,theta))
    # 
    y4 =  rpois(1, harzard4(x,theta))
    # 
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
    
   
    
    status_matrix[i,c1:c2] = x
    
  }
  
  
  
  ######################No regulation ####
  f_in1 = f_out[i,7:12]*P1/(P1 + P3) + f_out[i,13:18]*P1/(P1+P2)
  f_in2 = f_out[i,1:6]*P2/(P2+P3) + f_out[i,13:18]*P2/(P1+P2)
  f_in3 = f_out[i,1:6]*P3/(P2+P3) + f_out[i,7:12]*P3/(P1+P3)
   
  
  f_in[i,] = c(round(f_in1,digits = 0), round(f_in2,digits=0), round(f_in3,digits=0))
  ########
  update = status_matrix[i,]  +  f_in[i,] - f_out[i,]
  update[update<0.1]=0
  status_matrix[i,] = update
  
}  
return(status_matrix)
}


#####
datalist = list()
thetalist = list()
for(k in 1:n1){
  
  
  theta0 = as.numeric(thetas[k,])
  ###Stochastic
  mydata =  stofunc_travel(theta0)
  Confirmed1 = mydata[,3] + mydata[,4] + mydata[,5]
  Confirmed2 = mydata[,9] + mydata[,10] + mydata[,11]
  Confirmed3 = mydata[,15] + mydata[,16] + mydata[,17]
  minsus = min(mydata[,c(1,7,13)] )
  
  
  # #####
  if (minsus > .4&& 
      mydata[(d-28),2] > 1 && mydata[(d-28),8] > 1 && mydata[(d-28),14] > 1 &&
      mydata[d,4] > 2*mydata[d,5] &&  mydata[d,10] > 2*mydata[d,11] &&  mydata[d,16] > 2*mydata[d,17] &&
      mydata[d,5]< .3*Confirmed1[d] && mydata[d,11]< .3*Confirmed2[d] && mydata[d,17]< .3*Confirmed3[d] && 
      Confirmed1[d] < .5*P1 && Confirmed1[d] < .5*P2 && Confirmed3[d] < .5*P3 ) {
  #   
  #   
    thetalist[[k]] = theta0
    datalist[[k]] = mydata
  #   
  }
}

thetas_3travel = thetalist[lengths(thetalist)!=0]
datas_3travel = datalist[lengths(datalist)!=0]
travelout_3dat = traveloutdat



save(thetas_3travel, file="thetas_3travel.RData")
save(datas_3travel, file="datas_3travel.RData")
save(travelout_3dat, file="travelout_3dat.RData")
#########
