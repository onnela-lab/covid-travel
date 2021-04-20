
set.seed(2)
prior1 <- function(n){
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
n1=2000000
thetas = prior1(n1)
tmp = proc.time()


P =10^7
I1 = 15
A1 = 13
S1 = P-I1-A1
x1 = c(S1,I1,A1,0,0,0)# State corresponding S,I,A,R,D,Ru
d = 84
status_matrix1 = matrix(0,nrow =d ,ncol=6)
status_matrix1[1,] = x1

status_matrix2 = matrix(0,nrow =d ,ncol=6)
status_matrix2[1,] = x1

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


sto_generating_func1 = function(theta){
  
  
  for (i in 2:d){
    
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

#########
generating_func1 = function(theta){
  
  
  for (i in 2:d){
    
    x = status_matrix2[(i-1),]
    
    
    
    
    ##Generating Poisson values based on hazard functions
    y1 =  harzard1(x,theta)
    
    y2 =   harzard2(x,theta)
    
    
    y3 =  harzard3(x,theta)
    
    
    y4 =   harzard4(x,theta)
    
    
    
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
    status_matrix2[i,] = x
  }
  
  return(status_matrix2)
}

#####
datalist = list()
thetalist = list()
for(j in 1:n1){
  theta0 = as.numeric(thetas[j,])
  ###Stochastic
  mydata =  sto_generating_func1(theta0)
  Confirmed = mydata[,3] + mydata[,4] + mydata[,5]
  #####
  if (mydata[(d-28),1] > 0 && mydata[(d-28),2] > 0 && mydata[d,4] > 2*mydata[d,5] && mydata[d,5]> .01*Confirmed[d] && mydata[d,5]< .3*Confirmed[d] && Confirmed[d] < .5*P && min(theta0[-c(1,5)]) >.01 ) {
    
    
    thetalist[[j]] = theta0
    datalist[[j]] = mydata
    
  }
}

thetas_list = thetalist[lengths(thetalist)!=0]
datas_list = datalist[lengths(datalist)!=0]
save(thetas_list, file="thetas_list.RData")
save(datas_list, file="datas_list.RData")

