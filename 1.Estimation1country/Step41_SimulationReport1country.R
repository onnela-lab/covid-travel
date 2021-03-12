setwd("C:/Users/thl902/Desktop/CONETPublicationDec27/Simulation/Estimation1country")
library(ggplot2)
library(cowplot)
particles = 1000
prob=.01

library(tidyr)







####################Print Text output
library(xtable)
library(rmatio)

texoutfun1 = function(particles,prob,k){
  
  fname1 = paste("ZDistance2output_particles",particles,"_pacc",prob,".txt",sep="")
  mydat1 = as.matrix(read.table(fname1, header=T))
  mydat1 = mydat1[1:200,]
  

  fname4 = paste("ZDistance9output_particles",particles,"_pacc",prob,".txt",sep="")
  mydat4 = as.matrix(read.table(fname4, header=T))
  mydat4 = mydat4[1:200,]
  
  
  
  a1 =1+k*4
  b1 = (k+1)*4
  
  bias1 = mydat1[,a1:b1]*mydat1[,17:20]
 
  bias4 = mydat4[,a1:b1]*mydat4[,17:20]
  
  
  
  ######### BIAS#########
  # Mean and Sd each parameter
  H1a = rbind(apply(bias1,2,mean), apply(bias1,2,sd))
  # Mean and Sd of avarage
  h1 = apply(bias1,1,mean)
  H1b  = t(t(c(mean(h1), sd(h1))))
  H1 = cbind(H1b, H1a)
  
  #Relative and parameters
  H4a = rbind(apply(bias4,2,mean), apply(bias4,2,sd))
  h4 = apply(bias4,1,mean)
  H4b  = t(t(c(mean(h4), sd(h4))))
  H4 = cbind(H4b, H4a)
  
  H = rbind(H1,H4)
  
  
  colnames(H) = c("Ave", "alpha", "beta", "delta", "gamma")
  # Mean and Sd each parameter
  
  
  H = as.data.frame(H)
  tab <-xtable(H)
  
  digits(tab) <- rep(3,6)
  
  ##############
  
  
  
  
  print(tab,include.rownames=FALSE)
  
}


texoutfun1(1000,.01,1) #median, 5 to 8

##############Bias#########
texoutfun2 = function(particles,prob,k){
  
  fname1 = paste("ZDistance2output_particles",particles,"_pacc",prob,".txt",sep="")
  mydat1 = as.matrix(read.table(fname1, header=T))
  mydat1 = mydat1[1:200,]
  
  fname4 = paste("ZDistance9output_particles",particles,"_pacc",prob,".txt",sep="")
  mydat4 = as.matrix(read.table(fname4, header=T))
  mydat4 = mydat4[1:200,]
  a1 =1+k*4
  b1 = (k+1)*4
  #########RELATIVE BIAS#########
  # Mean and Sd each parameter
  
  
  
  H1a = rbind(apply(mydat1[,a1:b1],2,mean), apply(mydat1[,a1:b1],2,sd))
  # Mean and Sd of avarage
  h1 = apply(mydat1[,a1:b1],1,mean)
  H1b  = t(t(c(mean(h1), sd(h1))))
  H1 = cbind(H1b, H1a)
  ####Relative max
  
  #Relative and parameters
  H4a = rbind(apply(mydat4[,a1:b1],2,mean), apply(mydat4[,a1:b1],2,sd))
  h4 = apply(mydat4[,a1:b1],1,mean)
  H4b  = t(t(c(mean(h4), sd(h4))))
  H4 = cbind(H4b, H4a)
  
  H = rbind(H1,H4)
  
  
  colnames(H) = c("Ave", "alpha", "beta", "delta", "gamma")
  # Mean and Sd each parameter
  
  
  H = as.data.frame(H)
  tab <-xtable(H)
  
  digits(tab) <- rep(3,6)
  
  ##############
  
  
  
  
  print(tab,include.rownames=FALSE)
  
}

texoutfun1(1000,.01,1)# median, 5 to 8, bias

texoutfun2(1000,.01,1)# median, 5 to 8, relative bias

######IQR Range 

texoutfun2(1000,.01,2)

###IQR Covergae Rate

texoutfun2(1000, .01, 3)


#############BOXPLOTS#######
fname1 = paste("ZDistance2output_particles",particles,"_pacc",prob,".txt",sep="")
mydat1 = as.matrix(read.table(fname1, header=T))
mydat1 = mydat1[1:200,]


fname4 = paste("ZDistance9output_particles",particles,"_pacc",prob,".txt",sep="")
mydat4 = as.matrix(read.table(fname4, header=T))
mydat4 = mydat4[1:200,]

mydat1 = mydat1[,9:12]
mydat2 = mydat2[,9:12]
mydat3 = mydat3[,9:12]
mydat4 = mydat4[,9:12]


mydat1 = as.data.frame(mydat1)
mydat2 = as.data.frame(mydat1)
mydat3 = as.data.frame(mydat3)
mydat4 = as.data.frame(mydat4)






###
colnames(mydat2) = c("alpha","beta","delta","gamma")
long_mydat2 = mydat2 %>% gather(Parameter, Value, alpha:gamma)

head(long_mydat2)

p2a <-ggplot(long_mydat2, aes(x=Parameter, y=Value, fill=Parameter)) +
  geom_boxplot()+ ylab("L2 IQR Range")+ theme(
    axis.title.y = element_text(color="blue", size=10))+ theme(legend.position="none")
p2a
####
colnames(mydat3) = c("alpha","beta","delta","gamma")
long_mydat3 = mydat3 %>% gather(Parameter, Value, alpha:gamma)

head(long_mydat3)

p3a <-ggplot(long_mydat3, aes(x=Parameter, y=Value, fill=Parameter)) +
  geom_boxplot()+ ylab("L3 IQR Range") + theme(
    axis.title.y = element_text(color="blue", size=10))+ theme(legend.position="none")
p3a
#######
colnames(mydat1) = c("alpha","beta","delta","gamma")
long_mydat1 = mydat1 %>% gather(Parameter, Value, alpha:gamma)

head(long_mydat1)

p1a <-ggplot(long_mydat1, aes(x=Parameter, y=Value, fill=Parameter)) +
  geom_boxplot()+ ylab("Euclidean Distance")+ theme(
    axis.title.y = element_text(color="blue", size=14), legend.position = "none" )
p1a


colnames(mydat4) = c("alpha","beta","delta","gamma")
long_mydat4 = mydat4 %>% gather(Parameter, Value, alpha:gamma)

head(long_mydat4)
c(expression(alpha),expression(alpha),expression(delta),expression(gamma))


p4a <-ggplot(long_mydat4, aes(x=Parameter, y=Value, fill=Parameter)) +
  geom_boxplot()+ ylab("Proposed Distance")+ theme(
    axis.title.y = element_text(color="blue", size=14))
p4a

plot =  plot_grid(p1a,p4a)
title = ggdraw() + draw_label("Inter Quartile Range", color="blue", fontface ="bold",size=16) 

plot_grid(
  title, plot, ncol = 1, #rel height control the title height
  rel_heights = c(0.1, 1)
)
