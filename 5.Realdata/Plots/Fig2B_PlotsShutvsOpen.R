
##The code is created on March - 11, 2021
##Last edit: March 11 12 p.m.
##Purpose: Report the out put for different scenarios of traveling ######



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
library(ggsci)
tmptmp = proc.time()

load("coviddataJanJune20.Rdata")


load("Allinitialsforpredictionchangepoint.Rdata") 
numbercountries = length(coviddataJanJune20$population)
numberactualcountries = numbercountries - 1
durationprediction =  15
scenarios = 4
#Row 1: Travel
#Row 2 and 3, lower and upper percentile percent change in new cases
#Row 4 and 5, lower and upper percentile relative change in new cases
#Row 6  and 7, lower and upper percentile percent change in active confirmed
#Row 8 and 9, lower and upper percentile relative change in active confirmed

R0 =  rep(0,numberactualcountries)
for(countryconsider in 1:numberactualcountries){
  Rts = initials[[countryconsider]]
  Rt = Rts[2,2]/(Rts[2,3] + Rts[2,6])*Rts[1,1]/sum(Rts[1,])
  R0[countryconsider] = Rt
}

index1 = which(R0<=.9)
index2 = which(R0>.9 & R0<=1.1 )
index3 = which(R0>1.1)
indexall = 1:numberactualcountries

#################

matrix1 = matrix(0,numberactualcountries,scenarios) # matrix to save .975 new cases relative change values of each scenarios


for(countryconsider in 1:numberactualcountries){
  
  fname = paste('CIpandemicchangepointcountryconsider_',countryconsider,"_durationprediction_",durationprediction,".txt",sep="")
 
  if(file.exists(fname)){
    data = as.matrix(read.table(fname,header=T))
    matrix1[countryconsider,] =  data[9,]
    
  }#end for checking file.exists
  } #end loop for i



#################
mydata = as.data.frame(matrix1)
countries = coviddataJanJune20$country[indexall]
mydata$countries = countries
mydata$countrygroup = 1:numberactualcountries
colnames(mydata) = c("traffic2020","fullyopen","borderclose","averagecontrol","country","Group")
mydata$ratio = mydata$fullyopen/mydata$borderclose
mydata$Group[index1] = rep("Group1",length(index1))
mydata$Group[index2] = rep("Group2",length(index2))
mydata$Group[index3] = rep("Group3",length(index3))
####################################################

##################Reverse order figure
######################Change point area###############
sizetext =22
sizetext1 = 16
plt1 = ggplot(mydata, aes(x=borderclose, y= fullyopen, label = country)) +
  geom_point(size=.3) +  geom_text(aes(color=Group),size=3, nudge_x = 0.005, nudge_y = 0.005,check_overlap = F) + geom_abline(intercept = 0, slope = 1)+
  labs( x = "Border closure", y = "Fully open")+
  theme_classic()+
  theme(axis.text=element_text(size=sizetext),
    axis.title.x = element_text( size=sizetext, face="bold"),
        axis.title.y = element_text( size=sizetext, face="bold"), legend.position=c(.9, .22),
        legend.title = element_text( size=sizetext, 
                                     face="bold"),legend.text=element_text(size=sizetext1)
  )+ylim(0,3)+
  scale_color_tron()+
  scale_fill_tron()





plt2= ggplot(mydata1, aes(x=borderclose, y= fullyopen, label = country)) +
  geom_point(size=.3) +  geom_text(aes(color=Group),size=4) + geom_abline(intercept = 0, slope = 1)+
  labs( x = "", y = "")+
  theme_classic()+
  theme(axis.text=element_text(size=sizetext1),
    axis.title.x = element_text( size=sizetext1, face="bold"),
    axis.title.y = element_text( size=sizetext1, face="bold"), legend.position="none",
    legend.title = element_text( size=sizetext1, 
                                 face="bold"),legend.text=element_text(size=sizetext1)
  )+
  scale_color_tron()+
  scale_fill_tron()











plte = plt1+ annotation_custom(ggplotGrob(plt2), xmin = -0.1, xmax = 1.2, 
                               ymin = 1.2, ymax = 3)

fname1 = paste('Comparescenarioschangepoint',".pdf",sep="")
ggsave(fname1,plot=plte,scale=1)



