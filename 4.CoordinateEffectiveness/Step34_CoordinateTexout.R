
library(xtable)
policytime =14
inflation = 10/100

myfunction = function(countryconsider,scalefactor){

fname1 = paste('ZCoordinateeffectbenchmarkcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,".txt",sep="")
  
data =  read.table(fname1,header=T)
data = data*scalefactor
data = as.matrix(data)
data = data[-c(1,2,3),]
data = as.data.frame(data)
#index = c(2:9,12,13,10,11,16,17) # make change in new cases a head
#data1 = as.matrix(data[index,])
#data1 = as.data.frame(data1)
return(data)
}
#################
myfunctione = function(countryconsider,scalefactor){
  
  fname1 = paste('ZCoordinateeffectestimationcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,".txt",sep="")

  data =  read.table(fname1,header=T)
  data = data*scalefactor
  data = as.matrix(data)
  data = data[-c(1,2,3),]
  data = as.data.frame(data)
  #index = c(2:9,12,13,10,11,16,17) # make change in new cases a head
  #data1 = as.matrix(data[index,])
  #data1 = as.data.frame(data1)
  return(data)
}


#########
myfunction1 = function(index,digitsval,scalefactor){
a = index[1]
b =index[2]
data1 = (myfunction(a, scalefactor) +myfunction(b, scalefactor))/2
values <-xtable(data1)
digits(values) <- rep(digitsval,9)
print(values)
}
#########Group 1,2,3,4 with Rt range >1.1, 1<Rt<1.1, 0.9<Rt<1, Rt<.9
index1 = c(1,5)
index2 =c(2,6)
index3 = c(3,7)
index4 = c(4,8)
######
index= index1
t1 = myfunction(index[1],1)[c(9,10,23,24),]
t2 = myfunction(index[2],1)[c(9,10,23,24),]
e1 =  myfunctione(index[1],1)[c(9,10,23,24),]
e2 =  myfunctione(index[2],1)[c(9,10,23,24),]

t1
e1

t2
e2
###As the order we have here:
#row 1 and 2 are percent travel and expected travelers
#row3 and 4 lower and upper of percent unobserved I travel enter the country
#row 5and 6 lower and upper of percent unobserved travel enter after quarantine
#row 7 and 8 lower and upper percentile total A imported
#row 9 and 10 lower and upper percentile new cases
#row 11 and 12 lower and upper relative change in new active confirmed 
# row 13 and 14 lower and upper reduction relative R

#We return a matrix of 8 columns, each column corresponding to one policy
#15th and 16th row, lower and upper percentile  unobserved infected enter the country prequarantine
#17th and 18th row, lower and upper percentile unobserved infected enter the country post quarantine
#19th and 20th row, lower and upper percentile total active confirmed imported by travel
#21st and 22nd row, lower and upper percentile  change in new active confirmed
#23rd and 24th row, lower and upper percentile  change in new cases change
#25th and 26th row, lower and upper percentile  in Rt at the beginning
#27th and 28th row, lower and upper percentile reduction in Rt at the end



