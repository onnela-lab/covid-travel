#Row1-2: travelinbound_countryconsider, travelinbound_countryconsider_adjusted

#Row 3-4, lower and upper percentile number unobserved infected enter the country prequarantine
#Row 5-6, lower and upper percentile percent unobserved infected enter the country prequarantine

#Rows 7-8: lower and upper number unobserved infected enter the country post quarantine
#Row 9-10, lower and upper percentile percent unobserved infected enter the country post quarantine

#Row 11-12, lower and upper percentile total active confirmed imported
#Row 13-14 percent of active confirmed imported people


#Row 15-16, lower and upper percentile percent change in new cases
#Row 17-18, lower and upper percentile relative change in new cases

#Row 19-20, lower and upper percentile percent change in active confirmed
#Row 21-22, lower and upper percentile relative change in active confirmed


#Row 23-24, lower and upper percentile reduction in Rt
#Row 25-26, lower and upper percentile relative reduction in Rt
##Row27-28 are the percentage of inbound

library(xtable)
policytime =14
inflation = 10/100

numbersce = 8 


#################
myfunction = function(index){
  matrix1 = matrix(0,28,numbersce)
  t1= rep(0,length(index))
  for(i in 1:length(index)){
    countryconsider = index[i]
    fname = paste('ZCoordinateestimationcountryconsider_',countryconsider,".txt",sep="")
    
    if(file.exists(fname)){
      data = as.matrix(read.table(fname,header=T))
      ##reorderpoliciesforconsistentwithmaintext
      #index = c(1,2,3,4,6,8,5,7)
      #data= data[,index]
      matrix1 = matrix1 +data
      t1[i] = 1
    }#end for checking file.exists
  } #end loop for i
  
  matrix1 = matrix1/sum(t1)
  return(matrix1)
  
}

#########Group 1,2,3,4 with Rt range >1.1, 1<Rt<1.1, 0.9<Rt<1, Rt<.9
index4 = c(1,5)
index3 =c(2,3,6,7)
index2 = c(4,8)
index1 = 1:8
########### 
matrix = myfunction(index1)
H = as.data.frame(matrix)
scale1 = 100
scale2 = 100

for (i in 1:4){
  index = get(paste("index",i,sep=""))
  matrix = myfunction(index)
  H = as.data.frame(matrix)
  tab1index = c(17,18,21,22,15,16,19,20)
  tab2index = c(13,14,5,6,9,10,27,28)
  H1 = H[tab1index,]*scale1
  H2 = H[tab2index,]*scale2
  H2 = round(H2,digits=0)
  round(2.5,digits=0)
  # reportindex = c(1,2,3,4,6,8,5,7)
  # 
  # H1= H1[,reportindex]
  # H2 = H2[,reportindex]
  H1 = as.data.frame(H1)
  H2 = as.data.frame(H2)
 
  tab1 <- xtable(H1)
  tab2 <- xtable(H2)
  print(paste("This is 1st tab group",i))
  digits(tab1) <- rep(0,9)
  #print(tab1)
  print(paste("This is 2nd tab group",i))
  digits(tab2) <- rep(0,9)
  print(tab2)
  ############
  
}

