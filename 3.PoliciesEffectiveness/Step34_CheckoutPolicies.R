
policytime <- 14
inflation <- 10/100


myfunction =  function(countryconsider)

{
nrep = 350
h =  rep(0,nrep)
h1 = h
######################
for (iteration in 1:nrep){

tname = paste('Policieseffectestimationcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")

if(file.exists(tname)){
    h[iteration] = 1
    }
}
#####################
for (iteration in 1:nrep){
 tname1 = paste('Policieseffectbenchmarkcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")

if(file.exists(tname1)){
    h1[iteration] = 1
    }
}
#############################
index1 = which(h!=0)
index2 = which(h1!=0)
index = intersect(index1,index2)
index= sort(index)
return(index[1:200])  # get only 200 iterations for consistent,outputs made sure more than 200 available
}

##################################################
myfunction1 = function(countryconsider){
mymatrix = matrix(0,26,8) 
mymatrixa =  matrix(0,2,8)
mymatrix1 = mymatrix
mymatrix1a =  mymatrixa
index= myfunction(countryconsider)
################################
for (r in 1:length(index)){
iteration = index[r]
 tname = paste('Policieseffectestimationcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")

    tmp = as.matrix(read.table(tname,header=T))
    tmpa = tmp[1:2,]/tmp[1,1]
    #Distance 1
    mymatrix = mymatrix + tmp
    mymatrixa = mymatrixa + tmpa
  
}
mymatrix = mymatrix/length(index)
mymatrixa = mymatrixa/length(index)
mymatrixe = rbind(mymatrix,mymatrixa)
fname = paste("Zpoliciesestimationcountryconsider_",countryconsider,".txt",sep="")
write.table(mymatrixe, file = fname, row.names=F)



#####################
for (r in 1:length(index)){
iteration = index[r]
 tname1 = paste('Policieseffectbenchmarkcountryconsider_',countryconsider,"_durationpolicy_",policytime,"_inflation_",inflation,"_iteration_",iteration,".txt",sep="")

    tmp = as.matrix(read.table(tname1,header=T))
    tmpa = tmp[1:2,]/tmp[1,1]
    #Distance 1
    mymatrix1 = mymatrix1 + tmp
    mymatrix1a = mymatrix1a + tmpa

}
mymatrix1 = mymatrix1/length(index)
mymatrix1a = mymatrix1a/length(index)
mymatrixb = rbind(mymatrix1,mymatrix1a)
fname1 = paste("Zpoliciesbenchmarkcountryconsider_",countryconsider,".txt",sep="")
write.table(mymatrixb, file = fname1, row.names=F)


}
 
#####################################
for(i in 1:4){
myfunction1(i)}
