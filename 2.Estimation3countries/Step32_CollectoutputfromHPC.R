
args <- commandArgs(trailingOnly = TRUE)
particles  <- as.integer(args[1])
prob1 <- as.integer(args[2])

prob = prob1/1000
nrep = 300
h =h1 =  rep(NA,nrep)
mymatrix1 = matrix(NA,nrep,450) 

for (r in 1:nrep){
 
  tname = paste("Travelestimateoutput_particles",particles,"_pacc",prob,"_loop",r,".txt",sep="")


    if(file.exists(tname)){
    tmp = as.matrix(read.table(tname,header=T))
    #Distance 1
    mymatrix1[r,] = tmp
    h[r] = r
}
}



index = which(is.na(h)==1)


if(length(index)!=0){
  mymatrix1 = mymatrix1[-index,]
 
}

fname1 = paste("ZC3countriesDistanceoutput_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix1,file=fname1, row.names = F)





