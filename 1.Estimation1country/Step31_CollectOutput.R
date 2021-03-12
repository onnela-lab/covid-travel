particles = 1000
prob = .005
nrep = 59
h = rep(NA,300)
mymatrix1 = mymatrix2 = mymatrix3 = mymatrix4 = mymatrix5 = mymatrix6 = matrix(NA,300,21)
mymatrix7 = mymatrix8 = mymatrix9 = matrix(NA,300,21)


for (r in 0:nrep){
  tname =paste("Distanceoutput_particles",particles,"_pacc",prob,"_loop",r,".txt",sep="")

    if(file.exists(tname)){
    a = 1 + r*5
    b = 5 +r*5
    h[a:b] = a:b
    tmp = as.matrix(read.table(tname,header=T))
    #Distance 1
    mymatrix1[a:b,1:20] = tmp[1:5,]
    mymatrix1[a:b,21] = a:b
    #Distance 2
    mymatrix2[a:b,1:20] = tmp[6:10,]
    mymatrix2[a:b,21] = a:b
    #Distance 3
    mymatrix3[a:b,1:20] = tmp[11:15,]
    mymatrix3[a:b,21] = a:b
    #Distance 4
    mymatrix4[a:b,1:20] = tmp[16:20,]
    mymatrix4[a:b,21] = a:b
    #Distance 5
    mymatrix5[a:b,1:20] = tmp[21:25,]
    mymatrix5[a:b,21] = a:b
    #Distance 6
    mymatrix6[a:b,1:20] = tmp[26:30,]
    mymatrix6[a:b,21] = a:b
    #Distance 7
    mymatrix7[a:b,1:20] = tmp[31:35,]
    mymatrix7[a:b,21] = a:b
    #Distance 8
    mymatrix8[a:b,1:20] = tmp[36:40,]
    mymatrix8[a:b,21] = a:b
    #Distance 9
    mymatrix9[a:b,1:20] = tmp[41:45,]
    mymatrix9[a:b,21] = a:b
  }
  
}

index = which(is.na(h)==1)


if(length(index)!=0){
  mymatrix1 = mymatrix1[-index,]
  mymatrix2 = mymatrix2[-index,]
  mymatrix3 = mymatrix3[-index,]
  mymatrix4 = mymatrix4[-index,]
  mymatrix5 = mymatrix5[-index,]
  mymatrix6 = mymatrix6[-index,]
  mymatrix7 = mymatrix7[-index,]
  mymatrix8 = mymatrix8[-index,]
  mymatrix9 = mymatrix9[-index,]
  
}

fname1 = paste("ZDistance1output_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix1,file=fname1, row.names = F)


fname2 = paste("ZDistance2output_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix2,file=fname2, row.names = F)


fname3 = paste("ZDistance3output_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix3,file=fname3, row.names = F)


fname4 = paste("ZDistance4output_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix4,file=fname4, row.names = F)


fname5 = paste("ZDistance5output_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix5,file=fname5, row.names = F)


fname6 = paste("ZDistance6output_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix6,file=fname6, row.names = F)


fname7 = paste("ZDistance7output_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix7,file=fname7, row.names = F)


fname8 = paste("ZDistance8output_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix8,file=fname8, row.names = F)


fname9 = paste("ZDistance9output_particles",particles,"_pacc",prob,".txt",sep="")
write.table(mymatrix9,file=fname9, row.names = F)

