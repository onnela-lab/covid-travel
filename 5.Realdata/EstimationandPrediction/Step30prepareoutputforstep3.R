initials3 =  list()

for (countryconsider in 1:92){ 
fname1 = paste('Bestinitialandparametersforstep3_',countryconsider,".txt",sep="")
if(file.exists(fname1)){
tmp1 =  read.table(fname1, header =T)
initials3[[countryconsider]] = tmp1
} else{
initials3[[countryconsider]] = NA

}

}


fname = paste('InitialStep3',".Rdata",sep="")
save(initials3, file=fname)

