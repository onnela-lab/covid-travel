initials =  list()

for (countryconsider in 1:92){ 
fname1 = paste('Bestinitialandparametersforpredictionchangepointcountryconsider_',countryconsider,".txt",sep="")
if(file.exists(fname1)){
tmp1 =  read.table(fname1, header =T)
initials[[countryconsider]] = tmp1
} else{
initials[[countryconsider]] = NA

}

}

fname = paste('Allinitialsforpredictionchangepoint',".Rdata",sep="")
save(initials, file=fname)

