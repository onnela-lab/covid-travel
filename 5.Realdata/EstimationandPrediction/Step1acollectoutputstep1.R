JHKinitial =  list()

for (countryconsider in 1:92){ 
fname1 = paste('InitialandPosteriorprelim_JHKdata_countryconsider_',countryconsider,".txt",sep="")
tmp1 =  read.table(fname1, header =T)
JHKinitial[[countryconsider]] = tmp1

}


fname = paste('AllPreliminariesposteriorandinitialStep1',".Rdata",sep="")
save(JHKinitial, file=fname)

