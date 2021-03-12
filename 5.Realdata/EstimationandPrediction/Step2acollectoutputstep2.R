
##########################################################
JHKPosIniStep2 =  list()
thresholddeath = 0.05
for (countryconsider in 1:92){ 
fname1 = paste('GlobalPosteriorandInitial_JHKdata_countryconsider_',countryconsider,'_thresholddeath',thresholddeath,".txt",sep="")

if(file.exists(fname1)){
tmp1 =  read.table(fname1, header =T)
JHKPosIniStep2[[countryconsider]] = tmp1
} else{
JHKPosIniStep2[[countryconsider]] =NA
}
}

fname = paste('AllposteriorandinitialforPredictionStep2thresholddeath1',".Rdata",sep="")
save(JHKPosIniStep2, file=fname)

