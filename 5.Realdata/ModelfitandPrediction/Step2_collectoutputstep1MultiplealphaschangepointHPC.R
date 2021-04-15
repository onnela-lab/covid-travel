initialmat = matrix(0,93,6)
JHKposterior = list()
for (countryconsider in 1:92){ 

fname = paste('Posteriorandinitialfromstep1withmultiplealphachangepointcountry',countryconsider,".Rdata",sep="")
if(file.exists(fname)){
load(fname)
tmp1 =  output$initial
tmp2 = output$posteriorstep1
initialmat[countryconsider,] = tmp1
JHKposterior[[countryconsider]] = tmp2
}
}
#####Create a dummy initial for country 93



load("/n/holyscratch01/onnela_lab/ThienLe/Realdata/Modelfit/coviddataJanJune20.Rdata")

P = coviddataJanJune20$population #population all countries
worldpopulation = 7.8*10^9  #WORLD POPULATION
countrydummypop  = worldpopulation - sum(P)

initialmat[93,] = c(countrydummypop, rep(0,5))

fname1a = paste('AllPreliminariesposteriofromstep1changepoint',".Rdata",sep="")
save(JHKposterior, file=fname1a)
initialmat = as.matrix(initialmat)
fname1b = paste('Allinitialfromstep1changepoint',".txt",sep="")
write.table(initialmat, file=fname1b, row.names=F)

