rm(list=ls())
source('~/Documents/R/updateECLdatabase.R')
setwd('~/output_erai_925')
ECLyear('tracks_1987.dat','1987')
ECLyear('tracks_1987_wtmov0.dat','1987_wtmov0')
ECLyear('tracks_1987_wtmov1.dat','1987_wtmov1')
ECLyear('tracks_1987_pmem0.dat','1987_pmem0')
ECLyear('tracks_1987_pmem1.dat','1987_pmem1')
ECLyear('tracks_1987_c3_0.dat','1987_c3_0')
ECLyear('tracks_1987_c1_1.dat','1987_c1_1')
ECLyear('tracks_1987_c1_3.dat','1987_c1_3')
ECLyear('tracks_1987_c2_5.dat','1987_c2_5')
ECLyear('tracks_1987_c2_20.dat','1987_c2_20')
ECLyear('tracks_1987_c1_18.dat','1987_c1_18')
ECLyear('tracks_1987_new.dat','1987_new')
setwd('~/output_erai_925')
ECLyear('tracks_1987_1.dat','1987_cv1')
ECLyear('tracks_1987_05.dat','1987_cv05')


setwd('~/output_erai')
ECLyear('tracks_1987.dat','1987_surf')

##U/V winds for advstat - 200103 NCEP
rm(list=ls())
setwd('~/Documents/ECLs')
library("RNetCDF")
f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.surf.mean.nc')
time=var.get.nc(f1,'time')
lat=var.get.nc(f1,'lat')
lon=var.get.nc(f1,'lon')
time=as.Date((time/24-2),origin="0001-01-01")
I=which(time=="2001-03-01")
uwnd=var.get.nc(f1,'uwnd',c(1,1,I),c(144,73,1),unpack=T)
uwnd2=colMeans(uwnd)

advstat=matrix(0,71,4)
for(i in 1:71)
{
  advstat[i,1]=lat[73-i]
  advstat[i,2]=uwnd2[73-i]
}
read.table('advstat.U')->advstat2
advstat[,3]=advstat2[,2]
advstat[,4]=advstat[,3]-advstat[,2]

##V2 - all the different levels.
f1=open.nc('~/Documents/Data/NCEP Uwind/uwnd.mon.mean.nc')
uwnd=var.get.nc(f1,'uwnd',c(1,1,1,I),c(144,73,3,1),unpack=T)
advstat3=matrix(0,71,5)
advstat3[,1:2]=advstat[,1:2]
for(i in 1:71)
{
  advstat3[i,3]=mean(uwnd[,73-i,1])
  advstat3[i,4]=mean(uwnd[,73-i,2])
  advstat3[i,5]=mean(uwnd[,73-i,3])
}

plot(advstat3[,1],advstat3[,2],type='l',col='red')