###Just shoving this in here so don't do anywhere else
setwd('~/Documents/Data')
years=seq(1900,2012)
ESBrain=matrix(0,length(years),12)
for(i in 1:length(years))
  for(j in 1:12)
  {
    if(j<10) fin=paste('/media/Seagate Expansion Drive/monthly rainfall/',years[i],'0',j,'.grid',sep="")
    else fin=paste('/media/Seagate Expansion Drive/monthly rainfall/',years[i],j,'.grid',sep="")
    read.table(fin, sep="",skip=6,nrows=691)->data
    as.matrix(data)->data
    data[data<0]=NaN
    data<-data[nrow(data):1,]*escci$mask.escci
    ESBrain[i,j]=mean(data,na.rm=T)
  }
write.csv(ESBrain,file='ESBrain.csv')

rm(list=ls())
setwd('~/Documents/Data')
read.csv('ESBrain.csv',sep=';',header=T)->ESBrain
R<-R2<-matrix(NaN,113,3)
R[,1]<-R2[,1]<-ESBrain[,1]
R[,2]<-rowSums(ESBrain[,7:11])
R[1:112,3]<-rowSums(cbind(ESBrain[1:112,12:13],ESBrain[2:113,2:4]))
R2[,2]=R[,2]/sd(R[,2])
R2[,3]=R[,3]/sd(R[,3],na.rm=T)

read.table('~/Documents/Timeseries/n34.txt',header=T)->n34
E<-E2<-E3<-matrix(NaN,113,3)
E[,1]<-E2[,1]<-E3[,1]<-n34[,1]
E[,2]<-rowSums(n34[,7:11])
E[1:112,3]<-rowSums(cbind(n34[1:112,12:13],n34[2:113,2:4]))
E2[,2]=E[,2]/sd(E[,2])
E2[,3]=E[,3]/sd(E[,3],na.rm=T)

read.table('~/Documents/Timeseries/dmi.txt',header=T)->dmi
D<-D2<-matrix(NaN,113,3)
D[,1]<-D2[,1]<-dmi[,1]
D[,2]<-rowSums(dmi[,7:11])
D[1:112,3]<-rowSums(cbind(dmi[1:112,12:13],dmi[2:113,2:4]))
D2[,2]=D[,2]/sd(D[,2])
D2[,3]=D[,3]/sd(D[,3],na.rm=T)
E3[,2]=resid(lm(E2[,2] ~ D2[,2]))
E3[1:112,3]=resid(lm(E2[1:112,3] ~ D2[1:112,3]))

###Okay, now the combonents
corrs=matrix(0,3,2)
for(i in 2:3)
{
  corrs[1,i-1]=cor(R2[83:113,i],E2[83:113,i],use='complete.obs')
  corrs[2,i-1]=cor(R2[83:113,i],D2[83:113,i],use='complete.obs')
  corrs[3,i-1]=cor(R2[83:113,i],E3[83:113,i],use='complete.obs')
}

read.table('~/Documents/Timeseries/dmiJ.txt',header=T)->dmiJ
DJ<-DJ2<-EJ<-DJa<-matrix(NaN,31,3)
DJ[,1]<-DJ2[,1]<-EJ[,1]<-dmiJ[,1]
DJ[,2]<-rowSums(dmiJ[,7:11])
DJ[1:30,3]<-rowSums(cbind(dmiJ[1:30,12:13],dmiJ[2:31,2:4]))
DJ2[,2]=DJ[,2]/sd(DJ[,2])
DJ2[,3]=DJ[,3]/sd(DJ[,3],na.rm=T)
EJ[,2]=resid(lm(E2[83:113,2] ~ DJ2[,2]))
EJ[1:30,3]=resid(lm(E2[83:112,3] ~ DJ2[1:30,3]))
DJa[,2]=resid(lm(DJ2[,2] ~ E2[83:113,2]))
DJa[1:30,3]=resid(lm(DJ2[1:30,3] ~ E2[83:112,3]))
corrsJ=matrix(0,3,2)
for(i in 2:3)
{
  corrsJ[1,i-1]=cor(R2[83:113,i],E2[83:113,i],use='complete.obs',method='spearman')
  corrsJ[2,i-1]=cor(R2[83:113,i],DJ2[,i],use='complete.obs',method='spearman')
  corrsJ[3,i-1]=cor(R2[83:113,i],EJ[,i],use='complete.obs',method='spearman')
}

rm(list=ls())

setwd('~/Documents/Data')
years=seq(1982,2012)
cool=array(0,dim=c(691,886,length(years)))
for(i in 1:length(years))
  for(j in 5:10)
  {
    if(j<10) fin=paste('/media/Seagate Expansion Drive/monthly rainfall/',years[i],'0',j,'.grid',sep="")
    else fin=paste('/media/Seagate Expansion Drive/monthly rainfall/',years[i],j,'.grid',sep="")
    read.table(fin, sep="",skip=6,nrows=691)->data
    as.matrix(data)->data
    data[data<0]=NaN
    cool[,,i]=cool[,,i]+data[nrow(data):1,]
  }
warm=array(0,dim=c(691,886,length(years)-1))
for(i in 1:(length(years)-1))
{
  for(j in 11:12)
  {
    fin=paste('/media/Seagate Expansion Drive/monthly rainfall/',years[i],j,'.grid',sep="")
    read.table(fin, sep="",skip=6,nrows=691)->data
    as.matrix(data)->data
    data[data<0]=NaN
    warm[,,i]=warm[,,i]+data[nrow(data):1,]
  }
  for(j in 1:3)
  {
    fin=paste('/media/Seagate Expansion Drive/monthly rainfall/',years[i+1],'0',j,'.grid',sep="")
    read.table(fin, sep="",skip=6,nrows=691)->data
    as.matrix(data)->data
    data[data<0]=NaN
    warm[,,i]=warm[,,i]+data[nrow(data):1,]
  }
}

cool2=cool
warm2=warm
for(i in 1:691)
  for(j in 1:886)
  {
    cool2[i,j,]=cool[i,j,]/sd(cool[i,j,],na.rm=T)
    warm2[i,j,]=warm[i,j,]/sd(warm[i,j,],na.rm=T)
  }

corNw<-corNc<-corN2w<-corN2c<-corDw<-corDc<-matrix(0,691,886)
for(i in 1:691)
  for(j in 1:886)
  {
    corNw[i,j]=cor(warm2[i,j,],E2[83:112,3])
    corNc[i,j]=cor(cool2[i,j,],E2[83:113,2])
    corN2w[i,j]=cor(warm2[i,j,],EJ[1:30,3])
    corN2c[i,j]=cor(cool2[i,j,],EJ[,2])
    corDw[i,j]=cor(warm2[i,j,],DJ2[1:30,3])
    corDc[i,j]=cor(cool2[i,j,],DJ2[,2])    
  }

save(cool,warm,cool2,warm2,file='GDI_rev2_rain.RData')

source('~/Documents/R/corr.plot.R')
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
corr.plot(Useful$x,Useful$y,t(corNw*Useful$mask),mask,cm,'CorrW_rev2_N34_8011')
corr.plot(Useful$x,Useful$y,t(corN2w*Useful$mask),mask,cm,'CorrW_rev2_N34noDMI_8011')
corr.plot(Useful$x,Useful$y,t(corDw*Useful$mask),mask,cm,'CorrW_rev2_DMI_8011')
corr.plot(Useful$x,Useful$y,t(corNc*Useful$mask),mask,cm,'CorrC_rev2_N34_8012')
corr.plot(Useful$x,Useful$y,t(corN2c*Useful$mask),mask,cm,'CorrC_rev2_N34noDMI_8012')
corr.plot(Useful$x,Useful$y,t(corDc*Useful$mask),mask,cm,'CorrC_rev2_DMI_8012')
