setwd('~/Documents/Data')
years=seq(1900,2012)
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
save(cool,warm,file='AWAPrain.RData')

rm(list=ls())
setwd('~/Documents/Data')
years=seq(1900,2012)
read.csv('GDIlist.csv')->GDI
GDI=GDI[(GDI[,8]=="E" & GDI[,2]>=3 & GDI[,2]<=5),]

E<-W<-array(0,dim=c(691,886,length(years)))

for(i in 1:length(GDI[,1]))
{
  if(GDI[i,3]>=2000) fname<-paste('rainfall_2000-2008/rainfall-',GDI[i,3],'/r',GDI[i,1],'.txt',sep="") else fname<-paste('rainfall_',GDI[i,4],'0-',GDI[i,4],'9/rainfall-',GDI[i,3],'/r',GDI[i,1],'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->data
  as.matrix(data)->data
  data[data<0]=0
  data<-data[nrow(data):1,]
  E[,,GDI[i,9]]=(E[,,GDI[i,9]]+data)
  rm(data)
}

save(E,file="Erain_MAM.RData")

rm(list=ls())
setwd('G:/daily-rainfall')
read.csv('GDIlist.csv')->GDI
GDI=GDI[(GDI[,8]=="W" & GDI[,2]>=3 & GDI[,2]<=5),]

W<-array(0,dim=c(691,886,109))

for(i in 1:length(GDI[,1]))
{
  if(GDI[i,3]>=2000) fname<-paste('rainfall_2000-2008/rainfall-',GDI[i,3],'/r',GDI[i,1],'.txt',sep="") else fname<-paste('rainfall_',GDI[i,4],'0-',GDI[i,4],'9/rainfall-',GDI[i,3],'/r',GDI[i,1],'.txt',sep="")
  read.table(fname, sep="",skip=6,nrows=691)->data
  as.matrix(data)->data
  data[data<0]=0
  data<-data[nrow(data):1,]
  W[,,GDI[i,9]]=(W[,,GDI[i,9]]+data)
  rm(data)
}
save(W,file="Wrain_MAM.RData")



##Load count of missing values by year
years=seq(1950,2009)
missing=array(0,dim=c(691,886,length(years)))
for(i in 1:length(years))
{
  print(years[i])
  decade=floor(years[i]/10)
  dir=paste('/media/Seagate Expansion Drive/Data/daily rainfall/rainfall_',
            decade,'0-',decade,'9/rainfall-',years[i],"/",sep="")
  fnames=list.files(path=dir,pattern="*.txt")
  for(j in 1:length(fnames))
{
    fin=paste(dir,fnames[j],sep="")
    read.table(fin, sep="",skip=6,nrows=691)->data
    as.matrix(data)->data
    data=data[nrow(data):1,]
    test=matrix(0,691,886)
    I=which(data<0)
    test[I]=1
    missing[,,i]=missing[,,i]+test
  }
}