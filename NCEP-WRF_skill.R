rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/CMIP-WRF-ECLs/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)
wrf=c("R1","R2","R3")
namelist=rep("aaa",12)

events<-fixes<-list()

n=1
for(proj in c(100,240))
  for(r in 1:3)
  {
    if(proj==100) cv=1 else cv=1.35
    events[[n]]=read.csv(paste("outputUM/proj",proj,"/outputUM_ncep_WRF",wrf[r],"_50_rad2cv1/ECLevents_umelb_ncep_wrf",wrf[r],"_proj",proj,"_rad2cv",cv,"_9009.csv",sep=""))
    events[[n]]$Year=floor(events[[n]]$Date1/10000)
    events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
    
    fixes[[n]]=read.csv(paste("outputUM/proj",proj,"/outputUM_ncep_WRF",wrf[r],"_50_rad2cv1/ECLfixes_umelb_ncep_wrf",wrf[r],"_proj",proj,"_rad2cv",cv,"_9009.csv",sep=""))
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    n=n+1
  }

proj=100
eventsN=read.csv(paste("outputUM/proj",proj,"/outputUM_ncep_rad2cv1/ECLevents_umelb_ncep_proj",proj,"_rad2cv1_9009.csv",sep=""))
eventsN$Year=floor(eventsN$Date1/10000)
eventsN$Month=floor(eventsN$Date1/100)%%100

fixesN=read.csv(paste("outputUM/proj",proj,"/outputUM_ncep_rad2cv1/ECLfixes_umelb_ncep_proj",proj,"_rad2cv1_9009.csv",sep=""))
fixesN$Year=floor(fixesN$Date/10000)
fixesN$Month=floor(fixesN$Date/100)%%100
fixesN$Date2=as.POSIXct(paste(as.character(fixesN$Date),substr(fixesN$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

eventsE=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv",stringsAsFactors = F)
eventsE$Year=floor(eventsE$Date1/10000)
eventsE$Month=floor(eventsE$Date1/100)%%100
eventsE=eventsE[eventsE$Year>=1990 & eventsE$Year<=2009,]

fixesE=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv",stringsAsFactors = F)
fixesE$Year=floor(fixesE$Date/10000)
fixesE$Month=floor(fixesE$Date/100)%%100
fixesE$Date2=as.POSIXct(paste(as.character(fixesE$Date),substr(fixesE$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

###### Step 1 - Count by year/month

years=cbind(1990:2009,matrix(0,20,8))
colnames(years)<-c("Year","NCEP p100","ERAI p100","R1 p100","R2 p100","R3 p100","R1 p240","R2 p240","R3 p240")
months=cbind(1:12,matrix(0,12,8))
colnames(months)<-c("Month","NCEP p100","ERAI p100","R1 p100","R2 p100","R3 p100","R1 p240","R2 p240","R3 p240")

for(i in 1:20)
{
  I=which(eventsN$Year==years[i,1])
  years[i,2]=length(I)
  I=which(eventsE$Year==years[i,1])
  years[i,3]=length(I)
  for(j in 1:6)
  {
    I=which(events[[j]]$Year==years[i,1])
    years[i,j+3]=length(I)
  }
}

apply(years,2,mean)
for(i in 1:6) print(cor(years[,2],years[,i+2]))

for(i in 1:12)
{
  I=which(eventsN$Month==months[i,1])
  months[i,2]=length(I)
  I=which(eventsE$Month==months[i,1])
  months[i,3]=length(I)
  for(j in 1:6)
  {
    I=which(events[[j]]$Month==months[i,1])
    months[i,j+3]=length(I)
  }
}

months2=months
for(i in 1:12) months2[i,2:9]=100*months[i,2:9]/apply(months[,2:9],2,sum)

apply(months2[5:10,],2,sum)

plot(months2[,1],months2[,2],col=1,lwd=3,type="l",xlab="Month",ylab="% of ECLs",ylim=c(0,max(months2)))
for(i in 2:8) lines(months2[,1],months2[,i+1],col=i,lwd=3)


##########

yearsA=years
monthsA=months

for(i in 4:6) events[[i]]$CV2=events[[i]]$CVmax

thresh=22
thresh2=rep(0,6)

for(i in 1:6)
{
  data=events[[i]]
  b=order(data$CV2,decreasing=T)
  thresh2[i]=data$CV2[b[20*thresh]]
  if(is.na(thresh2[i])) thresh2[i]=min(data$CV2,na.rm=T)
  
  for(y in 1:20) yearsA[y,i+3]=length(which(data$Year==years[y,1] & data$CV2>=thresh2[i]))
  for(m in 1:12) monthsA[m,i+3]=length(which(data$Month==m & data$CV2>=thresh2[i]))
}


corrs=matrix(0,8,8)
for(i in 1:8)
  for(j in 1:8)
    corrs[i,j]=cor(yearsA[,i+1],yearsA[,j+1])
corrs[corrs==1]=NaN

months2=monthsA
for(i in 1:12) months2[i,2:9]=100*monthsA[i,2:9]/apply(monthsA[,2:9],2,sum)

apply(months2[5:10,],2,sum)

clist=c("black","grey","red","blue","purple")
plot(months2[,1],months2[,2],col=1,lwd=3,type="l",xlab="Month",ylab="% of ECLs",ylim=c(0,max(months2)))
for(i in 2:5) lines(months2[,1],months2[,i+1],col=clist[i],lwd=3)

###################### ECL matching?


matches=matrix(0,length(eventsN$ID),7)
for(i in 1:length(eventsN$ID))
{
  tmp=fixesN[(fixesN$ID==eventsN$ID[i] & fixesN$Location==1),]
  rn=range(tmp$Date2)
  
  I=which(fixesE$Date2<=rn[2]+(60*60*6) & fixesE$Date2>=rn[1]-(60*60*6) & fixesE$Location==1)
  if(length(I)>0) matches[i,1]=1
  
  for(j in 1:6)
  {
    J=which(events[[j]]$CV2>=thresh2[j])
    
    I=which(fixes[[j]]$Date2<=rn[2]+(60*60*6) & fixes[[j]]$Date2>=rn[1]-(60*60*6) & 
              fixes[[j]]$Location==1 & fixes[[j]]$ID%in%events[[j]]$ID[J])
    if(length(I)>0) matches[i,1+j]=1
  }
  
}

apply(matches,2,mean)


### Same version, different projection

for(n in 1:3)
{
  K=which(events[[n]]$CV2>=thresh2[n])
  matches=rep(0,length(K))
for(i in 1:length(K))
{
 tmp=fixes[[n]][(fixes[[n]]$ID==events[[n]]$ID[K[i]] & fixes[[n]]$Location==1),]
 rn=range(tmp$Date2)
  
  j=n+3
  J=which(events[[j]]$CV2>=thresh2[j])
  I=which(fixes[[j]]$Date2<=rn[2]+(60*60*6) & fixes[[j]]$Date2>=rn[1]-(60*60*6) & 
              fixes[[j]]$Location==1 & fixes[[j]]$ID%in%events[[j]]$ID[J])
    if(length(I)>0) matches[i]=1
}
print(mean(matches))
}



########## Locations

lat=seq(-40,-24,2)
lon=seq(148,160,2)

loc=array(0,c(9,7,8))

for(y in 1:length(lat))
  for(x in 1:length(lon))
  {
    I=which(fixesN$Lon>=lon[x]-1 & fixesN$Lon<lon[x]+1 & fixesN$Lat>=lat[y]-1 & fixesN$Lat<lat[y]+1 & fixesN$Location==1)
    loc[y,x,1]=length(I)
    I=which(fixesE$Lon>=lon[x]-1 & fixesE$Lon<lon[x]+1 & fixesE$Lat>=lat[y]-1 & fixesE$Lat<lat[y]+1 & fixesE$Location==1)
    loc[y,x,2]=length(I)
    
    for(n in 1:6)
    {
      J=which(events[[n]]$CV2>=thresh2[n])
      I=which(fixes[[n]]$Lon>=lon[x]-1 & fixes[[n]]$Lon<lon[x]+1 & fixes[[n]]$Lat>=lat[y]-1 & fixes[[n]]$Lat<lat[y]+1 & 
                fixes[[n]]$Location==1 & fixes[[n]]$ID%in%events[[n]]$ID[J] & fixes[[n]]$CV>=thresh2[n])
      loc[y,x,n+2]=length(I)
    }
    
  }

loc=loc/20

names=c("NCEP_p100","ERAI_p100","NCEP-R1_p100","NCEP-R2_p100","NCEP-R3_p100","NCEP-R1_p240","NCEP-R2_p240","NCEP-R3_p240")

library(maps)
ColorBar <- function(brks,cols,labels=NA)
{
  par(mar = c(3, 1, 3, 4), mgp = c(1, 1, 0), las = 1, cex = 1.2)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  if(is.na(labels[1])) labels=brks[seq(2, length(brks)-1)]
  
  axis(4, at = seq(1.5, length(brks) - 1.5), tick = TRUE, 
       labels = labels)
}
cm2=gray(seq(1,0.1,-0.15))
bb2=c(-0.5,0,1,2,3,4,5,100)

for(i in 1:8)
{
pdf(file=paste("outputUM/ECL_locations_",names[i],"_top22_d01.pdf",sep=""),width=6,height=4.5,pointsize=12)
layout(cbind(1,2),c(1,0.3))
par(mar=c(3,3,3,2))
image(lon,lat,t(loc[,,i]),xlab="",ylab="",breaks=bb2,col=cm2,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=names[i],cex.axis=1.5,cex.main=1.5)
map(xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),add=T,lwd=2)
ColorBar(bb2,cm2)
dev.off()
}

########## Whole CORDEX domain

fixesALL<-list()
fixesALL[[1]]=read.csv(paste("outputUM/proj",proj,"/outputUM_ncep_rad2cv1/ECLfixes_umelb_ncep_rad2cv1_9009_ALL.csv",sep=""))
fixesALL[[2]]=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_rad2cv1_9009_ALL.csv",stringsAsFactors = F)
n=3
for(proj in c(100,240))
  for(r in 1:3)
  {
    tmp=read.csv(paste("outputUM/proj",proj,"/outputUM_ncep_WRF",wrf[r],"_50_rad2cv1/ECLfixes_umelb_ncep_wrf",wrf[r],"_rad2cv1_9009_ALL.csv",sep=""))
    fixesALL[[n]]=tmp[tmp$CV>=thresh2[n-2],]
    n=n+1
  }

lat=seq(-50,0,2.5)
lon=seq(100,180,2.5)

loc2=array(0,c(21,33,8))

for(y in 1:length(lat))
  for(x in 1:length(lon))
    for(n in 1:8)
    {
      I=which(fixesALL[[n]]$Lon>=lon[x]-1.25 & fixesALL[[n]]$Lon<lon[x]+1.25 & fixesALL[[n]]$Lat>=lat[y]-1.25 & fixesALL[[n]]$Lat<lat[y]+1.25)
      loc2[y,x,n]=length(I)
    }

loc2=loc2/20

names=c("NCEP_p100","ERAI_p100","NCEP-R1_p100","NCEP-R2_p100","NCEP-R3_p100","NCEP-R1_p240","NCEP-R2_p240","NCEP-R3_p240")

for(i in 1:8)
{
  pdf(file=paste("outputUM/ECL_locations_CORDEX_",names[i],"_top22_d01_v2.pdf",sep=""),width=8,height=5,pointsize=12)
  layout(cbind(1,2),c(1,0.2))
  par(mar=c(3,3,3,2))
  image(lon,lat,t(loc2[,,i]),xlab="",ylab="",breaks=bb2,col=cm2,zlim=c(-Inf,Inf),main=names[i],cex.axis=1.5,cex.main=1.5,
        xlim=c(110,172.5),ylim=c(-45,-10))
  map(add=T,lwd=2)
  ColorBar(bb2,cm2)
  dev.off()
}

pdf(file=paste("outputUM/ECL_locations_CORDEX_NCEP-WRF_proj100_top22_d01_v2.pdf",sep=""),width=8,height=5,pointsize=12)
layout(cbind(1,2),c(1,0.2))
par(mar=c(3,3,3,2))
image(lon,lat,t(apply(loc2[,,3:5],c(1,2),mean)),xlab="",ylab="",breaks=bb2,
      col=cm2,zlim=c(-Inf,Inf),main="NCEP-WRF",cex.axis=1.5,cex.main=1.5,
      xlim=c(110,172.5),ylim=c(-45,-10))
map(add=T,lwd=2)
ColorBar(bb2,cm2)
dev.off()

####### Compare to merra

eventsM=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_merraN_rad2cv1_proj100_diff2/ECLevents_umelb_merra_rad2cv1_9009.csv")
eventsM=eventsM[eventsM$Date1>=19900000 & eventsM$Date2<=20100000,]
fixesM=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_merraN_rad2cv1_proj100_diff2/ECLfixes_umelb_merra_rad2cv1_9009_ALL.csv",stringsAsFactors = F)
fixesM=fixesM[fixesM$Date>=19900000 & fixesM$Date<=20100000,]

locM=array(0,c(21,33))

for(y in 1:length(lat))
  for(x in 1:length(lon))
  {
    I=which(fixesM$Lon>=lon[x]-1.25 & fixesM$Lon<lon[x]+1.25 & fixesM$Lat>=lat[y]-1.25 & fixesM$Lat<lat[y]+1.25)
      locM[y,x]=length(I)
}
locM=locM/20

pdf(file=paste("outputUM/ECL_locations_CORDEX_MERRA_proj100_d01_v2.pdf",sep=""),width=8,height=5,pointsize=12)
layout(cbind(1,2),c(1,0.2))
par(mar=c(3,3,3,2))
image(lon,lat,t(locM),xlab="",ylab="",breaks=bb2,
      col=cm2,zlim=c(-Inf,Inf),main="MERRA",cex.axis=1.5,cex.main=1.5,
      xlim=c(110,172.5),ylim=c(-45,-10))
map(add=T,lwd=2)
ColorBar(bb2,cm2)
dev.off()


###############
##########

## Histograms of ECL intensity, duration, etc.

fixesM=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_merraN_rad2cv1_proj100_diff2/ECLfixes_umelb_merra_rad2cv1_9009.csv",stringsAsFactors = F)
fixesM=fixesM[fixesM$Date>=19900000 & fixesM$Date<=20100000,]

cvthresh=c(seq(1,4,0.25),Inf)
cvhist=matrix(0,13,6)
rownames(cvhist)=cvthresh[1:13]
colnames(cvhist)=c("NCEP","ERAI","MERRA","R1","R2","R3")

for(i in 1:length(cvhist[,1]))
{
  cvhist[i,1]=length(which(eventsN$CV2>=cvthresh[i] & eventsN$CV2<cvthresh[i+1]))
  cvhist[i,2]=length(which(eventsE$CV2>=cvthresh[i] & eventsE$CV2<cvthresh[i+1]))
  cvhist[i,3]=length(which(eventsM$CV2>=cvthresh[i] & eventsM$CV2<cvthresh[i+1]))
  for(n in 1:3) cvhist[i,n+3]=length(which(events[[n]]$CV2>=cvthresh[i] & events[[n]]$CV2<cvthresh[i+1]))
}

clist=c("black","grey","blue","red","darkgreen","purple")
plot(cvthresh[1:13],cvhist[,1],lwd=2,type="l",
     xlab="Intensity",ylab="Number of events",ylim=c(0,250))
for(i in 2:6) lines(cvthresh[1:13],cvhist[,i],lwd=2,col=clist[i])

lenthresh=c(1,5,9,13,17,21,Inf)
lenhist=matrix(0,6,6)
rownames(lenhist)=lenthresh[1:6]
colnames(lenhist)=c("NCEP","ERAI","MERRA","R1","R2","R3")

for(i in 1:length(lenhist[,1]))
{
  lenhist[i,1]=length(which(eventsN$Length2>=lenthresh[i] & eventsN$Length2<lenthresh[i+1]))
  lenhist[i,2]=length(which(eventsE$Length2>=lenthresh[i] & eventsE$Length2<lenthresh[i+1]))
  lenhist[i,3]=length(which(eventsM$Length2>=lenthresh[i] & eventsM$Length2<lenthresh[i+1]))
  for(n in 1:3) lenhist[i,n+3]=length(which(events[[n]]$Length2>=lenthresh[i] & events[[n]]$Length2<lenthresh[i+1] & events[[n]]$CV2>=thresh2[n]))
}

###### That's odd - is it because these events move slower on average?

fixesN$Movement=NaN

library(sp)

for(i in 2:length(fixesN[,1]))
 if(fixesN$ID[i]==fixesN$ID[i-1])
   fixesN$Movement[i]=spDistsN1(as.matrix(cbind(fixesN$Lon[i],fixesN$Lat[i])),c(fixesN$Lon[i-1],fixesN$Lat[i-1]),longlat=T)

eventsN$Move2<-eventsN$Move<-rep(NaN,length(eventsN$ID))
for(i in 1:length(eventsN$ID))
{
  I=which(fixesN$ID==eventsN$ID[i])
  eventsN$Move[i]=mean(fixesN$Movement[I],na.rm=T)
  I=which(fixesN$ID==eventsN$ID[i] & fixesN$Location==1)
  eventsN$Move2[i]=mean(fixesN$Movement[I],na.rm=T)
}

fixesE$Movement=NaN
for(i in 2:length(fixesE[,1]))
  if(fixesE$ID[i]==fixesE$ID[i-1])
    fixesE$Movement[i]=spDistsN1(as.matrix(cbind(fixesE$Lon[i],fixesE$Lat[i])),c(fixesE$Lon[i-1],fixesE$Lat[i-1]),longlat=T)

eventsE$Move2<-eventsE$Move<-rep(NaN,length(eventsE$ID))
for(i in 1:length(eventsE$ID))
{
  I=which(fixesE$ID==eventsE$ID[i])
  eventsE$Move[i]=mean(fixesE$Movement[I],na.rm=T)
  I=which(fixesE$ID==eventsE$ID[i] & fixesE$Location==1)
  eventsE$Move2[i]=mean(fixesE$Movement[I],na.rm=T)
}

for(n in 1:6)
{
  fixes[[n]]$Movement=NaN
  for(i in 2:length(fixes[[n]][,1]))
    if(fixes[[n]]$ID[i]==fixes[[n]]$ID[i-1])
      fixes[[n]]$Movement[i]=spDistsN1(as.matrix(cbind(fixes[[n]]$Lon[i],fixes[[n]]$Lat[i])),c(fixes[[n]]$Lon[i-1],fixes[[n]]$Lat[i-1]),longlat=T)
    
    events[[n]]$Move2<-events[[n]]$Move<-rep(NaN,length(events[[n]]$ID))
    for(i in 1:length(events[[n]]$ID))
    {
      I=which(fixes[[n]]$ID==events[[n]]$ID[i])
      events[[n]]$Move[i]=mean(fixes[[n]]$Movement[I],na.rm=T)
      I=which(fixes[[n]]$ID==events[[n]]$ID[i] & fixes[[n]]$Location==1)
      events[[n]]$Move2[i]=mean(fixes[[n]]$Movement[I],na.rm=T)
    }
}

mthresh=c(seq(0,1200,60),Inf)
mhist=matrix(0,21,5)
rownames(mhist)=mthresh[1:21]
colnames(mhist)=c("NCEP","ERAI","R1","R2","R3")

for(i in 1:length(mhist[,1]))
{
  mhist[i,1]=length(which(eventsN$Move2>=mthresh[i] & eventsN$Move2<mthresh[i+1]))
  mhist[i,2]=length(which(eventsE$Move2>=mthresh[i] & eventsE$Move2<mthresh[i+1]))
  for(n in 1:3) mhist[i,n+2]=length(which(events[[n]]$Move2>=mthresh[i] & events[[n]]$Move2<mthresh[i+1] & events[[n]]$CV2>=thresh2[n]))
}

clist=c("black","grey","red","darkgreen","purple")
plot(mthresh[1:21],mhist[,1],lwd=2,type="l",
     xlab="Intensity",ylab="Number of events",ylim=c(0,250))
for(i in 2:5) lines(mthresh[1:21],mhist[,i],lwd=2,col=clist[i])




