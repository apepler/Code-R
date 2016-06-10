library(fields)
levels=seq(0,2,0.2)
cols=gray(seq(1,0.1,-0.1))
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0
mldb<-read.csv('~/Documents/ECLs/MLDB.csv',sep=";")
a=stack(mldb[,1:2])
ecl=unique(a[,1])
rain<-read.csv('~/Documents/ECLs/Impacts/morwenna.csv')
rain2=rain[which(rain[,6]==1),] 

##Extracts all fixes in Australia with CV>=0.25
ECLdatabase<-function(yearS,yearE,output)
{  
  ##Files with ECL track data for required years 
  year<-as.character(seq(yearS,yearE))
  fname1=paste('tracks_',year,'.dat',sep="")
  
  ##collates all the annual tracks into one file
  
  read.table(fname1[1], sep="",skip=1)->data
  I<-which(data[,7]<(-10) & data[,7]>(-50) & data[,6]>110 & data[,6]<180 & data[,9]>=0.25)
  fixes<-data[I,]
  
  for(i in 2:length(year)) ##Only takes broader Aus box to make faster
  {
    read.table(fname1[i], sep="",skip=1)->data
    I<-which(data[,7]<(-10) & data[,7]>(-50) & data[,6]>110 & data[,6]<180 & data[,9]>=0.25)
    fixes<-rbind(fixes,data[I,])
  }
  rm(data)
  
  Fixes<-data.frame(fixes[,3:9])
  names(Fixes)=c('Date','Time','Open','Lon','Lat','MSLP','CV') 
  outF=paste('ECLfixes_aus_',output,'.csv',sep="")
  write.csv(Fixes,outF)
}

##Extracts all fixes in Australia where event was > 2 fixes
ECLdatabase2<-function(yearS,yearE,output)
{
  year<-as.character(seq(yearS,yearE))
  fname1=paste('tracks_',year,'.dat',sep="")
  
  ##collates all the annual tracks into one file
  
  read.table(fname1[1], sep="",skip=1)->fixes
  I<-which(fixes[,7]<0)
  fixes<-fixes[I,]
  
  for(i in 2:length(year))
  {
    read.table(fname1[i], sep="",skip=1)->data
    I<-which(data[,7]<0)
    fixes<-rbind(fixes,data[I,])
  }
  rm(data)
  
  #Unique event numbers (starting with year)
  fixes[,1]<-fixes[,1]+floor(fixes[,3]/10000)*1000000
  
  ##Identify if within the ECL region
  fixes[,11]<-0
  I<-which(fixes[,7]<(-10) & fixes[,7]>(-50) & fixes[,6]>110 & fixes[,6]<180)
  fixes[I,11]=1
  
  ##Collates information for events
  x<-rle(fixes[,1])
  events<-cbind(x$values,x$lengths,matrix(data=0,nrow=length(x$values),ncol=8))
  
  ##Finds if a multi-fix event
  I<-which(events[,2]>1)
  y<-(fixes[,1] %in% events[I,1])
  
  ##Select only fixes from a multi-fix event, in aus region, CV>0.25
  J<-which(y==T & fixes[,11]==1 & fixes[,9]>=0.25)  
  Fixes<-data.frame(fixes[J,3:9])
  names(Fixes)=c('Date','Time','Open','Lon','Lat','MSLP','CV') 
  outF=paste('ECLfixes_aus_',output,'.csv',sep="")
  write.csv(Fixes,outF)
}

##Let's play with seasonality
ECLfig_seas<-function(data,name,smon,sname)
{
  mon=floor((data[,2] %% 10000)/100)
  if(smon[2]>smon[1]) data=data[mon>=smon[1] & mon<=smon[2],]
  else data=data[mon>=smon[1] | mon<=smon[2],]
  
  lats=seq(-49.5,-10.5,1)
  lons=seq(110.5,179.5,1)
  
  all<-matrix(0,length(lats),length(lons))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
    {
      I=which(data[,5]>=floor(lons[i]) & data[,5]<ceiling(lons[i]) & data[,6]>=floor(lats[j]) & data[,6]<ceiling(lats[j]))
      all[j,i]=length(I)
    }
  all=(all/19)
  I=which(all>1.95)
  all[I]=1.95
  
  tiff(file=paste("Figures/Locs_aus_",name,"_",sname,".tiff",sep=""), height=450, width=600)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=seq(0,2,0.2),col=cols,zlim=c(0,2))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("Figures/Locs_esb_",name,"_",sname,".tiff",sep=""), height=400, width=400)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=seq(0,2,0.2),col=cols,zlim=c(0,2),xlim=c(149,161),ylim=c(-41,-24))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

ECLfig_ann<-function(data,name)
{
  lats=seq(-49.5,-10.5,1)
  lons=seq(110.5,179.5,1)
  
  all<-matrix(0,length(lats),length(lons))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
    {
      I=which(data[,5]>=floor(lons[i]) & data[,5]<ceiling(lons[i]) & data[,6]>=floor(lats[j]) & data[,6]<ceiling(lats[j]))
      all[j,i]=length(I)
    }
  all=(all/19)
  I=which(all>2.95)
  all[I]=2.95
  tiff(file=paste("Figures/Locs_aus_",name,".tiff",sep=""), height=450, width=600)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=seq(0,3,0.3),col=cols,zlim=c(0,3))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("Figures/Locs_esb_",name,".tiff",sep=""), height=400, width=400)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=seq(0,3,0.3),col=cols,zlim=c(0,3),xlim=c(149,161),ylim=c(-41,-24))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

ECL_FA_seas<-function(data,name,smon,sname)
{
  data=data[data[,2]<20070000,]
  data[,9:10]=0
  data[,9]=data[,2] %in% ecl
  data[,10]=data[,2] %in% rain2[,3]
  
  mon=floor((data[,2] %% 10000)/100)
  if(smon[2]>smon[1]) data=data[mon>=smon[1] & mon<=smon[2],]
  else data=data[mon>=smon[1] | mon<=smon[2],]
  
  lats=seq(-40.5,-24,1.5)
  lons=seq(149.5,161,1.5)
  
  FA<-FAr<-matrix(0,length(lats),length(lons))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
    {
      I=which(data[,5]>=lons[i]-0.75 & data[,5]<lons[i]+0.75 & data[,6]>=lats[j]-0.75 & data[,6]<lats[j]+0.75)
      FA[j,i]=sum(data[I,9])/length(I) ##% of fixes w an MLDB low at 00Z on same/following day
      J=which(data[I,10]==1) ##Rain
      FAr[j,i]=sum(data[I[J],9])/length(J)
    }
  tiff(file=paste("Figures/FA_",name,"_",sname,".tiff",sep=""), height=400, width=400)
  image.plot(lons,lats,t(FA),xlab="",ylab="",breaks=seq(0,1,0.1),col=cols,zlim=c(0,1))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("Figures/FAr_",name,"_",sname,".tiff",sep=""), height=400, width=400)
  image.plot(lons,lats,t(FAr),xlab="",ylab="",breaks=seq(0,1,0.1),col=cols,zlim=c(0,1))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

ECL_rain_seas<-function(data,name,smon,sname)
{
  data[,9]=data[,2] %in% rain2[,3]
  
  mon=floor((data[,2] %% 10000)/100)
  if(smon[2]>smon[1]) data=data[mon>=smon[1] & mon<=smon[2],]
  else data=data[mon>=smon[1] | mon<=smon[2],]
  
  lats=seq(-40.5,-24,1.5)
  lons=seq(149.5,161,1.5)
  
  rr<-matrix(0,length(lats),length(lons))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
    {
      I=which(data[,5]>=lons[i]-0.75 & data[,5]<lons[i]+0.75 & data[,6]>=lats[j]-0.75 & data[,6]<lats[j]+0.75)
      rr[j,i]=sum(data[I,9])/length(I) ## of fixes w/ Morwenna rain on subs day
    }
  tiff(file=paste("Figures/Rprop_",name,"_",sname,".tiff",sep=""), height=400, width=400)
  image.plot(lons,lats,t(rr),xlab="",ylab="",breaks=seq(0,1,0.1),col=cols,zlim=c(0,1))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

##Returns the annual count of lows for a lat/lon grid & season given
##Assumes first column is date - if there's still a random column it must be removed in call
ECLcount_seas<-function(data,lats,lons,smon)
{
  mon=floor((data$Date %% 10000)/100)
  yy=floor(data$Date/10000) 
  dd=(lons[2]-lons[1])/2 ##This is the range for each lon/lat
  
  if(smon[2]>smon[1])
  {
    I=which(mon>=smon[1] & mon<=smon[2])
    data=data[I,]
    yy=yy[I]
    years=seq(min(yy),max(yy))
  }
  else 
  {
    I=which(mon>=smon[1] | mon<=smon[2])
    data=data[I,]
    yy=yy[I]
    J=which(mon[I]<6)
    yy[J]=yy[J]-1
    years=seq(min(yy)+1,max(yy)-1)
  }
  
  ann<-array(0,dim=c(length(lats),length(lons),length(years)))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
      for(k in 1:length(years))
      {
        I=which(data$Lon>=lons[i]-dd & data$Lon<lons[i]+dd & data$Lat>=lats[j]-dd & data$Lat<lats[j]+dd & yy==years[k] )
        ann[j,i,k]=length(I)
      }
  
  return(ann)
}

ECLfig_seas2<-function(data,name,smon,sname,ny)
{
  mon=floor((data[,2] %% 10000)/100)
  if(smon[2]>smon[1]) data=data[mon>=smon[1] & mon<=smon[2],]
  else data=data[mon>=smon[1] | mon<=smon[2],]
  
  lats=seq(-49.5,-10.5,1)
  lons=seq(110.5,179.5,1)
  
  all<-matrix(0,length(lats),length(lons))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
    {
      I=which(data[,5]>=floor(lons[i]) & data[,5]<ceiling(lons[i]) & data[,6]>=floor(lats[j]) & data[,6]<ceiling(lats[j]))
      all[j,i]=length(I)
    }
  all=(all/ny)
  if(smon[2]-smon[1]>10)
  {
    I=which(all>2.95)
    all[I]=2.95 
    bb=seq(0,3,0.3)
  }
  else
  {
    I=which(all>1.95)
    all[I]=1.95
    bb=seq(0,2,0.2)
  }
  
  tiff(file=paste("Figures/Locs_aus_",name,"_",sname,".tiff",sep=""), height=450, width=600)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("Figures/Locs_esb_",name,"_",sname,".tiff",sep=""), height=400, width=400)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)),xlim=c(149,161),ylim=c(-41,-24))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

##Uses 2 degree boxes instead
ECLfig_seas2a<-function(data,name,smon,sname,ny=19)
{
  mon=floor((data[,2] %% 10000)/100)
  if(smon[2]>smon[1]) data=data[mon>=smon[1] & mon<=smon[2],]
  else data=data[mon>=smon[1] | mon<=smon[2],]
  
  lats=seq(-48,-10,2)
  lons=seq(110,180,2)
  
  all<-matrix(0,length(lats),length(lons))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
    {
      I=which(data[,5]>=lons[i]-1 & data[,5]<lons[i]+1 & data[,6]>=lats[j]-1 & data[,6]<lats[j]+1)
      all[j,i]=length(I)
    }
  all=(all/ny)


  if(smon[2]-smon[1]>10)
  {
    I=which(all>9.95)
    all[I]=9.95 
    bb=seq(0,10)
  }
  else
  {
    I=which(all>4.95)
    all[I]=4.95 
    bb=seq(0,5,0.5)
  }
  
  tiff(file=paste("Figures/Locs_aus_",name,"_",sname,".tiff",sep=""), height=450, width=600)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("Figures/Locs_esb_",name,"_",sname,".tiff",sep=""), height=400, width=400)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)),xlim=c(149,161),ylim=c(-41,-24))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

ECLfig_seas2b<-function(data,name,smon,sname,ny,b)
{
  mon=floor((data[,2] %% 10000)/100)
  if(smon[2]>smon[1]) data=data[mon>=smon[1] & mon<=smon[2],]
  else data=data[mon>=smon[1] | mon<=smon[2],]
  
  lats=seq(-48,-10,2)
  lons=seq(110,180,2)
  
  all<-matrix(0,length(lats),length(lons))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
    {
      I=which(data[,5]>=lons[i]-1 & data[,5]<lons[i]+1 & data[,6]>=lats[j]-1 & data[,6]<lats[j]+1)
      all[j,i]=length(I)
    }
  all=(all/ny)
  
  
  if(smon[2]-smon[1]>10)
  {
    I=which(all>b-0.05)
    all[I]=b-0.05
    bb=seq(0,b,length.out=11)
  }
  else
  {
    I=which(all>((b/2)-0.05))
    all[I]=(b/2)-0.05
    bb=seq(0,b/2,length.out=11)
  }
  
  tiff(file=paste("Figures/Locs_aus_",name,"_",sname,".tiff",sep=""), height=450, width=600)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("Figures/Locs_esb_",name,"_",sname,".tiff",sep=""), height=400, width=400)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)),xlim=c(149,161),ylim=c(-41,-24))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

ECLfig_seas2c<-function(data,name,smon,sname,ny,b)
{
  mon=floor((data[,2] %% 10000)/100)
  if(smon[2]>smon[1]) data=data[mon>=smon[1] & mon<=smon[2],]
  else data=data[mon>=smon[1] | mon<=smon[2],]
  
  lats=seq(-48,-10,2)
  lons=seq(110,180,2)
  
  all<-matrix(0,length(lats),length(lons))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
    {
      I=which(data[,5]>=lons[i]-1 & data[,5]<lons[i]+1 & data[,6]>=lats[j]-1 & data[,6]<lats[j]+1)
      all[j,i]=length(unique(data[I,2]))
    }
  all=(all/ny)
  
  
  if(smon[2]-smon[1]>10)
  {
    I=which(all>b-0.05)
    all[I]=b-0.05
    bb=seq(0,b,length.out=11)
  }
  else
  {
    I=which(all>((b/2)-0.05))
    all[I]=(b/2)-0.05
    bb=seq(0,b/2,length.out=11)
  }
  
  tiff(file=paste("Figures/Locs_aus_",name,"_",sname,".tiff",sep=""), height=450, width=600)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
  tiff(file=paste("Figures/Locs_esb_",name,"_",sname,".tiff",sep=""), height=400, width=400)
  image.plot(lons,lats,t(all),xlab="",ylab="",breaks=bb,col=cols,zlim=c(min(bb),max(bb)),xlim=c(149,161),ylim=c(-41,-24))
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  dev.off()
}

ECLcount_all<-function(data,lats,lons)
{
  dd=(lons[2]-lons[1])/2 ##This is the range for each lon/lat
  ann<-array(0,dim=c(length(lats),length(lons)))
  for(i in 1:length(lons))
    for(j in 1:length(lats))
    {
        I=which(data$Lon>=lons[i]-dd & data$Lon<lons[i]+dd & data$Lat>=lats[j]-dd & data$Lat<lats[j]+dd)
        ann[j,i]=length(I)
    }  
  return(ann)
}
