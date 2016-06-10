####### Okay, taking bits from existing code for ECL typing.
### First, we need to load in wrf and restrict to just one year where a range of MLDB ECLs
rm(list=ls())
library(R.matlab)
library(RNetCDF)
library(abind)
library(sp)
library(geosphere)
library(abind)

setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing")

fixes<-events<-list()

for(r in 1:3) {
  fixes[[r]]=read.csv(paste("ECLfixes_d01_0708_R",r,"_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
  fixes[[r]]$Date2=as.POSIXct(paste(as.character(fixes[[r]]$Date),substr(fixes[[r]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  events[[r]]=read.csv(paste("ECLevents_d01_0708_R",r,"_rad2_p100_typing.csv",sep=""),stringsAsFactors=F)
}

##Comparing days? Useless

dates=seq.POSIXt(as.POSIXct("2007010100",format="%Y%m%d%H",tz="GMT"),as.POSIXct("2008123118",format="%Y%m%d%H",tz="GMT"),by="day")
dates2=data.frame(Date=dates,R1=rep("NA",length(dates)),R2=rep("NA",length(dates)),R3=rep("NA",length(dates)),stringsAsFactors=F)

for(i in 1:(length(dates)-1))
  for(r in 1:3)
  {
    I=which(fixes[[r]]$Date2>=dates[i] & fixes[[r]]$Date2<dates[i+1] & fixes[[r]]$Location==1)
    if(length(I)>0)
      {
      a=unique(fixes[[r]]$Type[I])
      if(length(a)>1) dates2[i,r+1]="Mixed" else dates2[i,r+1]=a
    }
  }

J=which(dates2[,2:4]=="NA",arr.ind=T)
j=sort(unique(J[,1]))
dates3=dates2[-j,]

### Comparing events, using ERAI as the comparator

erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)
erai_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)

erai=erai[,-c(1,2)]
erai_E=erai_E[,-c(1,2)]
erai$Date2=as.POSIXct(paste(as.character(erai$Date),substr(erai$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

erai_E$R3H<-erai_E$R2H<-erai_E$R1H<-erai_E$R3B<-erai_E$R2B<-erai_E$R1B<-rep("NA",length(erai_E[,1]))

library(sp)
for(i in 1:length(erai_E[,1]))
{
  I=which(erai$ID==erai_E$ID[i] & erai$Location==1)
  tmp=erai[I,]
  for(r in 1:3)
  {
    IDs=NaN
    for(x in 1:length(tmp[,1]))
    {
    J=which(fixes[[r]]$Date2>=tmp$Date2[x]-21600 & fixes[[r]]$Date2<=tmp$Date2[x]+21600)
    if(length(J)>0)
    {
    tmp2=fixes[[r]][J,]
    dist=spDistsN1(as.matrix(tmp2[,7:8]),as.numeric(tmp[x,6:7]),longlat=T)
    J=which(dist<=500)
    if(length(J)>0) IDs=c(IDs,unique(tmp2$ID[J]))
    }
    }
    IDs=sort(unique(IDs))
  if(length(IDs)==1)
  {
    I=which(events[[r]]$ID==IDs)
    erai_E[i,13+r]=events[[r]]$TypeSB[I]
    erai_E[i,16+r]=events[[r]]$HartType[I]
  } else if(length(IDs)>1)
  {
    I=which(events[[r]]$ID %in% IDs)
    SB=unique(events[[r]]$TypeSB[I])
    if(length(SB)>1) erai_E[i,13+r]="Many" else erai_E[i,13+r]=SB
    H=unique(events[[r]]$HartType[I])
    if(length(H)>1) erai_E[i,16+r]="Many" else erai_E[i,16+r]=H
  }
    
  }
}


######### What about lazy approach - count types for each version

####### Okay, taking bits from existing code for ECL typing.
### First, we need to load in wrf and restrict to just one year where a range of MLDB ECLs
rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing")

dom="d01"
cat="rad2_p100"

wnames=c("R1","R2","R3",
         "R1_notopo","R2_notopo","R3_notopo",
         "R1_BRAN","R2_BRAN","R3_BRAN",
         "R1_BRAN_noeac","R2_BRAN_noeac","R3_BRAN_noeac",
         "R2_BRAN_2eac","R2_BRAN_2eac","R3_BRAN_2eac")

Htypes=c("EC","SC","TC","Mixed")
Btypes<-c("ET","SSL","IT","CL")

Hart<-Browning<-matrix(0,16,4)
rownames(Hart)<-rownames(Browning)<-c(wnames,"ERAI")
colnames(Hart)<-Htypes
colnames(Browning)<-Btypes

for(w in 1:15)
{
  events=read.csv(paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_typing.csv",sep=""))
  for(k in 1:4)
  {
    I=which(events$HartType==Htypes[k])
    Hart[w,k]=length(I)
    I=which(events$TypeSB==Btypes[k])
    Browning[w,k]=length(I)
  }
}
w=16
for(k in 1:4)
{
  I=which(erai_E$HartType==Htypes[k])
  Hart[w,k]=length(I)
  I=which(erai_E$TypeSB==Btypes[k])
  Browning[w,k]=length(I)
}

#What about days associated with the events
for(w in 1:15)
{
  events=read.csv(paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_typing.csv",sep=""))
  fixes=read.csv(paste("ECLfixes_",dom,"_0708_",wnames[w],"_",cat,"_typing.csv",sep=""))
  for(k in 1:4)
  {
    I=which(events$HartType==Htypes[k])
    a=unique(events$ID[I])
    J=which(fixes$ID %in% a & fixes$Location==1)
    Hart[w,k]=length(unique(fixes$Date2[J]))
    
    I=which(events$TypeSB==Btypes[k])
    a=unique(events$ID[I])
    J=which(fixes$ID %in% a & fixes$Location==1)
    Browning[w,k]=length(unique(fixes$Date2[J]))
  }
}
w=16
for(k in 1:4)
{
  I=which(erai_E$HartType==Htypes[k])
  a=unique(erai_E$ID[I])
  J=which(erai$ID %in% a & erai$Location==1)
  Hart[w,k]=length(unique(erai$Date2[J]))
  
  I=which(erai_E$TypeSB==Btypes[k])
  a=unique(erai_E$ID[I])
  J=which(erai$ID %in% a & erai$Location==1)
  Browning[w,k]=length(unique(erai$Date2[J]))
}

Other=matrix(0,16,4)
rownames(Other)=c(wnames,"ERAI")
colnames(Other)=c("Formed","Entered","Intensified","Bomb")
for(w in 1:15)
{
  events=read.csv(paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_typing.csv",sep=""))
  for(k in 1:3) Other[w,k]=length(which(events$EnteredFormed==k))
  Other[w,4]=sum(events$Bomb)
}
w=16
for(k in 1:3) Other[w,k]=length(which(erai_E$EnteredFormed==k))
Other[w,4]=sum(erai_E$Bomb)