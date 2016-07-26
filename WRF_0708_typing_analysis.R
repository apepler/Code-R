rm(list=ls())
setwd("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing/")
library(ggplot2)
library(reshape2)
library(abind)
library(RNetCDF)

dom="d01"
cat="rad2_p100"

wnames=c("R1","R2","R3",
         "R1_notopo","R2_notopo","R3_notopo",
         "R1_BRAN","R2_BRAN","R3_BRAN",
         "R1_BRAN_noeac","R2_BRAN_noeac","R3_BRAN_noeac",
         "R2_BRAN_2eac","R2_BRAN_2eac","R3_BRAN_2eac")

events<-fixes<-events2<-fixes2<-list()

for(w in 1:length(wnames))
{
  dom="d01"
  
  fixes[[w]]=read.csv(paste("ECLfixes_",dom,"_0708_",wnames[w],"_",cat,"_typing_impacts.csv",sep=""),stringsAsFactors=F)
  fixes[[w]]$Date2=as.POSIXct(paste(as.character(fixes[[w]]$Date),substr(fixes[[w]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
  events[[w]]=read.csv(paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_typing_impacts.csv",sep=""))
  
  if(w<13)
  {
    dom="d02"
    
    fixes2[[w]]=read.csv(paste("ECLfixes_",dom,"_0708_",wnames[w],"_",cat,"_typing_impacts.csv",sep=""),stringsAsFactors=F)
    fixes2[[w]]$Date2=as.POSIXct(paste(as.character(fixes2[[w]]$Date),substr(fixes2[[w]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
    events2[[w]]=read.csv(paste("ECLevents_",dom,"_0708_",wnames[w],"_",cat,"_typing_impacts.csv",sep=""))
  }
}


statave<-matrix(0,15,18)
rownames(statave)=wnames
colnames(statave)=c("Count","Length2","MSLP2","CV2",
                    "EC","SC","TC","Mixed",
                    "MeanRain>=6","MaxRain>=50","MeanWind>=50km/h","MaxWind>=80km/h",
                    "Length>=5","CV>2","MSLP<=990","Bomb",
                    "Warm","Cold")

impcol=c(17:20,8,10) ## For 500, 22:25 for 250
xthresh=c(6,50,13.9,22.2,5,2)
types=c("EC","SC","TC","Mixed")

statave2=statave[1:12,]

for(w in 1:length(wnames))
{
  tmp=events[[w]]
  statave[w,1]=length(tmp[,1])
  statave[w,2:4]=apply(tmp[,8:10],2,mean,na.rm=T)
  for(x in 1:4) statave[w,x+4]=length(which(tmp$HartType==types[x]))
  for(x in 1:6) statave[w,x+8]=length(which(tmp[,impcol[x]]>=xthresh[x]))
  statave[w,15]=length(which(tmp$MSLP2<=990))
  statave[w,16]=sum(tmp$Bomb)
  
  month=floor(tmp$Date1/100)%%100
  I=which(month<=3 | month>=9)
  statave[w,17]=length(I)
  statave[w,18]=length(month)-length(I)
}

for(w in 1:12)
{
  tmp=events2[[w]]
  statave2[w,1]=length(tmp[,1])
  statave2[w,2:4]=apply(tmp[,8:10],2,mean,na.rm=T)
  for(x in 1:4) statave2[w,x+4]=length(which(tmp$HartType==types[x]))
  for(x in 1:6) statave2[w,x+8]=length(which(tmp[,impcol[x]]>=xthresh[x]))
  statave2[w,15]=length(which(tmp$MSLP2<=990))
  statave2[w,16]=sum(tmp$Bomb)
  
  month=floor(tmp$Date1/100)%%100
  I=which(month<=3 | month>=9)
  statave2[w,17]=length(I)
  statave2[w,18]=length(month)-length(I)
}


##So, now split further - Cool/EC, Cool/Other, Warm/EC, Warm/Other

statave3=matrix(NaN,15,10)

rownames(statave3)=wnames
colnames(statave3)=c("Count","Warm/EC","Warm/Other","Cool/EC","Cool/Other",
                    "Count2","Warm/EC 2","Warm/Other 2","Cool/EC 2","Cool/Other 2")

for(w in 1:length(wnames))
{
  tmp=events[[w]]
  statave3[w,1]=length(tmp[,1])
  
  month=floor(tmp$Date1/100)%%100
  I=which(month<=4 | month>=11)
  J=which(tmp$HartType[I]=="EC")
  statave3[w,2]=length(J)
  statave3[w,3]=length(I)-length(J)
  
  I=which(month>4 & month<11)
  J=which(tmp$HartType[I]=="EC")
  statave3[w,4]=length(J)
  statave3[w,5]=length(I)-length(J)
  
  if(w<=12)
  {
  tmp=events2[[w]]
  statave3[w,6]=length(tmp[,1])
  
  month=floor(tmp$Date1/100)%%100
  I=which(month<=4 | month>=11)
  J=which(tmp$HartType[I]=="EC")
  statave3[w,7]=length(J)
  statave3[w,8]=length(I)-length(J)
  
  I=which(month>4 & month<11)
  J=which(tmp$HartType[I]=="EC")
  statave3[w,9]=length(J)
  statave3[w,10]=length(I)-length(J)
  }
}

### Monthly count - both events & days

monthly<-array(NaN,c(12,15,4))
dimnames(monthly)[[2]]=wnames
dimnames(monthly)[[3]]=c("d01 events","d01 days","d02 events","d02 days")

for(w in 1:length(wnames))
{
  month=list()
  month[[1]]=floor(events[[w]]$Date1/100)%%100
  month[[2]]=floor(fixes[[w]]$Date[fixes[[w]]$Location==1]/100)%%100
  if(w<=12)
  {
  month[[3]]=floor(events2[[w]]$Date1/100)%%100
  month[[4]]=floor(fixes2[[w]]$Date[fixes2[[w]]$Location==1]/100)%%100
  for(m in 1:12)
    for(n in 1:4)
    {
      I=which(month[[n]]==m)
      monthly[m,w,n]=length(I)
    }
  } else {
    for(m in 1:12)
      for(n in 1:2)
      {
        I=which(month[[n]]==m)
        monthly[m,w,n]=length(I)
      }
  }
}

monthly2=abind(apply(monthly[3:5,,],c(2,3),sum),
               apply(monthly[6:8,,],c(2,3),sum),
               apply(monthly[9:11,,],c(2,3),sum),
               apply(monthly[c(12,1,2),,],c(2,3),sum),along=1)


###############
#################
### Play with typing

erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100_typing0708.csv")

I=which(erai$Location==1)
plot(erai$VL[I],erai$VU[I],xlim=c(-500,200),ylim=c(-500,200))
abline(h=0,col="red")
abline(v=0,col="red")

cor(erai$B[I],erai$VL[I])

plot(erai$VL[I],erai$B[I],xlim=c(-500,200),ylim=c(-50,200))
abline(h=0,col="red")
abline(v=0,col="red")

## How does this compare to WRF
I=which(erai$Location==1)
plot(erai$VL[I],erai$VU[I],xlim=c(-500,200),ylim=c(-500,200))

clist=c("red","green","blue")
for(i in 1:3)
{
  I=which(fixes[[i]]$Location==1)
  points(fixes[[i]]$VL[I],fixes[[i]]$VU[I],col=clist[i],pch=4)
}
abline(h=0,col="red")
abline(v=0,col="red")

I=which(erai$Location==1)
plot(erai$VL[I],erai$B[I],xlim=c(-500,200),ylim=c(-50,200))

clist=c("red","green","blue")
for(i in 1:3)
{
  I=which(fixes[[i]]$Location==1)
  points(fixes[[i]]$VL[I],fixes[[i]]$B[I],col=clist[i],pch=4)
}
abline(h=0,col="red")
abline(v=0,col="red")
