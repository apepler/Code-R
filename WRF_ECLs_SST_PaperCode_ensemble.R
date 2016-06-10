######### Intial stuff
rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007-06')
library(abind)
library("R.matlab")
library(fields)
library(map)
figdir="~/Documents/ECLs/WRFruns/0708/EACpaper/"

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}
source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))

cat=c("d01_p100_rad2cv1","d02_p100_rad2cv1")
dom2=c("d01","d02")
cat2=c("rad2_p100","rad2_p100")

SSTtype=c("Control","NoTopo","NoEAC","Both")
count<-array(NaN,c(3,3,5))

events<-fixes<-list()

c=1
n=1
  for(k in 1:3)
    for(r in 1:3)
  {
    dir=c(paste("ERAI_R",r,"_ensemble_BRAN/",cat[c],"/",sep=""),
          paste("ERAI_R",r,"_ensemble_BRAN_noeac/",cat[c],"/",sep=""),
          paste("ERAI_R",r,"_ensemble_BRAN_2eac/",cat[c],"/",sep=""))
    n2=1
  for(day in 27:31)
    for(hour in c("00"))
    {
      count[k,r,n2]=length(read.csv(paste(dir[k],"ECLevents_200705",day,hour,".csv",sep=""))[,1])
      if(n2==1)
        {
          
          events[[n]]=read.csv(paste(dir[k],"ECLevents_200705",day,hour,".csv",sep=""))
          fixes[[n]]=read.csv(paste(dir[k],"ECLfixes_200705",day,hour,".csv",sep=""))
        } else {
          
          events[[n]]=rbind(events[[k]],read.csv(paste(dir[k],"ECLevents_200705",day,hour,".csv",sep="")))
          fixes[[n]]=rbind(fixes[[k]],read.csv(paste(dir[k],"ECLfixes_200705",day,hour,".csv",sep="")))
        }
      n2=n2+1 
      }
      n=n+1    
    }

####### Also need ERAI

erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)
erai_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100_typing0708.csv",stringsAsFactors = F)

erai=erai[,-c(1,2)]
erai_E=erai_E[,-c(1,2)]
erai$Date2=as.POSIXct(paste(as.character(erai$Date),substr(erai$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

I=which(erai$Date>=20070601 & erai$Date<=20070630 & erai$Location==1)
a=unique(erai$ID[I])

erai_E=erai_E[erai_E$ID %in% a,]
erai=erai[erai$ID %in% a,]
eraiL=erai[erai$Location==1,]

####### For all the events, check if erai is matched 

erai_match<-array(0,c(5,3,3,5))
erai_match_cv<-array(NaN,c(5,3,3,5))
day=27:31

for(k in 1:3)
  for(r in 1:3)
  {
    dir=c(paste("ERAI_R",r,"_ensemble_BRAN/",cat[c],"/",sep=""),
          paste("ERAI_R",r,"_ensemble_BRAN_noeac/",cat[c],"/",sep=""),
          paste("ERAI_R",r,"_ensemble_BRAN_2eac/",cat[c],"/",sep=""))
    for(d in 1:5)
      {
          fixes=read.csv(paste(dir[k],"ECLfixes_200705",day[d],hour,".csv",sep=""))
          fixes$Date2=as.POSIXct(paste(as.character(fixes$Date),substr(fixes$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
          for(i in 1:5)
          {
            tmp=erai[erai$ID==erai_E$ID[i] & erai$Location==1,]
            rn=range(tmp$Date2)
            I=which(fixes$Date2<=rn[2]+(60*60*6) & fixes$Date2>=rn[1]-(60*60*6) & fixes$Location==1)
            if(length(I)>0)
              {
              erai_match[i,k,r,d]=1
              erai_match_cv[i,k,r,d]=max(fixes$CV[I])
          }

        }
    }
  }

a=apply(erai_match,c(1,3),sum)

erai_match_cv2=array(NaN,c(5,2,3,5))
for(i in 1:2)
  for(j in 1:3)
    for(k in 1:5)
      erai_match_cv2[,i,j,k]=erai_match_cv[,i+1,j,k]-erai_match_cv[,1,j,k]

######## ECL dates 
date=seq.POSIXt(from=as.POSIXct("2007060100",format="%Y%m%d%H",tz="GMT"),
                to=as.POSIXct("2007063018",format="%Y%m%d%H",tz="GMT"),
                by="6 hours")

erai_date=matrix(0,length(date),2)
for(i in 1:length(date))
{
  I=which(erai$Date2==date[i] & erai$Location==1)
  if(length(I)>0) 
  {
    erai_date[i,1]=1
    erai_date[i,2]=max(erai$CV[I])
  }
}
  
wrf_date<-wrf_cv<-array(0,c(length(date),3,3,5))

day=27:31

for(k in 1:3)
  for(r in 1:3)
  {
    dir=c(paste("ERAI_R",r,"_ensemble_BRAN/",cat[c],"/",sep=""),
          paste("ERAI_R",r,"_ensemble_BRAN_noeac/",cat[c],"/",sep=""),
          paste("ERAI_R",r,"_ensemble_BRAN_2eac/",cat[c],"/",sep=""))
    for(d in 1:5)
    {
      fixes=read.csv(paste(dir[k],"ECLfixes_200705",day[d],hour,".csv",sep=""))
      fixes$Date2=as.POSIXct(paste(as.character(fixes$Date),substr(fixes$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
      
      for(i in 1:120)
        {
      I=which(fixes$Date2==date[i] & fixes$Location==1)
      if(length(I)>0) 
      {
        wrf_date[i,k,r,d]=1
        wrf_cv[i,k,r,d]=max(fixes$CV[I])
      }
      }
      
    }
  }

###Amount of time?

a=apply(wrf_date,c(2,3,4),sum)
apply(a,1,mean)/4

t.test(as.vector(a[3,,]-a[1,,]))


wrf_date2=apply(wrf_date,c(1,2),sum)/15
wrf_cv[wrf_cv==0]=NaN
wrf_cv2=apply(wrf_cv,c(1,2),mean,na.rm=T)

erai_date[erai_date[,2]==0,2]=NaN
SSTtype=c("Control","NoEAC","2EAC")
clist=c("grey","blue","red")
pdf(file=paste(figdir,"DailyCV_Jun07_all_",cat[c],".pdf",sep=""),width=7,height=4)
plot(date,erai_date[,2],type="l",lwd=3,col="black",
     ylim=c(0,ceiling(max(cbind(wrf_cv2,erai_date[,2]),na.rm=T))),xlab="Date",ylab="CV",main="Mean ECL curvature across all runs")
for(k in 1:3) lines(date,wrf_cv2[,k],lwd=3,col=clist[k])
lines(date,erai_date[,2],lwd=3,col="black")
legend("topleft",legend=c("ERAI",SSTtype),col=c("black",clist),lwd=3,bty='n')
dev.off()

erai_date[erai_date[,1]==0,1]=NaN
pdf(file=paste(figdir,"ECLprop_Jun07_all_",cat[c],".pdf",sep=""),width=7,height=4)
plot(date,erai_date[,1]*100,type="l",lwd=3,col="black",
     ylim=c(0,100),xlab="Date",ylab="CV",main="Proportion of members with an ECL present")
for(k in 1:3) lines(date,wrf_date2[,k]*100,lwd=3,col=clist[k])
lines(date,erai_date[,1]*100,lwd=3,col="black")
legend("topleft",legend=c("ERAI",SSTtype),col=c("black",clist),lwd=3,bty='n')
dev.off()


############



for(k in 1:6)  fixes[[k]]$Date2=fixes[[k]]$Date+(as.numeric(fixes[[k]]$Time)-1)/4
dates=data.frame(Date=seq(20070601,20070630.75,0.25),matrix(0,120,6))
for(i in 1:120)
  for(k in 1:6)
  {
    I=which(fixes[[k]]$Date2==dates$Date[i] & fixes[[k]]$Location==1)
    if(length(I)>0) dates[i,k+1]=mean(fixes[[k]]$CV[I])
  } 

plot(dates$Date,dates[,2],type="l",lwd=3,col=clist[1],
     ylim=c(0,ceiling(max(dates[,2:5],na.rm=T))),xlab="Date",ylab="CV",main="Mean ECL curvature across all runs")
for(k in 2:4) lines(dates$Date,dates[,k+1],lwd=3,col=clist[k])
legend("topleft",legend=SSTtype,col=clist,lwd=3,bty='n')