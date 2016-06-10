rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/')
library(abind)

makePDF = function(data1,data2,xlabel,labloc="topleft") {
  a=density(data1,na.rm=T)
  b=density(data2,na.rm=T)
  
  lims=range(data1,data2,na.rm=T)
  if((lims[2]-lims[1])<10)
  {
    lims[1]=floor(lims[1])
    lims[2]=ceiling(lims[2])
  } else {
    lims[1]=floor(lims[1]/5)*5
    lims[2]=ceiling(lims[2]/5)*5
  }
  
  plot(a,col=rgb(0,0,1,1/4),xlim=lims,ylim=range(0,a$y,b$y),
       xlab=xlabel,ylab="Frequency",cex.main=1.2,
       main="2007-2008 ECL statistics")
  polygon(a,col=rgb(0,0,1,1/4),density=-1)
  polygon(b,col=rgb(1,0,0,1/4),density=-1)
  legend(labloc,legend=c("Control","NoTopo"),
         col=c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)),lwd=4,cex=1,bty="n")   
}

############

year=c(2007,2008)
domain=c("d01","d02")
cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5")

count<-countBRAN<-count_noeac<-count_2eac<-matrix(0,10,3)
sigN<-sig2<-matrix(0,10,4)

for(r in 1:3)
{
  n=1
  for(dom in domain)
    for(c in 1:5)
  {
    for(year in 2007:2008)
    {
      
      count[n,r]=count[n,r]+length(read.csv(paste("ECLevents_",dom,"_",year,"_R",r,"_",cat[c],".csv",sep=""))[,1])
      countBRAN[n,r]=countBRAN[n,r]+length(read.csv(paste("ECLevents_",dom,"_",year,"_R",r,"_BRAN_",cat[c],".csv",sep=""))[,1])
      count_noeac[n,r]=count_noeac[n,r]+length(read.csv(paste("ECLevents_",dom,"_",year,"_R",r,"_BRAN_noeac_",cat[c],".csv",sep=""))[,1])
      if(dom=="d01") count_2eac[n,r]=count_2eac[n,r]+length(read.csv(paste("ECLevents_",dom,"_",year,"_R",r,"_BRAN_2eac_",cat[c],".csv",sep=""))[,1])
    }
    n=n+1
    }}

diffNo=count_noeac/countBRAN
diff2=count_2eac[1:5,]/countBRAN[1:5,]

#### Different version based on ECLs p.a., so with thresh


cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240","rad2_p100_cv0.5")
cvthresh=c(80,70,60,50,40,30,20,10,NaN)
count<-array(NaN,c(10,5,3,8))


for(r in 1:3)
{
  n=1
  for(dom in c("d01","d02"))
    for(c in 1:5)
    {
      data=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                 read.csv(paste("ECLevents_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
      data2=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_BRAN_",cat[c],".csv",sep="")),
                  read.csv(paste("ECLevents_",dom,"_2008_R",r,"_BRAN_",cat[c],".csv",sep="")))
      data3=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")),
                  read.csv(paste("ECLevents_",dom,"_2008_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")))
      data5=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_notopo_",cat[c],".csv",sep="")),
                  read.csv(paste("ECLevents_",dom,"_2008_R",r,"_notopo_",cat[c],".csv",sep="")))
      if(dom=="d01")
      {
        data4=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_BRAN_2eac_",cat[c],".csv",sep="")),
                    read.csv(paste("ECLevents_",dom,"_2008_R",r,"_BRAN_2eac_",cat[c],".csv",sep="")))
      }
      
      a=order(data2$CV2,decreasing=T)
      for(t in 1:8)
        if(length(a)>=cvthresh[t])
      {
        b=data2$CV2[a[cvthresh[t]]]
        if(t==8) c=10 else c=data2$CV2[a[cvthresh[t+1]]]
        count[n,1,r,t]=length(which(data$CV2>=b & data$CV2<c))
        count[n,2,r,t]=length(which(data2$CV2>=b & data2$CV2<c))
        count[n,3,r,t]=length(which(data3$CV2>=b & data3$CV2<c))
        count[n,5,r,t]=length(which(data5$CV2>=b & data5$CV2<c))
        if(dom=="d01") count[n,4,r,t]=length(which(data4$CV2>=b & data4$CV2<c))
      }
      n=n+1 
    }
}

for(i in 2:8) print(mean(apply(count[c(1:2,4:7,9:10),4,,i:8],c(1,2),sum)/apply(count[c(1:2,4:7,9:10),2,,i:8],c(1,2),sum),na.rm=T))


## Make as boxplot like I did for NARCliM stuff?

change=100*(abind(count[,3,,2:8]/count[,2,,2:8],
                  count[,4,,2:8]/count[,2,,2:8],along=4)-1)

dimnames(change)[[2]]=c("R1","R2","R3")
dimnames(change)[[3]]=c(paste(cvthresh[2:7]/2,"-",cvthresh[3:8]/2),"<5")
dimnames(change)[[4]]=c("NoEAC","2EAC")
names(dimnames(change))=c("Source","WRF version","Threshold","Model")

pdf("~/Documents/ECLs/WRFruns/0708/Boxplot_ECLevents_change_cvthresh_noeac_2eac_vBRAN_d01_v2.pdf",width=9,height=4)
par(mar=c(4,4,1,1))
library(ggplot2)
library(reshape2)
tmp=melt(change[1:4,,,])
print(ggplot(tmp, aes(x = Threshold, y = value, fill = Model)) +
        geom_boxplot() +
        scale_fill_manual(values = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4)))  +
        scale_y_continuous(breaks=seq(-100, 150, 50)) +
        theme_bw() + ylab("Percentage change") + xlab("Number of events p.a. in control") +  geom_hline(yintercept = 0))
dev.off()


change=100*(count[1:4,3,,]/count[1:4,2,,]-1)
dimnames(change)[[2]]=c("R1","R2","R3")
dimnames(change)[[3]]=c(paste(cvthresh[1:7]/2,"-",cvthresh[2:8]/2),"<5")
names(dimnames(change))=c("Source","WRFversion","Threshold")

library(ggplot2)
library(reshape2)
tmp=melt(change[,,4:8])
print(ggplot(tmp, aes(x = Threshold, y = value, fill = WRFversion)) +
        geom_boxplot() +
        scale_fill_manual(values = c("red","blue","green"))  +
        scale_y_continuous(breaks=seq(-100, 150, 50)) +
        theme_bw() + ylab("Percentage change") + xlab("") +  geom_hline(yintercept = 0))


### Change in events for different event frequencies
## Function of resolution, wrf version, & choice of projection/radius

cvthresh=c(20,10,5,2)
cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")

cvBRAN<-cvNOEAC<-countBRAN<-countNOEAC<-array(0,c(8,3,4))

for(r in 1:3)
{
  n=1
  for(dom in c("d01","d02"))
    for(c in 1:4)
    {
      data2=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_BRAN_",cat[c],".csv",sep="")),
                  read.csv(paste("ECLevents_",dom,"_2008_R",r,"_BRAN_",cat[c],".csv",sep="")))
      data3=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")),
                  read.csv(paste("ECLevents_",dom,"_2008_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")))
      
      b1=order(data2$CV2,decreasing=T)
      b2=order(data3$CV2,decreasing=T)
      
      for(t in 1:4)
      {
        cvBRAN[n,r,t]=data2$CV2[b1[cvthresh[t]*2]]
        cvNOEAC[n,r,t]=data3$CV2[b2[cvthresh[t]*2]]
      }
      for(t in 1:4)
      {
        if(t<4)
        {
          countBRAN[n,r,t]=length(which(data2$CV2>=cvBRAN[n,r,t] & data2$CV2<cvBRAN[n,r,t+1]))
          countNOEAC[n,r,t]=length(which(data3$CV2>=cvBRAN[n,r,t] & data3$CV2<cvBRAN[n,r,t+1]))
        } else {
          countBRAN[n,r,t]=length(which(data2$CV2>=cvBRAN[n,r,t]))
          countNOEAC[n,r,t]=length(which(data3$CV2>=cvBRAN[n,r,t]))
        }
      }
      n=n+1 
    }
}




library(abind)
cvR2=abind(cv[1:4,2,],cvBRAN[1:4,2,],cv_noeac[1:4,2,],along=3)
dimnames(cvR2)[[1]]=cat
dimnames(cvR2)[[2]]=cvthresh
dimnames(cvR2)[[3]]=c("Default","BRAN","NoEAC")

### Okay, now loading individual datasets
cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
c=3
events<-eventsBRAN<-events_noeac<-events_2eac<-list()

n=1
for(dom in c("d01","d02"))
for(r in 1:3)
  {
    events[[n]]=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                      read.csv(paste("ECLevents_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
    events[[n]]$Year=floor(events[[n]]$Date1/10000)
    events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
    
    eventsBRAN[[n]]=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_BRAN_",cat[c],".csv",sep="")),
                      read.csv(paste("ECLevents_",dom,"_2008_R",r,"_BRAN_",cat[c],".csv",sep="")))
    eventsBRAN[[n]]$Year=floor(eventsBRAN[[n]]$Date1/10000)
    eventsBRAN[[n]]$Month=floor(eventsBRAN[[n]]$Date1/100)%%100
    
    events_noeac[[n]]=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")),
                          read.csv(paste("ECLevents_",dom,"_2008_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")))
    events_noeac[[n]]$Year=floor(events_noeac[[n]]$Date1/10000)
    events_noeac[[n]]$Month=floor(events_noeac[[n]]$Date1/100)%%100
    
    if(dom=="d01")
    {
    events_2eac[[n]]=rbind(read.csv(paste("ECLevents_",dom,"_2007_R",r,"_BRAN_2eac_",cat[c],".csv",sep="")),
                            read.csv(paste("ECLevents_",dom,"_2008_R",r,"_BRAN_2eac_",cat[c],".csv",sep="")))
    events_2eac[[n]]$Year=floor(events_2eac[[n]]$Date1/10000)
    events_2eac[[n]]$Month=floor(events_2eac[[n]]$Date1/100)%%100  
    }
    n=n+1     
}

count2<-matrix(0,6,3)
for(n in 1:6)
{
  count2[n,1]=length(which(eventsBRAN[[n]]$CV2>=2))
  count2[n,2]=length(which(events_noeac[[n]]$CV2>=2))
  if(n<4) count2[n,3]=length(which(events_2eac[[n]]$CV2>=2))
}


events_erai=read.csv("~/output/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv")
I=which(events_erai$Date2>=20070000 & events_erai$Date1<=20090000)
events_erai=events_erai[I,]

# for(n in 1:6) print(ks.test(eventsBRAN[[n]]$CV2,events_noeac[[n]]$CV2))
# 
# dens<-array(NaN,c(512,6,3))
# for(r in 1:6)
# {
#   a=density(events[[r]]$CV2,na.rm=T,from=1,to=5)
#   dens[,r,1]=a$y
#   a=density(eventsBRAN[[r]]$CV2,na.rm=T,from=1,to=5)
#   dens[,r,2]=a$y
#   a=density(events_noeac[[r]]$CV2,na.rm=T,from=1,to=5)
#   dens[,r,3]=a$y
# }

# plot(a$x,apply(dens[,1:3,1],1,mean),type="l",col="black",xlim=c(1,5),ylim=range(dens),
#      xlab="ECL intensity",ylab="Frequency",cex.main=1.2,lwd=3,main="")
# clist=c("black","blue","red")
# for(r in 2:3) lines(a$x,apply(dens[,1:3,r],1,mean),col=clist[r],lwd=3)
# legend("topright",legend=c("Default","BRAN SST","No EAC"),
#        col=clist,lwd=4,cex=1,bty="n")   

if(c<3) cvthresh=seq(1,2,0.125) else cvthresh=seq(1,4,0.25)
dens2<-array(0,c(12,6,4))
for(r in 1:6)
  for(i in 1:12)
  {
    I=which(events[[r]]$CV2>=cvthresh[i] & events[[r]]$CV2<cvthresh[i+1])
    dens2[i,r,1]=length(I)
    I=which(eventsBRAN[[r]]$CV2>=cvthresh[i] & eventsBRAN[[r]]$CV2<cvthresh[i+1])
    dens2[i,r,2]=length(I)
    I=which(events_noeac[[r]]$CV2>=cvthresh[i] & events_noeac[[r]]$CV2<cvthresh[i+1])
    dens2[i,r,3]=length(I)
    if(r<4)
    {
      I=which(events_2eac[[r]]$CV2>=cvthresh[i] & events_2eac[[r]]$CV2<cvthresh[i+1])
      dens2[i,r,4]=length(I)
    }
  }

densE=rep(0,12)
for(i in 1:12)
{
  I=which(events_erai$CV2>=cvthresh[i] & events_erai$CV2<cvthresh[i+1])
  densE[i]=length(I)
}

dimnames(dens2)[[1]]=cvthresh[1:12]
dimnames(dens2)[[2]]=c("R1 d01","R2 d01","R3 d01","R1 d02","R2 d02","R3 d02")

a=apply(dens2,c(2,3),sum)
mean(a[4:6,3]/a[4:6,1])

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_CVpdf_d01_SST_",cat[c],".pdf",sep=""),width=6,height=4,pointsize=12)
plot(cvthresh[1:12],apply(dens2[,1:3,1],1,mean),type="l",lwd=3,xlab="Intensity",ylab="Frequency",ylim=range(dens2))
lines(cvthresh[1:12],apply(dens2[,1:3,2],1,mean),lwd=3,col="grey")
lines(cvthresh[1:12],apply(dens2[,1:3,3],1,mean),lwd=3,col="blue")
lines(cvthresh[1:12],apply(dens2[,1:3,4],1,mean),lwd=3,col="red")
legend("topright",c("Control","BRAN","NoEAC","2EAC"),lwd=3,col=c("black","grey","blue","red"))
dev.off()

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_CVpdf_d01_SST_",cat[c],"_vERAI.pdf",sep=""),width=6,height=4,pointsize=12)
plot(cvthresh[1:12],densE,type="l",lwd=3,xlab="Intensity",ylab="Frequency",ylim=range(dens2))
lines(cvthresh[1:12],apply(dens2[,1:3,1],1,mean),lwd=3,col="red")
lines(cvthresh[1:12],apply(dens2[,1:3,2],1,mean),lwd=3,col="blue")
lines(cvthresh[1:12],apply(dens2[,4:6,1],1,mean),lwd=3,lty=2,col="blue")
lines(cvthresh[1:12],apply(dens2[,4:6,2],1,mean),lwd=3,lty=2,col="red")
legend("topright",c("ERAI","Control","BRAN","Control 10km","BRAN 10km"),lwd=3,col=c("black","blue","red","blue","red"))
dev.off()


dens3=array(0,c(12,6,4))
dens3[1,,]=dens2[1,,]
for(i in 2:12)
  dens3[i,,]=apply(dens2[1:i,,],c(2,3),sum)

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_CVcdf_d01_SST_",cat[c],".pdf",sep=""),width=6,height=4,pointsize=12)
plot(cvthresh[1:12],apply(dens3[,1:3,1],1,mean),type="l",lwd=3,xlab="Intensity",ylab="Cumulative total",ylim=range(dens3))
lines(cvthresh[1:12],apply(dens3[,1:3,2],1,mean),lwd=3,col="grey")
lines(cvthresh[1:12],apply(dens3[,1:3,3],1,mean),lwd=3,col="blue")
lines(cvthresh[1:12],apply(dens3[,1:3,4],1,mean),lwd=3,col="red")
legend("bottomright",c("Control","BRAN","NoEAC","2EAC"),lwd=3,col=c("black","grey","blue","red"))
dev.off()

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_CVpdf_change_d01_SST_",cat[c],".pdf",sep=""),width=6,height=4,pointsize=12)
plot(1,NA,xlim=c(1,max(cvthresh)),ylim=c(-max(cvthresh)*2,max(cvthresh)*2),xlab="Intensity",ylab="Frequency")
abline(h=0,col="red",lwd=2,lty=2)
for(r in 1:3) lines(cvthresh[1:12],dens2[,r,3]-dens2[,r,1],lwd=1,lty=2,col="blue")
for(r in 1:3) lines(cvthresh[1:12],dens2[,r,4]-dens2[,r,1],lwd=1,lty=2,col="red")
lines(cvthresh[1:12],apply(dens2[,1:3,3]-dens2[,1:3,1],1,mean),lwd=3,col="blue")
lines(cvthresh[1:12],apply(dens2[,1:3,4]-dens2[,1:3,1],1,mean),lwd=3,col="red")
legend("topright",c("NoEAC","2EAC"),lwd=3,col=c("blue","red"))
dev.off()

#######Trying for histogram

a=hist(events[[1]]$CV2,breaks=cvthresh,plot=F)
tmp=a
tmp$counts=apply(dens2[,,3]-dens2[,,2],1,mean)
plot(tmp,col=rgb(0,0,1,1/4),xlim=c(1,4),ylim=c(-4,4),
     xlab="Intensity (hPa/(deg.lat)^2)",ylab="Average change in events",main="")
tmp$counts=apply(dens2[,1:3,4]-dens2[,1:3,2],1,mean)
plot(tmp,col=rgb(1,0,0,1/4),add=T)
legend("topright",legend=c("NoEAC","2EAC"),
       col=c(rgb(0,0,1,1/4),col=rgb(1,0,0,1/4),"red"),lwd=8,cex=1,bty="n")   


### What about average length
lens=c(1:10,100)
len<-array(0,c(10,6,4))
for(i in 1:6)
for(j in 1:10)
{
  I=which(events[[i]]$Length2>=j & events[[i]]$Length2<lens[j+1])
  len[j,i,1]=length(I)
  I=which(eventsBRAN[[i]]$Length2>=j & eventsBRAN[[i]]$Length2<lens[j+1])
  len[j,i,2]=length(I)
  I=which(events_noeac[[i]]$Length2>=j & events_noeac[[i]]$Length2<lens[j+1])
  len[j,i,3]=length(I)
  if(i<4)
  {
    I=which(events_2eac[[i]]$Length2>=j & events_2eac[[i]]$Length2<lens[j+1])
    len[j,i,4]=length(I)
  }
}

########### ECL days

fixes<-fixesBRAN<-fixes_noeac<-fixes_2eac<-list()
c=3
n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    fixes[[n]]=rbind(read.csv(paste("ECLfixes_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                      read.csv(paste("ECLfixes_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    fixes[[n]]$Location2<-0
    I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
    fixes[[n]]$Location2[I]<-1
    
    fixesBRAN[[n]]=rbind(read.csv(paste("ECLfixes_",dom,"_2007_R",r,"_BRAN_",cat[c],".csv",sep="")),
                          read.csv(paste("ECLfixes_",dom,"_2008_R",r,"_BRAN_",cat[c],".csv",sep="")))
    fixesBRAN[[n]]$Year=floor(fixesBRAN[[n]]$Date/10000)
    fixesBRAN[[n]]$Month=floor(fixesBRAN[[n]]$Date/100)%%100
    fixesBRAN[[n]]$Location2<-0
    I<-which(fixesBRAN[[n]][,7]>=149 & fixesBRAN[[n]][,7]<=154 & fixesBRAN[[n]][,8]<(-37) & fixesBRAN[[n]][,8]>=-41)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=(149+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,7]<=(154+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,8]<(-31) & fixesBRAN[[n]][,8]>=-37)
    fixesBRAN[[n]]$Location2[I]<-1
    I<-which(fixesBRAN[[n]][,7]>=152 & fixesBRAN[[n]][,7]<=157 & fixesBRAN[[n]][,8]<=(-24) & fixesBRAN[[n]][,8]>=-31)
    fixesBRAN[[n]]$Location2[I]<-1
    
    fixes_noeac[[n]]=rbind(read.csv(paste("ECLfixes_",dom,"_2007_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")),
                            read.csv(paste("ECLfixes_",dom,"_2008_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")))
    fixes_noeac[[n]]$Year=floor(fixes_noeac[[n]]$Date/10000)
    fixes_noeac[[n]]$Month=floor(fixes_noeac[[n]]$Date/100)%%100
    fixes_noeac[[n]]$Location2<-0
    I<-which(fixes_noeac[[n]][,7]>=149 & fixes_noeac[[n]][,7]<=154 & fixes_noeac[[n]][,8]<(-37) & fixes_noeac[[n]][,8]>=-41)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=(149+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,7]<=(154+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,8]<(-31) & fixes_noeac[[n]][,8]>=-37)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=152 & fixes_noeac[[n]][,7]<=157 & fixes_noeac[[n]][,8]<=(-24) & fixes_noeac[[n]][,8]>=-31)
    fixes_noeac[[n]]$Location2[I]<-1
    
    if(dom=="d01")
    {
    fixes_2eac[[n]]=rbind(read.csv(paste("ECLfixes_",dom,"_2007_R",r,"_BRAN_2eac_",cat[c],".csv",sep="")),
                           read.csv(paste("ECLfixes_",dom,"_2008_R",r,"_BRAN_2eac_",cat[c],".csv",sep="")))
    fixes_2eac[[n]]$Year=floor(fixes_2eac[[n]]$Date/10000)
    fixes_2eac[[n]]$Month=floor(fixes_2eac[[n]]$Date/100)%%100
    fixes_2eac[[n]]$Location2<-0
    I<-which(fixes_2eac[[n]][,7]>=149 & fixes_2eac[[n]][,7]<=154 & fixes_2eac[[n]][,8]<(-37) & fixes_2eac[[n]][,8]>=-41)
    fixes_2eac[[n]]$Location2[I]<-1
    I<-which(fixes_2eac[[n]][,7]>=(149+(37+fixes_2eac[[n]][,8])/2) & fixes_2eac[[n]][,7]<=(154+(37+fixes_2eac[[n]][,8])/2) & fixes_2eac[[n]][,8]<(-31) & fixes_2eac[[n]][,8]>=-37)
    fixes_2eac[[n]]$Location2[I]<-1
    I<-which(fixes_2eac[[n]][,7]>=152 & fixes_2eac[[n]][,7]<=157 & fixes_2eac[[n]][,8]<=(-24) & fixes_2eac[[n]][,8]>=-31)
    fixes_2eac[[n]]$Location2[I]<-1
    }
    
    n=n+1     
  }

daylist=seq.Date(from=as.Date("2007-01-01"),to=as.Date("2008-12-31"),by=1)
daylist2=as.numeric(format.Date(daylist,format="%Y%m%d"))

dayCV<-array(NaN,c(length(daylist2),6,8))
dimnames(dayCV)[[2]]=c("R1","R2","R3","R1 d02","R2 d02","R3 d02")
dimnames(dayCV)[[3]]=c("Default","BRAN","No EAC","2EAC","Default close","BRAN close","No EAC close","2EAC close")
for(i in 1:length(daylist2))
for(j in 1:6)
{
  I=which(fixes[[j]]$Date==daylist2[i] & fixes[[j]]$Location==1)
  if(length(I)>0) dayCV[i,j,1]=max(fixes[[j]]$CV[I])
  I=which(fixes[[j]]$Date==daylist2[i] & fixes[[j]]$Location2==1)
  if(length(I)>0) dayCV[i,j,5]=max(fixes[[j]]$CV[I])
  
  I=which(fixesBRAN[[j]]$Date==daylist2[i] & fixesBRAN[[j]]$Location==1)
  if(length(I)>0) dayCV[i,j,2]=max(fixesBRAN[[j]]$CV[I])
  I=which(fixesBRAN[[j]]$Date==daylist2[i] & fixesBRAN[[j]]$Location2==1)
  if(length(I)>0) dayCV[i,j,6]=max(fixesBRAN[[j]]$CV[I])
  
  I=which(fixes_noeac[[j]]$Date==daylist2[i] & fixes_noeac[[j]]$Location==1)
  if(length(I)>0) dayCV[i,j,3]=max(fixes_noeac[[j]]$CV[I])
  I=which(fixes_noeac[[j]]$Date==daylist2[i] & fixes_noeac[[j]]$Location2==1)
  if(length(I)>0) dayCV[i,j,7]=max(fixes_noeac[[j]]$CV[I])
  
  if(j<4)
  {
    I=which(fixes_2eac[[j]]$Date==daylist2[i] & fixes_2eac[[j]]$Location==1)
    if(length(I)>0) dayCV[i,j,4]=max(fixes_noeac[[j]]$CV[I])
    I=which(fixes_2eac[[j]]$Date==daylist2[i] & fixes_2eac[[j]]$Location2==1)
    if(length(I)>0) dayCV[i,j,8]=max(fixes_noeac[[j]]$CV[I])
    
  }
}

fixes_erai=read.csv("~/output/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv")
I=which(fixes_erai$Date>=20070000 & fixes_erai$Date<=20090000)
fixes_erai=fixes_erai[I,]

dayCV_erai=rep(NaN,length(daylist))
for(i in 1:length(daylist2))
  {
    I=which(fixes_erai$Date==daylist2[i] & fixes_erai$Location==1)
    if(length(I)>0) dayCV_erai[i]=max(fixes_erai$CV[I])
}

dayECL=!is.na(dayCV)
dayECL2=apply(dayECL[,1:3,],c(1,3),sum)
dayECL_erai=!is.na(dayCV_erai)
dayECL3=cbind(dayECL_erai,dayECL2)
skillmatrix<-array(NaN,c(6,4,3))

for(i in 1:6)
  for(j in 1:4)
{
  H=length(which(dayECL_erai==1 & dayECL[,i,j]==1))
  M=length(which(dayECL_erai==1 & dayECL[,i,j]==0))
  FA=length(which(dayECL_erai==0 & dayECL[,i,j]==1))
  
  skillmatrix[i,j,1]=H/(H+M)
  skillmatrix[i,j,2]=FA/(H+FA)
  skillmatrix[i,j,3]=H/(H+M+FA)
  }

skillmatrix[4:6,4,]=NaN

CSI<-array(NaN,c(6,6,4))
for(i in 1:6)
  for(j in 1:6)
  for(k in 1:4)
  {
    H=length(which(dayECL[,i,k]==1 & dayECL[,j,k]==1))
    M=length(which(dayECL[,i,k]==1 & dayECL[,j,k]==0))
    FA=length(which(dayECL[,i,k]==0 & dayECL[,j,k]==1))
    CSI[i,j,k]=H/(H+M+FA)
  }
CSI[CSI>=1]=NaN
apply(CSI,3,mean,na.rm=T)

dimnames(skillmatrix)[[1]]=c("R1","R2","R3","R1 d02","R2 d02","R3 d02")
dimnames(skillmatrix)[[2]]=c("Default","BRAN","No EAC","2EAC")
dimnames(skillmatrix)[[3]]=c("Hit Rate","False Alarm Rate","CSI")
apply(skillmatrix,c(2,3),mean,na.rm=T)


plot(dayCV[,1,2],dayCV[,1,3],type="p",col="blue",pch=4,lwd=2,
     xlab="BRAN",ylab="No EAC",xlim=c(1,3.5),ylim=c(1,3.5),main="Daily max CV")
points(dayCV[,2,2],dayCV[,2,3],col="red",pch=4,lwd=2)
points(dayCV[,3,2],dayCV[,3,3],col="green",pch=4,lwd=2)
abline(a=0,b=1,col="black")
legend("topleft",c("R1","R2","R3"),col=c("blue","red","green"),pch=4)

apply(dayCV[,,3]-dayCV[,,2],2,mean,na.rm=T)

# lat=seq(-40,-25,2.5)
# lon=seq(145,160,2.5)
# loc<-array(0,c(7,7,6,4))
lat=seq(-40,-25,5)
lon=seq(145,160,5)
loc<-array(NaN,c(4,4,6,4))

for(i in 1:4)
  for(j in 1:4)
    for(k in 1:6)
  {
    I=which(fixes[[k]]$Lat>=lat[i]-2.5 & fixes[[k]]$Lat<lat[i]+2.5 & fixes[[k]]$Lon>=lon[j]-2.5 & fixes[[k]]$Lon<lon[j]+2.5 & fixes[[k]]$Location==1)
    loc[i,j,k,1]=length(I)
    I=which(fixesBRAN[[k]]$Lat>=lat[i]-2.5 & fixesBRAN[[k]]$Lat<lat[i]+2.5 & fixesBRAN[[k]]$Lon>=lon[j]-2.5 & fixesBRAN[[k]]$Lon<lon[j]+2.5 & fixesBRAN[[k]]$Location==1)
    loc[i,j,k,2]=length(I)
    I=which(fixes_noeac[[k]]$Lat>=lat[i]-2.5 & fixes_noeac[[k]]$Lat<lat[i]+2.5 & fixes_noeac[[k]]$Lon>=lon[j]-2.5 & fixes_noeac[[k]]$Lon<lon[j]+2.5 & fixes_noeac[[k]]$Location==1)
    loc[i,j,k,3]=length(I)
    if(k<4)
    {
    I=which(fixes_2eac[[k]]$Lat>=lat[i]-2.5 & fixes_2eac[[k]]$Lat<lat[i]+2.5 & fixes_2eac[[k]]$Lon>=lon[j]-2.5 & fixes_2eac[[k]]$Lon<lon[j]+2.5 & fixes_2eac[[k]]$Location==1)
    loc[i,j,k,4]=length(I)
    }
    }

loc_erai=matrix(0,7,7)
for(i in 1:7)
  for(j in 1:7)
    {
      I=which(fixes_erai$Lat>=lat[i]-1.25 & fixes_erai$Lat<lat[i]+1.25 & fixes_erai$Lon>=lon[j]-1.25 & fixes_erai$Lon<lon[j]+1.25)
      loc_erai[i,j]=length(I)
    }

cols=gray(seq(1,0.1,-0.15))
ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(1.5, length(brks) - 1.5, subsampleg), tick = TRUE, 
       labels = brks[seq(2, length(brks)-1, subsampleg)])
}

bb=c(-0.5,0,1,2,4,8,12,100)
names=c("Default","BRAN","No EAC","2EAC")
library("R.matlab")
readMat('~/Documents/GDI/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_d01_",cat[c],"_SST_vERAI.pdf",sep=""),width=11,height=4,pointsize=12)
layout(cbind(1,2,3,4),width=c(1,1,1,0.3))
image(lon,lat,t(loc_erai)/2,xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="ERAI",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat,t(apply(loc[,,1:6,1],c(1,2),mean))/2,xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="Control",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat,t(apply(loc[,,1:6,2],c(1,2),mean))/2,xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="BRAN",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(bb,cols)
dev.off()

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_d01_",cat[c],"_SST_v3.pdf",sep=""),width=8,height=8,pointsize=12)
layout(cbind(c(1,3),c(2,4),c(5,5)),width=c(1,1,0.3))
for(r in 1:4)
{
image(lon,lat,t(apply(loc[,,1:3,r],c(1,2),mean))/2,xlab="",ylab="",breaks=bb,col=cols,zlim=c(0,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=names[r],cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
}
ColorBar(bb,cols)
dev.off()

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
#bb2=c(-100,-5,-2,-1,0,1,2,5,100)
#bb2=c(-100,-10,-5,-2,-1,0,1,2,5,10,100)
bb2=c(-10000,-50,-25,-10,-5,0,5,10,25,50,10000)
cm=pal(10)

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_",cat[c],"_SST_change_vwrf_v3.pdf",sep=""),width=12,height=4,pointsize=12)
layout(cbind(1,2,3,4),c(1,1,1,0.3))
for(r in 1:3)
{
  image(lon,lat,100*t((loc[,,r,3]/loc[,,r,2])-1),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=paste("R",r,sep=""),cex.axis=1.5,cex.main=1.5)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
}
ColorBar(bb2,cm)
dev.off()

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_",cat[c],"_SST_change_vres_v3.pdf",sep=""),width=9,height=4,pointsize=12)
layout(cbind(1,2,3),c(1,1,0.3))
image(lon,lat,100*t(apply(((loc[,,1:3,3]/loc[,,1:3,2])-1),c(1,2),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="50 km",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat,100*t(apply(((loc[,,4:6,3]/loc[,,4:6,2])-1),c(1,2),mean)),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="10 km",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_",cat[c],"_SST_change_vbran_d01.pdf",sep=""),width=8,height=4,pointsize=12)
bb2=c(-10000,-5,-2.5,-1,-0.5,0,0.5,1,2.5,5,10000)
layout(cbind(1,2,3),c(1,1,0.3))
image(lon,lat,t(apply(loc[,,1:3,3]-loc[,,1:3,2],c(1,2),mean))/2,xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="No EAC",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat,t(apply(loc[,,1:3,4]-loc[,,1:3,2],c(1,2),mean))/2,xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="2EAC",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(bb2,cm)
dev.off()

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_",cat[c],"_SST_change_vbran_d01_locpc.pdf",sep=""),width=8,height=4,pointsize=12)
bb2=c(-10000,-40,-30,-20,-10,0,10,20,30,40,10000)
layout(cbind(1,2,3),c(1,1,0.3))
image(lon,lat,t(apply((loc[,,1:3,3]/loc[,,1:3,2])-1,c(1,2),mean))*100,xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="No EAC",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
image(lon,lat,t(apply((loc[,,1:3,4]/loc[,,1:3,2])-1,c(1,2),mean))*100,xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main="2EAC",cex.axis=1.5,cex.main=1.5)
contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
ColorBar(bb2,cm)
dev.off()


# Comparing location change for the two important radii
library("R.matlab")
readMat('~/Documents/Data/Useful.mat')->Useful
mask<-t(Useful$mask)
mask[is.na(mask)]=0

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
bb2=c(-100,-5,-2,-1,0,1,2,5,100)
cm=pal(8)
rad=c("500 km","500 km","200 km","200 km")
lat=seq(-40,-25,5)
lon=seq(145,160,5)

loc<-array(0,c(4,4,3,4,3))
names(dimnames(loc))<-c("Lat","Lon","R","Cat","Source")

for(c in 1:4)
{
  fixes<-fixesBRAN<-fixes_noeac<-list()
  n=1
  for(dom in c("d01"))
    for(r in 1:3)
    {
      fixes[[n]]=rbind(read.csv(paste("ECLfixes_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                       read.csv(paste("ECLfixes_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
      fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
      fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
      fixes[[n]]$Location2<-0
      I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
      fixes[[n]]$Location2[I]<-1
      I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
      fixes[[n]]$Location2[I]<-1
      I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
      fixes[[n]]$Location2[I]<-1
      
      fixesBRAN[[n]]=rbind(read.csv(paste("ECLfixes_",dom,"_2007_R",r,"_BRAN_",cat[c],".csv",sep="")),
                           read.csv(paste("ECLfixes_",dom,"_2008_R",r,"_BRAN_",cat[c],".csv",sep="")))
      fixesBRAN[[n]]$Year=floor(fixesBRAN[[n]]$Date/10000)
      fixesBRAN[[n]]$Month=floor(fixesBRAN[[n]]$Date/100)%%100
      fixesBRAN[[n]]$Location2<-0
      I<-which(fixesBRAN[[n]][,7]>=149 & fixesBRAN[[n]][,7]<=154 & fixesBRAN[[n]][,8]<(-37) & fixesBRAN[[n]][,8]>=-41)
      fixesBRAN[[n]]$Location2[I]<-1
      I<-which(fixesBRAN[[n]][,7]>=(149+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,7]<=(154+(37+fixesBRAN[[n]][,8])/2) & fixesBRAN[[n]][,8]<(-31) & fixesBRAN[[n]][,8]>=-37)
      fixesBRAN[[n]]$Location2[I]<-1
      I<-which(fixesBRAN[[n]][,7]>=152 & fixesBRAN[[n]][,7]<=157 & fixesBRAN[[n]][,8]<=(-24) & fixesBRAN[[n]][,8]>=-31)
      fixesBRAN[[n]]$Location2[I]<-1
      
      fixes_noeac[[n]]=rbind(read.csv(paste("ECLfixes_",dom,"_2007_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")),
                             read.csv(paste("ECLfixes_",dom,"_2008_R",r,"_BRAN_noeac_",cat[c],".csv",sep="")))
      fixes_noeac[[n]]$Year=floor(fixes_noeac[[n]]$Date/10000)
      fixes_noeac[[n]]$Month=floor(fixes_noeac[[n]]$Date/100)%%100
      fixes_noeac[[n]]$Location2<-0
      I<-which(fixes_noeac[[n]][,7]>=149 & fixes_noeac[[n]][,7]<=154 & fixes_noeac[[n]][,8]<(-37) & fixes_noeac[[n]][,8]>=-41)
      fixes_noeac[[n]]$Location2[I]<-1
      I<-which(fixes_noeac[[n]][,7]>=(149+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,7]<=(154+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,8]<(-31) & fixes_noeac[[n]][,8]>=-37)
      fixes_noeac[[n]]$Location2[I]<-1
      I<-which(fixes_noeac[[n]][,7]>=152 & fixes_noeac[[n]][,7]<=157 & fixes_noeac[[n]][,8]<=(-24) & fixes_noeac[[n]][,8]>=-31)
      fixes_noeac[[n]]$Location2[I]<-1
      
      for(i in 1:4)
        for(j in 1:4)
          {
            I=which(fixes[[n]]$Lat>=lat[i]-2.5 & fixes[[n]]$Lat<lat[i]+2.5 & fixes[[n]]$Lon>=lon[j]-2.5 & fixes[[n]]$Lon<lon[j]+2.5)
            loc[i,j,r,c,1]=length(I)
            I=which(fixesBRAN[[n]]$Lat>=lat[i]-2.5 & fixesBRAN[[n]]$Lat<lat[i]+2.5 & fixesBRAN[[n]]$Lon>=lon[j]-2.5 & fixesBRAN[[n]]$Lon<lon[j]+2.5)
            loc[i,j,r,c,2]=length(I)
            I=which(fixes_noeac[[n]]$Lat>=lat[i]-2.5 & fixes_noeac[[n]]$Lat<lat[i]+2.5 & fixes_noeac[[n]]$Lon>=lon[j]-2.5 & fixes_noeac[[n]]$Lon<lon[j]+2.5)
            loc[i,j,r,c,3]=length(I)
        }
      
      n=n+1
    }
}

change=apply(loc[,,,,3]-loc[,,,,2],c(1,2,4),mean)/2

pdf(file=paste("~/Documents/ECLs/WRFruns/0708/ECL_locations_0708_SST_change_byrad_v2.pdf",sep=""),width=9,height=4,pointsize=12)
layout(cbind(1,2,3),c(1,1,0.3))
for(i in c(1,3))
{
  image(lon,lat,t(change[,,i]),xlab="",ylab="",breaks=bb2,col=cm,zlim=c(-Inf,Inf),xlim=c(142.5,162.5),ylim=c(-42.5,-22.5),main=rad[i],cex.axis=1.5,cex.main=1.5)
  contour(Useful$x,Useful$y,mask,add=T,drawlabels=F)
  
}
ColorBar(bb2,cm)
dev.off()

######### Rain composites

rm(list=ls())
setwd('~/Documents/ECLs/WRFruns/0708/')

library(RNetCDF)
library(fields)
library(akima)
library("R.matlab")
readMat('~/Documents/GDI/Useful_ECL.mat')->Useful
readMat('~/Documents/GDI/Useful.mat')->Useful2
mask<-t(Useful$mask)
mask[is.na(mask)]=0
source('~/R/x86_64-pc-linux-gnu-library/2.15/gplots/R/rich.colors.R')
wrfv=c("R1","R2","R3")

cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
c=3
dir='/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/typing/'
events<-fixes<-comp<-list()
comp_d01<-array(0,c(21,21,3,4)) # X, Y, R, EAC type
comp_d02<-array(0,c(101,101,3,4))

tlist=c("","_BRAN","_BRAN_noeac","_BRAN_2eac")
name2=c("","","","_v2")
dir1=c(36,37,37,45)

n=1
for(dom in c("d01"))
  for(type in 1:4)
  for(r in 1:3)
  {
    events[[n]]=read.csv(paste(dir,"ECLevents_",dom,"_0708_R",r,tlist[type],"_",cat[c],name2[type],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    fixes[[n]]=read.csv(paste(dir,"ECLfixes_",dom,"_0708_R",r,tlist[type],"_",cat[c],name2[type],"_typing_impactsC.csv",sep=""),stringsAsFactors = F)
    
    fixes[[n]]$Location2<-0
    I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
    fixes[[n]]$Location2[I]<-1
    I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
    fixes[[n]]$Location2[I]<-1
    
    ### Rain stuff
    
    if(dom=="d01") a=open.nc(paste("/srv/ccrc/data",dir1[type],"/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",tlist[type],"/out/ECLrain_0708_",cat[c],name2[type],"_centred.nc",sep="")) else
      a=open.nc(paste("/srv/ccrc/data",dir1[type],"/z3478332/WRF/output/ERAI_R",r,"_nudging_default_2007",tlist[type],"/out/ECLrain_d02_0708_",cat[c],name2[type],".nc",sep=""))
    tmp=var.get.nc(a,"ECLrain")
    fixes[[n]]$MeanRain=apply(tmp,3,mean,na.rm=T)
    fixes[[n]]$MaxRain=apply(tmp,3,max,na.rm=T)
    I=which(fixes[[n]]$Location2==1)
    if(dom=="d01") {
      fixes[[n]]$MeanRainS=apply(tmp[,1:10,],3,mean,na.rm=T)
      fixes[[n]]$MeanRainN=apply(tmp[,12:21,],3,mean,na.rm=T)
      comp_d01[,,r,type]=apply(tmp[,,I],c(1,2),mean,na.rm=T) 
    } else {
      fixes[[n]]$MeanRainS=apply(tmp[,1:50,],3,mean,na.rm=T)
      fixes[[n]]$MeanRainN=apply(tmp[,52:101,],3,mean,na.rm=T)
      comp_d02[,,r,type]=apply(tmp[,,I],c(1,2),mean,na.rm=T)
    }
    #comp[[n]]=tmp
    n=n+1     
  }

source('~/Documents/R/color.palette.R')
pal <- color.palette(c("red","yellow","white","cyan","blue"), c(20,20,20,20))
#bb2=c(-100,-4,-3,-2,-1,0,1,2,3,4,100)
bb2=c(-100,-5,-2,-1,0,1,2,5,100)
layout(cbind(1,2,3))
par(mar=c(1,1,3,1))
for(i in 1:3) image(comp_d02[,,i,3]-comp_d02[,,i,2],col=pal(8),breaks=bb2,main=paste("RCM",i),axes=F)

image(apply(comp_d01[,,,4]-comp_d01[,,,2],c(1,2),mean),col=pal(8),breaks=bb2,axes=F)

ColorBar <- function(brks,cols,vert=T,subsampleg=1)
{
  par(mar = c(3, 1, 3, 3), mgp = c(1, 1, 0), las = 1, cex = 1)
  image(1, c(1:length(cols)), t(c(1:length(cols))), axes = FALSE, col = cols, 
        xlab = '', ylab = '')
  box()
  axis(4, at = seq(0.5, length(brks) - 0.5, subsampleg), tick = TRUE, 
       labels = brks[seq(1, length(brks), subsampleg)])
}
pdf(file="ECLrain_d01_centred.pdf",width=12,height=4)
layout(cbind(1,2,3,4),width=c(1,1,1,0.3))
image(apply(comp_d01[,,,3],c(1,2),mean),breaks=c(seq(0,14,1),100),col=tim.colors(15),main="NoEAC",axes=F,cex.main=2)
image(apply(comp_d01[,,,2],c(1,2),mean),breaks=c(seq(0,14,1),100),col=tim.colors(15),main="BRAN",axes=F,cex.main=2)
image(apply(comp_d01[,,,4],c(1,2),mean),breaks=c(seq(0,14,1),100),col=tim.colors(15),main="2EAC",axes=F,cex.main=2)
ColorBar(c(seq(0,14,1),100),tim.colors(15),subsampleg=2)
dev.off()

###########

##Okay, stuff about pdf of rain rates etc.

pdf_d02<-array(0,c(512,3,3,4)) # X, R, Type, c(FixMean,FixMax,EventMean,EventMax)
x_d02<-matrix(0,512,4)

n=1
for(t in 1:3)
   for(r in 1:3)
   {
     data=fixes[[n]]
     a=density(data$MeanRain,from=0,to=20)
     if(n==1) x_d02[,1]=a$x
     pdf_d02[,r,t,1]=a$y
     
     a=density(data$MaxRain,from=0,to=512)
     if(n==1) x_d02[,2]=a$x
     pdf_d02[,r,t,2]=a$y
     
     tmp1=data$ID+floor(data$Date/10000)*100
     tmp=unique(data$ID+floor(data$Date/10000)*100)
     data2=data.frame(ID=tmp,MeanRain=rep(0,length(tmp)),MaxRain=rep(0,length(tmp)))
     for(i in 1:length(tmp))
     {
       I=which(tmp1==tmp[i])
       data2[i,2]=max(data$MeanRain[I])
       data2[i,3]=max(data$MaxRain[I])}
     
     a=density(data2$MeanRain,from=0,to=20)
     if(n==1) x_d02[,3]=a$x
     pdf_d02[,r,t,3]=a$y
     
     a=density(data2$MaxRain,from=0,to=512)
     if(n==1) x_d02[,4]=a$x
     pdf_d02[,r,t,4]=a$y
     
     n=n+1
   }


clist=c("black","grey","blue")
plot(x_d02[,4],apply(pdf_d02[,,,4],1,max),
     type="l",col="white",lwd=3,xlab="Max Rainfall (mm)",ylab="Frequency",main="PDF of maximum 6hourly rainfall within 500 km of low")
for(i in 1:3)
  lines(x_d02[,4],apply(pdf_d02[,,i,4],1,mean),col=clist[i],lwd=3)
legend("topright",c("Control","BRAN","NoEAC"),lwd=3,col=clist)

stats<-array(0,c(3,3,3))
n=1
for(t in 1:3)
  for(r in 1:3)
  {
    data=fixes[[n]]
    tmp1=data$ID+floor(data$Date/10000)*100
    tmp=unique(data$ID+floor(data$Date/10000)*100)
    data2=data.frame(ID=tmp,MeanRain=rep(0,length(tmp)),MaxRain=rep(0,length(tmp)))
    for(i in 1:length(tmp))
    {
      I=which(tmp1==tmp[i])
      data2[i,2]=max(data$MeanRain[I])
      data2[i,3]=max(data$MaxRain[I])}
    
    stats[r,t,1]=length(tmp)
    I=which(data2$MeanRain>=6)
    stats[r,t,2]=length(I)
    I=which(data2$MaxRain>=50)
    stats[r,t,3]=length(I)
    n=n+1
  }

################## TEST: slp/gv composites


library(RNetCDF)
a=open.nc("/srv/ccrc/data36/z3478332/WRF/output/ERAI_R1_nudging_default_2007/out/ECLslp_0708_rad2_p100.nc")
ECLs=fixes[[1]]
tmp=var.get.nc(a,"ECL_slp")
I=which(ECLs$Location==1 & ECLs$CV>=2)
contour(apply(tmp[,,I],c(1,2),mean,na.rm=T))
