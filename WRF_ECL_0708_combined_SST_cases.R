###### New ECL days - using the cv0.5 runs - separating out if day by each threshold
rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/outputUM_wrf_2007_all/')
library(abind)

year=c(2007,2008)
domain=c("d01","d02")
cat=c("rad5_p100","rad5_p240","rad2_p100_cv0.5","rad2_p240")

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

dayECL=!is.na(dayCV)
dayECLa=apply(dayECL[,1:3,],c(1,3),sum)

dayECL_CV1=(!is.na(dayCV) & dayCV>=1)
dayECL_CV1a=apply(dayECL_CV1[,1:3,],c(1,3),sum)

##Matching w/in BRAN

missECL=array(0,c(length(daylist2),6,2))
for(i in 1:6)
{
  I=which(dayECL_CV1[,i,2]==1 & dayECL_CV1[,i,3]==0)
  missECL[I,i,1]=1
  I=which(dayECL_CV1[,i,2]==0 & dayECL_CV1[,i,3]==1)
  missECL[I,i,2]=1
}

missingdays=matrix(0,6,4)
for(i in 1:6)
{
I=which(dayECL_CV1[,i,2]==1 & dayECL_CV1[,i,3]==0)
a=dayCV[I,i,3]
missingdays[i,1]=length(a)
missingdays[i,2]=length(which(is.na(a)))

I=which(dayECL_CV1[,i,2]==0 & dayECL_CV1[,i,3]==1)
a=dayCV[I,i,2]
missingdays[i,3]=length(a)
missingdays[i,4]=length(which(is.na(a)))
}

clist=c("black","red","blue")
plot(dayCV[,1,2],dayCV[,1,3]-dayCV[,1,2],pch=4,col=clist[1],lwd=3,cex=1.5,
     xlab="Intensity in Control (hPa/deg.lat^2)",ylab="Difference in intensity between NoEAC and Control")
for(i in 2:3) points(dayCV[,i,2],dayCV[,i,3]-dayCV[,i,2],pch=4,col=clist[i],lwd=3,cex=1.5)

###########
## The other way - comparing events

