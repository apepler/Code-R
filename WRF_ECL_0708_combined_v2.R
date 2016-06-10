rm(list=ls())
setwd('/home/nfs/z3478332/output/outputUM_wrf_2007_all/')

year=c(2007,2008)
domain=c("d01","d02")
cat=c("rad5_p100","rad5_p240","rad2_p100","rad2_p240")
count<-countBRAN<-count_noeac<-matrix(0,4,3)

fixes<-fixesBRAN<-fixes_noeac<-fixes_notopo<-list()
c=3
n=1
for(dom in c("d01","d02"))
  for(r in 1:3)
  {
    fixes[[n]]=rbind(read.csv(paste("ECLfixes_",dom,"_2007_R",r,"_",cat[c],".csv",sep="")),
                     read.csv(paste("ECLfixes_",dom,"_2008_R",r,"_",cat[c],".csv",sep="")))
    fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
    fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
    fixes[[n]]$ID=fixes[[n]]$Year*1000+fixes[[n]]$ID
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
    fixesBRAN[[n]]$ID=fixesBRAN[[n]]$Year*1000+fixesBRAN[[n]]$ID
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
    fixes_noeac[[n]]$ID=fixes_noeac[[n]]$Year*1000+fixes_noeac[[n]]$ID
    fixes_noeac[[n]]$Location2<-0
    I<-which(fixes_noeac[[n]][,7]>=149 & fixes_noeac[[n]][,7]<=154 & fixes_noeac[[n]][,8]<(-37) & fixes_noeac[[n]][,8]>=-41)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=(149+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,7]<=(154+(37+fixes_noeac[[n]][,8])/2) & fixes_noeac[[n]][,8]<(-31) & fixes_noeac[[n]][,8]>=-37)
    fixes_noeac[[n]]$Location2[I]<-1
    I<-which(fixes_noeac[[n]][,7]>=152 & fixes_noeac[[n]][,7]<=157 & fixes_noeac[[n]][,8]<=(-24) & fixes_noeac[[n]][,8]>=-31)
    fixes_noeac[[n]]$Location2[I]<-1
    
    fixes_notopo[[n]]=rbind(read.csv(paste("ECLfixes_",dom,"_2007_R",r,"_notopo_",cat[c],".csv",sep="")),
                           read.csv(paste("ECLfixes_",dom,"_2008_R",r,"_notopo_",cat[c],".csv",sep="")))
    fixes_notopo[[n]]$Year=floor(fixes_notopo[[n]]$Date/10000)
    fixes_notopo[[n]]$Month=floor(fixes_notopo[[n]]$Date/100)%%100
    fixes_notopo[[n]]$ID=fixes_notopo[[n]]$Year*1000+fixes_notopo[[n]]$ID
    fixes_notopo[[n]]$Location2<-0
    I<-which(fixes_notopo[[n]][,7]>=149 & fixes_notopo[[n]][,7]<=154 & fixes_notopo[[n]][,8]<(-37) & fixes_notopo[[n]][,8]>=-41)
    fixes_notopo[[n]]$Location2[I]<-1
    I<-which(fixes_notopo[[n]][,7]>=(149+(37+fixes_notopo[[n]][,8])/2) & fixes_notopo[[n]][,7]<=(154+(37+fixes_notopo[[n]][,8])/2) & fixes_notopo[[n]][,8]<(-31) & fixes_notopo[[n]][,8]>=-37)
    fixes_notopo[[n]]$Location2[I]<-1
    I<-which(fixes_notopo[[n]][,7]>=152 & fixes_notopo[[n]][,7]<=157 & fixes_notopo[[n]][,8]<=(-24) & fixes_notopo[[n]][,8]>=-31)
    fixes_notopo[[n]]$Location2[I]<-1
    
    n=n+1     
  }

### Time spent in location

timelist=c(0:8,Inf)
loctime<-array(0,c(9,6,4,2))

for(i in 1:6)
  {
    IDs=unique(fixes[[i]]$ID)
    len=matrix(0,length(IDs),2)
    for(j in 1:length(IDs))
      {
      I=which(fixes[[i]]$ID==IDs[j])
      len[j,1]=sum(fixes[[i]]$Location[I])
      len[j,2]=sum(fixes[[i]]$Location2[I])
    }
    for(j in 1:9)
      for(k in 1:2)
        loctime[j,i,1,k]=length(which(len[,k]>=timelist[j] & len[,k]<timelist[j+1]))
    
    IDs=unique(fixes_notopo[[i]]$ID)
    len=matrix(0,length(IDs),2)
    for(j in 1:length(IDs))
    {
      I=which(fixes_notopo[[i]]$ID==IDs[j])
      len[j,1]=sum(fixes_notopo[[i]]$Location[I])
      len[j,2]=sum(fixes_notopo[[i]]$Location2[I])
    }
    for(j in 1:9)
      for(k in 1:2)
        loctime[j,i,2,k]=length(which(len[,k]>=timelist[j] & len[,k]<timelist[j+1]))
    
    IDs=unique(fixesBRAN[[i]]$ID)
    len=matrix(0,length(IDs),2)
    for(j in 1:length(IDs))
    {
      I=which(fixesBRAN[[i]]$ID==IDs[j])
      len[j,1]=sum(fixesBRAN[[i]]$Location[I])
      len[j,2]=sum(fixesBRAN[[i]]$Location2[I])
    }
    for(j in 1:9)
      for(k in 1:2)
        loctime[j,i,3,k]=length(which(len[,k]>=timelist[j] & len[,k]<timelist[j+1]))
    
    IDs=unique(fixes_noeac[[i]]$ID)
    len=matrix(0,length(IDs),2)
    for(j in 1:length(IDs))
    {
      I=which(fixes_noeac[[i]]$ID==IDs[j])
      len[j,1]=sum(fixes_noeac[[i]]$Location[I])
      len[j,2]=sum(fixes_noeac[[i]]$Location2[I])
    }
    for(j in 1:9)
      for(k in 1:2)
        loctime[j,i,4,k]=length(which(len[,k]>=timelist[j] & len[,k]<timelist[j+1]))
}

tmp=apply(loctime[,1:3,,],c(1,3,4),mean)
plot((1:8)*6,tmp[2:9,1,1],type="l",lwd=3,col="black",ylim=range(tmp[,,1]),xlab="Duration (hours)",ylab="Count")
for(i in 2:4) lines((1:8)*6,tmp[2:9,i,1],lwd=3,col=i)
legend("topright",legend=c("Control","NoTopo","BRAN","NoEAC"),lwd=3,col=1:4)

plot((1:8)*6,tmp[2:9,1,2],type="l",lwd=3,col="black",ylim=range(tmp[2:9,,2]),xlab="Duration (hours)",ylab="Count")
for(i in 2:4) lines((1:8)*6,tmp[2:9,i,2],lwd=3,col=i)
legend("topright",legend=c("Control","NoTopo","BRAN","NoEAC"),lwd=3,col=1:4)
