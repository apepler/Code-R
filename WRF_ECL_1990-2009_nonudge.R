rm(list=ls())
setwd('/srv/ccrc/data34/z3478332/ECLtracks/')
figdir="~/Documents/ECLs/WRFruns/0708/NoTopoPaper/"
source("~/Documents/R/ECL_functions.R")
dirs=c("ERA-nudge","ERA-nonudge","ERA-nonudge_notopo")
cat="p100_rad2cv1"
cat2="rad2_p100_cv1.0"

events<-fixes<-list()

for(n in 1:3)
{
  events[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,"/ECLevents_",dirs[n],"_",cat2,".csv",sep=""))
  events[[n]]$Year=floor(events[[n]]$Date1/10000)
  events[[n]]$Month=floor(events[[n]]$Date1/100)%%100
  events[[n]][events[[n]]==-Inf]=NaN
  
  fixes[[n]]<-read.csv(paste("outputUM_",dirs[n],"/",cat,"/ECLfixes_",dirs[n],"_",cat2,".csv",sep=""))
  fixes[[n]]$Year=floor(fixes[[n]]$Date/10000)
  fixes[[n]]$Month=floor(fixes[[n]]$Date/100)%%100
  fixes[[n]]$Date2=as.POSIXct(paste(as.character(fixes[[n]]$Date),substr(fixes[[n]]$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")
}

years=1980:2009
count=matrix(NaN,30,3)

for(i in 1:30)
  for(j in 1:3)
{
  I=which(events[[j]]$Year==years[i])
  count[i,j]=length(I)
  }
count[1:10,3]=NaN
apply(count[11:30,],2,mean,na.rm=T)

cor(count[11:30,2],count[11:30,3])
cor(count[11:30,1],count[11:30,3])

ks.test(events[[1]]$CV2,events[[3]]$CV2) # Not significantly different

cvthresh=c(seq(1,4,0.5),NaN)
cvcount=array(0,c(7,3))
for(x in 1:7)
  for(j in 1:3)
    cvcount[x,j]=length(which(events[[j]]$CV2>=cvthresh[x] & events[[j]]$CV2<cvthresh[x+1] &
                              events[[j]]$Year>=1990))

makePDF(events[[2]]$CV2[events[[2]]$Year>=1990],events[[3]]$CV2)

######### Matching

match=eventmatch(events[[2]][events[[2]]$Year>=1990,],fixes[[2]][fixes[[2]]$Year>=1990,],events[[3]],fixes[[3]])
match2=match[,6:8]-match[,1:3]
apply(match2,2,mean,na.rm=T)
plot(match[,1],match2[,1])
abline(h=0,col="red")

###### ERAI comp

erai_E=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv")
erai_E$Year=floor(erai_E$Date1/10000)
erai_E$Month=floor(erai_E$Date1/100) %% 100
erai=read.csv("/srv/ccrc/data34/z3478332/ECLtracks/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv")
erai$Year=floor(erai$Date/10000)
erai$Month=floor(erai$Date/100) %% 100
erai$Date2=as.POSIXct(paste(as.character(erai$Date),substr(erai$Time,1,2),sep=""),format="%Y%m%d%H",tz="GMT")

years=1990:2009
erai_ann=rep(0,20)
for(i in 1:20) erai_ann[i]=length(which(erai_E$Year==years[i]))

eracomp=matrix(0,4,3)
rownames(eracomp)=c("Cor","HR","FAR","CSI")
colnames(eracomp)=c("Nudge","No-nudge","No-nudge notopo")
for(i in 1:3) 
  {
  eracomp[1,i]=cor(count[11:30,i],erai_ann)
  eracomp[2:4,i]=CSI_days(erai[erai$Year>=1990,],fixes[[i]][fixes[[i]]$Year>=1990,])
}

###
