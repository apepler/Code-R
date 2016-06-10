########## Matching lows

fixesE<-read.csv("~/output/outputUM_erai_150_topo_rad2_proj100/ECLfixes_umelb_erai_150_topo_rad2_proj100.csv")
fixesE=fixesE[fixesE$Date>=20070000 & fixesE$Date<=20090000,]
eventsE<-read.csv("~/output/outputUM_erai_150_topo_rad2_proj100/ECLevents_umelb_erai_150_topo_rad2_proj100.csv")
eventsE=eventsE[eventsE$Date1>=20070000 & eventsE$Date1<=20090000,]
library(sp)

matchlows<-function(d1,d2,res)
{
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  ##So same time period & curvature threshold; only d1 is location-restricted
  d1=d1[y1>=ymin & y1<=ymax & d1$Location==1,]
  d2=d2[y2>=ymin & y2<=ymax,]
  
  match=matrix(NaN,length(d1[,1]))  
  for(i in 1:length(d1[,1])) 
  {
    I=which(d2$Date==d1$Date[i] & d2$Time==d1$Time[i])
    if(length(I)>0) match[i]=min(spDistsN1(as.matrix(cbind(d2$Lon[I],d2$Lat[I])),as.numeric(c(d1$Lon[i],d1$Lat[i])),longlat=TRUE))
  } 
  
  res=cbind(res,matrix(0,length(res),2))
  for(i in 1:length(res[,1])){
    res[i,2]=length(which(match<=res[i,1]))/length(match)
    res[i,3]=length(which(match<=res[i,1]))/sum(1-is.na(match))
  } 
  
  return(res)
}

matchevents<-function(d1,d2,res)
{
  ##Where different datasets in terms of curvature etc, but need to create location
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  d1=d1[y1>=ymin & y1<=ymax & d1$Location==1,]
  d2=d2[y2>=ymin & y2<=ymax,] 
  
  x<-rle(d1$ID)
  ev1<-data.frame(ID=x$values,Length=x$lengths,Location=rep(0,length(x$values)),MinMatch=rep(NaN,length(x$values)))
  for(i in 1:length(ev1[,1])) ev1[i,3]=max(d1$Location[d1$ID==ev1[i,1]])
  ev1=ev1[ev1$Location==1,]
  
  match1=matrix(NaN,length(d1[,1]))
  for(i in 1:length(d1[,1])) 
  {
    I=which(d2$Date==d1$Date[i] & d2$Time==d1$Time[i])
    if(length(I)>0) match1[i]=min(spDistsN1(as.matrix(d2[I,7:8]),as.numeric(d1[i,7:8]),longlat=TRUE))
  } 
  for(i in 1:length(ev1[,1])) ev1[i,4]=min(match1[d1$ID==ev1[i,1] & d1$Location==1,],na.rm=T)
  
  res=cbind(res,matrix(0,length(res),1))
  for(i in 1:length(res[,1])) res[i,2]=length(which(ev1[,4]<=res[i,1]))/length(ev1[,4]) 
  
  return(res)
}

matchF<-matchE<-cbind(c(10,50,100,250,500,1000),matrix(0,6,6))
for(i in 1:6)
{
  a=matchlows(fixesE,fixes[[i]],c(10,50,100,250,500,1000))
  matchF[,i+1]=a[,2]
  b=matchevents(fixesE,fixes[[i]],c(10,50,100,250,500,1000))
  matchE[,i+1]=b[,2]
}

### What about showing longitudinal distance?

matchlows_distance<-function(d1,d2)
{
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  ##So same time period & curvature threshold; only d1 is location-restricted
  d1=d1[y1>=ymin & y1<=ymax & d1$Location==1,]
  d2=d2[y2>=ymin & y2<=ymax,]
  
  match=matrix(NaN,length(d1[,1]),3) 
  
  match[,1]=d1$ID+floor(d1$Date/10000)*1000
  for(i in 1:length(d1[,1])) 
  {
    I=which(d2$Date==d1$Date[i] & d2$Time==d1$Time[i])
    if(length(I)>0)
    {
      J=(spDistsN1(as.matrix(cbind(d2$Lon[I],d2$Lat[I])),as.numeric(c(d1$Lon[i],d1$Lat[i])),longlat=TRUE))
      k=which(J==min(J))
      match[i,2]=min(J)
      match[i,3]=d2$Lon[I[k]]-d1$Lon[i]
    }
  } 
  return(match)
}

matchlows_CV<-function(d1,d2)
{
  y1=floor(d1$Date/10000)
  y2=floor(d2$Date/10000)
  ymin=max(min(y1),min(y2))
  ymax=min(max(y1),max(y2))
  
  ##So same time period & curvature threshold; only d1 is location-restricted
  d1=d1[y1>=ymin & y1<=ymax & d1$Location==1,]
  d2=d2[y2>=ymin & y2<=ymax,]
  
  match=matrix(NaN,length(d1[,1]),4) 
  
  match[,1]=d1$ID+floor(d1$Date/10000)*1000
  for(i in 1:length(d1[,1])) 
  {
    I=which(d2$Date==d1$Date[i] & d2$Time==d1$Time[i])
    if(length(I)>0)
    {
      J=(spDistsN1(as.matrix(cbind(d2$Lon[I],d2$Lat[I])),as.numeric(c(d1$Lon[i],d1$Lat[i])),longlat=TRUE))
      k=which(J==min(J))
      match[i,2]=min(J)
      match[i,3]=d1$CV[i]
      match[i,4]=d2$CV[I[k]]
    }
  } 
  return(match)
}

testmatch=matrix(0,6,6)
colnames(testmatch)=c("MatchedLows","MatchedEvents","LowAveCVdiff","EventAveCVdiff","% lows <-0.1","% lows >0.1")

res=250
for(i in 1:6)
{
  c=matchlows_CV(fixes[[i]],fixesNT[[i]])
  I=which(c[,2]<=res)
  testmatch[i,1]=length(I)/length(c[,1])
  testmatch[i,3]=mean(c[I,4]-c[I,3])
  J=which(c[I,4]-c[I,3]<=-0.1)
  testmatch[i,5]=length(J)/length(I)
  J=which(c[I,4]-c[I,3]>=0.1)
  testmatch[i,6]=length(J)/length(I)
  
  x<-unique(c[,1])
  d<-cbind(x,matrix(0,length(x),2))
  for(j in 1:length(x))
  {
    I=which(c[,1]==d[j,1])
    d[j,2]=min(c[I,2])
    J=which(c[I,2]<=res)
    d[j,3]=mean(c[I[J],4]-c[I[J],3])
  }
  
  I=which(d[,2]<=res)
  testmatch[i,2]=length(I)/length(d[,1])
  testmatch[i,4]=mean(d[I,3])
}


a=fixes[[3]][fixes[[3]]$Location==1,]
a=cbind(a,matchlows_distance(fixes[[3]],fixesNT[[3]]))
I=which(a[,18]<=1000)
a=a[I,]
x=rle(a[,17])

b<-cbind(x$values,x$lengths,matrix(0,length(x$values),9))
colnames(b)<-c("ID","Length","Date1","Year","Month","Lon","Lat","MSLP","CV","Radius","LonChange")
for(i in 1:length(x$values))
{
  I=which(a[,17]==b[i,1])
  b[i,3:5]=as.numeric(a[I[1],c(4,15,16)])
  b[i,6:11]=as.numeric(apply(a[I,c(7,8,9,10,12,19)],2,mean))
}

thresh=c(-Inf,-1,0,1,Inf)
stats<-matrix(0,4,11)
colnames(stats)=colnames(b)
for(i in 1:4)
{
  I=which(b[,11]>=thresh[i] & b[,11]<thresh[i+1])
  stats[i,]=apply(b[I,],2,mean)
}

###########
### What about for ECLs right along the coast?
#########

fixesE$Location2<-0
I<-which(fixesE[,7]>=149 & fixesE[,7]<=154 & fixesE[,8]<(-37) & fixesE[,8]>=-41)
fixesE$Location2[I]<-1
I<-which(fixesE[,7]>=(149+(37+fixesE[,8])/2) & fixesE[,7]<=(154+(37+fixesE[,8])/2) & fixesE[,8]<(-31) & fixesE[,8]>=-37)
fixesE$Location2[I]<-1
I<-which(fixesE[,7]>=152 & fixesE[,7]<=157 & fixesE[,8]<=(-24) & fixesE[,8]>=-31)
fixesE$Location2[I]<-1

for(n in 1:6)
{
  fixes[[n]]$Location2<-0
  I<-which(fixes[[n]][,7]>=149 & fixes[[n]][,7]<=154 & fixes[[n]][,8]<(-37) & fixes[[n]][,8]>=-41)
  fixes[[n]]$Location2[I]<-1
  I<-which(fixes[[n]][,7]>=(149+(37+fixes[[n]][,8])/2) & fixes[[n]][,7]<=(154+(37+fixes[[n]][,8])/2) & fixes[[n]][,8]<(-31) & fixes[[n]][,8]>=-37)
  fixes[[n]]$Location2[I]<-1
  I<-which(fixes[[n]][,7]>=152 & fixes[[n]][,7]<=157 & fixes[[n]][,8]<=(-24) & fixes[[n]][,8]>=-31)
  fixes[[n]]$Location2[I]<-1
  
  fixesNT[[n]]$Location2<-0
  I<-which(fixesNT[[n]][,7]>=149 & fixesNT[[n]][,7]<=154 & fixesNT[[n]][,8]<(-37) & fixesNT[[n]][,8]>=-41)
  fixesNT[[n]]$Location2[I]<-1
  I<-which(fixesNT[[n]][,7]>=(149+(37+fixesNT[[n]][,8])/2) & fixesNT[[n]][,7]<=(154+(37+fixesNT[[n]][,8])/2) & fixesNT[[n]][,8]<(-31) & fixesNT[[n]][,8]>=-37)
  fixesNT[[n]]$Location2[I]<-1
  I<-which(fixesNT[[n]][,7]>=152 & fixesNT[[n]][,7]<=157 & fixesNT[[n]][,8]<=(-24) & fixesNT[[n]][,8]>=-31)
  fixesNT[[n]]$Location2[I]<-1
}

testmatch=matrix(0,6,6)
colnames(testmatch)=c("MatchedLows","MatchedEvents","LowAveLonDiff","EventAveLonDiff","% lows <-1","% lows >1")

res=1000
for(i in 1:6)
{
  c=matchlows_distance(fixes[[i]],fixesNT[[i]])
  I=which(c[,2]<=res)
  testmatch[i,1]=length(I)/length(c[,1])
  testmatch[i,3]=mean(c[I,3])
  J=which(c[I,3]<=-1)
  testmatch[i,5]=length(J)/length(I)
  J=which(c[I,3]>=1)
  testmatch[i,6]=length(J)/length(I)
  
  x<-unique(c[,1])
  d<-cbind(x,matrix(0,length(x),2))
  for(j in 1:length(x))
  {
    I=which(c[,1]==d[j,1])
    d[j,2]=min(c[I,2])
    J=which(c[I,2]<=res)
    d[j,3]=mean(c[I[J],3])
  }
  
  I=which(d[,2]<=res)
  testmatch[i,2]=length(I)/length(d[,1])
  testmatch[i,4]=mean(d[I,3])
}

fixes2<-list()
for(i in 1:6)
{
  a=fixes[[i]][fixes[[i]]$Location2==1,]
  a=cbind(a,matchlows_distance(fixes[[i]],fixesNT[[i]]))
  I=which(a[,19]<=500)
  fixes2[[i]]=a[I,]
}

monthcount<-month1<-monthave<-matrix(0,3,6)

for(j in 1:6)
{
  I=which(fixes2[[j]]$Lat<=-38)
  monthcount[1,j]=length(I)
  if(length(I)>0) 
  {
    monthave[1,j]=mean(fixes2[[j]][I,20])
    J=which(fixes2[[j]][I,20]<=-1)
    month1[1,j]=length(J)/length(I)
  }
  
  I=which(fixes2[[j]]$Lat>(-38) & fixes2[[j]]$Lat<(-34))
  monthcount[2,j]=length(I)
  if(length(I)>0) 
  {
    monthave[2,j]=mean(fixes2[[j]][I,20])
    J=which(fixes2[[j]][I,20]<=-1)
    month1[2,j]=length(J)/length(I)
  }
  
  I=which(fixes2[[j]]$Lat>=-34)
  monthcount[3,j]=length(I)
  if(length(I)>0) 
  {
    monthave[3,j]=mean(fixes2[[j]][I,20])
    J=which(fixes2[[j]][I,20]<=-1)
    month1[3,j]=length(J)/length(I)
  }
}

freq<-matrix(0,4,7)

I=which(fixesE$Location==1)
freq[1,1]=length(I)
I=which(fixesE$Location2==1)
freq[2,1]=length(I)

for(i in 1:6)
{
  I=which(fixes[[i]]$Location==1)
  freq[1,1+i]=length(I)
  I=which(fixes[[i]]$Location2==1)
  freq[2,1+i]=length(I)
  I=which(fixesNT[[i]]$Location==1)
  freq[3,1+i]=length(I)
  I=which(fixesNT[[i]]$Location2==1)
  freq[4,1+i]=length(I)
  
}



######################
######################
#######################
## What makes a low more likely to be matched?

testmatch<-testcount<-matrix(NaN,200,3)

res=500
for(i in 1:3)
{
  c=matchlows_distance(fixes[[i]],fixesNT[[i]])
  I=which(c[,2]<=res)
  
  x<-unique(c[,1])
  for(j in 1:length(x))
  {
    I=which(c[,1]==x[j] & c[,2]<=res)
    testcount[j,i]=length(I)
    if(length(I)>0) testmatch[j,i]=mean(c[I,3])
  }
}

typecount<-typemean<-matrix(0,7,3)
rownames(typecount)<-rownames(typemean)<-c(types,"Bomb","Cool","Warm")

for(j in 1:3)
{
  for(i in 1:4)
  {
    I=which(tlist[,j]==types[i] & testcount[,j]>0)
    typecount[i,j]=length(I)
    typemean[i,j]<-mean(testmatch[I,j])
  }
  I=which(bomb[,j]==1 & testcount[,j]>0)
  typecount[5,j]=length(I)
  if(length(I)>0) typemean[5,j]=mean(testmatch[I,j])
  
  month=floor(events[[j]]$Date1/100)%%100
    I=which(month>=5 & month<=10 & testcount[1:length(month),j]>0)
    typecount[6,j]=length(I)
    if(length(I)>0) typemean[6,j]=mean(testmatch[I,j])
    I=which((month<5 | month>10) & testcount[1:length(month),j]>0)
    typecount[7,j]=length(I)
    if(length(I)>0) typemean[7,j]=mean(testmatch[I,j])
}

######## What about by season

