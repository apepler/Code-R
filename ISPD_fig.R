data=scan(file="/srv/ccrc/data34/z3478332/20CR/ispd_history_v3.00.txt",what="character",sep="$")
year1=as.integer(substr(data[2:length(data)],22,25))
year2=as.integer(substr(data[2:length(data)],35,38))
lat=as.numeric(substr(data[2:length(data)],48,53))
lon=as.numeric(substr(data[2:length(data)],55,60))

years=seq(1833,2010)

statcount=matrix(0,length(years),5)
colnames(statcount)=c("Year","Aust","Aust2","SH","World")
statcount[,1]=years

for(i in 1:length(years))
{
  I=which(year1<=years[i] & year2>=years[i])
  statcount[i,5]=length(I)
  J=which(lat[I]<0)
  statcount[i,4]=length(J)
  J=which(lat[I]<0 & lat[I]>(-60) & lon[I]>100 & lon[I]<180)
  statcount[i,3]=length(J)
  J=which(lat[I]<(-10) & lat[I]>(-45) & lon[I]>110 & lon[I]<155)
  statcount[i,2]=length(J)
}

plot(years,statcount[,5],log="y",type="l",
     xlab="Year",ylab="Number of stations",lwd=3,axes=F)
axis(1,at=seq(1820,2010,10))
axis(2,at=c(0,1,10,100,1000,10000))
for(i in 4:3) lines(years,statcount[,i],lwd=3,col=i)
legend("topleft",c("Global","S. Hem.","Aust. Reg."),col=c(1,4,3),lwd=3)