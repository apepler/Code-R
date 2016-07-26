mldb=read.csv("~/Documents/ECLs/mldb_events.csv",stringsAsFactors=F)
mldb$Year=floor(mldb$date/10000)
mldb$Month=floor(mldb$date/100)%%100

years=1970:2006
months=1:12
mon_ecls=matrix(0,37*12,4)

n=1
for(y in 1:37)
  for(m in 1:12)
  {
    mon_ecls[n,1]=years[y]*100+m
    I=which(mldb$Year==years[y] & mldb$Month==m)
    if(length(I)>0)
    {
      mon_ecls[n,2]=length(I)
      mon_ecls[n,3]=sum(mldb$sigimpactflag[I])
      mon_ecls[n,4]=sum(mldb$explosive[I])
    }
    n=n+1
  }

apply(mon_ecls,2,max)

mon_ecls2=mon_ecls
for(i in 2:length(mon_ecls[,1])) mon_ecls2[i,2:4]=mon_ecls[i,2:4]+mon_ecls[i-1,2:4]

I=which(mon_ecls2[,3]>=3)
J=unique(c(I,I-1))
sum(mon_ecls[J,3])
