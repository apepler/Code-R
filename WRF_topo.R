##Loading the input data
setwd("~/WRF")
library(RNetCDF)
f1=open.nc("wrfinput_2.5")
lat1=var.get.nc(f1,"XLAT")
lon1=var.get.nc(f1,"XLONG")
mask1=var.get.nc(f1,"LANDMASK")
topo1=var.get.nc(f1,"HGT")

f1=open.nc("wrfinput_0.5")
lat2=var.get.nc(f1,"XLAT")
lon2=var.get.nc(f1,"XLONG")
mask2=var.get.nc(f1,"LANDMASK")
topo2=var.get.nc(f1,"HGT")

topo2a=topo2[seq(3,215,5),seq(3,144,5)]
topo2b=matrix(0,43,28)
for(i in 1:43)
  for(j in 1:28)
    topo2b[i,j]=mean(topo2[((i*5)-4):(i*5),((j*5)-4):(j*5)])

library(fields)
image.plot(topo1)
image.plot(topo2a) #Centre point
image.plot(topo2b) #Average across grid

f1=open.nc('WRFtopo_regrid2.5.nc')
mask3=var.get.nc(f1,"LANDMASK")
topo3=var.get.nc(f1,"HGT")
image.plot(topo3)

#Same scale?
image.plot(topo1,breaks=seq(-50,1000,50),col=tim.colors(21)) 
image.plot(topo2b,breaks=seq(-50,1000,50),col=tim.colors(21)) 
image.plot(topo3,breaks=seq(-50,1000,50),col=tim.colors(21))

f1=open.nc('WRFtopo_regrid0.5.nc')
mask4=var.get.nc(f1,"LANDMASK")
topo4=var.get.nc(f1,"HGT")
image.plot(topo2)
image.plot(topo1)
image.plot(topo4)

topo5=topo4
topo5[mask2==0]=0
topo5[is.na(topo5)==1]=topo2[is.na(topo5)==1]
image.plot(topo5)

topo6=topo2
topo6[32:129,25:102]=topo5[32:129,25:102]
topo6[32:60,95:102]=topo2[32:60,95:102]
image.plot(topo6)

f1=open.nc('WRFtopo_regrid0.5.nc')
topo5a=var.get.nc(f1,'HGT_v1')
topo6a=var.get.nc(f1,'HGT_v2')
image.plot(topo5)
image.plot(topo5a)
image.plot(topo6)
image.plot(topo6a)







