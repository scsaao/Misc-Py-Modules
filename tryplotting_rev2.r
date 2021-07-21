require(Hmisc)
library(plotrix)

d=read.table('~/ASTROSAT_Ver/arf_xselect_anaEdited.log')
 
dt=d
dt[,4][which(dt[,9]==37)]=dt[,4][which(dt[,9]==37)]*(60/37)^2
dt[,4][which(dt[,9]==40)]=dt[,4][which(dt[,9]==40)]*(60/40)^2

plot(dt[,7],dt[,8],pch=1,cex=1.0,col=257,xlim=c(0,600),ylim=c(0,600),xlab="X [pix]",ylab="Y [pix]",main="Position of SXT pointing for vrious offset of BL Lac PKS 2155-304",cex.lab=1.5,cex.axis=1.5,axes=TRUE)
minor.tick(nx=10,ny=10)
draw.circle(300,300,300,lwd=2.5,col=256)
points(dt[,7],dt[,8],pch=1,cex=1.5,col=257)
dev.copy2pdf(file="PointingPositiononDetectorPlane.pdf")


#first Quadrant-----------
dt1=subset(dt,dt[,7]>350 & dt[,8]>250)

plot(dt1[,7],dt1[,8],pch=8,cex=1.5,col=257,xlim=c(0,600),ylim=c(0,600),xlab="X [pix]",ylab="Y [pix]",main="Position of SXT pointing for vrious offset of BL Lac PKS 2155-304",cex.lab=1.5,cex.axis=1.5)
minor.tick(nx=10,ny=10)
dev.copy2pdf(file="PointingPositiononDetectorPlaneFirstQuad.pdf")

#--Ang Dist from centre of dectector (300,300)
angdist1=sqrt((300-dt1[,7])^2+(300-dt1[,8])^2)*(4.2/60) 
bcksubcnt1=dt1[,4]-dt1[,6]*(60/100)^2

framecnt=dt1[,2]
framecnt_er=sqrt(dt1[,1])/dt1[,10]

#xlim=c(min(angdist1)-(max(angdist1)-min(angdist1))*0.1,max(angdist1)+(max(angdist1)-min(angdist1))*0.1)
ylim=c(1.5,7)
xlim=c(0,25)
plot(angdist1,framecnt,pch=19,cex=1.5,col=257,xlim=xlim,ylim=ylim,xlab=expression(paste(theta," [ arcmin ]")),ylab="cnt/s",main="Count Rates for First Quadrant",cex.lab=1.5,cex.axis=1.5)
segments(angdist1,framecnt-framecnt_er,angdist1,framecnt+framecnt_er,col=258)

minor.tick(nx=10,ny=10)
dev.copy2pdf(file="CountOffset_FirstQuadrant_fullFrame.pdf")


#----Second Quadrant----------

dt2=subset(dt,dt[,7]<350 & dt[,8]>250)

plot(dt2[,7],dt2[,8],pch=8,cex=1.5,col=257,xlim=c(0,600),ylim=c(0,600),xlab="X [pix]",ylab="Y [pix]",main="Position of SXT pointing for vrious offset of BL Lac PKS 2155-304",cex.lab=1.5,cex.axis=1.5)
minor.tick(nx=10,ny=10)
dev.copy2pdf(file="PointingPositiononDetectorPlaneSecondQuad.pdf")

#--Ang Dist from centre of dectector (300,300)
angdist2=sqrt((300-dt2[,7])^2+(300-dt2[,8])^2)*(4.2/60) 
bcksubcnt2=dt2[,4]-dt2[,6]*(60/100)^2
x=angdist2
y=bcksubcnt2

framecnt=dt2[,2]
framecnt_er=sqrt(dt2[,1])/dt2[,10]
y=framecnt
ey=framecnt_er

#xlim=c(min(x)-(max(x)-min(x))*0.1,max(x)+(max(x)-min(x))*0.1)
#ylim=c(0,7)
#xlim=c(0,25)

plot(x,y,pch=19,cex=1.5,col=257,xlim=xlim,ylim=ylim,xlab=expression(paste(theta," [ arcmin ]")),ylab="cnt/s",main="Count Rates for Second Quadrant",cex.lab=1.5,cex.axis=1.5)
segments(x,y-ey,x,y+ey,col=258)
minor.tick(nx=10,ny=10)
dev.copy2pdf(file="CountOffset_SecondQuadrant_fullFrame.pdf")


#----Third Quadrant----------

dt3=subset(dt,dt[,7]<350 & dt[,8]<250)

plot(dt3[,7],dt3[,8],pch=8,cex=1.5,col=257,xlim=c(0,600),ylim=c(0,600),xlab="X [pix]",ylab="Y [pix]",main="Position of SXT pointing for vrious offset of BL Lac PKS 2155-304",cex.lab=1.5,cex.axis=1.5)
minor.tick(nx=10,ny=10)
dev.copy2pdf(file="PointingPositiononDetectorPlaneThirdQuad.pdf")

#--Ang Dist from centre of dectector (300,300)
angdist3=sqrt((300-dt3[,7])^2+(300-dt3[,8])^2)*(4.2/60) 
bcksubcnt3=dt3[,4]-dt3[,6]*(60/100)^2
x=angdist3
y=bcksubcnt3

framecnt=dt3[,2]
framecnt_er=sqrt(dt3[,1])/dt3[,10]
y=framecnt
ey=framecnt_er

#xlim=c(min(x)-(max(x)-min(x))*0.1,max(x)+(max(x)-min(x))*0.1)
#ylim=c(0,2)
#xlim=c(0,25)

plot(x,y,pch=19,cex=1.5,col=257,xlim=xlim,ylim=ylim,xlab=expression(paste(theta," [ arcmin ]")),ylab="cnt/s",main="Count Rates for Third Quadrant",cex.lab=1.5,cex.axis=1.5)
segments(x,y-ey,x,y+ey,col=258)
minor.tick(nx=10,ny=10)
dev.copy2pdf(file="CountOffset_ThirdQuadrant_fullFrame.pdf")


#----Fourth Quadrant----------

dt4=subset(dt,dt[,7]>350 & dt[,8]<250)

plot(dt4[,7],dt4[,8],pch=8,cex=1.5,col=257,xlim=c(0,600),ylim=c(0,600),xlab="X [pix]",ylab="Y [pix]",main="Position of SXT pointing for vrious offset of BL Lac PKS 2155-304",cex.lab=1.5,cex.axis=1.5)
minor.tick(nx=10,ny=10)
dev.copy2pdf(file="PointingPositiononDetectorPlaneFourthQuad.pdf")

#--Ang Dist from centre of dectector (300,300)
angdist4=sqrt((300-dt4[,7])^2+(300-dt4[,8])^2)*(4.2/60) 
bcksubcnt4=dt4[,4]-dt4[,6]*(60/100)^2
x=angdist4
y=bcksubcnt4

framecnt=dt4[,2]
framecnt_er=sqrt(dt4[,1])/dt4[,10]
y=framecnt
ey=framecnt_er

#xlim=c(min(x)-(max(x)-min(x))*0.1,max(x)+(max(x)-min(x))*0.1)
#xlim=c(0,25)
#ylim=c(0,2)

plot(x,y,pch=19,cex=1.5,col=257,xlim=xlim,ylim=ylim,xlab=expression(paste(theta," [ arcmin ]")),ylab="cnt/s",main="Count Rates for Fourth Quadrant",cex.lab=1.5,cex.axis=1.5)
segments(x,y-ey,x,y+ey,col=258)
minor.tick(nx=10,ny=10)
dev.copy2pdf(file="CountOffset_FouthQuadrant_fullFrame.pdf")



