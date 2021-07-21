require(Hmisc)

plotfunc=function(x,ex,y,ey,...){
    plot(x,y,...)
    xcap=0.005*(max(x)-min(x))
    ycap=0.005*(max(y)-min(y))
    segments(x,y-ey,x,y+ey)
    segments(x-ex,y,x+ex,y)
    segments(x+xcap,y-ey,x-xcap,y-ey)
    segments(x+xcap,y+ey,x-xcap,y+ey)
    segments(x-ex,y+ycap,x-ex,y-ycap)
    segments(x+ex,y+ycap,x+ex,y-ycap)
}

pointsfunc=function(x,ex,y,ey,...){
    points(x,y,...)
    xcap=0.005*(max(x)-min(x))
    ycap=0.005*(max(y)-min(y))
    segments(x,y-ey,x,y+ey)
    segments(x-ex,y,x+ex,y)
    segments(x+xcap,y-ey,x-xcap,y-ey)
    segments(x+xcap,y+ey,x-xcap,y+ey)
    segments(x-ex,y+ycap,x-ex,y-ycap)
    segments(x+ex,y+ycap,x+ex,y-ycap)
}


binDataAve=function(vec,arm1_3_1st){
    for(i in 1:NROW(vec)){
        vecout=cbind(mean(arm1_3_1st[vec[i,1]:vec[i,2],11]),median(arm1_3_1st[vec[i,1]:vec[i,2],9]),sqrt(sum(arm1_3_1st[vec[i,1]:vec[i,2],7]^2))/NROW(arm1_3_1st[vec[i,1]:vec[i,2],7]),sd(arm1_3_1st[vec[i,1]:vec[i,2],9]),mean(arm1_3_1st[vec[i,1]:vec[i,2],2]),mean(arm1_3_1st[vec[i,1]:vec[i,2],3]),mean(arm1_3_1st[vec[i,1]:vec[i,2],5]))
        if(i==1){
            vecout1=vecout
            
        } else {
            vecout2=rbind(vecout1,vecout)
            vecout1=vecout2
        }
        
    }
    return(vecout2)
}


Source="PKS_2155-304"
range="0.3-8.0 keV"
Erange="0p3-8p0"

d=read.table('Final_OutParam0p3-8p0.txt')
dtfltrd=cbind(as.character(d[,1]),d[,11:13],d[,17:22])
angsep=sqrt((dtfltrd[,2]-300)^2+(dtfltrd[,3]-300)^2)
dtout=cbind(dtfltrd,angsep)[order(angsep),]
par(mar=c(4,5,3,1))
plot(dtout[,2],dtout[,3],pch=4,cex=1.3,cex.lab=1.5,cex.axis=1.5,xlab="RA [pix]",ylab="Dec [pix]",main="PKS 2155-3014 Pointing Details\n RA:21 58 52.065; Dec:-30 13 32.12",xlim=c(0,600),ylim=c(0,600))
minor.tick(nx=10,ny=10)

arm2_4=subset(dtout,dtout[,5]<57325.25)
arm1_3_red=subset(dtout,dtout[,2]<300 & dtout[,3]>300 & dtout[,5]>57325.25 )
arm2_4_rev=rbind(arm2_4,arm1_3_red)
arm2_4=arm2_4_rev

#arm1_3=subset(dtout,dtout[,5]>57325.25)
arm1_3=subset(dtout,dtout[,5]>57325.25 & as.character(dtout[,1])!="V56.0" & as.character(dtout[,1])!="V59.0_2nd")

arm2_4_1st=subset(arm2_4,arm2_4[,2]<300)
arm2_4_2nd=subset(arm2_4,arm2_4[,2]>300)

plotfunc(-arm2_4_1st[,11],arm2_4_1st[,4],arm2_4_1st[,9],arm2_4_1st[,10],pch=4,cex=1.3,cex.lab=1.5,cex.axis=1.5,xlab="Ang. Sep. [pix]",ylab="rate [c/s]",ylim=c(0,5),xlim=c(-300,300))
pointsfunc(arm2_4_2nd[,11],arm2_4_2nd[,4],arm2_4_2nd[,9],arm2_4_2nd[,10],pch=4,cex=1.3)

arm1_3_1st=subset(arm1_3,arm1_3[,2]<300)
arm1_3_2nd=subset(arm1_3,arm1_3[,2]>300)


plotfunc(-arm1_3_1st[,11],arm1_3_1st[,4],arm1_3_1st[,9],arm1_3_1st[,10],pch=4,cex=1.3,col=258,xlim=c(-300,300),ylim=c(0,5))
pointsfunc(arm1_3_2nd[,11],arm1_3_2nd[,4],arm1_3_2nd[,9],arm1_3_2nd[,10],pch=4,cex=1.3,ylim=c(5,0),col=258)

vec=cbind(c(1,6,12,14,23),c(5,11,13,22,25))
arm1_3_1stave=binDataAve(vec,arm1_3_1st)

vec=cbind(c(1,3,7,11,15),c(2,6,10,14,22))
arm1_3_2ndave=binDataAve(vec,arm1_3_2nd)


plotfunc(-arm1_3_1stave[,1],1.9,arm1_3_1stave[,2],arm1_3_1stave[,3],pch=4,cex=1.3,col=258,xlim=c(-300,300),ylim=c(0,5))
pointsfunc(arm1_3_2ndave[,1],1.9,arm1_3_2ndave[,2],arm1_3_2ndave[,3],pch=4,cex=1.3,ylim=c(5,0),col=258)

vec=cbind(c(1,5),c(4,8))
arm2_4_1stave=rbind(binDataAve(vec,arm2_4_1st),cbind(arm2_4_1st[9,11],arm2_4_1st[9,9],arm2_4_1st[9,7],sd(arm2_4_1st[9,9]),arm2_4_1st[9,2],arm2_4_1st[9,3],arm2_4_1st[9,5]) )

vec=cbind(c(1,7,12,19,24,26,29,31),c(6,11,18,23,25,28,30,38))
arm2_4_2ndave=binDataAve(vec,arm2_4_2nd)

write.table(rbind(arm2_4_1stave,arm2_4_2ndave),file=paste(Source,"_countMapFile_arm1-3_",Erange,".dat",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

write.table(rbind(arm1_3_1stave,arm1_3_2ndave),file=paste(Source,"_countMapFile_arm2-4_",Erange,".dat",sep=""),col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")


#---------------Plotting Count Map----------
plotfunc(-arm1_3_1stave[,1],1.9,arm1_3_1stave[,2],arm1_3_1stave[,3],pch=4,cex=1.3,col=258,xlim=c(-300,300),ylim=c(0,5.2),xlab=expression(paste(delta, " [ pix ] ")),ylab="Rate [c/s]",cex.lab=1.5,cex.axis=1.5,main="HBL PKS 2155-304 (z~0.116)")
pointsfunc(arm1_3_2ndave[,1],1.9,arm1_3_2ndave[,2],arm1_3_2ndave[,3],pch=4,cex=1.3,ylim=c(5,0),col=258)

pointsfunc(-arm2_4_1stave[,1],1.9,arm2_4_1stave[,2],arm2_4_1stave[,3],pch=5,cex=1.3,ylim=c(50,0),col=259)
pointsfunc(arm2_4_2ndave[,1],1.9,arm2_4_2ndave[,2],arm2_4_2ndave[,3],pch=5,cex=1.3,ylim=c(50,0),col=259)

abline(v=0.0,lty=6,lwd=2.3,col=256)
legend(x='topright',paste("E2-E1:", range),text.col=258,box.lty=4,box.lwd=1.5,box.col=256)
legend(x='topleft',c("Diag 1-3", "Diag 2-4"),pch=4:5,col=258:259,text.col=258:259,box.lty=4,box.lwd=1.5,box.col=256)

text(0,0.15,"(300,300)",col='red',cex=1.3)
minor.tick(nx=10,ny=10)

dev.copy2pdf(file=paste(Source,"_countMapFile_arm1-3_",Erange,".pdf",sep=""))
