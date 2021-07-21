import numpy as np
import os,sys
import time
import pysao, pyds9
from astropy.io import fits
import matplotlib.pyplot as plt
import pyfits
from centroidFindEdited import *
#global Source
Source="PKS2155-304"
CHANMIN=31
CHANMAX=801

"-------------------------------------------------"
" Sunil Chandra, TIFR, 14, February, 2016         "
" Before Using this tool please make sure that    "
" centroidFindEdited.py file is in same directory "
"-------------------------------------------------" 
#------Defining plotting Functions ------
def plottingNew(x,y,ey,lwd=1.3,size=[10.,8.75],submarzins=[0.1,0.9,0.25],xlab="TIME [s]",ylab='y(t)',xscale='',ylim=[0.0,7.0]):
    fig=plt.figure(figsize=np.array([size[0],size[1]],float))
    fig.subplots_adjust(left=float(submarzins[0]),right=float(submarzins[1]),hspace=float(submarzins[2]))
    #label_size=2.5;plt.rcParams['axes.titlesize']=5
    print "---------------\n", "Defining the window Size...\n", "The Size: "+str((size))+ "...\n", "If you wish to make changes do it in 'plottingNew' function and before the call of script...\n"
    print "Plotting The Results------"
    ax=fig.add_subplot(111)
    fig.suptitle((Source),fontsize=20)
    ax.errorbar(x,y,ey,fmt='o',lw=lwd,ecolor='red')
    ax.set_label_size=2.5
    ax.set_xlabel(xlab,fontsize=18)
    ax.set_ylabel(ylab,fontsize=18)
    ax.set_xlim(x.min()-(x.max()-x.min())*0.1,x.max()+(x.max()-x.min())*0.1)
    ax.set_ylim(ylim)
    fig.text(0.5,0.8,str(np.round(np.mean(y),2))+" +/- "+str(np.round(np.std(y),3))+" +/- "+str(np.round(np.sqrt(np.sum(ey**2))/len(ey),3)),horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,fontsize=16)
    print "The Device has made your plot...save it as per your choice..."
    plt.tick_params(axis='both',which='major',labelsize=12)
    plt.tick_params(axis='both',which='minor',labelsize=8)
    ax.axhline(y=np.mean(y),c="red",linewidth=0.9,zorder=1,linestyle='--')
#-------End of the plotting function definition--------


def xco_create(xco_path,outStem,regdir,CHANMIN,CHANMAX):
    for extr in ['src','bkg']:
        productstem=outStem+'_out_3p0-8p0_'+extr
        xco_outstr="{}_{}_out_3p0-8p0_scrpt.xco".format(outStem,extr)
        regfile="{}/{}_out{}.reg".format(regdir,outStem,extr)
        evedir="./"
        evefile=outStem+"_cl.evt"
        print outStem,xco_outstr,regfile,evefile
        fl_xco=open(os.path.join(xco_path,xco_outstr),'a')
        fl_xco.write("xsel\n")
        fl_xco.write("read events\n")
        fl_xco.write(("%s\n")%(evedir))
        fl_xco.write(("%s\n")%(evefile))
        fl_xco.write("yes\n")
        fl_xco.write("set xyname RAWX RAWY\n")
        fl_xco.write(("filter pha_cutoff %d %d \n")%(CHANMIN,CHANMAX))
        fl_xco.write(("filter region %s \n")%(regfile))
        fl_xco.write("extract all\n")
        fl_xco.write(("save all %s\n")%(productstem))
        fl_xco.write("quit\n")
        fl_xco.write("no")
        fl_xco.close()


def make_reg(reg_shape,reg_path,reg_outstr):
    global rawx, rawy, raw_r, raw_a, raw_b, area
    fl=open(os.path.join(reg_path,reg_outstr),'w')
    fl.write("# Region file format: DS9 version 4.1\n")
    fl.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    fl.write("image,\n")
    
    if reg_shape=='circle':
        rawx=float(raw_input('Enter raw X for region file'))
        rawy=float(raw_input('Enter raw Y for region file'))
        raw_r=float(raw_input('Enter radius of circle for region file'))
        area=np.pi*raw_r**2
        fl.write(("circle( %.5f, %.5f, %.5f )")%(rawx,rawy,raw_r))
        fl.close()
        return rawx, rawy, raw_r, area
    elif reg_shape=='ellipse':
        rawx=float(raw_input('Enter raw X for region file'))
        rawy=float(raw_input('Enter raw Y for region file'))
        raw_a=float(raw_input('Enter semi-major axis for region file'))
        raw_b=float(raw_input('Enter semi-minor axis for region file'))
        raw_angle=float(raw_input('Enter angle for major axis for region file'))
        area=np.pi*raw_a*raw_b
        #fl.write("ellipse("+rawx+","+rawy+","+raw_a+","+raw_b+","+raw_angle+")")
        fl.write(("ellipse( %.5f, %.5f, %.5f, %.5f, %.5f ) ")%(rawx,rawy,raw_a,raw_b,raw_angle))
        fl.close()
        return rawx, rawy, raw_a, raw_b, area

def getCoordFromFile(filename):
    fl=open(filename,'r')
    for str1 in fl:
        str2=str1.rstrip()
        str3=str2.split('(')
        if str3[0] == 'circle':
            str4 = (str3[1].split(')'))[0].split(',')
            srcarea=np.pi*float(str4[2])**2
            rawx,rawy=float(str4[0]),float(str4[1])
            return rawx,rawy,srcarea
        elif str3[0] == 'ellipse':
            str4 = (str3[1].split(')'))[0].split(',')
            srcarea=np.pi*float(str4[2])*float(str4[3])
            rawx,rawy=float(str4[0]),float(str4[1])
            #print rawx, rawy, srcarea
            return rawx,rawy,srcarea


#----------- Main Program Starts from here....
wdir=os.getcwd()
# Calling pyds9 for plotting device 
ds9=pyds9.DS9()
#--- Opening text file to write output...You may choose name as per your choice
fl_out=open(os.path.join(wdir,"Final_OutParam0p3-8p0.txt"),'a')
#----Input list file Name (path & input data files)....
fname="fileList.txt"
#-----Iterating over list of files...
count=0
fileIn=open(fname,'r')
for f in fileIn:
    dirstr=f.rstrip()
    #---Breaking paths in segments...
    spltddirstr=dirstr.split('/')
    print spltddirstr[0],spltddirstr[1],spltddirstr[len(spltddirstr)-1],spltddirstr[len(spltddirstr)-2],spltddirstr[len(spltddirstr)-3]
    #----Determining the data directory 
    dtdir=wdir+'/'+dirstr.split('AS1')[0]
    #----Generating the string for output data directory 
    outdirName=wdir+"/"+spltddirstr[0]+'/out'+spltddirstr[1]
    #----Verfying if output dir exists otherwise create ...
    if not os.path.exists(outdirName):
        print "Dir {} does not exists...\n So Creating Dir for You".format(outdirName)
        os.makedirs(outdirName)
    
    print "Output Data Directory:  {}".format(outdirName)
    print "Working Data Directory:  {}".format(dtdir)
    #-----Changing working dir to data directory 
    os.chdir(dtdir)
    #---Image file (*.img)
    fls = spltddirstr[-1]
    #---Reading the fits image file using astropy.io 
    hdulist=fits.open(fls)

#-------Reading Headers of img file......
    obsid=hdulist[0].header['OBS_ID']
    orbitid=(dtdir.split('level2_')[1]).split('/')[0]
    obsid=orbitid
    dateobs=hdulist[0].header['DATE-OBS']
    timeobs=hdulist[0].header['TIME-OBS']
    exptime=hdulist[0].header['EXPOSURE']
    mjd=hdulist[0].header['MJDREFI'] + ((hdulist[0].header['TSTART'] + hdulist[0].header['TSTOP'])/(2*86400.0))
    ra_pnt=hdulist[0].header['RA_PNT']
    dec_pnt=hdulist[0].header['DEC_PNT']
    hdulist.close()
#----------------------------Playing with Region files--------
    #----Sending image file to the ds9 for display
    ds9.set("file "+os.path.join(dtdir,fls))
    #----Input for region flag (1:circle & 2:ellipse)
    #
    #----Output Stem Structure 
    outStem=spltddirstr[-1].split('.img')[0]
    #----region file full ath
    reg_path=outdirName
    #----src region file Name  
    reg_outsrcstr="{}_outsrc.reg".format(outStem)
    print reg_path
    regfilename=(reg_path+'/'+reg_outsrcstr)
    if os.path.exists(regfilename):
    #----Running the tool to generate the region file and determing the centroid
        rawx,rawy,srcarea=getCoordFromFile(regfilename)
        centroid_xPos,centroid_yPos,centroid_unc=centroidFunc(fls,rawx,rawy,35.0)
        fl_out.write(("%s\t %s\t %s\t %.5f\t %.5f\t %.5f\t %.5f\t %.3f\t %.3f\t %.3f\t %.5f\t %.5f\t %.5f\t")%(obsid,dateobs,timeobs,mjd,exptime,ra_pnt,dec_pnt,rawx,rawy,srcarea,centroid_xPos,centroid_yPos,centroid_unc))
    else:
        print "Making the region files for you:---- \n Please provide some informations:"
        reg_flg=int(raw_input('Region Flag, 1=circle, 2=ellipse:'))
        if reg_flg==1:
            reg_shape="circle"
            rawx,rawy,raw_r,srcarea=make_reg(reg_shape,reg_path,reg_outsrcstr)
            centroid_xPos,centroid_yPos,centroid_unc=centroidFunc(fls,rawx,rawy,35.0)
            fl_out.write(("%s\t %s\t %s\t %.5f\t %.5f\t %.5f\t %.5f\t %.3f\t %.3f\t %.3f\t %.5f\t %.5f\t %.5f\t")%(obsid,dateobs,timeobs,mjd,exptime,ra_pnt,dec_pnt,rawx,rawy,srcarea,centroid_xPos,centroid_yPos,centroid_unc))

        elif reg_flg==2:
            reg_shape="ellipse"
            rawx,rawy,raw_a,raw_b,srcarea=make_reg(reg_shape,reg_path,reg_outsrcstr)
            centroid_xPos,centroid_yPos,centroid_unc=centroidFunc(fls,rawx,rawy,35.0)
            fl_out.write(("%s\t %s\t %s\t %.5f\t %.5f\t %.5f\t %.5f\t %.3f\t %.3f\t %.3f\t %.5f\t %.5f\t %.5f\t")%(obsid,dateobs,timeobs,mjd,exptime,ra_pnt,dec_pnt,rawx,rawy,srcarea,centroid_xPos,centroid_yPos,centroid_unc))
        else:
            print "Incompatible/Incorrect shape: \n use 1 (circle) or 2(ellipse) for correct one"
            os.exit()


    #print obsid,dateobs,timeobs,mjd,exptime,ra_pnt,dec_pnt,rawx,rawy,srcarea,centroid_xPos,centroid_yPos,centroid_unc

    print "Making BKG region files"
    reg_outbkgstr="{}_outbkg.reg".format(outStem)
    regfilenamebkg=reg_path+'/'+reg_outbkgstr
    if os.path.exists(regfilenamebkg):
        #----Running the tool to generate the region file and determing the centroid
        rawx,rawy,bkgarea=getCoordFromFile(regfilenamebkg)
        #centroid_xPos,centroid_yPos,centroid_unc=centroidFunc(fls,rawx,rawy,35.0)
        fl_out.write(("%.3f\t %.3f\t %.3f\t")%(rawx,rawy,bkgarea))

    else:
        print "Making Background region file for you....\n The shape is Circle"
        reg_shape='circle'
        rawx,rawy,raw_r,bkgarea=make_reg(reg_shape,reg_path,reg_outbkgstr)
        fl_out.write(("%.3f\t %.3f\t %.3f\t")%(rawx,rawy,bkgarea))

    print srcarea, bkgarea
#-------|------- Creating .XCO file ------------
    xco_path="./"
    xco_create(xco_path,outStem,reg_path,CHANMIN,CHANMAX)
    for extr in ['src', 'bkg']:
        os.system(('xselect @%s_%s_out_3p0-8p0_scrpt.xco')%(outStem,extr))

#------------- Working on lightcurves --------------
    f=fits.open(outStem+"_out_3p0-8p0_src.lc")
    #productstem=outStem+'_out_3p0-8p0_'+extr
    srclc=outStem+'_out_3p0-8p0_'+'src'
    bkglc=outStem+'_out_3p0-8p0_'+'bkg'
    outlc=outStem+'_out_3p0-8p0_'+'bkgCorrectedsrc'
    areafac=srcarea/bkgarea; newbin=30.0;print areafac, newbin
    nbint_auto=np.ceil((f[0].header['TSTOP']-f[0].header['TSTART'])/newbin)+2
    print areafac
    os.system(("lcmath %s.lc %s.lc %s.lc 1.0 %0.3f docor=yes no")%(srclc,bkglc,outlc,areafac))
    f.close()
#--------------Generating Rebbined Lightcurve -----------

    lcrvout=outlc.split('.')[0]
    giffileout=lcrvout+'.gif'
    os.system(("lcurve nser=1 plotfile=/Users/Sunil/Downloads/lcofile.pco dtnb=%f nbint=%.3f window='-' plot='yes' plotdev=%s/gif cfile1=%s.lc outfile=%s.flc  plotdnum=1")%(newbin,nbint_auto,giffileout,lcrvout,lcrvout))

    hdulist1=fits.open(lcrvout+'.flc')
    time=hdulist1[1].data['TIME']
    rate=hdulist1[1].data['RATE1']
    error=hdulist1[1].data['ERROR1']
    frac=hdulist1[1].data['FRACEXP']
    mjdlc=float(fits.open(lcrvout+'.lc')[0].header['MJD-OBS'])+(time/86400.0)
    hdulist1.close()
    
    lcouttxt=np.column_stack((mjdlc,rate,error))
    lcouttxtfile=np.savetxt(("%s_data.txt")%(outlc),lcouttxt,delimiter=" ")

    x=lcouttxt[:,0][lcouttxt[:,1]>0.0];y=lcouttxt[:,1][lcouttxt[:,1]>0.0];ey=lcouttxt[:,2][lcouttxt[:,1]>0.0]
    cnt_median=np.median(y)
    median_dev=np.sqrt( np.sum((y-cnt_median)**2)/(len(y)-1) )
    fl_out.write(("%0.7f %0.3f %0.3f %0.3f %0.3f %0.3f\n")%(x.mean(),np.mean(y),y.std(),np.sqrt(np.sum(ey**2))/len(ey),cnt_median,median_dev))
#fl_out.close()
#-------------------Plotting ------------
    plottingNew(x,y,ey,lwd=1.3,size=[10.,8.75],submarzins=[0.1,0.9,0.25],xlab="TIME [s]",ylab='y(t)',xscale='',ylim=[0.0,5.0])

    if not os.path.exists(os.path.join(wdir,"Results")):
        os.makedirs(os.path.join(wdir,"Results"))
        if not os.path.exists(os.path.join(wdir,"Results","VigResults")):
            os.makedirs(os.path.join(wdir,"Results","VigResults"))


    pltoutfileStem=spltddirstr[0]+'_'+spltddirstr[1]+"lcplot"
    pltoutStem=os.path.join(wdir,"Results","VigResults",pltoutfileStem)
    plt.savefig(pltoutStem+'.png')
    plt.savefig(pltoutStem+'.eps')
    plt.savefig(pltoutStem+'.pdf')
    #plt.show()
    plt.close()


# inp = raw_input('Enter to continue')
fl_out.close()
os.chdir(wdir)
#   break

