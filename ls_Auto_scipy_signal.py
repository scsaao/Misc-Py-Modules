#!/bin/python
import os; import time 
import pyfits; import numpy as np; import scipy as sp
from astroML.time_series import lomb_scargle_BIC, lomb_scargle_bootstrap
from scipy.signal import lombscargle
from astroML.time_series import lomb_scargle
from gatspy.periodic import LombScargle
from gatspy.periodic import LombScargleFast
import matplotlib.pyplot as plt


# use seaborn for plot styles
import seaborn; seaborn.set()

#------Defining Functions ------

def plottingNew(x,y,ey,periods,power,lwd=1.3,size=[10.,6.75],submarzins=[0.1,0.9,0.25],xlab1='TIME [days]', ylab1='y(t)',xlab2='$\omega$' + str('[days]'),ylab2='Power',xscale='log'):
    
    fig=plt.figure(figsize=np.array([size[0],size[1]],float))
    fig.subplots_adjust(left=float(submarzins[0]), right=float(submarzins[1]), hspace=float(submarzins[2]))
    #label_size=5;plt.rcParams['axes.titlesize']=18
    #------------------------------------------------------------
    # Plot the results
    print "--------------------------------------\n","Defining the Window Size... \n", "The Size : "+ str((size))+ "...\n", " The margins for Subplots (l,r,h): " + str((submarzins))+ "...\n", "If you wish to change, do it during the call of 'plottingNew' function...\n", "--------------------------------------\n"
    
    print "plotiing the results..."
    # First panel: the data
    ax = fig.add_subplot(211)
    ax.errorbar(x, y, ey, fmt='.k', lw=lwd, ecolor='gray')
    ax.set_xlabel(xlab1)
    ax.set_ylabel(ylab1)
    ax.set_xlim(x.min()-(x.max()-x.min())*0.1, x.max()+(x.max()-x.min())*0.1)
    
    # Second panel: the periodogram & significance levels
    ax1 = fig.add_subplot(212, xscale=xscale)
    ax1.plot(periods, power, '-', c='black', lw=1.5, zorder=1)
    #    ax1.plot([period[0], period[-1]], [sig1, sig1], ':', c='black')
    #    ax1.plot([period[0], period[-1]], [sig5, sig5], ':', c='black')
    
    #    ax1.annotate("", (0.3, 0.65), (0.3, 0.85), ha='center',
    #             arrowprops=dict(arrowstyle='->'))
    
    ax1.set_xlim(periods[0], periods[-1])
    ax1.set_ylim(power.min()-((power.max()-power.min())*0.05), power.max()+((power.max()-power.min())*0.05))
    
    ax1.set_xlabel(xlab2)
    ax1.set_ylabel(ylab2)

#-------End of the plotting function definition--------

instrument='sxt'
paramsIn=['TIME','Roll_RA','Roll_DEC','ANG_OFFSET']
filelist="fileListatServer"

with open(filelist) as input_file:
    for i, line in enumerate(input_file):
        filestring=line.split('/')
        dir1=filestring[len(filestring)-4];dir2=filestring[len(filestring)-2]
        #--------Creaking Output Directory Tree--------
        print "Creating Output Directory Tree...."
        if not os.path.exists(str(dir1)):
            os.mkdir(dir1);print "made a folder named:"+dir1
            os.chdir(dir1)
            if not os.path.exists(dir2):
                os.mkdir(dir2)
            os.chdir('..')
        print 'Current Working Diectory is: ', os.getcwd()
        
        #----------------------------------------
        print "Reading fits file using the path from a file: " + filelist
        #Reading Fits File from path 
        fileIN=pyfits.open(line.split('\n')[0])
        #----------------Reading Header and fetching file info-----------
        Source=fileIN[0].header['OBJECT']
        Instrument=fileIN[0].header['INSTRUME']
        tstart=float(fileIN[0].header['TSTARTI']) + float(fileIN[0].header['TSTARTF'])
        tstop=float(fileIN[0].header['TSTOPI']) + float(fileIN[0].header['TSTOPF'])
        exposure=tstop-tstart
        tstart_mjd=float(fileIN[0].header['MJDREFI']) + (tstart/86400.0)
        tstop_mjd=float(fileIN[0].header['MJDREFI']) + (tstop/86400.0)
        time_obs_mjd=float(fileIN[0].header['MJDREFI']) + ((tstop+tstart)/(2.0*86400.0))
        time_obs_UTC=fileIN[0].header['DATE-OBS']+'T'+fileIN[0].header['TIME-OBS']

        #---RA_NOM/DEC_NOM and/or RA_PNT/RA_PNT ---Nominal RA/Dec direction of the telescope [in Degree]
        ra_nom=float(fileIN[0].header['RA_NOM']);dec_nom=float(fileIN[0].header['DEC_NOM'])
        ra_pnt=float(fileIN[0].header['RA_REF']);dec_pnt=float(fileIN[0].header['DEC_REF'])

        #---Differences in RA, Dec, and diagonal Differences [in arcmin]-----
        delRA=ra_nom-ra_pnt;delDec=dec_nom-dec_pnt;delDiag=(np.sqrt(delRA**2 + delDec**2))*60.
        #----------------------------------------------------------------
        
        print "Reading the user input data file" + line.split('/')[-1]
        t_tmp=np.array(fileIN[1].data['TIME'],float)
        print "Number of data points used for this analysis: " + str(len(t_tmp))
        delT=[]
        print "estimating the smallest data-gap..."
        for u in range(len(t_tmp)-1):
            delT.append(t_tmp[u+1]-t_tmp[u])
        print "smallest data-gap in our time series: ",min(delT)
        print "Generating false errors (const=1) in case not given in input array: "
        dmag=np.ones(len(t_tmp))
        dy=np.zeros(len(t_tmp))
        # Choose a period grid
        print "Creating period vector based on the minum data gaps and number of data points..."
        periods = np.linspace(2*min(delT), int((t_tmp.max()-t_tmp.min())/3.0), 400000)
        ang_freqs = 2 * np.pi / periods

        #-----Deciding the Stem of OutputFile...
        OutStem=dir1+'_'+instrument+'_TSstudyOut_'   
        #----opening files to write output
        fileInfout =open(dir1+'/'+dir2+'/'+OutStem+'TS_TP.info','a')
        dtpow = [None]*(len(paramsIn)-1)
        statvector=[None]*(len(paramsIn)-1)

        #-----------------Wrting a file for Stat and Other header Info....
        fileInfout.write(('# CREATED/MODIFIED: %s\n')%(time.asctime()))
        fileInfout.write(('# Source: %s \t RA_POINT (Deg.): %0.3f \t DEC_POINT (Deg.): %0.3f\n')%(Source,ra_pnt,dec_pnt))
        fileInfout.write(('# Source: %s \t \t RA_Nom (Deg.): %0.3f \t DEC_Nom (Deg.): %0.3f\n')%(Source,ra_nom,dec_nom,))
        fileInfout.write(('# DATE-OBS (UTC): %s \tDATE-OBS (MJD): %0.3f\n')%(time_obs_UTC,time_obs_mjd))
        fileInfout.write(('# Exposure (s): %0.1f;\tTSTART(s): %0.1f (%0.1f) \tTSTOP(s): %0.1f (%0.1f)\n')%(exposure,tstart,tstart_mjd,tstop,tstop_mjd))
        fileInfout.write("#---------------------------------------------------")


        #--------------------------Main Loop For Analysis-----------------------        
        for j in range(len(paramsIn)-1):
            # compute the (unnormalized) periodogram
            # note pre-centering of y values!
            t=np.array(fileIN[1].data[paramsIn[0]],float)
            mag=np.array(fileIN[1].data[paramsIn[j + 1]],float)
            print "The TS being used is for: "+ str(paramsIn[j + 1])
            #-------------- Estimating the statistics on input data......
            stat=[len(t),round(max(t)-min(t),5),round(max(mag),5),round(min(mag),5),round(np.mean(mag),5),round(np.median(mag),5),round(np.std(mag),5)]  
            frmt=['N','Delta_T', 'MAX', 'MIN',  'MEAN', 'MEDIAN','Standard Dev.' ]
            #----Using Scipy.signal
            print "Performing the periodogram analysis using scipy.signal.lombscargle"
            power = lombscargle(t, mag - mag.mean(), ang_freqs)
            # normalize the power
            print "normalize the power"
            power *= 2 / (len(t) * mag.std() ** 2)

            #-------------Generating Output Vectors------
            if j == 0:
                statvector[j]=np.column_stack((frmt,stat))
                dtpow[j]=np.column_stack((periods,power))
                print "making two Dimentional array for output..."
            else:
                statvector[j]=np.column_stack((statvector[j-1],stat))
                dtpow[j]=np.column_stack((dtpow[j-1],power))
                print "Adding one more column to the output data matrix... "

            #----------Sketching Output on Graphics------ 
            plottingNew(x=t,y=mag,ey=dy,periods=periods,power=power,lwd=1.3,size=[10.,6.75],submarzins=[0.1,0.9,0.25],xlab1='TIME [days]', ylab1='y(t)',xlab2='$2 \pi / \omega$ ' + str(' [days]'),ylab2=r'$P_{LS} (\omega)$',xscale='log')
    
            print "Writing plots to files"
            plt.savefig(dir1+'/'+dir2+'/'+OutStem + str(paramsIn[j + 1]) + '.png')
            plt.savefig(dir1+'/'+dir2+'/'+OutStem + str(paramsIn[j + 1]) + '.eps')
            plt.savefig(dir1+'/'+dir2+'/'+OutStem + str(paramsIn[j + 1]) + '.pdf')
 
        print line.split('/')[1],len(line.split('/')),dir1,dir2

        #-----Writing The output vector to a file.....
        print "Writing Output to Files/..../"
        print dir1+'/'+dir2+'/'+OutStem+'TS_TP.out'
        fileformat=dir1+'/'+dir2+'/'+OutStem+'TS_TP.out'
        np.savetxt(fileformat,dtpow[j],delimiter=" ")
        np.savetxt(fileInfout,statvector[j],delimiter=" ",fmt='%15s')
        #fileOut.write(dtpow[j])
        #fileOut.close()
        fileIN.close()
print "{0} line(s) printed".format(i+1)
