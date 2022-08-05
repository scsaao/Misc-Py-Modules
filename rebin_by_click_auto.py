#!/bin/python

"""
 This script bins the data based on the clickes by the user.
 This is nearly a genralized script and hence can be used for any purpose.
 We do not claim that the resulting data plot will not be biased towards the clicks.
 Think and Judge your steps before clicking on the plot...
 Best Wishes!!!

updated on 03 May 2017 (01:15 AM) to take the input from command line...
You need to run > python lcbinmaker_auto.py -h
for informations about the input parameters. 

"""


def profile_plot(x=None, y=None, ex=None, ey=None, Source="Source", lwd=1.3,size=[10.,8.75], submarzins=[0.1,0.9,0.25], xlab="RAWX [pix]", ylab='Rate [c/s]', xscale='', yscale='', xlim=[0,600], outfile = "OutFIle.png"):

	import numpy as np
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	import seaborn; seaborn.set()

	fig=plt.figure(figsize=np.array([size[0],size[1]],float))
	fig.subplots_adjust(left=float(submarzins[0]),right=float(submarzins[1]),hspace=float(submarzins[2]))

	delt_ylim = (y.max() - y.min())*0.45
	ax=fig.add_subplot(111)
	ax.errorbar(x,y,ey,ex,fmt='ob',lw=lwd,ecolor='red')
	ax.set_label_size = 2.5
	ax.set_xlabel(xlab,fontsize=18)
	ax.set_ylabel(ylab,fontsize=18)

	ax.set_ylim([y.min() - delt_ylim, y.max() + delt_ylim])
	ax.set_xlim(xlim)
	seaborn.set_style("ticks")
	plt.tick_params(axis='both',which='major',labelsize=12)
	plt.tick_params(axis='both',which='minor',labelsize=8)

	plt.savefig('{}'.format(outfile))

	status = 0	

	print status

def rebin_data_using_plot(x=None, y=None, ex=None, ey=None, xlab=None, ylab=None, xlim=None):

	import numpy as np
	import matplotlib.pyplot as plt
	import seaborn; seaborn.set()

	seaborn.set_style("ticks")

	delt_ylim = (y.max() - y.min())*0.45	

	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.errorbar(x, y, ey, ex, "*r")
	#ax.errorbar(x1,y1,ex1,ey1,"ob")

	ax.set_label_size = 2.5
	ax.set_xlabel(xlab,fontsize=18)
	ax.set_ylabel(ylab,fontsize=18)

	ax.set_ylim([y.min() - delt_ylim, y.max() + delt_ylim])
	ax.set_xlim(xlim)
	
	
	coords = []
	def onclick(event):
		global ix, iy
		print 'button=%d, xdata=%f, ydata=%f'%(event.button, event.xdata, event.ydata)
		button, ix, iy = event.button, event.xdata, event.ydata
		#global coords
		
		if button == 3:
			fig.canvas.mpl_disconnect(cid)
			plt.close(event.canvas.figure)
		else:
			coords.append((ix, iy))
			print "Keep Going....Loop No. "
			
		return coords
	
	cid = fig.canvas.mpl_connect('button_press_event', onclick)
	plt.show()

	X_start = []
	X_stop = []
        #print len(coords)
	if (len(coords) % 2 == 0) & (len(coords) >0) :

		print "Length of Coords Selected are even so may be a correct choice"

		for i in range(len(coords))[::2] :
			X_start.append([coords[i][0]])
	
		for j in range(len(coords))[1::2] :
			X_stop.append([coords[j][0]])
	
		X_range = np.column_stack((X_start, X_stop))
	
	return X_range, len(X_range)



def do_binning_procedure2(lc_table, X_range, posstr = "TIME", ratestr = "RATE1", errorstr = "ERROR1", fracexpstr = "FRACEXP", size=[10.,8.75], submarzins=[0.1,0.9,0.25], xlim= None, binsize=120) :

    import numpy as np
    from astropy.table import Table, Column

    Total_cnts =  lc_table[ratestr] * binsize * lc_table[fracexpstr]
    #------------------------------------------
    tbl = Table()
    tbl.add_column(Column(Total_cnts, name = 'COUNTS'), index=4)
    tbl.add_column(Column(lc_table[posstr], name = 'TIME'), index=0)
    tbl.add_column(Column(lc_table[ratestr], name = 'RATE'), index=1)
    tbl.add_column(Column(lc_table[errorstr], name = 'ERROR'), index=2)
    tbl.add_column(Column(lc_table[fracexpstr], name = 'FRACEXP'), index=3)
    #------------------------------------------

    dt_mean=([])

    rawx_indx = np.array(tbl["TIME"],'float')

    for count in range(len(X_range)) :
         	
        subtbl = tbl[ (rawx_indx >= X_range[count,0]) & (rawx_indx <= X_range[count,1]) ]

        sdcnts = np.std(subtbl['COUNTS'])
        newtime, newcnt, newerror = np.nanmean(subtbl['TIME']), np.nansum(subtbl['COUNTS']), np.sqrt(np.nansum(subtbl['COUNTS']))
        delT = sum(binsize * subtbl[fracexpstr])
        #delT = (subtbl['TIME'][-1] - subtbl['TIME'][0])
        newcntr =  newcnt/delT
        newerror =  newerror/delT
     
        rateave = np.nanmean(subtbl['RATE'])
        ratestd = np.sqrt(sum(subtbl['ERROR']**2)/len(subtbl))
        dt_mean.append([newtime, newcntr, newerror, delT, newcnt, sdcnts/delT, rateave, ratestd])

    dt_mean = np.array(dt_mean, float)
    #------------------------------------------
    newtbl = Table()
    newtbl.add_column(Column(dt_mean[:,0], name = 'TIME'), index=0)
    newtbl.add_column(Column(dt_mean[:,1], name = 'RATE'), index=1)
    newtbl.add_column(Column(dt_mean[:,2], name = 'ERROR'), index=2)
    newtbl.add_column(Column(dt_mean[:,3], name = 'TBIN'), index=3)
    newtbl.add_column(Column(dt_mean[:,4], name = 'COUNTS'), index=4)
    newtbl.add_column(Column(dt_mean[:,5], name = 'STDERR'), index=5)
    newtbl.add_column(Column(dt_mean[:,6], name = 'RATEAVE'), index=6)
    newtbl.add_column(Column(dt_mean[:,7], name = 'RATESTD'), index=7)

    newtbl['TIME'].unit = 'MET'
    newtbl['RATE'].unit = 'c/s'
    newtbl['ERROR'].unit = 'c/s'
    newtbl['TBIN'].unit = 's'
    
    return tbl, newtbl 	 


def sxt_met2mjd(timevec) :
    mjdt = []
    for t in timevec :
        mjdt.append((t/86400.) + 55197.0)
    return mjdt

#-------------------------------------------------------------
def lcplotter(lcdt, exposure = 60, en = "0p3to7p0", outfile = "lcout.pdf", xlab = "xlab", ylab = "ylab", srcname = "source") :
    import numpy as np
    import os, pyfits
    from astropy.table import Table, Column
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    from matplotlib.ticker import MultipleLocator, FormatStrFormatter

    enflag = en.split("to")[0].split("p")[0] +"."+ en.split("to")[0].split("p")[1] + "-" + en.split("to")[1].split("p")[0] + "." +en.split("to")[1].split("p")[1]
 
    if len(lcdt) > 3 :
        lcdtflxm = np.nanmean(lcdt['RATE']); lcdtflsd = np.nanstd(lcdt['RATE'])
        lcdt = lcdt[lcdt['RATE'] >= (lcdtflxm - 5.5*lcdtflsd)]
        lcdt = lcdt[lcdt['RATE'] <= (lcdtflxm + 5.5*lcdtflsd)]
        lcdt['TIME'] = sxt_met2mjd(lcdt['TIME'])
        lcdt.add_column(Column(lcdt['TIME'], name = "MJD"), index = len(lcdt.columns))
        fig, ax = plt.subplots(1, sharex=True, figsize = [10,8])

        x = lcdt['TIME']-int(lcdt['TIME'][0]); y = lcdt['RATE']; ey = lcdt['ERROR']; ex = lcdt['TBIN']/(2*86400.)

        xlab = "{} - {}".format(xlab, int(lcdt['TIME'][0]))
        ylab = "{}".format(ylab)

        ax.errorbar(x, y, xerr = ex, yerr = ey, fmt="*k", markersize = 15, label = "{} keV".format(enflag))

        ylim = [y.min() - 0.3*(y.max() - y.min()),  y.max() + 0.3*(y.max() - y.min())]
        xlim = [x.min() - 0.1*(x.max() - x.min()),  x.max() + 0.2*(x.max() - x.min())]

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_xlabel(xlab, fontsize =20)
        ax.set_ylabel(ylab, fontsize =20)
        ax.set_title("SXT Lightcurve; Object : {}; Exp. : {} ks".format(srcname.upper(), round(exposure/1e3,2)), fontsize=20)
        
        ax.minorticks_on()
        #ax.tick_params(axis = 'X', which = 'major', labelsize = 12)
        #ax.tick_params(axis = 'X', which = 'minor', labelsize = 0)
        ax.tick_params(labelsize=20)
        ax.grid(which = 'major', alpha = 0.3, linewidth=1, color = 'k', linestyle = '--')
        ax.grid(which = 'minor', alpha = 0.2, linewidth=1, color = 'k', linestyle = '-.')

        ax.legend(loc=1, fontsize = 20)
        plt.savefig(outfile)
        fig.clf()
        plt.close()

        if os.path.exists(outfile) :
            status = "LCPlotter succeded in making lighcurve plots with name : \n {}".format(outfile)
        else :
            status = "Warning:: Something fissy happend could not write the curve plot"

    else :
        status = "Error:: The length of the lc data is significantly lower (=< 3), so skipping the plotting part"
    #plt.show()
    return status, lcdt


def main() :

    import numpy as np

    import matplotlib
    matplotlib.use('TKAgg')
    from astropy.table import Table, vstack
    import optparse, pyfits
    usage = "usage: %prog [options] "
    parser = optparse.OptionParser(usage)

    parser.add_option("-f", "--lcfile", dest = "lcfilename", help = "Input lc file name", default = None)

    parser.add_option("-o", "--outstr", dest = "outstr", help = "Addtional String for output", default = None)

    parser.add_option("-e", "--enstr", dest = "enstr", help = "String in format of 0p3to7p0 to convey energy band", default = "ep0to7p0")

    parser.add_option("-x", "--xstr", dest = "xstr", help = "String for first colname to be used", default = "TIME")

    parser.add_option("-y", "--ystr", dest = "ystr", help = "String for second colname to be used", default = "RATE")

    parser.add_option("-u", "--errstr", dest = "errstr", help = "String for err-axis colname to be used", default = "ERROR")

    parser.add_option("-v", "--fracexpstr", dest = "fracexpstr", help = "String the name fractiona exposures", default = "FRACEXP")

    parser.add_option("-l", "--xlab", dest = "xlab", help = "String for X-axis label", default = "TIME [MJD]")

    parser.add_option("-m", "--ylab", dest = "ylab", help = "String for Y-axis label", default = "RATE [c/s]")


    (options, args) = parser.parse_args()
    lcfilename = options.lcfilename
    outstr = options.outstr
    enstr, xlab, ylab = options.enstr, options.xlab, options.ylab 
    xstr, ystr, errstr, fracexpstr = options.xstr, options.ystr, options.errstr, options.fracexpstr

    print "===============================================\n"
    print "The input lc file name : {}\n".format(lcfilename)
    print "The string to be appended to the default output file names : {}\n".format(outstr)
    print "The string to denote the energy bands selected for lcfile {}\n".format(enstr)
    print "===============================================\n"

    f=pyfits.open(lcfilename)
    fhdr = f[0].header
    fhdu = f[1].data
    f.close()

    lc_table = fhdu[fhdu[fracexpstr]>0.75]
    binsize = float(fhdr['TIMEDEL'])
    srcname = fhdr['OBJECT']
    exposure = fhdr['EXPOSURE']
    #Source Name
    if len(srcname.split(" ")) > 1 :
        sname = srcname.split(" ")[0]
        for u in range(len(srcname.split(" "))-1) :
            sname = "{}{}".format(sname,srcname.split(" ")[u+1])
        srcname = sname
    #Output file name
    outfilename = "{}_{}_{}_avelc.png".format(srcname, enstr, outstr)
    outtblname = "{}_tbl.ipac".format(outfilename.split(".png")[0])

    thead, rhead, ehead, frachead = xstr, ystr, errstr, fracexpstr
    
    lc_table[thead] = lc_table[thead]+fhdr['TSTART']

    x, y, ey = lc_table[thead], lc_table[rhead], lc_table[ehead]
    xlim = [x.min()-0.1*(x.max()- x.min()), x.max()+0.1*(x.max()- x.min())]

    X_range, X_length = rebin_data_using_plot(x=x, y=y, ey=ey, ex=0, xlab = "X", ylab = "Y", xlim = xlim)

    tbl, newtbl = do_binning_procedure2(lc_table, X_range, posstr = thead, ratestr = rhead, errorstr = ehead, fracexpstr = frachead, size=[10.,8.75], submarzins=[0.1,0.9,0.25], xlim= xlim)

    stts, outtbl = lcplotter(newtbl, exposure = exposure, en = enstr, outfile = outfilename, xlab = xlab, ylab = ylab, srcname = srcname)


    outtbl.write(outtblname, format = "ascii.ipac")

    print "Please check files {} & {}, the output after averaging ".format(outfilename, outtblname)

if __name__ == "__main__":
    main()




