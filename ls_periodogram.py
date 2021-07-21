#-----------------------------------------------

import pyfits; import numpy as np; import scipy as sp
from astroML.time_series import lomb_scargle_BIC, lomb_scargle_bootstrap
from scipy.signal import lombscargle
from astroML.time_series import lomb_scargle
from gatspy.periodic import LombScargle
from gatspy.periodic import LombScargleFast
import matplotlib.pyplot as plt

# use seaborn for plot styles
import seaborn; seaborn.set()


dirStem='ASTROSAT_Ver'
paramsIn=['TIME','Roll_RA','Roll_DEC','ANG_OFFSET']

# Reading the user input data file
fileIN = pyfits.open('AS1P01_044T01_9000000072sxt_level1.mkf')
time=np.array(fileIN[1].data[:5000]['TIME'],float)
print "Number of data points used for this analysis: " + str(len(time))
delT=[]
print "estimating the smallest data-gap..."
for u in range(len(time)-1):
    delT.append(time[u+1]-time[u])

print "Generating false errors (const=1) in case not given in input array: "
dmag=np.ones(len(time))
dy=np.zeros(len(time))
# Choose a period grid
print "Making period vector based on the minum data gaps and number of data points...."
periods = np.linspace(2*min(delT), int((time.max()-time.min())/3.0), 400000)
ang_freqs = 2 * np.pi / periods

#----opening files to write output
fileout =open('OutFile1.out','w')

for i in range(len(paramsIn)-1):
    # compute the (unnormalized) periodogram
    # note pre-centering of y values!
    t=np.array(fileIN[1].data[:5000][paramsIn[0]],float)
    mag=np.array(fileIN[1].data[:5000][paramsIn[i + 1]],float)
    print "The TS being used is for: "+ str(paramsIn[i + 1])
    
    #----Using Scipy.signal
    print "Performing the periodogram analysis using GatSpy"
    power = lombscargle(t, mag - mag.mean(), ang_freqs)
    # normalize the power
    print "normalize the power"
    power *= 2 / (len(t) * mag.std() ** 2)
    
    #------------------------------------------------------------
    # Plot the results
    print "plotiing the results------"
    fig = plt.figure(figsize=(10, 6.75))
    fig.subplots_adjust(left=0.1, right=0.9, hspace=0.25)

    # First panel: the data
    ax = fig.add_subplot(211)
    ax.errorbar(t, mag, dy, fmt='.k', lw=1, ecolor='gray')
    ax.set_xlabel('time [s]')
    ax.set_ylabel(paramsIn[i + 1])
    ax.set_xlim(t.min()-(t.max()-t.min())*0.1, t.max()+(t.max()-t.min())*0.1)

    # Second panel: the periodogram & significance levels
    ax1 = fig.add_subplot(212, xscale='log')
    ax1.plot(periods, power, '-', c='black', lw=1.5, zorder=1)
#    ax1.plot([period[0], period[-1]], [sig1, sig1], ':', c='black')
#    ax1.plot([period[0], period[-1]], [sig5, sig5], ':', c='black')

#    ax1.annotate("", (0.3, 0.65), (0.3, 0.85), ha='center',
#             arrowprops=dict(arrowstyle='->'))

    ax1.set_xlim(periods[0], periods[-1])
    ax1.set_ylim(power.min()-((power.max()-power.min())*0.05), power.max()+((power.max()-power.min())*0.05))

    ax1.set_xlabel(r'period (s)')
    ax1.set_ylabel('power')
#    imageOuttxt=(dirStem)'TSstudyOut_'(paramsIn[i + 1])
    print "Writing plots to files"
    plt.savefig(str(dirStem) + 'TSstudyOut_' + str(paramsIn[i + 1]) + '.png')
    plt.savefig(str(dirStem) + 'TSstudyOut_' + str(paramsIn[i + 1]) + '.eps')
    plt.savefig(str(dirStem) + 'TSstudyOut_' + str(paramsIn[i + 1]) + '.pdf')

    #-----Using astroML
#    power_astroML = lomb_scargle(t, mag, dmag, ang_freqs)
#    D=lomb_scargle_bootstrap(t,mag,dmag,ang_freqs,generalized=True,N_bootstraps=1000,random_state=0)
#    sig1, sig5 = np.percentile(D, [99,95])
    
    #-----Using Gatspy LombScargle
#    model1 = LombScargle(fit_offset=True).fit(t, mag, dmag)
#    power1 = model1.score(periods)
    
    #-----Using Gatspy LombScargleFast
#    fmin = 1. / periods.max()
#    fmax = 1. / periods.min()
#    N = 10000
#    df = (fmax - fmin) / N

#    model2 = LombScargleFast().fit(t, mag, dmag)
#    power2 = model2.score_frequency_grid(fmin, df, N)
#    freqs2 = fmin + df * np.arange(N)

#    period, power = model.periodogram_auto(nyquist_factor=100)
#    print("period range: ({0}, {1})".format(period.min(), period.max()))
#    print("number of periods: {0}".format(len(period)))

#    dtpow=np.column_stack((periods,power,power_astroML,power1)v)
    dtpow = [None]*len(periods)
    if i == 0:
        dtpow[i]=np.column_stack((periods,power))
        print "making two Dimentional array for output..."
    else:
        dtpow[i]=np.column_stack((dtpow[i-1],power))
        print "Adding one more column to the output data matrix... "


#--------Writing Outputfiles....
print "Writing Output to Files/..../"
fileOut.write(dtpow[1])
fileOut.close()
fileIN.close()
