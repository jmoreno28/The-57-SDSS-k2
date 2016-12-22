import numpy as np
import math as math
import cmath as cmath
import psutil as psutil
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib import gridspec as gridspec
import argparse as argparse
import operator as operator
import warnings as warnings
import copy as copy
import time as time
import pdb
import os as os
import random

import kali.k2
import kali.s82
import kali.carma
import kali.util.mcmcviz as mcmcviz
from kali.util.mpl_settings import set_plot_params
import kali.util.triangle as triangle
import pickle


def doubleMADsfromMedian(y,thresh=3.5):
    # warning: this function does not check for NAs
    # nor does it address issues when 
    # more than 50% of your data have identical values
    m = np.median(y)
    abs_dev = np.abs(y - m)
    left_mad = np.median(abs_dev[y <= m])
    right_mad = np.median(abs_dev[y >= m])
    y_mad = left_mad * np.ones(len(y))
    y_mad[y > m] = right_mad
    modified_z_score = 0.6745 * abs_dev / y_mad
    modified_z_score[y == m] = 0
    return np.where(modified_z_score < thresh)[0]
    
# make new kepler flux
def kepFlux(sdsslc_r,sdsslc_g):
    sdss_g = sdsslc_g.y
    sdss_r = sdsslc_r.y
    sdss_t = sdsslc_g.t
    sdss_gerr = sdsslc_g.yerr
    
    if len(sdss_r) != len(sdss_g):
        print sdsslc_r.t[-1], sdsslc_g.t[-1]
        #catch missing g band measurements in the middle of the array 
        for i in range(0, len(sdss_t)):
            tolerance = 0.0002
            missed = np.isclose(sdsslc_r.t, sdss_t[i],tolerance)
            m = np.where(missed == True)[0]
            #print m, i
            if m != i:
                print m, i, len(sdss_t)
                sdss_g = np.insert(sdss_g, m ,0.)
                sdss_gerr = np.insert(sdss_gerr, m ,0.)
                sdss_t = np.insert(sdss_t, m ,sdsslc_r.t[m])
                print sdss_g[m], sdss_t[i], len(sdss_t),len(sdsslc_r.t)
            if len(sdss_r) == len(sdss_g):
                break
    c = sdss_r*0.8 + sdss_g*0.2
    fullr = np.where(sdss_g == 0.)[0]
    c[fullr] = sdss_r[fullr]
    print(c[fullr],sdss_r[fullr-1])
    c_err = np.zeros(len(sdss_t))
    for i in range (len(sdss_t)):
    	#c_err[i] = np.maximum(sdsslc_r.yerr[i], sdss_gerr[i])
    	
    	#compute error by adding in quadrature
    	thing1 = 0.8
    	thing2 = 0.2
    	c_err = np.sqrt( (thing1**2. * sdsslc_r.yerr**2.) + (thing2**2. * sdss_gerr**2.))

    c_t = sdss_t
    #return 0,0,0                   
    return c, c_err, c_t
    
def valOrBlank( i, x, size=2 ):
	if i >= len(x):
		return " "*size 
	else:
		value = x[i] 
		remainder = size - len( str(x[i] ) ) 
		outStr = str(value)
		if remainder > 0: 
			outStr += ( remainder*" ") 
		return outStr
def maxSize( listOfLists ):
	maxSizeOut = 0 
	for li in listOfLists :
		if len( li ) > maxSizeOut : 
			maxSizeOut = len(li ) 
	return maxSizeOut
def maxWidth( listOfLists ):
	maxWidthOut = 0 
	for li in listOfLists:
		for q in li : 
			if len( str(q ) ) > maxWidthOut:
				maxWidthOut = len( str(q) ) 
	return maxWidthOut 
#---------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-id', '--id', type = str, default = '220173631', help = r'EPIC ID')
parser.add_argument('-sdssid', '--sdssid', type = str, default = '005102.42-010244.4', help = r'EPIC ID')
parser.add_argument('-p', '--processing', type = str, default = 'k2sff', help = r'sap/pdcsap/k2sff/k2sc/k2varcat etc...')
parser.add_argument('-c', '--campaign', type = str, default = 'c08', help = r'Campaign')
parser.add_argument('-goid', '--goid', type = str, default = '', help = r'Guest Observer ID')
parser.add_argument('-gopi', '--gopi', type = str, default = '', help = r'Guest Observer PI')
parser.add_argument('-libcarmaChain', '--lC', type = str, default = 'libcarmaChain', help = r'libcarma Chain Filename')
parser.add_argument('-nsteps', '--nsteps', type = int, default = 1000, help = r'Number of steps per walker')
parser.add_argument('-nwalkers', '--nwalkers', type = int, default = 25*psutil.cpu_count(logical = True), help = r'Number of walkers')
parser.add_argument('-pMax', '--pMax', type = int, default = 4, help = r'Maximum C-AR order')
parser.add_argument('-pMin', '--pMin', type = int, default = 1, help = r'Minimum C-AR order')
parser.add_argument('-qMax', '--qMax', type = int, default = -1, help = r'Maximum C-MA order')
parser.add_argument('-qMin', '--qMin', type = int, default = -1, help = r'Minimum C-MA order')
parser.add_argument('--plot', dest = 'plot', action = 'store_true', help = r'Show plot?')
parser.add_argument('--no-plot', dest = 'plot', action = 'store_false', help = r'Do not show plot?')
parser.set_defaults(plot = True)
parser.add_argument('-minT', '--minTimescale', type = float, default = 2.0, help = r'Minimum allowed timescale = minTimescale*lc.dt')
parser.add_argument('-maxT', '--maxTimescale', type = float, default = 0.5, help = r'Maximum allowed timescale = maxTimescale*lc.T')
parser.add_argument('-maxS', '--maxSigma', type = float, default = 5.0, help = r'Maximum allowed sigma = maxSigma*var(lc)')
parser.add_argument('--stop', dest = 'stop', action = 'store_true', help = r'Stop at end?')
parser.add_argument('--no-stop', dest = 'stop', action = 'store_false', help = r'Do not stop at end?')
parser.set_defaults(stop = False)
parser.add_argument('--save', dest = 'save', action = 'store_true', help = r'Save files?')
parser.add_argument('--no-save', dest = 'save', action = 'store_false', help = r'Do not save files?')
parser.set_defaults(save = False)
parser.add_argument('--log10', dest = 'log10', action = 'store_true', help = r'Compute distances in log space?')
parser.add_argument('--no-log10', dest = 'log10', action = 'store_false', help = r'Do not compute distances in log space?')
parser.set_defaults(log10 = False)
parser.add_argument('--viewer', dest = 'viewer', action = 'store_true', help = r'Visualize MCMC walkers')
parser.add_argument('--no-viewer', dest = 'viewer', action = 'store_false', help = r'Do not visualize MCMC walkers')
parser.set_defaults(viewer = True)
args = parser.parse_args()

if (args.qMax >= args.pMax):
	raise ValueError('pMax must be greater than qMax')
if (args.qMax == -1):
	args.qMax = args.pMax - 1
if (args.qMin == -1):
	args.qMin = 0
if (args.pMin < 1):
	raise ValueError('pMin must be greater than or equal to 1')
if (args.qMin < 0):
	raise ValueError('qMin must be greater than or equal to 0')
	
lowerDaysLimit = 1.0 
upperDaysLimit = 3100.0 	
#------------------------------------------------------------------------
#Load target list ids



#Loading SDSS g and r band flux lightcurves: object #1
sdsslc_r = kali.s82.sdssLC(name = args.sdssid, band = 'r')
sdsslc_g = kali.s82.sdssLC(name = args.sdssid, band = 'g')


import seaborn as sns
c, c_err, c_t = kepFlux(sdsslc_r,sdsslc_g)
plt.scatter(c_t, c)
plt.savefig('TestCombineLC.png')
#plt.errorbar(sdsslc_r.t, c, yerr = c_err)
print c_err
print c
print sdsslc_r.y
plt.clf()
#sanity check: Is the lightcurve normalized? Answer : no
c_med = np.median(c)
c_norm = c/c_med
print("the average flux is %d",c_med)


full_lcy = c_norm
full_lcyerr = c_err/c_med
full_lct = c_t
w = np.where(full_lcy > 0)[0]
#compare time attibutes in sdss and k2 objects from kali
#print(len(sdsslc_g.t), len(k2lc.t),len(full_lct))
plt.plot(full_lct[w], full_lcy[w])
plt.savefig(args.sdssid+'combined.png')
plt.clf()
#match mask index
#mmatch errors and deal with outliers
y = full_lcy[w]
yt = full_lct[w]
z = doubleMADsfromMedian(y )
plt.plot(yt[z], y[z])
plt.savefig(args.sdssid+'g_r_noOutliers.png')


#reset K2LC
sdsslc_r.y = y

sdsslc_r.yerr = full_lcyerr[w]

sdsslc_r.t = yt

sdsslc_r.mask = np.zeros(len(sdsslc_r.y))
sdsslc_r.mask[z] = 1
sdsslc_r.cadence = np.arange(0,len(sdsslc_r.y))
plt.scatter(sdsslc_r.t , sdsslc_r.y*sdsslc_r.mask)
plt.savefig(args.sdssid+'nonzero_noOutliers.png')

#correct other timescale params 
#sigma = maxSigma*var(lc)
sdsslc_r.startT = 0.
sdsslc_r._dt = 100.0## Increment between epochs.
sdsslc_r._mindt = 1.0
sdsslc_r._maxdt = 300.0
sdsslc_r._T =sdsslc_r.t[-1] - sdsslc_r.t[0] ## Total duration of the light curve.
sdsslc_r._numCadences = len(sdsslc_r.y)
#self.cadence = 0.02



taskDict = dict()
DICDict= dict()
# if statement to check for pickle and load it
dataForResultsFile = [] 
for p in xrange(args.pMin, args.pMax + 1):

	for q in xrange(args.qMin, p):
		nt = kali.carma.CARMATask(p, q, nwalkers = args.nwalkers, nsteps = args.nsteps)

		print 'Starting libcarma fitting for p = %d and q = %d...'%(p, q)
		startLCARMA = time.time()
		nt.fit(sdsslc_r)
		stopLCARMA = time.time()
		timeLCARMA = stopLCARMA - startLCARMA
		print 'libcarma took %4.3f s = %4.3f min = %4.3f hrs'%(timeLCARMA, timeLCARMA/60.0, timeLCARMA/3600.0)

		Deviances = copy.copy(nt.LnPosterior[:,args.nsteps/2:]).reshape((-1))
		DIC = 0.5*math.pow(np.std(-2.0*Deviances),2.0) + np.mean(-2.0*Deviances)
		print 'C-ARMA(%d,%d) DIC: %+4.3e'%(p, q, DIC)
		DICDict['%d %d'%(p, q)] = DIC
		taskDict['%d %d'%(p, q)] = nt
		#save Theta vector for walker 10 at step 15
		#print r'Theta: ' + str(nt.Chain[:,10, 15])

		#print r'Rho: ' + str(nt.rootChain[:,10, 200:])
		irand = random.randint(0 , args.nwalkers -1)
		print r'Tau: ' + str(nt.timescaleChain[:,irand, -1])

		labelList = []
		for k in xrange(1, p + 1):
			labelList.append('a$_{%d}$'%(k))
		for u in xrange(0, p):
			labelList.append('b$_{%d}$'%(u))	
		print labelList
		figTitle = args.sdssid	
		#plot_res =  mcmcviz.vizTriangle(p, q, nt.Chain, labelList, figTitle+'%d_%d'%(p, q))
		iRandom = random.randint( 0, args.nwalkers -1  )
		#print r'Tau: ' + str(nt.timescaleChain[:, iRandom,  -1])
		timescalesRough = nt.timescaleChain[:, iRandom, -1] 
		pAllSame = [ ] 
		qAllSame = [ ] 
		timescales = [] 
		for t in timescalesRough:
			if t >= lowerDaysLimit and t <= upperDaysLimit : 
				timescales.append( t ) 	
		DICAllSame = [ ] 
		for tt in timescales : 
			DICAllSame.append( DIC ) 
			pAllSame.append( p ) 
			qAllSame.append(q ) 
		dataForResultsFile.append( pAllSame ) 
		dataForResultsFile.append( qAllSame ) 
		dataForResultsFile.append( timescales ) 
		dataForResultsFile.append( DICAllSame) 
		fname = str(figTitle+'%i-%i'%(p, q)+'TimescaleSamples')

		output = open(fname+'.pkl', 'wb')
		pickle.dump(np.array(taskDict['%i %i'%(p, q)].timescaleChain),output)
		output.close()
		output2 = open(fname+'.pkl', 'wb')
		pickle.dump(np.array(taskDict['%i %i'%(p, q)].Chain),output2)
		output2.close()

sortedDICVals = sorted(DICDict.items(), key = operator.itemgetter(1))
pBest = int(sortedDICVals[0][0].split()[0])
qBest = int(sortedDICVals[0][0].split()[1])

print(sortedDICVals[0][:])
print(sortedDICVals[1][:])

pNext = int(sortedDICVals[1][0].split()[0])
qNext = int(sortedDICVals[1][0].split()[1])
print(pNext,qNext)



print 'Best model is C-ARMA(%d,%d)'%(pBest, qBest)
lnPrior = nt.logPrior(sdsslc_r)
print 'The log prior is %e'%(lnPrior)
lnLikelihood = nt.logLikelihood(sdsslc_r)
print 'The log likelihood is %e'%(lnLikelihood)
lnPosterior = nt.logPosterior(sdsslc_r)

bestTask = taskDict['%d %d'%(pBest, qBest)]
Theta = bestTask.Chain[:,np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[0][0],np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[1][0]]
nt = kali.carma.CARMATask(pBest, qBest)
nt.set(sdsslc_r.dt, Theta)
nt.smooth(sdsslc_r)


#####PLOTTING
#FUNCTIONS IN CFFI/PYTHON
LCplot = sdsslc_r.plot(doShow = False)
LCplot.savefig(figTitle+'_LC.png')
plt.clf()
PSDplot = nt.plotpsd(LC = sdsslc_r, doShow = False)
PSDplot.savefig(figTitle+'_PSD.png')
plt.clf()
#FUNCTIONS IN LIBCARMA.PY
a = nt.plotacvf(LC = sdsslc_r)
a.savefig(figTitle + '_acvf.png')
plt.clf()
b = nt.plotacf(LC = sdsslc_r)
b.savefig(figTitle + '_acf.png')
plt.clf()
#c = nt.plotdacf()
#c.savefig(figTitle+'_test_dacf.png')
plt.clf()
d = nt.plotsf(LC = sdsslc_r)
d.savefig(figTitle+'_sf.png')
plt.clf()

#lagsEst, sfEst, sferrEst = LC.sf()
#lagsModel, sfModel = bestTask.sf(start = lagsEst[1], stop = lagsEst[-1], num = 5000, spacing = 'log')
#plt.clf()
#plt.loglog(lagsModel, sfModel, label = '$SF(\delta t)$ (model)', color = '#000000', zorder = 5)
#plt.errorbar(lagsEst, sfEst, sferrEst, label = '$SF(\delta t)$ (est)', fmt = 'o', capsize = 0, color = '#ff7f00', markeredgecolor = 'none', zorder = 0)
#plt.xlabel('$\log_{10}\delta t$')
#plt.ylabel('$\log_{10} SF$')
#plt.legend(loc = 2)
#plt.savefig(figTitle+'_sf.png')
#plt.clf()
#---------------------------------

#LagsEst, acfEst = LC.acvf()
#Lags, ACVF = bestTask.acvf(1.0e-4, 1.0e4, 10000, spacing = 'log')
#plt.plot(np.log10(Lags), ACVF, label = '$ACVF$ (model)', color = '#000000')
#plt.plot(np.log10(LagsEst), acfEst, label = '$ACVF(\delta t)$ (est)', fmt = 'o', capsize = 0, color = '#ff7f00', markeredgecolor = 'none')
#plt.xlabel(r'$\log_{10}\delta t$')
#plt.ylabel(r'$\log_{10} ACVF$')
#plt.legend(loc = 2)
#plt.savefig(figTitle+'acvf.png')
#plt.clf()


#----------------plotting timescale triangle

labelList = []
labelT = []
labelRoots = []
for k in xrange(1, pBest + 1):
	labelList.append('a$_{%d}$'%(k))
	labelT.append('$\tao_{%d}$'%(k))
	labelRoots.append('r$_{%d}$'%(k))
for u in xrange(0, pBest):
	labelList.append('b$_{%d}$'%(u))	
	labelT.append('$\tao_{MA%d}$'%(u))
	labelRoots.append('r$_{MA%d}$'%(u))

figTitle = args.sdssid	
#plot_res =  mcmcviz.vizTriangle(p, q, nt.Chain, labelList, figTitle+'%d_%d'%(p, q))
res = mcmcviz.vizWalkers(taskDict['%d %d'%(pBest, qBest)].timescaleChain, taskDict['%d %d'%(pBest, qBest)].LnPosterior, 0, labelT[0], 1, labelT[1])
plt.clf()
res_coeffs = mcmcviz.vizTriangle(pBest, qBest, taskDict['%d %d'%(pBest, qBest)].Chain, labelList , str(figTitle+'BestChain'+'%d-%d'%(pBest, qBest)))
plt.clf()
res_times = mcmcviz.vizTriangle(pBest, qBest, taskDict['%d %d'%(pBest, qBest)].timescaleChain, labelT, str(figTitle+'BestTimescales'+'%d-%d'%(pBest, qBest)))

plt.clf()

#reset labels list
labelList2 = []
labelT2 = []
labelRoots2 = []
for k in xrange(1, pNext + 1):
	labelList2.append('a$_{%d}$'%(k))
	labelT2.append('$tao_{%d}$'%(k))
	labelRoots2.append('r$_{%d}$'%(k))
for u in xrange(0, pNext):
	labelList2.append('b$_{%d}$'%(u))	
	labelT2.append('$tao_{MA%d}$'%(u))
	labelRoots2.append('r$_{MA%d}$'%(u))	
#vizTriangle(p, q, Chain, labelList, figTitle, doShow=False)
res_coeffsN = mcmcviz.vizTriangle(pBest, qBest, taskDict['%d %d'%(pNext, qNext)].Chain, labelList2 , str(figTitle+'Chain'+'%d-%d'%(pNext, qNext)),doShow=False)
res_timesN = mcmcviz.vizTriangle(pBest, qBest, taskDict['%d %d'%(pNext, qNext)].timescaleChain, labelT2, str(figTitle+'Timescales'+'%d-%d'%(pNext, qNext)),doShow=False)


#fileOut = open( "Results"  , "w" )
#mW = maxWidth( dataForResultsFile) 
#mS = maxSize( dataForResultsFile ) 
#fileOut.write( figTitle+" "+time.time()+"\n\n") 
#fileOut.write( "#p, q, timescales, DIC [repeated]\n")

#for i in range( mS ):
#	line = "" 
#	for li in dataForResultsFile : 
#		line += valOrBlank( i, li, mW ) + " "
#	fileOut.write( line +"\n") 
#fileOut.close() 







if args.stop:
	pdb.set_trace()



