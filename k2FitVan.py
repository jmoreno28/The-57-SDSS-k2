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
import k2
import libcarma as libcarma
import util.mcmcviz as mcmcviz
import sdss as sdss
from util.mpl_settings import set_plot_params
import util.triangle as triangle
import cPickle as pickle

lowerDaysLimit = 1.0 
upperDaysLimit = 79.0 


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


try: 
	os.environ['DISPLAY']
except KeyError as Err:
	warnings.warn('No display environment! Using matplotlib backend "Agg"')
	import matplotlib
	matplotlib.use('Agg')

try:
	import carmcmc as cmcmc
except ImportError:
	carma_pack = False
else:
	carma_pack = True

fhgt = 10
fwid = 16
set_plot_params(useTex = True)

parser = argparse.ArgumentParser()
parser.add_argument('-id', '--id', type = str, default = '205905563', help = r'EPIC ID')
parser.add_argument('-p', '--processing', type = str, default = 'sap', help = r'sap/pdcsap/k2sff/k2sc/k2varcat etc...')
parser.add_argument('-c', '--campaign', type = str, default = 'c03', help = r'Campaign')
parser.add_argument('-goid', '--goid', type = str, default = '', help = r'Guest Observer ID')
parser.add_argument('-gopi', '--gopi', type = str, default = '', help = r'Guest Observer PI')
parser.add_argument('-libcarmaChain', '--lC', type = str, default = 'libcarmaChain', help = r'libcarma Chain Filename')
parser.add_argument('-nsteps', '--nsteps', type = int, default = 250, help = r'Number of steps per walker')
parser.add_argument('-nwalkers', '--nwalkers', type = int, default = 25*psutil.cpu_count(logical = True), help = r'Number of walkers')
parser.add_argument('-pMax', '--pMax', type = int, default = 1, help = r'Maximum C-AR order')
parser.add_argument('-pMin', '--pMin', type = int, default = 1, help = r'Minimum C-AR order')
parser.add_argument('-qMax', '--qMax', type = int, default = -1, help = r'Maximum C-MA order')
parser.add_argument('-qMin', '--qMin', type = int, default = -1, help = r'Minimum C-MA order')
parser.add_argument('--plot', dest = 'plot', action = 'store_true', help = r'Show plot?')
parser.add_argument('--no-plot', dest = 'plot', action = 'store_false', help = r'Do not show plot?')
parser.set_defaults(plot = True)
parser.add_argument('-minT', '--minTimescale', type = float, default = 2.0, help = r'Minimum allowed timescale = minTimescale*lc.dt')
parser.add_argument('-maxT', '--maxTimescale', type = float, default = 0.5, help = r'Maximum allowed timescale = maxTimescale*lc.T')
parser.add_argument('-maxS', '--maxSigma', type = float, default = 2.0, help = r'Maximum allowed sigma = maxSigma*var(lc)')
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

LC = k2.k2LC(name = args.id, band = 'Kep', processing = args.processing, campaign = args.campaign, goid = args.goid, gopi = args.gopi)
#/Users/Jackster/Research/k2Analysis/k2c08AGN/Vanderberg



LC.minTimescale = args.minTimescale
LC.maxTimescale = args.maxTimescale
LC.maxSigma = args.maxSigma

taskDict = dict()
DICDict= dict()

dataForResultsFile = [] 
for p in xrange(args.pMin, args.pMax + 1):

	for q in xrange(args.qMin, p):
		nt = libcarma.basicTask(p, q, nwalkers = args.nwalkers, nsteps = args.nsteps)

		print 'Starting libcarma fitting for p = %d and q = %d...'%(p, q)
		startLCARMA = time.time()
		nt.fit(LC)
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
		#np.savetxt('testwalkers.txt', nt.Chain[:,:,100:])
		#print r'Rho: ' + str(nt.rootChain[:,10, 200:])
		irand = random.randint(0 , args.nwalkers -1)
		print r'Tau: ' + str(nt.timescaleChain[:,irand, -1])
		#def vizWalkers(Chain, LogPosterior, dim1, dim1Name, dim2, dim2Name):
		#res = mcmcviz.vizWalkers(taskDict['%d %d'%(2, 2)].Chain, taskDict['%d %d'%(2, 1)].LnPosterior, 0, 1,Name, dim2, dim2Name)
		#def vizTriangle(p, q, Chain, labelList, figTitle):
		labelList = []
		for k in xrange(1, p + 1):
			labelList.append('a$_{%d}$'%(k))
		for u in xrange(0, p):
			labelList.append('b$_{%d}$'%(u))	
		print labelList
		figTitle = args.id	
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
		

sortedDICVals = sorted(DICDict.items(), key = operator.itemgetter(1))
pBest = int(sortedDICVals[0][0].split()[0])
qBest = int(sortedDICVals[0][0].split()[1])

print(sortedDICVals[0][:])
print(sortedDICVals[1][:])

pNext = int(sortedDICVals[1][0].split()[0])
qNext = int(sortedDICVals[1][0].split()[1])
print(pNext,qNext)

print 'Best model is C-ARMA(%d,%d)'%(pBest, qBest)
lnPrior = nt.logPrior(LC)
print 'The log prior is %e'%(lnPrior)
lnLikelihood = nt.logLikelihood(LC)
print 'The log likelihood is %e'%(lnLikelihood)
lnPosterior = nt.logPosterior(LC)

bestTask = taskDict['%d %d'%(pBest, qBest)]
Theta = bestTask.Chain[:,np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[0][0],np.where(bestTask.LnPosterior == np.max(bestTask.LnPosterior))[1][0]]
nt = libcarma.basicTask(pBest, qBest)
nt.set(LC.dt, Theta)
nt.smooth(LC)


#####PLOTTING
#FUNCTIONS IN CFFI/PYTHON
LCplot = LC.plot(doShow = False)
LCplot.savefig(figTitle+'_LC.png')
plt.clf()
PSDplot = nt.plotpsd(LC = LC, doShow = False)
PSDplot.savefig(figTitle+'_PSD.png')
plt.clf()
#FUNCTIONS IN LIBCARMA.PY
a = nt.plotacvf(LC = LC)
a.savefig(figTitle + '_acvf.png')
plt.clf()
b = nt.plotacf(LC = LC)
b.savefig(figTitle + '_acf.png')
plt.clf()
#c = nt.plotdacf()
#c.savefig(figTitle+'_test_dacf.png')
plt.clf()
d = nt.plotsf(LC = LC)
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

figTitle = args.id	
#plot_res =  mcmcviz.vizTriangle(p, q, nt.Chain, labelList, figTitle+'%d_%d'%(p, q))
res = mcmcviz.vizWalkers(taskDict['%d %d'%(pBest, qBest)].timescaleChain, taskDict['%d %d'%(pBest, qBest)].LnPosterior, 0, labelT[0], 1, labelT[1])
plt.clf()
res_coeffs = mcmcviz.vizTriangle(pBest, qBest, taskDict['%d %d'%(pBest, qBest)].Chain, labelList , figTitle+'Chain'+'%d_%d'%(pBest, qBest))
plt.clf()
res_times = mcmcviz.vizTriangle(pBest, qBest, taskDict['%d %d'%(pBest, qBest)].timescaleChain, labelT, figTitle+'Timescales'+'%d_%d'%(pBest, qBest))
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

res_coeffsN = mcmcviz.vizTriangle(pBest, qBest, taskDict['%d %d'%(pNext, qNext)].Chain, labelList2 , figTitle+'Chain'+'%d_%d'%(pNext, qNext))
res_timesN = mcmcviz.vizTriangle(pBest, qBest, taskDict['%d %d'%(pNext, qNext)].timescaleChain, labelT2, figTitle+'Timescales'+'%d_%d'%(pNext, qNext))


fileOut = open( "Results"  , "w" )
mW = maxWidth( dataForResultsFile) 
mS = maxSize( dataForResultsFile ) 
fileOut.write( figTitle+" "+time.time()+"\n\n") 
fileOut.write( "#p, q, timescales, DIC [repeated]\n")

for i in range( mS ):
	line = "" 
	for li in dataForResultsFile : 
		line += valOrBlank( i, li, mW ) + " "
	fileOut.write( line +"\n") 
fileOut.close() 


#pickle.dump(taskDict, open( figTitle+'%d_%d'%(pBest, qBest)+".p", "wb" ) )






if args.stop:
	pdb.set_trace()
