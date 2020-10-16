#!/usr/bin/env python3
################################################
################################################

vmecFiles    = ['wistella_midscale' ,'NuhrenbergZille_1988_QHS','HSX_QHS_vacuum_ns201','n4qh.b4.a79a','Drevlak_qh_8_7','li383_1.4m_ns201','n3are_R7.75B5.7_hires','GarabedianQAS2_noCurrentOnAxis_ns201','estell_24_scaled','cfqs_freeBoundary_vacuum_hiRes','st_a34_i32v22_beta_35_scaledAUG_hires']
stellDesigns = ['WISTELL-A'         ,'NZ1988'                  ,'HSX'                 ,'KuQHS48'     ,'Drevlak'       ,'NCSX'            ,'ARIES-CS'             ,'QAS2'                                ,'ESTELL'          ,'CFQS'                          ,'Henneberg'   ]
etabar       = [0.791               ,0.157                     ,1.28                  ,0.147         ,0.0899          ,0.408             ,0.0740                 ,0.347                                 ,0.570             ,0.586                           ,0.302         ]
B0           = [2.54                ,0.205                     ,1.00                  ,1.20          ,3.97            ,1.55              ,5.69                   ,1.79                                  ,1.00              ,0.933                           ,2.41          ]
nzgridvec    = [300                 ,300                       ,300                   ,300           ,300             ,300               ,300                    ,300                                   ,300               ,300                             ,300           ]
nlambdavec   = [35                  ,35                        ,35                    ,35            ,35              ,35                ,35                     ,45                                    ,35                ,35                              ,35            ]
NFPtoInclude = [3                   ,2.5                       ,2.5                   ,2.5           ,2.5             ,2.5               ,2.5                    ,2.5                                   ,2.5               ,1.7                             ,1.7           ]
negridvec    = [12                  ,12                        ,12                    ,12            ,12              ,12                ,12                     ,12                                    ,12                ,12                              ,12            ]
deltaTvec    = [0.4                 ,0.4                       ,0.4                   ,0.4           ,0.4             ,0.4               ,0.4                    ,0.3                                   ,0.4               ,0.4                             ,0.4           ]
simTimevec   = [180                 ,180                       ,180                   ,180           ,180             ,180               ,180                    ,180                                   ,180               ,180                             ,180           ]
ngaussvec    = [3                   ,3                         ,3                     ,3             ,3               ,3                 ,3                      ,3                                     ,3                 ,3                               ,3             ]

tprimvec = [0,1,2,3,4,5,6,7,8] # -|grad(T)/T|
fprimvec = [0,1,2,3,4,5,6,7,8] # |grad(T)/T|\|grad(n)/n|
aky_min  = 0.01 # minimum ky to consider
aky_max  = 8.0 # maximum ky to consider
naky     = 20 # number of ky points
rr  = 0.01 # normalized flux psi/psi_a radial coordinate
stellsToRun  = [0,1,2,3,4,5,6,7,8,9,10] # select the indices of stellarators from stellDesigns to analyze. Start at 0
fractionToConsider = 0.7 # fraction of time fron the simulation period to consider

makeCleanVMEC2GS2 = 0 # 1 - make clean and then make VMEC2GS2, 0 - just make
ncores = 8 # number of cpu cores to run GS2
runGS2 = 1 # 1 -> run gs2; 0 -> don't run gs2 and only plot
plotFig = 0 # 1 -> save plots in pdf; 0 -> don't


vmecGS2interfaceFolder='/Users/rogeriojorge/Dropbox/PostDoc/NearAxisGyrokinetics/ConvergenceTests/VMEC_to_GS2'
runsPath='/Users/rogeriojorge/scratch/GS2critRuns/'
equilibriaFolder='/Users/rogeriojorge/Dropbox/PostDoc/NearAxisGyrokinetics/equilibria/'
gs2Path='/Users/rogeriojorge/local/gs2/bin'
wolframScript = "wolframscript"

import os
from os import path, fdopen, remove
from scipy.io import netcdf
import subprocess
from shutil import move, copymode, copyfile, rmtree
from tempfile import mkstemp
import matplotlib.pyplot as plt
import numpy as np

## User-defined functions from criticalFuncs.py
from criticalFuncs import replace, removeGS2, runGS2func, gammabyky

######################################
## Create executable for VMEC_2_GS2 ##
######################################
currentPath = os.getcwd()
os.chdir(vmecGS2interfaceFolder)
if makeCleanVMEC2GS2==1:
	print("Compiling VMEC to GS2")
	process = subprocess.call(['make','clean'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
process = subprocess.call(['make'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
copyfile("test_vmec_to_gs2_geometry_interface",currentPath+"/vmec2gs2")
copymode("test_vmec_to_gs2_geometry_interface",currentPath+"/vmec2gs2")
os.chdir(currentPath)
##############################

if not path.exists(runsPath):
	os.mkdir(runsPath)

for i in stellsToRun:
	# Select grid variables from input vectors
	stells = stellDesigns[i]
	nzgrid = nzgridvec[i]
	nlambda = nlambdavec[i]
	nfp = NFPtoInclude[i]
	ne = negridvec[i]
	dt = deltaTvec[i]
	nstep = int(simTimevec[i]/deltaTvec[i])
	ngauss = ngaussvec[i]
	etab = etabar[i]
	b0 = B0[i]
	equilibrium = equilibriaFolder+'wout_'+vmecFiles[i]+'.nc'

	# Create paths
	if not path.exists(runsPath+stells):
		os.mkdir(runsPath+stells)
	if not path.exists(runsPath+stells+'/r'+str(rr)):
		os.mkdir(runsPath+stells+'/r'+str(rr))
	GSpath = runsPath+stells+'/r'+str(rr)+'/'

	######################
	## Create VMEC grid ##
	######################
	gs2grid      = 'gs2_grid'+stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	gxgrid       = 'gx_grid'+stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	geom2math    = 'gridMath'+stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	f = netcdf.netcdf_file(equilibrium,'r',mmap=False)
	iotaVMEC = f.variables['iotaf'][()][1]
	nfpVMEC  = f.variables['nfp'][()]
	if not path.exists(GSpath+gs2grid):
		print(' Running VMEC_2_GS2')
		process = subprocess.call(['./vmec2gs2',equilibrium,GSpath+gs2grid,GSpath+geom2math,str(nfpVMEC*nfp/abs(iotaVMEC)),str(rr),GSpath+gxgrid,str(nzgrid),str(nlambda)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	########################
	## Get Near-Axis Grid ##
	########################
	gs2gridNA   = 'gs2_NAgrid'+stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	gxgridNA    = 'gx_NAgrid' +stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	sigmaFile   = runsPath+stells+'/Math'+stells+'.mx'
	bashCommand = wolframScript+" -noprompt -script solveSigma.m "+str(etab)+" "+str(b0)+" "+equilibrium+" "+stells+" "+GSpath+geom2math+" "+sigmaFile+" "+GSpath+gs2gridNA+" "+GSpath+gxgridNA
	if not path.exists(GSpath+gs2gridNA):
		print(' Running Mathematica')
		output = subprocess.call(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	#### Iterate over tprim and eta ####
	kysX  = np.zeros(shape=(len(tprimvec),len(fprimvec)))
	kysNA = np.zeros(shape=(len(tprimvec),len(fprimvec)))
	growthRatesX  = np.zeros(shape=(len(tprimvec),len(fprimvec)))
	growthRatesNA = np.zeros(shape=(len(tprimvec),len(fprimvec)))
	realFrequenciesX  = []
	realFrequenciesNA = []
	xdata = []
	ydata = []
	for countx, tprim in enumerate(tprimvec):
		for county, fprim in enumerate(fprimvec):
			print('Run with tprim='+str("%.2f" % tprim)+' and fprim='+str("%.2f" % fprim))
			#############
			## Run GS2 ##
			#############
			if not path.exists(GSpath+'gs2'):
				os.mkdir(GSpath+'gs2')
			copyfile(gs2Path+"/gs2",GSpath+"gs2/gs2")
			copymode(gs2Path+"/gs2",GSpath+"gs2/gs2")
			runName=stells+'tprim'+str("%.2f" % tprim)+'fprim'+str("%.2f" % fprim)+'.in'
			runNameNA=runName[:-3]+'NA.in'
			runGS2func(plotFig,currentPath,GSpath,runName,runNameNA,gs2grid,gxgrid,gs2gridNA,gxgridNA,stells,nzgrid,ne,dt,nstep,ngauss,tprim,fprim,aky_min,aky_max,naky,ncores,runGS2)
			## Compute ky, growth rate and frequency per ky
			kyX, kyNA, growthRateX, growthRateNA, realFrequencyX, realFrequencyNA = gammabyky(plotFig,currentPath,GSpath,runName,runNameNA,fractionToConsider)
			kysX[countx,county ]=kyX[growthRateX.index(max(growthRateX))]
			kysNA[countx,county]=kyNA[growthRateNA.index(max(growthRateNA))]
			realFrequenciesX[countx,county ]=realFrequencyX[growthRateX.index(max(growthRateX))]
			realFrequenciesNA[countx,county]=realFrequencyNA[growthRateNA.index(max(growthRateNA))]
			growthRatesX[countx,county] =np.max(growthRateX)
			growthRatesNA[countx,county]=np.max(growthRateNA)
			xdata.append(fprim)
			ydata.append(tprim)
	if plotFig==1:
		## Plot maximum growth rate and corresponding ky for each fprim and tprim
		print('Plotting max gamma and ky for each tprim and fprim')
		plotExtent=[min(fprimvec),max(fprimvec),min(tprimvec),max(tprimvec)]
		finalPlot(growthRatesX,growthRatesNA,plotExtent,'gamma',stells)
		finalPlot(realFrequenciesX,realFrequenciesNA,plotExtent,'omega',stells)
		finalPlot(kysX,kysNA,plotExtent,'ky',stells)

