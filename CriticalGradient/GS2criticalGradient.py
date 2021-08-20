#!/usr/bin/env python3
################################################
################################################

vmecFiles    = ['wistella_midscale' ,'NuhrenbergZille_1988_QHS','HSX_QHS_vacuum_ns201','n4qh.b4.a79a','Drevlak_qh_8_7','li383_1.4m_ns201','n3are_R7.75B5.7_hires','GarabedianQAS2_noCurrentOnAxis_ns201','estell_24_scaled','cfqs_freeBoundary_vacuum_hiRes','st_a34_i32v22_beta_35_scaledAUG_hires','LandremanSengupta2019_section5.4','ITER']
stellDesigns = ['WISTELL-A'         ,'NZ1988'                  ,'HSX'                 ,'KuQHS48'     ,'Drevlak'       ,'NCSX'            ,'ARIES-CS'             ,'QAS2'                                ,'ESTELL'          ,'CFQS'                          ,'Henneberg'                            ,'NAQS'                            ,'ITER']
etabar       = [0.791               ,0.1567                    ,1.28                  ,0.147         ,0.08914         ,0.408             ,0.0740                 ,0.347                                 ,0.570             ,0.586                           ,0.302                                  ,1.549                             ,0.115]
B0           = [2.54                ,0.2052                    ,1.00                  ,1.20          ,3.9736          ,1.55              ,5.69                   ,1.79                                  ,1.00              ,0.933                           ,2.41                                   ,1.001                             ,5.358]
nzgridvec    = [150                 ,162                       ,150                   ,150           ,182             ,150               ,150                    ,150                                   ,150               ,150                             ,150                                    ,100                               ,150]
nlambdavec   = [34                  ,26                        ,24                    ,24            ,28              ,24                ,24                     ,30                                    ,24                ,28                              ,32                                     ,24                                ,20]
NFPtoInclude = [3                   ,3                         ,4                     ,6             ,3               ,4                 ,6                      ,7                                     ,2                 ,2                               ,6                                      ,4                                 ,3]
negridvec    = [10                  ,10                        ,10                    ,10            ,10              ,10                ,10                     ,10                                    ,10                ,10                              ,10                                     ,10                                ,8]
deltaTvec    = [0.4                 ,0.4                       ,0.4                   ,0.4           ,0.4             ,0.4               ,0.4                    ,0.4                                   ,0.4               ,0.4                             ,0.4                                    ,0.4                               ,0.4]
simTimevec   = [70                  ,80                        ,60                    ,60            ,120             ,60                ,140                    ,100                                   ,70                ,60                              ,60                                     ,60                                ,50]
ngaussvec    = [3                   ,3                         ,3                     ,3             ,3               ,3                 ,3                      ,3                                     ,3                 ,3                               ,3                                      ,3                                 ,3]
QAorQH       = ['QH'                ,'QH'                      ,'QH'                  ,'QH'          ,'QH'            ,'QA'              ,'QA'                   ,'QA'                                  ,'QA'              ,'QA'                            ,'QA'                                   ,'QH'                              ,'A']
## 11 stellarators + 1 tokamak

tprimvec = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6]#[0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6] # -|grad(T)/T|
fprimvec = [0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5,6] # |grad(T)/T|\|grad(n)/n|
aky_min  = 0.01 # minimum ky to consider
aky_max  = 20 # maximum ky to consider
naky     = 40 # number of ky points
rr  = 0.01 # normalized flux psi/psi_a radial coordinate
stellsToRun  = [0,1,2,3,4,5,6,7,8,9,10,11] # select the indices of stellarators from stellDesigns to analyze. Start at 0
fractionToConsider = 0.6 # fraction of time fron the simulation period to consider

makeCleanVMEC2GS2 = 0 # 1 - make clean and then make VMEC2GS2, 0 - just make
ncores = 8 # number of cpu cores to run GS2
runGS2 = 1 # 1 -> run gs2; 0 -> don't run gs2 and only plot
plotFig = 1 # 1 -> save plots in pdf; 0 -> don't

vmecGS2interfaceFolder='/Users/rogeriojorge/Dropbox/PostDoc/NearAxisGyrokinetics/ConvergenceTests/VMEC_to_GS2'
runsPath='/Users/rogeriojorge/scratch/GS2critRuns/'
equilibriaFolder='/Users/rogeriojorge/Dropbox/PostDoc/NearAxisGyrokinetics/equilibria/'
gs2Path='/Users/rogeriojorge/local/gs2/bin'
sigmaScript='/Users/rogeriojorge/Dropbox/Postdoc/NearAxisGyrokinetics/solveSigma.py'
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
from criticalFuncs import replace, removeGS2, runGS2func, gammabyky, finalPlot

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
	print('Stellarator '+stells)
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
	#bashCommand = wolframScript+" -noprompt -script solveSigma.m "+str(etab)+" "+str(b0)+" "+equilibrium+" "+stells+" "+GSpath+geom2math+" "+sigmaFile+" "+GSpath+gs2gridNA+" "+GSpath+gxgridNA
	bashCommand = sigmaScript+" "+str(etab)+" "+str(b0)+" "+equilibrium+" "+stells+" "+GSpath+geom2math+" "+sigmaFile+" "+GSpath+gs2gridNA+" "+GSpath+gxgridNA
	if not path.exists(GSpath+gs2gridNA):
		print(' Running Near-Axis Grid')
		output = subprocess.call(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	#### Iterate over tprim and eta ####
	kysX  = np.zeros(shape=(len(tprimvec),len(fprimvec)))
	kysNA = np.zeros(shape=(len(tprimvec),len(fprimvec)))
	growthRatesX  = np.zeros(shape=(len(tprimvec),len(fprimvec)))
	growthRatesNA = np.zeros(shape=(len(tprimvec),len(fprimvec)))
	realFrequenciesX  = np.zeros(shape=(len(tprimvec),len(fprimvec)))
	realFrequenciesNA = np.zeros(shape=(len(tprimvec),len(fprimvec)))
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
			if np.max(growthRateX)<=0:
				kysX[countx,county ]= np.inf
				growthRatesX[countx,county] = 0
				realFrequenciesX[countx,county ]= np.inf
			else:
				kysX[countx,county ]=kyX[growthRateX.index(max(growthRateX))]
				growthRatesX[countx,county] =np.max(growthRateX)
				realFrequenciesX[countx,county ]=realFrequencyX[growthRateX.index(max(growthRateX))]
			if np.max(growthRateNA)<=0:
				kysNA[countx,county] = np.inf
				growthRatesNA[countx,county] = 0
				realFrequenciesNA[countx,county]= np.inf
			else:
				kysNA[countx,county]=kyNA[growthRateNA.index(max(growthRateNA))]
				growthRatesNA[countx,county]=np.max(growthRateNA)
				realFrequenciesNA[countx,county]=realFrequencyNA[growthRateNA.index(max(growthRateNA))]
			xdata.append(fprim)
			ydata.append(tprim)
	if plotFig==1:
		## Plot maximum growth rate and corresponding ky for each fprim and tprim
		print('Plotting max gamma and ky for each tprim and fprim')
		if not path.exists('Figures_r'+str(rr)):
			os.mkdir('Figures_r'+str(rr))
		plotExtent=[min(fprimvec),max(fprimvec),min(tprimvec),max(tprimvec)]
		finalPlot(growthRatesX,growthRatesNA,plotExtent,'gamma',stells,rr,aky_max)
		finalPlot(realFrequenciesX,realFrequenciesNA,plotExtent,'omega',stells,rr,aky_max)
		finalPlot(kysX,kysNA,plotExtent,'ky',stells,rr,aky_max)

