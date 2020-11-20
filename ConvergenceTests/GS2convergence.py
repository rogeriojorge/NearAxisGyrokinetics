#!/usr/bin/env python3
################################################
################################################

vmecFiles    = ['wistella_midscale' ,'NuhrenbergZille_1988_QHS','HSX_QHS_vacuum_ns201','n4qh.b4.a79a','Drevlak_qh_8_7','li383_1.4m_ns201','n3are_R7.75B5.7_hires','GarabedianQAS2_noCurrentOnAxis_ns201','estell_24_scaled','cfqs_freeBoundary_vacuum_hiRes','st_a34_i32v22_beta_35_scaledAUG_hires','LandremanSengupta2019_section5.4','ITER']
stellDesigns = ['WISTELL-A'         ,'NZ1988'                  ,'HSX'                 ,'KuQHS48'     ,'Drevlak'       ,'NCSX'            ,'ARIES-CS'             ,'QAS2'                                ,'ESTELL'          ,'CFQS'                          ,'Henneberg'                            ,'NAQS'                             ,'ITER']
etabar       = [0.791               ,0.157                     ,1.28                  ,0.147         ,0.0899          ,0.408             ,0.0740                 ,0.347                                 ,0.570             ,0.586                           ,0.302                                  ,1.549                              ,0.115]
B0           = [2.54                ,0.205                     ,1.00                  ,1.20          ,3.97            ,1.55              ,5.69                   ,1.79                                  ,1.00              ,0.933                           ,2.41                                   ,1.001                              ,5.358]
nzgridvec    = [300                 ,300                       ,300                   ,300           ,300             ,300               ,300                    ,350                                   ,300               ,300                             ,300                                    ,300                                ,80]
nlambdavec   = [40                  ,40                        ,40                    ,40            ,40              ,40                ,40                     ,45                                    ,40                ,40                              ,40                                     ,40                                 ,20]
NFPtoInclude = [4                   ,3                         ,3                     ,3             ,3               ,2                 ,2                      ,2                                     ,3                 ,2                               ,2                                      ,2                                  ,3]
negridvec    = [12                  ,12                        ,12                    ,12            ,12              ,12                ,12                     ,12                                    ,12                ,12                              ,12                                     ,12                                 ,12]
deltaTvec    = [0.4                 ,0.4                       ,0.4                   ,0.4           ,0.4             ,0.4               ,0.4                    ,0.3                                   ,0.4               ,0.4                             ,0.4                                    ,0.4                                ,0.3]
simTimevec   = [200                 ,200                       ,200                   ,200           ,200             ,200               ,200                    ,200                                   ,200               ,200                             ,200                                    ,200                                ,200]
ngaussvec    = [3                   ,3                         ,3                     ,3             ,3               ,3                 ,3                      ,3                                     ,3                 ,3                               ,3                                      ,3                                  ,3]
QAorQH       = ['QH'                ,'QH'                      ,'QH'                  ,'QH'          ,'QH'            ,'QA'              ,'QA'                   ,'QA'                                  ,'QA'              'QA'                             ,'QA'                                   ,'QA'                               ,'A']
## 11 stellarators + 1 tokamak

rr  = 0.3 # normalized flux psi/psi_a radial coordinate
stellsToRun  = [0,1,2,3,4,5,6,7,8,9,10] # select the indices of stellarators from stellDesigns to analyze. Start at 0
ncores = 8 # number of cpu cores to run GS2
runGS2 = 1 # 0 - only plot, no GS2, 1 - run GS2
makeCleanVMEC2GS2 = 1 # 1 - make clean and then make VMEC2GS2, 0 - just make

vmecGS2interfaceFolder='/Users/rogeriojorge/Dropbox/PostDoc/NearAxisGyrokinetics/ConvergenceTests/VMEC_to_GS2'
runsPath='/Users/rogeriojorge/scratch/GS2runs/'
equilibriaFolder='/Users/rogeriojorge/Dropbox/PostDoc/NearAxisGyrokinetics/equilibria/'
gs2Path='/Users/rogeriojorge/local/gs2/bin'
wolframScript = "wolframscript"

import os
from os import path, fdopen, remove
import subprocess
from shutil import move, copymode, copyfile, rmtree
from scipy.io import netcdf
from tempfile import mkstemp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties
import pickle

## User-defined functions from convergenceFuncs.py
from convergenceFuncs import runGS2func, replace, getgamma, eigenPlot, geomPlot, removeGS2, gammaPlot, finalGammaPlot, allGammaPlot

##############################

currentPath = os.getcwd()

print("Compiling VMEC to GS2")
os.chdir(vmecGS2interfaceFolder)
if makeCleanVMEC2GS2==1:
	process = subprocess.call(['make','clean'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
process = subprocess.call(['make'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
copyfile("test_vmec_to_gs2_geometry_interface",currentPath+"/vmec2gs2")
copymode("test_vmec_to_gs2_geometry_interface",currentPath+"/vmec2gs2")
os.chdir(currentPath)

##############################
## Variables to iterate: nzgrid, nlambda, nfp (field line length), negrid, deltat, nstep, ngauss
##############################
if not path.exists(runsPath):
	os.mkdir(runsPath)
if not path.exists('gammaPlots'):
#	rmtree('gammaPlots')
	os.mkdir('gammaPlots')
if not path.exists('omegaPlots'):
#	rmtree('omegaPlots')
	os.mkdir('omegaPlots')

gammaX  = []
gammaNA = []
omegaX  = []
omegaNA = []
fig, ax = plt.subplots()
with open('gamma.pickle','wb') as fid:
	pickle.dump(ax, fid)
fig, ax = plt.subplots()
with open('omega.pickle','wb') as fid:
	pickle.dump(ax, fid)
legendTXT=['base case','double ngauss','half dt','double nstep','double negrid','double nzgrid','double nlambda','double nfp']
for count, i in enumerate(stellsToRun):
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
	# Run GS2 and double resolution in each parameter
	print(stells+' base case')
	[gammaXbc, gammaNAbc] = runGS2func(stells,rr,nzgrid,nlambda,nfp,equilibrium,runsPath,etab,b0,ne,dt,nstep,ngauss,ncores,wolframScript,gs2Path,runGS2,currentPath)
	print(stells+' double ngauss')
	[gammaXngauss, gammaNAngauss] = runGS2func(stells,rr,nzgrid,nlambda,nfp,equilibrium,runsPath,etab,b0,ne,dt,nstep,2*ngauss,ncores,wolframScript,gs2Path,runGS2,currentPath)
	print(stells+' half dt')
	[gammaXdt, gammaNAdt] = runGS2func(stells,rr,nzgrid,nlambda,nfp,equilibrium,runsPath,etab,b0,ne,dt/2,2*nstep,ngauss,ncores,wolframScript,gs2Path,runGS2,currentPath)
	print(stells+' double nstep')
	[gammaXnstep, gammaNAnstep] = runGS2func(stells,rr,nzgrid,nlambda,nfp,equilibrium,runsPath,etab,b0,ne,dt,2*nstep,ngauss,ncores,wolframScript,gs2Path,runGS2,currentPath)
	print(stells+' double negrid')
	[gammaXnegrid, gammaNAnegrid] = runGS2func(stells,rr,nzgrid,nlambda,nfp,equilibrium,runsPath,etab,b0,2*ne,dt,nstep,ngauss,ncores,wolframScript,gs2Path,runGS2,currentPath)
	print(stells+' double nzgrid')
	[gammaXnzgrid, gammaNAnzgrid] = runGS2func(stells,rr,2*nzgrid,nlambda,nfp,equilibrium,runsPath,etab,b0,ne,dt,nstep,ngauss,ncores,wolframScript,gs2Path,runGS2,currentPath)
	print(stells+' double nlambda')
	[gammaXnlambda, gammaNAnlambda] = runGS2func(stells,rr,nzgrid,2*nlambda,nfp,equilibrium,runsPath,etab,b0,ne,dt,nstep,ngauss,ncores,wolframScript,gs2Path,runGS2,currentPath)
	print(stells+' double nfp')
	[gammaXnfp, gammaNAnfp] = runGS2func(stells,rr,2*nzgrid,nlambda,2*nfp,equilibrium,runsPath,etab,b0,ne,dt,nstep,ngauss,ncores,wolframScript,gs2Path,runGS2,currentPath)

	growthRateX  = [gammaXbc[0], gammaXngauss[0], gammaXdt[0], gammaXnstep[0], gammaXnegrid[0], gammaXnzgrid[0], gammaXnlambda[0], gammaXnfp[0]]
	growthRateNA = [gammaNAbc[0],gammaNAngauss[0],gammaNAdt[0],gammaNAnstep[0],gammaNAnegrid[0],gammaNAnzgrid[0],gammaNAnlambda[0],gammaNAnfp[0]]
	frequencyX   = [gammaXbc[1], gammaXngauss[1], gammaXdt[1], gammaXnstep[1], gammaXnegrid[1], gammaXnzgrid[1], gammaXnlambda[1], gammaXnfp[1]]
	frequencyNA  = [gammaNAbc[1],gammaNAngauss[1],gammaNAdt[1],gammaNAnstep[1],gammaNAnegrid[1],gammaNAnzgrid[1],gammaNAnlambda[1],gammaNAnfp[1]]

	gammaX.append(gammaXbc[0])
	gammaNA.append(gammaNAbc[0])
	omegaX.append(gammaXbc[1])
	omegaNA.append(gammaNAbc[1])

	## Plot growth rates
	if not path.exists('gamma'+'Plots/'+'gamma'+stells+'_r'+str(rr)+'.pdf'): gammaPlot(growthRateX,growthRateNA,'gamma',legendTXT,stells,rr,ngauss,dt,nstep,ne,nzgrid,nlambda,nfp)
	if not path.exists('omega'+'Plots/'+'omega'+stells+'_r'+str(rr)+'.pdf'): gammaPlot(frequencyX,frequencyNA,'omega',legendTXT,stells,rr,ngauss,dt,nstep,ne,nzgrid,nlambda,nfp)

	## Add points to figure with all stellarators and cases
	#gamma
	with open('gamma.pickle','rb') as fid:
		ax=pickle.load(fid)
	colors = cm.rainbow(np.linspace(0, 1, len(legendTXT)))
	for x, y, c, l in zip(growthRateNA, growthRateX, colors, legendTXT):
		if count==0:
			plt.scatter(x, y, color=c, label=l)
		else:
			plt.scatter(x, y, color=c)
	with open('gamma.pickle','wb') as fid:
		pickle.dump(ax, fid)
	plt.close()
	#omega
	with open('omega.pickle','rb') as fid:
		ax=pickle.load(fid)
	colors = cm.rainbow(np.linspace(0, 1, len(legendTXT)))
	for x, y, c, l in zip(frequencyNA, frequencyX, colors, legendTXT):
		if count==0:
			plt.scatter(x, y, color=c, label=l)
		else:
			plt.scatter(x, y, color=c)
	with open('omega.pickle','wb') as fid:
		pickle.dump(ax, fid)
	plt.close()

## Plot all stellarators together
print('Plotting final gamma and omega (X vs NA)')
#base case
finalGammaPlot(gammaNA,gammaX,'gamma',stellsToRun,stellDesigns,rr)
finalGammaPlot(omegaNA,omegaX,'omega',stellsToRun,stellDesigns,rr)
#all cases
allGammaPlot(gammaNA,gammaX,'gamma',stellsToRun,stellDesigns,rr)
allGammaPlot(omegaNA,omegaX,'omega',stellsToRun,stellDesigns,rr)
#move figures to folders
if not path.exists('Figures'):
	os.mkdir('Figures')
for f in glob.glob('*.pdf'): copyfile(f,'Figures/'+f)
for f in glob.glob('*.pdf'): remove(f)
for f in glob.glob('*.pickle'): remove(f)