from tempfile import mkstemp
from os import path, fdopen, remove
from shutil import move, copymode, copyfile
from scipy.io import netcdf
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import glob
import os
import subprocess
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties
import pickle

plotfontSize=20;figSize1=7.5;figSize2=4.0;legendfontSize=14;annotatefontSize=8
matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
matplotlib.rc('text', usetex=True);matplotlib.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"


###### RUN GS2 and diagnostics
def runGS2func(stells,rr,nzgrid,nlambda,nfp,equilibrium,runsPath,etab,b0,ne,dt,nstep,ngauss,ncores,wolframScript,gs2Path,runGS2,currentPath):
	# Create folders
	if not path.exists(runsPath+stells):
		os.mkdir(runsPath+stells)
	if not path.exists(runsPath+stells+'/r'+str(rr)):
		os.mkdir(runsPath+stells+'/r'+str(rr))
	if not path.exists(runsPath+stells+'/r'+str(rr)+'/nzgrid'+str(nzgrid)):
		os.mkdir(runsPath+stells+'/r'+str(rr)+'/nzgrid'+str(nzgrid))
	if not path.exists(runsPath+stells+'/r'+str(rr)+'/nzgrid'+str(nzgrid)+'/nlambda'+str(nlambda)):
		os.mkdir(runsPath+stells+'/r'+str(rr)+'/nzgrid'+str(nzgrid)+'/nlambda'+str(nlambda))
	if not path.exists(runsPath+stells+'/r'+str(rr)+'/nzgrid'+str(nzgrid)+'/nlambda'+str(nlambda)+'/nfp'+str(nfp)):
		os.mkdir(runsPath+stells+'/r'+str(rr)+'/nzgrid'+str(nzgrid)+'/nlambda'+str(nlambda)+'/nfp'+str(nfp))
	GSpath = runsPath+stells+'/r'+str(rr)+'/nzgrid'+str(nzgrid)+'/nlambda'+str(nlambda)+'/nfp'+str(nfp)+'/'
	#####################
	## Run VMEC_to_GS2 ##
	#####################
	gs2grid      = 'gs2_grid'+stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	gxgrid       = 'gx_grid'+stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	geom2math    = 'gridMath'+stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	f = netcdf.netcdf_file(equilibrium,'r',mmap=False)
	iotaVMEC = f.variables['iotaf'][()][1]
	nfpVMEC  = f.variables['nfp'][()]
	if not path.exists(GSpath+gs2grid):
		print(' Running VMEC_2_GS2')
		process = subprocess.call(['./vmec2gs2',equilibrium,GSpath+gs2grid,GSpath+geom2math,str(nfpVMEC*nfp/abs(iotaVMEC)),str(rr),GSpath+gxgrid,str(nzgrid),str(nlambda)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	#print('./vmec2gs2'+' '+equilibrium+' '+GSpath+gs2grid+' '+GSpath+geom2math+' '+str(nfpVMEC*nfp/abs(iotaVMEC))+' '+str(rr)+' '+GSpath+gxgrid+' '+str(nzgrid)+' '+str(nlambda))
	#########################
	## Get Near-Axis Grids ##
	#########################
	gs2gridNA   = 'gs2_NAgrid'+stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	gxgridNA    = 'gx_NAgrid' +stells+'r'+str(rr)+'nzgrid'+str(nzgrid)+'nlambda'+str(nlambda)+'nfp'+str(nfp)+'.out'
	sigmaFile   = runsPath+stells+'/Math'+stells+'.mx'
	bashCommand = wolframScript+" -noprompt -script solveSigma.m "+str(etab)+" "+str(b0)+" "+equilibrium+" "+stells+" "+GSpath+geom2math+" "+sigmaFile+" "+GSpath+gs2gridNA+" "+GSpath+gxgridNA
	if not path.exists(GSpath+gs2gridNA):
		print(' Running Mathematica')
		output = subprocess.call(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	######################################
	## Run GS2 to check for convergence ##
	######################################
	if not path.exists(GSpath+'gs2'):
		os.mkdir(GSpath+'gs2')
	copyfile(gs2Path+"/gs2",GSpath+"gs2/gs2")
	#print(gs2Path+"/gs2")
	#print(GSpath+"gs2/gs2")
	copymode(gs2Path+"/gs2",GSpath+"gs2/gs2")
	runName=stells+'ne'+str(ne)+'dt'+str(dt)+'nt'+str(nstep)+'ngauss'+str(ngauss)+'.in'
	copyfile("gs2Input.in",GSpath+"gs2"+"/"+runName)
	os.chdir(GSpath+"gs2")
	replace(runName,' negrid = 12',' negrid = '+str(ne))
	replace(runName,' delt = 0.2',' delt = '+str(dt))
	replace(runName,' nstep = 500',' nstep = '+str(int(nstep)))
	replace(runName,' ntheta = 64',' ntheta = '+str(2*nzgrid))
	replace(runName,' ngauss = 3',' ngauss = '+str(ngauss))
	# VMEC run
	replace(runName,' gridout_file = ',' gridout_file = "../'+gs2grid+'"')
	bashCommand = "mpirun -n "+str(ncores)+" ./gs2 "+runName
	if runGS2==1:
		if not path.exists(runName[:-3]+".out.nc"):
			print(' Running GS2 VMEC')
			output = subprocess.call(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Near-Axis run
	runNameNA=runName[:-3]+'NA.in'
	copyfile(runName,runNameNA)
	replace(runNameNA,' gridout_file = "../'+gs2grid+'"',' gridout_file = "../'+gs2gridNA+'"')
	bashCommand = "mpirun -n "+str(ncores)+" ./gs2 "+runNameNA
	if runGS2==1:
		if not path.exists(runNameNA[:-3]+".out.nc"):
			print(' Running GS2 Near-Axis')
			output = subprocess.call(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	# Remove spurious files
	removeGS2(os.getcwd())
	# Plot eigenfunctions
	if not path.exists(runName[:-3]+".out.nc_eigenphi.pdf"):
		print(' Plotting eigenfunctions')
		eigenPlot(runName[:-3]+".out.nc")
	if not path.exists(runNameNA[:-3]+".out.nc_eigenphi.pdf"):
		eigenPlot(runNameNA[:-3]+".out.nc")
	# Compare Near-Axis and VMEC's geometric coefficients
	if not path.exists('../'+stells+"_geom.pdf"):
		print(' Plotting geometric coefficients')
		geomPlot(stells,runName[:-3]+".out.nc",runNameNA[:-3]+".out.nc")
		move(stells+'_geom.pdf','../'+stells+'_geom.pdf')
	############################
	os.chdir(currentPath)
    ############################
	print(' Computing growth rate and frequency')
	return getgamma(GSpath+'gs2/'+runName[:-3]+".out.nc"), getgamma(GSpath+'gs2/'+runNameNA[:-3]+".out.nc")


## Function to replace text in files
def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    #Copy the file permissions from the old file to the new file
    copymode(file_path, abs_path)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)
    return 0

###### Function to remove spurious GS2 files
def removeGS2(Path):
    currDir=os.getcwd()
    os.chdir(Path)
    for f in glob.glob('*.amoments'): remove(f)
    for f in glob.glob('*.eigenfunc'): remove(f)
    for f in glob.glob('*.error'): remove(f)
    for f in glob.glob('*.fields'): remove(f)
    for f in glob.glob('*.g'): remove(f)
    for f in glob.glob('*.lpc'): remove(f)
    for f in glob.glob('*.mom2'): remove(f)
    for f in glob.glob('*.moments'): remove(f)
    for f in glob.glob('*.vres'): remove(f)
    for f in glob.glob('*.vres2'): remove(f)
    for f in glob.glob('*.exit_reason'): remove(f)
    for f in glob.glob('*.optim'): remove(f)
    for f in glob.glob('*.out'): remove(f)
    remove("gs2")
    os.chdir(currDir)
    return 0

###### Function to obtain growth rate and plot phi2
def getgamma(stellFile):
	fractionToConsider=0.35
	f = netcdf.netcdf_file(stellFile,'r',mmap=False)
	phi2 = np.log(f.variables['phi2'][()])
	t = f.variables['t'][()]
	startIndex = int(len(t)*(1-fractionToConsider))
	mask = np.isfinite(phi2)
	data_x = t[mask]
	data_y = phi2[mask]
	fit = np.polyfit(data_x[startIndex:], data_y[startIndex:], 1)
	poly = np.poly1d(fit)
	GrowthRate = fit[0]/2
	gamma  = np.mean(f.variables['omega'][()][startIndex:,0,0,1])
	omega  = np.mean(f.variables['omega'][()][startIndex:,0,0,0])
    #fitRes = np.poly1d(coeffs)
	if not path.exists(stellFile+'_phi2.pdf'):
		plt.figure(figsize=(figSize1,figSize2))
		##############
		plt.plot(t, phi2,'.', label=r'data - $\gamma_{GS2} = $'+str(gamma))
		plt.plot(t, poly(t),'-', label=r'fit - $\gamma = $'+str(GrowthRate))
		##############
		plt.legend(loc=0,fontsize=legendfontSize)
		plt.xlabel(r'$t$');plt.ylabel(r'$\ln |\phi|^2$')
		plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
		plt.savefig(stellFile+'_phi2.pdf', format='pdf')
		plt.close()
	return (GrowthRate,abs(omega))

###### Function to plot eigenfunctions
def eigenPlot(stellFile):
    f = netcdf.netcdf_file(stellFile,'r',mmap=False)
    y = f.variables['phi'][()]
    x = f.variables['theta'][()]
    plt.figure(figsize=(figSize1,figSize2))
    phiR0= y[0,0,int((len(x)-1)/2+1),0]
    phiI0= y[0,0,int((len(x)-1)/2+1),1]
    phi02= phiR0**2+phiI0**2
    phiR = (y[0,0,:,0]*phiR0+y[0,0,:,1]*phiI0)/phi02
    phiI = (y[0,0,:,1]*phiR0-y[0,0,:,0]*phiI0)/phi02
    ##############
    plt.plot(x, phiR, label=r'Re($\phi/\phi_0$)')
    plt.plot(x, phiI, label=r'Im($\phi/\phi_0$)')
    ##############
    plt.xlabel(r'$\theta$');plt.ylabel(r'$\phi$')
    plt.legend(loc="upper right")
    plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.93)
    plt.savefig(stellFile+'_eigenphi.pdf', format='pdf')
    plt.close()
    return 0

###### Function to plot geometry coefficients
def geomPlot(stells,stellFileX,stellFileNA):
    fX  = netcdf.netcdf_file(stellFileX,'r',mmap=False)
    fNA = netcdf.netcdf_file(stellFileNA,'r',mmap=False)
    theta      = fX.variables['theta'][()]
    matplotlib.rc('font', size=6)
    nrows=5; ncols=2;fig = plt.figure()
    ##
    plt.subplot(nrows, ncols, 1)
    plt.scatter(theta, fX.variables['gradpar'][()] , color='b', label='X', s=0.1)
    plt.scatter(theta, fNA.variables['gradpar'][()], color='r', label='NA', s=0.1)
    plt.xlabel(r'$\theta$');plt.ylabel(r'gradpar')
    ##
    plt.subplot(nrows, ncols, 2);
    plt.scatter(theta, fX.variables['bmag'][()] , color='b', label='X', s=0.1)
    plt.scatter(theta, fNA.variables['bmag'][()], color='r', label='NA', s=0.1)
    plt.xlabel(r'$\theta$');plt.ylabel(r'bmag')
    ##
    plt.subplot(nrows, ncols, 3);
    plt.scatter(theta, fX.variables['gds2'][()] , color='b', label='X', s=0.1)
    plt.scatter(theta, fNA.variables['gds2'][()], color='r', label='NA', s=0.1)
    plt.xlabel(r'$\theta$');plt.ylabel(r'gds2')
    ##
    plt.subplot(nrows, ncols, 4);
    plt.scatter(theta, fX.variables['gds21'][()] , color='b', label='X', s=0.1)
    plt.scatter(theta, fNA.variables['gds21'][()], color='r', label='NA', s=0.1)
    plt.xlabel(r'$\theta$');plt.ylabel(r'gds21')
    ##
    plt.subplot(nrows, ncols, 5);
    plt.scatter(theta, fX.variables['gds22'][()] , color='b', label='X', s=0.1)
    plt.scatter(theta, fNA.variables['gds22'][()], color='r', label='NA', s=0.1)
    plt.xlabel(r'$\theta$');plt.ylabel(r'gds22')
    ##
    plt.subplot(nrows, ncols, 6);
    plt.scatter(list(range(1, 1+len(fX.variables['lambda'][()]))),fX.variables['lambda'][()] , color='b', label='X', s=0.2)
    plt.scatter(list(range(1, 1+len(fNA.variables['lambda'][()]))),fNA.variables['lambda'][()], color='r', label='NA', s=0.2)
    plt.xlabel(r'');plt.ylabel(r'lambda')
    ##
    plt.subplot(nrows, ncols, 7);
    plt.scatter(theta, fX.variables['gbdrift'][()] , color='b', label='X', s=0.1)
    plt.scatter(theta, fNA.variables['gbdrift'][()], color='r', label='NA', s=0.1)
    plt.xlabel(r'$\theta$');plt.ylabel(r'gbdrift')
    ##
    plt.subplot(nrows, ncols, 8);
    plt.scatter(theta, fX.variables['gbdrift0'][()] , color='b', label='X', s=0.1)
    plt.scatter(theta, fNA.variables['gbdrift0'][()], color='r', label='NA', s=0.1)
    plt.xlabel(r'$\theta$');plt.ylabel(r'gbdrift0')
    ##
    plt.subplot(nrows, ncols, 9);
    plt.scatter(theta, fX.variables['cvdrift'][()] , color='b', label='X', s=0.1)
    plt.scatter(theta, fNA.variables['cvdrift'][()], color='r', label='NA', s=0.1)
    plt.xlabel(r'$\theta$');plt.ylabel(r'cvdrift')
    ##
    plt.subplot(nrows, ncols, 10);
    l1=plt.scatter(theta, fX.variables['cvdrift0'][()] , color='b', label='X', s=0.1)
    l2=plt.scatter(theta, fNA.variables['cvdrift0'][()], color='r', label='NA', s=0.1)
    plt.xlabel(r'$\theta$');plt.ylabel(r'cvdrift0')
    ##
    plt.subplots_adjust(left=0.08, bottom=0.08, right=0.98, top=0.97, wspace=0.27, hspace=0.3)
    fig.legend([l1,l2], ['X', 'NA'], loc = 'lower center', ncol=2)
    plt.savefig(stells+'_geom.pdf', format='pdf')
    plt.close()
    return 0

###### Function to plot convergence on growth rates and frequencies for each stellarator
def gammaPlot(growthRateX,growthRateNA,strLabel,legendTXT,stells,rr,ngauss,dt,nstep,ne,nzgrid,nlambda,nfp):
	plotfontSize=14;figSize1=5.7;figSize2=4.2;titleFontSize=16
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
	fig, ax = plt.subplots(figsize=(figSize1,figSize2))
	colors = cm.rainbow(np.linspace(0, 1, len(legendTXT)))
	for x, y, c, l in zip(growthRateNA, growthRateX, colors, legendTXT):
		plt.scatter(x, y, color=c, label=l)
	fontP = FontProperties()
	fontP.set_size('xx-small')
	legend=ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', prop=fontP)
	fig.suptitle(stells, fontsize=titleFontSize)
	ax.plot([0, 1], [0, 1], color='r', ls='--')
	plt.ylim(0.0,1.01*max(max(growthRateX),max(growthRateNA)))
	plt.xlim(0.0,1.01*max(max(growthRateX),max(growthRateNA)))
	plt.gca().set_aspect('equal', adjustable='box')
	if strLabel=='gamma':
		plt.xlabel(r'Near-Axis $\gamma$');plt.ylabel(r'VMEC $\gamma$')
	elif strLabel=='omega':
		plt.xlabel(r'Near-Axis $\omega$');plt.ylabel(r'VMEC $\omega$')
	plt.subplots_adjust(left=0.11, bottom=0.12, right=0.79, top=0.92)
	## Add text below legend with base case parameters
	offset = matplotlib.text.OffsetFrom(legend, (1.0, 0.0))
	infoString=r"\textbf{Base case}"+"\nngauss = "+str(ngauss)+"\ndeltat = "+str(dt)+"\nnstep = "+str(nstep)+"\nnegrid = "+str(ne)+"\nnzgrid = "+str(nzgrid)+"\nnlambda = "+str(nlambda)+"\nnfp = "+str(nfp)
	ax.annotate(infoString, xy=(0,0),size=10,
            xycoords='figure fraction', xytext=(-80,-20), textcoords=offset, 
            horizontalalignment='left', verticalalignment='top')
	fig.canvas.draw()
	## Save figure
	plt.savefig(strLabel+'Plots/'+strLabel+stells+'_r'+str(rr)+'.pdf', format='pdf')
	plt.close()
	return 0

###### Function to plot final converged growth rates and frequencies for all stellarators (base case)
def finalGammaPlot(gammaNA,gammaX,strLabel,stellsToRun,stellDesigns,rr):
	plotfontSize=16;annotatefontSize=8;figSize1=5.8;figSize2=5.2
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize);
	fig, ax = plt.subplots(figsize=(figSize1,figSize2))
	ax.scatter(gammaNA,gammaX)
	for count, i in enumerate(stellsToRun):
		ax.annotate(stellDesigns[i], (gammaNA[count],gammaX[count]), fontsize=annotatefontSize)
	ax.plot([0, 1], [0, 1], color='r', ls='--')
	if strLabel=='gamma':
		plt.xlabel(r'Near-Axis $\gamma$');plt.ylabel(r'VMEC $\gamma$')
		plt.ylim(0.95*min(min(gammaNA),min(gammaX)),1.06*max(max(gammaNA),max(gammaX)))
		plt.xlim(0.95*min(min(gammaNA),min(gammaX)),1.06*max(max(gammaNA),max(gammaX)))
	elif strLabel=='omega':
		plt.xlabel(r'Near-Axis $\omega$');plt.ylabel(r'VMEC $\omega$')
		plt.ylim(-0.005,1.1*max(max(gammaNA),max(gammaX)))
		plt.xlim(-0.005,1.1*max(max(gammaNA),max(gammaX)))
	plt.gca().set_aspect('equal', adjustable='box')
	plt.subplots_adjust(left=0.12, bottom=0.12, right=1.0, top=0.97)
	plt.savefig(strLabel+'Stells_r'+str(rr)+'.pdf', format='pdf')
	plt.close()

###### Function to plot final converged growth rates and frequencies for all stellarators (all cases)
def allGammaPlot(gammaNA,gammaX,strLabel,stellsToRun,stellDesigns,rr):
	with open(strLabel+'.pickle','rb') as fid:
		ax=pickle.load(fid)
	for count, i in enumerate(stellsToRun):
		ax.annotate(stellDesigns[i], (gammaNA[count],gammaX[count]), fontsize=annotatefontSize)
	ax.plot([0, 1], [0, 1], color='r', ls='--')
	if strLabel=='gamma':
		plt.xlabel(r'Near-Axis $\gamma$');plt.ylabel(r'VMEC $\gamma$')
		plt.ylim(0.95*min(min(gammaNA),min(gammaX)),1.06*max(max(gammaNA),max(gammaX)))
		plt.xlim(0.95*min(min(gammaNA),min(gammaX)),1.06*max(max(gammaNA),max(gammaX)))
	elif strLabel=='omega':
		plt.xlabel(r'Near-Axis $\omega$');plt.ylabel(r'VMEC $\omega$')
		plt.ylim(-0.005,1.1*max(max(gammaNA),max(gammaX)))
		plt.xlim(-0.005,1.1*max(max(gammaNA),max(gammaX)))
	plt.gca().set_aspect('equal', adjustable='box')
	fontP = FontProperties()
	fontP.set_size('xx-small')
	ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', prop=fontP)
	#fig.suptitle(stells, fontsize=titleFontSize)
	plt.subplots_adjust(left=-0.07, bottom=0.16, right=1.0, top=0.96)
	plt.savefig('all'+strLabel+'Stells_r'+str(rr)+'.pdf', format='pdf')
	plt.close()