from tempfile import mkstemp
from os import path, fdopen, remove
from shutil import move, copymode, copyfile, rmtree
import os
import glob
from scipy.io import netcdf
import subprocess
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from matplotlib.font_manager import FontProperties
from pylab import text

plotfontSize=20;figSize1=7.5;figSize2=4.0;legendfontSize=14;annotatefontSize=8
matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
matplotlib.rc('text', usetex=False);matplotlib.rcParams['text.latex.preamble']=r"\usepackage{amsmath}"

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

###### RUN GS2 and diagnostics
def runGS2func(plotFig,currentPath,GSpath,runName,runNameNA,gs2grid,gxgrid,gs2gridNA,gxgridNA,stells,nzgrid,ne,dt,nstep,ngauss,tprim,fprim,aky_min,aky_max,naky,ncores,runGS2):
    # Create VMEC input file
    copyfile("gs2Input.in",GSpath+"gs2"+"/"+runName)
    os.chdir(GSpath+"gs2")
    replace(runName,' negrid = 12',' negrid = '+str(ne))
    replace(runName,' delt = 0.2',' delt = '+str(dt))
    replace(runName,' nstep = 500',' nstep = '+str(int(nstep)))
    replace(runName,' ntheta = 64',' ntheta = '+str(2*nzgrid))
    replace(runName,' ngauss = 3',' ngauss = '+str(ngauss))
    replace(runName,' tprim = 6.9',' tprim = '+str(tprim))
    replace(runName,' fprim = 2.2',' fprim = '+str(fprim))
    replace(runName,' aky_min = 0.0',' aky_min = '+str(aky_min))
    replace(runName,' aky_max = 2.0',' aky_max = '+str(aky_max))
    replace(runName,' naky = 10',' naky = '+str(naky))
    # VMEC run
    replace(runName,' gridout_file = ',' gridout_file = "../'+gs2grid+'"')
    bashCommand = "mpirun -n "+str(ncores)+" ./gs2 "+runName
    if runGS2==1:
        if not path.exists(runName[:-3]+".out.nc"):
            print(' Running GS2 VMEC')
            output = subprocess.call(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # Near-Axis run
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
    if not path.exists(runName[:-3]+".out.nc_eigenphi.pdf") and plotFig==1:
        print("Not eigenplotting")
        #print(' Plotting eigenfunctions')
        #eigenPlot(runName[:-3]+".out.nc")
    if not path.exists(runNameNA[:-3]+".out.nc_eigenphi.pdf") and plotFig==1:
        print("Not eigenplotting")
        #eigenPlot(runNameNA[:-3]+".out.nc")
    # Compare Near-Axis and VMEC's geometric coefficients
    if not path.exists('../'+stells+"_geom.pdf") and plotFig==1:
        print(' Plotting geometric coefficients')
        geomPlot(stells,runName[:-3]+".out.nc",runNameNA[:-3]+".out.nc")
        move(stells+'_geom.pdf','../'+stells+'_geom.pdf')
    ############################
    os.chdir(currentPath)

###### Function to plot eigenfunctions
def eigenPlot(stellFile):
    f  = netcdf.netcdf_file(stellFile,'r',mmap=False)
    y  = f.variables['phi'][()]
    x  = f.variables['theta'][()]
    kyvec = f.variables['ky'][()]
    plt.figure(figsize=(8*figSize1,10*figSize2))
    #plt.subplots_adjust(hspace=1.0)
    for i,v in enumerate(range(len(kyvec))):
        v = v+1
        phiR0= y[i,0,int((len(x)-1)/2+1),0]
        phiI0= y[i,0,int((len(x)-1)/2+1),1]
        phi02= phiR0**2+phiI0**2
        phiR = (y[i,0,:,0]*phiR0+y[i,0,:,1]*phiI0)/phi02
        phiI = (y[i,0,:,1]*phiR0-y[i,0,:,0]*phiI0)/phi02
        ##############
        ax1 = plt.subplot(int(len(kyvec)/2),2,v)
        ax1.plot(x, phiR, label=r'Re($\phi/\phi_0$)')
        ax1.plot(x, phiI, label=r'Im($\phi/\phi_0$)')
        #plt.title(r'$\phi$')
        if i==1:
            fontP = FontProperties()
            #fontP.set_size('xx-small')
            legend=ax1.legend(bbox_to_anchor=(1.01, 1), loc='upper left', prop=fontP)
            #ax1.legend(loc="upper right",fontsize=legendfontSize)
        text(0.9, 0.7,"ky = "+str('%.2f' % kyvec[i]), ha='center', va='center', transform=ax1.transAxes,fontsize=legendfontSize)
        if v==int(len(kyvec))-1 or v==int(len(kyvec)):
            ax1.set_xlabel(r'$\theta$')
    plt.tight_layout()
    #plt.subplots_adjust(left=0.16, bottom=0.19, right=0.9, top=0.93)
    plt.savefig(stellFile+'_eigenphi.pdf', format='pdf', bbox_inches='tight')
    plt.close()
    return 0


##### Function to obtain gamma and omega for each ky
def gammabyky(plotFig,currentPath,GSpath,runName,runNameNA,fractionToConsider):
    os.chdir(GSpath+"gs2")
    # Compute growth rate:
    fX   = netcdf.netcdf_file(runName[:-3]+".out.nc",'r',mmap=False)
    fNA  = netcdf.netcdf_file(runNameNA[:-3]+".out.nc",'r',mmap=False)
    tX   = fX.variables['t'][()]
    tNA  = fNA.variables['t'][()]
    kyX  = fX.variables['ky'][()]
    kyNA = fNA.variables['ky'][()]
    phi2_by_kyX  = fX.variables['phi2_by_ky'][()]
    phi2_by_kyNA = fNA.variables['phi2_by_ky'][()]
    omegaX  = fX.variables['omega'][()]
    omegaNA = fNA.variables['omega'][()]
    fX.close()
    fNA.close()

    startIndexX  = int(len(tX)*(1-fractionToConsider))
    startIndexNA = int(len(tNA)*(1-fractionToConsider))

    growthRateX  = []
    growthRateNA = []
    ## assume that kyX=kyNA
    for i in range(len(kyX)):
        maskX  = np.isfinite(phi2_by_kyX[:,i])
        maskNA = np.isfinite(phi2_by_kyNA[:,i])
        data_xX = tX[maskX]
        data_yX = phi2_by_kyX[maskX,i]
        data_xNA = tNA[maskNA]
        data_yNA = phi2_by_kyNA[maskNA,i]
        fitX  = np.polyfit(data_xX[startIndexX:], np.log(data_yX[startIndexX:]), 1)
        fitNA = np.polyfit(data_xNA[startIndexNA:], np.log(data_yNA[startIndexNA:]), 1)
        thisGrowthRateX  = fitX[0]/2
        thisGrowthRateNA = fitNA[0]/2
        growthRateX.append(thisGrowthRateX)
        growthRateNA.append(thisGrowthRateNA)

    # Compute real frequency:
    realFreqVsTimeX  = []
    realFreqVsTimeNA = []
    realFrequencyX   = []
    realFrequencyNA  = []
    for i in range(len(kyX)):
        realFreqVsTimeX.append(omegaX[:,i,0,0])
        realFreqVsTimeNA.append(omegaNA[:,i,0,0])
        realFrequencyX.append(np.mean(realFreqVsTimeX[i][startIndexX:]))
        realFrequencyNA.append(np.mean(realFreqVsTimeNA[i][startIndexNA:]))


    if not path.exists(runName[:-3]+"_GammaOmegaKy.pdf") and plotFig==1:
        print(' Plotting gamma and omega per ky')

        numRows = 1
        numCols = 2

        plt.subplot(numRows, numCols, 1)
        plt.plot(kyX,growthRateNA,'.-',label=r'Near-Axis')
        plt.plot(kyNA,growthRateX,'.-',label=r'VMEC')
        plt.xlabel(r'$k_y$')
        plt.ylabel(r'$\gamma$')
        plt.xscale('log')
        #plt.rc('font', size=8)
        #plt.rc('axes', labelsize=8)
        #plt.rc('xtick', labelsize=8)
        plt.legend(frameon=False,prop=dict(size='xx-small'),loc=0)

        plt.subplot(numRows, numCols, 2)
        plt.plot(kyX,realFrequencyNA,'.-',label=r'Near-Axis')
        plt.plot(kyNA,realFrequencyX,'.-',label=r'VMEC')
        plt.xlabel(r'$k_y$')
        plt.ylabel(r'$\omega$')
        plt.xscale('log')
        plt.rc('font', size=8)
        plt.rc('axes', labelsize=8)
        plt.rc('xtick', labelsize=8)
        plt.legend(frameon=False,prop=dict(size=12),loc=0)

        plt.tight_layout()
        #plt.subplots_adjust(left=0.14, bottom=0.15, right=0.98, top=0.96)
        plt.savefig(runName[:-3]+"_GammaOmegaKy.pdf", format='pdf')
        plt.close()
        
    ############################
    os.chdir(currentPath)

    return kyX, kyNA, growthRateX, growthRateNA, realFrequencyX, realFrequencyNA

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

def finalPlot(dataX,dataNA,plotExtent,strLabel,stells,rr,aky_max):
    numRows = 1
    numCols = 2
    plt.figure(figsize=(11, 6))
    ax1=plt.subplot(numRows, numCols, 1)
    ax1.set_title("VMEC")
    im=plt.imshow(dataX, interpolation='hermite', origin='lower', extent=plotExtent, cmap='jet')
    clb=plt.colorbar(im,fraction=0.046, pad=0.04)
    #if strLabel=='gamma':
    #    plt.clim(0,1)
    #elif strLabel=='ky':
    #    plt.clim(0,aky_max)
    #else:
    plt.clim(np.min([np.min(dataX),np.min(dataNA)]),np.max([np.max(np.where(np.isinf(dataX),-np.Inf,dataX)),np.max(np.where(np.isinf(dataNA),-np.Inf,dataNA))])) 
    plt.xlabel(r'$a/L_n$', usetex=True)
    plt.ylabel(r'$a/L_T$', usetex=True)

    if strLabel=='gamma':
        clb.ax.set_title(r'$\gamma$', usetex=True)
    elif strLabel=='omega':
        clb.ax.set_title(r'$\omega$', usetex=True)
    elif strLabel=='ky':
        clb.ax.set_title(r'$k_y$', usetex=True)

    ax2=plt.subplot(numRows, numCols, 2)
    ax2.set_title("Near-Axis")
    im=plt.imshow(dataNA, interpolation='hermite', origin='lower', extent=plotExtent, cmap='jet')
    clb=plt.colorbar(im,fraction=0.046, pad=0.04)
    #if strLabel=='gamma':
    #    plt.clim(0,1)
    #elif strLabel=='ky':
    #    plt.clim(0,aky_max)
    #else:
    plt.clim(np.min([np.min(dataX),np.min(dataNA)]),np.max([np.max(np.where(np.isinf(dataX),-np.Inf,dataX)),np.max(np.where(np.isinf(dataNA),-np.Inf,dataNA))])) 
    plt.xlabel(r'$a/L_n$')
    plt.ylabel(r'$a/L_T$')

    plt.tight_layout()

    if strLabel=='gamma':
        clb.ax.set_title(r'$\gamma$', usetex=True)
    elif strLabel=='omega':
        clb.ax.set_title(r'$\omega$', usetex=True)
    elif strLabel=='ky':
        clb.ax.set_title(r'$k_y$', usetex=True)

    matplotlib.rc('font', size=16)
    plt.savefig('Figures_r'+str(rr)+'/'+stells+'_max'+strLabel+'.pdf', format='pdf', bbox_inches='tight')
    plt.close()


