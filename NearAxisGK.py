#!/usr/bin/env python3
################################################
################################################

number_of_field_periods_to_include = 8
normalizedfluxvec = [0.01]#,0.05,0.1,0.3]
N_phi = 200 #Number of points for B0, B1 and B2 calculation and plotting
plotSave = 0 #Spend time plotting results or not

stellDesigns=['WISTELL-A'         ,'NZ1988'                  ,'HSX'                 ,'KuQHS48'     ,'Drevlak'       ,'NCSX'            ,'ARIES-CS'             ,'QAS2'                                ,'ESTELL'          ,'CFQS'                          ,'Henneberg'   ]
etab        =[0.01                ,0.155                     ,1.33                  ,0.147         ,0.0861          ,0.403             ,0.0740                 ,0.341                                 ,0.563             ,0.569                           ,0.269         ]
Nrotations  =[4                   ,-6                        ,4                     ,4             ,-5              ,0                 ,0                      ,0                                     ,0                 ,0                               ,0             ]
vmecFiles   =['wistella_midscale' ,'NuhrenbergZille_1988_QHS','HSX_QHS_vacuum_ns201','n4qh.b4.a79a','Drevlak_qh_8_7','li383_1.4m_ns201','n3are_R7.75B5.7_hires','GarabedianQAS2_noCurrentOnAxis_ns201','estell_24_scaled','cfqs_freeBoundary_vacuum_hiRes','st_a34_i32v22_beta_35_scaledAUG_hires']

equilibriaFolder='equilibria/'
gs2gridsFolder='gs2grids/'
papergridsFolder='paperGrids/'
vmecGS2interfaceFolder='VMEC_to_GS2'
figuresFolder='Figures/'
MathDataFolder='data/'
toPaperFolder='toPaper/'

###### Start NearAxisGK
print("###### Near-Axis Gyrokinetics Interface ######")
print("Number of configurations = "+str(len(stellDesigns)))
import os
from os import path
import subprocess
from shutil import move, copymode, copyfile
import numpy as np
from scipy.io import netcdf
import sys
import warnings
import matplotlib.pyplot as plt
import matplotlib

if not sys.warnoptions:
    warnings.simplefilter("ignore")

###### Function to obtain B0 for each stellarator
def obtainB0B1B2(boozFile,nNormal):
	max_s_for_fit = 0.5
	f = netcdf.netcdf_file(boozFile,'r',mmap=False)
	bmnc = f.variables['bmnc_b'][()]
	ixm = f.variables['ixm_b'][()]
	ixn = f.variables['ixn_b'][()]
	jlist = f.variables['jlist'][()]
	ns = f.variables['ns_b'][()]
	nfp = f.variables['nfp_b'][()]
	Psi = f.variables['phi_b'][()]
	Psi_a = Psi[-1]
	iotaVMECt=f.variables['iota_b'][()][1]
	f.close()
	s_full = np.linspace(0,1,ns)
	ds = s_full[1] - s_full[0]
	#s_half = s_full[1:] - 0.5*ds
	s_half = s_full[jlist-1] - 0.5*ds
	mask = s_half < max_s_for_fit
	s_fine = np.linspace(0,1,400)
	sqrts_fine = s_fine
	phi = np.linspace(0,2*np.pi / nfp, N_phi)
	B0  = np.zeros(N_phi)
	B1s = np.zeros(N_phi)
	B1c = np.zeros(N_phi)
	B20 = np.zeros(N_phi)
	B2s = np.zeros(N_phi)
	B2c = np.zeros(N_phi)
	for jmn in range(len(ixm)):
		m = ixm[jmn]
		n = ixn[jmn] / nfp
		if m>2:
			continue
		if m==0:
        	# For m=0, fit a polynomial in s (not sqrt(s)) that does not need to go through the origin.
			degree = 3
			p = np.polyfit(s_half[mask], bmnc[mask,jmn], degree)
			B0 += p[-1] * np.cos(n*nfp*phi)
			B20 += p[-2] * np.cos(n*nfp*phi)
		if m==1:
        	# For m=1, fit a polynomial in sqrt(s) to an odd function
			x1 = np.sqrt(s_half[mask])
			y1 = bmnc[mask,jmn]
			x2 = np.concatenate((-x1,x1))
			y2 = np.concatenate((-y1,y1))
			degree = 5
			p = np.polyfit(x2,y2, degree)
			B1c += p[-2] * (np.sin(n*nfp*phi) * np.sin(nNormal*phi) + np.cos(n*nfp*phi) * np.cos(nNormal*phi))
			B1s += p[-2] * (np.sin(n*nfp*phi) * np.cos(nNormal*phi) - np.cos(n*nfp*phi) * np.sin(nNormal*phi))
			#B1c += p[-2] * np.cos(n*nfp*phi)
			#B1s += p[-2] * np.sin(n*nfp*phi)
		if m==2:
    		# For m=2, fit a polynomial in s (not sqrt(s)) that does need to go through the origin.
			x1 = s_half[mask]
			y1 = bmnc[mask,jmn]
			x2=x1
			y2=y1
			degree = 4
			p = np.polyfit(x2,y2, degree)
			B2c += p[-2] * (np.sin(n*nfp*phi) * np.sin(nNormal*phi) + np.cos(n*nfp*phi) * np.cos(nNormal*phi))
			B2s += p[-2] * (np.sin(n*nfp*phi) * np.cos(nNormal*phi) - np.cos(n*nfp*phi) * np.sin(nNormal*phi))
			#B2c += p[-2] * np.cos(n*nfp*phi)
			#B2s += p[-2] * np.sin(n*nfp*phi)
	# Convert expansion in sqrt(s) to an expansion in r
	BBar = np.mean(B0)
	sqrt_s_over_r = np.sqrt(np.pi * BBar / Psi_a)
	B1s *= sqrt_s_over_r
	B1c *= sqrt_s_over_r
	B20 *= sqrt_s_over_r*sqrt_s_over_r
	B2c *= sqrt_s_over_r*sqrt_s_over_r
	B2s *= sqrt_s_over_r*sqrt_s_over_r
	eta_bar = np.mean(B1c) / BBar
	return [BBar,eta_bar,B0,B1s,B1c,B20,B2c,B2s,phi,Psi_a,iotaVMECt,nfp]

print("Creating Bfield, equilibria and xbooz arrays")
equilibria=[]
booz=[]
BBarvec=[]
eta_barvec=[]
B0vec=[]
B1cvec=[]
B1svec=[]
B20vec=[]
B2cvec=[]
B2svec=[]
phiedge=[]
iotaVMEC=[]
nfpVMEC=[]
i=0
for stells in stellDesigns:
	equilibria.append(equilibriaFolder+'wout_'+vmecFiles[i]+'.nc')
	booz.append(equilibriaFolder+'boozmn_'+vmecFiles[i]+'.nc')
	BBart,eta_bart,B0t,B1st,B1ct,B20t,B2ct,B2st,phi,phiedget,iotat,nfpt = obtainB0B1B2(booz[i],Nrotations[i])
	BBarvec.append(BBart)
	eta_barvec.append(eta_bart)
	B0vec.append(B0t)
	B1svec.append(B1st)
	B1cvec.append(B1ct)
	B20vec.append(B20t)
	B2cvec.append(B2ct)
	B2svec.append(B2st)
	phiedge.append(phiedget)
	iotaVMEC.append(iotat)
	nfpVMEC.append(nfpt)
	i=i+1

phi=np.linspace(0,1,N_phi)

if plotSave==0:# path.exists(toPaperFolder+'B0QSEquilibria1.pdf'):
	None
else:
	print("Plotting Near-Axis Bfields for the designs")
	plotsIn1fig=5;plotfontSize=20;legendfontSize=14;figSize1=7.2;figSize2=4.0;
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize);
	plt.figure(figsize=(figSize1,figSize2));i=0;
	for stells in stellDesigns:
		if i<plotsIn1fig: plt.plot(phi, B0vec[i]/np.mean(B0vec[i]), label=stells);i=i+1;
	plt.legend(loc=0,fontsize=legendfontSize);plt.xlabel(r'$N_{fp} \varphi/2 \pi$');plt.ylabel(r'$B_0/\langle B_0 \rangle$');plt.ylim((0.985,1.02))
	plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
	plt.savefig(toPaperFolder+'B0QSEquilibria1.pdf', format='pdf')
	plt.figure(figsize=(figSize1,figSize2));i=0;
	for stells in stellDesigns:
		if i>=plotsIn1fig: plt.plot(phi, B0vec[i]/np.mean(B0vec[i]), label=stells); i=i+1;
		else: i=i+1;
	plt.legend(loc=0,fontsize=legendfontSize);plt.xlabel(r'$N_{fp} \varphi/2 \pi$');plt.ylabel(r'$B_0/\langle B_0 \rangle$');plt.ylim((0.985,1.02))
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
	plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
	plt.savefig(toPaperFolder+'B0QSEquilibria2.pdf', format='pdf')
	plt.figure(figsize=(figSize1,figSize2));i=0;
	for stells in stellDesigns:
		if i<plotsIn1fig: plt.plot(phi, B1cvec[i]/np.mean(B1cvec[i]), label=stells);i=i+1;
		else: i=i+1;
	plt.legend(loc=0,fontsize=legendfontSize);plt.xlabel(r'$N_{fp} \varphi/2 \pi$');plt.ylabel(r'$B_{1c}/\langle B_{1c} \rangle$');plt.ylim((0.7,1.21))
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
	plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
	plt.savefig(toPaperFolder+'B1cQSEquilibria1.pdf', format='pdf')
	plt.figure(figsize=(figSize1,figSize2));i=0;
	for stells in stellDesigns:
		if i>=plotsIn1fig: plt.plot(phi, B1cvec[i]/np.mean(B1cvec[i]), label=stells);i=i+1;
		else: i=i+1;
	plt.legend(loc=0,fontsize=legendfontSize);plt.xlabel(r'$N_{fp} \varphi/2 \pi$');plt.ylabel(r'$B_{1c}/\langle B_{1c} \rangle$');plt.ylim((0.7,1.21))
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
	plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
	plt.savefig(toPaperFolder+'B1cQSEquilibria2.pdf', format='pdf')
	plt.figure(figsize=(figSize1,figSize2));i=0;
	for stells in stellDesigns:
		if i<plotsIn1fig: plt.plot(phi, B1svec[i]/np.mean(B1cvec[i]), label=stells);i=i+1;
		else: i=i+1;
	plt.legend(loc=0,fontsize=legendfontSize);plt.xlabel(r'$N_{fp} \varphi/2 \pi$');plt.ylabel(r'$B_{1s}/\langle B_{1c} \rangle$');plt.ylim((-0.21,0.21))
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
	plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
	plt.savefig(toPaperFolder+'B1sQSEquilibria1.pdf', format='pdf')
	plt.figure(figsize=(figSize1,figSize2));i=0;
	for stells in stellDesigns:
		if i>=plotsIn1fig: plt.plot(phi, B1svec[i]/np.mean(B1cvec[i]), label=stells);i=i+1;
		else: i=i+1;
	plt.legend(loc=0,fontsize=legendfontSize);plt.xlabel(r'$N_{fp} \varphi/2 \pi$');plt.ylabel(r'$B_{1s}/\langle B_{1c} \rangle$');plt.ylim((-0.21,0.21))
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
	plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
	plt.savefig(toPaperFolder+'B1sQSEquilibria2.pdf', format='pdf')
	plt.figure(figsize=(figSize1,figSize2));i=0;
	for stells in stellDesigns:
		if i<plotsIn1fig: plt.plot(phi, B20vec[i]/np.mean(B20vec[i]), label=stells);i=i+1;
		else: i=i+1;
	plt.legend(loc=0,fontsize=legendfontSize);plt.xlabel(r'$N_{fp} \varphi/2 \pi$');plt.ylabel(r'$B_{20}/\langle B_{20} \rangle$');
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
	plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
	plt.savefig(toPaperFolder+'B20QSEquilibria1.pdf', format='pdf')
	plt.figure(figsize=(figSize1,figSize2));i=0;
	for stells in stellDesigns:
		if i>=plotsIn1fig: plt.plot(phi, B20vec[i]/np.mean(B20vec[i]), label=stells);i=i+1;
		else: i=i+1;
	plt.legend(loc=0,fontsize=legendfontSize);plt.xlabel(r'$N_{fp} \varphi/2 \pi$');plt.ylabel(r'$B_{20}/\langle B_{20} \rangle$');
	matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize)
	plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
	plt.savefig(toPaperFolder+'B20QSEquilibria2.pdf', format='pdf')

print("Compiling VMEC to GS2")
os.chdir(vmecGS2interfaceFolder)
process = subprocess.call(['make'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
copyfile("test_vmec_to_gs2_geometry_interface","../test_vmec_to_gs2_geometry_interface")
copymode("test_vmec_to_gs2_geometry_interface","../test_vmec_to_gs2_geometry_interface")
os.chdir('../')

if path.exists(gs2gridsFolder+'agrid'+stellDesigns[0]+'r'+str(normalizedfluxvec[0])+'.out'):
	None
else:
	print("Running VMEC_to_GS2 interface")
	for desired_normalized_toroidal_flux in normalizedfluxvec:
		print("  Normalized toroidal flux = "+str(desired_normalized_toroidal_flux))
		gs2_output=[]
		geometry_output=[]
		i=0
		for stells in stellDesigns:
			gs2_output.append(gs2gridsFolder+'grid'+stells+'r'+str(desired_normalized_toroidal_flux)+'.out')
			geometry_output.append(papergridsFolder+'gridMath'+stells+'r'+str(desired_normalized_toroidal_flux)+'.out')
			i=i+1

		###### Run VMEC_to_GS2 interface
		i=0
		for eq in equilibria:
			f = netcdf.netcdf_file(equilibria[i],'r',mmap=False)
			iota0 = f.variables['iotaf'][()][1]
			nfp   = f.variables['nfp'][()]
			process = subprocess.call(['./test_vmec_to_gs2_geometry_interface',eq,gs2_output[i],geometry_output[i],str(nfp*number_of_field_periods_to_include/abs(iota0)),str(desired_normalized_toroidal_flux)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			i=i+1

###### Mixed Near-Axis/VMEC coefficients
if path.exists(MathDataFolder+'aBcrossGradBdGradAlpha.txt'):
	None
else:
	print("Creating Mixed Near-Axis/VMEC coefficients")
	i=0
	BcrossGradBdGradAlpha=[]
	BcrossGradBdGradPsi=[]
	phiH=[]
	for stells in stellDesigns:
		#print("  Stellarator "+stellDesigns[i])
		B0temp=B0vec[i]
		B1ctemp=B1cvec[i]
		B1stemp=B1svec[i]
		if abs(iotaVMEC[i]-Nrotations[i])>4:
			phiMultiplier=2*np.int(nfpVMEC[i]*4*np.ceil(1/abs(iotaVMEC[i]-Nrotations[i])))
		elif abs(iotaVMEC[i]-Nrotations[i])>3:
			phiMultiplier=2*np.int(nfpVMEC[i]*6*np.ceil(1/abs(iotaVMEC[i]-Nrotations[i])))
		elif abs(iotaVMEC[i]-Nrotations[i])>2:
			phiMultiplier=2*np.int(nfpVMEC[i]*8*np.ceil(1/abs(iotaVMEC[i]-Nrotations[i])))
		elif abs(iotaVMEC[i]-Nrotations[i])>1:
			phiMultiplier=2*np.int(nfpVMEC[i]*8*np.ceil(1/abs(iotaVMEC[i]-Nrotations[i])))
		else:
			phiMultiplier=2*np.int(nfpVMEC[i]*10*np.ceil(1/abs(iotaVMEC[i]-Nrotations[i])))
		for kk in range(phiMultiplier-1):
			B0vec[i]=np.append(B0vec[i],B0temp)
			B1cvec[i]=np.append(B1cvec[i],B1ctemp)
			B1svec[i]=np.append(B1svec[i],B1stemp)
		phiH.append(np.linspace(-np.pi*phiMultiplier/2,np.pi*phiMultiplier/2, N_phi*phiMultiplier))
		vartheta=np.multiply(iotaVMEC[i]-Nrotations[i],phiH[i])

		B1cvec[i]=-B1cvec[i]
		
		#plotfontSize=20;legendfontSize=14;figSize1=7.5;figSize2=4.0;
		#matplotlib.rc('font', size=plotfontSize);matplotlib.rc('axes', titlesize=plotfontSize);
		#matplotlib.rc('text', usetex=True);matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
		#plt.figure(figsize=(figSize1,figSize2))
		for s in normalizedfluxvec:
			r=np.sqrt(2*phiedge[i]*s/BBarvec[i])
			Bmag=(B0vec[i]+np.multiply(r,np.multiply(B1cvec[i],np.cos(vartheta))+np.multiply(B1svec[i],np.sin(vartheta))))
			#plt.plot(phiH[i], Bmag, label=stells+' s='+str(s))
		#plt.legend(loc=0,fontsize=legendfontSize);plt.xlabel(r'$\varphi$');plt.ylabel(r'$B$');#plt.ylim((0.985,1.015))
		#plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
		#plt.savefig(figuresFolder+stells+'_BmagMixed.pdf', format='pdf')

		#plt.figure(figsize=(figSize1,figSize2))
		#BcrossGradBdGradAlpha.append(np.multiply(np.multiply(np.sqrt(BBarvec[i]/2),np.multiply(np.multiply(B0vec[i],B1cvec[i]),np.cos(vartheta))),np.power(np.reciprocal(B0vec[i]),3)))
		BcrossGradBdGradAlpha.append(np.multiply(np.multiply(np.multiply(B1cvec[i],np.cos(vartheta))+np.multiply(B1svec[i],np.sin(vartheta)),np.reciprocal(B0vec[i])),np.sqrt(1/(2*BBarvec[i]))))

		#plt.plot(phiH[i], BcrossGradBdGradAlpha[i])
		#plt.xlabel(r'$\varphi$');plt.ylabel(r'$\sqrt{\psi}\boldsymbol B \times \nabla B \cdot \nabla \alpha/B^3$');#plt.ylim((0.985,1.015))
		#plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
		#plt.savefig(figuresFolder+stells+'_BcrossGradBdGradAlphaMixed.pdf', format='pdf')
		
		#plt.figure(figsize=(figSize1,figSize2))
		#BcrossGradBdGradPsi.append(np.multiply(np.multiply(np.sqrt(2/BBarvec[i]),np.multiply(np.multiply(np.square(B0vec[i]),-B1cvec[i]),np.sin(vartheta))),np.power(np.reciprocal(B0vec[i]),3)))
		BcrossGradBdGradPsi.append(np.multiply(np.multiply(np.multiply(B1cvec[i],-np.sin(vartheta))-np.multiply(B1svec[i],np.cos(vartheta)),np.reciprocal(B0vec[i])),np.sqrt(2/(BBarvec[i]))))

		#plt.plot(phiH[i], BcrossGradBdGradPsi[i])
		#plt.xlabel(r'$\varphi$');plt.ylabel(r'$\sqrt{\psi}^{-1} \boldsymbol B \times \nabla B \cdot \nabla \psi/B^3$');#plt.ylim((0.985,1.015))
		#plt.subplots_adjust(left=0.16, bottom=0.19, right=0.98, top=0.97)
		#plt.savefig(figuresFolder+stells+'_BcrossGradBdGradPsiMixed.pdf', format='pdf')
		i=i+1

	with open(MathDataFolder+'phiH.txt','w+') as phif:
		with open(MathDataFolder+'BcrossGradBdGradAlpha.txt','w+') as BcrossGradBdGradAlphaf:
			with open(MathDataFolder+'BcrossGradBdGradPsi.txt','w+') as BcrossGradBdGradPsif:
				i=0
				#phif.write('{');BcrossGradBdGradAlphaf.write('{');BcrossGradBdGradPsif.write('{')
				for stells in stellDesigns:
					#phif.write('{');BcrossGradBdGradAlphaf.write('{');BcrossGradBdGradPsif.write('{')
					np.savetxt(phif, phiH[i], fmt='%s',newline=' ')
					np.savetxt(BcrossGradBdGradAlphaf, BcrossGradBdGradAlpha[i], fmt='%s',newline=' ')
					np.savetxt(BcrossGradBdGradPsif, BcrossGradBdGradPsi[i], fmt='%s',newline=' ')
					#phif.write('}');BcrossGradBdGradAlphaf.write('}');BcrossGradBdGradPsif.write('}')
					phif.write('\n');BcrossGradBdGradAlphaf.write('\n');BcrossGradBdGradPsif.write('\n')
					i=i+1
				#phif.write('}');BcrossGradBdGradAlphaf.write('}');BcrossGradBdGradPsif.write('}')

###### Near-Axis Comparison
print("Running Mathematica")
i=0
for eq in equilibria:
	bashCommand = "wolframscript -noprompt -script grid_NearAxis.wls "+stellDesigns[i]+" "+papergridsFolder+" "+eq+" "+booz[i]+" "+str(abs(eta_barvec[i]))+" "+figuresFolder+" "+gs2gridsFolder+" "+MathDataFolder+" "+toPaperFolder+" "+str(plotSave)+" "+str(BBarvec[i])+" "+str(i+1)+" "+str(len(normalizedfluxvec))
	for rr in normalizedfluxvec:
		bashCommand = bashCommand+" "+str(rr)
	print(stellDesigns[i])
	print("etabar = "+str(abs(eta_barvec[i])))
	print("B0 = "+str(BBarvec[i]))
	print("Working...", end='', flush=True)
	#print(bashCommand)
	#exit()
	output = subprocess.call(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	print(" Success!")
	i=i+1

#pkill -9 WolframKernel; pkill -9 WolframScript; pkill -9 Mathematica








