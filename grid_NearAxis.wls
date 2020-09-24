#!/usr/bin/env wolframscript
(* ::Package:: *)

LaunchKernels[];

ParallelNeeds["PacletManager`"]
PacletInstall["MaTeX-1.7.6.paclet",IgnoreVersion->True];
<< MaTeX`;
SetOptions[$Output, FormatType->OutputForm];
currentPath=DirectoryName @ $InputFileName;
SetDirectory[currentPath];
pi = N[Pi,10];

nperiod=1
drhodpsi=1
rmaj=1
kxfac=1

stell            = $ScriptCommandLine[[2]];
papergridsFolder = $ScriptCommandLine[[3]];
VMECfileIn       = $ScriptCommandLine[[4]];
BOOZfileIn       = $ScriptCommandLine[[5]];
etabar           = ToExpression@$ScriptCommandLine[[6]];
figuresFolder    = $ScriptCommandLine[[7]];
gs2gridsFolder   = $ScriptCommandLine[[8]];
MathDataFolder   = $ScriptCommandLine[[9]];
toPaperFolder    = $ScriptCommandLine[[10]];
plotSave         = ToExpression@$ScriptCommandLine[[11]];
B0               = ToExpression@$ScriptCommandLine[[12]];
(*B0 = Import[BOOZfileIn, {"Datasets", "bmnc_b"}][[1, 1]];*)
stellN           = ToExpression@$ScriptCommandLine[[13]];
nradius          = ToExpression@$ScriptCommandLine[[14]];
radiusvec        = Table[ToExpression@$ScriptCommandLine[[14+ii]],{ii,1,nradius}];
saveStellFile    = MathDataFolder<>stell<>"NA_MathData.mx";

If[FileExistsQ[saveStellFile],
	DumpGet[saveStellFile];Aminor = Import[VMECfileIn, {"Datasets", "Aminor_p"}];,
	ntgrid = 128; nlambda[j] = 50; nperiod = 1; drhodpsi = 1.0; rmaj = 1;
	shat[j] = 1.0; kxfac = 1; q = 1; sigma0 = 0; ns = 856; iota0 = 0.1;
	nz = 2*ntgrid; x = 0.08; y = 0.0; alpha0 = pi;
	raxis = Import[VMECfileIn, {"Datasets", "raxis_cc"}];
	zaxis = Import[VMECfileIn, {"Datasets", "zaxis_cs"}];
	R0 = raxis[[1]]; eR = raxis[[2 ;; 5]]; eZ = zaxis[[2 ;; 5]];
	NFP = Import[VMECfileIn, {"Datasets", "nfp"}];
	Aminor = Import[VMECfileIn, {"Datasets", "Aminor_p"}];
	phiEDGE = Last[Import[VMECfileIn, {"Datasets", "phi"}]]/(2 pi);
	(*B0 = Import[BOOZfileIn, {"Datasets", "bmnc_b"}][[1, 1]];*)
	I2oB0 = 0.;
	nIts = 10; nLineS = 4; tol = 10^-6;
	sigmaIn = 0.1 Cos[NFP t];
	nzMinMax = 0.0; (*number of 2pi additions on theta*)
	(*Axis Equations*)
	R[t_] = R0 + Sum[eR[[n]]*Cos[n NFP t], {n, 1, Length[eR]}];
	Z[t_] = Sum[eZ[[n]]*Sin[n NFP t], {n, 1, Length[eZ]}];
	closedCurve[t_] = {R[t] Cos[t], R[t] Sin[t], Z[t]};
	curvFunc[t_]    = Chop[ComplexExpand[Norm[Cross[closedCurve'[t], closedCurve''[t]]]]/ComplexExpand[Norm[closedCurve'[t]]^3], 10^-10];
	torsFunc[t_]    = Chop[Dot[Cross[closedCurve'[t], closedCurve''[t]], closedCurve'''[t]]/ComplexExpand[Norm[Cross[closedCurve'[t], closedCurve''[t]]]]^2, 10^-10];
	sprimeFunc[t_]  = Chop[ComplexExpand[Norm[closedCurve'[t]]], 10^-10] // Quiet;
	(*Fast Numerical Integrator*)
	GaussLegendreQuadrature[f_, {x_, a_, b_}, n_Integer: 10, prec_: MachinePrecision] := 
		Module[{nodes, weights},
			{nodes, weights} = Most[NIntegrate`GaussRuleData[n, prec]];
			(b - a) weights.Map[Function[x, f], Rescale[nodes, {0, 1}, {a, b}]]
		];
	(*Axis Length*)
	Laxis = GaussLegendreQuadrature[sprimeFunc[t], {t, 0, 2 pi}, 250, 20];
	(*Number of Rotations of Normal Vector*)
	b0[t_] = Chop[closedCurve'[t]/sprimeFunc[t]];
	k0[t_] = Chop[1/curvFunc[t] b0'[t]/sprimeFunc[t]];
	quadrant = 0; nNormal = 0; diffQuadrant = 0; qnew = 0;
	Table[normalx = k0[t][[1]]*Cos[t] + k0[t][[2]]*Sin[t];
	  normaly = k0[t][[3]];
	  If[Abs[normalx] > Abs[normaly], If[normalx > 0, qnew = 2, qnew = 4],
	    If[normaly > 0, qnew = 1, qnew = 3]];
	  If[quadrant == 0, quadrant = qnew, 
	   diffQuadrant = qnew - quadrant;];
	  If[Abs[diffQuadrant] == 3, diffQuadrant = diffQuadrant/3]; 
	  nNormal = nNormal + diffQuadrant/2; 
	  quadrant = qnew;, {t, 0, 2*Pi, 2*Pi/NFP/15}];
	(*Convert VMEC phi to axis phi*)
	phiAxis=Compile[{{phi,_Real}},(2*pi/Laxis)*GaussLegendreQuadrature[sprimeFunc[t],{t,0,phi},180,12]];
	phiA=Interpolation[Table[Evaluate@{t,phiAxis[t]},{t,-40*2.*pi,40*2.*pi,(40.*pi)/1500.}]]//Quiet;
	(*Compiled sigma equation and its derivatives*)
	F0[Iota_]:=(2 Pi)/Laxis (Iota-nNormal);
	sigmap[t_,Iota_,Sigma_,etab_]:=-F0[Iota] (1+Sigma^2+etab^4/curvFunc[t]^4)-(2 (etab^2) )/curvFunc[t]^2 (torsFunc[t]-I2oB0);
	fsigma=Compile[{{t,_Real},{Iota,_Real},{Sigma,_Real},{etab,_Real}},Evaluate[-sprimeFunc[t]sigmap[t,Iota,Sigma,etab]]]//Quiet;
	dfdsigma=Compile[{{t,_Real},{Iota,_Real},{Sigma,_Real},{etab,_Real}},Evaluate[-sprimeFunc[t]D[sigmap[t,Iota,Sigma,etab],Sigma]]]//Quiet;
	dfdiota0=Compile[{{t,_Real},{Iota,_Real},{Sigma,_Real},{etab,_Real}},Evaluate[-sprimeFunc[t]D[sigmap[t,Iota,Sigma,etab],Iota]]]//Quiet;
	(*Residual matrix*)
	Clear[qsresidualsigma];
	qsresidualsigma[Dmat_,sigmaVec_,iota_,etab_,ns_,sigma0_]:=Module[
		{delta=(2 pi)/ns },
		fsigmaMat=Table[fsigma[delta i,iota,sigmaVec[[i]],etab],{i,1,ns}];
		dfdsigmaMat=Table[dfdsigma[delta i,iota,sigmaVec[[i]],etab],{i,1,ns}];
		dfdiota0Mat=Table[dfdiota0[delta i,iota,sigmaVec[[i]],etab],{i,1,ns}];
		f1Mat={
			Dmat.sigmaVec+fsigmaMat,
			sigmaVec[[ns]]-sigma0
		}//Flatten;
		myArrayFlatten=Flatten/@Flatten[#,{{1,3}}]&;
		Jacf1Mat=myArrayFlatten@{
			{Dmat+DiagonalMatrix[dfdsigmaMat],dfdiota0Mat}
		};
		Jacf1Mat=Append[Jacf1Mat,{ConstantArray[0, ns-1],1,0}//Flatten];
		{f1Mat,Jacf1Mat}
	];
	(*Newton's Method*)
	Clear[solveSigma,Dmat,res];
	solveSigma[nIts_,tol_,nLineS_,sigmaIn_,iota00_,etab_,ns_,sigma0_]:=Module[{deltaold=0 },
		iotanew=iota00;
		delta=(2 pi)/ns;
		sigmaVec=Table[sigmaIn/.t->i delta,{i,1,ns}];
		(*Derivative Matrix*)
		column1=Prepend[ Table[ 0.5 (-1)^n Cot[n/2 delta],{n,1,ns-1}],0];
		column2=RotateRight[Reverse[column1]];
		Dmat=ToeplitzMatrix[column1,column2];
		For[i=1,i<=nIts,i++,
			{f1,Jacf1}=qsresidualsigma[Dmat,sigmaVec,iotanew,etab,ns,sigma0];
			deltay=-LinearSolve[Jacf1,f1]//Quiet;
			res[i]=Norm[deltay];
			If[i>5&&(Abs[Norm[deltay]-Norm[deltaold]])/Sqrt[ns]<tol,Break[]];
			For[j=1,j<nLineS,j++,
				sigmatemp=sigmaVec+deltay[[1;;ns]];
				iota0temp=iotanew+deltay[[ns+1]];
				{f2,Jacf2}=qsresidualsigma[Dmat,sigmatemp,iota0temp,etab,ns,sigma0];
				If[Abs[Norm[f2]]<Abs[Norm[f1]]/1.2,Break[],deltay=deltay/2];
			];
			sigmaVec=sigmatemp;
			iotanew=iota0temp;
			deltaold=deltay;
		];
		{sigmaVec,iotanew}
	];
	(*Solve Sigma Equation*)
	{sigmaSol, iota} = solveSigma[nIts, tol, nLineS, sigmaIn, iota0, etabar, ns, sigma0];
	sigmatemp = ListInterpolation[Flatten@sigmaSol, {0, 2 pi}, Method -> "Spline", InterpolationOrder -> 3];
	sigma[y_] := sigmatemp[Mod[y, 2 pi]];
	Print["iota = "<>ToString[AccountingForm[Chop[iota,10^(-7)],7]]];
	DumpSave[saveStellFile, {sigma,iota,curvFunc,torsFunc,sprimeFunc,etabar,B0,nNormal,Aminor,phiEDGE,sigmatemp,Laxis,phiA}];
];

PlotOptions[tit_] = {PlotRange -> All, AxesLabel -> {Style["\[Phi]", 20, FontFamily -> "Latin Modern Math"], 
	     			Style["", 20, FontFamily -> "Latin Modern Math"]}, ImageSize -> Medium, 
	   				LabelStyle -> Directive[16, FontFamily -> "Latin Modern Math"], 
	   				PlotStyle -> Thickness[0.01], PlotLabel -> Style[tit, 20, FontFamily -> "Latin Modern Math"]};
PlotOptions2[tit_] = {PlotRange -> All, AxesLabel -> {Style["\[Phi]", 20, FontFamily -> "Latin Modern Math"], 
	     			Style[tit, 20, FontFamily -> "Latin Modern Math"]}, ImageSize -> Medium, 
	   				LabelStyle -> Directive[16, FontFamily -> "Latin Modern Math"], PlotStyle -> Thickness[0.01]};
textMagnification = 1.7;
Table[
	(*<< MaTeX`;*)
	fileIn            = papergridsFolder<>"gridMath"<>stell<>"r"<>ToString@radiusvec[[j]]<>".out";
	(*Grid Parameters*)
	tgridVMEC[j]    	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "tgrid("][[1]]]],StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	gbdriftVMEC[j]  	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "gbdrift("][[1]]]], StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	gbdrift0VMEC[j] 	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "gbdrift0("][[1]]]], StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	cvdriftVMEC[j]  	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "cvdrift("][[1]]]], StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	cvdrift0VMEC[j] 	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "cvdrift0("][[1]]]], StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	gradparVMEC[j]  	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "gradpar("][[1]]]], StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	bmagVMEC[j]     	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "bmag("][[1]]]], StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	gds2VMEC[j]     	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "gds2("][[1]]]], StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	gds21VMEC[j]    	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "gds21("][[1]]]],StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^","E-" :> "*^-"}];
	gds22VMEC[j]    	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "gds22("][[1]]]],StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	lambdaVMEC[j]    	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "lambda("][[1]]]],StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
	shat[j]         	  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "shat"][[1]]]], StringSplit[{}[[1]]] // Quiet], {"E+" :> "*^", "E-" :> "*^-"}][[2]];
	normalizedtorFlux[j]  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "normalized_flux"][[1]]]], StringSplit[{}[[1]]] // Quiet], {"E+" :> "*^", "E-" :> "*^-"}][[2]];
	alphaVMEC[j]          = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "alpha"][[1]]]],StringSplit[{}[[1]]] // Quiet], {"E+" :> "*^", "E-" :> "*^-"}][[2]];
	nlambda[j]            = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "nlambda"][[1]]]],StringSplit[{}[[1]]] // Quiet], {"E+" :> "*^", "E-" :> "*^-"}][[2]];
	(*Geometry Definitions*)
	rVMEC[j] 			  = -Sqrt[((2 phiEDGE normalizedtorFlux[j])/B0)];
	Phi[theta_,j] 	      = (theta - alphaVMEC[j])/(iota - nNormal);(*phiA[(theta-alphaVMEC[j])/(iota-nNormal)];*)
	bmagNew[theta_,j]     = (Aminor^2 B0 (1 + rVMEC[j] etabar Cos[theta]))/(2 phiEDGE);
	gradparNew[theta_,j]  = ((2 Aminor pi (1 + rVMEC[j] etabar Cos[theta]))/Laxis)/(sprimeFunc[(alphaVMEC[j]-theta)/(iota-nNormal)]*2*pi/Laxis);
	gds2New[theta_,j]     = ((Aminor^2) B0 )/(2 phiEDGE) ((etabar^2 Cos[theta]^2)/curvFunc[Phi[theta,j]]^2 + (curvFunc[Phi[theta,j]]^2 (Sin[theta] + Cos[theta] sigma[Phi[theta,j]])^2)/etabar^2);
	gds21New[theta_,j]    = -1/(2 phiEDGE) Aminor^2  shat[j] ((B0 etabar^2 Cos[theta] Sin[theta])/curvFunc[Phi[theta,j]]^2 + 1/etabar^2 B0 curvFunc[Phi[theta,j]]^2 (Sin[theta] + Cos[theta] sigma[Phi[theta,j]]) (-Cos[theta] + Sin[theta] sigma[Phi[theta,j]]));
	gds22New[theta_,j]    = (Aminor^2 B0 shat[j]^2 (etabar^4 Sin[theta]^2 + curvFunc[Phi[theta,j]]^4 (Cos[theta] - Sin[theta] sigma[Phi[theta,j]])^2))/(2 phiEDGE etabar^2 curvFunc[Phi[theta,j]]^2);
	gbdriftNew[theta_,j]  = (2 Sqrt[2] etabar Cos[theta])/Sqrt[B0/phiEDGE] (1 - 0 2 rVMEC[j] etabar Cos[theta]); 
	cvdriftNew[theta_,j]  = gbdriftNew[theta,j];
	gbdrift0New[theta_,j] = -2 Sqrt[2] Sqrt[phiEDGE/B0] shat[j] etabar Sin[theta] (1 - 0 2 rVMEC[j] etabar Cos[theta]); 
	cvdrift0New[theta_,j] = gbdrift0New[theta,j];
	zminNA[j]		    = x/.Last[FindMaximum[bmagNew[x,j],{x,Min[tgridVMEC[j]]}]];
	zmaxNA[j]		    = x/.Last[FindMaximum[bmagNew[x,j],{x,Max[tgridVMEC[j]]}]];
	tgridNA[j] = Table[i,{i,zminNA[j],zmaxNA[j],(zmaxNA[j]-zminNA[j])/(Length@tgridVMEC[j]-1)}];
	zmin = Min[tgridVMEC[j]];
	zmax = Max[tgridVMEC[j]];
	lambdamin=((2 phiEDGE)/(Aminor^2 B0))/(1+Abs[rVMEC[j]*etabar]);
	lambdamax=((2 phiEDGE)/(Aminor^2 B0))/(1-Abs[rVMEC[j]*etabar]);
	lambdaNA=Table[i,{i,lambdamin,lambdamax,(lambdamax-lambdamin)/(nlambda[j]-1)}];
	(*Comparison of GS2 quantities*)
	(*
	plotGrid1 = Grid[{
	   {Show[Plot[gds2New[   z (iota - nNormal),j], {z, zmin, zmax}, PlotStyle -> {Red, Dashed}], ListLinePlot[Transpose@{tgridVMEC[j], gds2VMEC[j]}],  AxesOrigin -> {0, 0.6 Min[Quiet@FindMinimum[gds2New[z,j], {z, 0.6}][[1]], Min[gds2VMEC[j]]]}, Evaluate@PlotOptions["gds2"]],
	    Show[Plot[gds21New[  z (iota - nNormal),j], {z, zmin, zmax}, PlotStyle -> {Red, Dashed}], ListLinePlot[Transpose@{tgridVMEC[j], gds21VMEC[j]}], AxesOrigin -> {0, 1.2 Min[Quiet@FindMinimum[gds21New[z,j], {z, 0.1}][[1]], Min[gds21VMEC[j]]]}, Evaluate@PlotOptions["gds21"]],
	    Show[Plot[gds22New[  z (iota - nNormal),j], {z, zmin, zmax}, PlotStyle -> {Red, Dashed},  PlotLegends -> LineLegend[{Style["Near-Axis", 20, FontFamily -> "Latin Modern Math"]}]], ListLinePlot[Transpose@{tgridVMEC[j], gds22VMEC[j]}, PlotLegends -> LineLegend[{Style["VMEC", 20, FontFamily -> "Latin Modern Math"]}]], AxesOrigin -> {0, Min[0 Quiet@FindMinimum[gds22New[z], {z, 0.6}][[1]], Min[gds22VMEC[j]]]}, Evaluate@PlotOptions["gds22"]]},
	   {Show[Plot[gbdriftNew[z (iota - nNormal),j], {z, zmin, zmax}, PlotStyle -> {Red, Dashed}], ListLinePlot[Transpose@{tgridVMEC[j], gbdriftVMEC[j]}], AxesOrigin -> {0, Min[gbdriftNew[zmin (iota - nNormal),j], Min[gbdriftVMEC[j]]]}, Evaluate@PlotOptions["gbdrift"]],
	    Show[Plot[cvdriftNew[z (iota - nNormal),j], {z, zmin, zmax}, PlotStyle -> {Red, Dashed}], ListLinePlot[Transpose@{tgridVMEC[j], cvdriftVMEC[j]}], AxesOrigin -> {0, Min[cvdriftNew[zmin (iota - nNormal),j], Min[cvdriftVMEC[j]]]}, Evaluate@PlotOptions["cvdrift"]],
	    Show[ListLinePlot[lambdaVMEC[j],PlotLegends -> LineLegend[{Style["VMEC",20,FontFamily -> "Latin Modern Math"]}]],ListLinePlot[lambdaNA,PlotStyle -> {Red,Dashed},PlotLegends -> LineLegend[{Style["Near-Axis",20,FontFamily -> "Latin Modern Math"]}]],Evaluate@PlotOptions["lambda"]]},
	   {Show[ListLinePlot[Transpose@{tgridVMEC[j], gbdrift0VMEC[j]}], Plot[gbdrift0New[z (iota - nNormal),j], {z, zmin, zmax}, PlotStyle -> {Red, Dashed}], AxesOrigin -> {0, 0}, Evaluate@PlotOptions["gddrift0"]],
	    Show[ListLinePlot[Transpose@{tgridVMEC[j], cvdrift0VMEC[j]}], Plot[cvdrift0New[z (iota - nNormal),j], {z, zmin, zmax}, PlotStyle -> {Red, Dashed}], AxesOrigin -> {0, 0}, Evaluate@PlotOptions["cvdrift0"]],
	    Show[Plot[bmagNew[   z (iota - nNormal),j], {z, zmin, zmax}, PlotStyle -> {Red, Dashed}, PlotLegends -> LineLegend[{Style["Near-Axis", 20, FontFamily -> "Latin Modern Math"]}]], ListLinePlot[Transpose@{tgridVMEC[j], bmagVMEC[j]}, PlotLegends -> LineLegend[{Style["VMEC", 20, FontFamily -> "Latin Modern Math"]}]](*,AxesOrigin -> {0,Min[0.96bmagNew[zmin (iota-nNormal)],Min[bmagVMEC[j]]]}*), Evaluate@PlotOptions["bmag"]]}
	   }];
	Export[figuresFolder<>stell<>"r"<>ToString@normalizedtorFlux[j]<>"_comparison.pdf", plotGrid1];
	Export[MathDataFolder<>stell<>"r"<>ToString@normalizedtorFlux[j]<>"_comparison.mx", "AllowRasterization" -> True];
	*)
	(*Comparison of geometric quantities*)
	BstrengthNA[theta_,j]					 = bmagNew[theta,j]*(2 phiEDGE)/Aminor^2;
	BstrengthGK[j]          				 = bmagVMEC[j]*(2 phiEDGE)/Aminor^2;
	gradAlpha2r2NA[theta_,j]  		         = (B0/2)*((2 phiEDGE)/((Aminor^2) B0 ))*gds2New[theta,j];
	gradAlpha2r2GK[j]          		         = (B0/2)*((2 phiEDGE)/((Aminor^2) B0 ))*gds2VMEC[j];
	gradAlphaPsiNA[theta_,j]  		         = ((2 phiEDGE)/(Aminor^2 shat[j]))*gds21New[theta,j];
	gradAlphaPsiGK[j]          		         = ((2 phiEDGE)/(Aminor^2 shat[j]))*gds21VMEC[j];
	gradPsi2r2NA[theta_,j]    		         = ((2 B0 phiEDGE)/(Aminor^2 shat[j]^2))*gds22New[theta,j]*2/B0;
	gradPsi2r2GK[j]            	             = ((2 B0 phiEDGE)/(Aminor^2 shat[j]^2))*gds22VMEC[j]*2/B0;
	BcrossGradBdotGradAlphaOB3NA[theta_,j]   = Sqrt[B0/2]*gbdriftNew[theta,j]/(2 Sqrt[2] Sqrt[B0 phiEDGE]);
	BcrossGradBdotGradAlphaOB3GK[j]		     = Sqrt[B0/2]*gbdriftVMEC[j]/(2 Sqrt[2] Sqrt[B0 phiEDGE]);
	BcrossGradBdotGradPsiOB3NA[theta_,j]     = Sqrt[2/B0]*gbdrift0New[theta,j]*Sqrt[B0/phiEDGE]/(2*Sqrt[2]*shat[j]);
	BcrossGradBdotGradPsiOB3GK[j]		     = Sqrt[2/B0]*gbdrift0VMEC[j]*Sqrt[B0/phiEDGE]/(2*Sqrt[2]*shat[j]);

	(*Output to GS2 grid*)
	ToStringA[x_]:=ToString[FortranForm@N[1. Chop[x,10^-16],10]];
	nzNA=Length[tgridNA[j]];
	ntgrid=Floor[nzNA/2];
	paramTheta[i_]:=tgridNA[j][[i]]*(iota-nNormal);
	paramZ[z_]:=tgridNA[j][[i]];
	text="nlambda\n"<>ToString[nlambda[j]]<>"\nlambda";
	For[i=1,i<=nlambda[j],i++,text=text<>"\n"<>ToString[lambdaNA[[i]]]];
	text=text<>"\nntgrid nperiod ntheta drhodpsi rmaj shat kxfac q";
	text=text<>"\n"<>ToString@ntgrid<>" "<>ToString@nperiod<>" "<>ToString@nzNA<>" "<>ToString@drhodpsi<>" "<>ToString@rmaj<>" "<>ToString@shat[j]<>" "<>ToString@kxfac<>" "<>ToStringA[1/iota];
	text=text<>"\ngbdrift gradpar grho tgrid";
	For[i=1,i<=nzNA,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[gbdriftNew[zz,j]]<>" "<>ToStringA[gradparNew[zz,j]]<>" "<>ToStringA[1]<>" "<>ToStringA[paramZ[i]]];
	text=text<>"\ncvdrift gds2 bmag tgrid";
	For[i=1,i<=nzNA,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[cvdriftNew[zz,j]]<>" "<>ToStringA[gds2New[zz,j]]<>" "<>ToStringA[bmagNew[zz,j]]<>" "<>ToStringA[paramZ[i]]];
	text=text<>"\ngds21 gds22 tgrid";
	For[i=1,i<=nzNA,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[gds21New[zz,j]]<>" "<>ToStringA[gds22New[zz,j]]<>" "<>ToStringA[paramZ[i]]];
	text=text<>"\ncvdrift0 gbdrift0 tgrid";
	For[i=1,i<=nzNA,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[cvdrift0New[zz,j]]<>" "<>ToStringA[gbdrift0New[zz,j]]<>" "<>ToStringA[paramZ[i]]];
	text=text<>"\nRplot Rprime tgrid";
	For[i=1,i<=nzNA,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[0.0]<>" "<>ToStringA[0.0]<>" "<>ToStringA[paramZ[i]]];
	text=text<>"\nZplot Zprime tgrid";
	For[i=1,i<=nzNA,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[0.0]<>" "<>ToStringA[0.0]<>" "<>ToStringA[paramZ[i]]];
	text=text<>"\naplot aprime tgrid";
	For[i=1,i<=nzNA,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[0.0]<>" "<>ToStringA[0.0]<>" "<>ToStringA[paramZ[i]]];
	Export[gs2gridsFolder<>"grid"<>stell<>"r"<>ToString@normalizedtorFlux[j]<>"NA.out",text,"String",CharacterEncoding->"ASCII"];
	
	(*Export Everything and Exit*)
	(*
	DumpSave[MathDataFolder<>"/MathData"<>stell<>"r"<>ToString@normalizedtorFlux[j]<>".mx", "Global`"];
	*)
,{j,1,nradius}];

If[plotSave==1,

phiHT						 = Import[MathDataFolder<>"phiH.txt","Table"][[stellN]];
BcrossGradBdGradAlphaMixedT  = Import[MathDataFolder<>"BcrossGradBdGradAlpha.txt","Table"][[stellN]];
BcrossGradBdGradPsiMixedT    = Import[MathDataFolder<>"BcrossGradBdGradPsi.txt","Table"][[stellN]];

BcrossGradBdGradAlphaMixed = ListInterpolation[BcrossGradBdGradAlphaMixedT, {{Min[phiHT], Max[phiHT]}}];
BcrossGradBdGradPsiMixed   = ListInterpolation[BcrossGradBdGradPsiMixedT,   {{Min[phiHT], Max[phiHT]}}];

(* Plot Figures for Paper*)
legNA=Table[Row[{"Near-Axis s = ", N @ radiusvec[[j]]}],{j,1,nradius}];
legVMEC=Table[Row[{stell<>" s = ", N @ radiusvec[[j]]}],{j,1,nradius}];
thetaBooz[z_,j_]=alphaVMEC[j]+(iota - nNormal)*phiA[z];

rplotBmag=Show[
			   Plot[Evaluate[Table[BstrengthNA[thetaBooz[z,j],j],{j,1,nradius}]], {z, zminNA[1], zmaxNA[1],}(*, PlotLegends->legNA*), PlotStyle->{Dashed,Dashed,Dashed,Dashed,Dashed,Dashed}],
			   ListLinePlot[Evaluate[Table[Transpose@{tgridVMEC[j], BstrengthGK[j]},{j,1,nradius}]](*,PlotLegends->legVMEC*)]
		  ,Evaluate@PlotOptions2[MaTeX["B",Magnification->textMagnification]]];
rplotgradAlpha2r2=Show[
					   ListLinePlot[Evaluate[Table[Transpose@{tgridVMEC[j], gradAlpha2r2GK[j]},{j,1,nradius}]](*,PlotLegends->legVMEC*)],
					   Plot[Evaluate[Table[gradAlpha2r2NA[thetaBooz[z,j],j],{j,1,1}]], {z, zminNA[1],, zmaxNA[1],},PlotStyle->Black(*,PlotLegends->{"Near-Axis"}*)]
	      		  , Evaluate@PlotOptions2[MaTeX["\\psi |\\nabla \\alpha|^2",Magnification->textMagnification]]];
rplotgradPsi2r2=Show[
			         ListLinePlot[Evaluate[Table[Transpose@{tgridVMEC[j], gradPsi2r2GK[j]},{j,1,nradius}]](*,PlotLegends->legVMEC*)],
			         Plot[Evaluate[Table[gradPsi2r2NA[thetaBooz[z,j],j],{j,1,1}]], {z, zminNA[1],, zmaxNA[1],},PlotRange->All,PlotStyle->Black(*,PlotLegends->{"Near-Axis"}*)]
		  		, Evaluate@PlotOptions2[MaTeX["|\\nabla \\psi|^2/\psi",Magnification->textMagnification]]];
rplotgradAlphaPsi=Show[
			  		   ListLinePlot[Evaluate[Table[Transpose@{tgridVMEC[j], gradAlphaPsiGK[j]},{j,1,nradius}]](*,PlotLegends->legVMEC*)],
			  		   Plot[Evaluate[Table[gradAlphaPsiNA[thetaBooz[z,j],j],{j,1,1}]], {z, zminNA[1],, zmaxNA[1],},PlotStyle->Black(*,PlotLegends->{"Near-Axis"}*)]
		  		  , Evaluate@PlotOptions2[MaTeX["\\nabla \\psi \\cdot \\nabla \\alpha",Magnification->textMagnification]]];
rplotBcrossGradBdotGradAlphaOB3=Show[
									 ListLinePlot[Evaluate[Table[Transpose@{tgridVMEC[j], BcrossGradBdotGradAlphaOB3GK[j]},{j,1,nradius}]](*,PlotLegends->legVMEC*)],
									 Plot[Evaluate[Table[BcrossGradBdotGradAlphaOB3NA[thetaBooz[z,j],j],{j,1,1}]], {z, zminNA[1],, zmaxNA[1],},PlotStyle->Black(*,PlotLegends->{"Near-Axis"}*)],
									 Plot[BcrossGradBdGradAlphaMixed[phiA[z]], {z, zminNA[1],, zmaxNA[1],},PlotStyle->{Magenta, DotDashed}(*,PlotLegends->{"Near-Axis"}*)]
			  					, Evaluate@PlotOptions2[MaTeX["{\\sqrt \\psi} \\mathbf B \\times \\nabla B \\cdot \\nabla \\alpha/B^3",Magnification->textMagnification]]];
rplotBcrossGradBdotGradPsiOB3=Show[
								   ListLinePlot[Evaluate[Table[Transpose@{tgridVMEC[j], BcrossGradBdotGradPsiOB3GK[j]},{j,1,nradius}]](*,PlotLegends->legVMEC*)],
								   Plot[Evaluate[Table[BcrossGradBdotGradPsiOB3NA[thetaBooz[z,j],j],{j,1,1}]], {z, zminNA[1],, zmaxNA[1],},PlotStyle->Black(*,PlotLegends->{"Near-Axis"}*)],
								   Plot[BcrossGradBdGradPsiMixed[phiA[z]], {z, zminNA[1],, zmaxNA[1],},PlotStyle->{Magenta,DotDashed}(*,PlotLegends->{"Near-Axis"}*)]
			  				  , Evaluate@PlotOptions2[MaTeX["{\\sqrt \\psi}^{-1}\\mathbf B \\times \\nabla B \\cdot \\nabla \\psi/B^3",Magnification->textMagnification]]];

legFont[tit_] = Style[tit, 19, FontFamily -> "Latin Modern Math"];
paperLeg = LineLegend[
  {Directive[ColorData[97, "ColorList"][[1]], Dashed], Directive[ColorData[97, "ColorList"][[2]], Dashed],
  Directive[ColorData[97, "ColorList"][[3]], Dashed], Directive[ColorData[97, "ColorList"][[4]], Dashed],
  ColorData[97, "ColorList"][[1]], ColorData[97, "ColorList"][[2]], ColorData[97, "ColorList"][[3]],
  ColorData[97, "ColorList"][[4]], Black, Directive[Magenta, DotDashed]}
  , {legFont["Near-Axis s=0.01"], legFont["s=0.05"], legFont["s=0.1"],
     legFont["s=0.3"], legFont[stell <> " s=0.01"], legFont["s=0.05"],
     legFont["s=0.1"], legFont["s=0.3"],  legFont["Near-Axis"],  legFont["Mixed"]}
  , LegendLayout -> {"Row", 3}
 ];

paperPlot=Grid[{
	{Style[stell, 20, FontFamily -> "Latin Modern Math"], SpanFromLeft},
	{rplotBmag,rplotgradAlpha2r2},
	{rplotgradPsi2r2,rplotgradAlphaPsi},
	{rplotBcrossGradBdotGradAlphaOB3,rplotBcrossGradBdotGradPsiOB3},
	{paperLeg,SpanFromLeft}
	}];

Export[toPaperFolder<>stell<>"_plot.pdf",paperPlot];
(*Export[toPaperFolder<>stell<>"_plot.mx" ,paperPlot, "AllowRasterization" -> True];*)
]

CloseKernels[];Exit[];
