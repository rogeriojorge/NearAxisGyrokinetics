#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[Directory[]];
pi = N[Pi,10];
mu0 = 4*pi*10^(-7)

etabar        = ToExpression[$ScriptCommandLine[[2]]];
B0            = ToExpression@$ScriptCommandLine[[3]];
VMECfileIn    = $ScriptCommandLine[[4]];
stell         = $ScriptCommandLine[[5]];
fileIn        = $ScriptCommandLine[[6]];
saveStellFile = $ScriptCommandLine[[7]];
gs2gridNA     = $ScriptCommandLine[[8]];
gxgridNA      = $ScriptCommandLine[[9]];

nperiod  = 1;
drhodpsi = 1;
rmaj     = 1;
kxfac    = 1;

If[FileExistsQ[saveStellFile],DumpGet[saveStellFile];,
	raxis = Import[VMECfileIn, {"Datasets", "raxis_cc"}];
	zaxis = Import[VMECfileIn, {"Datasets", "zaxis_cs"}];
	R0 = raxis[[1]]; eR = If[Length[raxis] == 1, 0, raxis[[2 ;; 6]]]; eZ = If[Length[zaxis] == 1, 0, zaxis[[2 ;; 6]]];
	NFP = Import[VMECfileIn, {"Datasets", "nfp"}];

	sigma0 = 0;
	ns = 856;
	nIts = 10;
	nLineS = 4;
	tol = 10^-6;
	sigmaIn = 0.1*Cos[NFP t];
	iota0 = 0.1;
	
	phiEDGE = Abs[Last[Import[VMECfileIn, {"Datasets", "phi"}]]];
	ctor = Import[VMECfileIn, {"Datasets", "ctor"}];
	I2oB0 = (mu0*ctor)/(2*phiEDGE);

	(************)

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
	
	(*Convert VMEC phi to axis phi*)
	phiAxis=Compile[{{phi,_Real}},(2*pi/Laxis)*GaussLegendreQuadrature[sprimeFunc[t],{t,0,phi},180,12]];
	phiA=Interpolation[Table[Evaluate@{t,phiAxis[t]},{t,-40*2.*pi,40*2.*pi,(40.*pi)/1500.}]]//Quiet;

	DumpSave[saveStellFile, {sigma,iota,curvFunc,torsFunc,sprimeFunc,etabar,B0,nNormal,Aminor,phiEDGE,sigmatemp,Laxis,phiA}];
];

(**********************)
(*Output to GS2 and GX*)
(**********************)

(*Normalization Parameters*)
Aminor = Abs[Import[VMECfileIn, {"Datasets", "Aminor_p"}]];
phiEDGE = Abs[Last[Import[VMECfileIn, {"Datasets", "phi"}]]/(2 pi)];

(*Grid Parameters*)
tgridVMEC    	   = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "tgrid("][[1]]]],StringSplit[{}[[1]]] // Quiet][[4 ;;]], {"E+" :> "*^", "E-" :> "*^-"}];
shat         	   = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "shat"][[1]]]], StringSplit[{}[[1]]] // Quiet], {"E+" :> "*^", "E-" :> "*^-"}][[2]];
normalizedtorFlux  = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "normalized_flux"][[1]]]], StringSplit[{}[[1]]] // Quiet], {"E+" :> "*^", "E-" :> "*^-"}][[2]];
alphaVMEC          = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "alpha"][[1]]]],StringSplit[{}[[1]]] // Quiet], {"E+" :> "*^", "E-" :> "*^-"}][[2]];
nlambda            = ToExpression@StringReplace[DeleteCases[Flatten[StringSplit[FindList[fileIn, "nlambda"][[1]]]],StringSplit[{}[[1]]] // Quiet], {"E+" :> "*^", "E-" :> "*^-"}][[2]];

(*Geometry Definitions*)

rVMEC 			    = -Sqrt[((2 phiEDGE normalizedtorFlux)/B0)];
Phi[theta_] 	    = (theta - alphaVMEC)/(iota - nNormal);
bmagNew[theta_]     = (Aminor^2 B0 (1 + rVMEC etabar Cos[theta]))/(2 phiEDGE);
gradparNew[theta_]  = Sign[Last[Import[VMECfileIn, {"Datasets", "phi"}]]]*((2 Aminor pi (1 + rVMEC etabar Cos[theta]))/Laxis)/(sprimeFunc[(alphaVMEC-theta)/(iota-nNormal)]*2*pi/Laxis);
gds2New[theta_]     = ((Aminor^2) B0 )/(2 phiEDGE) ((etabar^2 Cos[theta]^2)/curvFunc[Phi[theta]]^2 + (curvFunc[Phi[theta]]^2 (Sin[theta] + Cos[theta] sigma[Phi[theta]])^2)/etabar^2);
gds21New[theta_]    = -Sign[Last[Import[VMECfileIn, {"Datasets", "phi"}]]]/(2 phiEDGE) Aminor^2  shat ((B0 etabar^2 Cos[theta] Sin[theta])/curvFunc[Phi[theta]]^2 + 1/etabar^2 B0 curvFunc[Phi[theta]]^2 (Sin[theta] + Cos[theta] sigma[Phi[theta]]) (-Cos[theta] + Sin[theta] sigma[Phi[theta]]));
gds22New[theta_]    = (Aminor^2 B0 shat^2 (etabar^4 Sin[theta]^2 + curvFunc[Phi[theta]]^4 (Cos[theta] - Sin[theta] sigma[Phi[theta]])^2))/(2 phiEDGE etabar^2 curvFunc[Phi[theta]]^2);
gbdriftNew[theta_]  = Sign[Last[Import[VMECfileIn, {"Datasets", "phi"}]]]*(2 Sqrt[2] etabar Cos[theta])/Sqrt[B0/phiEDGE] (1 - 0 2 rVMEC etabar Cos[theta]); 
cvdriftNew[theta_]  = gbdriftNew[theta];
gbdrift0New[theta_] = -2 Sqrt[2] Sqrt[phiEDGE/B0] shat etabar Sin[theta] (1 - 0 2 rVMEC etabar Cos[theta]); 
cvdrift0New[theta_] = gbdrift0New[theta];
lambdamin=((2 phiEDGE)/(Aminor^2 B0))/(1+Abs[rVMEC*etabar]);
lambdamax=((2 phiEDGE)/(Aminor^2 B0))/(1-Abs[rVMEC*etabar]);
lambda=Table[i,{i,lambdamin,lambdamax,(lambdamax-lambdamin)/(nlambda-1)}];

(*Output to GS2 grid*)
ToStringA[x_]:=ToString[FortranForm@N[1. Chop[x,10^-16],10]];
nz=Length[tgridVMEC];
ntgrid=Floor[nz/2];
paramTheta[i_]:=tgridVMEC[[i]]*(iota-nNormal);
paramZ[z_]:=tgridVMEC[[i]];
text="nlambda\n"<>ToString[nlambda]<>"\nlambda";
For[i=1,i<=nlambda,i++,text=text<>"\n"<>ToString[lambda[[i]]]];
text=text<>"\nntgrid nperiod ntheta drhodpsi rmaj shat kxfac q";
text=text<>"\n"<>ToString@ntgrid<>" "<>ToString@nperiod<>" "<>ToString[nz-1]<>" "<>ToString@drhodpsi<>" "<>ToString@rmaj<>" "<>ToString@shat<>" "<>ToString@kxfac<>" "<>ToStringA[1/iota];
text=text<>"\ngbdrift gradpar grho tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[gbdriftNew[zz]]<>" "<>ToStringA[gradparNew[zz]]<>" "<>ToStringA[1.0]<>" "<>ToStringA[paramZ[i]]];
text=text<>"\ncvdrift gds2 bmag tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[cvdriftNew[zz]]<>" "<>ToStringA[gds2New[zz]]<>" "<>ToStringA[bmagNew[zz]]<>" "<>ToStringA[paramZ[i]]];
text=text<>"\ngds21 gds22 tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[gds21New[zz]]<>" "<>ToStringA[gds22New[zz]]<>" "<>ToStringA[paramZ[i]]];
text=text<>"\ncvdrift0 gbdrift0 tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[cvdrift0New[zz]]<>" "<>ToStringA[gbdrift0New[zz]]<>" "<>ToStringA[paramZ[i]]];
text=text<>"\nRplot Rprime tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[0.0]<>" "<>ToStringA[0.0]<>" "<>ToStringA[paramZ[i]]];
text=text<>"\nZplot Zprime tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[0.0]<>" "<>ToStringA[0.0]<>" "<>ToStringA[paramZ[i]]];
text=text<>"\naplot aprime tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[0.0]<>" "<>ToStringA[0.0]<>" "<>ToStringA[paramZ[i]]];
Export[gs2gridNA,text,"String",CharacterEncoding->"ASCII"];

(*Output to GX grid*)
text="ntgrid nperiod ntheta drhodpsi rmaj shat kxfac q scale";
text=text<>"\n"<>ToString@ntgrid<>" "<>ToString@nperiod<>" "<>ToString[nz-1]<>" "<>ToString@drhodpsi<>" "<>ToString@rmaj<>" "<>ToString@shat<>" "<>ToString@kxfac<>" "<>ToStringA[1/iota]<>" 1.0";
text=text<>"\ngbdrift gradpar grho tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[gbdriftNew[zz]]<>" "<>ToStringA[gradparNew[zz]]<>" "<>ToStringA[1]<>" "<>ToStringA[paramZ[i]]];
text=text<>"\ncvdrift gds2 bmag tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[cvdriftNew[zz]]<>" "<>ToStringA[gds2New[zz]]<>" "<>ToStringA[bmagNew[zz]]<>" "<>ToStringA[paramZ[i]]];
text=text<>"\ngds21 gds22 tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[gds21New[zz]]<>" "<>ToStringA[gds22New[zz]]<>" "<>ToStringA[paramZ[i]]];
text=text<>"\ncvdrift0 gbdrift0 tgrid";
For[i=1,i<=nz,i++,zz=paramTheta[i];text=text<>"\n"<>ToStringA[cvdrift0New[zz]]<>" "<>ToStringA[gbdrift0New[zz]]<>" "<>ToStringA[paramZ[i]]];
Export[gxgridNA,text,"String",CharacterEncoding->"ASCII"];