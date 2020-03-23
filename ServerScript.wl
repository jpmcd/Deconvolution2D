(* ::Package:: *)

(*SetDirectory[NotebookDirectory[]];*)

Import["TwoDimBumpWave.wl"]
$ProcessorCount

(*Set parameters and kernel function definitions*)
numGridPoints = 400; 
GridWidth = 10.;
GridSpacing = GridWidth/numGridPoints;
numGridIntervals = 20;
gridSizes = {.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7};
epsG = 10^(-10);
tailG = 10;
KG[x_,y_] = Exp[-(x^2+y^2)/2];
LG[K_,x_,y_] = (x^2+y^2-1)*K[x,y];

(*Compute function tables in parallel*)
Print[AbsoluteTiming[funListGauss = makeFunctions[GridSpacing,numGridPoints,numGridIntervals,gridSizes,KG,LG,epsG,tailG];]];

(*Correct for Interval Arithmetic of gradient at origin*)
temp = fixGradientDivisionByZero[funListGauss,4,5,8];
funListGauss[[4]] = temp;

(*Save tables*)
Do[Export[ToString[StringForm["./flg_20190401_``.mat",i]],funListGauss[[All,All,All,i]],"MAT"],{i,Dimensions[funListGauss][[4]]}];
Quit[]
