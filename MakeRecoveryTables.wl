(* ::Package:: *)

SetDirectory[NotebookDirectory[]];
Import["TwoDimBumpWave.wl"]


(*Set parameters and kernel definitions.*)
NumSep = 400; 
GridWidth = 10.;
GridSep = GridWidth/NumSep;
epsG = 10^(-10);
numGridIntervals = 20;
DeltaList = Range[3.75,4.60,0.05];
GridSizes = {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.89};
KG[x_,y_] = Exp[-(x^2+y^2)/2];


(*Import sub-segment to cell distances. Distances are computed in script segmentdistances.m with MATLAB CVX package.*)
dists = importMatrix["partition_dists_20191108_npartitions_100.mat"];


(*Import previously built envelopes.*)
FLS = Import["~/Data/Supremums_20190531.mat"];
FLM = Import["~/Data/Monotones_20190831.mat"];

(*To build monotone and supremum envelopes uncomment and run following cells.*)
(*FLG = compileTables[];*)
(*FLM = compileFunListMonotone[GridSep,NumSep,FLG,epsG];*)
(*FLS = compileFunListSupremum[GridSep,NumSep,FLG];*)


(*Compute norms, coefficients, and recovery success.*)
TableF[gridIndex_,fIndex_,d_] := FLM[[gridIndex,fIndex,getIndex[d,FLM[[gridIndex,15]]]]]
matnorms = makeMatrixNorms[FLM,DeltaList,epsG];
norms = makeNormTables[FLM,DeltaList,epsG];
coeffs = makeAlphaBetaTables[FLM,DeltaList,epsG];
makeQTables[FLM,FLS,DeltaList,dists,epsG];
makeExtendedTables[FLM,FLS,DeltaList,dists,epsG];


(*Data to show that Q is bounded for points beyond Delta.*)
makeQTablesForDistantPoints[FLM,DeltaList,epsG];
