(* ::Package:: *)

BeginPackage["TwoDimBumpWave`"]
EPS::usage="10^(-5) constant used in calculations."

computeBumpWave::usage="Construct functions to algebraically compute and simplify bump, waves and derivatives, and Hessian's largest eigenvalues centered at the origin."
computeEigenvalues::usage="Construct eigenvalue functions. Merged into computeBumpWave."

makeSampleTable::usage="Builds rows of a function table, fixing the function evaluation position, grid spacing and varying the position of the spike between its three nearest samples. The maximum function value over the spike positions is returned."
makeFunctionTable::usage="Builds a three dimensional table array to store upper bounds of bumps, waves and derivatives across a square partition of [-10,10]x[-10,10] and different choices for grid spacing. makeSampleTable is called for each entry of the table."
makeSampleTableSigned::usage="Similar to makeSampleTable but extends range to negative coordinates."
makeFunctionTableSigned::usage="Similar to makeFunctionTable but calling makeSampleTableSigned, extending range to negative coordinates."

makeMonotone::usage="Given (tab,eps) where tab is a function table and eps is a number, the resulting table arr is monotonized so that arr[[i,j]] >= max(arr[[k,j]],eps) for any i<=k.  Used to obtain monotonic upper bounds on non-negative functions."
makeMaxDistanceTable::usage="Determine max distances from corners of squares from a two-dimensional square partition."
makeMinDistanceTable::usage="Determine min distances from corners of squares from a two-dimensional square partition."

makeFun::usage="Calls makeFunctionTable."
makeAbsFun::usage="Calls makeFunctionTable using abs(f) and then calls makeMonotone on the result using eps."
makeSignFun::usage="Similar to makeAbsFun but calls makeFunctionTableSigned."

makeDistFunTables::usage="Wrapper to create monotonized function arrays. The last array provides a distance argument."
makePartition::usage="Take a square partition of [-gridSep*numSep]x[-gridSep*numSep] in the plane and construct a flattened array of the sorted, nearest and furthest distances of each square from the origin."

makeFunctionsBW::usage="Wrapper to construct list of function tables for bumps, waves using makeFun and makeAbsFun.  
Expects bwlist to be of the form output by computeBumpWave.  Returns the following list of functions where d denotes differentiation and A denotes absolute-monotonized: (AB,AW,AS,DB,ADB,ADW,ADS,DDB,ADDB,ADDW,ADDS)"
makeFunctions::usage="Given (gammaWidth,DeltaWidth,numPerGamma,numDelta,numPerDelta,kappa,K,eps,tail,factor), calls makeFunctionsBW with computeBumpWave[K]"
makeEigenvalues::usage="Similar to makeFunctions but for eigenvalue functions."
fixGradientDivisionByZero::usage="Gradient envelope calculation at the origin for bump functions."

getNearbySpikes::usage="Lists the coordinates of the nearest points in the nearest cells of the first sector of the grid. There are six sectors by hexagonal symmetry."
spikeDistances::usage="Wrapper to get norms of points returned by getNearbySpikes."

getIndex::usage="Lists the coordinates of the nearest points in the nearest cells of the first sector of the grid. There are six sectors by hexagonal symmetry."
getOverlapping::usage="Get the indices of a partition that overlap with an interval."
buildSupremum::usage="Build a supremum for a function table f defined on a partition."

computeNormMat::usage="Compute matrix norms by summing functions using distances from spikes (worst-case distances). Infinity norm bounds are returned in block matrix form of identity minus bump-wave interpolation matrix, {{I-B,W1,W2},{D1B,I-D1W1,D1W2},{D2B,D2W1,I-D2W2}}."
computeAlphaBeta::usage="Compute bump/wave coefficients and bump coefficient lower bound."
computeBounds::usage="Computes the bounds on Q, DQ, DDQ over the interval [0,Delta/2]. Assumes a single length or range for grid spacing. Thus, this should be computed for each choice of grid spacing.
Spike at the origin is assumed. Returns $Failed if computeAlphaBeta fails."

DDQIntegral::usage="Calculate an integral across subintervals of [0,Delta/2] to determine where the following expression is negative:
  \\int_{0}^{t}v^TH(x_0+sv)v (t-s) ds."
DQIntegral::usage="Calculate an integral across subintervals of [0,Delta/2] to determine where the following expression is negative:
  \\int_{u_1}^{t}\\nabla Q(x_0+sv)\\cdot v ds."

makeQTables::usage="Uses integral method to determine if Q is bounded. Output of computeBounds DQ, and DDQ are integrated over segments spanning distance 0 to Delta/2 from spike at the origin.
Boundedness of Q is determined over each segment. Tables for recovery and lengths for boundedness for DDQ and DQ are saved."
makeQTablesForDistantPoints::usage="Determines bound on Q for points farther than Delta from each spike."
makeExtendedTables::usage=""

testQRegionBound::usage="OLD: testQRegionBound doesn't use integral method and is less accurate. Tests whether the dual combination Q is properly bounded. Determines if the necessary conditions properly overlap. Guarantees that Q is strongly concave near the spike at origin with sign +1 and |Q|<1 within Delta/2 from the origin."

importMatrix::usage="File I/O method for opening data files."
getPartitionDistancesToSpikes::usage="Helper function."
importAndFlatten::usage="Helper function."
getNearFunctionTables::usage="Imports bump Hessian eigenvalue and directional derivative envelope arrays."
testSuccess::usage="OLD: makeQTables makes DAT files for recovery results. testQRegionBound doesn't use integral method and is less accurate. Check for successful recovery."
makeRecoveryTables::usage="OLD: testQRegionBound doesn't use integral method and is less accurate. Organizes recovery results as table for different minimum separation and grid spacing choices."

compileTables::usage="Import unflattened (i.e. 4 dimensional) function tables."
compileTablesEV::usage="Import unflattened (i.e. 4 dimensional) function tables for eigenvalues."
compileFunListMonotone::usage="Create flattened, monotonized function arrays from an unflattened table."
compileFunListSupremum::usage="Create flattened, non-monotonized supremum tables for gradient and eigenvalue functions from function table."

plotSupremum::usage="Plot non-monotonized supremum envelopes. Note buildSupremum can take several minutes."
plotSupremumFromFile::usage="Plot non-monotonized supremum envelopes that are already built and imported."
plotMonotone::usage="Plot monotonized function envelopes."
plotMonotoneFromFile::usage="Plot envelopes that are already monotized and imported."
plotBoundList::usage="Using computeBounds output, plot dual combination and derivatives Q, DQ and DDQ over [0,Delta/2]."
makeMatrixNorms::usage="Organizes block interpolation matrix norms as table for different minimum separation and grid spacing choices to be saved for external plotting."
makeNormTables::usage="Organizes Schur complement norms as table for different minimum separation and grid spacing choices to be saved for external plotting."
makeAlphaBetaTables::usage="Organizes coefficient norms and bump coefficient minimum as table for different minimum separation and grid spacing choices to be saved for external plotting."
lowerBoundSuccess::usage="Checks lower bound (Q > -1) for a list of spike separations as needed for recovery."

plotTable::usage="Plots tables of booleans with given tick lists. Ticks are rounded to the grid."
plotValueTable::usage="Plots tables of values with given tick lists. Ticks are rounded to the grid."


Begin["`Private`"]


(*Constant*)
EPS = 10^(-5);


P[x_,y_] = (x^2+y^2-1);


(*Construct functions to algebraically compute and simplify bump,
  waves and derivatives centered at the origin. The three sample
  points that make up each bump/wave are:
  s_1 = (-x,-y)
  s_2 = (s-x,-y)
  s_3 = (-x,s-y)

  K = kernel function
  D1K = first partial of K
  D2K = second partial of K
  L = maximum eigenvalue of Hessian of K

  For each function definition, the arguments are:
  tx = x-coordinate of location of function evaluation
  ty = y-coordinate of location of function evaluation
  x = x-offset of origin-centered spike from lower left sample
  y = y-offset of spike from lower left sample
  s = sample gridpoint separation distance
*)
computeBumpWave[K_,D1K_,D2K_,L_] := 
	Module[{Minv,B,W1,W2,d1B,d1W1,d1W2,d2B,dBdotN,d2W1,d2W2,LB,LBA,LW1,LW2},
		Minv[x_,y_,s_] = Inverse[{{K[-x,-y],K[s-x,-y],K[-x,s-y]},{-D1K[-x,-y],-D1K[s-x,-y],-D1K[-x,s-y]},{-D2K[-x,-y],-D2K[s-x,-y],-D2K[-x,s-y]}}];
		Print["Simplifying Bump"];
		B[tx_,ty_,x_,y_,s_] = FullSimplify[({{K[-x-tx,-y-ty],K[s-x-tx,-y-ty],K[-x-tx,s-y-ty]}}.Minv[x,y,s].{{1},{0},{0}})[[1,1]]];
		Print["Simplifying Wave1"]; (*Subbing 0 in third coordinate by inverse coefficient matrix formula and choices for s_1, s_2, s_3.*)
		W1[tx_,ty_,x_,y_,s_] = FullSimplify[({{K[-x-tx,-y-ty],K[s-x-tx,-y-ty],0}}.Minv[x,y,s].{{0},{1},{0}})[[1,1]]];
		Print["Simplifying Wave2"]; (*0 in second coordinate for the same reason.*)
		W2[tx_,ty_,x_,y_,s_] = FullSimplify[({{K[-x-tx,-y-ty],0,K[-x-tx,s-y-ty]}}.Minv[x,y,s].{{0},{0},{1}})[[1,1]]];
		Print["Simplifying DxBump"];
		d1B[tx_,ty_,x_,y_,s_] = FullSimplify[D[B[tx,ty,x,y,s],tx]];
		Print["Simplifying DxWave1"];
		d1W1[tx_,ty_,x_,y_,s_] = FullSimplify[D[W1[tx,ty,x,y,s],tx]];
		Print["Simplifying DxWave2"];
		d1W2[tx_,ty_,x_,y_,s_] = FullSimplify[D[W2[tx,ty,x,y,s],tx]];
		Print["Simplifying DyBump"];
		d2B[tx_,ty_,x_,y_,s_] = FullSimplify[D[B[tx,ty,x,y,s],ty]];
		Print["Simplifying DyWave1"];
		d2W1[tx_,ty_,x_,y_,s_] = FullSimplify[D[W1[tx,ty,x,y,s],ty]];
		Print["Simplifying DyWave2"];
		d2W2[tx_,ty_,x_,y_,s_] = FullSimplify[D[W2[tx,ty,x,y,s],ty]];
		Print["Simplifying DBumpDotN"];
		dBdotN[tx_,ty_,x_,y_,s_] = FullSimplify[(tx*d1B[tx,ty,x,y,s]+ty*d2B[tx,ty,x,y,s])/Sqrt[tx^2+ty^2]];(*Note this causes INF at origin, must be fixed. Handled by fixGradientDivisionByZero below.*)
		Print["Simplifying Bump Eigenvalue"];
		LB[tx_,ty_,x_,y_,s_] = FullSimplify[({{L[K,-x-tx,-y-ty],L[K,s-x-tx,-y-ty],L[K,-x-tx,s-y-ty]}}.Minv[x,y,s].{{1},{0},{0}})[[1,1]]];
		Print["Simplifying Absolute Bump Eigenvalue"];
		LBA[tx_,ty_,x_,y_,s_] = FullSimplify[{K[-x-tx,-y-ty]*Max[P[-x-tx,-y-ty],1],K[s-x-tx,-y-ty]*Max[P[s-x-tx,-y-ty],1],K[-x-tx,s-y-ty]*Max[P[-x-tx,s-y-ty],1]}.Minv[x,y,s][[;;,1]]];
		Print["Simplifying Wave1 Eigenvalue"];
		LW1[tx_,ty_,x_,y_,s_] = FullSimplify[K[-x-tx,-y-ty]*Max[P[-x-tx,-y-ty],1]*Abs[Minv[x,y,s][[1,2]]]+K[s-x-tx,-y-ty]*Max[P[s-x-tx,-y-ty],1]*Abs[Minv[x,y,s][[2,2]]]];
		Print["Simplifying Wave2 Eigenvalue"];
		LW2[tx_,ty_,x_,y_,s_] = FullSimplify[K[-x-tx,-y-ty]*Max[P[-x-tx,-y-ty],1]*Abs[Minv[x,y,s][[1,3]]]+K[-x-tx,s-y-ty]*Max[P[-x-tx,s-y-ty],1]*Abs[Minv[x,y,s][[3,3]]]];
		{B,W1,W2,d1B,d1W1,d1W2,d2B,d2W1,d2W2,dBdotN,LB,LBA,LW1,LW2}
	]
computeBumpWave[K_,L_] := 
	Module[{DX,DY},
		DX[x_,y_] = D[K[x,y],x];
		DY[x_,y_] = D[K[x,y],y];
		computeBumpWave[K,DX,DY,L]
	]
computeEigenvalues[K_,D1K_,D2K_,L_] :=
	Module[{Minv,LBA,LW1,LW2},
		Minv[x_,y_,s_] = Inverse[{{K[-x,-y],K[s-x,-y],K[-x,s-y]},{-D1K[-x,-y],-D1K[s-x,-y],-D1K[-x,s-y]},{-D2K[-x,-y],-D2K[s-x,-y],-D2K[-x,s-y]}}];
		Print[DateString[]];
		Print["Simplifying Absolute Bump Eigenvalue"];
		LBA[tx_,ty_,x_,y_,s_] = FullSimplify[{K[-x-tx,-y-ty]*Max[P[-x-tx,-y-ty],1],K[s-x-tx,-y-ty]*Max[P[s-x-tx,-y-ty],1],K[-x-tx,s-y-ty]*Max[P[-x-tx,s-y-ty],1]}.Minv[x,y,s][[;;,1]]];
		Print[DateString[]];
		Print["Simplifying Wave1 Eigenvalue"];
		LW1[tx_,ty_,x_,y_,s_] = FullSimplify[K[-x-tx,-y-ty]*Max[P[-x-tx,-y-ty],1]*Abs[Minv[x,y,s][[1,2]]]+K[s-x-tx,-y-ty]*Max[P[s-x-tx,-y-ty],1]*Abs[Minv[x,y,s][[2,2]]]];
		Print[DateString[]];
		Print["Simplifying Wave2 Eigenvalue"];
		LW2[tx_,ty_,x_,y_,s_] = FullSimplify[K[-x-tx,-y-ty]*Max[P[-x-tx,-y-ty],1]*Abs[Minv[x,y,s][[1,3]]]+K[-x-tx,s-y-ty]*Max[P[-x-tx,s-y-ty],1]*Abs[Minv[x,y,s][[3,3]]]];
		Print[DateString[]];
		{LBA,LW1,LW2}
	]
computeEigenvalues[K_,L_] := 
	Module[{DX,DY},
		DX[x_,y_] = D[K[x,y],x];
		DY[x_,y_] = D[K[x,y],y];
		computeEigenvalues[K,DX,DY,L]
	]


(*Sub-routine to build rows of function tables.
  The inner loop iterates through the range of
  offsets (x,y) of the origin-centered spike from s_1,
  evaluating the function for every position of a
  spike between its three nearest samples.
  The maximum of this array is returned, since
  a max over the possible sample points for the
  chosen function is desired.

  tx = x-coordinate for function evaluation
  ty = y-coordinate for function evaluation
  numGridIntervals = number of intervals to partition location of spike offset (with respect to samples s1, s2, s3)
  sLow = minimum distance (in sigma) for intervals between sample grid separation
  sHigh = maximum distance (in sigma) for intervals between sample grid separation
  f = function (bump, wave, or derivatives)
  
  kLeft, kRight, kLow and kHigh are defined in terms of sHigh, where kLeft's lower limit is 0 and kRight's higher limit is sHigh/2. Similar for kLow and kHigh.
*)
makeSampleTable[tx_,ty_,numGridIntervals_,sLow_,sHigh_,f_] :=
	Module[{ret,w,v,kLeft,kRight,kLow,kHigh},
		w = (sHigh/2)/numGridIntervals;
		ret = Table[-Infinity,{k,1,numGridIntervals},{l,1,numGridIntervals}]; (*Result array*)
		Do[
			{kLeft,kRight} = {(k-1)*w,k*w}; (*Select interval of x-coordinate of spike offset*)
			{kLow,kHigh} = {(l-1)*w,l*w}; (*Select interval of y-coordinate of spike offset*)
			ret[[k,l]] = Max[f[tx,ty,Interval[{kLeft,kRight}],Interval[{kLow,kHigh}],Interval[{sLow,sHigh}]]]; (*Max of f over spike offset and grid spacing intervals*)
			,{k,1,numGridIntervals}
			,{l,1,numGridIntervals}
		];
		Max[ret] (*Maximize since envelope is upper bound on function over choices of spike location*)
	]
makeSampleTableSigned[tx_,ty_,numGridIntervals_,sLow_,sHigh_,f_] :=
	Module[{ret,w,v,kLeft,kRight,kLow,kHigh,m1,m2},
		w = (sHigh/2)/numGridIntervals;
		ret = Table[-Infinity,{k,1,numGridIntervals},{l,1,numGridIntervals}]; (*Result array*)
		Do[
			{kLeft,kRight} = {(k-1)*w,k*w}; (*Select interval of x-coordinate of spike offset*)
			{kLow,kHigh} = {(l-1)*w,l*w}; (*Select interval of y-coordinate of spike offset*)
			m1 = Max[f[tx,ty,Interval[{kLeft,kRight}],Interval[{kLow,kHigh}],Interval[{sLow,sHigh}]]];
			m2 = Max[f[tx,ty,-Interval[{kLeft,kRight}],-Interval[{kLow,kHigh}],-Interval[{sLow,sHigh}]]]; (*Same function but with spike in third quadrant, *)
			ret[[k,l]] = Max[m1,m2]; (*Max of f over spike offset and grid spacing intervals*)
			,{k,1,numGridIntervals}
			,{l,1,numGridIntervals}
		];
		Max[ret] (*Maximize since envelope is upper bound on function over choices of spike location*)
	]


(*Routine to build the function tables. Inner loop
  fixes the location of the function's argument;
  i.e. position where the function is evaluated away from spikes and samples,
  The loop also selects the different sample grid spacing ranges.
  The range of the function argument is [-gridSep*numSep,gridSep*numSep] in both dimensions,
  broken up into squares of side length gridSep.
  
  gridSep = separation distance for the function grid (not grid spacing)
  numSep = number of function grid points in each direction
  numGridIntervals = number of intervals to partition grid spacing for spike offset location
  gridSizes = grid spacing, i.e. |s1-s2| and |s1-s3|
  f = function (bump, wave, or derivatives)
  eps = bound on tail function height
  tail = distance beyond which function is marginal than eps
*)
makeFunctionTable[gridSep_,numSep_,numGridIntervals_,gridSizes_,f_,eps_,tail_] :=
	Module[{w},
		w = gridSep;
		Print["Making Function Table"];
		Print[DateString[]];
		ParallelTable[
			Module[{DLeft,DRight,DLow,DHigh},
				{DLeft,DRight} = {w*tx-EPS,w*(tx+1)+EPS};
				{DLow,DHigh} = {w*ty-EPS,w*(ty+1)+EPS};
				If[Min[DLeft^2,DRight^2]+Min[DLow^2,DHigh^2]>=tail^2,
					eps,
					makeSampleTable[Interval[{DLeft,DRight}],Interval[{DLow,DHigh}],numGridIntervals,gridSizes[[s]],gridSizes[[s+1]],f]
				]
			]
			,{tx,-numSep,numSep}(*Selects horizontal interval on function grid*)
			,{ty,-numSep,numSep}(*Select vertical interval on function grid*)
			,{s,Length[gridSizes]-1}
			,DistributedContexts -> Automatic
		]
	]
makeFunctionTableSigned[gridSep_,numSep_,numGridIntervals_,gridSizes_,f_,eps_,tail_] :=
	Module[{w},
		w = gridSep;
		Print["Making Function Table"];
		Print[DateString[]];
		ParallelTable[
			Module[{DLeft,DRight,DLow,DHigh},
				{DLeft,DRight} = {w*tx-EPS,w*(tx+1)+EPS};
				{DLow,DHigh} = {w*ty-EPS,w*(ty+1)+EPS};
				If[Min[DLeft^2,DRight^2]+Min[DLow^2,DHigh^2]>=tail^2,
					eps,
					makeSampleTableSigned[Interval[{DLeft,DRight}],Interval[{DLow,DHigh}],numGridIntervals,gridSizes[[s]],gridSizes[[s+1]],f]
				]
			]
			,{tx,-numSep,numSep}(*Selects horizontal interval on function grid*)
			,{ty,-numSep,numSep}(*Select vertical interval on function grid*)
			,{s,Length[gridSizes]-1}
			,DistributedContexts -> Automatic
		]
	]


(*Routines to monotonize a function table.*)
makeMonotone[tab_,eps_] := Reverse[FoldList[Max[#1,#2,eps]&,Reverse[tab]]](*Make a monotonically decreasing table from a table*)
(*Determine min/max distances from corners of boxes from a two-dimensional grid.

  gridSep = function grid separation distance
  numSep = number of function grid points in each direction 
  F = function handle
*)
makeFDistanceTable[gridSep_,numSep_,F_] :=
	Module[{w},
		w = gridSep;
		Table[
			Module[{DLeft, DRight, DLow, DHigh},
				{DLeft, DRight} = {w*tx-EPS, w*(tx+1)+EPS};
				{DLow, DHigh} = {w*ty-EPS, w*(ty+1)+EPS};
				Sqrt[F[DLeft^2, DRight^2]+F[DLow^2, DHigh^2]]
			]
			,{tx,-numSep,numSep}
			,{ty,-numSep,numSep}
		]
	]
makeMaxDistanceTable[gridSep_,numSep_] := makeFDistanceTable[gridSep,numSep,Max]
makeMinDistanceTable[gridSep_,numSep_] := makeFDistanceTable[gridSep,numSep,Min]
(*Wrappers to build regular, or absolute monotonized function tables*)
makeFun[gridSep_,numSep_,numGridIntervals_,gridSizes_,f_,eps_,tail_] := makeFunctionTable[gridSep,numSep,numGridIntervals,gridSizes,f,eps,tail] 
makeAbsFun[gridSep_,numSep_,numGridIntervals_,gridSizes_,f_,eps_,tail_] := makeFunctionTable[gridSep,numSep,numGridIntervals,gridSizes,Abs[f[#1,#2,#3,#4,#5]]&,eps,tail]
makeSignFun[gridSep_,numSep_,numGridIntervals_,gridSizes_,f_,eps_,tail_] := makeFunctionTableSigned[gridSep,numSep,numGridIntervals,gridSizes,Abs[f[#1,#2,#3,#4,#5]]&,eps,tail]
(*Create monotonized function arrays for each table in funTables.
  This builds a one dimensional, monotonized array of bump/wave
  envelopes from two dimensional function tables from
  makeFunctionTable and a flattened distance table from
  makeMax/MinDistanceTable.
  
  gridSep = function grid separation distance
  numSep = number of function grid points in each direction
  funTables = table of two dimensional arrays of function values over grid
  eps = lower/tail value for function height
*)
makeDistFunTables[gridSep_,numSep_,funTables_,eps_] := 
	Module[{distTable,flatDistTable,ordering,sortedDistTable,sortedFunTables,monotoneFunTables},
		Print["Making Tables..."];
		distTable = makeMaxDistanceTable[gridSep,numSep];
		flatDistTable = Flatten[distTable];
		ordering = Ordering[flatDistTable];
		sortedDistTable = flatDistTable[[ordering]];
		sortedFunTables = Table[Flatten[tab][[ordering]],{tab,funTables}];
		monotoneFunTables = Table[makeMonotone[tab,eps],{tab,sortedFunTables}];
		Append[monotoneFunTables,sortedDistTable] (*The first 14 tables list function values while the last gives the distance for each entry.*)
	]
(*Take a square partition of [-gridSep*numSep]x[-gridSep*numSep]
  in the plane and construct a flattened array of the sorted,
  nearest and furthest distances of each square from the origin.

  gridSep = function grid separation distance
  numSep = number of function grid points in each direction 
*)
makePartition[gridSep_,numSep_] :=
	Module[{maxes,mins},
		maxes = makeMaxDistanceTable[gridSep,numSep];
		mins = makeMinDistanceTable[gridSep,numSep];
		Sort[DeleteDuplicates[Flatten[{mins,maxes}]]]
	]


(*Constructs list of function tables for the bumps, waves, and derivatives. bwlist should be of form output by computeBumpWave.*)
makeFunctionsBW[gridSep_,numSep_,numGridIntervals_,gridSizes_,bwlist_,eps_,tail_] :=
	Module[{B,W1,W2,d1B,d1W1,d1W2,d2B,d2W1,d2W2,dBdotN,LB,LBA,LW1,LW2},
		{B,W1,W2,d1B,d1W1,d1W2,d2B,d2W1,d2W2,dBdotN,LB,LBA,LW1,LW2} = bwlist;
		{
			makeFun[gridSep,numSep,numGridIntervals,gridSizes,B,eps,tail],
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,W1,eps,tail],
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,W2,eps,tail],
			makeFun[gridSep,numSep,numGridIntervals,gridSizes,dBdotN,eps,tail],
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,d1B,eps,tail],
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,d1W1,eps,tail],
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,d1W2,eps,tail],
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,d2B,eps,tail],
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,d2W1,eps,tail],
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,d2W2,eps,tail],
			makeFun[gridSep,numSep,numGridIntervals,gridSizes,LB,eps,tail],
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,LBA,eps,tail],
			makeSignFun[gridSep,numSep,numGridIntervals,gridSizes,LW1,eps,tail],
			makeSignFun[gridSep,numSep,numGridIntervals,gridSizes,LW2,eps,tail]
		}
	]
makeFunctions[gridSep_,numSep_,numGridIntervals_,gridSizes_,K_,L_,eps_,tail_] := makeFunctionsBW[gridSep,numSep,numGridIntervals,gridSizes,computeBumpWave[K,L],eps,tail]
makeEigenvalues[gridSep_,numSep_,numGridIntervals_,gridSizes_,K_,L_,eps_,tail_] :=
	Module[{LBA,LW1,LW2},
		{LBA,LW1,LW2} = computeEigenvalues[K,L];
		{
			makeAbsFun[gridSep,numSep,numGridIntervals,gridSizes,LBA,eps,tail],
			makeSignFun[gridSep,numSep,numGridIntervals,gridSizes,LW1,eps,tail],
			makeSignFun[gridSep,numSep,numGridIntervals,gridSizes,LW2,eps,tail]
		}
	]
(*Gradient envelope calculation at the origin for bump functions.*)
fixGradientDivisionByZero[funTable_,dBi_,d1Bi_,d2Bi_] :=
	Module[{d1B,d2B,dimensions},
		d1B = funTable[[d1Bi]];
		d2B = funTable[[d2Bi]];
		dimensions = Dimensions[funTable];
		Table[
			If[Abs[funTable[[dBi,i,j,k]]]==Infinity,Sqrt[d1B[[i,j,k]]^2+d2B[[i,j,k]]^2],funTable[[dBi,i,j,k]]],
			{i,dimensions[[2]]},
			{j,dimensions[[3]]},
			{k,dimensions[[4]]}
		]
	]


(*Lists the coordinates of the nearest points in the nearest cells of the
  first sector of the grid. There are six sectors by hexagonal symmetry.*)
getNearbySpikes[] :=
	Module[{s},
		s = Sqrt[3];
		{{1,0},
		{3*s/4,0},{s/2,1/2},
		{5*s/4,0},{s,1/2},{3*s/4,5/4},
		{7*s/4,0},{3*s/2,1/2},{5*s/4,5/4},{s,2},
		{9*s/4,0},{2*s,1/2},{7*s/4,5/4},{3*s/2,2},{5*s/4,11/4},
		{11*s/4,0},{5*s/2,1/2},{9*s/4,5/4},{2*s,2},{7*s/4,11/4},{3*s/2,7/2},
		{13*s/4,0},{3*s,1/2},{11*s/4,5/4},{5*s/2,2},{9*s/4,11/4},{2*s,7/2},{7*s/4,17/4},
		{15*s/4,0},{7*s/2,1/2},{13*s/4,5/4},{3*s,2},{11*s/4,11/4},{5*s/2,7/2},{9*s/4,17/4},{2*s,10/2}}
	]
spikeDistances[] :=
	Table[Norm[x],{x,getNearbySpikes[]}]


(*Routines to compute infinity norm upper bounds on the entries of the
  bump-wave interpolation matrix.
*)
(*For a sorted table get the index for the largest entry
  in the table less than or equal the value d.

  d = specified value
  table = sorted table
*)
getIndex[d_,table_] :=
	Module[{start,mid,end},
		start=1;
		end=Length[table];
		If[d<=table[[start]],Return[start],Null];
		If[table[[end]]<=d,Return[end],Null];
		While[end-start>1,
			mid = Floor[(start+end)/2];
			If[table[[mid]]<=d,start=mid,end=mid];
		];
		start
	]
(*Get the indices of a partition that overlap with an interval.

  left = left endpoint of interval
  right = right endpoint of interval
  partition = sorted table of partition values 
*)
getOverlapping[left_,right_,partition_] :=
	Module[{start,end},
		start = getIndex[left,partition];
		(*NOTE: it seems like the +1 can be dropped from the next line.*)
		end = Min[getIndex[right,partition]+1,Length[partition]];
		{start,end}
	]
(*Build a supremum for a function table f defined on a partition.

  minDistTable = flattened table of distances of nearest corners of grid squares
  maxDistTable = flattened table of distances of farthest corners of grid squares
  f = flattened table giving function defined on the same grid
  partition = array of points giving the intervals for the supremum
*)
buildSupremum[minDistTable_,maxDistTable_,f_,partition_] :=
	Module[{supremum,start,end,left,right,fval},
		supremum=Table[-Infinity,Length[partition]];
		Do[
			If[Mod[i,100000]==0,Print[i];];
			fval=f[[i]];
			{left,right}={minDistTable[[i]],maxDistTable[[i]]};
			{start,end}=getOverlapping[left,right,partition];
			(*Print[StringForm["i=``,num=``,start=``,end=``,fval=``",i,end-start,start,end,fval]];*)
			supremum[[start;;end]]=Map[Max[#1,fval]&,supremum[[start;;end]]];
			,{i,Length[minDistTable]}
		];
		Print["Finished supremum"];
		supremum
	]
(*Compute matrix norms by summing functions using distances from spikes (worst-case distances). Norms are returned in block matrix form of identity minus bump-wave interpolation matrix.

  funList = table containing monotonized bump, wave and derivative functions (for one set of grid sizes)
  Delta = minimum spike separation, to scale distances
  dists = unscaled distances of spike cells from origin
  distTable = table of distances for the function arrays
  eps = bound on contribution from all very distant spikes
*)
computeNormMat[funList_,Delta_,dists_,distTable_,eps_] :=
	Block[{ab,aw1,aw2,ad1b,ad1w1,ad1w2,ad2b,adb,ad2w1,ad2w2,inds,IB,D1B,D2B,W1,ID1W1,D2W1,W2,D1W2,ID2W2},
		ab = funList[[1]];
		aw1 = funList[[2]];
		aw2 = funList[[3]];
		ad1b = funList[[5]];
		ad1w1 = funList[[6]];
		ad1w2 = funList[[7]];
		ad2b = funList[[8]];
		ad2w1 = funList[[9]];
		ad2w2 = funList[[10]];
		inds = Map[getIndex[#,distTable]&,Delta*dists];
		IB = 6*Sum[ab[[i]],{i,inds}]+eps;
		D1B = 6*Sum[ad1b[[i]],{i,inds}]+eps; (*Sums are trimmed using tail bounds on the bumps, waves, and derivatives*)
		D2B = 6*Sum[ad2b[[i]],{i,inds}]+eps;
		W1 = 6*Sum[aw1[[i]],{i,inds}]+eps;
		ID1W1 = 6*Sum[ad1w1[[i]],{i,inds}]+eps;
		D2W1 = 6*Sum[ad2w1[[i]],{i,inds}]+eps;
		W2 = 6*Sum[aw2[[i]],{i,inds}]+eps;
		D1W2 = 6*Sum[ad1w2[[i]],{i,inds}]+eps;
		ID2W2 = 6*Sum[ad2w2[[i]],{i,inds}]+eps;
		{{IB,W1,W2},{D1B,ID1W1,D1W2},{D2B,D2W1,ID2W2}}
	]
(*Routines to compute alpha, beta bounds.

  funList = table containing monotonized bump, wave and derivative functions (for one set of grid sizes)
  Delta = minimum spike separation, to scale distances
  distTable = table of distances for the function arrays
  eps = bound on contribution to norms from all very distant spikes
*)
computeAlphaBeta[funList_,Delta_,distTable_,eps_] :=
	Block[{IB,D1B,D2B,W1,ID1W1,D2W1,W2,D1W2,ID2W2,D2W2Inv,IS1,S1Inv,S2,IS3,S3Inv,alpha,beta,alpha1},
		{{IB,W1,W2},{D1B,ID1W1,D1W2},{D2B,D2W1,ID2W2}} = computeNormMat[funList,Delta,spikeDistances[],distTable,eps];
		If[ID2W2 >= 1-EPS,Return[$Failed],Null];
		D2W2Inv = 1/(1-ID2W2);
		IS1 = ID1W1+D1W2*D2W2Inv*D2W1; (*Infinity norm bound on I-S1 where S1 is the Schur complement of D2W2*)
		If[IS1 >= 1-EPS,Return[$Failed],Null];
		S1Inv = 1/(1-IS1); (*Infinity norm bound on inverse of S1*)
		S2 = D1B+D1W2*D2W2Inv*D2B;
		IS3 = IB+W1*S1Inv*S2+W2*D2W2Inv*(D2W1*S1Inv*S2+D2B); (*Infinity norm bound on I-S3 where S3 is the Schur complement of DW*)
		If[IS3 >= 1-EPS,Return[$Failed],Null];
		S3Inv = 1/(1-IS3); (*Infinity norm bound on inverse of S3*)
		alpha = S3Inv;
		beta = S1Inv*S2*S3Inv;
		alpha1 = 1-S3Inv*IS3;
		{alpha,beta,alpha1}
	]


(*Computes the bounds on Q, DQ, DDQ over the interval [0,Delta/2].
  Assumes a single length or range for grid spacing. Thus, this
  should be computed for each choice of grid spacing. Spike at the
  origin is assumed. Returns $Failed if computeAlphaBeta fails.

  funList = table containing monotonized bump, wave and derivative functions
  supFunList = table containing non-monotonized function supremums
  Delta = minimum spike separation, to scale distances
  dists = array containing distances from segments on x-axis to each spike cell, already scaled by Delta
  eps = bound on contribution to norms from all very distant spikes
  
  RETURNS:
  Q = upper bound on dual combination over interval [0,Delta], divided into numSegments sub-intervals
  DQ = upper bound on dual combination first derivative in radially outward direction
  DDQ = upper bound on max eigenvalue of Hessian of dual combination
  QLowerBound = lower bound on dual combination
*)
computeBounds[funList_,supFunList_,Delta_,dists_,eps_] :=
	Block[{distTable,distTable2,numSegments,numCells,coeffs,alpha,beta,alpha1,endpoints,indexTable,indexRanges,farIndexTable,farIndexTableExtended,
				ab,aw1,aw2,dbdotn,ad1b,ad1w1,ad1w2,ad2b,adb,ad2w1,ad2w2,lb,alb,alw1,alw2,db,ddb,
				B,W,FarB,FarW,DB,DW,FarDB,FarDW,DDB,DDW,FarDDB,FarDDW,Q,DQ,DDQ,QLowerBound
		},
		{ab,aw1,aw2,dbdotn,ad1b,ad1w1,ad1w2,ad2b,ad2w1,ad2w2,lb,alb,alw1,alw2,distTable} = funList;
		{db,ddb,distTable2} = supFunList;
		coeffs = computeAlphaBeta[funList,Delta,distTable,eps];
		If[coeffs===$Failed,Return[$Failed],{alpha,beta,alpha1}=coeffs];
		numSegments = Dimensions[dists][[1]];
		numCells = Dimensions[dists][[2]];
		endpoints = Table[i*Delta/(numSegments),{i,0,numSegments}];
		(*The ranges for indices for near bump's first and second derivatives*)
		indexRanges = Table[{getIndex[endpoints[[i]],distTable2],getIndex[endpoints[[i+1]],distTable2]+1},{i,numSegments}];
		(*Print["Making near index table..."];*)
		indexTable = Table[getIndex[d,distTable],{d,endpoints}];
		(*Print["Making far index table..."];*)
		farIndexTable = ParallelTable[getIndex[dists[[i,j]],distTable],{i,numSegments},{j,numCells},DistributedContexts -> Automatic];
		(*Max[dists,endpoints] takes the max of either the interval to cell distance or the interval's distance to nearest spike, since all spikes must be farther than the nearest one.*)
		farIndexTableExtended = ParallelTable[getIndex[Max[dists[[i,j]],endpoints[[i]]],distTable],{i,numSegments},{j,numCells},DistributedContexts -> Automatic];
		(*Print["Summing Tables"];*)
		B = Table[alpha*Max[ab[[indexTable[[i]];;indexTable[[i+1]]]]],{i,numSegments}];
		W = Table[beta*(aw1[[indexTable[[i]]]]+aw2[[indexTable[[i]]]]),{i,numSegments}];
		FarB = Table[alpha*Sum[ab[[farIndexTableExtended[[i,j]]]],{j,numCells}],{i,numSegments}];
		FarW = Table[beta*Sum[aw1[[farIndexTableExtended[[i,j]]]]+aw2[[farIndexTableExtended[[i,j]]]],{j,numCells}],{i,numSegments}];
		DB = Table[Max[alpha1*Max[db[[indexRanges[[i,1]];;indexRanges[[i,2]]]]],alpha*Max[db[[indexRanges[[i,1]];;indexRanges[[i,2]]]]]],{i,numSegments}];
		DW = Table[beta*(Sqrt[ad1w1[[indexTable[[i]]]]^2+ad2w1[[indexTable[[i]]]]^2]+Sqrt[ad1w2[[indexTable[[i]]]]^2+ad2w2[[indexTable[[i]]]]^2]),{i,numSegments}];
		FarDB = Table[alpha*Sum[Sqrt[ad1b[[farIndexTable[[i,j]]]]^2+ad2b[[farIndexTable[[i,j]]]]^2],{j,numCells}],{i,numSegments}];
		FarDW = Table[beta*Sum[Sqrt[ad1w1[[farIndexTable[[i,j]]]]^2+ad2w1[[farIndexTable[[i,j]]]]^2]+Sqrt[ad1w2[[farIndexTable[[i,j]]]]^2+ad2w2[[farIndexTable[[i,j]]]]^2],{j,numCells}],{i,numSegments}];
		DDB = Table[Max[alpha1*Max[ddb[[indexRanges[[i,1]];;indexRanges[[i,2]]]]],alpha*Max[ddb[[indexRanges[[i,1]];;indexRanges[[i,2]]]]]],{i,numSegments}];
		DDW = Table[beta*(alw1[[indexTable[[i]]]]+alw2[[indexTable[[i]]]]),{i,numSegments}];
		FarDDB = Table[alpha*Sum[alb[[farIndexTable[[i,j]]]],{j,numCells}],{i,numSegments}];
		FarDDW = Table[beta*Sum[alw1[[farIndexTable[[i,j]]]]+alw2[[farIndexTable[[i,j]]]],{j,numCells}],{i,numSegments}];
		Q = B+W+FarB+FarW;
		DQ = DB+DW+FarDB+FarDW;
		DDQ = DDB+DDW+FarDDB+FarDDW;
		QLowerBound = -W-FarB-FarW;
		{Q,DQ,DDQ,QLowerBound}
	]


(*Integration routines to prove boundedness of dual combination Q.*)
(*Calculate an integral across subintervals of [0,Delta/2]
  to determine where the following expression is negative:
  
  \int_{0}^{t}v^TH(x_0+sv)v (t-s) ds
  
  DDQ = upper bound on max eigenvalue of Hessian of dual combination Q
*)
DDQIntegral[DDQ_] :=
	Module[{endpoints,integralTable,numSegments},
		numSegments = Length[DDQ];
		endpoints = Table[i/(numSegments),{i,0,numSegments}];
		(*Note endpoints table indices must be indexed starting from 1 not 0.*)
		(*integralTable = Table[Sum[Integrate[DDQ[[i]]*(endpoints[[j+1]]-x),{x,endpoints[[i]],endpoints[[i+1]]}],{i,1,j}],{j,numSegments}];*)
		integralTable = Table[Sum[DDQ[[i]]*(endpoints[[j+1]]/(numSegments)-(endpoints[[i+1]]^2-endpoints[[i]]^2)/2),{i,1,j}],{j,numSegments}];(*This is the calculated integral from the expression*)
		integralTable
	]
(*Calculate an integral across subintervals of [0,Delta/2]
  to determine where the following expression is negative:

  \int_{u_1}^{t}\nabla Q(x_0+sv)\cdot v ds
  
  DQ = upper bound on dual combination first derivative in radially outward direction
*)
DQIntegral[DQ_] :=
	Module[{endpoints,integralTable,numSegments},
		numSegments = Length[DQ];
		integralTable = Table[Sum[DQ[[i]],{i,k}],{k,1,numSegments}];
		integralTable
	]


(*Uses integral method to determine if Q is bounded.
  The output of computeBounds DQ, and DDQ are 
  integrated over segments spanning distance 0 to Delta/2
  from spike at the origin. Boundedness of Q is 
  determined over each segment. Tables for recovery and
  lengths for boundedness for DDQ and DQ are saved.
  
  monotoneFunList = table containing monotonized bump, wave and derivative functions
  supFunList = table containing non-monotonized function supremums
  DeltaList = array of values to use for minimum spike separations
  absoluteDists = array containing distances from segments on x-axis to each spike cell, NOT yet scaled by Delta
  eps = marginal term, used to bound contributions from distant spikes or for error tolerance purposes
*)
makeQTables[monotoneFunList_,supFunList_,DeltaList_,absoluteDists_,eps_]:=
	Module[{ends},
		ends = Table[
			Module[{QList,Q,DQ,DDQ,QLowerBound,DQI,DDQI,u1,u2,u3,numSegments,dx,j,DQTest,DDQTest,isQBound,X},
				Print["i: ",i,"; Delta: ",Delta];
				QList = computeBounds[monotoneFunList[[i]],supFunList[[i]],Delta,Delta*absoluteDists,eps];
				If[QList===$Failed,Return[{0,0,0},Module]];
				{Q,DQ,DDQ,QLowerBound} = QList;
				numSegments = Length[Q];
				dx = 1/(numSegments);
				DDQI = DDQIntegral[DDQ];
				DQI = DQIntegral[DQ];
				u1 = 0;
				u2 = 0;
				u3 = 0;
				isQBound = True;
				DQTest = True;
				DDQTest = True;
				For[j=1,j<=Length[Q],j++,
					If[DDQTest&&(DDQI[[j]]<0),u1=j;, (*If DDQ integral was negative and remains negative, update u1*)
						DDQTest=False; (*Else update that DDQ integral is no longer negative and check...*)
						If[DQTest&&u1>0&&(Min[Table[DQI[[j]]-DQI[[k]],{k,1,u1}]]<0),u2=j;, (*... if DQ integral is negative*)
							DQTest=False;
							If[Q[[j]]>=1,isQBound=False;];
						];
					];
				];
				If[isQBound,u3=1;,u3=0;];
				{u1*dx*1.,u2*dx*1.,u3}
			]
			,{Delta,DeltaList}
			,{i,Dimensions[monotoneFunList][[1]]}
		];
		Export["/home/mcdonald/Dropbox/Granda/20180507paper/asy/dat/Deltas.dat",DeltaList,"TSV"]; (*For proper alignment in a plot, an extra Delta value might need to be added to this array.*)
		Export["/home/mcdonald/Dropbox/Granda/20180507paper/asy/dat/DDQ.dat",ends[[;;,;;,1]],"TSV"];
		Export["/home/mcdonald/Dropbox/Granda/20180507paper/asy/dat/DQ.dat",ends[[;;,;;,2]],"TSV"];
		Export["/home/mcdonald/Dropbox/Granda/20180507paper/asy/dat/Q.dat",ends[[;;,;;,3]],"TSV"];
	]
(*Determine if Q is small enough for points farther than
  Delta. In that case, we can use the matrix norms for I-B,
  W1 and W2 and the coefficient bounds to get a simple
  upper bound.
*)
makeQTablesForDistantPoints[monotoneFunList_,DeltaList_,eps_]:=
	Module[{results},
		results = Table[
			Module[{funList,distTable,IB,D1B,D2B,W1,ID1W1,D2W1,W2,D1W2,ID2W2,coeffs,alpha,beta,alpha1},
				funList = monotoneFunList[[i]];
				distTable = funList[[15]];
				{{IB,W1,W2},{D1B,ID1W1,D1W2},{D2B,D2W1,ID2W2}} = computeNormMat[funList,Delta,spikeDistances[],distTable,eps];
				coeffs = computeAlphaBeta[funList,Delta,distTable,eps];
				If[coeffs===$Failed,Return[100,Module],{alpha,beta,alpha1}=coeffs];
				alpha*IB+beta*W1+beta*W2
			]
			,{Delta,DeltaList}
			,{i,Dimensions[monotoneFunList][[1]]}
		];
		Export["/home/mcdonald/Dropbox/Granda/20180507paper/asy/dat/QDistantPoints.dat",results,"TSV"];
	]
makeExtendedTables[monotoneFunList_,supFunList_,DeltaList_,absoluteDists_,eps_]:=
	Module[{ends},
		ends = Table[
			Module[{QList,Q,DQ,DDQ,QLowerBound,numSegments,j,isQBound},
				Print["i: ",i,"; Delta: ",Delta];
				QList = computeBounds[monotoneFunList[[i]],supFunList[[i]],Delta,Delta*absoluteDists,eps];
				If[QList===$Failed,Return[0,Module]];
				{Q,DQ,DDQ,QLowerBound} = QList;
				isQBound = True;
				For[j=50,j<=Length[Q],j++,If[Q[[j]]>=1,isQBound=False;]];
				If[isQBound,1,0]
			]
			,{Delta,DeltaList}
			,{i,Dimensions[monotoneFunList][[1]]}
		];
		Export["/home/mcdonald/Dropbox/Granda/20180507paper/asy/dat/QExt.dat",ends,"TSV"];
	]


(*Tests whether the dual combination Q is properly bounded.
  testQRegionBound determines if the necessary conditions
  properly overlap.
  
  Q = upper bound on dual combination over interval [0,Delta/2], divided into numSegments sub-intervals
  DQ = upper bound on dual combination first derivative in radially outward direction
  DDQ = upper bound on max eigenvalue of Hessian of dual combination
*)
testQRegionBound[Q_,DQ_,DDQ_] :=
	Module[{isDQNeg,isDDQNeg,isQBound,i},
		isQBound = False;
		isDQNeg = False;
		isDDQNeg = True;
		For[i = 1, i <= Length[Q], i++,
			If[isDDQNeg && DDQ[[i]] < 0, isDDQNeg = True, isDDQNeg = False];
			If[(isDDQNeg || isDQNeg) && DQ[[i]] < 0, isDQNeg = True, isDQNeg = False];
			If[(isDDQNeg || isDQNeg || isQBound) && Q[[i]] < 1, isQBound=True, isQBound=False];
			If[!isDDQNeg && !isDQNeg && !isQBound, isQBound=False; Break[], Null];
		];
		isQBound
	]


(*Some I/O methods and useful wrapper functions.*)
importMatrix[matfile_] := Import[matfile][[1]]
getPartitionDistancesToSpikes[distFileString_,Delta_] := Delta*importMatrix[distFileString]
importAndFlatten[matfile_] := Flatten[importMatrix[matfile]]
getNearFunctionTables[partitionfile_,supremumDBfile_,supremumDDBfile_] :=
	Module[{partition,supremumDB,supremumDDB},
		partition = importAndFlatten[partitionfile];
		supremumDB = importAndFlatten[supremumDBfile];
		supremumDDB = importAndFlatten[supremumDDBfile];
		{supremumDB,supremumDDB,partition}
	]
(*Check for successful recovery. OLD: makeQTables makes DAT files for recovery
  results. testQRegionBound doesn't use integral method and is less accurate.

  Delta = minimum spike separation
  gridSep = function grid separation distance
  numSep = number of function grid points in each direction
  functionTable = unflattened table of bump, wave and derivative functions
  distFileString = filename of distances from x-axis intervals to spike cells
  eps = error tolerance bound
*)
testSuccess[Delta_,gridSep_,numSep_,functionTable_,distFileString_,eps_] :=
	Module[{dists,numPartitions,distTable,distfuntables,partition,supremumDB,supremumLB,neardistancefunctiontables,intervals,Q,DQ,DDQ,QLowerBound},
		dists = getPartitionDistancesToSpikes[distFileString,Delta];
		distfuntables = makeDistFunTables[gridSep,numSep,functionTable,eps];
		neardistancefunctiontables = getNearFunctionTables["partition.mat","DBdotNsupremum.mat","LBsupremum.mat"];
		{Q,DQ,DDQ,QLowerBound}=computeBounds[distfuntables,neardistancefunctiontables,1,Delta,dists,eps];
		testQRegionBound[Q,DQ,DDQ]
	]
(*OLD: testQRegionBound doesn't use integral method and
  is less accurate. Organizes recovery results as table
  for different minimum separation and grid spacing choices.

  monotoneFunList = table containing monotonized bump, wave and derivative functions
  supFunList = table containing non-monotonized function supremums
  DeltaList = array of values to use for minimum spike separations
  absoluteDists = array containing distances from segments on x-axis to each spike cell, NOT yet scaled by Delta
  eps = marginal term, used to bound contributions from distant spikes or for error tolerance purposes
*)
makeRecoveryTables[monotoneFunList_,supFunList_,DeltaList_,absoluteDists_,eps_]:=
	Table[
		Module[{QList,Q,DQ,DDQ,DQI,DDQI,QLowerBound},
			Print["i: ",i,"; Delta: ",Delta];
			QList = computeBounds[monotoneFunList[[i]],supFunList[[i]],Delta,Delta*absoluteDists,eps];
			If[QList===$Failed,False,{Q,DQ,DDQ,QLowerBound}=QList;DDQI=DDQIntegral[DDQ];DQI=DQIntegral[DQ];testQRegionBound[Q,DQI,DDQI]]
		]
		,{i,Dimensions[monotoneFunList][[1]]}
		,{Delta,DeltaList}
	]


(*Import unflattened (i.e. 4 dimensional) function tables.*)
compileTables[]:=
	Table[
		Module[{filename},
			filename = ToString[StringForm["~/Data/flg_20190531_``.mat",i]];
			Import[filename]
		]
		,{i,16}
	]
compileTablesEV[]:=
	Table[
		Module[{filename},
			filename = ToString[StringForm["~/Data/ev_20190829_``.mat",i]];
			Import[filename]
		]
		,{i,16}
	]
(*Create flattened, monotonized function arrays from an unflattened table.

  gridSep = function grid separation distance
  numSep = number of function grid points in each direction
  funList = unflattened table of bump, wave and derivative functions
  eps = error tolerance bound
*)
compileFunListMonotone[gridSep_,numSep_,funList_,eps_]:=Table[makeDistFunTables[gridSep,numSep,fList,eps],{fList,funList}]
(*Create flattened, non-monotonized supremum tables for gradient and eigenvalue functions from function table.

  gridSep = function grid separation distance
  numSep = number of function grid points in each direction
  funList = unflattened table of bump, wave and derivative functions
*)
compileFunListSupremum[gridSep_,numSep_,funList_]:=
	Module[{minDistTable,maxDistTable,partition},
		maxDistTable = makeMaxDistanceTable[gridSep,numSep];
		minDistTable = makeMinDistanceTable[gridSep,numSep];
		partition = Sort[DeleteDuplicates[Flatten[{minDistTable,maxDistTable}]]];
		Print[DateString[]];
		Table[
			Module[{sup,DBDN,LB},
				DBDN = buildSupremum[Flatten[minDistTable],Flatten[maxDistTable],Flatten[fList[[4]]],partition];
				LB = buildSupremum[Flatten[minDistTable],Flatten[maxDistTable],Flatten[fList[[11]]],partition];
				Print[DateString[]];
				{DBDN,LB,partition}
			]
			,{fList,funList}
		]
	]


(*Plot non-monotonized supremum envelopes. Note buildSupremum can take several minutes.

  gridSep = function grid separation distance
  numSep = number of function grid points in each direction
  functionTable = unflattened table of bump, wave and derivative functions
*)
plotSupremum[gridSep_,numSep_,functionTable_]:=
	Module[{sup,minDistTable,maxDistTable,partition},
		maxDistTable = makeMaxDistanceTable[gridSep,numSep];
		minDistTable = makeMinDistanceTable[gridSep,numSep];
		partition = Sort[DeleteDuplicates[Flatten[{minDistTable,maxDistTable}]]];
		sup = buildSupremum[Flatten[minDistTable],Flatten[maxDistTable],Flatten[functionTable],partition];
		ListLinePlot[Table[{partition[[i]],sup[[i]]},{i,Length[partition]}],PlotRange->All,ImageSize->Full]
	]
(*Plot non-monotonized supremum envelopes that are already built and imported.

  gridSep = function grid separation distance
  numSep = number of function grid points in each direction
  functionTable = unflattened table of bump, wave and derivative functions
  gridSizeIndex = index determining range of sample separations
  supIndex = index choosing supremum of either gradient (DBDN) or eigenvalue (L)
*)
plotSupremumFromFile[gridSep_,numSep_,functionTable_,gridSizeIndex_,supIndex_]:=
	Module[{sup,partition},
		partition = functionTable[[gridSizeIndex,3]];
		sup = functionTable[[gridSizeIndex,supIndex]];
		ListLinePlot[Table[{partition[[i]],sup[[i]]},{i,Length[partition]}],PlotRange->All,ImageSize->Large]
	]
(*Monotonize then plot function envelopes.

  gridSep = function grid separation distance
  numSep = number of function grid points in each direction
  functionTable = unflattened 2D table of bump, wave or derivative function
  eps = error tolerance bound
*)
plotMonotone[gridSep_,numSep_,functionTable_,eps_]:=
	Module[{table,partition},
		{table,partition} = makeDistFunTables[gridSep,numSep,{functionTable},eps];
		ListLinePlot[Table[{partition[[i]],table[[i]]},{i,Length[partition]}],PlotRange->All,ImageSize->Full]
	]
(*Plot envelopes that are already monotized and imported.

  gridSep = function grid separation distance
  numSep = number of function grid points in each direction
  functionTable = unflattened table of bump, wave and derivative functions
  gridSizeIndex = index determining range of sample separations
  supIndex = index choosing supremum of either gradient (DBDN) or eigenvalue (L)
*)
plotMonotoneFromFile[gridSep_,numSep_,functionTable_,gridSizeIndex_,supIndex_]:=
	Module[{sup,partition},
		partition = functionTable[[gridSizeIndex,15]];
		sup = functionTable[[gridSizeIndex,supIndex]];
		ListLinePlot[Table[{partition[[i]],sup[[i]]},{i,Length[partition]}],PlotRange->All,ImageSize->Large]
	]
(*Plot dual combination and derivatives Q, DQ and DDQ over [0,Delta/2]

  monotoneFunList = table containing monotonized bump, wave and derivative functions (of ALL sample grid sizes)
  supFunList = table containing non-monotonized function supremums
  gridSizeIndex = index determining range of sample separations
  Delta = minimum spike separation
  absoluteDists = array containing distances from segments on x-axis to each spike cell, NOT yet scaled by Delta
  eps = error tolerance bound
*)
plotBoundList[monotoneFunList_,supFunList_,gridSizeIndex_,Delta_,absoluteDists_,eps_] :=
	Module[{Q,DQ,DDQ,QLowerBound,X,line},
		{Q,DQ,DDQ,QLowerBound}=computeBounds[monotoneFunList[[gridSizeIndex]],supFunList[[gridSizeIndex]],Delta,Delta*absoluteDists,eps];
		X = Range[.01,.5,.01];
		ListLinePlot[{Table[{X[[i]],Q[[i]]},{i,Length[X]}],Table[{X[[i]],DQ[[i]]},{i,Length[X]}],Table[{X[[i]],DDQ[[i]]},{i,Length[X]}]},PlotRange->All,PlotLegends->{"Q","DQ","DDQ"},GridLines->{{},{1}}]
		(*ListLinePlot[{Table[{X[[i]],Q[[i]]},{i,Length[X]}],Table[{X[[i]],DQ[[i]]},{i,Length[X]}]},PlotRange->All,PlotLegends->{"Q","DQ"},GridLines->{{},{1}}]*)
	]


(*Routines to organize tables so they can be saved
  in formats that are easily imported by asymptote
  files for plotting.

  monotoneFunList = table containing monotonized bump, wave and derivative functions
  DeltaList = array of values to use for minimum spike separations
  eps = marginal term, used to bound contributions from distant spikes or for error tolerance purposes
*)
makeMatrixNorms[monotoneFunList_,DeltaList_,eps_]:=
	Module[{tab},
		tab = Table[
			Module[{IB,D1B,D2B,W1,ID1W1,D2W1,W2,D1W2,ID2W2},
				{{IB,W1,W2},{D1B,ID1W1,D1W2},{D2B,D2W1,ID2W2}} = computeNormMat[fList,Delta,spikeDistances[],fList[[15]],eps];
				{IB,W1,W2,D1B,ID1W1,D1W2,D2B,D2W1,ID2W2}
			]
			,{fList,monotoneFunList}
			,{Delta,DeltaList}
		];
		Transpose[tab,{2,3,1}]
	]
makeNormTables[monotoneFunList_,DeltaList_,eps_]:=
	Module[{tab,ID2W2vals,IS1vals,IS3vals},
		tab = Table[
			Module[{IB,D1B,D2B,W1,ID1W1,D2W1,W2,D1W2,ID2W2,D2W2Inv,IS1,S1Inv,S2,IS3},
				{{IB,W1,W2},{D1B,ID1W1,D1W2},{D2B,D2W1,ID2W2}} = computeNormMat[fList,Delta,spikeDistances[],fList[[15]],eps];
				D2W2Inv = 1/(1-ID2W2);
				IS1 = ID1W1+D1W2*D2W2Inv*D2W1; (*Infinity norm bound on I-S1 where S1 is the Schur complement of D2W2*)
				S1Inv = 1/(1-IS1); (*Infinity norm bound on inverse of S1*)
				S2 = D1B+D1W2*D2W2Inv*D2B;
				IS3 = IB+W1*S1Inv*S2+W2*D2W2Inv*(D2W1*S1Inv*S2+D2B); (*Infinity norm bound on I-S3 where S3 is the Schur complement of DW*)
				If[ID2W2 >= 1-EPS,{IS1,IS3}={100000,100000},If[IS1 >= 1-EPS,IS3=100000]];
				{ID2W2,IS1,IS3}
			]
			,{fList,monotoneFunList}
			,{Delta,DeltaList}
		];
		ID2W2vals = Table[tab[[i,j,1]],{i,Dimensions[monotoneFunList][[1]]},{j,Length[DeltaList]}];
		IS1vals = Table[tab[[i,j,2]],{i,Dimensions[monotoneFunList][[1]]},{j,Length[DeltaList]}];
		IS3vals = Table[tab[[i,j,3]],{i,Dimensions[monotoneFunList][[1]]},{j,Length[DeltaList]}];
		{ID2W2vals,IS1vals,IS3vals}
	]
makeAlphaBetaTables[monotoneFunList_,DeltaList_,eps_]:=
	Module[{tab,alphaVals,betaVals,alpha1Vals},
		tab = Table[
			Module[{v},
				v = computeAlphaBeta[fList,Delta,fList[[15]],eps];
				If[v===$Failed,{-10000,-10000,10000},v]
			]
			,{fList,monotoneFunList}
			,{Delta,DeltaList}
		];
		alphaVals = Table[tab[[i,j,1]],{i,Dimensions[monotoneFunList][[1]]},{j,Length[DeltaList]}];
		betaVals = Table[tab[[i,j,2]],{i,Dimensions[monotoneFunList][[1]]},{j,Length[DeltaList]}];
		alpha1Vals = Table[1-tab[[i,j,3]],{i,Dimensions[monotoneFunList][[1]]},{j,Length[DeltaList]}];
		{alphaVals,betaVals,alpha1Vals}
	]


(*Check the lower limit for spike separation such
  that Q > -1. Needed for boundedness proof.

  monotoneFunList = table containing monotonized bump, wave and derivative functions
  supFunList = table containing non-monotonized function supremums
  DeltaList = array of values to use for minimum spike separations
  absoluteDists = array containing distances from segments on x-axis to each spike cell, NOT yet scaled by Delta
  eps = marginal term, used to bound contributions from distant spikes or for error tolerance purposes
*)
lowerBoundSuccess[monotoneFunList_,supFunList_,DeltaList_,absoluteDists_,eps_]:=
	Module[{tab},
		tab = Table[
			Module[{QList,QLowerBound,success},
				Print["i: ",i,"; Delta: ",Delta];
				success = -1;
				QList = computeBounds[monotoneFunList[[i]],supFunList[[i]],Delta,Delta*absoluteDists,eps];
				If[QList===$Failed,success=0];
				QLowerBound = QList[[4]];
				If[Max[QLowerBound]<1,success=1,success=0];
				success
			]
			,{Delta,DeltaList}
			,{i,Dimensions[monotoneFunList][[1]]}
		];
		Export["/home/mcdonald/Dropbox/Granda/20180507paper/asy/dat/QLowerBound.dat",tab,"TSV"];
	]


(*Plots tables of booleans with given tick lists.

  table = table values to plot
  gridSizes = grid spacing
  DeltaList = array of values to use for minimum spike separations
  K = kernel function
*)
plotTable[table_,gridSizes_,DeltaList_,K_] :=
	Module[{m,xTicks,yTicks,t},
		yTicks = Table[{k,DeltaList[[k]]},{k,1,Length[DeltaList]}];
		xTicks = Table[{k,Rotate[Text[gridSizes[[k]]],Pi/2]},{k,1,Length[gridSizes]-1}];
		m = MatrixPlot[Table[If[table[[j,i]],1,0],{i,1,Length[DeltaList]},{j,1,Length[gridSizes]-1}],FrameTicks->{yTicks,xTicks},Mesh->All];
		Show[m,AxesLabel->{None,None},FrameLabel->{{"\[CapitalDelta] (in units of \[Sigma])",None},{"\[Gamma] (in units of \[Sigma])",None}},PlotLabel->"K[t]="TraditionalForm[K[Global`t]],LabelStyle->{GrayLevel[0]}]
	]
(*Plots tables of values with given tick lists.*)
plotValueTable[table_,gridSizes_,DeltaList_,K_] :=
	Module[{m,xTicks,yTicks,t},
		yTicks = Table[{k,DeltaList[[k]]},{k,1,Length[DeltaList]}];
		xTicks = Table[{k,Rotate[Text[gridSizes[[k]]],Pi/2]},{k,1,Length[gridSizes]}];
		m = MatrixPlot[Table[table[[j,i]],{i,1,Length[DeltaList]},{j,1,Length[gridSizes]-1}],FrameTicks->{yTicks,xTicks},Mesh->All,PlotLegends -> Placed[BarLegend[Automatic], {After, Center}]];
		Show[m,AxesLabel->{None,None},FrameLabel->{{"\[CapitalDelta] (in units of \[Sigma])",None},{"\[Gamma] (in units of \[Sigma])",None}},PlotLabel->"K[t]="TraditionalForm[K[Global`t]],LabelStyle->{GrayLevel[0]}]
	]


End[]
EndPackage[]
