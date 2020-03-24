# Deconvolution2D

The necessary packages to reconstruct the proof are listed and described below:

|File | Description|
|-----|-----|
|TwoDimBumpWave.wl | Library defining functions needed for computational proof.|
|ServerScript.wl   | Script computes and saves function tables using Interval Arithmetic, capable of running on a compute cluster.|                       
|run-math.sbatch | Script to run ServerScript.wl using Slurm scheduler.|                        
|segmentdistances.m | MATLAB script to compute distances between hexagon cells and segment of length Delta, where Delta is the minimum separation between  spikes.|                        
|MakeRecoveryTables.wl | Script computing monotonized and supremum function envelopes (tables). Pre-computed tables can be loaded. Also computes tables for matrix norms, bump wave coefficients, and recovery results.|

To create the function envelopes, the Mathematica script `ServerScript.wl` must be run.
This script computes supremums for bump and wave functions for a selection of grid spacing intervals, defined therein, using Interval Arithmetic.
We include a script `run-math.sbatch` which calls `ServerScript.wl` and provides options for a cluster managed by Slurm job scheduler.

The MATLAB script `segmentdistances.m` computes distances between a set of hexagonal cells surrounding the origin and N intervals of equal lengths that partition the segment on the horizontal axis from the origin to (1,0).
Its output is a file with name formatted as `partition_dists_YYYYMMDD_npartitions_N.mat` with date and N substituted.
By default we set `N=100`.

After running `segmentdistances.m`, recovery results can be obtained by running `MakeRecoveryTables.wl`.
First this script computes envelopes, both monotonized and non-monotonic as needed, for the bump and wave functions using the output of `ServerScript.wl`.
Since this computation is time-consuming these envelopes are saved and can be imported subsequently in future uses by uncommenting the relevant lines.
With these envelopes, recovery results are produced for a list of minimum separation sizes.
