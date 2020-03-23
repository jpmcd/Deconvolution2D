# Deconvolution2D

The necessary packages to reconstruct the proof are located in
`deconv_2D_local/deconvolution_proof/2D`:

TwoDimBumpWave.wl     - Library containing function definitions
                        needed for computational proof
                        
ServerScript.wl       - Script computes and saves function envelopes
                        using Interval Arithmetic, capable of running
                        on a compute cluster
                        
run-math.sbatch       - Slurm script to run ServerScript.wl using NYU
                        Compute Cluster
                        
segmentdistances.m    - MATLAB script to compute distances between
                        Delta/2 length segment and hexagon cells
                        
MakeRecoveryTables.nb - Notebook to compute monotonized and supremum
                        function tables. Pre-computed tables can be
                        loaded. Also computes tables for matrix norms,
                        bump wave coefficients, and recovery stats,
                        and saves in DAT files.
                        
NearestSpikes.nb      - Notebook to compute nearest spike positions
                        from distant points. Computes dual combination
                        over arrangements and bounds for minimum
                        separation and grid size choices.
