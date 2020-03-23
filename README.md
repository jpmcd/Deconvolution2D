# Deconvolution2D

The necessary packages to reconstruct the proof are located in
`deconv_2D_local/deconvolution_proof/2D`:

TwoDimBumpWave.wl     - Library defining functions needed for
                        computational proof.
                        
ServerScript.wl       - Script computes and saves function tables
                        using Interval Arithmetic, capable of running
                        on a compute cluster.
                        
run-math.sbatch       - Slurm script to run ServerScript.wl using NYU
                        Compute Cluster.
                        
segmentdistances.m    - MATLAB script to compute distances between
                        hexagon cells and segment of length Delta,
                        where Delta is the minimum separation between
                        spikes.
                        
MakeRecoveryTables.wl - Script computing monotonized and supremum
                        function envelopes (tables). Pre-computed
                        tables can be loaded. Also computes tables for
                        matrix norms, bump wave coefficients, and
                        recovery stats, and saves in DAT files.
