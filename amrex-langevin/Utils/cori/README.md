# Running on Cori

This directory has sample SLURM scripts for running on Cori Haswell and Cori KNL nodes.

Compiling for Cori KNL requires the `craype-mic-knl` module loaded,
so starting from the default modules:

`module swap craype-haswell craype-mic-knl`

To use the default intel compilers, use the `COMP=intel` flag to `make`.

Or to load the GNU compilers:

`module swap PrgEnv-intel/6.0.5 PrgEnv-gnu/6.0.5`
