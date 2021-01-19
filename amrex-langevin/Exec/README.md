# Compiling

In this Exec directory, the AMReX version of NRRB can be compiled like:

```
$ make -j [N] AMREX_HOME=[path to the amrex repository directory]
```

Where [N] sets the number of parallel make threads and the path to AMReX
is the path to the top level of the amrex git repository, e.g. /home/user/amrex.

It is sometimes useful to `make realclean AMREX_HOME=[path to the amrex repository directory]`

This will delete the executable and all the object files built along the way
so it is sometimes helpful if we want to just compile all over again "from scratch"
while debugging.

Even though this is a pretty short `make` rule to delete the temporary build files,
`AMREX_HOME` still needs to be defined since all the make rules are defined
within AMReX.

If you are having issues that can't immediately be identified, pull the latest version of amrex and try again!

Passing `USE_MPI=TRUE` to `make` will compile the code with MPI enabled.

```
$ make -j [N] AMREX_HOME=[path to the amrex repository directory] USE_MPI=TRUE
```

# Running

The AMReX example can be run by executing the "main[...].ex" executable file
and passing an inputs file as an argument, e.g.:

```
$ main3d.gnu.ex inputs
```

Any of the values in `inputs` can be overriden on the command line like this:

```
$ main3d.gnu.ex inputs nrrb.l=0.5
```

# Using HDF5

The simple way to tell the AMReX build system to enable HDF5 is to pass
`USE_HDF5=TRUE` on the `make` line.

Then when running the code, set `use_hdf5=true` either in the inputs file or as
a command line argument for the executable.

Before compiling, there are some minor setup steps at NERSC vs on a local workstation:

## Parallel HDF5 on Cori at NERSC

On Cori, before compiling, we need to do the following:

```
$ module swap PrgEnv-intel PrgEnv-gnu
$ module load cray-hdf5-parallel
```

Then we can compile like, e.g.:

```
$ make -j 8 USE_HDF5=TRUE USE_MPI=TRUE USE_OMP=TRUE
```

And then if we want to quickly test this in an interactive session:

```
$ salloc -N 1 -C haswell -q interactive -t 00:30:00
$ srun -N 1 -n 4 ./main3d.gnu.haswell.MPI.OMP.ex inputs use_hdf5=true
```

## HDF5 in Local environment

On a local workstation, it can be necessary to tell the build system where to
find the HDF5 header files and library files.

The AMReX build system looks for a path stored in the `HDF5_HOME` environment
variable and in that path expects to find subdirectories `include` and `lib`
containing the headers and library files.

So on my workstation, I had to compile like this:

```
$ make -j USE_HDF5=TRUE HDF5_HOME=/home/dewillcox/local/hdf5_1.10.4_mpich_noc
```

# Verification

## Plotfiles

`Langevin-Compare.py` can be used to write out cell values from an AMReX plotfile
and also calculate the infinity or L2 norms between cell values in an AMReX
plotfile and a lattice save file from the original NRRB code.

It relies on yt (https://yt-project.org/) and was tested with Python 3.

`Langevin-Compare.py -h` will print the flags and what they do, for example:

This will compare the values in AMReX plotfile `plt00000` right after initialization
to the equivalent lattice save file:

```
$ python3 Langevin-Compare.py plt00000 -c ../../field_configs/v4_mu_-0.100_w_0.100_Nx_21_nL_0_field_config.txt

yt : [INFO     ] 2019-10-09 23:01:52,931 Parameters: current_time              = 0.0
yt : [INFO     ] 2019-10-09 23:01:52,931 Parameters: domain_dimensions         = [21 21 80]
yt : [INFO     ] 2019-10-09 23:01:52,931 Parameters: domain_left_edge          = [0. 0. 0.]
yt : [INFO     ] 2019-10-09 23:01:52,931 Parameters: domain_right_edge         = [21. 21. 80.]
Comparing plt00000 to ../../field_configs/v4_mu_-0.100_w_0.100_Nx_21_nL_0_field_config.txt ...
  field     Abs Inf Norm   Rel Inf Norm   Abs L2 Norm    Rel L2 Norm
 phi_1_Re   4.999956e-07   4.881816e-06   5.164217e-05   1.717776e-04
 phi_1_Im   4.999989e-07   4.975708e-06   5.142428e-05   1.719121e-04
 phi_2_Re   4.999866e-07   4.892570e-06   5.133029e-05   1.705982e-04
 phi_2_Im   4.999988e-07   4.943755e-06   5.136080e-05   1.716011e-04
```

That shows the absolute and relative infinity and L2 norms between the cell values.

Those are nonzero because the original NRRB code uses 6 significant figures for the lattice save.

Using the `-csf` flag, we can say how many figures we want to compare,
and the values in the AMReX plotfile will be rounded to that many S.F. for calculating the norms:

```
$ python3 Langevin-Compare.py plt00000 -c ../../field_configs/v4_mu_-0.100_w_0.100_Nx_21_nL_0_field_config.txt -csf 6

yt : [INFO     ] 2019-10-09 23:22:48,360 Parameters: current_time              = 0.0
yt : [INFO     ] 2019-10-09 23:22:48,361 Parameters: domain_dimensions         = [21 21 80]
yt : [INFO     ] 2019-10-09 23:22:48,361 Parameters: domain_left_edge          = [0. 0. 0.]
yt : [INFO     ] 2019-10-09 23:22:48,361 Parameters: domain_right_edge         = [21. 21. 80.]
Comparing plt00000 to ../../field_configs/v4_mu_-0.100_w_0.100_Nx_21_nL_0_field_config.txt ...
 (Field values in AMReX plotfile rounded to 6 significant figures for the comparison)
  field     Abs Inf Norm   Rel Inf Norm   Abs L2 Norm    Rel L2 Norm
 phi_1_Re   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
 phi_1_Im   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
 phi_2_Re   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
 phi_2_Im   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
```

## Observables

`observables-compare.py` takes an AMReX observable file (`-a` flag) and an original observable file to compare (`-c` flag)
and prints out the absolute and relative error norms for each observable, e.g.:

```
$ python3 observables-compare.py -a test.log -c ../../logfile_D_2_Nx_21_Nt_80_dt_0.05_nL_10_eps_0.01_m_1.0_wtr_0.10_w_0.1_l_0.100_mu_-0.100.log
Comparing test.log to ../../logfile_D_2_Nx_21_Nt_80_dt_0.05_nL_10_eps_0.01_m_1.0_wtr_0.10_w_0.1_l_0.100_mu_-0.100.log ...
  Observable       Abs Inf Norm   Rel Inf Norm   Abs L2 Norm    Rel L2 Norm
 Re[phi^{*}phi]    0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
 Im[phi^{*}phi]    0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
 Re[<n>]           0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
 Im[<n>]           0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
 Re[<Lz>]          3.189200e-02   3.008104e+00   5.126316e-02   4.326049e+00
 Im[<Lz>]          2.984680e-02   4.727911e-01   7.208607e-02   1.210221e+00
 Re[<S>]           7.300000e-03   1.217695e-03   1.203818e-02   1.253478e-03
```

## MPI

To compare the AMReX version with MPI enabled with the serial version of the original code,
the random numbers generated in the initialization and advance are replaced with 1.0.

The reason for this is that if we were to set a constant seed for the random number generators
each MPI rank would use that seed for its part of the domain. We would have to map the boxes
owned by each MPI rank to the parts of the domain in the original code and reset the
random number generator there separately for each sub-domain.

Since the arithmetic operations to do the advance don't change depending on
exactly which random numbers are generated, replacing the random numbers with a constant
makes this comparison very easy.

Passing the `USE_TEST_CONSTANT_RNG=TRUE` flag to `make` for both the AMReX and original versions
of the code will replace all random numbers with 1.0.

For the AMReX version with MPI, I tested with the inputs in `inputs_mpi_test` with 4 MPI ranks where
the domain is divided into boxes no larger than `8x8x8` by setting `max_grid_size = 8`.

Here are the field norms for the 10th plotfiles for each version of the code,
showing we get exactly the same field values in each case:

```
$ python3 Langevin-Compare.py mpi-4-constrng-grid8/plt00010 -c ../../field_configs/v4_mu_-0.100_w_0.100_Nx_21_nL_10_field_config.txt -csf 6
yt : [INFO     ] 2019-10-11 17:09:04,702 Parameters: current_time              = 0.09999999999999999
yt : [INFO     ] 2019-10-11 17:09:04,702 Parameters: domain_dimensions         = [21 21 80]
yt : [INFO     ] 2019-10-11 17:09:04,702 Parameters: domain_left_edge          = [0. 0. 0.]
yt : [INFO     ] 2019-10-11 17:09:04,702 Parameters: domain_right_edge         = [21. 21. 80.]
Comparing mpi-4-constrng-grid8/plt00010 to ../../field_configs/v4_mu_-0.100_w_0.100_Nx_21_nL_10_field_config.txt ...
 (Field values in AMReX plotfile rounded to 6 significant figures for the comparison)
   field     Abs Inf Norm   Rel Inf Norm   Abs L2 Norm    Rel L2 Norm
  phi_1_Re   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
  phi_1_Im   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
  phi_2_Re   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
  phi_2_Im   0.000000e+00   0.000000e+00   0.000000e+00   0.000000e+00
```

# Using PSC Bridges-2 GPU nodes

Each node has 8 NVIDIA V100 GPUs. Just like on Cori GPU, we use 1 MPI rank per GPU.

## Compiling

List of modules to load:

```
Currently Loaded Modules:
  1) anaconda3/2020.07   2) gcc/10.2.0   3) cuda/11.1.1   4) openmpi/3.1.6-gcc10.2.0   5) allocations/1.0   6) hdf5/1.10.7
```

Compile for GPUs with:

```
make -j 8 USE_MPI=TRUE USE_CUDA=TRUE USE_HDF5=TRUE
```

## Running

Uses the `gpu_visible.sh` and `bridges2.MPI.CUDA.slurm` scripts in `nrrb/amrex-langevin/Utils/bridges-2`

Both scripts are commented for details. The "GPU visible" script makes it
so each MPI rank can only see 1 GPU. The SLURM script is for submitting batch
jobs to the queue.

### Interactive Job

To get an interactive job on, e.g. 2 GPU nodes for 30 minutes:

```
interact -p GPU -t 00:30:00 -N 2 --ntasks-per-node=8 -n 16 --gres=gpu:8
```

The `--ntasks-per-node` and `--gres=gpu:8` options would generally stay the
same to give you access to all 8 GPUs on each node.

To run on a different number of nodes, change `-N [number-of-nodes]` and then
change `-n [number-of-nodes * ntasks-per-node]`.

Running in the interactive session:

```
mpirun ./gpu_visible.sh ./main3d.gnu.MPI.CUDA.ex inputs
```

### Batch Job

Use a SLURM script like `bridges2.MPI.CUDA.slurm` and pass all the job options to the
`sbatch` command like this:

```
sbatch -p GPU -t 00:30:00 -N 1 --ntasks-per-node=8 --gres=gpu:8 bridges2.MPI.CUDA.slurm
```

Requires the `gpu_visible.sh` script in the current directory.

To modify the job options:

- Change walltime & number of nodes:
-- `-t 00:30:00`: walltime in HH:MM:SS
-- `-N 1`: number of nodes

- Always use these options for the GPU nodes
-- `-p GPU`: request GPU nodes
-- `--ntasks-per-node=8`: 8 MPI ranks per node, 1 for each GPU
-- `--gres=gpu:8`: use all 8 GPUs on each node
