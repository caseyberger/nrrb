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
