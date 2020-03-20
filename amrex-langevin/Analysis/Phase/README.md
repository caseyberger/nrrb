# Field phase analysis

These jupyter notebooks demonstrate how the fixed phase initialization works
and different ways of numerically calculating the field phase from the real
components of the field.

This will generate the plotfile for the `Langevin-Phase-Unit-Circulation.ipynb` notebook:

```
./main3d.gnu.ex inputs_mpi_omp_test nrrb.problem_type=1 nsteps=0
```

This will generate the plotfile for the `Langevin-Phase-Circulation-Cancellation.ipynb` notebook:

```
./main3d.gnu.ex inputs_mpi_omp_test nrrb.problem_type=2 nsteps=0
```

This will generate the plotfile for the `Langevin-Phase-Negative-Unit-Circulation.ipynb` notebook:

```
./main3d.gnu.ex inputs_mpi_omp_test nrrb.problem_type=3 nsteps=0
```

In these analysis notebooks, those plotfiles were manually organized into
`problem_1`, `problem_2`, and `problem_3` directories.

The circulation from each of those plotfiles is shown in the `Langevin-Circulation-Initialization.ipynb` notebook.

These notebooks demonstrate a way of calculating phase and circulation that
gives 1, 0, and -1 for the respective initializations of `nrrb.problem_type=1`,
`nrrb.problem_type=2`, and `nrrb.problem_type=3`.
