#include "Langevin.H"

using namespace amrex;

void WritePlotfile(int Langevin_step, const Real Langevin_time,
                   const MultiFab& lattice, const Vector<std::string>& component_names,
                   const MultiFab& lattice_aux, const Vector<std::string>& auxiliary_names,
                   const Geometry& geom, const NRRBParameters& nrrb)
{
    const std::string& pltfile = amrex::Concatenate("plt",Langevin_step,5);

    // consolidate field variable names and auxiliary variable names into a single vector
    Vector<std::string> plot_vars;
    for (auto s : component_names) plot_vars.push_back(s);
    for (auto s : auxiliary_names) plot_vars.push_back(s);

    // consolidate phi and auxiliary variables into one multifab for I/O with 0 ghost cells
    MultiFab lattice_plot(lattice.boxArray(), lattice.DistributionMap(),
                          lattice.nComp() + AIdx::NAux, 0);

    // copy all 4 phi variables to the first nComp() = 4 components of lattice_plot
    MultiFab::Copy(lattice_plot, lattice, 0, 0, lattice.nComp(), 0);

    // copy all auxiliary variables to the following components of lattice_plot
    MultiFab::Copy(lattice_plot, lattice_aux, 0, lattice.nComp(), AIdx::NAux, 0);

    // write lattice_plot to a plotfile
    if (nrrb.use_hdf5) {
        WriteSingleLevelPlotfileHDF5(pltfile, lattice_plot, plot_vars, geom, Langevin_time, 0);
    } else {
        WriteSingleLevelPlotfile(pltfile, lattice_plot, plot_vars, geom, Langevin_time, 0);
    }
}