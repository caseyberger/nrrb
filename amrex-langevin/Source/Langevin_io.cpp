#include "Langevin.H"

using namespace amrex;

std::string MakeOutputFilename(const amrex::Geometry& geom, const NRRBParameters& nrrb, const int& nsteps) {
    // Domain geometry
    const auto domain_box = geom.Domain();
    const int length_x = domain_box.length(0);
    const int length_t = domain_box.length(AMREX_SPACEDIM-1);

    // Construct the output file suffix string using runtime parameters
    std::ostringstream output_stream;
    output_stream << "D_" << AMREX_SPACEDIM-1 << "_Nx_" << length_x << "_Nt_" << length_t;
    output_stream << "_dt_" << nrrb.dtau << "_nL_" << nsteps << "_eps_" << nrrb.eps;
    output_stream << "_m_" << nrrb.m << "_wtr_" <<nrrb.w_t;
    output_stream << "_wz_" << nrrb.w << "_l_" << nrrb.l << "_mu_" << nrrb.mu;

    // Construct the output file name
    std::string output_file = "output_" + output_stream.str();
    return output_file;
}

void CreateOutputFile(const std::string& filename) {
    if (ParallelDescriptor::IOProcessor())
    {
        using namespace ClassyHDF;

        // Create the HDF5 output file and raise an error if the file already exists,
        // to prevent unintentionally overwriting simulation data.
        try {
            File outputFile(filename + ".h5", FileMode::exists_is_error);
        } catch (...) {
            amrex::Error("File initialization error -- it may already exist");
        }
    }
}

void WritePlotfile(int Langevin_step, const Real Langevin_time,
                   const MultiFab& lattice, const Vector<std::string>& component_names,
                   const MultiFab& lattice_aux, const Vector<std::string>& auxiliary_names,
                   const Geometry& geom, const NRRBParameters& nrrb, const std::string& output_file)
{
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
#ifdef AMREX_USE_HDF5
    const Vector<std::string> write_group = {"Lattice", "step_" + std::to_string(Langevin_step)};
    WriteSingleLevelPlotfileHDF5(output_file, lattice_plot, plot_vars, geom, Langevin_time, Langevin_step, write_group);
#else
    const std::string pltfile = amrex::Concatenate("plt",Langevin_step,7);
    WriteSingleLevelPlotfile(pltfile, lattice_plot, plot_vars, geom, Langevin_time, Langevin_step);
#endif
}