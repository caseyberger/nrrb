#include <iomanip>
#include "AMReX_Print.H"
#include "Observables.H"
#include "Observables_Kernel.H"

using namespace amrex;

/*
Computes and prints Real and Imaginary parts of observables

don't forget to modify mu, m, w, wtr, and l by dtau if they appear in observables:
mu = dtau*mu; m = m/dtau; w = dtau*w; wtr = dtau* wtr; l = dtau*l;
*/
//void Equal_Time_Correlators(double *** Lattice, int size, int Nx, int Nt, std::string logfilename);

Observables::Observables(const amrex::Geometry& geom, const NRRBParameters& nrrb, const int& nsteps, bool restart)
{
	const auto domain_box = geom.Domain();
	const int length_x = domain_box.length(0);
	const int length_t = domain_box.length(AMREX_SPACEDIM-1);

	// Construct the logfile suffix string using runtime parameters
	std::ostringstream logfile_stream;
	logfile_stream << "D_" << AMREX_SPACEDIM-1 << "_Nx_" << length_x << "_Nt_" << length_t;
	logfile_stream << "_dt_" << nrrb.dtau << "_nL_" << nsteps << "_eps_" << nrrb.eps;
	logfile_stream << "_m_" << nrrb.m << "_wtr_" <<nrrb.w_t;
	logfile_stream << "_wz_" << nrrb.w << "_l_" << nrrb.l << "_mu_" << nrrb.mu << ".log";
	logfile_suffix = logfile_stream.str();

	// Construct the log file names
	observable_log_file = "logfile_" + logfile_suffix;

	// Initialize circulation radii and logfile names
	circulation.emplace_back(nrrb.circulation_radius_1, logfile_suffix);
	circulation.emplace_back(nrrb.circulation_radius_2, logfile_suffix);

	Print() << "logfile name = logfile_" << logfile_suffix << std::endl;

    // If we are not restarting from a checkpoint then initialize the log files
    // otherwise, assume the log files are already present.
    if (!restart)
    {
        initialize_files(geom);
    }
}

void Observables::initialize_files(const amrex::Geometry& geom)
{
    if (ParallelDescriptor::IOProcessor())
    {
	    std::ofstream obsFile;

		// Write observables file header
	    obsFile.open(observable_log_file, std::fstream::trunc);
        obsFile << "#step  ";
        obsFile << "Re[phi^{*}phi]      ";
        obsFile << "Im[phi^{*}phi]      ";
        obsFile << "Re[<n>]             ";
        obsFile << "Im[<n>]             ";
        obsFile << "Re[<Lz>]            ";
        obsFile << "Im[<Lz>]            ";
        obsFile << "Re[<S>]             ";
        obsFile << "Im[<S>]             ";
        obsFile << "dt (sec)" << std::endl;
        obsFile.close();

		// Write circulation log files
		for (auto& circ : circulation) circ.init_file(geom);
    }
}

void Observables::update(const int nL, const amrex::MultiFab& Lattice, const amrex::GeometryData& geom,
                         const NRRBParameters& nrrb_parm)
{
    const int Ncomp = Lattice.nComp();

    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_operations;
    ReduceData<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    // for (MFIter mfi(Lattice, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    for (MFIter mfi(Lattice, false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        const Array4<const Real>& L_obs = Lattice.array(mfi);

        reduce_operations.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (const int i, const int j, const int t) -> ReduceTuple
                {
                    // Evaluate observable contributions for this cell (i,j,t)
                    const auto observables = local_observables(i, j, t, L_obs, geom, nrrb_parm);

                    return {observables[Obs::PhiSqRe],
                            observables[Obs::PhiSqIm],
                            observables[Obs::DensRe],
                            observables[Obs::DensIm],
                            observables[Obs::LzRe],
                            observables[Obs::LzIm],
                            observables[Obs::SRe],
                            observables[Obs::SIm],
							observables[Obs::Circ1],
							observables[Obs::Circ2]};
                });
    }

    ReduceTuple reduced_observables = reduce_data.value();

    // MPI reduction to the IO Processor
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::PhiSqRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::PhiSqIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::DensRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::DensIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::LzRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::LzIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::SRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::SIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::Circ1>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::Circ2>(reduced_observables), IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
	    std::ofstream obsFile;

		// Write reduced observables
	    obsFile.open(observable_log_file, std::fstream::app);
        obsFile << std::setw(6)  << std::left << nL << ' ';
        obsFile << std::setw(19) << std::left << amrex::get<Obs::PhiSqRe>(reduced_observables) << ' ';
        obsFile << std::setw(19) << std::left << amrex::get<Obs::PhiSqIm>(reduced_observables) << ' ';
        obsFile << std::setw(19) << std::left << amrex::get<Obs::DensRe>(reduced_observables) << ' ';
        obsFile << std::setw(19) << std::left << amrex::get<Obs::DensIm>(reduced_observables) << ' ';
        obsFile << std::setw(19) << std::left << amrex::get<Obs::LzRe>(reduced_observables) << ' ';
        obsFile << std::setw(19) << std::left << amrex::get<Obs::LzIm>(reduced_observables) << ' ';
        obsFile << std::setw(19) << std::left << amrex::get<Obs::SRe>(reduced_observables) << ' ';
        obsFile << std::setw(19) << std::left << amrex::get<Obs::SIm>(reduced_observables) << ' ';
        obsFile << std::setw(11) << std::left << 0.0 << std::endl;
        obsFile.close();

		circulation[0].set_circulation(amrex::get<Obs::Circ1>(reduced_observables));
		circulation[1].set_circulation(amrex::get<Obs::Circ2>(reduced_observables));

		// Write reduced circulation
		for (auto& circ : circulation) circ.write();
    }
}
