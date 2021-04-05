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

Observables::Observables(const amrex::Geometry& geom, const amrex::DistributionMapping& dm, const amrex::BoxArray& ba, const NRRBParameters& nrrb, const int& nsteps)
{
	const auto domain_box = geom.Domain();
	const int length_x = domain_box.length(0);
	const int length_t = domain_box.length(AMREX_SPACEDIM-1);

	// Construct the logfile suffix string using runtime parameters
	std::ostringstream logfile_stream;
	logfile_stream << "D_" << AMREX_SPACEDIM-1 << "_Nx_" << length_x << "_Nt_" << length_t;
	logfile_stream << "_dt_" << nrrb.dtau << "_nL_" << nsteps << "_eps_" << nrrb.eps;
	logfile_stream << "_m_" << nrrb.m << "_wtr_" <<nrrb.w_t;
	logfile_stream << "_wz_" << nrrb.w << "_l_" << nrrb.l << "_mu_" << nrrb.mu;
    std::string logfile_suffix = logfile_stream.str();

	// Construct the log file names
	observable_log_file = "logfile_" + logfile_suffix + ".h5";

	// Initialize circulation radii and logfile names
	circulation.emplace_back(nrrb.circulation_radius_1);
	circulation.emplace_back(nrrb.circulation_radius_2);

	Print() << "logfile name = " << observable_log_file << std::endl;

	initialize_files(geom, nrrb, nsteps);

    // Setup our density profile particle container
    density_profile.Setup(geom, dm, ba, nrrb.profile_max_grid_size);
}

void Observables::initialize_files(const amrex::Geometry& geom, const NRRBParameters& nrrb, const int& nsteps)
{
    if (ParallelDescriptor::IOProcessor())
    {
        using namespace ClassyHDF;

        // Make an HDF5 observable log file
        File obsFile(observable_log_file, FileMode::trunc);

        // Write simulation parameters as attributes
        const auto domain_box = geom.Domain();
        const int length_x = domain_box.length(0);
        const int length_t = domain_box.length(AMREX_SPACEDIM-1);

        obsFile.attr<int>("SpatialDim", AMREX_SPACEDIM-1);
        obsFile.attr<int>("Nx", length_x);
        obsFile.attr<int>("Nt", length_t);
        obsFile.attr<int>("nL", nsteps);
        obsFile.attr<double>("dt", nrrb.dtau);
        obsFile.attr<double>("eps", nrrb.eps);
        obsFile.attr<double>("mass", nrrb.m);
        obsFile.attr<double>("w_trap", nrrb.w_t);
        obsFile.attr<double>("w_z", nrrb.w);
        obsFile.attr<double>("l", nrrb.l);
        obsFile.attr<double>("mu", nrrb.mu);

		// Create datasets in the observable file
        obsFile.create_dataset<int>("Steps");
        obsFile.create_dataset<Real>("LangevinTimes");

        auto make_obs_re_im = [&] (const std::string& obs) {
            Group ogroup = obsFile.get_group(obs);
            ogroup.create_dataset<Real>("Re");
            ogroup.create_dataset<Real>("Im");
        };

        make_obs_re_im("PhiSq");
        make_obs_re_im("Dens");
        make_obs_re_im("Lz");
        make_obs_re_im("S");
        make_obs_re_im("KE");
        make_obs_re_im("Vtr");
        make_obs_re_im("Vint");
        make_obs_re_im("S_tau");
        make_obs_re_im("S_del");
        make_obs_re_im("S_trap");
        make_obs_re_im("S_w");
        make_obs_re_im("S_int");

		// Create circulation datasets
		for (auto& circ : circulation) {
            circ.init_file(observable_log_file, geom);
        }
    }
}

void Observables::update(const int nL, const Real Ltime, const amrex::MultiFab& Lattice,
                         const amrex::GeometryData& geom, const NRRBParameters& nrrb_parm)
{
    const int Ncomp = Lattice.nComp();

    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, 
              ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum,ReduceOpSum, ReduceOpSum,ReduceOpSum, ReduceOpSum,
              ReduceOpSum, ReduceOpSum,ReduceOpSum, ReduceOpSum,ReduceOpSum, ReduceOpSum,ReduceOpSum, ReduceOpSum,
              ReduceOpSum, ReduceOpSum> reduce_operations;
    ReduceData<Real, Real, Real, Real, Real, Real, Real, Real, Real, Real,Real, Real, Real, 
               Real,Real, Real, Real, Real,Real, Real, Real, Real,Real, Real, Real, Real> reduce_data(reduce_operations);
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
                            observables[Obs::KERe],
                            observables[Obs::KEIm],
                            observables[Obs::VtrRe],
                            observables[Obs::VtrIm],
                            observables[Obs::VintRe],
                            observables[Obs::VintIm],
                            observables[Obs::StauRe],
                            observables[Obs::StauIm],
                            observables[Obs::SdelRe],
                            observables[Obs::SdelIm],
                            observables[Obs::StrapRe],
                            observables[Obs::StrapIm],
                            observables[Obs::SwRe],
                            observables[Obs::SwIm],
                            observables[Obs::SintRe],
                            observables[Obs::SintIm],
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
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::KERe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::KEIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::VtrRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::VtrIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::VintRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::VintIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::StauRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::StauIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::SdelRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::SdelIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::StrapRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::StrapIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::SwRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::SwIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::SintRe>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::SintIm>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::Circ1>(reduced_observables), IOProc);
    ParallelDescriptor::ReduceRealSum(amrex::get<Obs::Circ2>(reduced_observables), IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        using namespace ClassyHDF;

        // Write reduced observables
        File obsFile(observable_log_file);

        obsFile.append(Data<int>("Steps", {nL}));
        obsFile.append(Data<Real>("LangevinTimes", {Ltime}));

        auto append_obs_re_im = [&] (const std::string& obs,
                                     const Real& oRe, const Real& oIm) {
            Group ogroup = obsFile.get_group(obs);
            ogroup.append(Data<Real>("Re", {oRe}));
            ogroup.append(Data<Real>("Im", {oIm}));
        };

        append_obs_re_im("PhiSq" , amrex::get<Obs::PhiSqRe>(reduced_observables),
                                   amrex::get<Obs::PhiSqIm>(reduced_observables));
        append_obs_re_im("Dens"  , amrex::get<Obs::DensRe>(reduced_observables),
                                   amrex::get<Obs::DensIm>(reduced_observables));
        append_obs_re_im("Lz"    , amrex::get<Obs::LzRe>(reduced_observables),
                                   amrex::get<Obs::LzIm>(reduced_observables));
        append_obs_re_im("S"     , amrex::get<Obs::SRe>(reduced_observables),
                                   amrex::get<Obs::SIm>(reduced_observables));
        append_obs_re_im("KE"    , amrex::get<Obs::KERe>(reduced_observables),
                                   amrex::get<Obs::KEIm>(reduced_observables));
        append_obs_re_im("Vtr"   , amrex::get<Obs::VtrRe>(reduced_observables),
                                   amrex::get<Obs::VtrIm>(reduced_observables));
        append_obs_re_im("Vint"  , amrex::get<Obs::VintRe>(reduced_observables),
                                   amrex::get<Obs::VintIm>(reduced_observables));
        append_obs_re_im("S_tau" , amrex::get<Obs::StauRe>(reduced_observables),
                                   amrex::get<Obs::StauIm>(reduced_observables));
        append_obs_re_im("S_del" , amrex::get<Obs::SdelRe>(reduced_observables),
                                   amrex::get<Obs::SdelIm>(reduced_observables));
        append_obs_re_im("S_trap", amrex::get<Obs::StrapRe>(reduced_observables),
                                   amrex::get<Obs::StrapIm>(reduced_observables));
        append_obs_re_im("S_w"   , amrex::get<Obs::SwRe>(reduced_observables),
                                   amrex::get<Obs::SwIm>(reduced_observables));
        append_obs_re_im("S_int" , amrex::get<Obs::SintRe>(reduced_observables),
                                   amrex::get<Obs::SintIm>(reduced_observables));

		circulation[0].set_circulation(amrex::get<Obs::Circ1>(reduced_observables));
		circulation[1].set_circulation(amrex::get<Obs::Circ2>(reduced_observables));

		// Write reduced circulation
		for (auto& circ : circulation) {
            circ.write(observable_log_file);
        }
    }

    // Calculate and save the density profile
    density_profile.AccumulateProfile(Lattice, nrrb_parm);
    density_profile.Write(nL, Ltime, nrrb_parm);
}
