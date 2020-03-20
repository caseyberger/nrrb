#include <iomanip>
#include "AMReX_Print.H"
#include "Observables.H"

using namespace amrex;

/*
Computes and prints Real and Imaginary parts of observables

don't forget to modify mu, m, w, wtr, and l by dtau if they appear in observables:
mu = dtau*mu; m = m/dtau; w = dtau*w; wtr = dtau* wtr; l = dtau*l;
*/
//void Equal_Time_Correlators(double *** Lattice, int size, int Nx, int Nt, std::string logfilename);

Observables::Observables(const amrex::Geometry& geom, const NRRBParameters& nrrb, const int& nsteps)
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

	initialize_files(geom);
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

Real Circulation(amrex::Array4<const amrex::Real> const& Lattice, const amrex::Box& box,
				 const amrex::GeometryData& geom, int radius) {
	// Box and Domain geometry
	const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);

	// For this to work properly, the following must be true:
	// - the domain sizes in x and y must be odd
	// - the center of the domain must satisfy the following:
    const int x_center = 0.5 * (domain_lo.x + domain_hi.x);
    const int y_center = 0.5 * (domain_lo.y + domain_hi.y);

	// Loop left/right edges in the x-dimension
	const int i_left = x_center - radius;
    const int i_right = x_center + radius;

	// Loop bottom/top edges in the y-dimension
    const int j_bottom = y_center - radius;
    const int j_top = y_center + radius;

	// Number of points in the t-dimension in the domain
	const long Nt = domain_box.length(AMREX_SPACEDIM-1);

	// Set t coordinate where we want to calculate circulation
	const int t = Nt/2;

	// Return 0.0 for this box if t is not inside it
	if (t < lo.z || t > hi.z)
	{
		return 0.0;
	}

	// Initialize circulation sum for this part of the domain
	Real circ_sum = 0.0;

	auto Theta = [&](int i, int j, int k) -> Real {
		// tan(theta) = (phi_1_Im + phi_2_Re)/(phi_1_Re - phi_2_Im);
		Real theta = std::atan2(Lattice(i,j,k,Field(1,C::Im)) + Lattice(i,j,k,Field(2,C::Re)),
					            Lattice(i,j,k,Field(1,C::Re)) - Lattice(i,j,k,Field(2,C::Im)));

		// return theta in the range [0, 2*pi]
		if (theta < 0.0)
			theta = theta + Constants::TwoPi;

		return theta;
	};

    auto DeltaTheta = [&](int ilp, int jlp, int klp, int il, int jl, int kl) -> Real {
        // Returns Theta(ilp, jlp, klp) - Theta(il, jl, kl) = theta_l+1 - theta_l
        Real dtheta = Theta(ilp, jlp, klp) - Theta(il, jl, kl);

        // Put dtheta into the range [-pi, 0] or [0, pi]
        if (dtheta < -Constants::Pi)
            dtheta +=  Constants::TwoPi;
        else if (dtheta > Constants::Pi)
            dtheta += -Constants::TwoPi;

        return dtheta;
    };

	// We are summing contributions from theta_t_l+1 - theta_t_l
	// where l denotes (x, y) for a loop lattice site
	// and l+1 denotes (x, y) for the next lattice site on the loop
	// chosen by a loop traversal in the -z direction using the right-hand-rule.
	//
	// The loop starts from the point S, around the center of the domain C:
	//
	// ^ y
	// |      ^------->
	// |      |       |
	// |      |   C   |
	// |      |       |
	// |      S<------v
	// |
	// .-----------> x

	// Loop over y-dimension to add contributions from left/right edges
	// of the loop contained in this box.
	//
	// This loop includes corners, so we do not double count them in the next loop over x
	for (int j = lo.y; j <= hi.y; ++j) {
		// if left cell at this y is within the box, add its phase difference
		if (i_left >= lo.x && i_left <= hi.x) {
			if (j >= j_bottom && j < j_top) {
				circ_sum += DeltaTheta(i_left, j+1, t, i_left, j, t);
			}
			else if (j == j_top) {
				circ_sum += DeltaTheta(i_left+1, j, t, i_left, j, t);
			}
		}

		// if right cell at this y is within the box, add its phase difference
		if (i_right >= lo.x && i_right <= hi.x) {
			if (j > j_bottom && j <= j_top) {
				circ_sum += DeltaTheta(i_right, j-1, t, i_right, j, t);
			}
			else if (j == j_bottom) {
				circ_sum += DeltaTheta(i_right-1, j, t, i_right, j, t);
			}
		}
	}

	// Loop over x-dimension to add contributions from bottom/top edges
	// of the loop contained in this box.
	for (int i = lo.x; i <= hi.x; ++i) {
		// In the following `if` statement, we do not include corners
		// (i.e. we use `>` and `<` instead of `>=` and `<=`)
		// because the above loop over y already accounted for them.
		if (i > i_left && i < i_right) {
			// if top cell at this x is within the box, add its phase difference
			if (j_top >= lo.y && j_top <= hi.y) {
				circ_sum += DeltaTheta(i+1, j_top, t, i, j_top, t);
			}

			// if bottom cell at this x is within the box, add its phase difference
			if (j_bottom >= lo.y && j_bottom <= hi.y) {
				circ_sum += DeltaTheta(i-1, j_bottom, t, i, j_bottom, t);
			}
		}
	}

	// adding the circulation for one loop and one box to the total circulation for this loop
	Real circ = circ_sum / Constants::TwoPi;
	return circ;
}
