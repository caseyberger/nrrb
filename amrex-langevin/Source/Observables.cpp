#include <iomanip>
#include "AMReX_Print.H"
#include "Langevin.H"

using namespace amrex;

/*
Computes and prints Real and Imaginary parts of observables

don't forget to modify mu, m, w, wtr, and l by dtau if they appear in observables:
mu = dtau*mu; m = m/dtau; w = dtau*w; wtr = dtau* wtr; l = dtau*l;
*/
//void Equal_Time_Correlators(double *** Lattice, int size, int Nx, int Nt, std::string logfilename);
Real S_tau_Re(int i,int j,int t,int a, Real mu, amrex::Array4<const amrex::Real> const& Lattice);
Real S_tau_Im(int i,int j,int t,int a, Real mu, amrex::Array4<const amrex::Real> const& Lattice);
Real S_del_Re(int i,int j,int t,int a, Real m, amrex::Array4<const amrex::Real> const& Lattice,const amrex::GeometryData& geom);
Real S_del_Im(int i,int j,int t,int a, Real m, amrex::Array4<const amrex::Real> const& Lattice,const amrex::GeometryData& geom);
Real S_trap_Re(int i,int j,int t,int a, Real w_t, const Real r2, amrex::Array4<const amrex::Real> const& Lattice);
Real S_trap_Im(int i,int j,int t,int a, Real w_t, const Real r2, amrex::Array4<const amrex::Real> const& Lattice);
Real S_w_Re(int i,int j,int t,int a, Real w, const Real x,const Real y, amrex::Array4<const amrex::Real> const& Lattice);
Real S_w_Im(int i,int j,int t,int a, Real w, const Real x,const Real y, amrex::Array4<const amrex::Real> const& Lattice);
Real S_int_Re(int i,int j,int t,int a, Real l, amrex::Array4<const amrex::Real> const& Lattice);
Real S_int_Im(int i,int j,int t,int a, Real l, amrex::Array4<const amrex::Real> const& Lattice);
void Circulation(amrex::Array4<const amrex::Real> const& Lattice, const amrex::Box& box, const amrex::GeometryData& geom,
			long int Nt, int radius, std::string filename);

void initialize_observables(const std::string observable_log_file)
{
    // Write observables header file
    if (ParallelDescriptor::IOProcessor())
    {
	    std::ofstream obsFile;
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
    }
}

void update_observables(const int nL, const amrex::MultiFab& Lattice, const amrex::GeometryData& geom,
                        const NRRBParameters& nrrb_parm, const std::string observable_log_file)
{
    const int Ncomp = Lattice.nComp();

    ReduceOps<ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum, ReduceOpSum> reduce_operations;
    ReduceData<Real, Real, Real, Real, Real, Real, Real, Real> reduce_data(reduce_operations);
    using ReduceTuple = typename decltype(reduce_data)::Type;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(Lattice); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        const Array4<const Real>& L_obs = Lattice.array(mfi);

        reduce_operations.eval(bx, reduce_data,
                [=] (Box const& bx) -> ReduceTuple
                {
                    const auto observables = compute_observables(bx, Ncomp, L_obs, geom, nrrb_parm);
                    return {observables[Obs::PhiSqRe],
                            observables[Obs::PhiSqIm],
                            observables[Obs::DensRe],
                            observables[Obs::DensIm],
                            observables[Obs::LzRe],
                            observables[Obs::LzIm],
                            observables[Obs::SRe],
                            observables[Obs::SIm]};
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

    // Write reduced observables
    if (ParallelDescriptor::IOProcessor())
    {
	    std::ofstream obsFile;
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
    }
}

amrex::Vector<amrex::Real> compute_observables(const amrex::Box& box, const int Ncomp,
                                               const amrex::Array4<const amrex::Real>& Lattice,
                                               const amrex::GeometryData& geom,
                                               const NRRBParameters nrrb_parm){
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
	long int Nt = box.length(2);

	//Initialize observable variables
    amrex::Vector<amrex::Real> observables(Obs::NumObservables, 0.0);
	Real dens_Re = 0;
	Real dens_Im = 0;
	Real phisq_Re = 0.;
	Real phisq_Im = 0.;
	Real Lz_Re = 0.;
	Real Lz_Im = 0.;
    Real S_Re = 0.;
    Real S_Im = 0.;

	//modify parameters by dtau
	const Real mu = nrrb_parm.dtau * nrrb_parm.mu;
	const Real m = nrrb_parm.m / nrrb_parm.dtau;
	const Real w_t = nrrb_parm.dtau * nrrb_parm.w_t;
	const Real w = nrrb_parm.dtau * nrrb_parm.w;
	const Real l = nrrb_parm.dtau * nrrb_parm.l;

	//define spacing parameters (x, y r^2, etc)
	const auto domain_xlo = geom.ProbLo();
    const auto domain_xhi = geom.ProbHi();

    const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);
	long int domain_volume = domain_box.volume();

    const Real x_center = 0.5 * (domain_xlo[0] + domain_xhi[0]);
    const Real y_center = 0.5 * (domain_xlo[1] + domain_xhi[1]);
    const Real dx_cell = geom.CellSize(0);
    const Real dy_cell = geom.CellSize(1);

	//loop over lattice
	for (int j = lo.y; j <= hi.y; ++j){
		for (int i = lo.x; i <= hi.x; ++i){
            const Real x_cell = domain_xlo[0] + dx_cell * (i + 0.5 - domain_lo.x);
            const Real y_cell = domain_xlo[1] + dy_cell * (j + 0.5 - domain_lo.y);
            const Real r2 = std::pow(x_cell - x_center, 2) + std::pow(y_cell - y_center, 2);

            //density profiles are averaged over time loop only
			Real dp_Re = 0.;
			Real dp_Im = 0.;
            for (int t = lo.z; t <= hi.z; ++t){
                for (int a = 1; a <= 2; a++){
                    //Field modulus squared
                    phisq_Re += 0.5 * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
                    phisq_Re += -0.5 * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
                    phisq_Im += Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
                    for (int b = 1; b <= 2; b++){
						//Density
						dp_Re += delta(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
						dp_Re += -delta(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
						dp_Re += -epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
						dp_Re += -epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
						dp_Im += delta(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
						dp_Im += delta(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
						dp_Im += epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
						dp_Im += -epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
					}//loop over b = 1,2
				}//loop over a = 1,2
				//Angular momentum
				//#if (AMREX_SPACEDIM == 3)
				if (AMREX_SPACEDIM == 3){
					//define x and y
					const Real x = x_cell - x_center;
	    			const Real y = y_cell - y_center;
					for (int a = 1; a <= 2; a++){//loop over a
						//real sum over a only
						Lz_Re += 0.5 * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(a,C::Im));
						Lz_Re += 0.5 * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(a,C::Re));
						Lz_Re += -0.5 * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(a,C::Im));
						Lz_Re += -0.5 * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(a,C::Re));
						Lz_Re += y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
						Lz_Re += -x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
						//imaginary sum over a only
						Lz_Im += -0.5 * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(a,C::Re));
						Lz_Im += 0.5 * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(a,C::Im));
						Lz_Im += 0.5 * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(a,C::Re));
						Lz_Im += -0.5 * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(a,C::Im));
						Lz_Im += -0.5 * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
						Lz_Im += 0.5 * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
						Lz_Im += 0.5 * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
						Lz_Im += -0.5 * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
						for (int b = 1; b <= 2; b++){//loop over b
							//real sum over a and b
							Lz_Re += 0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Re));
							Lz_Re += -0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Im));
							Lz_Re += 0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Im));
							Lz_Re += -0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Re));
							//imaginary sum over a and b
							Lz_Im += -0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Im));
							Lz_Im += -0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Re));
							Lz_Im += 0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Im));
							Lz_Im += 0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Re));
						}//loop over b for Lz
					}//loop over a for Lz
				}//if dim+1 = 3
				//#endif
				//ACTION
				for (int a = 1; a <= 2; a++){
					//S_tau
					S_Re += S_tau_Re(i,j,t,a,mu,Lattice);
					S_Im += S_tau_Im(i,j,t,a,mu,Lattice);
					//S_del
					S_Re += S_del_Re(i,j,t,a,m,Lattice,geom);
					S_Im += S_del_Im(i,j,t,a,m,Lattice,geom);
					//S_trap
					S_Re += S_trap_Re(i,j,t,a,w_t,r2,Lattice);
					S_Im += S_trap_Im(i,j,t,a,w_t,r2,Lattice);
					//S_w
					//#if (AMREX_SPACEDIM == 3)
					if (AMREX_SPACEDIM == 3){
						const Real x = x_cell - x_center;
	    				const Real y = y_cell - y_center;
						S_Re += S_w_Re(i,j,t,a,w,x,y,Lattice);
						S_Im += S_w_Im(i,j,t,a,w,x,y,Lattice);
					}
					//#endif
					//S_int
					S_Re += S_int_Re(i,j,t,a,l,Lattice);
					S_Im += S_int_Im(i,j,t,a,l,Lattice);
				}//loop over a (for the Action)
            }//loop over t
			dp_Re = 0.5 * exp(mu) * dp_Re;
			dp_Im = 0.5 * exp(mu) * dp_Im;
			//add local density to avg density variables
			dens_Re += dp_Re;
			dens_Im += dp_Im;
		}//loop over x
	}//loop over y

	//std::cout << "Phi^{*}phi successfully computed: " << phisq[0] << " + i" << phisq[1] << std::endl;
	//Equal_Time_Correlators(Lattice, size, Nx, Nt, filename);

	//std::cout << "Lz successfully computed: Lz = " << Lz[0] << " + i" << Lz[1] << std::endl;
	//compute circulation
	/* long int Nx = box.length(0); */
	/* Circulation(Lattice, box, geom, Nt, Nx/4, filename); */
	/* Circulation(Lattice, box, geom, Nt, Nx/2 - 1, filename); */

    observables[Obs::PhiSqRe] = phisq_Re / domain_volume;
    observables[Obs::PhiSqIm] = phisq_Im / domain_volume;
    observables[Obs::DensRe] = dens_Re / domain_volume;
    observables[Obs::DensIm] = dens_Im / domain_volume;
    observables[Obs::LzRe] = Lz_Re / domain_volume;
    observables[Obs::LzIm] = Lz_Im / domain_volume;
    observables[Obs::SRe] = S_Re / domain_volume;
    observables[Obs::SIm] = S_Im / domain_volume;

    return observables;
}

double S_tau_Re(int i,int j,int t,int a, double mu,amrex::Array4<const amrex::Real> const& Lattice){
	double S_Re = 0.;
	S_Re += 0.5 * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Re)) ;
	S_Re += -0.5 * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Im)) ;
	S_Re += -0.5 * exp(mu) * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) ;
	S_Re += 0.5 * exp(mu) * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) ;
	for (int b=1; b<=2; b++){
		S_Re += 0.5 * exp(mu) * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im)) ;
		S_Re += 0.5 * exp(mu) * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re)) ;
	}
	return S_Re;
}

double S_tau_Im(int i,int j,int t,int a,double mu, amrex::Array4<const amrex::Real> const& Lattice){
	double S_Im = 0.;
	S_Im += Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Im)) ;
	S_Im +=-0.5 * exp(mu) * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) ;
	S_Im +=-0.5 * exp(mu) * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) ;
	for (int b=1; b<=2; b++){
		S_Im += -0.5 * exp(mu) * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) ;
		S_Im += 0.5 * exp(mu) * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) ;
	}
	return S_Im;
}

double S_del_Re(int i,int j,int t,int a,double m, amrex::Array4<const amrex::Real> const& Lattice,const amrex::GeometryData& geom){
	const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);
    double S_Re = 0.;
	S_Re += 0.5*m*(AMREX_SPACEDIM -1)*Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Re)) ;
	S_Re += -0.5*m*(AMREX_SPACEDIM -1)*Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Im)) ;
	for (int d = 1; d <= AMREX_SPACEDIM-1; d++) { //loop over adjacent spatial sites!
		//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
		// -(1/2m)\phi_{a,i+1}^{R} - (1/2m)\phi_{a,ii1}^{R}
        if (d == 1)
        {
            if (!(i == domain_lo.x || i == domain_hi.x))
            {
				S_Re += -0.25 * m * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i+1,j,t,Field(a,C::Re));
				S_Re += 0.25 * m * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i+1,j,t,Field(a,C::Im));
				S_Re += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t,Field(a,C::Re));
				S_Re += 0.25 * m * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i-1,j,t,Field(a,C::Im));;
				for (int b=1; b<=2; b++){
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i+1,j,t,Field(b,C::Im));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i+1,j,t,Field(b,C::Re));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i-1,j,t,Field(b,C::Im));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i-1,j,t,Field(b,C::Re));
				}//loop over b
            }//checking x boundaries
        }//loop over x dim
        else if (d == 2)
        {
            if (!(j == domain_lo.y || j == domain_hi.y))
            {
				S_Re += -0.25 * m * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i,j+1,t,Field(a,C::Re));
				S_Re += 0.25 * m * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i,j+1,t,Field(a,C::Im));
				S_Re += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t,Field(a,C::Re));
				S_Re += 0.25 * m * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i,j-1,t,Field(a,C::Im));;
				for (int b=1; b<=2; b++){
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i,j+1,t,Field(b,C::Im));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i,j+1,t,Field(b,C::Re));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i,j-1,t,Field(b,C::Im));
					S_Re += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i,j-1,t,Field(b,C::Re));
				}//loop over b
            }//checking y boundaries
        }//loop over y dim
	}//derivative loop
	return S_Re;
}

double S_del_Im(int i,int j,int t,int a, double m, amrex::Array4<const amrex::Real> const& Lattice,const amrex::GeometryData& geom){
	const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);
    double S_Im = 0.;
	S_Im += m*(AMREX_SPACEDIM -1)*Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Im)) ;
	for (int d = 1; d <= AMREX_SPACEDIM-1; d++) { //loop over adjacent spatial sites!
		//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
		// -(1/2m)\phi_{a,i+1}^{R} - (1/2m)\phi_{a,ii1}^{R}
		if (d == 1)
        {
            if (!(i == domain_lo.x || i == domain_hi.x))
            {
            	S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i+1,j,t,Field(a,C::Im));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i+1,j,t,Field(a,C::Re));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t,Field(a,C::Im));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t,Field(a,C::Re));
				for (int b=1; b<=2; b++){
					S_Im += -0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i+1,j,t,Field(b,C::Re));
					S_Im += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i+1,j,t,Field(b,C::Im));
					S_Im += -0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i-1,j,t,Field(b,C::Re));
					S_Im += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i-1,j,t,Field(b,C::Im));
				}//loop over b
            }//checking x boundaries
        }//loop over x dim
        else if (d == 2)
        {
            if (!(j == domain_lo.y || j == domain_hi.y))
            {
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i,j+1,t,Field(a,C::Im));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i,j+1,t,Field(a,C::Re));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t,Field(a,C::Im));
				S_Im += -0.25 * m * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t,Field(a,C::Re));
				for (int b=1; b<=2; b++){
					S_Im += -0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i,j+1,t,Field(b,C::Re));
					S_Im += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i,j+1,t,Field(b,C::Im));
					S_Im += -0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re))  * Lattice(i,j-1,t,Field(b,C::Re));
					S_Im += 0.25 * m * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im))  * Lattice(i,j-1,t,Field(b,C::Im));
				}//loop over b
            }//checking y boundaries
        }//loop over y dim
	}//derivative loop
	return S_Im;
}

double S_trap_Re(int i,int j,int t,int a, double w_t,const Real r2, amrex::Array4<const amrex::Real> const& Lattice){
	double S_Re = 0.;
	S_Re += 0.25 * w_t * w_t * r2 * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) ;
	S_Re += -0.25 * w_t * w_t * r2 * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) ;
	for (int b=1; b<=2; b++){
		S_Re += -0.25 * w_t * w_t * r2 * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) ;
		S_Re += -0.25 * w_t * w_t * r2 * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) ;
	}//loop over b
	return S_Re;
}

double S_trap_Im(int i,int j,int t,int a,double w_t,const Real r2, amrex::Array4<const amrex::Real> const& Lattice){
	double S_Im = 0.;
	S_Im += 0.25 * w_t * w_t * r2 * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) ;
	S_Im += 0.25 * w_t * w_t * r2 * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) ;
	for (int b=1; b<=2; b++){
		S_Im += 0.25 * w_t * w_t * r2 * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) ;
		S_Im += -0.25 * w_t * w_t * r2 * epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) ;
	}//loop over b
	return S_Im;
}

double S_w_Re(int i,int j,int t,int a,double w,const Real x,const Real y, amrex::Array4<const amrex::Real> const& Lattice){
	double S_Re = 0.;
	for (int b = 1; b<=2; b++){
		S_Re += 0.5 * w * epsilon(a,b) * (x - y) * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) ;
		S_Re += -0.5 * w * epsilon(a,b) * (x - y) * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) ;
		S_Re += -0.5 * w * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j-1,t-1,Field(b,C::Re)) ;
		S_Re += 0.5 * w * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j-1,t-1,Field(b,C::Im)) ;
		S_Re += 0.5 * w * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Re));
		S_Re += -0.5 * w * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i-1,j,t-1,Field(b,C::Im)) ;
		S_Re += 0.5 * w * (x - y) * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) ;
		S_Re += 0.5 * w * (x - y) * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) ;
		S_Re += -0.5 * w * x * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j-1,t-1,Field(b,C::Im)) ;
		S_Re += -0.5 * w * x * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j-1,t-1,Field(b,C::Re)) ;
		S_Re += 0.5* w * y * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i-1,j,t-1,Field(b,C::Im));
		S_Re += 0.5* w * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Re));
	}

	return S_Re;
}

double S_w_Im(int i,int j,int t,int a,double w,const Real x,const Real y, amrex::Array4<const amrex::Real> const& Lattice){
	double S_Im = 0.;
	for (int b = 1; b<=2; b++){
		S_Im += 0.5 * w * epsilon(a,b) * (x-y) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
		S_Im += 0.5 * w * epsilon(a,b) * (x-y) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Im += -0.5 * w * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Im));
		S_Im += -0.5 * w * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j-1,t-1,Field(b,C::Re));
		S_Im += 0.5 * w * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Im));
		S_Im += 0.5 * w * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Re));
		S_Im += 0.5 * w * (y-x) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
		S_Im += -0.5 * w * (y-x) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
		S_Im += 0.5 * w * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Re));
		S_Im += -0.5 * w * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Im));
		S_Im += -0.5 * w * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Re));
		S_Im += 0.5 * w * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Im));
	}
	return S_Im;
}
							
double S_int_Re(int i,int j,int t,int a,double l, amrex::Array4<const amrex::Real> const& Lattice){
	double S_Re = 0.;
	for (int b=1; b<=2; b++){
		S_Re += 0.25 * l * Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re));
		S_Re += -0.25 * l * Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re));
		S_Re += 0.25*l*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) -Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Re += 0.5*l*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im))  - Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Re += 0.5*l*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) -Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Re += -0.5*l*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) +Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Re += -0.5*l*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) +Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Re += l*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Re += 0.25*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Im))  - Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Im)) );
		S_Re += 0.25*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Re))  - Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Re)) );
		S_Re += 0.25*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im))  - Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Re += 0.25*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re))  - Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Re += 0.5*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Re))  - Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Im)) );
		S_Re += 0.5*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im))  - Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) );
	}//loop over b
	return S_Re;
}	

double S_int_Im(int i,int j,int t,int a,double l, amrex::Array4<const amrex::Real> const& Lattice){
	double S_Im = 0.;
	for (int b=1; b<=2; b++){
		//S_int
		S_Im += 0.5*l*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) -Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Im += 0.5*l*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) -Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.5*l*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) -Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.5*l*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) -Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Im += 0.5*l*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) -Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.5*l*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) -Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.25*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re))  - Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) );
		S_Im += 0.25*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Re))  - Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Re)) );
		S_Im += 0.25*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Im))  - Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t,Field(b,C::Im)) );
		S_Im += 0.25*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Im)) -Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Re)) );
		S_Im += 0.5*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Re))  + Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t,Field(b,C::Im)) );
		S_Im += -0.5*l*epsilon(a,b)*(Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Re)) *Lattice(i,j,t-1,Field(b,C::Im))  + Lattice(i,j,t,Field(a,C::Im)) *Lattice(i,j,t,Field(a,C::Re)) *Lattice(i,j,t-1,Field(a,C::Im)) *Lattice(i,j,t-1,Field(b,C::Re)) );
	}//loop over b
	return S_Im;
}					


/*void Equal_Time_Correlators(double *** Lattice, int size, int Nx, int Nt, std::string logfilename){
	//compute the equal time correlator for each point on the lattice to each other point
	bool new_lat = true;
	int t = Nt/2;
	std::ofstream Gfile;
	std::string G_filename = "G_"+logfilename.substr(8);
	Gfile.open(G_filename, std::fstream::app);
	std::ofstream G2file;
	std::string G2_filename = "G2_"+logfilename.substr(8);
	G2file.open(G2_filename, std::fstream::app);
	//std::cout << "t = " << t << std::endl;
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){	
			double GijRe = 0.;
			double GijIm = 0.;
			double G2ijRe = 0.;
			double G2ijIm = 0.;
			for (int a = 1; a <= 2; a++){
				std::vector<double> phia_i = phi_a(Lattice, i, t, a,new_lat);
				for (int b = 1; b <=2; b++){
					std::vector<double> phib_j = phi_a(Lattice, j, t, b, new_lat);
					GijRe += 0.5*delta(a,b)*phia_i[0]*phib_j[0]-0.5*delta(a,b)*phia_i[1]*phib_j[1];
					GijRe += -0.5*epsilon(a,b)*phia_i[0]*phib_j[1]-0.5*epsilon(a,b)*phia_i[1]*phib_j[0];
					GijIm += 0.5*delta(a,b)*phia_i[0]*phib_j[1]+0.5*delta(a,b)*phia_i[1]*phib_j[0];
					GijIm += 0.5*epsilon(a,b)*phia_i[0]*phib_j[0]-0.5*epsilon(a,b)*phia_i[1]*phib_j[1];
					G2ijRe += 0.25*phia_i[0]*phia_i[0]*phib_j[0]*phib_j[0];
					G2ijRe -= 0.25*phia_i[0]*phia_i[0]*phib_j[1]*phib_j[1];
					G2ijRe -= 0.25*phia_i[1]*phia_i[1]*phib_j[0]*phib_j[0];
					G2ijRe += 0.25*phia_i[1]*phia_i[1]*phib_j[1]*phib_j[1];
					G2ijRe += phia_i[0]*phia_i[1]*phib_j[0]*phib_j[1];
					G2ijIm += 0.5*phia_i[0]*phia_i[0]*phib_j[0]*phib_j[1];
					G2ijIm -= 0.5*phia_i[1]*phia_i[1]*phib_j[0]*phib_j[1];
					G2ijIm += 0.5*phia_i[0]*phia_i[1]*phib_j[0]*phib_j[0];
					G2ijIm -= 0.5*phia_i[0]*phia_i[1]*phib_j[1]*phib_j[1];
					}//loop over b
				}//loop over a
			if (GijIm >= 0){
				Gfile << "(" << GijRe << "+" << GijIm << "j)" << ",";
			}
			else{
				Gfile << "(" << GijRe << GijIm << "j)" << ",";
			}
			if (G2ijIm >= 0){
				G2file << "(" << G2ijRe << "+" << G2ijIm << "j)" << ",";
			}
			else{
				G2file << "(" << G2ijRe << G2ijIm << "j)" << ",";
			}
		}//loop over x' (j)
	}//loop over x (i)
	Gfile << std::endl;
	G2file << std::endl;
	Gfile.close();
	G2file.close();
}//updated 2.1.19
*/

//NEW VERSION
void Circulation(amrex::Array4<const amrex::Real> const& Lattice, const amrex::Box& box, const amrex::GeometryData& geom, 
			long int Nt, int radius, std::string filename){
	//find the total circulation over the lattice
	//open circulation file
	std::ofstream circ_file;
	std::string r_string = std::to_string(radius);
	std::string circ_filename = "Circ_loop_"+r_string+"_"+filename.substr(8);
	circ_file.open(circ_filename, std::fstream::app);

	const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    const auto domain_xlo = geom.ProbLo();
    const auto domain_xhi = geom.ProbHi();
    const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);
    const Real x_center = 0.5 * (domain_xlo[0] + domain_xhi[0]);
    const Real y_center = 0.5 * (domain_xlo[1] + domain_xhi[1]);
    const Real dx_cell = geom.CellSize(0);
    const Real dy_cell = geom.CellSize(1);

	//Set coordinates
	int t = Nt/2;
	const Real x_left = x_center - radius;
    const Real y_bottom = y_center - radius;
    const Real x_right = x_center + radius;
    const Real y_top = y_center + radius;

	//Print loop center and radius to logfile
	circ_file << "center = (" << x_center << ","<< y_center << ") and loop radius = " << radius ;
	circ_file << std::endl;

	//Initialize sum over theta
	double theta_sum = 0.0;

	//std::cout << "center = (" << x_center << ","<< y_center << ") and loop radius = " << radius << std::endl;
	//std::cout << "{(i,j)}s = ";
	//figure out how to pick out i, j for that distance from the center
	for (int j = lo.y; j <= hi.y; ++j){
		const Real y = domain_xlo[1] + dy_cell * (j + 0.5 - domain_lo.y);
		if (y >= y_bottom && y <= y_top){
			int i_left = (x_left - domain_xlo[0]) / dx_cell - 0.5 + domain_lo.x;
			int i_right = (x_right - domain_xlo[0]) / dx_cell - 0.5 + domain_lo.x;
			//check loop indices
			//std::cout << "(" << i_left << "," << j << "), ";
			//std::cout << "(" << i_right << "," << j << "), ";
			//theta = (phi_1_Im + phi_2_Re)/(phi_1_Re - phi_2_Im);
			theta_sum += Lattice(i_left,j,t,Field(1,C::Im)) / (Lattice(i_left,j,t,Field(1,C::Re)) - Lattice(i_left,j,t,Field(2,C::Im)));
			theta_sum += Lattice(i_left,j,t,Field(2,C::Re)) / (Lattice(i_left,j,t,Field(1,C::Re)) - Lattice(i_left,j,t,Field(2,C::Im)));
			theta_sum += Lattice(i_right,j,t,Field(1,C::Im)) / (Lattice(i_right,j,t,Field(1,C::Re)) - Lattice(i_right,j,t,Field(2,C::Im)));
			theta_sum += Lattice(i_right,j,t,Field(2,C::Re)) / (Lattice(i_right,j,t,Field(1,C::Re)) - Lattice(i_right,j,t,Field(2,C::Im)));
		}
	}
	for (int i = lo.x; i <= hi.x; ++i){
		const Real x = domain_xlo[0] + dx_cell * (i + 0.5 - domain_lo.x);
		if (x >= x_left && x <= x_right){
			int j_top = (y_top - domain_xlo[1]) / dy_cell - 0.5 + domain_lo.y;
			int j_bottom = (y_bottom - domain_xlo[1]) / dy_cell - 0.5 + domain_lo.y;
			//check loop indices
			//std::cout << "(" << i << "," << j_top << "), ";
			//std::cout << "(" << i << "," << j_bottom << "), ";
			//theta = (phi_1_Im + phi_2_Re)/(phi_1_Re - phi_2_Im);
			theta_sum += Lattice(i,j_top,t,Field(1,C::Im)) / (Lattice(i,j_top,t,Field(1,C::Re)) - Lattice(i,j_top,t,Field(2,C::Im)));
			theta_sum += Lattice(i,j_top,t,Field(2,C::Re)) / (Lattice(i,j_top,t,Field(1,C::Re)) - Lattice(i,j_top,t,Field(2,C::Im)));
			theta_sum += Lattice(i,j_bottom,t,Field(1,C::Im)) / (Lattice(i,j_bottom,t,Field(1,C::Re)) - Lattice(i,j_bottom,t,Field(2,C::Im)));
			theta_sum += Lattice(i,j_bottom,t,Field(2,C::Re)) / (Lattice(i,j_bottom,t,Field(1,C::Re)) - Lattice(i,j_bottom,t,Field(2,C::Im)));
		}
	}
	//write results to circulation file
	circ_file << theta_sum/(8.*atan(1.)) << ","; //adding the circulation for one loop to the total circulation
	circ_file << std::endl;
	circ_file.close();
}