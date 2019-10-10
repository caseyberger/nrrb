#include "Langevin.H"
#include "AMReX_Print.H"
#include <iomanip>

using namespace amrex;

/*
Computes and prints Real and Imaginary parts of observables

don't forget to modify mu, m, w, wtr, and l by dtau if they appear in observables: 
mu = dtau*mu; m = m/dtau; w = dtau*w; wtr = dtau* wtr; l = dtau*l;
*/
//void Equal_Time_Correlators(double *** Lattice, int size, int Nx, int Nt, std::string logfilename);
void Circulation(amrex::Array4<amrex::Real> const& Lattice, const amrex::Box& box, const amrex::GeometryData& geom, 
			long int Nt, int radius, std::string filename);

void compute_observables(double m, double l, double w, double w_t, double dtau, double mu, int nL,
						double delta_t, std::string filename,
						const amrex::Box& box, const int Ncomp, 
                        amrex::Array4<amrex::Real> const& Lattice,
                        const amrex::GeometryData& geom){
	//n is number of steps in Langevin time
	std::ofstream logfile;
	logfile.open(filename, std::fstream::app);

    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

	long int volume = box.volume();
	Print() << "in Observables" << std::endl;
	Print() << "volume = " << volume << std::endl;
	long int Nt = box.length(2);
	Print() << "Nt = " << Nt << std::endl;

	//create logfile for density profiles
	std::fstream dpfile;
	std::string dp_filename = "dp_"+filename.substr(8);
	dpfile.open(dp_filename, std::fstream::app);

	//Initialize observable variables
	double densRe = 0;
	double densIm = 0;
	double phisqRe = 0.;
	double phisqIm = 0.;
	double LzRe = 0.;
	double LzIm = 0.;
    double S_Re = 0.;
    double S_Im = 0.;

	//modify parameters by dtau
	mu = dtau*mu;
	m = m/dtau;
	w_t = dtau*w_t;
	w = dtau*w;
	l = dtau*l;

	//define spacing parameters (x, y r^2, etc)
	const auto domain_xlo = geom.ProbLo();
    const auto domain_xhi = geom.ProbHi();
    const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);
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
			double dpRe = 0.;
			double dpIm = 0.;
            for (int t = lo.z; t <= hi.z; ++t){
                for (int a = 1; a <=2; a++){
                    //Field modulus squared
                    phisqRe += 0.5 * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
                    phisqRe += -0.5 * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
                    phisqIm += Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));

                    for (int b = 1; b <=2; b++){
						//Density
						dpRe += delta(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
						dpRe += -delta(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
						dpRe += -epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
						dpRe += -epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
						dpIm += delta(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Re));
						dpIm += delta(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Im));
						dpIm += epsilon(a,b) * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t-1,Field(b,C::Re));
						dpIm += -epsilon(a,b) * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t-1,Field(b,C::Im));
					}//loop over b = 1,2
				}//loop over a = 1,2

			//Angular momentum
			#if (AMREX_SPACEDIM == 3)
				//define x and y
				const Real x = x_cell - x_center;
    			const Real y = y_cell - y_center;
				for (int a = 1; a <=2; a++){//loop over a
					//real sum over a only
					LzRe += 0.5 * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(a,C::Im));
					LzRe += 0.5 * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(a,C::Re));
					LzRe += -0.5 * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(a,C::Im));
					LzRe += -0.5 * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(a,C::Re));
					LzRe += y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
					LzRe += -x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Im));
					//imaginary sum over a only
					LzIm += -0.5 * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(a,C::Re));
					LzIm += 0.5 * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(a,C::Im ));
					LzIm += 0.5 * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(a,C::Re));
					LzIm += -0.5 * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(a,C::Im));
					LzIm += -0.5 * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
					LzIm += 0.5 * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
					LzIm += 0.5 * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j,t,Field(a,C::Re));
					LzIm += -0.5 * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j,t,Field(a,C::Im));
					for (int b = 1; b <=2; b++){//loop over b
						//real sum over a and b
						LzRe += 0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Re)); 
						LzRe += -0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Im));
						LzRe += 0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im))* Lattice(i-1,j,t-1,Field(b,C::Im));
						LzRe += -0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re))* Lattice(i-1,j,t-1,Field(b,C::Re));
						//imaginary sum over a and b
						LzIm += -0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i,j-1,t-1,Field(b,C::Im));
						LzIm += -0.5 * epsilon(a,b) * x * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i,j-1,t-1,Field(b,C::Re));
						LzIm += 0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Re)) * Lattice(i-1,j,t-1,Field(b,C::Im));
						LzIm += 0.5 * epsilon(a,b) * y * Lattice(i,j,t,Field(a,C::Im)) * Lattice(i-1,j,t-1,Field(b,C::Re));
					}//loop over b for Lz
				}//loop over a for Lz
			#endif
				//ACTION
				/*
				//Stau
				for (int a =1; a<=2; a++){
					//S_tau 
					S_Re += 0.5 * phia[0]*phia[0];
					S_Re += -0.5 * phia[1]*phia[1];
					S_Re += -0.5 * exp(mu) * phia[0]*phia_mt[0];
					S_Re += 0.5 * exp(mu) * phia[1]*phia_mt[1];
					S_Im += phia[0]*phia[1];
					S_Im +=-0.5 * exp(mu) * phia[0]*phia_mt[1];
					S_Im +=-0.5 * exp(mu) * phia[1]*phia_mt[0];
					//S_del
					S_Re += 0.5*m*(AMREX_SPACEDIM -1)*phia[0]*phia[0];
					S_Re += -0.5*m*(AMREX_SPACEDIM -1)*phia[1]*phia[1];
					S_Im += m*(AMREX_SPACEDIM -1)*phia[0]*phia[1];

					for (int d = 1; d <= AMREX_SPACEDIM-1; d++) { //loop over adjacent spatial sites!
						//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
						// -(1/2m)\phi_{a,i+1}^{R} - (1/2m)\phi_{a,ii1}^{R}
				        if (d == 1)
				        {
				            if (!(i == domain_lo.x || i == domain_hi.x))
				            {
								S_Re += -0.25 * m * phia[0] * Lattice(i+1,j,t,Field(a,C::Re));
								S_Re += 0.25 * m * phia[1] * Lattice(i+1,j,t,Field(a,C::Im));
								S_Re += -0.25 * m * phia[0]* Lattice(i-1,j,t,Field(a,C::Re));
								S_Re += 0.25 * m * phia[1] * Lattice(i-1,j,t,Field(a,C::Im));;
								S_Im += -0.25 * m * phia[0] * Lattice(i+1,j,t,Field(a,C::Im));
								S_Im += -0.25 * m * phia[1] * Lattice(i+1,j,t,Field(a,C::Re));
								S_Im += -0.25 * m * phia[0]* Lattice(i-1,j,t,Field(a,C::Im));
								 S_Im += -0.25 * m * phia[1]* Lattice(i-1,j,t,Field(a,C::Re));
								for (int b=1; b<=2; b++){
									std::vector<double> phib_pi = phi_a(Lattice,i_p,t,b,new_lat);//phi_{b,i+1} --> i = x, y, z
									std::vector<double> phib_mi = phi_a(Lattice,i_n,t,b,new_lat);//phi_{b,i-1} --> i = x, y, z
									S_Re += 0.25 * m * epsilon(a,b) * phia[0] * Lattice(i+1,j,t,Field(b,C::Im));
									S_Re += 0.25 * m * epsilon(a,b) * phia[1] * Lattice(i+1,j,t,Field(b,C::Re));
									S_Re += 0.25 * m * epsilon(a,b) * phia[0] * Lattice(i-1,j,t,Field(b,C::Im));
									S_Re += 0.25 * m * epsilon(a,b) * phia[1] * Lattice(i-1,j,t,Field(b,C::Re));
									S_Im += -0.25 * m * epsilon(a,b) * phia[0] * Lattice(i+1,j,t,Field(b,C::Re));
									S_Im += 0.25 * m * epsilon(a,b) * phia[1] * Lattice(i+1,j,t,Field(b,C::Im));
									S_Im += -0.25 * m * epsilon(a,b) * phia[0] * Lattice(i-1,j,t,Field(b,C::Re));
									S_Im += 0.25 * m * epsilon(a,b) * phia[1] * Lattice(i-1,j,t,Field(b,C::Im));
								}//loop over b
				            }
				        } 
				        else if (d == 2)
				        {
				            if (!(j == domain_lo.y || j == domain_hi.y))
				            {
								S_Re += -0.25 * m * phia[0] * Lattice(i,j+1,t,Field(a,C::Re));
								S_Re += 0.25 * m * phia[1] * Lattice(i,j+1,t,Field(a,C::Im));
								S_Re += -0.25 * m * phia[0]* Lattice(i,j-1,t,Field(a,C::Re));
								S_Re += 0.25 * m * phia[1] * Lattice(i,j-1,t,Field(a,C::Im));;
								S_Im += -0.25 * m * phia[0] * Lattice(i,j+1,t,Field(a,C::Im));
								S_Im += -0.25 * m * phia[1] * Lattice(i,j+1,t,Field(a,C::Re));
								S_Im += -0.25 * m * phia[0]* Lattice(i,j-1,t,Field(a,C::Im));
								 S_Im += -0.25 * m * phia[1]* Lattice(i,j-1,t,Field(a,C::Re));
								for (int b=1; b<=2; b++){
									std::vector<double> phib_pi = phi_a(Lattice,i_p,t,b,new_lat);//phi_{b,i+1} --> i = x, y, z
									std::vector<double> phib_mi = phi_a(Lattice,i_n,t,b,new_lat);//phi_{b,i-1} --> i = x, y, z
									S_Re += 0.25 * m * epsilon(a,b) * phia[0] * Lattice(i,j+1,t,Field(b,C::Im));
									S_Re += 0.25 * m * epsilon(a,b) * phia[1] * Lattice(i,j+1,t,Field(b,C::Re));
									S_Re += 0.25 * m * epsilon(a,b) * phia[0] * Lattice(i,j-1,t,Field(b,C::Im));
									S_Re += 0.25 * m * epsilon(a,b) * phia[1] * Lattice(i,j-1,t,Field(b,C::Re));
									S_Im += -0.25 * m * epsilon(a,b) * phia[0] * Lattice(i,j+1,t,Field(b,C::Re));
									S_Im += 0.25 * m * epsilon(a,b) * phia[1] * Lattice(i,j+1,t,Field(b,C::Im));
									S_Im += -0.25 * m * epsilon(a,b) * phia[0] * Lattice(i,j-1,t,Field(b,C::Re));
									S_Im += 0.25 * m * epsilon(a,b) * phia[1] * Lattice(i,j-1,t,Field(b,C::Im));
								}//loop over b
				            }//checking that it doesn't go over the lattice in y direction
				        }//derivative in y direction
					}//derivative loop

					//S_trap
					S_Re += 0.25 * w_t * w_t * r2 * phia[0]*phia_mt[0];
					S_Re += -0.25 * w_t * w_t * r2 * phia[1]*phia_mt[1];
					S_Im += 0.25 * w_t * w_t * r2 * phia[0]*phia_mt[1];
					S_Im += 0.25 * w_t * w_t * r2 * phia[1]*phia_mt[0];

					for (int b=1; b<=2; b++){
						//S_tau
						S_Re += 0.5 * exp(mu) * epsilon(a,b) * phia[0]*phib_mt[1];
						S_Re += 0.5 * exp(mu) * epsilon(a,b) * phia[1]*phib_mt[0];
						S_Im += -0.5 * exp(mu) * epsilon(a,b) * phia[0]*phib_mt[0];
						S_Im += 0.5 * exp(mu) * epsilon(a,b) * phia[1]*phib_mt[1];
						//S_del part done in loop over dim
						//S_trap
						S_Re += -0.25 * w_t * w_t * r2 * epsilon(a,b) * phia[0]*phib_mt[1];
						S_Re += -0.25 * w_t * w_t * r2 * epsilon(a,b) * phia[1]*phib_mt[0];
						S_Im += 0.25 * w_t * w_t * r2 * epsilon(a,b) * phia[0]*phib_mt[0];
						S_Im += -0.25 * w_t * w_t * r2 * epsilon(a,b) * phia[1]*phib_mt[1];
						//S_w
						#if (AMREX_SPACEDIM == 3)
							S_Re += 0.5 * w * epsilon(a,b) * (x - y) * phia[0]*phib_mt[0];
							S_Re += -0.5 * w * epsilon(a,b) * (x - y) * phia[1]*phib_mt[1];
							S_Re += -0.5 * w * epsilon(a,b) * x * phia[0]*phib_my[0];
							S_Re += 0.5 * w * epsilon(a,b) * x * phia[1]*phib_my[1];
							S_Re += 0.5 * w * epsilon(a,b) * y * phia[0]*phib_mx[0];
							S_Re += -0.5 * w * epsilon(a,b) * y * phia[1]*phib_mx[1];
							S_Re += 0.5 * w * (x - y) * phia[0]*phib_mt[1];
							S_Re += 0.5 * w * (x - y) * phia[1]*phib_mt[0];
							S_Re += -0.5 * w * x * phia[0]*phib_my[1];
							S_Re += -0.5 * w * x * phia[1]*phib_my[0];
							//KEEP FORMATTING THE ACTION BELOW HERE
							S_Re += 0.5*w*y*(phia[0]*phib_mx[1]+phia[1]*phib_mx[0]);
							S_Im += 0.5*w*epsilon(a,b)*(x-y)*(phia[0]*phib_mt[1]+phia[1]*phib_mt[0]);
							S_Im += -0.5*w*epsilon(a,b)*x*(phia[0]*phib_my[1]+phia[1]*phib_my[0]);
							S_Im += 0.5*w*epsilon(a,b)*y*(phia[0]*phib_mx[1]+phia[1]*phib_mx[0]);
							S_Im += 0.5*w*(y-x)*(phia[0]*phib_mt[0]-phia[1]*phib_mt[1]);
							S_Im += 0.5*w*x*(phia[0]*phib_my[0]-phia[1]*phib_my[1]);
							S_Im += -0.5*w*y*(phia[0]*phib_mx[0]-phia[1]*phib_mx[1]);
						#endif
						//S_int
						S_Re += 0.25*l*(phia[1]*phia[1]*phib_mt[0]*phib_mt[0]-phia[0]*phia[0]*phib_mt[0]*phib_mt[0]);
						S_Re += 0.25*l*(phia[0]*phia[0]*phib_mt[1]*phib_mt[1]-phia[1]*phia[1]*phib_mt[1]*phib_mt[1]);
						S_Re += 0.5*l*(phia[1]*phia_mt[1]*phib[1]*phib_mt[1] - phia[1]*phia_mt[1]*phib[0]*phib_mt[0]);
						S_Re += 0.5*l*(phia[0]*phia_mt[0]*phib[0]*phib_mt[0]-phia[0]*phia_mt[0]*phib[1]*phib_mt[1]);
						S_Re += -0.5*l*(phia[0]*phia_mt[1]*phib[0]*phib_mt[1]+phia[1]*phia_mt[0]*phib[1]*phib_mt[0]);
						S_Re += -0.5*l*(phia[0]*phia_mt[1]*phib[1]*phib_mt[0]+phia[1]*phia_mt[0]*phib[0]*phib_mt[1]);
						S_Re += l*(phia[0]*phia[1]*phib_mt[0]*phib_mt[1]);
						S_Re += 0.25*l*epsilon(a,b)*(phia[0]*phia_mt[0]*phia_mt[0]*phib[1] - phia[0]*phia_mt[1]*phia_mt[1]*phib[1]);
						S_Re += 0.25*l*epsilon(a,b)*(phia[1]*phia_mt[0]*phia_mt[0]*phib[0] - phia[1]*phia_mt[1]*phia_mt[1]*phib[0]);
						S_Re += 0.25*l*epsilon(a,b)*(phia[1]*phia[1]*phia_mt[0]*phib_mt[1] - phia[0]*phia[0]*phia_mt[0]*phib_mt[1]);
						S_Re += 0.25*l*epsilon(a,b)*(phia[1]*phia[1]*phia_mt[1]*phib_mt[0] - phia[0]*phia[0]*phia_mt[1]*phib_mt[0]);
						S_Re += 0.5*l*epsilon(a,b)*(phia[0]*phia_mt[1]*phia_mt[0]*phib[0] - phia[1]*phia_mt[1]*phia_mt[0]*phib[1]);
						S_Re += 0.5*l*epsilon(a,b)*(phia[1]*phia_mt[0]*phib[1]*phib_mt[1] - phia[1]*phia_mt[0]*phib[0]*phib_mt[0]);
						S_Im += 0.5*l*(phia[0]*phia[1]*phib_mt[1]*phib_mt[1]-phia[0]*phia[1]*phib_mt[0]*phib_mt[0]);
						S_Im += 0.5*l*(phia[1]*phia[1]*phib_mt[0]*phib_mt[1]-phia[0]*phia[0]*phib_mt[0]*phib_mt[1]);
						S_Im += 0.5*l*(phia[1]*phia_mt[0]*phib[0]*phib_mt[0]-phia[0]*phia_mt[1]*phib[1]*phib_mt[1]);
						S_Im += 0.5*l*(phia[0]*phia_mt[0]*phib[0]*phib_mt[1]-phia[1]*phia_mt[1]*phib[1]*phib_mt[0]);
						S_Im += 0.5*l*(phia[0]*phia_mt[1]*phib[0]*phib_mt[0]-phia[1]*phia_mt[0]*phib[1]*phib_mt[1]);
						S_Im += 0.5*l*(phia[0]*phia_mt[0]*phib[1]*phib_mt[0]-phia[1]*phia_mt[1]*phib[0]*phib_mt[1]);
						S_Im += 0.25*l*epsilon(a,b)*(phia[0]*phia[0]*phia_mt[0]*phib_mt[0] - phia[0]*phia[0]*phia_mt[1]*phib_mt[1]);
						S_Im += 0.25*l*epsilon(a,b)*(phia[0]*phia_mt[1]*phia_mt[1]*phib[0] - phia[0]*phia_mt[0]*phia_mt[0]*phib[0]);
						S_Im += 0.25*l*epsilon(a,b)*(phia[1]*phia_mt[0]*phia_mt[0]*phib[1] - phia[1]*phia_mt[1]*phia_mt[1]*phib[1]);
						S_Im += 0.25*l*epsilon(a,b)*(phia[1]*phia[1]*phia_mt[1]*phib_mt[1]-phia[1]*phia[1]*phia_mt[0]*phib_mt[0]);
						S_Im += 0.5*l*epsilon(a,b)*(phia[1]*phia_mt[1]*phia_mt[0]*phib[0] + phia[0]*phia_mt[1]*phia_mt[0]*phib[1]);
						S_Im += -0.5*l*epsilon(a,b)*(phia[1]*phia[0]*phia_mt[0]*phib_mt[1] + phia[1]*phia[0]*phia_mt[1]*phib_mt[0]);
					}//loop over b
				}//loop over a (for the Action) */
            }//loop over t
			dpRe = 0.5*exp(mu)*dpRe;
			dpIm = 0.5*exp(mu)*dpIm;
			//add local density to avg density variables
			densRe += dpRe;
			densIm += dpIm;
			//write dp to logfile
			if (dpIm < 0){
				dpfile << "(" << dpRe/(1.0*Nt) << dpIm/(1.0*Nt) << "j)" << ",";
				}
			else{
				dpfile << "(" << dpRe/(1.0*Nt) << "+" << dpIm/(1.0*Nt) << "j)" << ",";
				}		
		}//loop over x
	}//loop over y
	dpfile << std::endl;
	dpfile.close();

	//std::cout << "Phi^{*}phi successfully computed: " << phisq[0] << " + i" << phisq[1] << std::endl;
	//Equal_Time_Correlators(Lattice, size, Nx, Nt, filename);
	
	//std::cout << "Lz successfully computed: Lz = " << Lz[0] << " + i" << Lz[1] << std::endl;
	//compute circulation
	long int Nx = box.length(0);
	Circulation(Lattice, box, geom, Nt, Nx/4, filename);
	Circulation(Lattice, box, geom, Nt, Nx/2 - 1, filename);

	//save values to file
	logfile << std::setw(6) << std::left << nL << ' ';
	logfile << std::setw(19) << std::left << phisqRe/volume << ' ';
	logfile << std::setw(19) << std::left << phisqIm/volume << ' ';
	logfile << std::setw(19) << std::left << densRe/volume << ' ';
	logfile << std::setw(19) << std::left << densIm/volume << ' ';	
	logfile << std::setw(19) << std::left << LzRe/volume << ' ';
	logfile << std::setw(19) << std::left << LzIm/volume << ' ';
	logfile << std::setw(19) << std::left << S_Re/volume << ' ';	
	logfile << std::setw(19) << std::left << S_Im/volume << ' ';	
	logfile << std::setw(11) << std::left << delta_t << std::endl;		

	logfile.close();
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
void Circulation(amrex::Array4<amrex::Real> const& Lattice, const amrex::Box& box, const amrex::GeometryData& geom, 
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

	//figure out how to pick out i, j for that distance from the center
	for (int j = lo.y; j <= hi.y; ++j){
		const Real y = domain_xlo[1] + dy_cell * (j + 0.5 - domain_lo.y);
		if (y >= y_bottom && y <= y_top){
			int i_left = (x_left - domain_xlo[0]) / dx_cell - 0.5 + domain_lo.x;
			int i_right = (x_right - domain_xlo[0]) / dx_cell - 0.5 + domain_lo.x;
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
