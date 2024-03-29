#ifndef _KA_RE_K_H
#define _KA_RE_K_H

#include "Langevin.H"

using namespace amrex;

/*
Defines the function that evolves the lattice in Langevin time, as well as the real and imaginary
drift functions, K_a.

Includes definition of an antisymmetric tensor, epsilon, and functions that modify the location
on the lattice to one positive or negative step in the direction specified.
*/

//COMPUTING THE REAL DRIFT FUNCTION
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real K_a_Re(Real m, Real l,Real w, Real w_t, int a, Real dtau, Real mu, amrex::Array4<amrex::Real> const& Lattice, const amrex::GeometryData& geom, int i, int j, int t){
	const auto domain_xlo = geom.ProbLo();
	const auto domain_xhi = geom.ProbHi();
	const auto domain_box = geom.Domain();
	const auto domain_lo = amrex::lbound(domain_box);
	const auto domain_hi = amrex::ubound(domain_box);
	const Real x_center = 0.5 * (domain_xlo[0] + domain_xhi[0]);
	const Real y_center = 0.5 * (domain_xlo[1] + domain_xhi[1]);
	const Real dx_cell = geom.CellSize(0);
	const Real dy_cell = geom.CellSize(1);
	const Real x_cell = domain_xlo[0] + dx_cell * (i + 0.5 - domain_lo.x);
	const Real y_cell = domain_xlo[1] + dy_cell * (j + 0.5 - domain_lo.y);

	//make dimensionless radius squared
	const Real r2_cell = (std::pow(x_cell - x_center, 2) + std::pow(y_cell - y_center, 2)) / (dx_cell * dy_cell);

	Real Ka = 0.0;
	//doing -K_{a}^{R} and then returning -Ka at the end for simplicity
	//modify mu, m, w, wtr, and l by dtau
	mu = dtau*mu;
	m = dx_cell * dy_cell * m / dtau;
	w_t = dtau*w_t;
	w = dtau*w;
	l = dtau*l;

	//CHEMICAL POTENTIAL + TIME DERIVATIVE PART OF Ka (calculated for all values of mu)
	//\phi_{a,r}^{R} - e^{\mu}/2 (\phi_{a,r-t}^{R} + \phi_{a,r+t}^{R})
	//no special conditions here because we have periodic boundary conditions in time
	Ka += Lattice(i,j,t,Field(a,C::Re));
	Ka += -0.5 * std::exp(mu) * Lattice(i,j,t-1,Field(a,C::Re));
	Ka += -0.5 * std::exp(mu) * Lattice(i,j,t+1,Field(a,C::Re));

	for (int b=1; b<=2; b++){
		//\eps_{ab} e^{\mu}/2 (\phi_{b,r-t}^{I}-\phi_{b,r+t}^{I})
		Ka += epsilon(a,b) * 0.5 * std::exp(mu) * Lattice(i,j,t-1,Field(b,C::Im));
		Ka -= epsilon(a,b) * 0.5 * std::exp(mu) * Lattice(i,j,t+1,Field(b,C::Im));
	}//last checked on 10.07.19

	//SPATIAL DERIVATIVE PART OF Ka (calculated for all cases w d>0)
	//To test 0+1 dimensional case, just skip this loop
	//(d/m)\phi_{a}^{R}
	Ka += (AMREX_SPACEDIM-1)*Lattice(i,j,t,Field(a,C::Re))/m;
	for (int d = 1; d <= AMREX_SPACEDIM-1; d++) { //loop over adjacent spatial sites!
		if (d == 1)
		{
			Ka += -0.5 * Lattice(i+1,j,t,Field(a,C::Re)) / m;
			Ka += -0.5 * Lattice(i-1,j,t,Field(a,C::Re)) / m;
		}
		else if (d == 2)
		{
			Ka += -0.5 * Lattice(i,j+1,t,Field(a,C::Re)) / m;
			Ka += -0.5 * Lattice(i,j-1,t,Field(a,C::Re)) / m;
		}
	}//updated on 10.07.19 - part mixing a and b removed per notes

	//ROTATIONAL AND TRAPPING PARTS OF Ka - ONLY WORKS IN 2D AND IF OMEGA > 0.
	//TRAPPING PIECE
	if (!(w_t == 0.)){
		// 0.25*m*w_t^{2}*r2*(\phi_{a,r+t}^{R}+\phi_{a,r-t}^{R})
		Ka += 0.25*m*w_t*w_t*r2_cell * Lattice(i,j,t+1,Field(a,C::Re));
		Ka += 0.25*m*w_t*w_t*r2_cell * Lattice(i,j,t-1,Field(a,C::Re));

		for (int b=1; b<=2; b++){
			//TRAPPING PIECE
			// -0.25*m*w_t^{2} r2*\eps_{ab}(\phi_{b,r+t}^{I}+\phi_{b,r-t}^{I})
			Ka += -0.25*m*w_t*w_t*r2_cell*epsilon(a,b) * Lattice(i,j,t+1,Field(b,C::Im));
			Ka += -0.25*m*w_t*w_t*r2_cell*epsilon(a,b) * Lattice(i,j,t-1,Field(b,C::Im));
		}
	}//only sum if w_t != 0

	// ROTATIONAL PIECE
#if (AMREX_SPACEDIM == 3)
	if (!(w == 0.)){
		const Real x = (x_cell - x_center) / dx_cell; //make dimensionless
		const Real y = (y_cell - y_center) / dy_cell; //make dimensionless

		// w (x-y)(\phi_{a,r+t}^{I}+\phi_{a,r-t}^{I})
		Ka += 0.5*w*x * Lattice(i,j,t+1,Field(a, C::Im));
		Ka += -0.5*w*y * Lattice(i,j,t+1,Field(a, C::Im));
		Ka += 0.5*w*x * Lattice(i,j,t-1,Field(a, C::Im));
		Ka += -0.5*w*y * Lattice(i,j,t-1,Field(a, C::Im));

		// (w/2) y (\phi_{a,r+x+t}^{I} + \phi_{a,r-x-t}^{I}) - (w/2) x (\phi_{a,r+y+t}^{I} + \phi_{a,r-y-t}^{I})
		Ka += 0.5*w*y * Lattice(i+1,j,t+1,Field(a,C::Im));
		Ka += 0.5*w*y * Lattice(i-1,j,t-1,Field(a,C::Im));
		Ka += -0.5*w*x * Lattice(i,j+1,t+1,Field(a,C::Im));
		Ka += -0.5*w*x * Lattice(i,j-1,t-1,Field(a,C::Im));

		for (int b=1; b<=2; b++){
			//\eps_{ab} (w/2) (y (\phi_{b,r-x-t}^{R} + \phi_{b,r+x+t}^{R})- x (\phi_{b,r-y-t}^{R}+\phi_{b,r+y+t}^{R}))
			Ka += -0.5*epsilon(a,b)*w*x * Lattice(i,j+1,t+1,Field(b,C::Re));
			Ka += -0.5*epsilon(a,b)*w*x * Lattice(i,j-1,t-1,Field(b,C::Re));

			Ka += 0.5*epsilon(a,b)*w*y * Lattice(i-1,j,t-1,Field(b,C::Re));
			Ka += 0.5*epsilon(a,b)*w*y * Lattice(i+1,j,t+1,Field(b,C::Re));

			//\eps_{ab} (w/2) (x-y)(\phi_{b,r+t}^{R}+\phi_{b,r-t}^{R})
			Ka += 0.5*epsilon(a,b)*w*x * Lattice(i,j,t+1,Field(b,C::Re));
			Ka += 0.5*epsilon(a,b)*w*x * Lattice(i,j,t-1,Field(b,C::Re));

			Ka += -0.5*epsilon(a,b)*w*y * Lattice(i,j,t+1,Field(b,C::Re));
			Ka += -0.5*epsilon(a,b)*w*y * Lattice(i,j,t-1,Field(b,C::Re));
		}
	}//only sum if wz != 0
//last updated 10.07.19
#endif

	//INTERACTION PART OF Ka
	if (!(l == 0.)){
		for (int b=1; b<=2; b++){
			Ka += l*Lattice(i,j,t,Field(b,C::Re))*Lattice(i,j,t-1,Field(a,C::Re))*Lattice(i,j,t-1,Field(b,C::Re));
			Ka += -l*Lattice(i,j,t,Field(b,C::Re))*Lattice(i,j,t-1,Field(a,C::Im))*Lattice(i,j,t-1,Field(b,C::Im));

			Ka += l*Lattice(i,j,t,Field(b,C::Re))*Lattice(i,j,t+1,Field(a,C::Re))*Lattice(i,j,t+1,Field(b,C::Re));
			Ka += -l*Lattice(i,j,t,Field(b,C::Re))*Lattice(i,j,t+1,Field(a,C::Im))*Lattice(i,j,t+1,Field(b,C::Im));

			Ka += -1.0*l*Lattice(i,j,t,Field(b,C::Im))*Lattice(i,j,t-1,Field(a,C::Re))*Lattice(i,j,t-1,Field(b,C::Im));
			Ka += -l*Lattice(i,j,t,Field(b,C::Im))*Lattice(i,j,t-1,Field(a,C::Im))*Lattice(i,j,t-1,Field(b,C::Re));

			Ka += -1.0*l*Lattice(i,j,t,Field(b,C::Im))*Lattice(i,j,t+1,Field(a,C::Re))*Lattice(i,j,t+1,Field(b,C::Im));
			Ka += -l*Lattice(i,j,t,Field(b,C::Im))*Lattice(i,j,t+1,Field(a,C::Im))*Lattice(i,j,t+1,Field(b,C::Re));

			Ka += 0.5*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t-1,Field(b,C::Im))*Lattice(i,j,t-1,Field(b,C::Im));
			Ka += -0.5*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t-1,Field(b,C::Re))*Lattice(i,j,t-1,Field(b,C::Re));

			Ka += 0.5*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t+1,Field(b,C::Im))*Lattice(i,j,t+1,Field(b,C::Im));
			Ka += -0.5*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t+1,Field(b,C::Re))*Lattice(i,j,t+1,Field(b,C::Re));

			Ka += l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t-1,Field(b,C::Im))*Lattice(i,j,t-1,Field(b,C::Re));
			Ka += l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t+1,Field(b,C::Im))*Lattice(i,j,t+1,Field(b,C::Re));

			Ka += -1.*epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t-1,Field(a,C::Re))*Lattice(i,j,t-1,Field(b,C::Im));
			Ka += -epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t-1,Field(a,C::Im))*Lattice(i,j,t-1,Field(b,C::Re));

			Ka += epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t+1,Field(a,C::Re))*Lattice(i,j,t+1,Field(b,C::Im));
			Ka += epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t+1,Field(a,C::Im))*Lattice(i,j,t+1,Field(b,C::Re));

			Ka += epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t-1,Field(a,C::Im))*Lattice(i,j,t-1,Field(b,C::Im));
			Ka += -epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t-1,Field(a,C::Re))*Lattice(i,j,t-1,Field(b,C::Re));

			Ka += epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t+1,Field(a,C::Re))*Lattice(i,j,t+1,Field(b,C::Re));
			Ka += -epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t+1,Field(a,C::Im))*Lattice(i,j,t+1,Field(b,C::Im));

			Ka += 0.5*epsilon(a,b)*l*Lattice(i,j,t,Field(b,C::Im))*Lattice(i,j,t-1,Field(a,C::Re))*Lattice(i,j,t-1,Field(a,C::Re));
			Ka += -0.5*epsilon(a,b)*l*Lattice(i,j,t,Field(b,C::Im))*Lattice(i,j,t-1,Field(a,C::Im))*Lattice(i,j,t-1,Field(a,C::Im));

			Ka += 0.5*epsilon(a,b)*l*Lattice(i,j,t,Field(b,C::Im))*Lattice(i,j,t+1,Field(a,C::Im))*Lattice(i,j,t+1,Field(a,C::Im));
			Ka += -0.5*epsilon(a,b)*l*Lattice(i,j,t,Field(b,C::Im))*Lattice(i,j,t+1,Field(a,C::Re))*Lattice(i,j,t+1,Field(a,C::Re));

			Ka += epsilon(a,b)*l*Lattice(i,j,t,Field(b,C::Re))*Lattice(i,j,t-1,Field(a,C::Im))*Lattice(i,j,t-1,Field(a,C::Re));
			Ka += -epsilon(a,b)*l*Lattice(i,j,t,Field(b,C::Re))*Lattice(i,j,t+1,Field(a,C::Im))*Lattice(i,j,t+1,Field(a,C::Re));

			Ka += 0.5*epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t-1,Field(a,C::Im))*Lattice(i,j,t-1,Field(a,C::Im));
			Ka += -0.5*epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t-1,Field(a,C::Re))*Lattice(i,j,t-1,Field(a,C::Re));

			Ka += 0.5*epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t+1,Field(a,C::Re))*Lattice(i,j,t+1,Field(a,C::Re));
			Ka += -0.5*epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Im))*Lattice(i,j,t+1,Field(a,C::Im))*Lattice(i,j,t+1,Field(a,C::Im));

			Ka += epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t+1,Field(a,C::Im))*Lattice(i,j,t+1,Field(a,C::Re));
			Ka += -epsilon(a,b)*l*Lattice(i,j,t,Field(a,C::Re))*Lattice(i,j,t-1,Field(a,C::Im))*Lattice(i,j,t-1,Field(a,C::Re));
		}//last update 10.07.19 to make a gauge interaction
	}//only sum if l!=0
	return -1. * dx_cell * dy_cell * Ka;
}

#endif
