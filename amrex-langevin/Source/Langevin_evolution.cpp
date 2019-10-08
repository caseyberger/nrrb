#include "Langevin.H"

using namespace amrex;

/*
   Defines the function that evolves the lattice in Langevin time, as well as the real and imaginary
   drift functions, K_a.

   Includes definition of an antisymmetric tensor, epsilon, and functions that modify the location
   on the lattice to one positive or negative step in the direction specified.
*/

void Langevin_evolution(Real m, Real l, Real w, Real w_t, Real dtau, Real mu, Real eps,
                        const amrex::Box& box, const int Ncomp,
                        amrex::Array4<amrex::Real> const& Lattice_old,
                        amrex::Array4<amrex::Real> const& Lattice_new,
                        const amrex::GeometryData& geom)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < Ncomp; ++n) {
        for (int t = lo.z; t <= hi.z; ++t) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    double eta_1 = RandomNormal(0.0, sqrt(2.0));
                    double eta_2 = RandomNormal(0.0, sqrt(2.0));

                    //phi1_Re
                    Lattice_new(i,j,t,0) = Lattice_old(i,j,t,0) + eps * K_a_Re(m,l,w,w_t,1,dtau,mu,Lattice_old,geom,i,j,t) + sqrt(eps) * eta_1;

                    //phi1_Im
                    Lattice_new(i,j,t,1) = Lattice_old(i,j,t,1) + eps * K_a_Im(m,l,w,w_t,1,dtau,mu,Lattice_old,geom,i,j,t);

                    //phi2_Re
                    Lattice_new(i,j,t,2) = Lattice_old(i,j,t,2) + eps * K_a_Re(m,l,w,w_t,2,dtau,mu,Lattice_old,geom,i,j,t) + sqrt(eps) * eta_2;

                    //phi2+Im
                    Lattice_new(i,j,t,3) = Lattice_old(i,j,t,3) + eps * K_a_Im(m,l,w,w_t,2,dtau,mu,Lattice_old,geom,i,j,t);
                }
            }
        }
    }
}
