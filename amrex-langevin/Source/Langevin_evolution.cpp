#include "Langevin.H"
#include "Ka_Re_Kernel.H"
#include "Ka_Im_Kernel.H"

using namespace amrex;

/*
   Defines the function that evolves the lattice in Langevin time, as well as the real and imaginary
   drift functions, K_a.

   Includes definition of an antisymmetric tensor, epsilon, and functions that modify the location
   on the lattice to one positive or negative step in the direction specified.
*/

void Langevin_auxiliary(Real m, Real l, Real w,
                        Real w_t, Real dtau, Real mu,
                        const amrex::Box& box, const int Ncomp,
                        amrex::Array4<amrex::Real> const& Lattice_old,
                        amrex::Array4<amrex::Real> const& Lattice_aux,
                        const amrex::GeometryData& geom)
{
    ParallelFor(box,
    [=] AMREX_GPU_DEVICE (int i, int j, int t) noexcept
    {
#ifdef TEST_CONSTANT_RNG
        Real eta_1 = 1.0;
        Real eta_2 = 1.0;
#else
#ifdef TEST_UNIFORM_RNG
        Real eta_1 = 2.0 * Random() - 1.0;
        Real eta_1 = 2.0 * Random() - 1.0;
#else
        Real eta_1 = RandomNormal(0.0, sqrt(2.0));
        Real eta_2 = RandomNormal(0.0, sqrt(2.0));
#endif
#endif

        // K_1_Re
        Lattice_aux(i,j,t,AIdx::K_1_Re) = K_a_Re(m,l,w,w_t,1,dtau,mu,Lattice_old,geom,i,j,t);

        // K_1_Im
        Lattice_aux(i,j,t,AIdx::K_1_Im) = K_a_Im(m,l,w,w_t,1,dtau,mu,Lattice_old,geom,i,j,t);

        // K_2_Re
        Lattice_aux(i,j,t,AIdx::K_2_Re) = K_a_Re(m,l,w,w_t,2,dtau,mu,Lattice_old,geom,i,j,t);

        // K_2_Im
        Lattice_aux(i,j,t,AIdx::K_2_Im) = K_a_Im(m,l,w,w_t,2,dtau,mu,Lattice_old,geom,i,j,t);

        // eta_1
        Lattice_aux(i,j,t,AIdx::eta_1) = eta_1;

        // eta_2
        Lattice_aux(i,j,t,AIdx::eta_2) = eta_2;
    });
}

void Langevin_evolution(Real eps, const amrex::Box& box,
                        amrex::Array4<amrex::Real> const& Lattice_old,
                        amrex::Array4<amrex::Real> const& Lattice_aux,
                        amrex::Array4<amrex::Real> const& Lattice_new)
{
    ParallelFor(box,
    [=] AMREX_GPU_DEVICE (int i, int j, int t) noexcept
    {
        // Advance phi_1_Re
        Lattice_new(i,j,t,Field(1,C::Re)) = (Lattice_old(i,j,t,Field(1,C::Re)) +
                                             eps * Lattice_aux(i,j,t,AIdx::K_1_Re)  +
                                             std::sqrt(eps) * Lattice_aux(i,j,t,AIdx::eta_1));

        // Advance phi_1_Im
        Lattice_new(i,j,t,Field(1,C::Im)) = (Lattice_old(i,j,t,Field(1,C::Im)) +
                                             eps * Lattice_aux(i,j,t,AIdx::K_1_Im));

        // Advance phi_2_Re
        Lattice_new(i,j,t,Field(2,C::Re)) = (Lattice_old(i,j,t,Field(2,C::Re)) +
                                             eps * Lattice_aux(i,j,t,AIdx::K_2_Re) +
                                             std::sqrt(eps) * Lattice_aux(i,j,t,AIdx::eta_2));

        // Advance phi_2_Im
        Lattice_new(i,j,t,Field(2,C::Im)) = (Lattice_old(i,j,t,Field(2,C::Im)) +
                                             eps * Lattice_aux(i,j,t,AIdx::K_2_Im));
    });
}
