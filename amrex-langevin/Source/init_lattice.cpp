#include <AMReX_Box.H>
#include <AMReX_Array4.H>
#include <AMReX_Utility.H>

void init_lattice (const amrex::Box& box, const int Ncomp, amrex::Array4<amrex::Real> const& Lattice)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for (int n = 0; n < Ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    // Get uniformly distributed random number in [0,1)
                    // by calling amrex::Random() and scale it to [-1, 1)
                    Lattice(i, j, k, n) = 2.0 * amrex::Random() - 1.0;
                }
            }
        }
    }
}