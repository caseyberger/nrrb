#include "Langevin.H"

using namespace amrex;

/*
    Initializes the domain ...
*/

void Langevin_initialization(MultiFab& lattice, const Geometry& geom, const NRRBParameters& parm)
{
    // Initialize lattice using an MFIter (MultiFab Iterator)
    // This loops over the array data corresponding to boxes owned by this MPI rank.

    // We markup this MFIter loop with OpenMP to distribute tiles to threads.
    // The AMReX random number generator interface is threadsafe.
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(lattice, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // This gets the index bounding box corresponding to the current MFIter object mfi.
        const Box& bx = mfi.tilebox();

        // This gets an Array4, a light wrapper for the underlying data that mfi points to.
        // The Array4 object provides accessor functions so it can be treated like a 4-D array
        // with dimensions (x, y, z, component).
        Array4<Real> const& Lattice = lattice.array(mfi);

        if (parm.problem_type == 1)
        {
            initialize_fixed_circulation(bx, geom.data(), Lattice);
        }
        else if (parm.problem_type == 0)
        {
            ParallelFor(bx, lattice.nComp(), [=](int i, int j, int k, int n) {
#ifdef TEST_CONSTANT_RNG
                Lattice(i, j, k, n) = 1.0;
#else
                Lattice(i, j, k, n) = 2.0 * Random() - 1.0;
#endif
            });
        }
        else
        {
            amrex::Error("Unrecognized problem_type in inputs.");
        }
    }
}


void initialize_fixed_circulation(const amrex::Box& box, const amrex::GeometryData& geom,
                                  amrex::Array4<amrex::Real> const& Lattice)
{
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

	// We are setting the phase difference theta_t_l+1 - theta_t_l
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

    for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                // Figure out the radius of our loop
                const int delta_x = std::abs(i - x_center);
                const int delta_y = std::abs(j - y_center);
                const int radius  = std::max(delta_x, delta_y);

                int loop_path_length = 0;
                bool found_location = false;
                // Figure out how many lattice sites we are from the start
                // of the loop at its lower left corner.
                //
                // First check if we are on the left edge ...
                if (i == x_center - radius) {
                    // set our path length to our distance from the point S.
                    loop_path_length = j - (y_center - radius);
                    found_location = true;
                } else {
                    // if we aren't on the left edge, add the length of the left edge
                    // to our path length
                    loop_path_length += 2 * radius;
                }

                if (!found_location && j == y_center + radius) {
                    // add our position along the top edge to our path length
                    loop_path_length += i - (x_center - radius);
                    found_location = true;
                } else if (!found_location) {
                    // if we aren't on the top edge, add the length of the top edge
                    // to our path length
                    loop_path_length += 2 * radius;
                }

                if (!found_location && i == x_center + radius) {
                    // add our position along the right edge to our path length
                    loop_path_length += (y_center + radius) - j;
                    found_location = true;
                } else if (!found_location) {
                    // if we aren't on the right edge, add the length of the right edge
                    // to our path length
                    loop_path_length += 2 * radius;
                }

                if (!found_location && j == y_center - radius) {
                    // add our position along the bottom edge to our path length
                    loop_path_length += (x_center + radius) - i;
                    found_location = true;
                } else if (!found_location) {
                    // if we aren't on the bottom edge, add the length of the bottom edge
                    // to our path length
                    loop_path_length += 2 * radius;
                }

                // As a safety check, assert that we have found the location of our cell
                // along its corresponding loop.
                AMREX_ASSERT(found_location);

                // Total number of lattice sites along this loop (not counting S twice):
                const int total_loop_length = 4 * (2 * radius);

                if (total_loop_length == 0) {
                    Lattice(i,j,k,Field(1,C::Im)) = 0.0;
                    Lattice(i,j,k,Field(2,C::Im)) = 0.0;
                    Lattice(i,j,k,Field(1,C::Re)) = 0.0;
                    Lattice(i,j,k,Field(2,C::Re)) = 0.0;
                } else {
                    // Set the field phase at this lattice site so that the phase at S is 0
                    // After total_loop_length number of phase increments from S we return back to S.
                    // The phase at S should be the same as when we started (modulo 2*pi),
                    // so, e.g. Theta(S + total_loop_length) = Theta(S) + 2 * pi
                    // so the phase increment from one loop site to the next is (2 * pi) / total_loop_length
                    const Real phase_increment = Constants::TwoPi / total_loop_length;

                    // since we are setting the phase at S to 0.0, this is just the total phase in the cell
                    const Real phase = phase_increment * loop_path_length;

                    // sets the field components so the total field is just |phi|*exp(I * theta).
                    // since tan(theta) = (phi_1_Im + phi_2_Re)/(phi_1_Re - phi_2_Im),
                    // if we set phi_1_Re = 1, phi_2_Im = phi_1_Im = 0, then phi_2_Re = tan(theta)
                    Lattice(i,j,k,Field(1,C::Im)) = 0.0;
                    Lattice(i,j,k,Field(2,C::Im)) = 0.0;
                    // setting the signs this way means atan2 can put the angle in the appropriate quadrant
                    Lattice(i,j,k,Field(1,C::Re)) = (phase > Constants::HalfPi && phase < Constants::ThreeHalvesPi) ? -1.0 : 1.0;
                    const Real abs_tan_phase = std::abs(std::tan(phase));
                    Lattice(i,j,k,Field(2,C::Re)) = (phase > Constants::Pi && phase < Constants::TwoPi) ? -abs_tan_phase : abs_tan_phase;
                }
            }
        }
    }
}