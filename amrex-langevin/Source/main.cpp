
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

void main_main    ();

void init_lattice (const amrex::Box&, const int, amrex::Array4<amrex::Real> const&);

int main (int argc, char* argv[])
{
    // AMReX's Initialize will parse argc and argv to read runtime parameters
    // and initialize components of AMReX, including the random number generator.
    amrex::Initialize(argc,argv);

    main_main();

    amrex::Finalize();
    return 0;
}

void main_main ()
{
    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = amrex::second();

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, nsteps, plot_int;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

    // inputs parameters (these have been read in by amrex::Initialize already)
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size",max_grid_size);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be writtenq
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Default nsteps to 10, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        pp.queryarr("is_periodic", is_periodic);
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                         {AMREX_D_DECL( 1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array
    int Nghost = 1;

    // Ncomp = number of components for each array
    int Ncomp  = 4;
    const Vector<std::string> component_names = {"phi1", "phi2", "phi3", "phi4"};

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // we allocate two lattice multifabs; one will store the old state, the other the new.
    MultiFab lattice_old(ba, dm, Ncomp, Nghost);
    MultiFab lattice_new(ba, dm, Ncomp, Nghost);

    // Initialize lattice_old using an MFIter (MultiFab Iterator)
    // This loops over the array data corresponding to boxes owned by this MPI rank.
    // We could markup this loop with OpenMP if we like.
    for ( MFIter mfi(lattice_old); mfi.isValid(); ++mfi )
    {
        // This gets the index bounding box corresponding to the current MFIter object mfi.
        const Box& bx = mfi.validbox();

        // This gets an Array4, a light wrapper for the underlying data that mfi points to.
        // The Array4 object provides accessor functions so it can be treated like a 4-D array
        // with dimensions (x, y, z, component).
        Array4<Real> const& L_old = lattice_old.array(mfi);

        // Two options for setting the values in L_old:
        // (1) Pass the Box and Array4 to our custom function to manually loop through indices
        // (2) Use amrex::ParallelFor to loop through the box indices with a lambda function we provide

        // (1) pass the Box and Array4 to a function that loops over indices
        init_lattice(bx, Ncomp, L_old);

        // (2) amrex::ParallelFor() lets us do this in just 3 lines
        ParallelFor(bx, Ncomp, [=](int i, int j, int k, int n) {
            L_old(i, j, k, n) = 2.0 * Random() - 1.0;
        });
    }

    // We initialized the interior of the domain (above, we got bx by calling MFIter::validbox())
    // so now we should fill the ghost cells using our periodic boundary conditions in the Geometry object geom.
    lattice_old.FillBoundary(geom.periodicity());

    // AMReX also provides high-level MultiFab operations like Copy
    // Here we copy Ncomp components from lattice_old to lattice_new
    // starting at component index 0 in each and with Nghost ghost cells.
    MultiFab::Copy(lattice_new, lattice_old, 0, 0, Ncomp, Nghost);

    // time = starting time in the simulation
    Real time = 0.0;

    // To check our initialization, we can write out plotfiles ...
    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
        WriteSingleLevelPlotfile(pltfile, lattice_new, component_names, geom, time, 0);
    }

    // e.g. compute the time step
    const Real* dx = geom.CellSize();
    Real dt = 0.9*dx[0]*dx[0] / (2.0*AMREX_SPACEDIM);

    for (int n = 1; n <= nsteps; ++n)
    {
        MultiFab::Copy(lattice_old, lattice_new, 0, 0, Ncomp, 0);

        // e.g. new_lattice = old_lattice + dt * (something)
        // ADVANCE VARIABLES HERE, for now just copy old to new
        MultiFab::Copy(lattice_new, lattice_old, 0, 0, Ncomp, 0);

        time = time + dt;

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && n%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,5);
            WriteSingleLevelPlotfile(pltfile, lattice_new, component_names, geom, time, 0);
        }
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    // over all processors
    Real stop_time = amrex::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
