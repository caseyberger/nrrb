#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>
#include "Langevin.H"

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

    // AMREX_SPACEDIM: number of dimensions (space + time)
    int max_grid_size, nsteps, plot_int;
    Vector<int> n_cell(AMREX_SPACEDIM, 1);

    // Periodicity and Boundary Conditions
    // Defaults to External Dirichlet (0.0) in space, periodic in time
    Vector<int> is_periodic(AMREX_SPACEDIM, 0);
    Vector<int> domain_lo_bc_types(AMREX_SPACEDIM, BCType::ext_dir);
    Vector<int> domain_hi_bc_types(AMREX_SPACEDIM, BCType::ext_dir);
    // Set periodicity and BCs for time
    is_periodic[AMREX_SPACEDIM-1] = 1;
    domain_lo_bc_types[AMREX_SPACEDIM-1] = BCType::int_dir;
    domain_hi_bc_types[AMREX_SPACEDIM-1] = BCType::int_dir;

    struct NRRBParameters {
        Real m;
        Real l;
        Real w;
        Real w_t;
        Real dtau;
        Real mu;
        Real eps;
        int seed_init;
        int seed_run;
    };

    NRRBParameters nrrb_parm;
    nrrb_parm.m = 1.0;
    nrrb_parm.l = 0.0;
    nrrb_parm.w = 0.0;
    nrrb_parm.w_t = 0.0;
    nrrb_parm.dtau = 0.0;
    nrrb_parm.mu = 0.0;
    nrrb_parm.eps = 0.0;
    nrrb_parm.seed_init = -1;
    nrrb_parm.seed_run = -1;

    // inputs parameters (these have been read in by amrex::Initialize already)
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of
        //   a square (or cubic) domain.
        pp.getarr("n_cell",n_cell);

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
        pp.queryarr("domain_lo_bc_types", domain_lo_bc_types);
        pp.queryarr("domain_hi_bc_types", domain_hi_bc_types);

        ParmParse pp_nrrb("nrrb");
        pp_nrrb.get("m", nrrb_parm.m);
        pp_nrrb.get("l", nrrb_parm.l);
        pp_nrrb.get("w", nrrb_parm.w);
        pp_nrrb.get("w_t", nrrb_parm.w_t);
        pp_nrrb.get("dtau", nrrb_parm.dtau);
        pp_nrrb.get("mu", nrrb_parm.mu);
        pp_nrrb.get("eps", nrrb_parm.eps);
        pp_nrrb.query("seed_init", nrrb_parm.seed_init);
        pp_nrrb.query("seed_run", nrrb_parm.seed_run);
    }

    // if we set a random seed to use, then reinitialize the random number generator with it
    if (nrrb_parm.seed_init != -1)
    {
        Print() << "got seed_init = " << nrrb_parm.seed_init << std::endl;
        amrex::ResetRandomSeed(nrrb_parm.seed_init);
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell[0]-1, n_cell[1]-1, n_cell[2]-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

        // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL( 0.0, 0.0, 0.0)},
                         {AMREX_D_DECL( n_cell[0], n_cell[1], n_cell[2])});

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

    Vector<BCRec> lattice_bc(Ncomp);
    for (int n = 0; n < Ncomp; ++n)
    {
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
        {
            // is_periodic overrides inputs in domain_(lo/hi)_bc_type
            if (geom.isPeriodic(i))
            {
                lattice_bc[n].setLo(i, BCType::int_dir);
                lattice_bc[n].setHi(i, BCType::int_dir);
            }
            else
            {
                lattice_bc[n].setLo(i, domain_lo_bc_types[i]);
                lattice_bc[n].setHi(i, domain_hi_bc_types[i]);
            }
        }
    }

    // Initialize lattice_old using an MFIter (MultiFab Iterator)
    // This loops over the array data corresponding to boxes owned by this MPI rank.

    // We could markup this loop with OpenMP if we like.
    // The AMReX random number generator interface is threadsafe.
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(lattice_old); mfi.isValid(); ++mfi )
    {
        // This gets the index bounding box corresponding to the current MFIter object mfi.
        const Box& bx = mfi.validbox();

        // This gets an Array4, a light wrapper for the underlying data that mfi points to.
        // The Array4 object provides accessor functions so it can be treated like a 4-D array
        // with dimensions (x, y, z, component).
        Array4<Real> const& L_old = lattice_old.array(mfi);

        ParallelFor(bx, Ncomp, [=](int i, int j, int k, int n) {
            L_old(i, j, k, n) = 2.0 * Random() - 1.0;
        });
    }

    // We initialized the interior of the domain (above, we got bx by calling MFIter::validbox())
    // so now we should fill the ghost cells using our periodic boundary conditions in the Geometry object geom.
    lattice_old.FillBoundary(geom.periodicity());
    FillDomainBoundary(lattice_old, geom, lattice_bc);

    // AMReX also provides high-level MultiFab operations like Copy
    // Here we copy Ncomp components from lattice_old to lattice_new
    // starting at component index 0 in each and with Nghost ghost cells.
    MultiFab::Copy(lattice_new, lattice_old, 0, 0, Ncomp, Nghost);

    amrex::Print() << "Fields initialized" << std::endl;

    // time = starting time in the simulation
    Real Ltime = 0.0;

    // To check our initialization, we can write out plotfiles ...
    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
    {
        int n = 0;
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
        WriteSingleLevelPlotfile(pltfile, lattice_new, component_names, geom, Ltime, 0);
    }

    for (int n = 1; n <= nsteps; ++n)
    {
        // if we set a random seed to use, then reinitialize the random number generator with it
        if (nrrb_parm.seed_run != -1)
        {
            Print() << "got seed_run = " << nrrb_parm.seed_run << std::endl;
            amrex::ResetRandomSeed(nrrb_parm.seed_run);
        }

        MultiFab::Copy(lattice_old, lattice_new, 0, 0, Ncomp, Nghost);

        // make this return source term
        // do saxpy
        // Advance lattice
        for ( MFIter mfi(lattice_old); mfi.isValid(); ++mfi )
        {
            // This gets the index bounding box corresponding to the current MFIter object mfi.
            const Box& bx = mfi.validbox();

            // This gets an Array4, a light wrapper for the underlying data that mfi points to.
            // The Array4 object provides accessor functions so it can be treated like a 4-D array
            // with dimensions (x, y, z, component).
            Array4<Real> const& L_old = lattice_old.array(mfi);
            Array4<Real> const& L_new = lattice_new.array(mfi);

            Langevin_evolution(nrrb_parm.m, nrrb_parm.l, nrrb_parm.w, nrrb_parm.w_t,
                               nrrb_parm.dtau, nrrb_parm.mu, nrrb_parm.eps,
                               bx, Ncomp, L_old, L_new, geom.data());
        }

        // fill ghost cells
        lattice_new.FillBoundary(geom.periodicity());
        FillDomainBoundary(lattice_new, geom, lattice_bc);

        // Calculate observables WITH THE NEW LATTICE
        /*
        if (n % autocorrelation_step == 0) {
        for ( MFIter mfi(lattice_new); mfi.isValid(); ++mfi )
        {
            // This gets the index bounding box corresponding to the current MFIter object mfi.
            const Box& bx = mfi.validbox();
// This gets an Array4, a light wrapper for the underlying data that mfi points to.  The Array4 object provides accessor functions so it can be treated like a 4-D array with dimensions (x, y, z, component).
            Array4<Real> const& L_new = lattice_new.array(mfi);

            // Calculate_observables(m, l, w, w_t, dtau, mu, eps, lattice_old, lattice_new, geom);
        }
        }
        */

        Ltime = Ltime + nrrb_parm.eps;

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && n%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,5);
            WriteSingleLevelPlotfile(pltfile, lattice_new, component_names, geom, Ltime, 0);
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
