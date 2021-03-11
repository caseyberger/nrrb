#include "Langevin.H"
#include "Observables.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    // AMReX's Initialize will parse argc and argv to read runtime parameters
    // and initialize components of AMReX, including the random number generator.
    amrex::Initialize(argc,argv);

    langevin_main();

    amrex::Finalize();
    return 0;
}

void langevin_main()
{
    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = amrex::second();

    // run_fom is the figure of merit for the run
    // we define as total cells advanced per microsecond
    Real run_fom = 0.0;

    // AMREX_SPACEDIM: number of dimensions (space + time)
    int max_grid_size, nsteps, plot_int;

    // only plot after taking this many steps
    int plot_after_step = -1;

    // only plot after this much langevin time has elapsed
    Real plot_after_time = -1.0;

    // number of lattice sites in each dimension of the domain
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

    NRRBParameters nrrb_parm;

    int autocorrelation_step = 1;

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
        //  If plot_int < 0 then no plot files will be written
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Only plot after taking this many langevin steps
        pp.query("plot_after_step", plot_after_step);

        // Only plot after this much langevin time has elapsed
        pp.query("plot_after_time", plot_after_time);

        // Use HDF5 for writing plotfiles?
        pp.query("use_hdf5", nrrb_parm.use_hdf5);

        // Get the max grid size to use for the density profile
        // By default, set it to the same as the CL domain max grid size
        int _profile_max_grid_size = -1;
        pp.query("profile_max_grid_size", _profile_max_grid_size);
        nrrb_parm.profile_max_grid_size = (_profile_max_grid_size == -1) ? max_grid_size : _profile_max_grid_size;

        // Default nsteps to 10, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);
        pp.query("autocorrelation_step", autocorrelation_step);

        // Periodicity and boundary conditions
        pp.queryarr("is_periodic", is_periodic);
        pp.queryarr("domain_lo_bc_types", domain_lo_bc_types);
        pp.queryarr("domain_hi_bc_types", domain_hi_bc_types);

        // Problem-specific parameters
        ParmParse pp_nrrb("nrrb");
        pp_nrrb.get("m", nrrb_parm.m);
        pp_nrrb.get("l", nrrb_parm.l);
        pp_nrrb.get("w", nrrb_parm.w);
        pp_nrrb.get("w_t", nrrb_parm.w_t);
        pp_nrrb.get("dtau", nrrb_parm.dtau);
        pp_nrrb.get("mu", nrrb_parm.mu);
        pp_nrrb.get("eps", nrrb_parm.eps);

        pp_nrrb.get("circulation_radius_1", nrrb_parm.circulation_radius_1);
        pp_nrrb.get("circulation_radius_2", nrrb_parm.circulation_radius_2);

        pp_nrrb.query("seed_init", nrrb_parm.seed_init);
        pp_nrrb.query("seed_run", nrrb_parm.seed_run);

        // Problem type: 0 = default, 1 = fixed phase for circulation testing
        pp_nrrb.query("problem_type", nrrb_parm.problem_type);
    }

#ifdef TEST_SEED_RNG
    // if we set a random seed to use, then reinitialize the random number generator with it
    Print() << "Resetting random seed using seed_init = " << nrrb_parm.seed_init << std::endl;
    amrex::ResetRandomSeed(nrrb_parm.seed_init);
#endif

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

        // This defines the physical box extents in each direction.
        RealBox real_box({AMREX_D_DECL( 0.0, 0.0, 0.0)},
                         {AMREX_D_DECL( static_cast<Real>(n_cell[0]),
                                        static_cast<Real>(n_cell[1]),
                                        static_cast<Real>(n_cell[2]))});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array
    int Nghost = 1;

    // Ncomp = number of components for each array
    int Ncomp  = 4;
    const Vector<std::string> component_names = {"phi_1_Re", "phi_1_Im", "phi_2_Re", "phi_2_Im"};

    // Naux = number of auxiliary components for constructing the Langevin time derivatives of phi
    const Vector<std::string> auxiliary_names = {"K_1_Re", "K_1_Im", "K_2_Re", "K_2_Im", "eta_1", "eta_2"};

    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // We allocate two lattice multifabs; one will store the old state, the other the new.
    MultiFab lattice_old(ba, dm, Ncomp, Nghost);
    MultiFab lattice_new(ba, dm, Ncomp, Nghost);

    // Also allocate a multifab for storing auxiliary variables: the drift functions K and RNG eta
    MultiFab lattice_aux(ba, dm, AIdx::NAux, 0);
    lattice_aux.setVal(0.0);

    // Create an observables object for updating and writing observable log files
    Observables observables(geom, dm, ba, nrrb_parm, nsteps);

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

    // Initialize the valid regions in the domain interior
    Langevin_initialization(lattice_old, geom, nrrb_parm);

    // We initialized the valid regions in the interior of the domain (and not the ghost cells).
    // so now we fill the ghost cells using our boundary conditions in the Geometry object geom.
    lattice_old.FillBoundary(geom.periodicity());
    FillDomainBoundary(lattice_old, geom, lattice_bc);

    // AMReX also provides high-level MultiFab operations like Copy
    // Here we copy Ncomp components from lattice_old to lattice_new
    // starting at component index 0 in each and with Nghost ghost cells.
    MultiFab::Copy(lattice_new, lattice_old, 0, 0, Ncomp, Nghost);

    amrex::Print() << "Fields initialized" << std::endl;

    // time = starting time in the simulation
    Real Ltime = 0.0;

    // Write observables after initialization for reference
    observables.update(0, lattice_new, geom.data(), nrrb_parm);

    // To check our initialization, we also write a plotfile
    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
    {
        WritePlotfile(0, Ltime, lattice_new, component_names, lattice_aux, auxiliary_names, geom, nrrb_parm);
    }

    // init_time is the current time post-initialization
    Real init_time = amrex::second();

    for (int n = 1; n <= nsteps; ++n)
    {
#ifdef TEST_SEED_RNG
        // if we set a random seed to use, then reinitialize the random number generator with it
        Print() << "Resetting random seed using seed_run = " << nrrb_parm.seed_run << std::endl;
        amrex::ResetRandomSeed(nrrb_parm.seed_run);
#endif

        MultiFab::Copy(lattice_old, lattice_new, 0, 0, Ncomp, Nghost);

        // Advance lattice

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(lattice_old, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // This gets the index bounding box corresponding to the current MFIter object mfi.
            const Box& bx = mfi.tilebox();

            // This gets an Array4, a light wrapper for the underlying data that mfi points to.
            // The Array4 object provides accessor functions so it can be treated like a 4-D array
            // with dimensions (x, y, z, component).
            Array4<Real> const& L_old = lattice_old.array(mfi);
            Array4<Real> const& L_aux = lattice_aux.array(mfi);
            Array4<Real> const& L_new = lattice_new.array(mfi);

            // First, calculate the auxiliary quantities we need to construct the Langevin right hand side
            Langevin_auxiliary(nrrb_parm.m, nrrb_parm.l, nrrb_parm.w,
                               nrrb_parm.w_t, nrrb_parm.dtau, nrrb_parm.mu,
                               bx, Ncomp, L_old, L_aux, geom.data());

            // Then update the lattice using the Langevin right hand side terms and the stepsize eps
            Langevin_evolution(nrrb_parm.eps, bx, L_old, L_aux, L_new);
        }

        // Fill ghost cells
        lattice_new.FillBoundary(geom.periodicity());
        FillDomainBoundary(lattice_new, geom, lattice_bc);

        // Calculate observables
        if (n % autocorrelation_step == 0)
        {
            observables.update(n, lattice_new, geom.data(), nrrb_parm);
        }

        Ltime = Ltime + nrrb_parm.eps;

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data every plot_int depending on the optional
        // step or time where we turn on plotting.
        if (plot_int > 0 && n % plot_int == 0 &&
            n >= plot_after_step && Ltime >= plot_after_time)
        {
            WritePlotfile(n, Ltime, lattice_new, component_names, lattice_aux, auxiliary_names, geom, nrrb_parm);
        }

        // Update the figure of merit with the number of cells advanced
        // We can use BoxArray.numPts() since this is single-level
        run_fom += ba.numPts();
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    // over all processors
    Real stop_time = amrex::second();
    Real total_time = stop_time - strt_time;
    Real advance_time = stop_time - init_time;

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(total_time,IOProc);
    ParallelDescriptor::ReduceRealMax(advance_time,IOProc);

    // Update the figure of merit, dividing cells advanced by time to advance in microseconds.
    run_fom = run_fom / advance_time / 1.e6;

    // Tell the I/O Processor to write out the total runtime,
    // time to advance solution (w/o initialization), and figure of merit.
    Print() << "Run time = " << total_time << std::endl;
    Print() << "Run time w/o initialization = " << advance_time << std::endl;
    Print() << "  Average number of cells advanced per microsecond = " << std::fixed << std::setprecision(3) << run_fom << std::endl;
}

std::string generate_filename(std::string inputs[], int size){
    std::string str_D = inputs[0];
    std::string str_Nx = inputs[4];
    std::string str_Nt = inputs[5];
    std::string str_dt = inputs[6];
    std::string str_nL = inputs[7];
    std::string str_wtr = inputs[8];
    std::string str_eps = inputs[9];    
    std::string str_m = inputs[3];
    std::string str_w = inputs[1];
    std::string str_l = inputs[2];
    std::string str_mu = inputs[10];
    //std::string str_mu = std::to_string(inputs[10]);
    //mu_stream << std::fixed << std::setprecision(3) << inputs[10]; //truncate mu for filename
    //std::string str_mu = mu_stream.str();
    std::string filename = "logfile_D_"+str_D+"_Nx_"+str_Nx+"_Nt_"+str_Nt+"_dt_"+str_dt+"_nL_"+str_nL+"_eps_"+str_eps+"_m_"+str_m+"_wtr_"+str_wtr+"_w_"+str_w+"_l_"+str_l+"_mu_"+str_mu+".log";
    //std::cout << filename << std::endl;
    return filename;
}
