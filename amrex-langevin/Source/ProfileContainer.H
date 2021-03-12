#ifndef PROFILE_CONTAINER_H_
#define PROFILE_CONTAINER_H_

#include <AMReX_Particles.H>
#include <AMReX_ParticleMesh.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <Observables_Kernel.H>

using namespace amrex;

// Indexes into particle real struct data for the quantities to create
// a 2D profile after averaging along the time axis.
struct PIdx
{
    enum {DensRe = 0, DensIm, NumProfile};
};


class ProfileContainer
    : public amrex::ParticleContainer<PIdx::NumProfile, 0, 0, 0>
{

    // Typedef the base particle container, particle iterator and the particle type
    typedef amrex::ParticleContainer<PIdx::NumProfile, 0, 0, 0> BasePC;
    typedef amrex::ParIter<BasePC::NStructReal, 0, 0, 0> ParIter;
    typedef amrex::Particle<BasePC::NStructReal, 0> ParticleType;

    // Geometry, distribution mapping, and box array for the lattice and profile
    amrex::ParGDB gdb_lattice, gdb_profile;

    // Profile mesh data
    MultiFab profile;

    // 1/N_tau, where N_tau is the number of time points in the CL domain
    Real inv_Ntau;

    // Names of profile components
    Vector<std::string> profile_component_names;

public:

    void Setup(const Geometry& lattice_geom,
               const DistributionMapping& lattice_dmap,
               const BoxArray& lattice_ba,
               const int profile_max_grid_size)
    {
        Define(lattice_geom, lattice_dmap, lattice_ba);

        gdb_lattice = amrex::ParGDB(lattice_geom, lattice_dmap, lattice_ba);

        // create the profile geometry, box array, and distribution mapping
        const auto lattice_gdata = lattice_geom.data();
        const auto lattice_domain = lattice_gdata.Domain();
        const auto lattice_lo = amrex::lbound(lattice_domain);
        const auto lattice_hi = amrex::ubound(lattice_domain);
        inv_Ntau = 1.0_rt / static_cast<Real>(lattice_domain.length(2));

        const IntVect profile_lo(AMREX_D_DECL(lattice_lo.x, lattice_lo.y, 1));
        const IntVect profile_hi(AMREX_D_DECL(lattice_hi.x, lattice_hi.y, 1));
        const Box profile_dbox(profile_lo, profile_hi);

        BoxArray profile_ba(profile_dbox);
        profile_ba.maxSize(profile_max_grid_size);

        const auto lattice_xlo = lattice_gdata.ProbLo();
        const auto lattice_xhi = lattice_gdata.ProbHi();

        RealBox profile_real_box({AMREX_D_DECL(lattice_xlo[0], lattice_xlo[1], 0.0)},
                                 {AMREX_D_DECL(lattice_xhi[0], lattice_xhi[1], 1.0)});

        Geometry profile_geom(profile_dbox, &profile_real_box, CoordSys::cartesian, lattice_geom.isPeriodicArray().data());

        DistributionMapping profile_dmap(profile_ba);

        // define the profile MultiFab and save its geometry, box array, and distribution mapping
        profile.define(profile_ba, profile_dmap, PIdx::NumProfile, 0);

        gdb_profile = amrex::ParGDB(profile_geom, profile_dmap, profile_ba);

        // define the component names we're going to profile
        profile_component_names.push_back("DensityProfile_Re");
        profile_component_names.push_back("DensityProfile_Im");
    }


    void AccumulateProfile(const MultiFab& lattice, const NRRBParameters& nrrb_parm)
    {
        // Set lattice GDB
        BasePC::m_gdb_object = gdb_lattice;

        // Clear particles
        clearParticles();

        // Initialize one particle per location in x-y for each box
        // And accumulate density profile contributions
        CreateParticleSums(lattice, nrrb_parm);

        // Set profile GDB
        BasePC::m_gdb_object = gdb_profile;

        // Redistribute particles to the profile
        Redistribute();

        // Clear the profile before accumulating particle sums
        profile.setVal(0.0);

        // Average particle sums to the mesh
        AverageToMesh();
    }


    MultiFab const& GetMeshProfile()
    {
        // Return a const reference to the profile multifab
        return profile;
    }

    void Write(const int Langevin_step, const Real Langevin_time, const NRRBParameters& nrrb_parm)
    {
        // Write the density profile MultiFab to a plotfile at this step
        const std::string& density_pltfile = amrex::Concatenate("density_profile_plt",Langevin_step,5);
        const auto profile_geom = gdb_profile.Geom(0);
        if (nrrb_parm.use_hdf5) {
            WriteSingleLevelPlotfileHDF5(density_pltfile, profile, profile_component_names, profile_geom, Langevin_time, Langevin_step);
        } else {
            WriteSingleLevelPlotfile(density_pltfile, profile, profile_component_names, profile_geom, Langevin_time, Langevin_step);
        }
    }


    void CreateParticleSums(const MultiFab& lattice, const NRRBParameters& nrrb_parm)
    {
        const int lev = 0;

        for (MFIter mfi(lattice, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& tile_box = mfi.tilebox();

            const auto lo = amrex::lbound(tile_box);
            const auto hi = amrex::ubound(tile_box);

            const int length_x = hi.x - lo.x + 1;
            const int length_y = hi.y - lo.y + 1;

            auto& particles = GetParticles(lev);

            auto grid = mfi.index();
            auto tile = mfi.LocalTileIndex();
            auto& particle_tile = particles[std::make_pair(grid, tile)];

            const int num_particles = length_x * length_y;
            particle_tile.resize(num_particles);

            long new_pid;
            #ifdef _OPENMP
            #pragma omp critical
            #endif
            {
                new_pid = ParticleType::NextID();
                ParticleType::NextID(new_pid + num_particles);
            }

            ParticleType* pstruct = particle_tile.GetArrayOfStructs()().data();
            int procID = ParallelDescriptor::MyProc();

            auto const& L_array = lattice.array(mfi);

            amrex::ParallelFor(num_particles,
            [=] AMREX_GPU_DEVICE (int i) noexcept {
                ParticleType& p = pstruct[i];

                p.id()  = new_pid + i;
                p.cpu() = procID;

                int cell_idx_x = i % length_x + lo.x;
                int cell_idx_y = i / length_x + lo.y;

                // Store the x-y cell coordinates of the particle as its
                // x-y position coordinates. The exact position of the particle
                // in the cell doesn't matter since we're simply using the particle
                // to hold a sum for these x-y indices.
                p.pos(0) = cell_idx_x;
                p.pos(1) = cell_idx_y;
                p.pos(2) = 0.5;

                // accumulate box contents along the time dimension
                p.rdata(PIdx::DensRe) = 0;
                p.rdata(PIdx::DensIm) = 0;

                for (int k = lo.z; k <= hi.z; ++k) {
                    Real dens_Re, dens_Im;
                    local_density(cell_idx_x, cell_idx_y, k, L_array, nrrb_parm, dens_Re, dens_Im);
                    p.rdata(PIdx::DensRe) += dens_Re;
                    p.rdata(PIdx::DensIm) += dens_Im;
                }
            });
        }
    }


    void AverageToMesh()
    {
        amrex::ParticleToMesh(*this, profile, 0,
        [=] AMREX_GPU_DEVICE (const ParticleType& p,
                              amrex::Array4<amrex::Real> const& parr)
        {
            // Get the cell index where the particle is located
            int i = amrex::Math::round(p.pos(0));
            int j = amrex::Math::round(p.pos(1));
            int k = amrex::Math::round(p.pos(2));

            // Sum particle density into its corresponding profile mesh cell
            amrex::Gpu::Atomic::AddNoRet(&parr(i, j, k, PIdx::DensRe), p.rdata(PIdx::DensRe));
            amrex::Gpu::Atomic::AddNoRet(&parr(i, j, k, PIdx::DensIm), p.rdata(PIdx::DensIm));
        });

        // Since we're averaging in tau, divide the sums over tau by the number of points in tau
        profile.mult(inv_Ntau, 0, profile.nComp());
    }

};

#endif