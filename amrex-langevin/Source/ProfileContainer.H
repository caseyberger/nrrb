#ifndef PROFILE_CONTAINER_H_
#define PROFILE_CONTAINER_H_

#include <AMReX_Particles.H>
#include <AMReX_ParticleMesh.H>
#include <AMReX_Geometry.H>

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

    amrex::ParGDB gdb_lattice, gdb_profile;

public:

    ProfileContainer::
    ProfileContainer(const Geometry& lattice_geom,
                     const DistributionMapping& lattice_dmap,
                     const BoxArray& lattice_ba,
                     const Geometry& profile_geom,
                     const DistributionMapping& profile_dmap,
                     const BoxArray& profile_ba)
                     : ParticleContainer<PIdx::NumProfile, 0, 0, 0>(lattice_geom, lattice_dmap, lattice_ba)
    {
        gdb_lattice = amrex::ParGDB(lattice_geom, lattice_dmap, lattice_ba);
        gdb_profile = amrex::ParGDB(profile_geom, profile_dmap, profile_ba);
    }

    void AccumulateProfile(const MultiFab& state, MultiFab& profile)
    {
        // Set lattice GDB
        BasePC::m_gdb_object = gdb_lattice;

        // Clear particles
        clearParticles();

        // Initialize one particle per location in x-y for each box
        // And accumulate density profile contributions
        AccumulateParticles(state);

        // Set profile GDB
        BasePC::m_gdb_object = gdb_profile;

        // Redistribute particles using the profile geometry,
        // box array, and distribution mapping.
        Redistribute();

        // Accumulate particle averages to the profile mesh
        AccumulateMesh(profile);
    }


    void SaveMeshProfile(const MultiFab& profile)
    {
        // Save the mesh profile to an output file
    }


    void AccumulateParticles(const MultiFab& state)
    {
        const int lev = 0;

        for (MFIter mfi(state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
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

            amrex::ParallelFor(num_particles,
            [=] AMREX_GPU_DEVICE (int i) noexcept {
                ParticleType& p = pstruct[i];

                p.id()  = new_pid + i;
                p.cpu() = procID;

                p.pos(0) = i % length_x;
                p.pos(1) = i / length_x;
                p.pos(2) = 0;

                // accumulate box contents along the time dimension
                for (int k = lo.z; k <= hi.z; ++k) {
                    p.rdata(PIdx::DensRe) += 0;
                    p.rdata(PIdx::DensIm) += 0;
                }
            });
        }
    }


    void AccumulateMesh(amrex::MultiFab& profile)
    {
        const auto plo = geom.ProbLoArray();
        const auto dxi = geom.InvCellSizeArray();

        const auto dbox = geom.Domain();
        const auto dlo  = amrex::lbound(dbox);

        amrex::ParticleToMesh(*this, profile, 0,
        [=] AMREX_GPU_DEVICE (const ParticleType& p,
                              amrex::Array4<amrex::Real> const& sarr)
        {
            // Get particle offset from physical domain lo
            amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0];
            amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1];
            amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2];

            // Get the cell index where the particle is located
            int i = amrex::Math::floor(lx) + dlo.x;
            int j = amrex::Math::floor(ly) + dlo.y;
            int k = amrex::Math::floor(lz) + dlo.z;

            // Sum particle data into the mesh cell where it's located
            for (int n = 0; n < Ncomp; ++n) {
                sarr(i, j, k, n) = p.rdata(PIdx::NCoord + n);
            }
        });
    }

};

#endif
