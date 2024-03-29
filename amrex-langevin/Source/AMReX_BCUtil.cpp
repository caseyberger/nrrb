#include <AMReX_BCUtil.H>
#include <AMReX_PhysBCFunct.H>

namespace amrex
{

namespace {

void dummy_cpu_fill_extdir (Box const& bx, Array4<Real> const& dest,
                            const int dcomp, const int numcomp,
                            GeometryData const& geom, const Real time,
                            const BCRec* bcr, const int bcomp,
                            const int orig_comp)
{
    // set external Dirichlet (BCType::ext_dir) BCs to 0.0
    const auto& domain = geom.Domain();
    const auto& dom_lo = amrex::lbound(domain);
    const auto& dom_hi = amrex::ubound(domain);

    ParallelFor(bx, numcomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int mcomp) {
        int n = dcomp + mcomp;

        if ((i < dom_lo.x && bcr[n].lo(0) == BCType::ext_dir) ||
            (i > dom_hi.x && bcr[n].hi(0) == BCType::ext_dir) ||

            (j < dom_lo.y && bcr[n].lo(1) == BCType::ext_dir) ||
            (j > dom_hi.y && bcr[n].hi(1) == BCType::ext_dir) ||

            (k < dom_lo.z && bcr[n].lo(2) == BCType::ext_dir) ||
            (k > dom_hi.z && bcr[n].hi(2) == BCType::ext_dir))
        {
                dest(i, j, k, n) = 0.0;
        }
    });
}

struct dummy_gpu_fill_extdir
{
    AMREX_GPU_DEVICE
    void operator() (const IntVect& iv, Array4<Real> const& dest,
                     const int dcomp, const int numcomp,
                     GeometryData const& geom, const Real time,
                     const BCRec* bcr, const int bcomp,
                     const int orig_comp) const
        {
            // set external Dirichlet (BCType::ext_dir) BCs to 0.0
            const auto& domain = geom.Domain();
            const auto& dom_lo = amrex::lbound(domain);
            const auto& dom_hi = amrex::ubound(domain);

            const auto cell = iv.dim3();
            const int i = cell.x;
            const int j = cell.y;
            const int k = cell.z;

            for (int n = dcomp; n < dcomp + numcomp; ++n) {
                if ((i < dom_lo.x && bcr[n].lo(0) == BCType::ext_dir) ||
                    (i > dom_hi.x && bcr[n].hi(0) == BCType::ext_dir) ||

                    (j < dom_lo.y && bcr[n].lo(1) == BCType::ext_dir) ||
                    (j > dom_hi.y && bcr[n].hi(1) == BCType::ext_dir) ||

                    (k < dom_lo.z && bcr[n].lo(2) == BCType::ext_dir) ||
                    (k > dom_hi.z && bcr[n].hi(2) == BCType::ext_dir))
                {
                        dest(i, j, k, n) = 0.0;
                }
            }
        }
};

}

void FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Vector<BCRec>& bc)
{
    if (geom.isAllPeriodic()) return;
    if (phi.nGrow() == 0) return;

    AMREX_ALWAYS_ASSERT(phi.ixType().cellCentered());

#if !(defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA) && defined(AMREX_GPU_PRAGMA_NO_HOST))
    if (Gpu::inLaunchRegion())
    {
#endif
        GpuBndryFuncFab<dummy_gpu_fill_extdir> gpu_bndry_func(dummy_gpu_fill_extdir{});
        PhysBCFunct<GpuBndryFuncFab<dummy_gpu_fill_extdir> > physbcf
            (geom, bc, gpu_bndry_func);
        physbcf.FillBoundary(phi, 0, phi.nComp(), phi.nGrowVect(), 0.0, 0);
#if !(defined(AMREX_USE_CUDA) && defined(AMREX_USE_GPU_PRAGMA) && defined(AMREX_GPU_PRAGMA_NO_HOST))
    }
    else
    {
        CpuBndryFuncFab cpu_bndry_func(dummy_cpu_fill_extdir);;
        PhysBCFunct<CpuBndryFuncFab> physbcf(geom, bc, cpu_bndry_func);
        physbcf.FillBoundary(phi, 0, phi.nComp(), phi.nGrowVect(), 0.0, 0);
    }
#endif
}

}
