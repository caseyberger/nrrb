#ifndef LANGEVIN_H
#define LANGEVIN_H

#include <sstream>
#include <string>
#include <math.h>
#include <AMReX_REAL.H>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_Array4.H>
#include <AMReX_Utility.H>
#include <AMReX_Tuple.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>

enum class C { Re = 0, Im };

namespace Obs {
enum ObservableIndices {PhiSqRe = 0, PhiSqIm, DensRe, DensIm,
                        LzRe, LzIm, SRe, SIm, NumObservables};
};

struct NRRBParameters {
    amrex::Real m;
    amrex::Real l;
    amrex::Real w;
    amrex::Real w_t;
    amrex::Real dtau;
    amrex::Real mu;
    amrex::Real eps;
    int seed_init;
    int seed_run;

    NRRBParameters()
    {
        m = 1.0;
        l = 0.0;
        w = 0.0;
        w_t = 0.0;
        dtau = 0.0;
        mu = 0.0;
        eps = 0.0;
        seed_init = -1;
        seed_run = -1;
    }
};

AMREX_INLINE
int Field(int a, C complex_type) {
    if (a == 1) {
        if (complex_type == C::Re)
            return 0;
        else
            return 1;
    }
    else if (a == 2) {
        if (complex_type == C::Re)
            return 2;
        else
            return 3;
    }
    else {
        amrex::Error("Field() received a != 1 or 2.");
        return -1;
    }
}

//ANTI-SYMMETRIC TENSOR
AMREX_INLINE
amrex::Real epsilon(int a, int b){
	if (a == b){return 0.;}
	else if (a == 1 and b ==2){return 1.;}
	else if (a ==2 and b == 1){return -1.;}
	else {return 0.;}//added 10.30.17
}

//DIRAC DELTA
AMREX_INLINE
amrex::Real delta(int a, int b){
    if (a == b){return 1.;}
    else {return 0.;}
}

void langevin_main();

void Langevin_evolution(amrex::Real m, amrex::Real l, amrex::Real w, amrex::Real w_t, amrex::Real dtau, amrex::Real mu, amrex::Real eps,
                        const amrex::Box& box, const int Ncomp,
                        amrex::Array4<amrex::Real> const& Lattice_old,
                        amrex::Array4<amrex::Real> const& Lattice_new,
                        const amrex::GeometryData& geom);

amrex::Real K_a_Re(amrex::Real m, amrex::Real l,amrex::Real w, amrex::Real w_t, int a, amrex::Real dtau, amrex::Real mu, amrex::Array4<amrex::Real> const& Lattice, const amrex::GeometryData& geom, int i, int j, int t);

amrex::Real K_a_Im(amrex::Real m, amrex::Real l,amrex::Real w, amrex::Real w_t, int a, amrex::Real dtau, amrex::Real mu, amrex::Array4<amrex::Real> const& Lattice, const amrex::GeometryData& geom, int i, int j, int t);

void initialize_observables(const std::string);

void update_observables(const int, const amrex::MultiFab&, const amrex::GeometryData &,
                        const NRRBParameters&, const std::string);

amrex::Vector<amrex::Real> compute_observables(const amrex::Box& box, const int Ncomp,
                                               const amrex::Array4<const amrex::Real>& Lattice,
                                               const amrex::GeometryData& geom,
                                               const NRRBParameters nrrb_parm);
#endif
