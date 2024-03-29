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
#include <ClassyHDF.H>

enum class C { Re = 0, Im };

namespace AIdx {
    enum AuxiliaryIndexes {K_1_Re = 0, K_1_Im, K_2_Re, K_2_Im, eta_1, eta_2, NAux};
}

namespace Constants {
    constexpr amrex::Real Pi = 3.1415926535897932384e0;
    constexpr amrex::Real TwoPi = 2.0 * Pi;
    constexpr amrex::Real HalfPi = 0.5 * Pi;
    constexpr amrex::Real ThreeHalvesPi = 1.5 * Pi;
}

struct NRRBParameters {
    amrex::Real m;
    amrex::Real l;
    amrex::Real w;
    amrex::Real w_t;
    amrex::Real dtau;
    amrex::Real mu;
    amrex::Real eps;

    int circulation_radius_1;
    int circulation_radius_2;

    int seed_init;
    int seed_run;

    int problem_type;
    int profile_max_grid_size;

    NRRBParameters()
    {
        m = 1.0;
        l = 0.0;
        w = 0.0;
        w_t = 0.0;
        dtau = 0.0;
        mu = 0.0;
        eps = 0.0;

        circulation_radius_1 = 2;
        circulation_radius_2 = 9;

        seed_init = -1;
        seed_run = -1;

        // problem_type = 0 is the default setup
        // other problem types may be, e.g. specialized test cases
        problem_type = 0;

        // when we read in the NRRB parameters, we default
        // profile_max_grid_size to the CL domain max_grid_size
        // if the user doesn't set it manually.
        profile_max_grid_size = -1;
    }
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real epsilon(int a, int b){
	if (a == b){return 0.;}
	else if (a == 1 and b ==2){return 1.;}
	else if (a ==2 and b == 1){return -1.;}
	else {return 0.;}//added 10.30.17
}

//DIRAC DELTA
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real delta(int a, int b){
    if (a == b){return 1.;}
    else {return 0.;}
}

void langevin_main();

void Langevin_initialization(amrex::MultiFab& lattice, const amrex::Geometry& geom, const NRRBParameters& parm);

void initialize_fixed_circulation(const amrex::Box& box, const amrex::GeometryData& geom,
                                  amrex::Array4<amrex::Real> const& Lattice, const int problem_type);

void Langevin_auxiliary(amrex::Real m, amrex::Real l, amrex::Real w,
                        amrex::Real w_t, amrex::Real dtau, amrex::Real mu,
                        const amrex::Box& box, const int Ncomp,
                        amrex::Array4<amrex::Real> const& Lattice_old,
                        amrex::Array4<amrex::Real> const& Lattice_aux,
                        const amrex::GeometryData& geom);

void Langevin_evolution(amrex::Real eps, const amrex::Box& box,
                        amrex::Array4<amrex::Real> const& Lattice_old,
                        amrex::Array4<amrex::Real> const& Lattice_aux,
                        amrex::Array4<amrex::Real> const& Lattice_new);

std::string MakeOutputFilename(const amrex::Geometry& geom, const NRRBParameters& nrrb, const int& nsteps);

void CreateOutputFile(const std::string& filename);

void WritePlotfile(int Langevin_step, const amrex::Real Langevin_time,
                   const amrex::MultiFab& lattice, const amrex::Vector<std::string>& component_names,
                   const amrex::MultiFab& lattice_aux, const amrex::Vector<std::string>& auxiliary_names,
                   const amrex::Geometry& geom, const NRRBParameters& nrrb, const std::string& output_file);

#endif
