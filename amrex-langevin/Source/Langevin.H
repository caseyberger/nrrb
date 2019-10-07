#include <math.h>

#ifndef LANGEVIN_H
#define LANGEVIN_H

enum class Complex { Re = 0, Im };

AMREX_INLINE
int Field(int a, Complex complex_type) {
    // Don: assert a = 1 or 2
    if (a == 1) {
        if (complex_type == Complex::Re)
            return 0;
        else
            return 1;
    }
    else if (a == 2) {
        if (complex_type == Complex::Re)
            return 2;
        else
            return 3;
    }
    else return -1;
}

//ANTI-SYMMETRIC TENSOR
AMREX_INLINE
amrex::Real epsilon(int a, int b){
	if (a == b){return 0.;}
	else if (a == 1 and b ==2){return 1.;}
	else if (a ==2 and b == 1){return -1.;}
	else {return 0.;}//added 10.30.17
}

void Langevin_evolution(amrex::Real m, amrex::Real l, amrex::Real w, amrex::Real w_t, amrex::Real dtau, amrex::Real mu, amrex::Real eps,
        const amrex::Box& box,
        amrex::Array4<amrex::Real> const& Lattice_old,
                        amrex::Array4<amrex::Real> const& Lattice_new,
                        const amrex::GeometryData& geom);

amrex::Real K_a_Re(amrex::Real m, amrex::Real l,amrex::Real w, amrex::Real w_t, int a, amrex::Real dtau, amrex::Real mu, amrex::Array4<amrex::Real> const& Lattice, const amrex::GeometryData& geom, int i, int j, int t);

amrex::Real K_a_Im(amrex::Real m, amrex::Real l,amrex::Real w, amrex::Real w_t, int a, amrex::Real dtau, amrex::Real mu, amrex::Array4<amrex::Real> const& Lattice, const amrex::GeometryData& geom, int i, int j, int t);
#endif
