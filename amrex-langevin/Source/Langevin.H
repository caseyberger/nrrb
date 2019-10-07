#include <math.h>

#ifndef LANGEVIN_H
#define LANGEVIN_H

void Langevin_evolution(amrex::Real m, amrex::Real l, amrex::Real w, amrex::Real w_t, amrex::Real dtau, amrex::Real mu, amrex::Real eps,
                        const amrex::Box& box,
                        amrex::Array4<amrex::Real> const& Lattice_old,
                        amrex::Array4<amrex::Real> const& Lattice_new,
                        const amrex::GeometryData& geom);

int field_comp(int a, bool is_real);

amrex::Real epsilon(int a, int b);

amrex::Real K_a_Re(amrex::Real m, amrex::Real l, amrex::Real w, amrex::Real w_t, int a, amrex::Real dtau, amrex::Real *** Lattice, int dim, int Nx, int Nt, int i, int t, amrex::Real mu);

amrex::Real K_a_Im(amrex::Real m, amrex::Real l, amrex::Real w, amrex::Real w_t, int a, amrex::Real dtau, amrex::Real *** Lattice, int dim, int Nx, int Nt, int i, int t, amrex::Real mu);
#endif
