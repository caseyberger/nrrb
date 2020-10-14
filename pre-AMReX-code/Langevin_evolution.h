#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <string>
#include <sstream>
#include <cstdlib>
#include <omp.h>

#ifndef LANGEVIN_EVOLUTION_H
#define LANGEVIN_EVOLUTION_H

void Langevin_evolution(double m, double l, double w, double w_t, double dtau, double *** Lattice, int size, int dim, int Nx, int Nt, double mu, double eps);
double epsilon(int a, int b);
std::vector<int> positive_step(int dim, int dir, int i, int t, int Nx, int Nt);
std::vector<int> negative_step(int dim, int dir, int i, int t, int Nx, int Nt);
std::vector<double> phi_a(double *** Lattice, int i, int t, int a, bool new_lat);
double K_a_Re(double m, double l, double w, double w_t, int a, double dtau, double *** Lattice, int dim, int Nx, int Nt, int i, int t, double mu);
double K_a_Im(double m, double l, double w, double w_t, int a, double dtau, double *** Lattice, int dim, int Nx, int Nt, int i, int t, double mu);
void test_step_functions(double *** Lattice, int size, int dim, int Nx, int Nt);
#endif