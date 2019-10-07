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
#include <vector>
#include <cstdlib>
#include "Langevin_evolution.h"
#include <omp.h>
/*
Defines the function that evolves the lattice in Langevin time, as well as the real and imaginary
drift functions, K_a.

Includes definition of an antisymmetric tensor, epsilon, and functions that modify the location 
on the lattice to one positive or negative step in the direction specified.
*/

void Langevin_evolution(double m, double l, double w, double w_t, double dtau, double *** Lattice, int size, int dim, int Nx, int Nt, double mu, double eps){
	//std::default_random_engine generator;
	double random_num_gen_time = 0.;
	double Langevin_calc_time = 0.;
	double Lattice_update_time = 0.;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::normal_distribution<double> noise1(0.,sqrt(2.));
	std::normal_distribution<double> noise2(0.,sqrt(2.));
	//evolve the lattice using the old lattice copy
	#pragma omp parallel
	#pragma omp for collapse(2)
	for(int i = 0; i<size; i++){
		for (int t = 0; t<Nt;t++){
			clock_t rn_0 = clock();
			double eta_1 = noise1(mt);
			double eta_2 = noise2(mt);
			clock_t rn_f = clock();
			random_num_gen_time += float(rn_f - rn_0)/CLOCKS_PER_SEC;
			clock_t CL_0 = clock();
			//phi1_Re
			Lattice[i][t][0] = Lattice[i][t][4]+eps*K_a_Re(m,l,w,w_t,1,dtau,Lattice,dim,Nx,Nt,i,t,mu) + sqrt(eps)*eta_1;
			//phi1_Im
			Lattice[i][t][1] = Lattice[i][t][5]+eps*K_a_Im(m,l,w,w_t,1,dtau,Lattice,dim,Nx,Nt,i,t,mu);
			//phi2_Re
			Lattice[i][t][2] = Lattice[i][t][6]+eps*K_a_Re(m,l,w,w_t,2,dtau,Lattice,dim,Nx,Nt,i,t,mu) + sqrt(eps)*eta_2;
			//phi2+Im
			Lattice[i][t][3] = Lattice[i][t][7]+eps*K_a_Im(m,l,w,w_t,2,dtau,Lattice,dim,Nx,Nt,i,t,mu);
			clock_t CL_f = clock();
			Langevin_calc_time += float(CL_f - CL_0)/CLOCKS_PER_SEC;
		}
	}
	//update the old values to the new values now that the evolution is complete
	clock_t up_0 = clock();
	double phi_1_old = 0.;
	double phi_1_new = 0.;
	double phi_2_old = 0.;
	double phi_2_new = 0.;
	for(int i = 0; i<size; i++){
		for (int t = 0; t<Nt;t++){
			phi_1_old += Lattice[i][t][4]*Lattice[i][t][4] + Lattice[i][t][5]*Lattice[i][t][5];
			phi_1_new += Lattice[i][t][0]*Lattice[i][t][0] + Lattice[i][t][1]*Lattice[i][t][1];
			phi_2_old += Lattice[i][t][6]*Lattice[i][t][6] + Lattice[i][t][7]*Lattice[i][t][7];
			phi_2_new += Lattice[i][t][2]*Lattice[i][t][2] + Lattice[i][t][3]*Lattice[i][t][3];
			Lattice[i][t][4] = Lattice[i][t][0];
			Lattice[i][t][5] = Lattice[i][t][1];
			Lattice[i][t][6] = Lattice[i][t][2];
			Lattice[i][t][7] = Lattice[i][t][3];
		}
	}
	clock_t up_f = clock();
	Lattice_update_time += float(up_f - up_0)/CLOCKS_PER_SEC;
	//std::cout << "dim = " << dim << " mu = " << mu << " Nx = " << Nx << " Nt = " << Nt << " dtau = " << dtau << " lambda " << l << std::endl;
	//std::cout << "|phi_{1}|^{2} old = " << phi_1_old << "   |phi_{1}|^{2} new = " << phi_1_new << std::endl;
	//std::cout << "|phi_{2}|^{2} old = " << phi_2_old << "   |phi_{2}|^{2} new = " << phi_2_new << std::endl;
	std::cout << "Time spent generating random numbers = " << random_num_gen_time << std::endl;
	std::cout << "Time spent calculating Langevin integral step = " << Langevin_calc_time << std::endl;
	std::cout << "Time spent updating lattice at end = " << Lattice_update_time << std::endl;
}