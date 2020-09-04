#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <random>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <string>
#include <sstream>
#include <cstdlib>

#include "Langevin_evolution.h"

#ifndef OBSERVABLES_H
#define OBSERVABLES_H

void compute_observables(double m, double l, double w, double w_t, double dtau, double *** Lattice, int size, int dim, int Nx, int Nt, double mu, int n, double delta_t, std::string filename);

#endif