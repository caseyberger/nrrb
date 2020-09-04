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

using namespace std;


#ifndef LATTICE_SAVE_H
#define LATTICE_SAVE_H
void lattice_save(double *** Lattice, int size, int dim, int Nx, int Nt, double mu, double w, int n);
#endif