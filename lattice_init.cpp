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
#include "lattice_init.h"
/*
Dynamically generates a lattice with spacial, temporal, and phi elements.
At each point in space (r) and time (t), four random numbers between -1 and 1
are generated to initialize the 4 phi vectors.
*/

void lattice_init(double *** Lattice, int size, int time_size){
#ifdef TEST_SEED_RNG
	int seed = 8134;
	std::mt19937 generator(seed);
    std::uniform_real_distribution<double> distribution(-1.,1.);
#else
	std::random_device rd;
	std::mt19937 generator(rd());
	std::uniform_real_distribution<double> distribution(-1.0,1.0);
#endif
	//assign each field a random value, uniformly distributed between -1 and 1
	for (int k = 0; k<4;k++){
		for (int j = 0; j<time_size;j++){
			for (int i = 0; i<size;i++){
				int Nx = sqrt(size);
				int x = i%Nx;
				int y = ((i - x)/Nx)%Nx;
				if (k == 0){std::cout << "(" << x << "," << y << ") ";}
#ifdef TEST_CONSTANT_RNG
				double r = 1.0;
#else
				double r = distribution(generator);
#endif
				Lattice[i][j][k] = r; //time-evolved fields - for future reference
				Lattice[i][j][k+4] = r; //old fields - for future reference
			}
		}
	}
	std::cout << "Fields initialized" << std::endl;
}