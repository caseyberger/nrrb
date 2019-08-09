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

/*
Writes the values at each point in the lattice to a file.

This needs to be updated for the dimensionality!
*/

void lattice_save(double *** Lattice, int size, int dim, int Nx, int Nt, double mu, double w, int n){
	//declare lattice size parameters
	//prep files
	std::stringstream mu_stream, w_stream;
	mu_stream << std::fixed << std::setprecision(3) << mu; //truncate mu for filename
	std::string str_mu = mu_stream.str();
	w_stream << std::fixed << std::setprecision(3) << w; //truncate mu for filename
	std::string str_w = w_stream.str();
	std::string str_N = std::to_string(Nx);
	std::string fname = "field_configs/v4_mu_" + str_mu + "_w_"+str_w+"_N_"+str_N+"_field_config.txt";//output filename
	std::ofstream fout; //output stream
	fout.open(fname, std::ios::out|std::ios::app); //open the file
	//check if files are open
    if (fout.is_open()){
    	fout << "Langevin iteration " << n << std::endl;
    	fout << std::setw(10)<<"coords"; //x,y,z,t for 3d
		fout << std::setw(12) << "phi_1^{R}" << std::setw(12) << "phi_1^{I}";
		fout << std::setw(12) << "phi_2^{R}" << std::setw(12) << "phi_2^{I}" << std::endl;
		for (int i = 0; i<size; i++){
			for(int j = 0; j<Nt; j++){
				if (dim == 1){
					int t = j;
					int x = i%Nx;
					fout << std::setw(6) << "{" << x << "," << t << "}";
					for (int k = 0; k < 4; k++){
						fout << std::setw(12) << Lattice[i][j][k];
					}
					fout << std::endl;
				}
				else if (dim == 2){
					int t = j;
					int x = i%Nx;
					int y = ((i - x)/Nx)%Nx;
					fout << std::setw(4) << "{" << x << "," << y << "," << t << "}";
					for (int k = 0; k < 4; k++){
						fout << std::setw(12) << Lattice[i][j][k];
					}
					fout << std::endl;
				}
				else if (dim ==3){
					int t = j;
					int x = i%Nx;
					int y = ((i - x)/Nx)%Nx;
					int z = ((i - x - Nx*y)/(Nx*Nx))%(Nx);
					fout << std::setw(2) << "{" << x << "," << y << "," << z << "," << t << "}";
					for (int k = 0; k < 4; k++){
						fout << std::setw(12) << Lattice[i][j][k];
					}
					fout << std::endl;
				}
			}
		}
		fout.close();
	}
	else{
		std::cerr << "Unable to open file " << std::endl;
        exit(10);
	}
}