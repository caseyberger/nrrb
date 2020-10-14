/*
Non-relativistic Bosons at Finite Chemical Potential and Zero Temperature with Rotating Potential
Casey Berger
Last edited: 10/02/19
Version 4.2

This loops over mu, wtrap, AND lambda now - so I can leave it running longer

The chemical potential, trap, angular momentum, and interaction are all implemented as
external gauge fields.

This version computes observables every 100th step 

*/

//standard headers
#include <iostream>
#include <vector>
#include <string>
#include <iomanip> //setw
#include <sstream> //stringstream stream
#include <fstream> 
#include <cmath>  //pow
#include <omp.h>
#include <chrono>

//custom headers
#include "test.h"
#include "lattice_init.h"
#include "lattice_save.h"
#include "Langevin_evolution.h"
#include "Observables.h"

using namespace std;

//function declaration
std::string generate_filename(std::string inputs[], int size);
std::vector<double> param_list(string param_string);

int main(int argc, char *argv[]) {
	clock_t all_t0 = clock(); //
	//set number of steps between samples:
	int acf_step = 1;
	int lattice_save_step = 1;

	//read in parameters
	string str, filename;
	vector<double> mu_vals, wtrap_vals, l_vals;
	string inputs [11] = {"dim","omega","lambda","m","Nx","Nt","dtau","nL", "w_t", "eps","mu"};//read in keywords for parameters
	
	int dim, Nx, Nt, vol, nL;
	double w, m, dtau, eps;
	if (argc != 2){ //exits if input file is not given
		cerr << "Usage: ./v1.4 input_CLB.txt"<<endl << "Exiting program" << endl;
		exit(10);
	}
	else {
		ifstream input_file(argv[1]);
		if (!input_file.is_open()){
			cerr << "input file cannot be opened";
			exit(10);
		}
		else{
			int count = 0;
			while (count < 11) {
				while (getline(input_file, str)) {
					size_t found = str.find(inputs[count]);
					size_t start;
					if (found != string::npos) {
						start = str.find_last_of(' ');
						inputs[count] = str.substr(start + 1);
						count++;
					}
				}
			}
			dim = stod(inputs[0]);
			w = stod(inputs[1]);
			string l_str = inputs[2];
			m = stod(inputs[3]);
			Nx = stoi(inputs[4]);
			vol = pow(Nx,dim);//spatial volume!
			Nt = stoi(inputs[5]);
			dtau = stod(inputs[6]);
			nL = stoi(inputs[7]);
			string w_str = inputs[8];
			eps = stod(inputs[9]);
			string mu_str = inputs[10];
			mu_vals = param_list(mu_str);
			wtrap_vals = param_list(w_str);
			l_vals = param_list(l_str);
			//cout << "parameters acquired" <<endl;
		}
		if (dim != 2 and w !=0.){
			cerr << "Rotation can only be implemented in a 2 dimensional system. Dim = " << dim << std::endl;
			exit(10);
		}
		//GENERATE LOGFILE AND START CLOCK
		//loop over lambdas
		int n_l = l_vals.size();
		for (int i_l = 0; i_l < n_l; i_l++){
			double l = l_vals[i_l];
			//loop over wtraps
			int n_w = wtrap_vals.size();
			for (int i_w = 0; i_w < n_w; i_w++){
				double w_t = wtrap_vals[i_w];
				//loop over chemical potentials 
				int n_mu = mu_vals.size();
				for (int i_mu = 0; i_mu < n_mu; i_mu++){
					//cout << "looping over mu" << endl;
					double mu = mu_vals[i_mu];
					stringstream mu_stream,l_stream,wtr_stream;
					mu_stream << std::fixed << std::setprecision(3) << mu; //truncate mu for filename
					std::string str_mu = mu_stream.str();
					l_stream << std::fixed << std::setprecision(3) << l; //truncate lambda for filename
					std::string str_l = l_stream.str();
					wtr_stream << std::fixed << std::setprecision(2) << w_t; //truncate wtr for filename
					std::string str_wtr = wtr_stream.str();
					inputs[2] = str_l;
					inputs[8] = str_wtr;
					inputs[10] = str_mu;
					filename = generate_filename(inputs, 10);
					ofstream logfile;
					logfile.open(filename,std::fstream::app);
					//cout << "Logfile created" <<endl;
					logfile << "Starting clock" << endl;
					clock_t CPU_t0 = clock();
					std::chrono::time_point<std::chrono::system_clock> t0 = std::chrono::system_clock::now();
					logfile << "Initializing fields" << endl;
					logfile << setw(7) << left << "#step ";
					logfile << setw(20) << left << "Re[phi^{*}phi] ";
					logfile << setw(20) << left << "Im[phi^{*}phi] ";
					logfile << setw(20) << left << "Re[<n>] ";
					logfile << setw(20) << left << "Im[<n>] ";
					logfile << setw(20) << left << "Re[<Lz>] ";
					logfile << setw(20) << left << "Im[<Lz>] ";
					logfile << setw(20) << left << "Re[<S>] ";
					logfile << setw(20) << left << "Im[<S>] ";
					logfile << setw(12) << left << "dt (sec) " << endl;
					logfile.close();

					//END GENERATE LOGFILE AND START CLOCK
					//Initalize the lattice - dynamically allocate the memory for the lattice
					double *** Lattice = new double**[vol];
					for(int i = 0; i < vol; i++){
						Lattice[i] = new double*[Nt];
					}
					//allocation - 4 field values, plus space to duplicate (old and new Lattice sites)
					for(int i = 0; i < vol; i++){
						for (int j = 0; j<Nt; j++){
							Lattice[i][j] = new double[8];
						}
					}
					lattice_init(Lattice, vol, Nt);//initialize phi everywhere
					//test_step_functions(Lattice, vol, dim, Nx, Nt);//print out step in every direction from every lattice site
					lattice_save(Lattice, vol, dim, Nx, Nt, mu, w, 0); //save the initial configuration

					//print time for initialization
					clock_t init_time = clock();
					clock_t init_time_2 = init_time - all_t0;
					double total_init_time = float(init_time_2)/CLOCKS_PER_SEC;
					std::cout << "time to initialize lattice and logfile = " << total_init_time << std::endl;

					clock_t Langevin_t0 = clock();

					std::chrono::time_point<std::chrono::system_clock> ti = std::chrono::system_clock::now();
					for(int k = 1; k <= nL; k++){
						std::cout <<"k = " <<  k << std::endl;
						//use Langevin equations to evolve your fields
						Langevin_evolution(m, l, w, w_t, dtau, Lattice, vol, dim, Nx, Nt, mu, eps);
						std::chrono::time_point<std::chrono::system_clock> tk = std::chrono::system_clock::now();
						std::chrono::duration<double> delta_t = tk - ti;
						if (k%acf_step == 0){
							compute_observables(m, l, w, w_t, dtau, Lattice, vol, dim, Nx, Nt, mu, k, delta_t.count(),filename); //compute the observables at each 100th step
						}
						if (k%lattice_save_step == 0){
							lattice_save(Lattice, vol, dim, Nx, Nt, mu, w, k); //save the configuration
						}//save the lattice configuration every so often
						ti = std::chrono::system_clock::now();
					}//Langevin loop
					clock_t Langevin_tf = clock();
					clock_t Langevin_time = Langevin_tf - Langevin_t0;
					//print time for Langevin evol
					double total_Langevin_time = float(Langevin_time)/CLOCKS_PER_SEC;
					std::cout << "time to to full Langevin evolution = " << total_Langevin_time << std::endl;

					delete [] Lattice; //de-allocate the memory
					clock_t CPU_tf = clock();
					std::chrono::time_point<std::chrono::system_clock> tf = std::chrono::system_clock::now();
					clock_t CPU_dt = CPU_tf - CPU_t0;
					double CPU_time = float(CPU_dt)/CLOCKS_PER_SEC;
					std::chrono::duration<double> elapsed_seconds = tf - t0;
					logfile.open(filename, std::fstream::app);
					logfile << "Total running time was " << elapsed_seconds.count() << " seconds, and total CPU time was " << CPU_time << " seconds" << endl;
					logfile.close();
				}//loop over mu
			}//loop over wtrap
		}//loop over lambda
	}//if statement checking for input file
	return 0;
}

//function definition
std::string generate_filename(std::string inputs[], int size){
	std::string str_D = inputs[0];
	std::string str_Nx = inputs[4];
	std::string str_Nt = inputs[5];
	std::string str_dt = inputs[6];
	std::string str_nL = inputs[7];
	std::string str_wtr = inputs[8];
	std::string str_eps = inputs[9];	
	std::string str_m = inputs[3];
	std::string str_w = inputs[1];
	std::string str_l = inputs[2];
	std::string str_mu = inputs[10];
	//std::string str_mu = std::to_string(inputs[10]);
	//mu_stream << std::fixed << std::setprecision(3) << inputs[10]; //truncate mu for filename
	//std::string str_mu = mu_stream.str();
	std::string filename = "logfile_D_"+str_D+"_Nx_"+str_Nx+"_Nt_"+str_Nt+"_dt_"+str_dt+"_nL_"+str_nL+"_eps_"+str_eps+"_m_"+str_m+"_wtr_"+str_wtr+"_w_"+str_w+"_l_"+str_l+"_mu_"+str_mu+".log";
	//std::cout << filename << std::endl;
	return filename;
}

std::vector<double> param_list(string param_string){
	size_t pos = 0;
	std::vector<double> param_vals;
	param_string.erase(0,1);
	while ((pos = param_string.find(',')) != std::string::npos) {
		double pval = stod(param_string.substr(0, pos));
		param_vals.push_back(pval);
		param_string.erase(0, pos + 1);
	}
	param_vals.push_back(stod(param_string.substr(0,param_string.length()-1)));
	return param_vals;
}
