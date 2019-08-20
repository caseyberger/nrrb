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
#include <chrono>

#include "Observables.h"
#include "Langevin_evolution.h"
/*
Computes and prints Real and Imaginary parts of observables

don't forget to modify mu, m, w, wtr, and l by dtau if they appear in observables: 
mu = dtau*mu; m = m/dtau; w = dtau*w; wtr = dtau* wtr; l = dtau*l;
*/

double delta(int a, int b);
std::vector<double> Density(double *** Lattice, int size, int dim, int Nx, int Nt, double dtau, double mu);
void DensityProfiles(double *** Lattice, int size, int dim, int Nx, int Nt, double dtau, double mu, std::string logfilename);
std::vector<double> Field_Modulus_Squared(double *** Lattice, int size, int Nx, int Nt);
void Equal_Time_Correlators(double *** Lattice, int size, int Nx, int Nt, std::string logfilename);
std::vector<double> Average_Lz(double *** Lattice, int size, int Nx, int Nt, int dim);
double Theta(double *** Lattice, int i, int t);
bool loop_is_on_lattice(int Nx, int Nt, int x, int y, int length);
void Circulation(double *** Lattice, int size, int Nx, int Nt, int dim, int length, std::string logfilename);
std::vector<double> Action(double *** Lattice, int size, int Nx, int Nt, int dim, double dtau, double m, double mu, double w, double w_t, double l);

void compute_observables(double m, double l, double w, double w_t, double dtau, double *** Lattice, int size, int dim, int Nx, int Nt, double mu, int n, double delta_t, std::string filename){
	//n is number of steps in Langevin time
	std::ofstream logfile;
	logfile.open(filename, std::fstream::app);
	if (n >= 0){//maybe unnecessary?
		double phi_1_Re = 0.0;
		double phi_1_Im = 0.0;
		double phi_2_Re = 0.0;
		double phi_2_Im = 0.0;
		for (int i = 0; i < size; i++){
			for (int t = 0; t < Nt; t++){
				phi_1_Re += Lattice[i][t][0];
				phi_1_Im += Lattice[i][t][1];
				phi_2_Re += Lattice[i][t][2];
				phi_2_Im += Lattice[i][t][3];
			}
		}
		//double volume = 1.0*pow(Nx+1,dim)*Nt;
		double volume = 1.0*pow(Nx,dim)*Nt;
		//std::cout << "Volume successfully computed: V = " << volume << std::endl;
		std::vector<double> density= Density(Lattice, size, dim, Nx, Nt, dtau, mu);
		double dRe = density[0];
		double dIm = density[1];
		//std::cout << "Density successfully computed: density = " << dRe << " + i" << dIm << std::endl;
		DensityProfiles(Lattice, size, dim, Nx, Nt, dtau, mu, filename);
		std::vector<double> phisq = Field_Modulus_Squared(Lattice, size, Nx, Nt);
		//std::cout << "Phi^{*}phi successfully computed: " << phisq[0] << " + i" << phisq[1] << std::endl;
		Equal_Time_Correlators(Lattice, size, Nx, Nt, filename);
		std::vector<double> Lz = Average_Lz(Lattice,size,Nx,Nt,dim);
		//std::cout << "Lz successfully computed: Lz = " << Lz[0] << " + i" << Lz[1] << std::endl;
		Circulation(Lattice,size,Nx,Nt,dim,2,filename);
		Circulation(Lattice,size,Nx,Nt,dim,4,filename);
		Circulation(Lattice,size,Nx,Nt,dim,Nx-2,filename);
		//std::cout << "Circulation successfully computed" << std::endl;
		std::vector<double> Action_S = Action(Lattice, size, Nx, Nt, dim, dtau, m, mu, w, w_t, l);
		//save values to file
		logfile << std::setw(6) << std::left << n << ' ';
		logfile << std::setw(19) << std::left << phisq[0]/volume << ' ';
		logfile << std::setw(19) << std::left << phisq[1]/volume << ' ';
		logfile << std::setw(19) << std::left << dRe/volume << ' ';
		logfile << std::setw(19) << std::left << dIm/volume << ' ';	
		logfile << std::setw(19) << std::left << Lz[0]/volume << ' ';
		logfile << std::setw(19) << std::left << Lz[1]/volume << ' ';
		logfile << std::setw(19) << std::left << Action_S[0]/volume << ' ';	
		logfile << std::setw(19) << std::left << Action_S[1]/volume << ' ';	
		logfile << std::setw(11) << std::left << delta_t << std::endl;		
	}
	logfile.close();
}

double delta(int a, int b){
	if (a == b){return 1.;}
	else {return 0.;}
}

std::vector<double> Density(double *** Lattice, int size, int dim, int Nx, int Nt, double dtau, double mu){
	std::vector<double> d;
	//modify mu by dtau
	mu = dtau*mu;
	bool new_lat = true;
	d.push_back(0.);
	d.push_back(0.);
	double densRe = 0.0;
	double densIm = 0.0;
	for (int i = 0; i< size; i++){
		for (int t=0; t<Nt; t++){
			std::vector<int> vn = negative_step(dim, 4, i, t, Nx, Nt);
			int tn = vn[1];
			for (int a = 1; a <=2; a++){
				std::vector<double> phia = phi_a(Lattice,i,t,a,new_lat);
				for (int b = 1; b <=2; b++){	
					std::vector<double> phib_mt = phi_a(Lattice,i,tn,b,new_lat);
					densRe += delta(a,b)*phia[0]*phib_mt[0]- delta(a,b)*phia[1]*phib_mt[1];
					densRe += -1.0*epsilon(a,b)*phia[1]*phib_mt[0] - epsilon(a,b)*phia[0]*phib_mt[1];
					densIm += delta(a,b)*phia[1]*phib_mt[0]+delta(a,b)*phia[0]*phib_mt[1];
					densIm += epsilon(a,b)*phia[0]*phib_mt[0]- epsilon(a,b)*phia[1]*phib_mt[1];
				}//loop over b = 1,2
			}//loop over a = 1,2
		}//loop over time
	}//loop over space
	d[0] = 0.5*exp(mu)*densRe;
	d[1] = 0.5*exp(mu)*densIm;
	//this will be normalized by the volume when it is written to the logfile
	return d;
}//checked 7.24.19

void DensityProfiles(double *** Lattice, int size, int dim, int Nx, int Nt, double dtau, double mu, std::string logfilename){
	std::fstream dpfile;
	std::string dp_filename = "dp_"+logfilename.substr(8);
	dpfile.open(dp_filename, std::fstream::app);
	//modify mu by dtau
	mu = dtau*mu;
	bool new_lat = true;
	for (int i = 0; i< size; i++){
		double dpRe = 0.;
		double dpIm = 0.;
		for (int t=0; t<Nt; t++){
			std::vector<int> vn = negative_step(dim, 4, i, t, Nx, Nt);
			int tn = vn[1];
			for (int a = 1; a <=2; a++){
				std::vector<double> phia = phi_a(Lattice,i,t,a,new_lat);
				for (int b = 1; b <=2; b++){	
					std::vector<double> phib_mt = phi_a(Lattice,i,tn,b,new_lat);
					dpRe += delta(a,b)*phia[0]*phib_mt[0]- delta(a,b)*phia[1]*phib_mt[1];
					dpRe += -1.0*epsilon(a,b)*phia[1]*phib_mt[0] - epsilon(a,b)*phia[0]*phib_mt[1];
					dpIm += delta(a,b)*phia[1]*phib_mt[0]+delta(a,b)*phia[0]*phib_mt[1];
					dpIm += epsilon(a,b)*phia[0]*phib_mt[0]- epsilon(a,b)*phia[1]*phib_mt[1];
				}//loop over b = 1,2
			}//loop over a = 1,2
		}//loop over time
		dpRe = 0.5*exp(mu)*dpRe/(1.0*Nt);
		dpIm = 0.5*exp(mu)*dpIm/(1.0*Nt);
		if (dpIm > 0){
			dpfile << "(" << dpRe << "+" << dpIm << "j)" << ",";
		}
		else{
			dpfile << "(" << dpRe << dpIm << "j)" << ",";
		}		
	}//loop over space
	dpfile << std::endl;
	dpfile.close();
}//checked 2.1.19

std::vector<double> Field_Modulus_Squared(double *** Lattice, int size, int Nx, int Nt){
	//sum the field modulus squared over every lattice point 
	std::vector<double> phisq;
	bool new_lat = true;
	phisq.push_back(0.0);
	phisq.push_back(0.0);
	double phisqRe = 0.0;
	double phisqIm = 0.0;
	for (int i = 0; i < size; i++){
		for (int t = 0; t < Nt; t++){
			for (int a = 1; a <= 2; a++){
				std::vector<double> phia = phi_a(Lattice, i, t, a,new_lat);
				phisqRe += 0.5*phia[0]*phia[0] - 0.5*phia[1]*phia[1];
				phisqIm += phia[0]*phia[1];
			}
		}
	}
	phisq[0] = phisqRe;
	phisq[1] = phisqIm;
	return phisq;
}//checked 2.1.19

void Equal_Time_Correlators(double *** Lattice, int size, int Nx, int Nt, std::string logfilename){
	//compute the equal time correlator for each point on the lattice to each other point
	bool new_lat = true;
	int t = Nt/2;
	std::ofstream Gfile;
	std::string G_filename = "G_"+logfilename.substr(8);
	Gfile.open(G_filename, std::fstream::app);
	std::ofstream G2file;
	std::string G2_filename = "G2_"+logfilename.substr(8);
	G2file.open(G2_filename, std::fstream::app);
	//std::cout << "t = " << t << std::endl;
	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){	
			double GijRe = 0.;
			double GijIm = 0.;
			double G2ijRe = 0.;
			double G2ijIm = 0.;
			for (int a = 1; a <= 2; a++){
				std::vector<double> phia_i = phi_a(Lattice, i, t, a,new_lat);
				for (int b = 1; b <=2; b++){
					std::vector<double> phib_j = phi_a(Lattice, j, t, b, new_lat);
					GijRe += 0.5*delta(a,b)*phia_i[0]*phib_j[0]-0.5*delta(a,b)*phia_i[1]*phib_j[1];
					GijRe += -0.5*epsilon(a,b)*phia_i[0]*phib_j[1]-0.5*epsilon(a,b)*phia_i[1]*phib_j[0];
					GijIm += 0.5*delta(a,b)*phia_i[0]*phib_j[1]+0.5*delta(a,b)*phia_i[1]*phib_j[0];
					GijIm += 0.5*epsilon(a,b)*phia_i[0]*phib_j[0]-0.5*epsilon(a,b)*phia_i[1]*phib_j[1];
					G2ijRe += 0.25*phia_i[0]*phia_i[0]*phib_j[0]*phib_j[0];
					G2ijRe -= 0.25*phia_i[0]*phia_i[0]*phib_j[1]*phib_j[1];
					G2ijRe -= 0.25*phia_i[1]*phia_i[1]*phib_j[0]*phib_j[0];
					G2ijRe += 0.25*phia_i[1]*phia_i[1]*phib_j[1]*phib_j[1];
					G2ijRe += phia_i[0]*phia_i[1]*phib_j[0]*phib_j[1];
					G2ijIm += 0.5*phia_i[0]*phia_i[0]*phib_j[0]*phib_j[1];
					G2ijIm -= 0.5*phia_i[1]*phia_i[1]*phib_j[0]*phib_j[1];
					G2ijIm += 0.5*phia_i[0]*phia_i[1]*phib_j[0]*phib_j[0];
					G2ijIm -= 0.5*phia_i[0]*phia_i[1]*phib_j[1]*phib_j[1];
					}//loop over b
				}//loop over a
			if (GijIm >= 0){
				Gfile << "(" << GijRe << "+" << GijIm << "j)" << ",";
			}
			else{
				Gfile << "(" << GijRe << GijIm << "j)" << ",";
			}
			if (G2ijIm >= 0){
				G2file << "(" << G2ijRe << "+" << G2ijIm << "j)" << ",";
			}
			else{
				G2file << "(" << G2ijRe << G2ijIm << "j)" << ",";
			}
		}//loop over x' (j)
	}//loop over x (i)
	Gfile << std::endl;
	G2file << std::endl;
	Gfile.close();
	G2file.close();
}//updated 2.1.19

std::vector<double> Average_Lz(double *** Lattice, int size, int Nx, int Nt, int dim){
	if (dim == 2){
		//sum the angular momentum over every lattice point
		std::vector<double> Lz;
		//std::cout << "Vector for Lz created" << std::endl;
		bool new_lat = true;
		Lz.push_back(0.);
		Lz.push_back(0.);
		//std::cout << "Lz initialized to 0" << std::endl;
		double LzRe = 0.;
		double LzIm = 0.;
		for (int t=0; t<Nt; t++){
			for (int i = 0; i< size; i++){
				int xint = i%(Nx); //x value based on (i,t) lattice numbering
				std::vector<int> vx = negative_step(dim, 1, i, t, Nx, Nt);
				int xn = vx[0];//define a step back in x direction
				//std::cout << "x, x-1 =" << xint << " " << xn%Nx << std::endl;
				int yint = ((i - xint)/Nx)%Nx; //y value based on (i,t) lattice numbering
				std::vector<int> vy = negative_step(dim, 2, i, t, Nx, Nt);
				int yn= vy[0];//define a step back in y direction
				//std::cout << "y, y-1 ="  << yint << " " << ((yn-xint)/Nx)%Nx << std::endl;
				double x = 1.*xint - 0.5*(Nx+1); //shift x so our rotation is about the center of the lattice
				double y = 1.*yint - 0.5*(Nx+1); //shift y so our rotation is about the center of the lattice
				//std::cout << "x-Nx/2, y-Nx/2 = "  << x << "," << y << std::endl;
				for (int a = 1; a <=2; a++){//loop over a
					std::vector<double> phia = phi_a(Lattice,i,t,a,new_lat);
					std::vector<double> phia_x = phi_a(Lattice,xn,t,a,new_lat);
					std::vector<double> phia_y = phi_a(Lattice,yn,t,a,new_lat);
					//real sum over a only
					LzRe += 0.5*x*phia[0]*phia_y[1]+0.5*x*phia[1]*phia_y[0];
					LzRe += -0.5*y*phia[0]*phia_x[1]-0.5*y*phia[1]*phia_x[0];
					LzRe += 1.0*y*phia[0]*phia[1] - 1.0*x*phia[0]*phia[1];
					//imaginary sum over a only
					LzIm += -0.5*x*phia[0]*phia_y[0]+0.5*x*phia[1]*phia_y[1];
					LzIm += 0.5*y*phia[0]*phia_x[0]-0.5*y*phia[1]*phia_x[1];
					LzIm += -0.5*y*phia[0]*phia[0]+0.5*y*phia[1]*phia[1];
					LzIm += 0.5*x*phia[0]*phia[0]-0.5*x*phia[1]*phia[1];
					for (int b = 1; b <=2; b++){//loop over b
						std::vector<double> phib_x = phi_a(Lattice,xn,t,b,new_lat);
						std::vector<double> phib_y = phi_a(Lattice,yn,t,b,new_lat);
						//real sum over a and b
						LzRe += 0.5*epsilon(a,b)*x*phia[0]*phib_y[0] - 0.5*epsilon(a,b)*x*phia[1]*phib_y[1];
						LzRe += 0.5*epsilon(a,b)*y*phia[1]*phib_x[1]- 0.5*epsilon(a,b)*y*phia[0]*phib_x[0];
						//imaginary sum over a and b
						LzIm += -0.5*epsilon(a,b)*x*phia[0]*phib_y[1]-0.5*epsilon(a,b)*x*phia[1]*phib_y[0];
						LzIm += 0.5*epsilon(a,b)*y*phia[0]*phib_x[1]+0.5*epsilon(a,b)*y*phia[1]*phib_x[0];
					}//loop over b
				}//loop over a
			}//loop over space
		}//loop over time
		Lz[0] = LzRe;
		Lz[1] = LzIm;
		return Lz;
	}
	else{
		std::vector<double> Lz;
		//std::cout << "Vector for Lz created" << std::endl;
		Lz.push_back(0.);
		Lz.push_back(0.);
		return Lz;
	}
}

double Theta(double *** Lattice, int i, int t){
	//compute the phase at point x,t
	double theta = 0.0;
	bool new_lat = true;
	std::vector<double> phi_1= phi_a(Lattice,i,t,1,new_lat);
	double phi_1_Re = phi_1[0];
	double phi_1_Im = phi_1[1];
	std::vector<double> phi_2 = phi_a(Lattice,i,t,2,new_lat);
	double phi_2_Re = phi_2[0];
	double phi_2_Im = phi_2[1];
	theta = (phi_1_Im + phi_2_Re)/(phi_1_Re - phi_2_Im);
	return theta;
}

bool loop_is_on_lattice(int Nx, int Nt, int x, int y, int length){
	//check that the full l x l loop is contained on the lattice
	//int x = i%(Nx);
	//int y = ((i - x)/Nx)%Nx;
	//std::cout << "x,y = " << x << ", " << y << " and Nx = " << Nx << std::endl;
	//std::cout << "length of loop is " << length << std::endl;
	if (x + length < Nx){
		if (y + length < Nx){
			//std::cout<< "loop is contained on lattice" << std::endl;
			return true;
		}
		else{
			//std::cout<< "loop is not contained on lattice" << std::endl;
			return false;
		}
	}
	else{
		//std::cout<< "loop is not contained on lattice" << std::endl;
		return false;
	}
}

/*void Circulation(double *** Lattice, int size, int Nx, int Nt, int dim, int length,std::string logfilename){
	//find the total circulation over the lattice
	if (dim == 2){
		std::ofstream circ_file;
		std::string circ_filename = "Circ_"+logfilename.substr(8);
		circ_file.open(circ_filename, std::fstream::app);
		//bool new_lat = true;
		for (int i = 0; i< size; i++){
			double loop = 0.0;
			int t = Nt/2;
			if (loop_is_on_lattice(Nx, Nt, i, length)){
				int directions[4] = {2,1,-2,-1};
				std::vector<int> x = {i,t};
				for (int d=0; d<4; d++){
					//std::cout << "Loop segment in direction " << directions[d] << std::endl;
					std::vector<int> xplusj;
					if (directions[d] > 0){
						int dir = directions[d];
						xplusj = positive_step(dim, dir, x[0], t, Nx, Nt);
					}
					else{
						int dir = abs(directions[d]);
						xplusj = negative_step(dim, dir, x[0], t, Nx, Nt);
					}
					//std::cout << "x = " << x[0] << ", x+j = " << xplusj[0] << std::endl;
					double theta = Theta(Lattice,x[0],t);
					double theta_j;
					if (xplusj[0] == 9999){
						theta_j = 0.;
					}
					else{
						theta_j = Theta(Lattice,xplusj[0],t);
					}
					loop += theta_j - theta;
					for (int j=1; j < length; j++){
						if (directions[d] > 0){
							int dir = directions[d];
							x = positive_step(dim, dir, x[0], t, Nx, Nt);//move x over 1 along j
							xplusj = positive_step(dim, dir, x[0], t, Nx, Nt);//move x+j accordingly
						}
						else{
							int dir = abs(directions[d]);
							x = negative_step(dim, dir, x[0], t, Nx, Nt);//move x over 1 along j
							xplusj = negative_step(dim, dir, x[0], t, Nx, Nt);//move x+j accordingly
						}
						//std::cout << "x = " << x[0] << ", x+j = " << xplusj[0] << std::endl;
						theta = Theta(Lattice,x[0],t);
						if (xplusj[0] == 9999){
							theta_j = 0.;
							//std::cout << theta_j << std::endl;
						}
						else{
							theta_j = Theta(Lattice,xplusj[0],t);
						}
						loop += theta_j - theta;
						//std::cout << "theta_{x+j} = " << theta_j << ' ' << "theta_{x} = " << theta << std::endl;
					}
				}//loop over direction (1,-2,-1,2)
				//std::cout << "circulation around one loop = " << loop/(8.*atan(1.)) << std::endl;
				circ_file <<loop/(8.*atan(1.)) << ","; //adding the circulation for one loop to the total circulation
			}//checking loop is contained within lattice
		}//loop over space (i) 
		circ_file << std::endl;
		circ_file.close();
	}//do nothing if we are not in two dimensions
}//OLD VERSION*/


//NEW VERSION
void Circulation(double *** Lattice, int size, int Nx, int Nt, int dim, int length,std::string logfilename){
	//find the total circulation over the lattice
	if (dim == 2){
		std::ofstream circ_file;
		string length_str = to_string(length);
		std::string circ_filename = "Circ_loop_"+length_str+"_"+logfilename.substr(8);
		circ_file.open(circ_filename, std::fstream::app);
		int center = Nx/2;
		int x = center - length/2;
		int y = center - length/2;
		int t = Nt/2;
		circ_file << "center = (" << center << ","<< center << ") and start = (" << x << "," << y << ")";
		circ_file << std::endl;
		double loop = 0.0;
		if (loop_is_on_lattice(Nx, Nt, x, y, length)){
			int directions[4] = {2,1,-2,-1};	
			int i = x + Nx*y;		
			std::vector<int> xvec = {i,t};
			for (int d=0; d<4; d++){
				//std::cout << "Loop segment in direction " << directions[d] << std::endl;
				std::vector<int> xplusj;
				if (directions[d] > 0){
					int dir = directions[d];
					xplusj = positive_step(dim, dir, xvec[0], t, Nx, Nt);
				}
				else{
					int dir = abs(directions[d]);
					xplusj = negative_step(dim, dir, xvec[0], t, Nx, Nt);
				}
				//std::cout << "x = " << x[0] << ", x+j = " << xplusj[0] << std::endl;
				double theta = Theta(Lattice,xvec[0],t);
				double theta_j;
				if (xplusj[0] == 9999){
					theta_j = 0.;
				}
				else{
					theta_j = Theta(Lattice,xplusj[0],t);
				}
				loop += theta_j - theta;
				for (int j=1; j < length; j++){
					if (directions[d] > 0){
						int dir = directions[d];
						xvec = positive_step(dim, dir, xvec[0], t, Nx, Nt);//move x over 1 along j
						xplusj = positive_step(dim, dir, xvec[0], t, Nx, Nt);//move x+j accordingly
					}
					else{
						int dir = abs(directions[d]);
						xvec = negative_step(dim, dir, xvec[0], t, Nx, Nt);//move x over 1 along j
						xplusj = negative_step(dim, dir, xvec[0], t, Nx, Nt);//move x+j accordingly
					}
					//std::cout << "x = " << x[0] << ", x+j = " << xplusj[0] << std::endl;
					theta = Theta(Lattice,xvec[0],t);
					if (xplusj[0] == 9999){
						theta_j = 0.;
						//std::cout << theta_j << std::endl;
					}
					else{
						theta_j = Theta(Lattice,xplusj[0],t);
					}
					loop += theta_j - theta;
					//std::cout << "theta_{x+j} = " << theta_j << ' ' << "theta_{x} = " << theta << std::endl;
				}
			}//loop over direction (1,-2,-1,2)
			//std::cout << "circulation around one loop = " << loop/(8.*atan(1.)) << std::endl;
			circ_file << loop/(8.*atan(1.)) << ","; //adding the circulation for one loop to the total circulation
		}//checking loop is contained within lattice
		circ_file << std::endl;
		circ_file.close();
	}//do nothing if we are not in two dimensions
}

std::vector<double> Action(double *** Lattice, int size, int Nx, int Nt, int dim, double dtau, double m, double mu, double w, double w_t, double l){
	//same lattice time, separation of length r - how to do this in multiple directions? What about diagonals??
	std::vector<double> S_vec;
	//modify mu, m, w, and l by dtau
	mu = dtau*mu;
	m = m/dtau;
	w = dtau*w;
	l = dtau*l;
	double wt2 = 1.0*w_t*w_t;
	bool new_lat = true;
	S_vec.push_back(0.);
	S_vec.push_back(0.);
	double S_Re = 0.;
	double S_Im = 0.;
	for (int t=0; t<Nt; t++){
		for (int i = 0; i< size; i++){
			std::vector<int> vn = negative_step(dim, 4, i, t, Nx, Nt);
			std::vector<int> vp = positive_step(dim, 4, i, t, Nx, Nt);
			std::vector<int> vpx = positive_step(dim, 1, i, t, Nx, Nt);
			std::vector<int> vnx = negative_step(dim, 1, i, t, Nx, Nt);
			std::vector<int> vpy = positive_step(dim, 2, i, t, Nx, Nt);
			std::vector<int> vny = negative_step(dim, 2, i, t, Nx, Nt);
			int tn = vn[1];
			int tp = vp[1];
			//Stau
			for (int a =1; a<=2; a++){
				std::vector<double> phia = phi_a(Lattice,i,t,a,new_lat);
				std::vector<double> phia_mt = phi_a(Lattice,i,tn,a,new_lat);
				std::vector<double> phia_px = phi_a(Lattice,vpx[0],tp,a,new_lat);//phi_{a,r+x+t}
				std::vector<double> phia_mx = phi_a(Lattice,vnx[0],tn,a,new_lat);//phi_{a,r-x-t}
				std::vector<double> phia_py = phi_a(Lattice,vpy[0],tp,a,new_lat);//phi_{a,r+y+t}
				std::vector<double> phia_my = phi_a(Lattice,vny[0],tn,a,new_lat);//phi_{a,r-y-t}
				//S_tau 
				S_Re += 0.5*(phia[0]*phia[0] - phia[1]*phia[1]);
				S_Re += -0.5*exp(mu)*(phia[0]*phia_mt[0] - phia[1]*phia_mt[1]);
				S_Im += phia[0]*phia[1];
				S_Im +=-0.5*exp(mu)*(phia[0]*phia_mt[1] + phia[1]*phia_mt[0]);
				//S_del
				S_Re += 0.5*m*dim*(phia[0]*phia[0] - phia[1]*phia[1]);
				S_Im += m*dim*(phia[0]*phia[1]);
				for (int d = 1; d <= dim; d++){//loop over adjacent spatial sites!
					std::vector<int> vp = positive_step(dim, d, i, t, Nx, Nt);
					std::vector<int> vn = negative_step(dim, d, i, t, Nx, Nt);
					int i_p = vp[0];
					int i_n = vn[0];
					std::vector<double> phia_pi = phi_a(Lattice,i_p,t,a,new_lat);//phi_{a,i+1} --> i = x, y, z
					std::vector<double> phia_mi = phi_a(Lattice,i_n,t,a,new_lat);//phi_{a,i-1} --> i = x, y, z
					//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
					if (!(i_p == 9999. || i_n == 9999.)){
						S_Re += -0.25*m*(phia[0]*phia_pi[0]-phia[1]*phia_pi[1]);
						S_Re += -0.25*m*(phia[0]*phia_mi[0]-phia[1]*phia_mi[1]);
						S_Im += -0.25*m*(phia[0]*phia_pi[1] + phia[1]*phia_pi[0]);
						S_Im += -0.25*m*(phia[0]*phia_mi[1] + phia[1]*phia_mi[0]);
						for (int b=1; b<=2; b++){
							std::vector<double> phib_pi = phi_a(Lattice,i_p,t,b,new_lat);//phi_{b,i+1} --> i = x, y, z
							std::vector<double> phib_mi = phi_a(Lattice,i_n,t,b,new_lat);//phi_{b,i-1} --> i = x, y, z
							S_Re += 0.25*m*epsilon(a,b)*(phia[0]*phib_pi[1]+phia[1]*phib_pi[0]);
							S_Re += 0.25*m*epsilon(a,b)*(phia[0]*phib_mi[1]+phia[1]*phib_mi[0]);
							S_Im += -0.25*m*epsilon(a,b)*(phia[0]*phib_pi[0]-phia[1]*phib_pi[1]);
							S_Im += -0.25*m*epsilon(a,b)*(phia[0]*phib_mi[0]-phia[1]*phib_mi[1]);

						}//loop over b
					}//checking that the derivative doesn't go over the edge
				}//loop over adjacent spatial sites, end S_del part
				//S_trap
				if (dim ==2){
					int xint = i%(Nx);
					int yint = ((i - xint)/Nx)%Nx;
					//shift x, y by (Nx-1)/2
					double x = 1.0*xint - 0.5*(Nx-1);
					double y = 1.0*yint - 0.5*(Nx-1);
					double r2 = x*x+y*y;
					S_Re += 0.25*wt2*r2*(phia[0]*phia_mt[0] - phia[1]*phia_mt[1]);
					S_Im += 0.25*wt2*r2*(phia[0]*phia_mt[1] + phia[1]*phia_mt[0]);
				}//end S_trap
				for (int b=1; b<=2; b++){
					//S_tau
					std::vector<double> phib = phi_a(Lattice,i,t,b,new_lat);
					std::vector<double> phib_mt = phi_a(Lattice,i,tn,b,new_lat);
					std::vector<double> phib_my = phi_a(Lattice,vny[0],tn,b,new_lat);
					std::vector<double> phib_mx = phi_a(Lattice,vnx[0],tn,b,new_lat);
					S_Re += 0.5*exp(mu)*epsilon(a,b)*(phia[0]*phib_mt[1]+phia[1]*phib_mt[0]);
					S_Im += -0.5*exp(mu)*epsilon(a,b)*(phia[0]*phib_mt[0] - phia[1]*phib_mt[1]);
					//S_del part done in loop over dim
					//S_trap and S_w
					if (dim ==2){
						int xint = i%(Nx);
						int yint = ((i - xint)/Nx)%Nx;
						//shift x, y by (Nx-1)/2
						double x = 1.0*xint - 0.5*(Nx-1);
						double y = 1.0*yint - 0.5*(Nx-1);
						double r2 = x*x+y*y;
						//S_trap
						S_Re += -0.25*wt2*r2*epsilon(a,b)*(phia[0]*phib_mt[1] + phia[1]*phib_mt[0]);
						S_Im += 0.25*wt2*r2*epsilon(a,b)*(phia[0]*phib_mt[0] - phia[1]*phib_mt[1]);
						//S_w
						S_Re += 0.5*w*epsilon(a,b)*(x-y)*(phia[0]*phib_mt[0]-phia[1]*phib_mt[1]);
						S_Re += -0.5*w*epsilon(a,b)*x*(phia[0]*phib_my[0]-phia[1]*phib_my[1]);
						S_Re += 0.5*w*epsilon(a,b)*y*(phia[0]*phib_mx[0]-phia[1]*phib_mx[1]);
						S_Re += 0.5*w*(x-y)*(phia[0]*phib_mt[1]+phia[1]*phib_mt[0]);
						S_Re += -0.5*w*x*(phia[0]*phib_my[1]+phia[1]*phib_my[0]);
						S_Re += 0.5*w*y*(phia[0]*phib_mx[1]+phia[1]*phib_mx[0]);
						S_Im += 0.5*w*epsilon(a,b)*(x-y)*(phia[0]*phib_mt[1]+phia[1]*phib_mt[0]);
						S_Im += -0.5*w*epsilon(a,b)*x*(phia[0]*phib_my[1]+phia[1]*phib_my[0]);
						S_Im += 0.5*w*epsilon(a,b)*y*(phia[0]*phib_mx[1]+phia[1]*phib_mx[0]);
						S_Im += 0.5*w*(y-x)*(phia[0]*phib_mt[0]-phia[1]*phib_mt[1]);
						S_Im += 0.5*w*x*(phia[0]*phib_my[0]-phia[1]*phib_my[1]);
						S_Im += -0.5*w*y*(phia[0]*phib_mx[0]-phia[1]*phib_mx[1]);
					}//end S_trap and S_w part
					//S_int
					S_Re += 0.25*l*(phia[1]*phia[1]*phib_mt[0]*phib_mt[0]-phia[0]*phia[0]*phib_mt[0]*phib_mt[0]);
					S_Re += 0.25*l*(phia[0]*phia[0]*phib_mt[1]*phib_mt[1]-phia[1]*phia[1]*phib_mt[1]*phib_mt[1]);
					S_Re += 0.5*l*(phia[1]*phia_mt[1]*phib[1]*phib_mt[1] - phia[1]*phia_mt[1]*phib[0]*phib_mt[0]);
					S_Re += 0.5*l*(phia[0]*phia_mt[0]*phib[0]*phib_mt[0]-phia[0]*phia_mt[0]*phib[1]*phib_mt[1]);
					S_Re += -0.5*l*(phia[0]*phia_mt[1]*phib[0]*phib_mt[1]+phia[1]*phia_mt[0]*phib[1]*phib_mt[0]);
					S_Re += -0.5*l*(phia[0]*phia_mt[1]*phib[1]*phib_mt[0]+phia[1]*phia_mt[0]*phib[0]*phib_mt[1]);
					S_Re += l*(phia[0]*phia[1]*phib_mt[0]*phib_mt[1]);
					S_Re += 0.25*l*epsilon(a,b)*(phia[0]*phia_mt[0]*phia_mt[0]*phib[1] - phia[0]*phia_mt[1]*phia_mt[1]*phib[1]);
					S_Re += 0.25*l*epsilon(a,b)*(phia[1]*phia_mt[0]*phia_mt[0]*phib[0] - phia[1]*phia_mt[1]*phia_mt[1]*phib[0]);
					S_Re += 0.25*l*epsilon(a,b)*(phia[1]*phia[1]*phia_mt[0]*phib_mt[1] - phia[0]*phia[0]*phia_mt[0]*phib_mt[1]);
					S_Re += 0.25*l*epsilon(a,b)*(phia[1]*phia[1]*phia_mt[1]*phib_mt[0] - phia[0]*phia[0]*phia_mt[1]*phib_mt[0]);
					S_Re += 0.5*l*epsilon(a,b)*(phia[0]*phia_mt[1]*phia_mt[0]*phib[0] - phia[1]*phia_mt[1]*phia_mt[0]*phib[1]);
					S_Re += 0.5*l*epsilon(a,b)*(phia[1]*phia_mt[0]*phib[1]*phib_mt[1] - phia[1]*phia_mt[0]*phib[0]*phib_mt[0]);
					S_Im += 0.5*l*(phia[0]*phia[1]*phib_mt[1]*phib_mt[1]-phia[0]*phia[1]*phib_mt[0]*phib_mt[0]);
					S_Im += 0.5*l*(phia[1]*phia[1]*phib_mt[0]*phib_mt[1]-phia[0]*phia[0]*phib_mt[0]*phib_mt[1]);
					S_Im += 0.5*l*(phia[1]*phia_mt[0]*phib[0]*phib_mt[0]-phia[0]*phia_mt[1]*phib[1]*phib_mt[1]);
					S_Im += 0.5*l*(phia[0]*phia_mt[0]*phib[0]*phib_mt[1]-phia[1]*phia_mt[1]*phib[1]*phib_mt[0]);
					S_Im += 0.5*l*(phia[0]*phia_mt[1]*phib[0]*phib_mt[0]-phia[1]*phia_mt[0]*phib[1]*phib_mt[1]);
					S_Im += 0.5*l*(phia[0]*phia_mt[0]*phib[1]*phib_mt[0]-phia[1]*phia_mt[1]*phib[0]*phib_mt[1]);
					S_Im += 0.25*l*epsilon(a,b)*(phia[0]*phia[0]*phia_mt[0]*phib_mt[0] - phia[0]*phia[0]*phia_mt[1]*phib_mt[1]);
					S_Im += 0.25*l*epsilon(a,b)*(phia[0]*phia_mt[1]*phia_mt[1]*phib[0] - phia[0]*phia_mt[0]*phia_mt[0]*phib[0]);
					S_Im += 0.25*l*epsilon(a,b)*(phia[1]*phia_mt[0]*phia_mt[0]*phib[1] - phia[1]*phia_mt[1]*phia_mt[1]*phib[1]);
					S_Im += 0.25*l*epsilon(a,b)*(phia[1]*phia[1]*phia_mt[1]*phib_mt[1]-phia[1]*phia[1]*phia_mt[0]*phib_mt[0]);
					S_Im += 0.5*l*epsilon(a,b)*(phia[1]*phia_mt[1]*phia_mt[0]*phib[0] + phia[0]*phia_mt[1]*phia_mt[0]*phib[1]);
					S_Im += -0.5*l*epsilon(a,b)*(phia[1]*phia[0]*phia_mt[0]*phib_mt[1] + phia[1]*phia[0]*phia_mt[1]*phib_mt[0]);
				}//loop over b
			}//loop over a
		}//loop over space
	}//loop over time
	S_vec[0] = S_Re;
	S_vec[1] = S_Im;
	return S_vec;
}