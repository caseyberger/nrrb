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
//function declarations:
double K_a_Re(double m, double l, double w, double w_t,int a, double dtau, double *** Lattice, int dim, int Nx, int Nt, int i, int t, double mu);
double K_a_Im(double m, double l, double w, double w_t,int a, double dtau, double *** Lattice, int dim, int Nx, int Nt, int i, int t, double mu);
double epsilon(int a, int b);
std::vector<int> positive_step(int dim, int dir, int i, int t, int Nx, int Nt);
std::vector<int> negative_step(int dim, int dir, int i, int t, int Nx, int Nt);
std::vector<double> phi_a(double *** Lattice, int i, int t, int a, bool new_lat);
void test_step_functions(double *** Lattice, int size, int dim, int Nx, int Nt);

void Langevin_evolution(double m, double l, double w, double w_t, double dtau, double *** Lattice, int size, int dim, int Nx, int Nt, double mu, double eps){
	//std::default_random_engine generator;
	double random_num_gen_time = 0.;
	double Langevin_calc_time = 0.;
	double Lattice_update_time = 0.;
	//std::random_device rd;
	//std::mt19937 mt(rd());
	int seed = 61;
	std::mt19937 generator(seed);
	//evolve the lattice using the old lattice copy
	//#pragma omp parallel
	//#pragma omp for collapse(2)
	//std::normal_distribution<double> noise1(0.,sqrt(2.));
	//std::uniform_real_distribution<double> noise1(0.,1.);
	//noise1.reset();
	//for (int num = 0; num < 10; num++){	
	//	double eta_1 = noise1(generator);
	//	std::cout << eta_1 << std::endl;
	//}

	for(int i = 0; i<size; i++){
		for (int t = 0; t<Nt;t++){
			clock_t rn_0 = clock();
			//std::normal_distribution<double> noise1(0.,sqrt(2.));
			//std::normal_distribution<double> noise2(0.,sqrt(2.));
			std::uniform_real_distribution<double> noise1(-1.,1.);
			std::uniform_real_distribution<double> noise2(-1.,1.);
			double eta_1 = noise1(generator);
			double eta_2 = noise2(generator);
			clock_t rn_f = clock();
			random_num_gen_time += float(rn_f - rn_0)/CLOCKS_PER_SEC;
			clock_t CL_0 = clock();
			//phi1_Re
			Lattice[i][t][0] = Lattice[i][t][4] + eps*K_a_Re(m,l,w,w_t,1,dtau,Lattice,dim,Nx,Nt,i,t,mu);// + sqrt(eps)*eta_1;
			//Lattice[i][t][0] = eta_1;

			//phi1_Im
			Lattice[i][t][1] = Lattice[i][t][5] + eps*K_a_Im(m,l,w,w_t,1,dtau,Lattice,dim,Nx,Nt,i,t,mu);
			//Lattice[i][t][1] = 0.;

			//phi2_Re
			Lattice[i][t][2] = Lattice[i][t][6] + eps*K_a_Re(m,l,w,w_t,2,dtau,Lattice,dim,Nx,Nt,i,t,mu);// + sqrt(eps)*eta_2;
			//Lattice[i][t][2] = eta_2;

			//phi2_Im
			Lattice[i][t][3] = Lattice[i][t][7] + eps*K_a_Im(m,l,w,w_t,2,dtau,Lattice,dim,Nx,Nt,i,t,mu);
			//Lattice[i][t][3] = 0.;

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

//COMPUTING THE REAL DRIFT FUNCTION
double K_a_Re(double m, double l,double w, double w_t, int a, double dtau, double *** Lattice, int dim, int Nx, int Nt, int i, int t, double mu){
	double Ka = 0.0;
	//doing -K_{a}^{R} and then returning -Ka at the end for simplicity
	//modify mu, m, w, wtr, and l by dtau
	mu = dtau*mu;
	m = m/dtau;
	w_t = dtau*w_t;
	w = dtau*w;
	l = dtau*l;
	//define bool new_lat = true if we are using the new lattice 
	//and bool old_lat = false if we are using the old one
	bool old_lat = false; 
	std::vector<int> vp = positive_step(dim, 4, i, t, Nx, Nt);
	std::vector<int> vn = negative_step(dim, 4, i, t, Nx, Nt);
	int t_p = vp[1];
	int t_n = vn[1];
	std::vector<double> phia = phi_a(Lattice,i,t,a,old_lat);//phi_a
	std::vector<double> phia_pt = phi_a(Lattice,i,t_p,a,old_lat);//phi_{a,t+1}
	std::vector<double> phia_mt = phi_a(Lattice,i,t_n,a,old_lat); //phi_{a,t-1}
	
	//CHEMICAL POTENTIAL + TIME DERIVATIVE PART OF Ka
	//\phi_{a,r}^{R} - e^{\mu}/2 (\phi_{a,r-t}^{R} + \phi_{a,r+t}^{R})
	//no special conditions here because we have periodic boundary conditions in time
	Ka += phia[0]- 0.5*exp(mu)*phia_mt[0] - 0.5*exp(mu)*phia_pt[0];
	for (int b=1; b<=2; b++){
		std::vector<double> phib_mt = phi_a(Lattice,i,t_n,b,old_lat); //phi_{b,t-1}
		std::vector<double> phib_pt = phi_a(Lattice,i,t_p,b,old_lat); //phi_{b,t+1}
		//\eps_{ab} e^{\mu}/2 (\phi_{b,r-t}^{I}-\phi_{b,r+t}^{I})
		Ka += epsilon(a,b)*0.5*exp(mu)*phib_mt[1];
		Ka -= epsilon(a,b)*0.5*exp(mu)*phib_pt[1];
	}//last checked on 01.31.19
	
	//SPATIAL DERIVATIVE PART OF Ka
	//To test 0+1 dimensional case, just skip this loop
	//(d/m)\phi_{a}^{R} 
	Ka += 1.*dim*phia[0]/m;
	for (int d = 1; d <= dim; d++){//loop over adjacent spatial sites!
		std::vector<int> vp = positive_step(dim, d, i, t, Nx, Nt);
		std::vector<int> vn = negative_step(dim, d, i, t, Nx, Nt);
		int i_p = vp[0];
		int i_n = vn[0];
		std::vector<double> phia_pi = phi_a(Lattice,i_p,t,a,old_lat);//phi_{a,i+1} --> i = x, y, z
		std::vector<double> phia_mi = phi_a(Lattice,i_n,t,a,old_lat);//phi_{a,i-1} --> i = x, y, z
		// -(1/2m)\phi_{a,i+1}^{R} - (1/2m)\phi_{a,ii1}^{R}
		//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
		if (!(i_p == 9999. || i_n == 9999.)){
				Ka += -0.5*phia_pi[0]/m -0.5*phia_mi[0]/m;
			}
	}//updated on 10.14.18 - part mixing a and b removed per notes
	//updated on 1.31.19
	
	//ROTATIONAL AND TRAPPING PARTS OF Ka - ONLY WORKS IN 2D AND IF OMEGA > 0.
	if (dim ==2){
		int xint = i%(Nx);
		int yint = ((i - xint)/Nx)%Nx;
		//shift x, y by (Nx-1)/2
		double x = 1.0*xint - 0.5*(Nx-1);
		double y = 1.0*yint - 0.5*(Nx-1);
		double r2 = x*x+y*y;
		//if (t == 40){
		//	std::cout << "{x,y} = {" << xint << "," << yint << "}, "<< "(x - x_c) = " << x << ", (y - y_c) = " << y << std::endl;
		//}
		//double r2 = 1.;
		std::vector<int> vpx = positive_step(dim, 1, i, t, Nx, Nt);
		std::vector<int> vnx = negative_step(dim, 1, i, t, Nx, Nt);
		std::vector<int> vpy = positive_step(dim, 2, i, t, Nx, Nt);
		std::vector<int> vny = negative_step(dim, 2, i, t, Nx, Nt);
		std::vector<double> phia_px = phi_a(Lattice,vpx[0],t_p,a,old_lat);//phi_{a,r+x+t}
		std::vector<double> phia_mx = phi_a(Lattice,vnx[0],t_n,a,old_lat);//phi_{a,r-x-t}
		std::vector<double> phia_py = phi_a(Lattice,vpy[0],t_p,a,old_lat);//phi_{a,r+y+t}
		std::vector<double> phia_my = phi_a(Lattice,vny[0],t_n,a,old_lat);//phi_{a,r-y-t}
		//TRAPPING PIECE
		// 0.25*m*w_t^{2}*r2*(\phi_{a,r+t}^{R}+\phi_{a,r+t}^{R})
		Ka += 0.25*m*w_t*w_t*r2*phia_pt[0]+0.25*m*w_t*w_t*r2*phia_mt[0];
		//ROTATIONAL PIECE
		// w (x-y)(\phi_{a,r+t}^{I}+\phi_{a,r-t}^{I})
		Ka += 0.5*w*x*phia_pt[1] - 0.5*w*y*phia_pt[1] + 0.5*w*x*phia_mt[1] - 0.5*w*y*phia_mt[1];
		// (w/2) y (\phi_{a,r+x}^{I} + \phi_{a,r-x}^{I}) - (w/2) x (\phi_{a,r+y}^{I} + \phi_{a,r-y}^{I}) 
		//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
		if (!(vpx[0] == 9999. || vnx[0] == 9999.)){
			Ka += 0.5*w*y*phia_px[1] + 0.5*w*y*phia_mx[1];
		}
		if (!(vpy[0] == 9999. || vny[0] == 9999.)){
			Ka += -0.5*w*x*phia_py[1] - 0.5*w*x*phia_my[1];
		}
		for (int b=1; b<=2; b++){
			std::vector<double> phib_mx = phi_a(Lattice,vnx[0],t_n,b,old_lat);//phi_{b,r-x-t}
			std::vector<double> phib_my= phi_a(Lattice,vny[0],t_n,b,old_lat);//phi_{b,r-y-t}
			std::vector<double> phib_px = phi_a(Lattice,vpx[0],t_p,b,old_lat);//phi_{b,r+x+t}
			std::vector<double> phib_py= phi_a(Lattice,vpy[0],t_p,b,old_lat);//phi_{b,r+y+t}
			std::vector<double> phib_mt = phi_a(Lattice,i,t_n,b,old_lat); //phi_{b,t-1}
			std::vector<double> phib_pt = phi_a(Lattice,i,t_p,b,old_lat); //phi_{b,t+1}
			//\eps_{ab} (w/2) (y (\phi_{b,r-x}^{R} + \phi_{b,r+x}^{R})- x (\phi_{b,r-y}^{R}+\phi_{b,r+y}^{R}))
			Ka += -0.5*epsilon(a,b)*w*x*phib_py[0] - 0.5*epsilon(a,b)*w*x*phib_my[0];
			Ka += 0.5*epsilon(a,b)*w*y*phib_mx[0] + 0.5*epsilon(a,b)*w*y*phib_px[0];
			//\eps_{ab} (w/2) (x-y)(\phi_{b,r+t}^{R}+\phi_{b,r-t}^{R})
			Ka += 0.5*epsilon(a,b)*w*x*phib_pt[0] +0.5*epsilon(a,b)*w*x*phib_mt[0];
			Ka += -0.5*epsilon(a,b)*w*y*phib_pt[0] -0.5*epsilon(a,b)*w*y*phib_mt[0];
			//TRAPPING PIECE
			// -0.25*m*w_t^{2} r2*\eps_{ab}(\phi_{b,r+t}^{I}+\phi_{b,r-t}^{I})
			Ka += -0.25*m*w_t*w_t*r2*epsilon(a,b)*phib_pt[1] -0.25*m*w_t*w_t*r2*epsilon(a,b)*phib_mt[1];
		}
	}
	//last updated 2.20.19

	//INTERACTION PART OF Ka
	for (int b=1; b<=2; b++){
		std::vector<double> phib = phi_a(Lattice,i,t,b,old_lat);//\phi_{b}
		std::vector<double> phib_mt = phi_a(Lattice,i,t_n,b,old_lat);//phi_{b,t-1}
		std::vector<double> phib_pt = phi_a(Lattice,i,t_p,b,old_lat);//phi_{b,t+1}
		Ka += l*phib[0]*phia_mt[0]*phib_mt[0] - l*phib[0]*phia_mt[1]*phib_mt[1];
		Ka += l*phib[0]*phia_pt[0]*phib_pt[0] - l*phib[0]*phia_pt[1]*phib_pt[1];
		Ka += -1.0*l*phib[1]*phia_mt[0]*phib_mt[1] - l*phib[1]*phia_mt[1]*phib_mt[0];
		Ka += -1.0*l*phib[1]*phia_pt[0]*phib_pt[1] - l*phib[1]*phia_pt[1]*phib_pt[0];
		Ka += 0.5*l*phia[0]*phib_mt[1]*phib_mt[1]-0.5*l*phia[0]*phib_mt[0]*phib_mt[0];
		Ka += 0.5*l*phia[0]*phib_pt[1]*phib_pt[1]-0.5*l*phia[0]*phib_pt[0]*phib_pt[0];
		Ka += l*phia[1]*phib_mt[1]*phib_mt[0]+l*phia[1]*phib_pt[1]*phib_pt[0];
		Ka += -1.*epsilon(a,b)*l*phia[0]*phia_mt[0]*phib_mt[1]-epsilon(a,b)*l*phia[0]*phia_mt[1]*phib_mt[0];
		Ka += epsilon(a,b)*l*phia[0]*phia_pt[0]*phib_pt[1]+epsilon(a,b)*l*phia[0]*phia_pt[1]*phib_pt[0];
		Ka += epsilon(a,b)*l*phia[1]*phia_mt[1]*phib_mt[1]-epsilon(a,b)*l*phia[1]*phia_mt[0]*phib_mt[0];
		Ka += epsilon(a,b)*l*phia[1]*phia_pt[0]*phib_pt[0]-epsilon(a,b)*l*phia[1]*phia_pt[1]*phib_pt[1];
		Ka += 0.5*epsilon(a,b)*l*phib[1]*phia_mt[0]*phia_mt[0]-0.5*epsilon(a,b)*l*phib[1]*phia_mt[1]*phia_mt[1];
		Ka += 0.5*epsilon(a,b)*l*phib[1]*phia_pt[1]*phia_pt[1]-0.5*epsilon(a,b)*l*phib[1]*phia_pt[0]*phia_pt[0];
		Ka += epsilon(a,b)*l*phib[0]*phia_mt[1]*phia_mt[0]-epsilon(a,b)*l*phib[0]*phia_pt[1]*phia_pt[0];
		Ka += 0.5*epsilon(a,b)*l*phia[1]*phia_mt[1]*phia_mt[1]-0.5*epsilon(a,b)*l*phia[1]*phia_mt[0]*phia_mt[0];
		Ka += 0.5*epsilon(a,b)*l*phia[1]*phia_pt[0]*phia_pt[0]-0.5*epsilon(a,b)*l*phia[1]*phia_pt[1]*phia_pt[1];
		Ka += epsilon(a,b)*l*phia[0]*phia_pt[1]*phia_pt[0]-epsilon(a,b)*l*phia[0]*phia_mt[1]*phia_mt[0];
	}//last update 2.22.19 to make a gauge interaction
	/*
	std::stringstream mu_stream, w_stream, dtau_stream;
	mu_stream << std::fixed << std::setprecision(3) << mu; //truncate mu for filename
	std::string str_mu = mu_stream.str();
	w_stream << std::fixed << std::setprecision(3) << w; //truncate w for filename
	std::string str_w = w_stream.str();
	dtau_stream << std::fixed << std::setprecision(3) << dtau; //truncate dtau for filename
	std::string str_dtau = dtau_stream.str();
	std::string str_N = std::to_string(Nx);
	std::string str_Nt = std::to_string(Nt);
	std::string fname = "KaRe_mu_" + str_mu + "_w_"+str_w+"_Nx_"+str_N+"_Nt_"+str_Nt+"_dt_"+str_dtau+".txt";//output filename
	std::ofstream fout; //output stream
	fout.open(fname, std::ios::out|std::ios::app); //open the file
	//check if files are open
    if (fout.is_open()){
    	fout << std::setw(3) << i;
    	fout << std::setw(10) <<  Ka << std::endl;
    	fout.close();
    }
    else{
		std::cerr << "Unable to open file " << std::endl;
        exit(10);
	}*/
	return -1.*Ka;
}


//COMPUTING THE IMAGINARY DRIFT FUNCTION
double K_a_Im(double m, double l, double w, double w_t, int a, double dtau, double *** Lattice, int dim, int Nx, int Nt, int i, int t, double mu){
	double Ka = 0.0;
	//doing -K_{a}^{I} and then returning -Ka at the end for simplicity
	//modify mu, m, w, wtr, and l by dtau
	mu = dtau*mu;
	m = m/dtau;
	w = dtau*w;
	w_t = dtau*w_t;
	l = dtau*l;
	//define bool new_lat = true if we are using the new lattice 
	//and bool old_lat = false if we are using the old one
	bool old_lat = false; 
	std::vector<int> vp = positive_step(dim, 4, i, t, Nx, Nt);
	std::vector<int> vn = negative_step(dim, 4, i, t, Nx, Nt);
	int t_p = vp[1];
	int t_n = vn[1];
	std::vector<double> phia = phi_a(Lattice,i,t,a,old_lat);//phi_a
	std::vector<double> phia_pt = phi_a(Lattice,i,t_p,a,old_lat);//phi_{a,t+1}
	std::vector<double> phia_mt = phi_a(Lattice,i,t_n,a,old_lat); //phi_{a,t-1}
	//CHEMICAL POTENTIAL + TIME DERIVATIVE PART OF Ka
	//\phi_{a,r}^{I} - e^{\mu}/2 (\phi_{a,r+t}^{I} + \phi_{a,r-t}^{I})
	//no special boundary concerns because we have periodic boundary conditions in time
	Ka += phia[1] - 0.5*exp(mu)*phia_mt[1] - 0.5*exp(mu)*phia_pt[1];
	for (int b=1; b<=2; b++){
		std::vector<double> phib_mt = phi_a(Lattice,i,t_n,b,old_lat); //phi_{b,t-1}
		std::vector<double> phib_pt = phi_a(Lattice,i,t_p,b,old_lat); //phi_{b,t+1}
		// -\eps_{ab} e^{\mu}/2 ( \phi_{b,r-t}^{R} - \phi_{b,r+t}^{R})
		Ka -= 0.5*epsilon(a,b)*exp(mu)*phib_mt[0];
		Ka += 0.5*epsilon(a,b)*exp(mu)*phib_pt[0];
	}//checked on 10.14.18
	//checked again on 1.31.19

	//SPATIAL DERIVATIVE PART OF Ka
	//To test 0+1 dimensional case, just skip this loop
	//(d/m)\phi_{a}^{I} 
	Ka += 1.0*dim*phia[1]/m;
	for (int d = 1; d <= dim; d++){//loop over adjacent spatial sites!
		std::vector<int> vp = positive_step(dim, d, i, t, Nx, Nt);
		std::vector<int> vn = negative_step(dim, d, i, t, Nx, Nt);
		int i_p = vp[0];
		int i_n = vn[0];
		std::vector<double> phia_pi = phi_a(Lattice,i_p,t,a,old_lat);//phi_{a,i+1} --> i = x, y, z
		std::vector<double> phia_mi = phi_a(Lattice,i_n,t,a,old_lat);//phi_{a,i-1} --> i = x, y, z
		//-(1/2m)(\phi_{a,i+1}^{I} +\phi_{a,i-1}^{I})
		//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
		if (!(i_p == 9999. || i_n == 9999.)){
			Ka += -0.5*phia_pi[1]/m - 0.5*phia_mi[1]/m;
		}
	}//updated on 10.14.18 - part mixing a and b removed per notes
	//updated on 1.31.19

	//ROTATIONAL AND TRAPPING PART OF Ka - ONLY WORKS IN 2D AND IF OMEGA > 0.
	if (dim ==2){
		int xint = i%(Nx);
		int yint = ((i - xint)/Nx)%Nx;
		//shift x, y by (Nx-1)/2
		double x = 1.0*xint - 0.5*(Nx-1);
		double y = 1.0*yint - 0.5*(Nx-1);
		double r2 = x*x+y*y;
		//double r2 = 1.;
		std::vector<int> vpx = positive_step(dim, 1, i, t, Nx, Nt);
		std::vector<int> vnx = negative_step(dim, 1, i, t, Nx, Nt);
		std::vector<int> vpy = positive_step(dim, 2, i, t, Nx, Nt);
		std::vector<int> vny = negative_step(dim, 2, i, t, Nx, Nt);
		std::vector<double> phia_px = phi_a(Lattice,vpx[0],t_p,a,old_lat);//phi_{a,r+x+t}
		std::vector<double> phia_mx = phi_a(Lattice,vnx[0],t_n,a,old_lat);//phi_{a,r-x-t}
		std::vector<double> phia_py = phi_a(Lattice,vpy[0],t_p,a,old_lat);//phi_{a,r+y+t}
		std::vector<double> phia_my = phi_a(Lattice,vny[0],t_n,a,old_lat);//phi_{a,r-y-t}
		//TRAPPING PIECE
		// 0.25*m*w_t^{2}*r2*(\phi_{a,r+t}^{I}+\phi_{a,r+t}^{I})
		Ka += 0.25*m*w_t*w_t*r2*phia_pt[1]+0.25*m*w_t*w_t*r2*phia_mt[1];
		//ROTATIONAL PIECE
		// w(y-x)\phi_{a}^{R}
		Ka += 0.5*w*y*phia_pt[0] - 0.5*w*x*phia_pt[0] + 0.5*w*y*phia_mt[0] - 0.5*w*x*phia_mt[0];
		// (w/2) x (\phi_{a,r+y}^{R} + \phi_{a,r-y}^{R})-(w/2) y (\phi_{a,r+x}^{R} + \phi_{a,r-x}^{R})
		//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
		if (!(vpx[0] == 9999. || vnx[0] == 9999.)){
			Ka += -0.5*w*y*phia_px[0] - 0.5*w*y*phia_mx[0];
		}
		if (!(vpy[0] == 9999. || vny[0] == 9999.)){
			Ka +=  0.5*w*x*phia_py[0] + 0.5*w*x*phia_my[0];
		}
		for (int b=1; b<=2; b++){
			std::vector<double> phib_mx = phi_a(Lattice,vnx[0],t_n,b,old_lat);//phi_{b,r-x-t}
			std::vector<double> phib_my= phi_a(Lattice,vny[0],t_n,b,old_lat);//phi_{b,r-y-t}
			std::vector<double> phib_px = phi_a(Lattice,vpx[0],t_p,b,old_lat);//phi_{b,r+x+t}
			std::vector<double> phib_py= phi_a(Lattice,vpy[0],t_p,b,old_lat);//phi_{b,r+y+t}
			std::vector<double> phib_mt = phi_a(Lattice,i,t_n,b,old_lat); //phi_{b,t-1}
			std::vector<double> phib_pt = phi_a(Lattice,i,t_p,b,old_lat); //phi_{b,t+1}
		//\eps_{ab} (w/2) (y (\phi_{b,r-x}^{I}+\phi_{b,r+x}^{I}) - x (\phi_{b,r-y}^{I}+\phi_{b,r+y}^{I}) )
			Ka += -0.5*epsilon(a,b)*w*x*phib_py[1] - 0.5*epsilon(a,b)*w*x*phib_my[1];
			Ka += 0.5*epsilon(a,b)*w*y*phib_mx[1] + 0.5*epsilon(a,b)*w*y*phib_px[1];
			//\eps_{ab} (w/2) (x-y)(\phi_{b,r+t}^{I}+\phi_{b,r-t}^{I})
			Ka += 0.5*epsilon(a,b)*w*x*phib_pt[1] +0.5*epsilon(a,b)*w*x*phib_mt[1];
			Ka += -0.5*epsilon(a,b)*w*y*phib_pt[1] -0.5*epsilon(a,b)*w*y*phib_mt[1];
			//TRAPPING PIECE
			// 0.25*w_t^{2}*r2*\eps_{ab}(\phi_{b,r+t}^{R}+\phi_{b,r+t}^{R})
			Ka += 0.25*m*w_t*w_t*r2*epsilon(a,b)*phib_pt[0] + 0.25*m*w_t*w_t*r2*epsilon(a,b)*phib_mt[0];
		}
	}
	//last updated 2.20.19 
	
	//INTERACTION PART OF Ka
	for (int b=1; b<=2; b++){
		std::vector<double> phib = phi_a(Lattice,i,t,b,old_lat);//\phi_{b}
		std::vector<double> phib_mt = phi_a(Lattice,i,t_n,b,old_lat);//phi_{b,t-1}
		std::vector<double> phib_pt = phi_a(Lattice,i,t_p,b,old_lat);//phi_{b,t+1}
		Ka += l*phib[0]*phia_mt[0]*phib_mt[1] + l*phib[0]*phia_mt[1]*phib_mt[0];
		Ka += l*phib[0]*phia_pt[0]*phib_pt[1] + l*phib[0]*phia_pt[1]*phib_pt[0];
		Ka += -1.0*l*phib[1]*phia_mt[1]*phib_mt[1] + l*phib[1]*phia_mt[0]*phib_mt[0];
		
		Ka += -1.0*l*phib[1]*phia_pt[1]*phib_pt[1] + l*phib[1]*phia_pt[0]*phib_pt[0];
		
		Ka += 0.5*l*phia[1]*phib_mt[1]*phib_mt[1]-0.5*l*phia[1]*phib_mt[0]*phib_mt[0];
		
		Ka += 0.5*l*phia[1]*phib_pt[1]*phib_pt[1]-0.5*l*phia[1]*phib_pt[0]*phib_pt[0];
		
		Ka += -1.0*l*phia[0]*phib_mt[1]*phib_mt[0]-l*phia[0]*phib_pt[1]*phib_pt[0];
		
		Ka += -1.*epsilon(a,b)*l*phia[1]*phia_mt[0]*phib_mt[1]-epsilon(a,b)*l*phia[1]*phia_mt[1]*phib_mt[0];
		
		Ka += epsilon(a,b)*l*phia[1]*phia_pt[0]*phib_pt[1]+epsilon(a,b)*l*phia[1]*phia_pt[1]*phib_pt[0];
		
		Ka += epsilon(a,b)*l*phia[0]*phia_mt[0]*phib_mt[0]-epsilon(a,b)*l*phia[0]*phia_mt[1]*phib_mt[1];
		
		Ka += epsilon(a,b)*l*phia[0]*phia_pt[1]*phib_pt[1]-epsilon(a,b)*l*phia[0]*phia_pt[0]*phib_pt[0];
		
		Ka += 0.5*epsilon(a,b)*l*phib[0]*phia_mt[1]*phia_mt[1]-0.5*epsilon(a,b)*l*phib[0]*phia_mt[0]*phia_mt[0];
		
		Ka += 0.5*epsilon(a,b)*l*phib[0]*phia_pt[0]*phia_pt[0]-0.5*epsilon(a,b)*l*phib[0]*phia_pt[1]*phia_pt[1];
		
		Ka += epsilon(a,b)*l*phib[1]*phia_mt[1]*phia_mt[0]-epsilon(a,b)*l*phib[1]*phia_pt[1]*phia_pt[0];
		
		Ka += 0.5*epsilon(a,b)*l*phia[0]*phia_mt[0]*phia_mt[0]-0.5*epsilon(a,b)*l*phia[0]*phia_mt[1]*phia_mt[1];
		
		Ka += 0.5*epsilon(a,b)*l*phia[0]*phia_pt[1]*phia_pt[1]-0.5*epsilon(a,b)*l*phia[0]*phia_pt[0]*phia_pt[0];
		
		Ka += epsilon(a,b)*l*phia[1]*phia_pt[1]*phia_pt[0]-epsilon(a,b)*l*phia[1]*phia_mt[1]*phia_mt[0];
	}//last update 2.22.19 to make a gauge interaction
	/*
	std::stringstream mu_stream, w_stream, dtau_stream;
	mu_stream << std::fixed << std::setprecision(3) << mu; //truncate mu for filename
	std::string str_mu = mu_stream.str();
	w_stream << std::fixed << std::setprecision(3) << w; //truncate w for filename
	std::string str_w = w_stream.str();
	dtau_stream << std::fixed << std::setprecision(3) << dtau; //truncate dtau for filename
	std::string str_dtau = dtau_stream.str();
	std::string str_N = std::to_string(Nx);
	std::string str_Nt = std::to_string(Nt);
	std::string fname = "KaIm_mu_" + str_mu + "_w_"+str_w+"_Nx_"+str_N+"_Nt_"+str_Nt+"_dt_"+str_dtau+".txt";//output filename
	std::ofstream fout; //output stream
	fout.open(fname, std::ios::out|std::ios::app); //open the file
	//check if files are open
    if (fout.is_open()){
    	fout << std::setw(3) << i;
    	fout << std::setw(10) <<  Ka << std::endl;
    	fout.close();
    }
    else{
		std::cerr << "Unable to open file " << std::endl;
        exit(10);
	}*/
	return -1.*Ka;
}

//ANTI-SYMMETRIC TENSOR
double epsilon(int a, int b){
	if (a == b){return 0;}
	else if (a == 1 and b ==2){return 1.;}
	else if (a ==2 and b == 1){return -1.;}
	else {return 0.;}//added 10.30.17
}

//FUNCTION TO MOVE ONE STEP FORWARD IN ONE DIRECTION
std::vector<int> positive_step(int dim, int dir, int i, int t, int Nx, int Nt){
	//dir ranges from 1 to 4 (1 = x, 2 = y, 3 = z, 4 = t)
	//TRY making j = 9999 when we go outside the lattice. 
	//Then have phi_a check for that and return 0 if it's found
	//This means we can't go larger than 21^3 or 99^2 for the spatial lattice
	int j=i; //new position one step in the positive direction in space
	int k=t; //new position one step in the positive direction in time
	std::vector<int> v;
	if (dim == 1){
		int x = i;
		//std::cout << " positive_step, dir = " << dir << " {x,t} = " << x << ' ' << t << std::endl;
		if (dir == 4){
			if (t == Nt-1){k = 0;}
			else {k = t + 1;}
		}
		else if (dir ==1){
			if(x == Nx-1){j = 9999;}
			else{j = i+1;}
		}
		else {
			std::cerr << "Invalid direction " << dir << " for {x,t} = " << x << ' ' << t << std::endl;
			exit(10);
		}
		i = j; //update i
		t = k; //update t
		x = j;
		//std::cout << " {x,t}_new = " << x << ' ' << t << std::endl << std::endl;
	}//close loop for single spatial dimension
	if (dim == 2){
		int x = i%(Nx);
		int y = ((i - x)/Nx)%Nx;
		//std::cout << " positive_step dir = " << dir<< " {x,y,t} = " << x << ' '<< y << ' ' << t << std::endl;
		if (dir == 4){
			j = i;
			if (t == Nt-1){k = 0;}
			else {k = t + 1;}
		}
		else if (dir ==1){
			k = t;
			if(x == Nx-1){j = 9999;}
			else{j = i+ 1;}
		}
		else if (dir ==2){
			k = t;
			if(y == Nx-1){j = 9999;}
			else {j = i + Nx;}
			}
		else {
			std::cerr << "Invalid direction " << dir << " for {x,y,t} = " << x << ' '<< y << ' ' << t << std::endl;
			exit(10);
		}
		i = j;//update i
		t = k;//update t
		if (j == 9999){
			x = 9999;
			y = 9999;
		}
		else{
			x = j%(Nx);
			y = ((j - x)/Nx)%Nx;
		}
		//std::cout << " {x,y,t}_new = " << x << ' '<< y << ' ' << t << std::endl << std::endl;
	}
	if (dim ==3){
		int x = i%(Nx);
		int y = ((i - x)/Nx)%Nx;
		int z = ((i - x - (Nx*y))/(Nx*Nx))%(Nx);
		//std::cout << " positive_step, dir = " << dir<< " {x,y,z,t} = " << x << ' '<< y << ' ' << z << ' ' << t << std::endl;
		if (dir == 4){
			j = i;
			if (t == Nt-1){k = 0;}
			else {k = t + 1;}
		}
		else if (dir ==1){
			k = t;
			if(x == Nx-1){j=9999;}
			else{j = i+ 1;}
		}
		else if (dir ==2){
			k = t;
			if(y == Nx-1){j =9999;}
			else {j = i + Nx;}
			}
		else if (dir ==3){
			k = t;
			if(z == Nx-1){j = 9999;}
			else {j = i + Nx*Nx;}
			}
		else {
			std::cerr << "Invalid direction " << dir << " for {x,y,t} = " << x << ' '<< y << ' ' << t << std::endl;
			exit(10);
		}
		i = j;//update i
		t = k;//update t
		if (j == 9999){
			x = 9999;
			y = 9999;
			z = 9999;
		}
		else{
			x = j%(Nx);
			y = ((j - x)/Nx)%Nx;
			z = ((j - x - Nx*y)/(Nx*Nx))%(Nx);
		}
		//std::cout << " {x,y,z,t}_new = " << x << ' '<< y << ' ' << z << ' ' << t << std::endl;
	}
	v.push_back(j);
	v.push_back(k);
	return v;
}

//FUNCTION TO MOVE ONE STEP BACKWARD IN ONE DIRECTION
std::vector<int> negative_step(int dim, int dir, int i, int t, int Nx, int Nt){
	//dir ranges from 1 to 4 (1 = x, 2 = y, 3 = z, 4 = t)
	int j=i; //new position one step in the positive direction in space
	int k=t; //new position one step in the positive direction in time
	std::vector<int> v;
	if (dim == 1){
		int x = i;
		//std::cout << " negative_step, dir = " << dir << " {x,t} = " << x << ' ' << t << std::endl;
		if (dir == 4){
			if (t == 0){k = Nt-1;}
			else {k = t - 1;}
		}
		else if (dir ==1){
			if(x == 0){j = 9999;}
			else{j = i-1;}
		}
		else {
			std::cerr << "Invalid direction " << dir << " for {x,t} = " << x << ' ' << t << std::endl;
			exit(10);
		}
		i = j; //update i
		t = k; //update t
		x = j;
		//std::cout << " {x,t}_new = " << x << ' ' << k << std::endl << std::endl;
	}//close loop for single spatial dimension
	if (dim == 2){
		int x = i%(Nx);
		int y = ((i - x)/Nx)%Nx;
		//std::cout << " negative_step, dir = " << dir<< " {x,y,t} = " << x << ' '<< y << ' ' << t << std::endl;
		if (dir == 4){
			j = i;
			if (t == 0){k = Nt-1;}
			else {k = t - 1;}
		}
		else if (dir ==1){
			k = t;
			if(x == 0){j = 9999;}
			else{j = i - 1;}
		}
		else if (dir ==2){
			k = t;
			if(y == 0){j = 9999;}
			else {j = i - Nx;}
			}
		else {
			std::cerr << "Invalid direction " << dir << " for {x,y,t} = " << x << ' '<< y << ' ' << t << std::endl;
			exit(10);
		}
		i = j;//update i
		t = k;//update t
		if (j == 9999){
			x = 9999;
			y = 9999;
		}
		x = j%(Nx);
		y = ((j - x)/Nx)%Nx;
		//std::cout << " {x,y,t}_new = " << x << ' '<< y << ' ' << t << std::endl << std::endl;
	}
	if (dim ==3){
		//int jmax = Nx*Nx*(Nx-1);
		int x = i%(Nx);
		int y = ((i - x)/Nx)%Nx;
		int z = ((i - x - (Nx*y))/(Nx*Nx))%(Nx);
		//std::cout << " negative_step, dir = " << dir<< " {x,y,z,t} = " << x << ' '<< y << ' ' << z << ' ' << t << std::endl;
		if (dir == 4){
			j = i;
			if (t == 0){k = Nt-1;}
			else {k = t - 1;}
		}
		else if (dir ==1){
			k = t;
			if(x == 0){j = 9999;}
			else{j = i - 1;}
		}
		else if (dir ==2){
			k = t;
			if(y == 0){j = 9999;}
			else {j = i - Nx;}
			}
		else if (dir ==3){
			k = t;
			if(z == 0){j = 9999;}
			else {j = i - Nx*Nx;}
		}
		else{
			std::cerr << "Invalid direction " << dir << " for {x,y,t} = " << x << ' '<< y << ' ' << t << std::endl;
			exit(10);
		}
		i = j;//update i
		t = k;//update t
		if (j == 9999){
			x = 9999;
			y = 9999;
			z = 9999;
		}
		x = j%(Nx);
		y = ((j - x)/Nx)%Nx;
		z = ((j - x - Nx*y)/(Nx*Nx))%(Nx);
		//std::cout << " {x,y,z,t}_new = " << x << ' '<< y << ' ' << z << ' ' << t << std::endl;
	}
	v.push_back(j);
	v.push_back(k);
	return v;
}

//FUNCTION TO GET FIELD FROM LATTICE
std::vector<double> phi_a(double *** Lattice, int i, int t, int a, bool new_lat){
	std::vector<double> phia;
	if (a == 1 and new_lat == true){
		if (i == 9999){//if we are looking outside the spatial lattice, the field = 0.
			phia.push_back(0.);
			phia.push_back(0.);
		}
		else{
			phia.push_back(Lattice[i][t][0]);
			phia.push_back(Lattice[i][t][1]);
		}
	}
	else if (a == 2 and new_lat == true){
		if (i == 9999){
			phia.push_back(0.);
			phia.push_back(0.);
		}
		else{
			phia.push_back(Lattice[i][t][2]);
			phia.push_back(Lattice[i][t][3]);
		}
	}
	else if (a == 1 and new_lat == false){
		if (i == 9999){
			phia.push_back(0.);
			phia.push_back(0.);
		}
		else{
			phia.push_back(Lattice[i][t][4]);
			phia.push_back(Lattice[i][t][5]);
		}
	}
	else if (a == 2 and new_lat == false){
		if (i == 9999){
			phia.push_back(0.);
			phia.push_back(0.);
		}
		else{
			phia.push_back(Lattice[i][t][6]);
			phia.push_back(Lattice[i][t][7]);
		}
	}
	else {
		std::cerr << "Invalid index " << a << " for phi" << std::endl;
		exit(10);
	}
	return phia;
}

void test_step_functions(double *** Lattice, int size, int dim, int Nx, int Nt){
	//prep file
	std::stringstream stream;
	stream << std::fixed << std::setprecision(1) << Nx; //truncate mu for filename
	std::string str_Nx = stream.str();
	std::string fname = "lattice_sites.txt";//output filename
	std::ofstream fout; //output stream
	fout.open(fname, std::ios::out|std::ios::app); //open the file
	//check if files are open
    if (fout.is_open()){
    	//touch every site on the lattice and print its location, 
    	//plus a step forward and back in each direction
    	fout << "Checking forward and backwards step functions for lattice of Nx = " << str_Nx << std::endl;
		int count = 0;
		for (int i = 0; i < size; i++){
			for (int t = 0; t < Nt; t++){
				if (dim == 3){
					int x = i%(Nx);
					int y = ((i - x)/Nx)%Nx;
					int z = ((i - x - (Nx*y))/(Nx*Nx))%(Nx);
					//initialize x, y, and z
					//print i,t and x, y, z, t
					fout << "{i,t} = " << "{" << i << "," << t << "}" << std::endl;
					fout << "{x,y,z,t} = " << "{" << x << "," << y << "," << z << "," << t << "}" << std::endl;
					//std::vector<int> negative_step(int dim, int dir, int i, int t, int Nx, int Nt)
					for (int d = 1; d <= dim+1; d++){
						std::vector<int> vp = positive_step(dim,d,i,t,Nx,Nt);
						int ip = vp[0];
						int xp,yp,zp;
						if (ip == 9999){
							xp = 9999;
							yp = 9999;
							zp = 9999;
						}
						else{
							xp = ip%(Nx);
							yp = ((ip - xp)/Nx)%Nx;
							zp = ((ip - xp - (Nx*yp))/(Nx*Nx))%(Nx);
						}
						int tp = vp[1];
						count++;
						fout << "Forward step in direction "<< d << std::endl;
						fout << "{x,y,z,t} = " << "{" << xp << "," << yp << "," << zp << "," << tp << "}" << std::endl;
						std::vector<double> phi_1 = phi_a(Lattice, ip, tp, 1, true);
						std::vector<double> phi_2 = phi_a(Lattice, ip, tp, 2, true);
						fout << "phi_1 = " << phi_1[0] << "+ i*"<< phi_1[1] << "phi_2 = " << phi_2[0] << "+ i*"<< phi_2[1] << std::endl;
						std::vector<int> vn = negative_step(dim,d,i,t,Nx,Nt);
						int in = vn[0];
						int xn, yn, zn;
						if (in == 9999){
							xn = 9999;
							yn = 9999;
							zn = 9999;
						}
						else{
							xn = in%(Nx);
							yn = ((in - xn)/Nx)%Nx;
							zn = ((in - xn - (Nx*yn))/(Nx*Nx))%(Nx);
						}
						int tn = vn[1];
						fout << "Backward step in direction "<< d << std::endl;
						fout << "{x,y,z,t} = " << "{" << xn << "," << yn << "," << zn << "," << tn << "}" << std::endl;
						phi_1 = phi_a(Lattice, in, tn, 1, true);
						phi_2 = phi_a(Lattice, in, tn, 2, true);
						fout << "phi_1 = " << phi_1[0] << "+ i*"<< phi_1[1] << "phi_2 = " << phi_2[0] << "+ i*"<< phi_2[1] << std::endl;
						count++;
					}
					fout << std::endl;
				}
				else if (dim == 2){
					int x = i%(Nx);
					int y = ((i - x)/Nx)%Nx;
					//initialize x, y, and z
					//print i,t and x, y, z, t
					fout << "{i,t} = " << "{" << i << "," << t << "}" << std::endl;
					fout << "{x,y,t} = " << "{" << x << "," << y << "," << t << "}" << std::endl;
					//std::vector<int> negative_step(int dim, int dir, int i, int t, int Nx, int Nt)
					for (int d = 1; d <= dim; d++){
						std::vector<int> vp = positive_step(dim,d,i,t,Nx,Nt);
						int ip = vp[0];
						int xp, yp;
						if (ip == 9999){
							xp = 9999;
							yp = 9999;
						}
						else{
							xp = ip%(Nx);
							yp = ((ip - xp)/Nx)%Nx;
						}
						int tp = vp[1];
						count++;
						fout << "Forward step in direction "<< d << std::endl;
						fout << "{x,y,z,t} = " << "{" << xp << "," << yp << "," << tp << "}" << std::endl;
						std::vector<int> vn = negative_step(dim,d,i,t,Nx,Nt);
						int in = vn[0];
						int xn, yn;
						if (in == 9999){
							xn = 9999;
							yn = 9999;
						}
						else{
							xn = in%(Nx);
							yn = ((in - xn)/Nx)%Nx;
						}
						int tn = vn[1];
						fout << "Backward step in direction "<< d << std::endl;
						fout << "{x,y,t} = " << "{" << xn << "," << yn << "," << tn << "}" << std::endl;
						count++;
					}
					std::vector<int> vp = positive_step(dim,4,i,t,Nx,Nt);
					int ip = vp[0];
					int xp = ip%(Nx);
					int yp = ((ip - xp)/Nx)%Nx;
					int tp = vp[1];
					count++;
					fout << "Forward step in direction "<< 4 << std::endl;
					fout << "{x,y,t} = " << "{" << xp << "," << yp <<  "," << tp << "}" << std::endl;
					std::vector<int> vn = negative_step(dim,4,i,t,Nx,Nt);
					int in = vn[0];
					int xn = in%(Nx);
					int yn = ((in - xn)/Nx)%Nx;
					int tn = vn[1];
					fout << "Backward step in direction "<< 4 << std::endl;
					fout << "{x,y,t} = " << "{" << xn << "," << yn <<  "," << tn << "}" << std::endl;
					count++;
					fout << std::endl;
				}
				else if (dim == 1){
					int x = i;
					//initialize x, y, and z
					//print i,t and x, y, z, t
					fout << "{i,t} = " << "{" << i << "," << t << "}" << std::endl;
					fout << "{x,t} = " << "{" << x << "," << t << "}" << std::endl;
					//std::vector<int> negative_step(int dim, int dir, int i, int t, int Nx, int Nt)
					for (int d = 1; d <= dim; d++){
						std::vector<int> vp = positive_step(dim,d,i,t,Nx,Nt);
						int ip = vp[0];
						int xp = ip;
						int tp = vp[1];
						count++;
						fout << "Forward step in direction "<< d << std::endl;
						fout << "{x,t} = " << "{" << xp << "," << tp << "}" << std::endl;
						std::vector<int> vn = negative_step(dim,d,i,t,Nx,Nt);
						int in = vn[0];
						int xn = in;
						int tn = vn[1];
						fout << "Backward step in direction "<< d << std::endl;
						fout << "{x,y,t} = " << "{" << xn << "," << tn << "}" << std::endl;
						count++;
					}
					std::vector<int> vp = positive_step(dim,4,i,t,Nx,Nt);
					int ip = vp[0];
					int xp = ip;
					int tp = vp[1];
					count++;
					fout << "Forward step in direction "<< 4 << std::endl;
					fout << "{x,t} = " << "{" << xp <<  "," << tp << "}" << std::endl;
					std::vector<int> vn = negative_step(dim,4,i,t,Nx,Nt);
					int in = vn[0];
					int xn = in;
					int tn = vn[1];
					fout << "Backward step in direction "<< 4 << std::endl;
					fout << "{x,t} = " << "{" << xn <<  "," << tn << "}" << std::endl;
					count++;
					fout << std::endl;
				}
				else{
					fout << "Dimensions = " << dim << std::endl;
				}
			}
		}
		fout << "Check: " << 2*dim + 2 << " x lattice spacetime volume = " << 2*(dim+1)*size*Nt << ", and count = " << count << std::endl;
		fout.close();
	}
	else{
		std::cerr << "Unable to open file " << std::endl;
        exit(10);
	}
}