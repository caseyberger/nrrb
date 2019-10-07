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
#include "Langevin.H"
#include <omp.h>
/*
Defines the function that evolves the lattice in Langevin time, as well as the real and imaginary
drift functions, K_a.

Includes definition of an antisymmetric tensor, epsilon, and functions that modify the location
on the lattice to one positive or negative step in the direction specified.
*/

/* enum ComplexComponents { Real = 0, Imag }; */

//COMPUTING THE REAL DRIFT FUNCTION
Real K_a_Re(Real m, Real l,Real w, Real w_t, int a, Real dtau, Real mu, amrex::Array4<amrex::Real> const& Lattice_old, const amrex::GeometryData& geom, int i, int j, int t){
    const auto domain_xlo = geom.ProbLo();
    const auto domain_xhi = geom.ProbHi();
    const auto domain_box = geom.Domain();
    const auto domain_lo = amrex::lbound(domain_box);
    const auto domain_hi = amrex::ubound(domain_box);
    const Real x_center = 0.5 * (domain_xlo[0] + domain_xhi[0]);
    const Real y_center = 0.5 * (domain_xlo[1] + domain_xhi[1]);
    const Real dx_cell = geom.CellSize(0);
    const Real dy_cell = geom.CellSize(1);
    const Real x_cell = domain_xlo[0] + dx_cell * (i + 0.5 - domain_lo[0]);
    const Real y_cell = domain_xlo[1] + dy_cell * (j + 0.5 - domain_lo[1]);
    const Real r2_cell = (x_cell - x_center)**2 + (y_cell - y_center)**2;

	Real Ka = 0.0;
	//doing -K_{a}^{R} and then returning -Ka at the end for simplicity
	//modify mu, m, w, wtr, and l by dtau
	mu = dtau*mu;
	m = m/dtau;
	w_t = dtau*w_t;
	w = dtau*w;
	l = dtau*l;

	//CHEMICAL POTENTIAL + TIME DERIVATIVE PART OF Ka
	//\phi_{a,r}^{R} - e^{\mu}/2 (\phi_{a,r-t}^{R} + \phi_{a,r+t}^{R})
	//no special conditions here because we have periodic boundary conditions in time
	Ka += Lattice_old(i,j,t,field_comp(a,true));
	Ka += - 0.5 * exp(mu) * Lattice_old(i,j,t-1,field_comp(a,true));
	Ka += - 0.5 * exp(mu) * Lattice_old(i,j,t+1,field_comp(a,true));

	for (int b=1; b<=2; b++){
		//\eps_{ab} e^{\mu}/2 (\phi_{b,r-t}^{I}-\phi_{b,r+t}^{I})
		Ka += epsilon(a,b) * 0.5 * exp(mu) * Lattice_old(i,j,t-1,field_comp(b,false));
		Ka -= epsilon(a,b) * 0.5 * exp(mu) * Lattice_old(i,j,t+1,field_comp(b,false));
	}//last checked on 10.07.19

	//SPATIAL DERIVATIVE PART OF Ka
	//To test 0+1 dimensional case, just skip this loop
	//(d/m)\phi_{a}^{R}
	Ka += 1.*dim*Lattice(i,j,t,field_comp(a,true))/m;
	for (int d = 1; d <= AMREX_SPACEDIM-1; d++) { //loop over adjacent spatial sites!
		//if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
		// -(1/2m)\phi_{a,i+1}^{R} - (1/2m)\phi_{a,ii1}^{R}
        if (d == 1)
        {
            if (!(i == domain_lo[0] || i == domain_hi[0]))
            {
				Ka += -0.5 * Lattice_old(i+1,j,t,field_comp(a,true)) / m;
				Ka += -0.5 * Lattice_old(i-1,j,t,field_comp(a,true)) / m;
            }
        } else if (d == 2)
        {
            if (!(j == domain_lo[1] || j == domain_hi[1]))
            {
				Ka += -0.5 * Lattice_old(i,j+1,t,field_comp(a,true)) / m;
				Ka += -0.5 * Lattice_old(i,j-1,t,field_comp(a,true)) / m;
            }
        }
	}//updated on 10.07.19 - part mixing a and b removed per notes

	//ROTATIONAL AND TRAPPING PARTS OF Ka - ONLY WORKS IN 2D AND IF OMEGA > 0.
	//TRAPPING PIECE

	// 0.25*m*w_t^{2}*r2*(\phi_{a,r+t}^{R}+\phi_{a,r-t}^{R})
    Ka += 0.25*m*w_t*w_t*r2_cell * Lattice_old(i,j,t+1,field_comp(a,true));
    Ka += 0.25*m*w_t*w_t*r2_cell * Lattice_old(i,j,t-1,field_comp(a,true));

    for (int b=1; b<=2; b++){
        //TRAPPING PIECE
        // -0.25*m*w_t^{2} r2*\eps_{ab}(\phi_{b,r+t}^{I}+\phi_{b,r-t}^{I})
        Ka += -0.25*m*w_t*w_t*r2_cell*epsilon(a,b) * Lattice_old(i,j,t+1,field_comp(b,false));
        Ka += -0.25*m*w_t*w_t*r2_cell*epsilon(a,b) * Lattice_old(i,j,t-1,field_comp(b,false));
    }

    // ROTATIONAL PIECE
#if (AMREX_SPACEDIM == 3)
    const Real x = x_cell - x_center;
    const Real y = y_cell - y_center;

	// w (x-y)(\phi_{a,r+t}^{I}+\phi_{a,r-t}^{I})
	Ka += 0.5*w*x * Lattice_old(i,j,t+1,field_comp(a, false));
	Ka += -0.5*w*y * Lattice_old(i,j,t+1,field_comp(a, false));
	Ka += 0.5*w*x * Lattice_old(i,j,t-1,field_comp(a, false));
	Ka += -0.5*w*y * Lattice_old(i,j,t-1,field_comp(a, false));

	// (w/2) y (\phi_{a,r+x+t}^{I} + \phi_{a,r-x-t}^{I}) - (w/2) x (\phi_{a,r+y+t}^{I} + \phi_{a,r-y-t}^{I})
    //if r is on the border, this term is zero because it's a derivative of \phi_{a,r}\phi_{a,r+i}^{R}
    if (!(i == domain_lo[0] || i == domain_hi[0]))
    {
	    Ka += 0.5*w*y * Lattice_old(i+1,j,t+1,field_comp(a,false));
	    Ka += 0.5*w*y * Lattice_old(i-1,j,t-1,field_comp(a,false));
    }

    if (!(j == domain_lo[1] || j == domain_hi[1]))
    {
	    Ka += -0.5*w*x * Lattice_old(i,j+1,t+1,field_comp(a,false));
	    Ka += -0.5*w*x * Lattice_old(i,j-1,t-1,field_comp(a,false));
    }

    for (int b=1; b<=2; b++){
        //\eps_{ab} (w/2) (y (\phi_{b,r-x-t}^{R} + \phi_{b,r+x+t}^{R})- x (\phi_{b,r-y-t}^{R}+\phi_{b,r+y+t}^{R}))
        Ka += -0.5*epsilon(a,b)*w*x * Lattice_old(i,j+1,t+1,field_comp(b,true));
        Ka += - 0.5*epsilon(a,b)*w*x * Lattice_old(i,j-1,t-1,field_comp(b,true));

        Ka += 0.5*epsilon(a,b)*w*y * Lattice_old(i-1,j,t-1,field_comp(b,true));
        Ka += 0.5*epsilon(a,b)*w*y * Lattice_old(i+1,j,t+1,field_comp(b,true));

        //\eps_{ab} (w/2) (x-y)(\phi_{b,r+t}^{R}+\phi_{b,r-t}^{R})
        Ka += 0.5*epsilon(a,b)*w*x * Lattice_old(i,j,t+1,field_comp(b,true));
        Ka += 0.5*epsilon(a,b)*w*x * Lattice_old(i,j,t-1,field_comp(b,true));

        Ka += -0.5*epsilon(a,b)*w*y * Lattice_old(i,j,t+1,field_comp(b,true));
        Ka += -0.5*epsilon(a,b)*w*y * Lattice_old(i,j,t-1,field_comp(b,true));
    }
}
//last updated 10.07.19
#endif

	//INTERACTION PART OF Ka
	for (int b=1; b<=2; b++){
		Ka += l*Lattice(i,j,t,field_comp(b,true))*Lattice(i,j,t-1,field_comp(a,true))*Lattice(i,j,t-1,field_comp(b,true));
		Ka += - l*Lattice(i,j,t,field_comp(b,true))*Lattice(i,j,t-1,field_comp(a,false))*Lattice(i,j,t-1,field_comp(b,false));

		Ka += l*Lattice(i,j,t,field_comp(b,true))*Lattice(i,j,t+1,field_comp(a,true))*Lattice(i,j,t+1,field_comp(b,true));
		Ka += - l*Lattice(i,j,t,field_comp(b,true))*Lattice(i,j,t+1,field_comp(a,false))*Lattice(i,j,t+1,field_comp(b,false));

		Ka += -1.0*l*Lattice(i,j,t,field_comp(b,false))*Lattice(i,j,t-1,field_comp(a,true))*Lattice(i,j,t-1,field_comp(b,false));
		Ka += - l*Lattice(i,j,t,field_comp(b,false))*Lattice(i,j,t-1,field_comp(a,false))*Lattice(i,j,t-1,field_comp(b,true));

		Ka += -1.0*l*Lattice(i,j,t,field_comp(b,false))*Lattice(i,j,t+1,field_comp(a,true))*Lattice(i,j,t+1,field_comp(b,false));
		Ka += - l*Lattice(i,j,t,field_comp(b,false))*Lattice(i,j,t+1,field_comp(a,false))*Lattice(i,j,t+1,field_comp(b,true));

		Ka += 0.5*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t-1,field_comp(b,false))*Lattice(i,j,t-1,field_comp(b,false));
		Ka += -0.5*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t-1,field_comp(b,true))*Lattice(i,j,t-1,field_comp(b,true));

		Ka += 0.5*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t+1,field_comp(b,false))*Lattice(i,j,t+1,field_comp(b,false));
		Ka += -0.5*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t+1,field_comp(b,true))*Lattice(i,j,t+1,field_comp(b,true));

		Ka += l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t-1,field_comp(b,false))*Lattice(i,j,t-1,field_comp(b,true));
		Ka += l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t+1,field_comp(b,false))*Lattice(i,j,t+1,field_comp(b,true));

		Ka += -1.*epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t-1,field_comp(a,true))*Lattice(i,j,t-1,field_comp(b,false));
		Ka += -epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t-1,field_comp(a,false))*Lattice(i,j,t-1,field_comp(b,true));

		Ka += epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t+1,field_comp(a,true))*Lattice(i,j,t+1,field_comp(b,false));
		Ka += epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t+1,field_comp(a,false))*Lattice(i,j,t+1,field_comp(b,true));

		Ka += epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t-1,field_comp(a,false))*Lattice(i,j,t-1,field_comp(b,false));
		Ka += -epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t-1,field_comp(a,true))*Lattice(i,j,t-1,field_comp(b,true));

		Ka += epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t+1,field_comp(a,true))*Lattice(i,j,t+1,field_comp(b,true));
		Ka += -epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t+1,field_comp(a,false))*Lattice(i,j,t+1,field_comp(b,false));

		Ka += 0.5*epsilon(a,b)*l*Lattice(i,j,t,field_comp(b,false))*Lattice(i,j,t-1,field_comp(a,true))*Lattice(i,j,t-1,field_comp(a,true));
		Ka += -0.5*epsilon(a,b)*l*Lattice(i,j,t,field_comp(b,false))*Lattice(i,j,t-1,field_comp(a,false))*Lattice(i,j,t-1,field_comp(a,false));

		Ka += 0.5*epsilon(a,b)*l*Lattice(i,j,t,field_comp(b,false))*Lattice(i,j,t+1,field_comp(a,false))*Lattice(i,j,t+1,field_comp(a,false));
		Ka += -0.5*epsilon(a,b)*l*Lattice(i,j,t,field_comp(b,false))*Lattice(i,j,t+1,field_comp(a,true))*Lattice(i,j,t+1,field_comp(a,true));

		Ka += epsilon(a,b)*l*Lattice(i,j,t,field_comp(b,true))*Lattice(i,j,t-1,field_comp(a,false))*Lattice(i,j,t-1,field_comp(a,true));
		Ka += -epsilon(a,b)*l*Lattice(i,j,t,field_comp(b,true))*Lattice(i,j,t+1,field_comp(a,false))*Lattice(i,j,t+1,field_comp(a,true));

		Ka += 0.5*epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t-1,field_comp(a,false))*Lattice(i,j,t-1,field_comp(a,false));
		Ka += -0.5*epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t-1,field_comp(a,true))*Lattice(i,j,t-1,field_comp(a,true));

		Ka += 0.5*epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t+1,field_comp(a,true))*Lattice(i,j,t+1,field_comp(a,true));
		Ka += -0.5*epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,false))*Lattice(i,j,t+1,field_comp(a,false))*Lattice(i,j,t+1,field_comp(a,false));

		Ka += epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t+1,field_comp(a,false))*Lattice(i,j,t+1,field_comp(a,true));
		Ka += -epsilon(a,b)*l*Lattice(i,j,t,field_comp(a,true))*Lattice(i,j,t-1,field_comp(a,false))*Lattice(i,j,t-1,field_comp(a,true));
	}//last update 10.07.19 to make a gauge interaction
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
