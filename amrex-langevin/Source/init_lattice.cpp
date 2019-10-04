
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

void init_lattice (MultiFab lattice_mf)
{

    for ( MFIter mfi(lattice_mf); mfi.isValid(); ++mfi )
    {

      const Box& bx = mfi.validbox();

// void lattice_init(double *** Lattice, int size, int time_size)
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_real_distribution<double> distribution(-1.0,1.0);
        //assign each field a random value, uniformly distributed between -1 and 1

        for (int i = 0; i<size;i++){

                int Nx = sqrt(size);
                int x = i%Nx;
                int y = ((i - x)/Nx)%Nx;
                std::cout << "(" << x << "," << y << ") ";

                for (int j = 0; j<time_size;j++){

                        for (int k = 0; k<4;k++){

                                double r = distribution(mt);
                                //double r = 1.;
//                              Lattice[i][j][k] = r; //time-evolved fields - for future reference
//                              Lattice[i][j][k+4] = r; //old fields - for future reference
                        }
                }
        }
    }
    std::cout << "Fields initialized" << std::endl;
}
