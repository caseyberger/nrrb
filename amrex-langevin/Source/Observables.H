#include "Langevin.H"

namespace Obs {
enum ObservableIndices {PhiSqRe = 0, PhiSqIm, DensRe, DensIm,
                        LzRe, LzIm, SRe, SIm, Theta1, Theta2,
                        NumObservables};
};

namespace {
struct RadCirc
{
    amrex::Real radius;
    amrex::Real circulation;
    std::string logfile;

    RadCirc(amrex::Real circ_radius, std::string logfile_suffix)
    {
        // Initialize struct data
        radius = circ_radius;
        circulation = 0.0;
        logfile = "Circ_loop_" + std::to_string(radius) + "_" + logfile_suffix;
    }

    void init_file(const amrex::Geometry& geom)
    {
        // Get domain center location
        const auto domain_xlo = geom.ProbLo();
        const auto domain_xhi = geom.ProbHi();
        const Real x_center = 0.5 * (domain_xlo[0] + domain_xhi[0]);
        const Real y_center = 0.5 * (domain_xlo[1] + domain_xhi[1]);

        // Create the circulation logfile
        std::ofstream obsFile;
        obsFile.open(logfile, std::fstream::trunc);

        // Print loop center and radius to logfile
        obsFile << "center = (" << x_center << ","<< y_center << ") and loop radius = " << radius;
        obsFile << std::endl;
        obsFile.close();
    }

    void write()
    {
        std::ofstream obsFile;
        obsFile.open(logfile, std::fstream::app);
        obsFile << circulation << ",\n";
        obsFile.close();
    }
};
};

class Observables
{
private:
	std::string logfile_suffix;
	std::string observable_log_file;
    amrex::Vector<RadCirc> circulation;

public:
    Observables(const amrex::Geometry& geom, const NRRBParameters& nrrb, const int& nsteps);
    void initialize_files(const amrex::Geometry& geom);
    void update(const int nL, const amrex::MultiFab& Lattice, 
                const amrex::GeometryData& geom, const NRRBParameters& nrrb_parm);
};