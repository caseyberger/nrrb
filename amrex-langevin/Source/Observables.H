#ifndef _OBSERVABLES_H
#define _OBSERVABLES_H

#include "Langevin.H"

namespace Obs {
enum ObservableIndices {PhiSqRe = 0, PhiSqIm, DensRe, DensIm,
                        LzRe, LzIm, SRe, SIm, Circ1, Circ2,
                        NumObservables};
};

struct RadCirc
{
    int radius;
    amrex::Real circulation;
    std::string logfile;

    RadCirc(int circ_radius, std::string logfile_suffix)
    {
        // Initialize struct data
        radius = circ_radius;
        circulation = 0.0;

        logfile = "Circ_loop_" + std::to_string(radius) + "_" + logfile_suffix;
    }

    void init_file(const amrex::Geometry& geom)
    {
        // Get domain center location
        const auto domain_box = geom.Domain();
        const auto domain_lo = amrex::lbound(domain_box);
        const auto domain_hi = amrex::ubound(domain_box);

        // For this to work properly, the following must be true:
        // - the domain sizes in x and y must be odd
        // - the center of the domain must satisfy the following:
        const int x_center = 0.5 * (domain_lo.x + domain_hi.x);
        const int y_center = 0.5 * (domain_lo.y + domain_hi.y);

        // Create the circulation logfile
        std::ofstream obsFile;
        obsFile.open(logfile, std::fstream::trunc);

        // Print loop center and radius to logfile
        obsFile << "center = (" << x_center << ","<< y_center << ") and loop radius = " << radius;
        obsFile << std::endl;
        obsFile.close();
    }

    void set_circulation(const amrex::Real loop_circ)
    {
        circulation = loop_circ;
    }

    void write()
    {
        std::ofstream obsFile;
        obsFile.open(logfile, std::fstream::app);
        obsFile << circulation << ",\n";
        obsFile.close();
    }
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

#endif
