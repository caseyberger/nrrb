#ifndef _OBSERVABLES_H
#define _OBSERVABLES_H

#include "Langevin.H"
#include "ProfileContainer.H"

struct RadCirc
{
    int radius;
    std::string circ_rad_name;
    amrex::Real circulation;

    RadCirc(int circ_radius)
    {
        // Initialize struct data
        radius = circ_radius;
        circulation = 0.0;
        circ_rad_name = "radius=" + std::to_string(radius);
    }

    void init_file(const std::string& logfile, const amrex::Geometry& geom)
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

        // Create the circulation datasets
        using namespace ClassyHDF;

        File obsFile(logfile);
        Group gcirc = obsFile.get_group("Circulation");
        // Todo: add loop center as an attribute of the circulation group
        gcirc.create_dataset<amrex::Real>(circ_rad_name);
    }

    void set_circulation(const amrex::Real loop_circ)
    {
        circulation = loop_circ;
    }

    void write(const std::string& logfile)
    {
        using namespace ClassyHDF;

        // append the circulation to its dataset in the observables file
        File obsFile(logfile);
        Group gcirc = obsFile.get_group("Circulation");
        gcirc.append(Data<amrex::Real>(circ_rad_name, {circulation}));
    }
};

class Observables
{
private:
	std::string observable_log_file;
    amrex::Vector<RadCirc> circulation;
    ProfileContainer density_profile;

public:
    Observables(const amrex::Geometry& geom, const amrex::DistributionMapping& dm, const amrex::BoxArray& ba, const NRRBParameters& nrrb, const int& nsteps);
    void initialize_files(const amrex::Geometry& geom, const NRRBParameters& nrrb, const int& nsteps);
    void update(const int nL, const amrex::Real Ltime, const amrex::MultiFab& Lattice,
                const amrex::GeometryData& geom, const NRRBParameters& nrrb_parm);
};

#endif
