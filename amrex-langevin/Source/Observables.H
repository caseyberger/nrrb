#include "Langevin.H"

namespace Obs {
enum ObservableIndices {PhiSqRe = 0, PhiSqIm, DensRe, DensIm,
                        LzRe, LzIm, SRe, SIm, NumObservables};
};
class Observables
{
private:
	std::string logfile_suffix;
	std::string observable_log_file;
    std::map<int, std::string> circulation_radii_logs;

public:
    Observables(const amrex::GeometryData& geom, const NRRBParameters& nrrb, const int& nsteps);
    void initialize_files(const amrex::GeometryData& geom);
    void update(const int nL, const amrex::MultiFab& Lattice, 
                const amrex::GeometryData& geom, const NRRBParameters& nrrb_parm);
};