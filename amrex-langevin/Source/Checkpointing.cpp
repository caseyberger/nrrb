#include <Langevin.H>

using namespace amrex;

namespace {
    // utility to skip to next line in Header
    void GotoNextLine (std::istream& is)
    {
        constexpr std::streamsize bl_ignore_max { 100000 };
        is.ignore(bl_ignore_max, '\n');
    }
}

void WriteCheckpointFile (const MultiFab& lattice, const int step, const Real time)
{
    /* Saves the lattice data, step, time, and timestep dt in a checkpoint file. */

    // chk00010            write a checkpoint file with this root directory
    // chk00010/Header     this contains information you need to save (e.g., finest_level, t_new, etc.) and also
    //                     the BoxArrays at each level
    // chk00010/Level_0/
    // chk00010/Level_1/
    // etc.                these subdirectories will hold the MultiFab data at each level of refinement

    // checkpoint file name, e.g., chk00010
    const std::string chk_file = "chk";
    const std::string& checkpointname = amrex::Concatenate(chk_file,step,7);

    amrex::Print() << "Writing checkpoint " << checkpointname << "\n";

    const int finest_level = 0; // we are currently single-level
    const int nlevels = finest_level+1;

    // ---- prebuild a hierarchy of directories
    // ---- dirName is built first.  if dirName exists, it is renamed.  then build
    // ---- dirName/subDirPrefix_0 .. dirName/subDirPrefix_nlevels-1
    // ---- if callBarrier is true, call ParallelDescriptor::Barrier()
    // ---- after all directories are built
    // ---- ParallelDescriptor::IOProcessor() creates the directories
    amrex::PreBuildDirectorHierarchy(checkpointname, "Level_", nlevels, true);

    // write Header file
    if (ParallelDescriptor::IOProcessor()) {

        std::string HeaderFileName(checkpointname + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;
        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        HeaderFile.open(HeaderFileName.c_str(), std::ofstream::out   |
                                                std::ofstream::trunc |
                                                std::ofstream::binary);
        if( ! HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);

        // write out title line
        HeaderFile << "Checkpoint file for AMReX-CL\n";

        // write out finest_level
        HeaderFile << finest_level << "\n";

        // write the number of components
        HeaderFile << lattice.nComp() << "\n";

        // write the number of ghost cells
        HeaderFile << lattice.nGrow() << "\n";

        // write out step number
        HeaderFile << step << " ";
        HeaderFile << "\n";

        // write out time
        HeaderFile << time << " ";
        HeaderFile << "\n";

        // write the BoxArray
        lattice.boxArray().writeOn(HeaderFile);
        HeaderFile << '\n';

        // write the number of OpenMP threads we are using
        HeaderFile << OpenMP::get_max_threads() << '\n';

        // write the random number generation state
        amrex::SaveRandomState(HeaderFile);
    }

    // write the MultiFab data to, e.g., chk00010/Level_0/
    const int lev = 0;
    VisMF::Write(lattice,
                 amrex::MultiFabFileFullPrefix(lev, checkpointname, "Level_", "Cell"));

}


void ReadCheckpointFile (const std::string restart_chkfile, MultiFab& lattice_old, MultiFab& lattice_new, int& step, Real& time)
{
    /* Defines lattice_old and lattice_new, and initializes lattice_new
       and the step, time, and timestep dt. */

    amrex::Print() << "Restart from checkpoint " << restart_chkfile << "\n";

    // Header
    std::string File(restart_chkfile + "/Header");

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);

    std::string line, word;
    int finest_level;

    int chk_ncomp, chk_nghost;

    // read in title line
    std::getline(is, line);

    // read in finest_level
    is >> finest_level;
    GotoNextLine(is);

    // read in number of components & assert they are the same as here
    is >> chk_ncomp;
    GotoNextLine(is);

    // read in number of ghost cells & assert they are the same as here
    is >> chk_nghost;
    GotoNextLine(is);

    // read in "new" step
    is >> word;
    step = std::stoi(word);

    // read in "new" time
    is >> word;
    time = std::stod(word);

    const int lev = 0;

    // read in level 'lev' BoxArray from Header
    BoxArray ba;
    ba.readFrom(is);
    GotoNextLine(is);

    // read the number of OpenMP threads from the previous run
    int previous_omp_nthreads = 1;
    is >> previous_omp_nthreads;
    GotoNextLine(is);

    // read the random number generation state
    amrex::RestoreRandomState(is, previous_omp_nthreads, step);

    // create a distribution mapping
    DistributionMapping dm { ba, ParallelDescriptor::NProcs() };

    // build MultiFab data
    lattice_old.define(ba, dm, chk_ncomp, chk_nghost);
    lattice_new.define(ba, dm, chk_ncomp, chk_nghost);

    // read in the MultiFab data to initialize "new" time lattice
    VisMF::Read(lattice_new,
                amrex::MultiFabFileFullPrefix(lev, restart_chkfile, "Level_", "Cell"));

}