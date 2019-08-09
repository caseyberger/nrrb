#############################################################################
#
# Makefile v4.1 of NRRB
#
#
# Author:  Casey E Berger
# Date:    07/30/2019 (version 4.1)
#
# Description:
# ------------
# GNU make is expected to use the Makefile. Other versions of makes
# may or may not work.
#
#
# For use with a profiler, compile with
# make-profiler-libraries --lib-type=static
#
#
#===========================================================================

CXX = g++
#CXXFLAGS = -Wall -O3 -g -std=c++0x -fopenmp
CXXFLAGS = -Wall -g -std=c++0x --coverage
TARGET = v4.1
#LFLAGS = -L/opt/ddt/lib/64 -Wl,--undefined=malloc -ldmalloc -Wl,--allow-multipledefinition -pgc++libs

# ****************************************************
# Targets needed to bring the executable up to date

v4.1: v4.1.o test.o lattice_init.o lattice_save.o Langevin_evolution.o Observables.o
	$(CXX) $(CXXFLAGS) -o v4.1 v4.1.o test.o lattice_init.o lattice_save.o Langevin_evolution.o Observables.o

v4.1.o: v4.1.cpp test.h lattice_init.h lattice_save.h Langevin_evolution.h Observables.h
	$(CXX) $(CXXFLAGS) -c v4.1.cpp

# Dependencies
test.o: test.cpp test.h
	$(CXX) $(CXXFLAGS) -c test.cpp

lattice_init.o: lattice_init.cpp lattice_init.h
	$(CXX) $(CXXFLAGS) -c lattice_init.cpp

lattice_save.o: lattice_save.cpp lattice_save.h
	$(CXX) $(CXXFLAGS) -c lattice_save.cpp

Langevin_evolution.o: Langevin_evolution.cpp Langevin_evolution.h
	$(CXX) $(CXXFLAGS) -c Langevin_evolution.cpp

Observables.o: Observables.cpp Observables.h Langevin_evolution.h
	$(CXX) $(CXXFLAGS) -c Observables.cpp

clean:
	rm *.o