#############################################################################
#
# Makefile v4.2 of NRRB
#
#
# Author:  Casey E Berger
# Date:    08/20/2019 (version 4.2)
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

USE_OMP ?= FALSE

USE_TEST_CONSTANT_RNG ?= FALSE
USE_TEST_SEED_RNG ?= FALSE
USE_TEST_UNIFORM_RNG ?= FALSE

ifeq ($(USE_TEST_CONSTANT_RNG), TRUE)
	DEFINES += -DTEST_CONSTANT_RNG
endif

ifeq ($(USE_TEST_SEED_RNG), TRUE)
	DEFINES += -DTEST_SEED_RNG
endif

ifeq ($(USE_TEST_UNIFORM_RNG), TRUE)
	DEFINES += -DTEST_UNIFORM_RNG
endif

CXX = g++

ifeq ($(USE_OMP), TRUE)
	CXXFLAGS = -Wall -O3 -g -std=c++0x -fopenmp
else
	CXXFLAGS = -Wall -g -std=c++0x --coverage
endif

CXXFLAGS += $(DEFINES)

TARGET = v4.2
#LFLAGS = -L/opt/ddt/lib/64 -Wl,--undefined=malloc -ldmalloc -Wl,--allow-multipledefinition -pgc++libs

# ****************************************************
# Targets needed to bring the executable up to date

v4.2: v4.2.o test.o lattice_init.o lattice_save.o Langevin_evolution.o Observables.o
	$(CXX) $(CXXFLAGS) -o v4.2 v4.2.o test.o lattice_init.o lattice_save.o Langevin_evolution.o Observables.o

v4.2.o: v4.2.cpp test.h lattice_init.h lattice_save.h Langevin_evolution.h Observables.h
	$(CXX) $(CXXFLAGS) -c v4.2.cpp

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