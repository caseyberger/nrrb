DIM ?= 3

TOP := $(CL_HOME)

EBASE := ComplexLangevin

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

ifeq ($(USE_TEST_CONSTANT_RNG), TRUE)
	DEFINES += -DTEST_CONSTANT_RNG
endif

ifeq ($(USE_TEST_SEED_RNG), TRUE)
	DEFINES += -DTEST_SEED_RNG
endif

ifeq ($(USE_TEST_UNIFORM_RNG), TRUE)
	DEFINES += -DTEST_UNIFORM_RNG
endif

Bdirs   := Source
Bpack   += $(foreach dir, $(Bdirs), $(TOP)/$(dir)/Make.package)
Blocs   += $(foreach dir, $(Bdirs), $(TOP)/$(dir))

CLASSYHDF_SRC = $(TOP)/Submodules/ClassyHDF/Source
Bpack += $(CLASSYHDF_SRC)/Make.package
Blocs += $(CLASSYHDF_SRC)

include $(Bpack)

INCLUDE_LOCATIONS += $(Blocs)
VPATH_LOCATIONS   += $(Blocs)

Pdirs             := Base Particle
Ppack             += $(foreach dir, $(Pdirs), $(AMREX_HOME)/Src/$(dir)/Make.package)

include $(Ppack)

all: $(executable)
	@echo SUCCESS!

include $(AMREX_HOME)/Tools/GNUMake/Make.rules
