# This makefile makes fortran code in present directory and subdirectories 
# specified in variable $(SUBDIRS). Command line options include, aside from
# standard make options, the variable $(OPT), which specifies level of 
# optimization:
#     + OPT = O (default) optimizes code.
#     + OPT = g creates debugging executable
#     + OPT = p creates profiling executable
# and the fortran compiler variable $(FC):
#     + FC = f90 (default) uses Absoft f90 
#     + FC = lf95 uses Lahey lf95.
# When PETSc is employed, the optimization level should be specified by
# the variable $(BOPT).
#
# A call example that uses lf95 in debugging mode is:
#
#        make FC=lf95 OPT=g
#
# A call example that employs PETSc is:
#
#        make BOPT=g petsc

# System-dependent variables

PETSC_DIR=/usr/local/petsc-2.2.0
PETSC_ARCH=linux_absoft_
HDF5_HOME=/usr/local

# Petsc include

ifneq ($(strip $(BOPT)),)
  PETSCTRUE = true
  include ${PETSC_DIR}/bmake/common/base
endif

#Define compiler flags

FC = f90

# Flags for Absoft f90
ifeq ($(FC),f90)
  OPTIMIZATION = -O2 -cpu:host
#  DEBUG        = -g -et -Rb -Rp
  DEBUG        = -g -ggdb
  PROFILE      = -P
  STATIC       = -s
  MODPATH      = -p
  ADDMODPATH   = -p
  FFLAGS       = -w
  LIBS         = -llapack_f90 -lblas_f90 -lU77 -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz
  VERBOSE      = -v
  GMODPATH  = -p$(HDF5_HOME)/lib

  ifdef PETSCTRUE
    CPPFLAGS += -Dabsoft_ -Dpetsc -DNVAR=8 -I$(HDF5_HOME)/include
    LIBS      = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz \
                -I$(MPI_HOME)/include -L$(MPI_HOME)/lib -lmpichf90 -lmpich -lpmpich 
    FC        = $(MPI_HOME)/bin/mpif90
  endif
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  OPTIMIZATION = -O
#  DEBUG        = -g --chkglobal
  DEBUG        = -g --chk ase
  PROFILE      =
  STATIC       = 
  MODPATH      = -M
  ADDMODPATH   = -I
  FFLAGS      +=
  LIBS         = -llapackmt -lblasmt
#  LIBS         = -llapack -lblas -lf2c
  VERBOSE      = --verbose

  ifdef PETSCTRUE
    CPPFLAGS  += -DNVAR=8 -Dlahey 
    FFLAGS    += -DNVAR=8 --ml cdecl
    PETSC_ARCH = linux_lahey
    FC        = $(MPI_HOME)/bin/mpif90
  endif
endif

# Flags for Intel ifort
ifeq ($(FC),ifort)
  OPTIMIZATION = -O2 -mp -axW
#  DEBUG = -g
  DEBUG        = -g -check -traceback
  PROFILE      = -p
  STATIC       =
  MODPATH      = -I
  ADDMODPATH   = -I
  FFLAGS      += -vec_report0
  LIBS         = -llapack -lblas -lf2c
#  LIBS         = -llapack_intel -lblas_intel
  VERBOSE      = -v

  ifdef PETSCTRUE
    CPPFLAGS  += -DNVAR=8
    FFLAGS    += -DNVAR=8
    PETSC_ARCH = linux_intel
    FC        = $(MPI_HOME)/bin/mpif90
  endif
endif

#Assemble compiler options

OPT = O 

ifneq (,$(findstring O,$(OPT)))
FFLAGS += $(OPTIMIZATION)
endif
ifneq (,$(findstring g,$(OPT)))
FFLAGS += $(DEBUG)
endif
ifneq (,$(findstring p,$(OPT)))
FFLAGS += $(PROFILE)
endif
ifneq (,$(findstring s,$(OPT)))
FFLAGS += $(STATIC)
endif
ifneq (,$(findstring v,$(OPT)))
FFLAGS += $(VERBOSE)
endif

#Define linker flags

LDFLAGS = 

#Define relevant directories

SUBDIRS = src plot

#Export required variables

export FC FFLAGS CPPFLAGS MODPATH ADDMODPATH GMODPATH LIBS HDF5_HOME \
       PETSC_DIR PETSC_ARCH PETSCTRUE

#Define targets

.PHONY: pixie3d pixplot distclean $(SUBDIRS)

pixie3d: src

pixplot: plot

$(SUBDIRS):
	$(MAKE) -e -C $@ 

petsc:
	$(MAKE) -e -C src petsc

distclean: 
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir distclean;  done

