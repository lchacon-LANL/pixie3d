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

PETSC_DIR =/usr/local/petsc-2.2.0
HDF5_HOME =/usr/local

# Petsc include

ifdef BOPT
  PETSCTARGET = petsc
  HDF5 = true
endif

#Define compiler flags

FC = f90

# Flags for Absoft f90
ifeq ($(FC),f90)
  OPTIMIZATION = -O2 -cpu:host
  DEBUG        = -g -et -Rb -Rp -Rc
#  DEBUG        = -g -en
  PROFILE      = -P
  STATIC       = -s
  MODFLAG      = -p
  ADDMODFLAG   = -p
  FFLAGS       = -w
  LIBS        += -llapack_f90 -lblas_f90 -lU77
  VERBOSE      = -v

  ifdef BOPT
    PETSC_ARCH = linux_absoft_
    CPPFLAGS  += -Dabsoft_ 
    include ${PETSC_DIR}/bmake/common/base
    FC         = $(MPI_HOME)/bin/mpif90
  endif
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  OPTIMIZATION = -O
#  DEBUG        = -g --chkglobal
  DEBUG        = -g --chk aseu --warn --f95
#  DEBUG        = -g --f95
  PROFILE      =
  STATIC       = 
  MODFLAG      = -M
  ADDMODFLAG   = -I
  FFLAGS      += -X9
#  LIBS         = -llapackmt -lblasmt
  LIBS         = -llapack -lblas -lf2c
  VERBOSE      = --verbose

  ifdef BOPT
    CPPFLAGS  += -Dlahey 
    FFLAGS    += --ml cdecl
    PETSC_ARCH = linux_lahey
    include ${PETSC_DIR}/bmake/common/base
    FC         = $(MPI_HOME)/bin/mpif90
  endif
endif

# Flags for Intel ifort
ifeq ($(FC),ifort)
  OPTIMIZATION = -O2 -mp -axW
  DEBUG = -g
#  DEBUG        = -g -check -traceback
  PROFILE      = -p
  STATIC       =
  MODFLAG      = -I
  ADDMODFLAG   = -I
  FFLAGS      += -vec_report0 -w
  LIBS         = -llapack -lblas -lf2c
#  LIBS         = -llapack_intel -lblas_intel
  VERBOSE      = -v

  ifdef BOPT
    CPPFLAGS  += -DNVAR=8
    FFLAGS    += -DNVAR=8
    PETSC_ARCH = linux_intel
    include ${PETSC_DIR}/bmake/common/base
    FC         = $(MPI_HOME)/bin/mpif90
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

MODPATH = $(MODFLAG).

# HDF5 setup

ifdef HDF5
  LIBS     += -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz 
  CPPFLAGS += -Dhdf5 -I$(HDF5_HOME)/include
  MODPATH  += $(ADDMODFLAG)$(HDF5_HOME)/lib
endif

# Petsc setup

ifdef BOPT
  CPPFLAGS += -Dpetsc -DNVAR=8
  MODPATH  += $(ADDMODFLAG)$(MPI_HOME)/include
endif

#Export required variables

export FC FFLAGS CPPFLAGS MODFLAG ADDMODFLAG MODPATH LIBS HDF5_HOME \
       BOPT PETSC_DIR PETSC_ARCH

#Define targets

.PHONY: pixie3d pixplot distclean petsc $(SUBDIRS)

pixie3d: src

pixplot: plot

$(SUBDIRS):
	$(MAKE) -e -C $@ $(PETSCTARGET)

distclean: 
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir distclean;  done

