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
HDF5_HOME =/usr/local/hdf5/parallel/f90_

ifndef HOST
   HOST = `hostname`
endif

#Define compiler flags

FC = f90

# Flags for Absoft f90
ifeq ($(FC),f90)
  OPTIMIZATION = -O3 -cpu:host
#  DEBUG        = -g -et -Rb -Rp -Rc
  DEBUG        = -g -en -et -trap=DIVBYZERO,INVALID
  PROFILE      = -P
  STATIC       = -s
  MODFLAG      = -p
  ADDMODFLAG   = -p
  VERBOSE      = -v
  RELEASE      = absoft

  ifdef BOPT
    PETSC_ARCH = linux_absoft_
    CPPFLAGS  += -Dabsoft_ 
    include ${PETSC_DIR}/bmake/common/base
    FC         = $(MPI_HOME)/bin/mpif90
    RELEASE      = absoft_
  else
    FFLAGS       = -w -YEXT_NAMES=LCS -YEXT_SFX=_ -YCFRL=1
    LDFLAGS      = -lU77 
    ifeq ($(HOST),nip.lanl.gov)
#      LIBS      += -llapack -lblas -lg2c 
      LIBS      += -llapack_f90_ -lblas_f90_
    else
      LIBS      += -llapack -lblas
    endif
  endif

  HDF5 = true
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  OPTIMIZATION = -O
#  DEBUG        = -g --chkglobal --warn --f95
  DEBUG        = -g --chk aseu --warn --f95
#  DEBUG        = -g --f95
  PROFILE      =
  STATIC       = 
  MODFLAG      = -M
  ADDMODFLAG   = -I
  VERBOSE      = --verbose
  RELEASE      = lahey

  ifdef BOPT
    CPPFLAGS  += -Dlahey 
    FFLAGS    += --ml cdecl
    PETSC_ARCH = linux_lahey
    include ${PETSC_DIR}/bmake/common/base
    FC         = $(MPI_HOME)/bin/mpif90
  else
    FFLAGS    += -X9
    LDFLAGS    =
    ifeq ($(HOST),nip.lanl.gov)
      LIBS      += -llapack -lblas -lg2c
    else
      LIBS      += -llapack -lblas
    endif
#    LIBS       = -llapackmt -lblasmt
  endif
endif

# Flags for Intel ifort
ifeq ($(FC),ifort)
  OPTIMIZATION = -O2 -mp -axW
  DEBUG = -g -check all -traceback
#  DEBUG        = -g -check -traceback
  PROFILE      = -p
  STATIC       =
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  RELEASE      = intel

  ifdef BOPT
    CPPFLAGS  += -DNVAR=8
    FFLAGS    += -DNVAR=8
    PETSC_ARCH = linux_intel
    include ${PETSC_DIR}/bmake/common/base
    FC         = $(MPI_HOME)/bin/mpif90
  else
    FFLAGS      += -vec_report0 -w
    LDFLAGS      = 
    ifeq ($(HOST),nip.lanl.gov)
      LIBS      += -llapack -lblas -lg2c
    else
      LIBS      += -llapack -lblas
    endif
#    LIBS         = -llapack_intel -lblas_intel
  endif
endif

# Flags for Intel ifort
ifeq ($(FC),g95)
  OPTIMIZATION = -O2
  DEBUG = -g
#  DEBUG        = -g -check -traceback
  PROFILE      = -pg
  STATIC       =
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  RELEASE      = g95

  ifdef BOPT
    CPPFLAGS  += -DNVAR=8
    FFLAGS    += -DNVAR=8
    PETSC_ARCH = linux_intel
    include ${PETSC_DIR}/bmake/common/base
    FC         = $(MPI_HOME)/bin/mpif90
  else
    FFLAGS      += 
    ifeq ($(HOST),nip.lanl.gov)
      LIBS      += -llapack -lblas -lg2c
    else
      LIBS      += -llapack -lblas
    endif
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

#Define relevant directories

SUBDIRS = src plot

MODPATH = $(MODFLAG).

# HDF5 setup

ifdef HDF5
  LIBS     += -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz 
  CPPFLAGS += -Dhdf5 -I$(HDF5_HOME)/include
  MODPATH  += $(ADDMODFLAG)$(HDF5_HOME)/lib
endif

# PIXIE setup

ifdef VECPOT
  CPPFLAGS += -Dvec_pot
  ifdef VECPOT
    TARGET = code_a
  endif
endif

# VMEC setup

ifdef VMEC
  LIBS     += -L../contrib/vmec/lib -lstell
  CPPFLAGS += -Dvmec
  MODPATH  += $(ADDMODFLAG)../contrib/vmec/lib
endif

# Petsc setup

ifdef BOPT
  ifdef VECPOT
    TARGET = petsc_a
  else
    TARGET = petsc
  endif
  CPPFLAGS += -Dpetsc -DNVAR=8 -I$(PETSC_DIR)/include -I${PETSC_DIR}/bmake/$(PETSC_ARCH) -I$(MPI_HOME)/include
#  MODPATH  += $(ADDMODFLAG)$(MPI_HOME)/include
endif

# XDRAW setup

ifdef NOBC
  CPPFLAGS += -Dnobc
endif

#Export required variables

export FC FFLAGS CPPFLAGS MODFLAG ADDMODFLAG MODPATH LIBS LDFLAGS HDF5_HOME \
       BOPT PETSC_DIR PETSC_ARCH VECPOT VMEC

#Define targets

.PHONY: pixie3d pixplot distclean petsc all contrib contrib_clean \
        vmec vmec_clean setup serial-tests serial-tests-rebuild $(SUBDIRS)

all: contrib src plot

pixie3d: src

pixplot: plot

$(SUBDIRS):
	$(MAKE) -e -C $@ $(TARGET)

distclean: contrib_clean
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir distclean;  done

setup: ;
	-cd plot ; ln -s ../src/Makefile
	-for subdir in common plot tests/serial ; do \
		$(MAKE) -C $$subdir setup;  done

serial-tests: ;
	$(MAKE) -C tests/serial test

serial-tests-rebuild: ;
	$(MAKE) -C tests/serial rebuild

# CONTRIBUTED SOFTWARE

contrib: vmec

contrib_clean: vmec_clean

vmec:
ifdef VMEC
	$(MAKE) -e -C contrib/vmec/LIBSTELL release
endif

vmec_clean:
ifdef VMEC
	$(MAKE) -e -C contrib/vmec/LIBSTELL/Release -f makelibstell clean
endif
