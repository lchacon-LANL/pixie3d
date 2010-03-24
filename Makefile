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

# Defaults

FC  = gfortran
OPT = O 

PETSC_DIR ?=/usr/local/petsc-2.3.3
HDF5_HOME ?=/usr/local/hdf5/parallel/mpich2_f90_
HDF5_LIBS ?= -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5

PREPROC = -D

HDF5 = t

#FPA = t

LIBS = -llapack -lblas

#CPPFLAGS += -DRFX

# System-dependent variables

-include common/make/make.mach.inc

#Define compiler flags

-include common/make/make.comp.inc

#Define relevant directories

BINDIR = $(PWD)/bin

SUBDIRS = src plot

MODPATH = $(MODFLAG).

# PIXIE setup

REL1=0
REL2=9
CPPFLAGS += $(PREPROC)REL1=$(REL1) $(PREPROC)REL2=$(REL2)

ifdef VECPOT
  CPPFLAGS += $(PREPROC)vec_pot
  TARGET = code_a
endif

ifdef PER_BC_SYNC
  ifndef BOPT
    CPPFLAGS += $(PREPROC)PER_BC_SYNC
  endif
endif

# HDF5 setup

ifeq ($(HDF5),t)
  CONTRIBLIBS = $(HDF5_LIBS)
  CPPFLAGS   += $(PREPROC)hdf5 $(HDF5_INC)
  MODPATH    += $(ADDMODFLAG).$(HDF5_MOD)
endif

# VMEC setup

ifeq ($(VMEC),t)
  VMEC_DIR     = contrib/vmec/LIBSTELL
  ifeq ($(NETCDF),t)
    CPPFLAGS   += $(PREPROC)NETCDF 
  endif

  CONTRIBLIBS += -L$(PWD)/contrib/vmec/lib -lstell $(NETCDF_LIBS)
  CPPFLAGS    += $(PREPROC)vmec $(NETCDF_INC)
  MODPATH     += $(ADDMODFLAG)$(PWD)/contrib/vmec/lib
endif

# ARPACK setup

ifdef ARPACK
  CONTRIBLIBS += $(ARPACK_LIBS)
  CPPFLAGS    += $(PREPROC)arpack
endif

# FPA setup

ifdef FPA
  CONTRIBLIBS += -L$(PWD)/common/contrib/fpa/lib -lfpa
  CPPFLAGS    += $(PREPROC)FPA
  MODPATH     += $(ADDMODFLAG)$(PWD)/common/contrib/fpa/lib
endif

# LSODE setup

CONTRIBLIBS += -L$(PWD)/common/contrib/lsode -llsode

# SLATEC setup

CONTRIBLIBS += -L$(PWD)/common/contrib/slatec/lib -lslatec
MODPATH     += $(ADDMODFLAG)$(PWD)/common/contrib/slatec/lib

# PETSC setup

ifdef BOPT
  include ${PETSC_DIR}/conf/base

  ifdef VECPOT
    TARGET = petsc_a
  else
    TARGET = petsc
  endif

  CPPFLAGS += $(PREPROC)petsc $(PREPROC)NVAR=8 -I$(PETSC_DIR)/include -I${PETSC_DIR}/$(PETSC_ARCH)/include

  ifdef PETSC_C
    CPPFLAGS += $(PREPROC)petsc_c
    SNES_OPT = -snes_mf
  endif
endif

#SAMRAI setup

ifdef SAMR
   include ${SAMRAI}/config/Makefile.config

   PDIM = 3
   OBJECT=${SAMRAI}
   CXXFLAGS_EXTRA += -DNDIM=$(PDIM)

   TARGET = samrai
   CPPFLAGS += -Dsamrai

ifdef VECPOT
  CPPFLAGS += -Dvec_pot
endif

   CPPFLAGS_EXTRA += -I${AMRUTILITIES_HOME}/include
   CPPFLAGS_EXTRA += -I${SAMRSOLVERS_HOME}/include
   CXXFLAGS_EXTRA += -I${AMRUTILITIES_HOME}/include
   CXXFLAGS_EXTRA += -I${SAMRSOLVERS_HOME}/include
   LDFLAGS_EXTRA += -L${AMRUTILITIES_HOME}/lib   
   LDLIBS_EXTRA += ${AMRUTILITIES_HOME}/lib/libAMRUtils3d.a
   LDLIBS_EXTRA += ${SAMRSOLVERS_HOME}/lib/liblinearops3d.a
   LDLIBS_EXTRA += ${SAMRSOLVERS_HOME}/lib/libmlsolvers3d.a
   LDLIBS_EXTRA += ${SAMRSOLVERS_HOME}/lib/libpreconditionerbase3d.a
   LDFLAGS_EXTRA += -L${SAMRSOLVERS_HOME}/lib   
   LDLIBS_EXTRA += ${SAMRSOLVERS_HOME}/lib/libtimeintegrators3d.a
   #LDLIBS_EXTRA += -lAMRUtils${PDIM}d
endif

#Export required variables

export FC FFLAGS CPPFLAGS MODFLAG ADDMODFLAG MODPATH LIBS LDFLAGS HDF5_HOME \
       H5LIBS MPI_HOME BOPT PETSC_DIR PETSC_ARCH VECPOT VMEC ARPACK FPA SNES_OPT \
       BINDIR CONTRIBLIBS MPIEXEC FLINKER PETSC_SNES_LIB SAMR CPPFLAGS_EXTRA \
       CXXFLAGS_EXTRA LDFLAGS_EXTRA LDLIBS_EXTRA LIBSAMRAI3D LIBSAMRAI

#Define targets

.PHONY: pixie3d pixplot distclean petsc all contrib contrib_clean \
        vmec vmec_clean setup serial-tests-b serial-tests-a \
        rebuild-serial-tests-a rebuild-serial-tests-b \
        parallel-tests-b parallel-tests-a \
        rebuild-parallel-tests-a rebuild-parallel-tests-b $(SUBDIRS)

all: contrib src plot

pixie3d: contrib src

pixplot: contrib plot

$(SUBDIRS):
	$(MAKE) -e -C $@ $(TARGET)

# SETUP

setup:
	-cd plot ; ln -s ../src/Makefile
	-tar xzf contrib.tgz
	-for subdir in common plot tests/serial tests/parallel ; do \
		$(MAKE) -C $$subdir setup;  done

# SERIAL TESTS

all-serial-tests: serial-tests-a serial-tests-b

rebuild-all-serial-tests: rebuild-serial-tests-a rebuild-serial-tests-b

serial-tests-a: ;
	$(MAKE) -e -C tests/serial test-a

serial-tests-b: ;
	$(MAKE) -e -C tests/serial test-b

rebuild-serial-tests-a: ;
	$(MAKE) -e -C tests/serial rebuild-a

rebuild-serial-tests-b: ;
	$(MAKE) -e -C tests/serial rebuild-b

# PARALLEL TESTS

all-parallel-tests: parallel-tests-a parallel-tests-b

rebuild-all-parallel-tests: rebuild-parallel-tests-a rebuild-parallel-tests-b

parallel-tests-a: ;
ifdef BOPT
	$(MAKE) -e -C tests/parallel test-a
else
	-@echo "Please, specify BOPT in make command"
endif

parallel-tests-b: ;
ifdef BOPT
	$(MAKE) -e -C tests/parallel test-b
else
	-@echo "Please, specify BOPT in make command"
endif

rebuild-parallel-tests-a: ;
ifdef BOPT
	$(MAKE) -e -C tests/parallel rebuild-a
else
	-@echo "Please, specify BOPT in make command"
endif

rebuild-parallel-tests-b: ;
ifdef BOPT
	$(MAKE) -e -C tests/parallel rebuild-b
else
	-@echo "Please, specify BOPT in make command"
endif

# CONTRIBUTED LIBRARIES

contrib:
	$(MAKE) -e -C common contrib
ifeq ($(VMEC),t)
	$(MAKE) -e -C $(VMEC_DIR) release INC_PATH=$(NETCDF_INC)
endif

contrib_clean:
	$(MAKE) -e -C common contrib_clean
ifeq ($(VMEC),t)
	$(MAKE) -e -C $(VMEC_DIR)/Release -f makelibstell clean
endif

# CLEAN ALL

allclean: contrib_clean distclean

distclean:
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir distclean;  done
