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

FC  = f90
OPT = O 

PETSC_DIR ?=/usr/local/petsc-2.3.3
HDF5_HOME ?=/usr/local/hdf5/parallel/mpich2_f90_
HDF5_LIBS ?= -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5

PREPROC = -D

HDF5 = t

LIBS = -llapack -lblas

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

# PETSC setup

ifdef BOPT
  include ${PETSC_DIR}/bmake/common/base
endif

# HDF5 setup

ifeq ($(HDF5),t)
  CONTRIBLIBS = $(HDF5_LIBS)
  CPPFLAGS   += $(PREPROC)hdf5 -I$(HDF5_HOME)/include
  MODPATH    += $(ADDMODFLAG)$(HDF5_HOME)/lib
endif

# VMEC setup

ifdef VMEC
  ifdef NETCDF
    CPPFLAGS   += $(PREPROC)NETCDF 
  endif

  CONTRIBLIBS += -L../contrib/vmec/lib -lstell $(NETCDF_LIBS)
  CPPFLAGS    += $(PREPROC)vmec $(NETCDF_INC)
  MODPATH     += $(ADDMODFLAG)../contrib/vmec/lib
endif

# ARPACK setup

ifdef ARPACK
  CONTRIBLIBS += $(ARPACK_LIBS)
  CPPFLAGS    += $(PREPROC)arpack
endif

# Petsc setup

ifdef BOPT
  ifdef VECPOT
    TARGET = petsc_a
  else
    TARGET = petsc
  endif

  CPPFLAGS += $(PREPROC)petsc $(PREPROC)NVAR=8 -I$(PETSC_DIR)/include -I${PETSC_DIR}/bmake/$(PETSC_ARCH)

  ifdef PETSC_C
    CPPFLAGS += $(PREPROC)petsc_c
    SNES_OPT = -snes_mf
  endif
endif

#Export required variables

export FC FFLAGS CPPFLAGS MODFLAG ADDMODFLAG MODPATH LIBS LDFLAGS HDF5_HOME \
       H5LIBS MPI_HOME BOPT PETSC_DIR PETSC_ARCH VECPOT VMEC ARPACK SNES_OPT \
       BINDIR CONTRIBLIBS MPIEXEC

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

# contrib: vmec ftracer arpack

# contrib_clean: vmec_clean ftracer_clean arpack_clean

# contrib_setup: ftracer_setup

contrib: vmec arpack

contrib_clean: vmec_clean arpack_clean

vmec:
ifdef VMEC
	$(MAKE) -e -C contrib/vmec/LIBSTELL release INC_PATH=$(NETCDF_INC)
endif

vmec_clean:
ifdef VMEC
	$(MAKE) -e -C contrib/vmec/LIBSTELL/Release -f makelibstell clean
endif

# ftracer:
# 	$(MAKE) -e -C contrib/field_tracer

# ftracer_clean:
# 	$(MAKE) -e -C contrib/field_tracer distclean

# ftracer_setup:
# 	$(MAKE) -e -C contrib/field_tracer setup

arpack:
ifdef ARPACK
	$(MAKE) -e -C contrib/arpack PLAT=$(FC) home=$(PWD)/contrib/arpack lib
ifdef BOPT
	$(MAKE) -e -C contrib/arpack PLAT=$(FC) home=$(PWD)/contrib/arpack plib
endif
endif

arpack_clean:
ifdef ARPACK
	$(MAKE) -e -C contrib/arpack PLAT=$(FC) home=$(PWD)/contrib/arpack clean
endif

# CLEAN ALL

distclean: contrib_clean
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir distclean;  done
