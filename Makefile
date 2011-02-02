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
#     + FC = gfortran uses Gfortran.
#     + ... (others documented in common/make/make.comp.inc)
# When PETSc is employed, the optimization level should be specified by
# the variable $(BOPT).
#
# A call example that uses lf95 in debugging mode is:
#
#        make FC=lf95 OPT=g
#
# A call example that employs PETSc is:
#
#        make BOPT=g pixie3d

# System-dependent variables

-include common/make/make.mach.inc

#Define compiler flags

-include common/make/make.comp.inc

#Define relevant directories

BINDIR = $(PWD)/bin

SUBDIRS = src plot

MODPATH = $(MODFLAG).

# PIXIE setup

REL1=1
REL2=5
CPPFLAGS += $(PREPROC)REL1=$(REL1) $(PREPROC)REL2=$(REL2)

ifdef VECPOT
  CPPFLAGS += $(PREPROC)vec_pot
  TARGET = code_a
endif

ifdef PER_BC_SYNC
#  ifndef BOPT
    CPPFLAGS += $(PREPROC)PER_BC_SYNC
#  endif
endif

# ADIOS setup

ifeq ($(ADIOS),t)
  CONTRIBLIBS += $(ADIOS_LIBS)
  CPPFLAGS   += $(PREPROC)adios -I$(ADIOS_HOME)/include
#  MODPATH    += $(ADDMODFLAG)$(ADIOS_HOME)/include
endif

# HDF5 setup

ifeq ($(HDF5),t)
  CONTRIBLIBS += $(HDF5_LIBS) 
  CPPFLAGS    += $(PREPROC)hdf5 $(PREPROC)H5_USE_16_API $(HDF5_INC)
  MODPATH     += $(ADDMODFLAG)$(HDF5_MOD)
endif

# VMEC setup

ifeq ($(VMEC),t)
  VMEC_DIR     = contrib/vmec/LIBSTELL

  CONTRIBLIBS += -L$(PWD)/contrib/vmec/lib -lstell
  CPPFLAGS    += $(PREPROC)vmec
  MODPATH     += $(ADDMODFLAG)$(PWD)/contrib/vmec/lib

  ifeq ($(NETCDF),t)
    CPPFLAGS   += $(PREPROC)NETCDF $(NETCDF_INC)
    CONTRIBLIBS += $(NETCDF_LIBS)
  endif
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

   SUBDIRS += common/driver-samrai

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
   LDFLAGS_EXTRA += -L${SAMRSOLVERS_HOME}/lib   
   LDLIBS_EXTRA += ${SAMRSOLVERS_HOME}/lib/libtimeintegrators3d.a
   #LDLIBS_EXTRA += -lAMRUtils${PDIM}d
endif

#Export required variables
export FC FFLAGS CPPFLAGS MODFLAG ADDMODFLAG MODPATH LIBS LDFLAGS \
       H5LIBS MPI_HOME BOPT PETSC_DIR PETSC_ARCH VECPOT VMEC ARPACK FPA SNES_OPT \
       BINDIR CONTRIBLIBS MPIEXEC FLINKER PETSC_SNES_LIB CPPFLAGS_EXTRA \
       CXXFLAGS_EXTRA LDFLAGS_EXTRA LDLIBS_EXTRA SAMR SAMRAI LIBSAMRAI3D LIBSAMRAI \
       AMRUTILITIES_HOME SAMRSOLVERS_HOME

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
#	@echo "-- Machine = $(MACHINE)"
#	@echo "-- Do make in $@ with FC=$(FC)"
#	@echo "-- HDF5_HOME=$(HDF5_HOME)"
#	@echo "-- HDF5_LIBS=$(HDF5_LIBS)"
#	@echo "-- MODPATH=$(MODPATH)"
	$(MAKE) -e -C $@ $(TARGET)


# SETUP

setup:
	-cd plot ; ln -s ../src/Makefile
	-tar xzf contrib.tgz
	-for subdir in common plot tests/serial tests/parallel ; do \
		$(MAKE) -C $$subdir setup;  done

# TESTS

tests: tests-a tests-b

rebuild-tests: rebuild-tests-a rebuild-tests-b

tests-a: ;
ifdef BOPT
	$(MAKE) -e -C tests/parallel test-a
else
	$(MAKE) -e -C tests/serial test-a
endif

tests-b: ;
ifdef SAMR
	$(MAKE) -e -C tests/samrai test-b
else
ifdef BOPT
	$(MAKE) -e -C tests/parallel test-b
else
	$(MAKE) -e -C tests/serial test-b
endif
endif

rebuild-tests-a: ;
ifdef BOPT
	$(MAKE) -e -C tests/parallel rebuild-a
else
	$(MAKE) -e -C tests/serial rebuild-a
endif

rebuild-tests-b: ;
ifdef BOPT
	$(MAKE) -e -C tests/parallel rebuild-b
else
	$(MAKE) -e -C tests/serial rebuild-b
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


