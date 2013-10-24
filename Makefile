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

FPA=t
FLUX=f

COMMONDIR =$(PWD)/common
CONTRIBDIR=$(PWD)/contrib
BINDIR    =$(PWD)/bin

# System-dependent variables

-include $(COMMONDIR)/make/make.mach.inc

# Define compiler flags

-include $(COMMONDIR)/make/make.comp.inc

# External packages configuration

-include $(COMMONDIR)/make/make.lib.inc

# PIXIE3D setup

SUBDIRS = src plot

REL1=2
REL2=0
CPPFLAGS += $(PREPROC)REL1=$(REL1) $(PREPROC)REL2=$(REL2)

ifeq ($(PIT),t)
  PIT_DIR   = $(CONTRIBDIR)/parareal
  CPPFLAGS += -Dpit
endif

ifeq ($(FLUX),t)
  CPPFLAGS += -Dflux_rhs
endif

ifeq ($(VECPOT),t)
  CPPFLAGS += $(PREPROC)vec_pot
  TARGET = code_a
endif

ifeq ($(PER_BC_SYNC),t)
#  ifndef BOPT
    CPPFLAGS += $(PREPROC)PER_BC_SYNC
#  endif
endif

VMEC=t
ifeq ($(VMEC),t)
  VMEC_DIR     = $(CONTRIBDIR)/vmec/LIBSTELL

  CONTRIBLIBS += -L$(CONTRIBDIR)/vmec/lib -lstell
  CPPFLAGS    += $(PREPROC)vmec
  MODPATH     += $(ADDMODFLAG)$(CONTRIBDIR)/vmec/lib

  ifeq ($(NETCDF),t)
    CPPFLAGS   += $(PREPROC)NETCDF $(NETCDF_INC)
    CONTRIBLIBS += $(NETCDF_LIBS)
  endif
endif

ifdef BOPT
  CPPFLAGS += $(PREPROC)NVAR=8

  ifdef VECPOT
    TARGET = petsc_a
  else
    TARGET = petsc
  endif
endif

ifeq ($(SAMR),t)
  TARGET = samrai
  CPPFLAGS += -Dflux_rhs
endif

#Export required variables

export FC FFLAGS CPPFLAGS MODFLAG ADDMODFLAG MODPATH LIBS LDFLAGS \
       H5LIBS MPI_HOME BOPT PETSC_DIR PETSC_ARCH VECPOT VMEC ARPACK FPA SNES_OPT \
       BINDIR CONTRIBLIBS MPIEXEC FLINKER PETSC_SNES_LIB CPPFLAGS_EXTRA \
       CXXFLAGS_EXTRA LDFLAGS_EXTRA LDLIBS_EXTRA SAMR SAMRAI LIBSAMRAI3D LIBSAMRAI \
       AMRUTILITIES_HOME SAMRSOLVERS_HOME HOST

#Define targets

.PHONY: pixie3d pixplot distclean petsc all contrib contrib_clean \
        vmec vmec_clean setup parareal \
        serial-tests-b serial-tests-a \
        rebuild-serial-tests-a rebuild-serial-tests-b \
        parallel-tests-b parallel-tests-a \
        rebuild-parallel-tests-a rebuild-parallel-tests-b $(SUBDIRS)

all: contrib parareal src plot

pixie3d: contrib src 

pixplot: contrib plot

$(SUBDIRS):
#	@echo "-- Machine = $(HOST)"
#	@echo "-- Do make in $@ with FC=$(FC)"
#	@echo "-- HDF5_HOME=$(HDF5_HOME)"
#	@echo "-- HDF5_LIBS=$(HDF5_LIBS)"
#	@echo "-- MODPATH=$(MODPATH)"
	$(MAKE) --no-print-directory -e -C $@ $(TARGET)

# SETUP

setup:
	-@cd plot ; ln -s -f ../src/Makefile
	-@tar xzf contrib.tgz
	-@for subdir in $(COMMONDIR) plot tests/serial tests/parallel tests/samrai; do \
		$(MAKE) -C $$subdir setup;  done

# TESTS

tests: tests-a tests-b

rebuild-tests: rebuild-tests-a rebuild-tests-b

tests-a: ;
ifdef SAMR
	$(MAKE) --no-print-directory -e -C tests/samrai test-a
else
ifdef BOPT
	$(MAKE) --no-print-directory -e -C tests/parallel test-a
else
	$(MAKE) --no-print-directory -e -C tests/serial test-a
endif
endif

tests-b: ;
ifdef SAMR
	$(MAKE) --no-print-directory -e -C tests/samrai test-b
else
ifdef BOPT
	$(MAKE) --no-print-directory -e -C tests/parallel test-b
else
	$(MAKE) --no-print-directory -e -C tests/serial test-b
endif
endif

rebuild-tests-a: ;
ifdef SAMR
	$(MAKE) --no-print-directory -e -C tests/samrai rebuild-a
else
ifdef BOPT
	$(MAKE) --no-print-directory -e -C tests/parallel rebuild-a
else
	$(MAKE) --no-print-directory -e -C tests/serial rebuild-a
endif
endif

rebuild-tests-b: ;
ifdef SAMR
	$(MAKE) --no-print-directory -e -C tests/samrai rebuild-b
else
ifdef BOPT
	$(MAKE) --no-print-directory -e -C tests/parallel rebuild-b
else
	$(MAKE) --no-print-directory -e -C tests/serial rebuild-b
endif
endif

# COMMON CONTRIBUTED LIBRARIES

contrib:
	$(MAKE) --no-print-directory -e -C $(COMMONDIR) contrib
ifeq ($(VMEC),t)
	$(MAKE) --no-print-directory -e -C $(VMEC_DIR) release INC_PATH=$(NETCDF_INC)
endif
ifeq ($(PIT),t)
	$(MAKE) --no-print-directory -e -C $(PIT_DIR) all	
endif

contrib_clean:
	$(MAKE) --no-print-directory -e -C $(COMMONDIR) contrib_clean
ifeq ($(VMEC),t)
	$(MAKE) --no-print-directory -e -C $(VMEC_DIR)/Release -f makelibstell clean
endif
ifeq ($(PIT),t)
	$(MAKE) --no-print-directory -e -C $(PIT_DIR) clean	
endif

# CLEAN ALL

allclean: contrib_clean distclean

distclean:
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) --no-print-directory -C $$subdir distclean;  done


