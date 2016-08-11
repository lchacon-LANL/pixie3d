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

FLUX=f
PER_BC_SYNC=t
VMEC=t
COARSE_MG = t

COMMONDIR =$(PWD)/common
CONTRIBDIR=$(PWD)/contrib
BINDIR    =$(PWD)/bin

# System-dependent variables

-include $(COMMONDIR)/make/make.mach.inc

# Define compiler flags

-include $(COMMONDIR)/make/make.comp.inc

# External packages configuration

-include $(COMMONDIR)/make/make.lib.inc

# FRAMEWORK setup

# PIXIE3D setup

SUBDIRS = eq src plot

REL1=3
REL2=3.8
CPPFLAGS += $(PREPROC)REL1=$(REL1) $(PREPROC)REL2=$(REL2)

ifeq ($(COARSE_MG),t)
  CPPFLAGS += $(PREPROC)coarse_MG
endif

ifeq ($(PIT),t)
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

  ifeq ($(VECPOT),t)
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
       H5LIBS MPI_HOME BOPT PETSC_DIR PETSC_ARCH VECPOT VMEC ARPACK SNES_OPT \
       BINDIR CONTRIBLIBS MPIEXEC FLINKER PETSC_SNES_LIB CPPFLAGS_EXTRA \
       CXXFLAGS_EXTRA LDFLAGS_EXTRA LDLIBS_EXTRA SAMR SAMRAI LIBSAMRAI3D LIBSAMRAI \
       AMRUTILITIES_HOME SAMRSOLVERS_HOME HOST PIT PREPROC

#Define targets

.PHONY: pixie3d pixplot distclean petsc all contrib contrib_clean \
        vmec vmec_clean setup tests rebuild-tests \
        $(SUBDIRS)

all: contrib $(SUBDIRS)

pixie3d: contrib src 

pixplot: contrib plot

$(SUBDIRS):
	$(MAKE) --no-print-directory -e -C $@ $(TARGET)

# SETUP

setup: contrib_setup
	-@cd plot ; ln -s -f ../src/Makefile
	-@for subdir in $(COMMONDIR) plot eq tests/serial tests/parallel ; do \
		$(MAKE) -C $$subdir setup;  done

# TESTS

tests: ;
ifdef VECPOT
ifdef SAMR
	$(MAKE) --no-print-directory -e -C tests/samrai test-a
else
ifdef BOPT
	$(MAKE) --no-print-directory -e -C tests/parallel test-a
else
	$(MAKE) --no-print-directory -e -C tests/serial test-a
endif
endif
else
ifdef SAMR
	$(MAKE) --no-print-directory -e -C tests/samrai test-b
else
ifdef BOPT
	$(MAKE) --no-print-directory -e -C tests/parallel test-b
else
	$(MAKE) --no-print-directory -e -C tests/serial test-b
endif
endif
endif

rebuild-tests: ;
ifdef VECPOT
ifdef SAMR
	$(MAKE) --no-print-directory -e -C tests/samrai rebuild-a
else
ifdef BOPT
	$(MAKE) --no-print-directory -e -C tests/parallel rebuild-a
else
	$(MAKE) --no-print-directory -e -C tests/serial rebuild-a
endif
endif
else
ifdef SAMR
	$(MAKE) --no-print-directory -e -C tests/samrai rebuild-b
else
ifdef BOPT
	$(MAKE) --no-print-directory -e -C tests/parallel rebuild-b
else
	$(MAKE) --no-print-directory -e -C tests/serial rebuild-b
endif
endif
endif

# COMMON CONTRIBUTED LIBRARIES

contrib:
ifeq ($(VMEC),t)
	$(MAKE) --no-print-directory -e -C $(VMEC_DIR) release INC_PATH=$(NETCDF_INC)
endif
	$(MAKE) --no-print-directory -e -C $(COMMONDIR) contrib

contrib_clean:
ifeq ($(VMEC),t)
	@$(MAKE) -e -C $(VMEC_DIR)/Release -f makelibstell clean
endif
	@$(MAKE) -e -C $(COMMONDIR) contrib_clean

contrib_setup: ;
	-@tar xzf contrib.tgz

contrib_pack: ;
	-@rm -f contrib.tgz > /dev/null
	-@tar czf contrib.tgz contrib

# CLEAN ALL

allclean: contrib_clean distclean

distclean:
	-@for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir distclean;  done


