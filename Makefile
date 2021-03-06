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
COARSE_MG=t
USER_DT=t
RFX=f

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

ifeq ($(USER_DT),t)
  CPPFLAGS += $(PREPROC)USER_PROVIDED_DT
endif

ifeq ($(COARSE_MG),t)
  CPPFLAGS += $(PREPROC)coarse_MG
endif

ifeq ($(PER_BC_SYNC),t)
#  ifndef BOPT
    CPPFLAGS += $(PREPROC)PER_BC_SYNC
#  endif
endif

ifeq ($(PIT),t)
  CPPFLAGS += -Dpit
endif

ifdef BOPT
  CPPFLAGS += $(PREPROC)NVAR=8
endif

# PIXIE3D setup

GIT_BRANCH=$(shell git rev-parse --abbrev-ref HEAD)
GIT_COMMIT_HASH=$(shell git rev-parse --short HEAD)
GIT_COMMIT_TAG=$(shell git describe --tags `git rev-list --tags --max-count=1`)
GIT_COMMIT_DATE=$(shell git show -s --format=%cd --date=format:'%Y-%m-%d' $(GIT_COMMIT_HASH))

REL1=$(GIT_COMMIT_TAG)-$(GIT_COMMIT_HASH)

CPPFLAGS += $(PREPROC)REL1=\"$(REL1)\" $(PREPROC)GIT_BRANCH=\"$(GIT_BRANCH)\"
CPPFLAGS += $(PREPROC)GIT_COMMIT_HASH=\"$(GIT_COMMIT_HASH)\"
CPPFLAGS += $(PREPROC)GIT_COMMIT_DATE=\"$(GIT_COMMIT_DATE)\"

CPPFLAGS += $(PREPROC)COMPILE_HOST=\"$(MACHINE)\"

SUBDIRS = eq src plot

ifeq ($(FLUX),t)
  CPPFLAGS += -Dflux_rhs
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

export VMEC

ifeq ($(VECPOT),t)
  CPPFLAGS += $(PREPROC)vec_pot
  ifdef BOPT
    TARGET = petsc_a
  else
    TARGET = code_a
  endif
else
  ifdef BOPT
    TARGET = petsc
  else
    TARGET = code
  endif
endif

export VECPOT

ifeq ($(SAMR),t)
  TARGET = samrai
  CPPFLAGS += $(PREPROC)flux_rhs
endif

ifeq ($(RFX),t)
  CPPFLAGS += $(PREPROC)RFX
endif

#Export required variables

ifdef NOCOMMON
  COMMON = 
else
  COMMON = common
endif

export CPPFLAGS MODPATH CONTRIBLIBS BINDIR PIT REL1 COMMON #MPI_HOME

#Define targets

.PHONY: pixie3d pixplot srcclean allclean distclean petsc all contrib contrib_clean \
        vmec vmec_clean setup tests rebuild-tests testclean $(SUBDIRS)

all: $(SUBDIRS)

pixie3d: src 

pixplot: plot

$(SUBDIRS): contrib
	$(MAKE) --no-print-directory -e -C $@ $(TARGET)

# SETUP

setup: contrib_setup
	-@cd plot ; ln -s -f ../src/Makefile
	-@for subdir in $(COMMONDIR) plot eq tests/serial tests/parallel ; do \
		$(MAKE) -C $$subdir setup;  done

# TESTS

tests: ;
ifeq ($(VECPOT),t)
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
ifneq ($(VECPOT),t)
ifdef BOPT
	$(MAKE) --no-print-directory -e -C eq/test_parallel test
else
	$(MAKE) --no-print-directory -e -C eq/test_serial test
endif
endif

rebuild-tests: ;
ifeq ($(VECPOT),t)
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
ifneq ($(VECPOT),t)
ifdef BOPT
	$(MAKE) --no-print-directory -e -C eq/test_parallel rebuild
else
	$(MAKE) --no-print-directory -e -C eq/test_serial rebuild
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

srcclean:
	-@for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

allclean: contrib_clean distclean

distclean:
	-@for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir distclean;  done

testclean:
	-@for subdir in tests/serial tests/parallel ; do \
		$(MAKE) -C $$subdir clean;  done



