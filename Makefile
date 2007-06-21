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

# System-dependent variables

ifndef HOST
   HOST = `hostname`
endif

ifeq ($(HOST),ra4.igi.cnr.it)
  PETSC_DIR =/ricercatori/ft/petsc-2.2.0
  HDF5_HOME =
  CPPFLAGS += -DRFX
else ifeq ($(HOST),moscato.lanl.gov)
  PETSC_DIR =/usr/local/petsc-2.3.3
  HDF5_HOME =/usr/local/hdf5/parallel/mpich2_f90_
  MPI_HOME  =/usr/local/mpich2-1.0.5/f90_
else ifeq ($(HOST),cayenne1)
  MPI_HOME  =/packages/mpi/mpich2-1.0.5p4-gcc-4.1-pgi-7.0-debug
  OPT       = g
  FC        = pgf95
  PETSC_C   = t
  FFLAGS    += -i8
else
  PETSC_DIR =/usr/local/petsc-2.2.0
  HDF5_HOME =/usr/local/hdf5/parallel/f90_
endif

#Define compiler flags

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
  CPPFLAGS    += -Dabsoft_ 

  ifdef BOPT
    PETSC_ARCH = linux_absoft_
    include ${PETSC_DIR}/bmake/common/base
#    FC         = $(MPI_HOME)/bin/mpif90
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

  HDF5 = t
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  OPTIMIZATION = -O
#  DEBUG        = -g --chkglobal --warn --f95
  DEBUG        = -g --chk ase --warn --f95
#  DEBUG        = -g --f95
  PROFILE      =
  STATIC       = 
  MODFLAG      = -M
  ADDMODFLAG   = -I
  VERBOSE      = --verbose
  CPPFLAGS    += -Dlf95 

  ifdef BOPT
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
  DEBUG        = -g -check all -traceback
#  DEBUG        = -g -check -traceback
  PROFILE      = -p
  STATIC       =
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += -Difort

  ifdef BOPT
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
  CPPFLAGS    += -Dg95

  ifdef BOPT
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


# Flags for Portland Group f90
ifeq ($(FC),pgf95)
  OPTIMIZATION = -fastsse -Mipa=fast
  DEBUG        = -g -Mbounds -Mchkptr -Ktrap=fp
#  DEBUG        = -g
  PROFILE      = -pg
  STATIC       = 
  MODFLAG      = -module 
  ADDMODFLAG   = -module 
  VERBOSE      = -v
  CPPFLAGS    += -Dpgf90

  ifdef BOPT
#    PETSC_ARCH = linux64_gcc_pgf90
    include ${PETSC_DIR}/bmake/common/base
    FC         = $(MPI_HOME)/bin/mpif90
  else
#    FFLAGS       = -Minform=severe
    LDFLAGS      = 
    ifeq ($(HOST),nip.lanl.gov)
      LIBS      += -llapack -lblas -lg2c 
    else
      LIBS      += -llapack -lblas
    endif
  endif

endif

#Assemble compiler options

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

# PIXIE setup

REL1=0
REL2=9
CPPFLAGS += -DREL1=$(REL1) -DREL2=$(REL2)

ifdef VECPOT
  CPPFLAGS += -Dvec_pot
  ifdef VECPOT
    TARGET = code_a
  endif
endif

# HDF5 setup

ifdef HDF5
  H5LIBS    = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5
  CPPFLAGS += -Dhdf5 -I$(HDF5_HOME)/include
  MODPATH  += $(ADDMODFLAG)$(HDF5_HOME)/lib
endif

# VMEC setup

ifdef VMEC
  VMECLIBS  = -L../contrib/vmec/lib -lstell
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

  ifdef PETSC_C
    CPPFLAGS += -Dpetsc_c
  endif
endif

#Export required variables

export FC FFLAGS CPPFLAGS MODFLAG ADDMODFLAG MODPATH LIBS LDFLAGS HDF5_HOME \
       H5LIBS MPI_HOME BOPT PETSC_DIR PETSC_ARCH VECPOT VMEC VMECLIBS

#Define targets

.PHONY: pixie3d pixplot distclean petsc all contrib contrib_clean \
        vmec vmec_clean setup serial-tests-b serial-tests-a \
        rebuild-serial-tests-a rebuild-serial-tests-b \
        parallel-tests-b parallel-tests-a \
        rebuild-parallel-tests-a rebuild-parallel-tests-b $(SUBDIRS)

all: contrib src plot

pixie3d: src

pixplot: plot

$(SUBDIRS):
	$(MAKE) -e -C $@ $(TARGET)

# SETUP

setup: ;
	-cd plot ; ln -s ../src/Makefile
	-for subdir in common plot tests/serial tests/parallel ; do \
		$(MAKE) -C $$subdir setup;  done

# SERIAL TESTS

all-serial-tests: serial-tests-a serial-tests-b

serial-tests-a: ;
	$(MAKE) -C tests/serial test-a

serial-tests-b: ;
	$(MAKE) -C tests/serial test-b

rebuild-serial-tests-a: ;
	$(MAKE) -C tests/serial rebuild-a

rebuild-serial-tests-b: ;
	$(MAKE) -C tests/serial rebuild-b

# PARALLEL TESTS

all-parallel-tests: parallel-tests-a parallel-tests-b

parallel-tests-a: ;
	$(MAKE) -C tests/parallel test-a

parallel-tests-b: ;
	$(MAKE) -C tests/parallel test-b

rebuild-parallel-tests-a: ;
	$(MAKE) -C tests/parallel rebuild-a

rebuild-parallel-tests-b: ;
	$(MAKE) -C tests/parallel rebuild-b

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

# CLEAN ALL

distclean: contrib_clean
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir distclean;  done
