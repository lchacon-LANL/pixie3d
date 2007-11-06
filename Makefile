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

PETSC_DIR =/usr/local/petsc-2.3.3
HDF5_HOME =/usr/local/hdf5/parallel/mpich2_f90_
MPI_HOME  =/usr/local/mpich2-1.0.5/f90_

PREPROC = -D

# System-dependent variables

ifndef HOST
  HOST = `hostname`
endif

ifeq ($(HOST),ra4)
  PETSC_DIR =/ricercatori/ft/petsc-2.2.0
  HDF5_HOME =
  CPPFLAGS += -DRFX
endif

ifeq ($(HOST),cayenne1)
  MPI_HOME  =/packages/mpi/mpich2-1.0.5p4-gcc-4.1-pgi-7.0-debug
  OPT       = g
  FC        = pgf95
  PETSC_C   = t
  BOPT      = O
#  FFLAGS    += -i8
endif

ifeq ($(HOST),gongora.lanl.gov)
  PETSC_DIR =/usr/local/petsc-2.2.0
  HDF5_HOME =/usr/local/hdf5/parallel/f90_
  MPI_HOME  =/usr/local/mpich-1.2.5.2/f90_
endif

ifeq ($(HOST),bassi)
  FC = xlf95
endif

ifeq ($(HOST),franklin)
  FC = ftn
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
  CPPFLAGS    += $(PREPROC)absoft_ 

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
  CPPFLAGS    += $(PREPROC)lf95 

  ifdef BOPT
    FFLAGS    += --ml cdecl
    PETSC_ARCH = linux_lahey
    include ${PETSC_DIR}/bmake/common/base
    override FC = $(MPI_HOME)/bin/mpif90
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
  CPPFLAGS    += $(PREPROC)ifort

  ifdef BOPT
    PETSC_ARCH = linux_intel
    include ${PETSC_DIR}/bmake/common/base
    override FC = $(MPI_HOME)/bin/mpif90
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
  CPPFLAGS    += $(PREPROC)g95

  ifdef BOPT
    PETSC_ARCH = linux_intel
    include ${PETSC_DIR}/bmake/common/base
    override FC = $(MPI_HOME)/bin/mpif90
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
  CPPFLAGS    += $(PREPROC)pgf90

  ifdef BOPT
#    PETSC_ARCH = linux64_gcc_pgf90
    include ${PETSC_DIR}/bmake/common/base
    override  FC = $(MPI_HOME)/bin/mpif90
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

# Flags for IBM xlf95  (NERSC's Bassi)
ifeq ($(FC),xlf95)
  OPTIMIZATION = -O3 -qstrict -qarch=pwr3 -qtune=pwr3
#  DEBUG        = -g -C
  DEBUG        = -g
  PROFILE      = -P
  STATIC       = -s
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
#  FFLAGS       = -qmaxmem=-1 -bmaxstack:0x22000000 -qsave -blpdata -qsmallstack -qfixed=72
  FFLAGS       = -qfixed=72

  PREPROC      = -WF,-D
  CPPFLAGS    += $(PREPROC)xlf

  ifdef BOPT
    include ${PETSC_DIR}/bmake/common/base
    override FC = mpxlf95_r $(FFLAGS)
  else
    LDFLAGS      = 
    LIBS      += $(LAPACK) -lblas
  endif

  HDF5 = f
endif

# Flags for NERSC's Franklin ftn (wrapper for pgf90)
ifeq ($(FC),ftn)
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
    include ${PETSC_DIR}/bmake/common/base
    override  FC = ftn $(FFLAGS)
  else
#    FFLAGS       = -Minform=severe
    LDFLAGS      =
#    LIBS      += $(LAPACK) -lblas
  endif

  HDF5 = f
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

BINDIR = $(PWD)/bin

SUBDIRS = src plot

MODPATH = $(MODFLAG).

# PIXIE setup

REL1=0
REL2=9
CPPFLAGS += $(PREPROC)REL1=$(REL1) $(PREPROC)REL2=$(REL2)

ifdef VECPOT
  CPPFLAGS += $(PREPROC)vec_pot
  ifdef VECPOT
    TARGET = code_a
  endif
endif

# HDF5 setup

ifeq ($(HDF5),t)
  H5LIBS    = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5
  CPPFLAGS += $(PREPROC)hdf5 -I$(HDF5_HOME)/include
  MODPATH  += $(ADDMODFLAG)$(HDF5_HOME)/lib
endif

# VMEC setup

ifdef VMEC
  VMECLIBS  = -L../contrib/vmec/lib -lstell
  CPPFLAGS += $(PREPROC)vmec
  MODPATH  += $(ADDMODFLAG)../contrib/vmec/lib
endif

# Petsc setup

ifdef BOPT
  ifdef VECPOT
    TARGET = petsc_a
  else
    TARGET = petsc
  endif

  CPPFLAGS += $(PREPROC)petsc $(PREPROC)NVAR=8 -I$(PETSC_DIR)/include -I${PETSC_DIR}/bmake/$(PETSC_ARCH) -I$(MPI_HOME)/include

  ifdef PETSC_C
    CPPFLAGS += $(PREPROC)petsc_c
    SNES_OPT = -snes_mf
  endif
endif

#Export required variables

export FC FFLAGS CPPFLAGS MODFLAG ADDMODFLAG MODPATH LIBS LDFLAGS HDF5_HOME \
       H5LIBS MPI_HOME BOPT PETSC_DIR PETSC_ARCH VECPOT VMEC VMECLIBS SNES_OPT \
       BINDIR

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

rebuild-all-serial-tests: rebuild-serial-tests-a rebuild-serial-tests-b

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

rebuild-all-parallel-tests: rebuild-parallel-tests-a rebuild-parallel-tests-b

parallel-tests-a: ;
	$(MAKE) -C tests/parallel test-a

parallel-tests-b: ;
	$(MAKE) -C tests/parallel test-b

rebuild-parallel-tests-a: ;
	$(MAKE) -C tests/parallel rebuild-a

rebuild-parallel-tests-b: ;
	$(MAKE) -C tests/parallel rebuild-b

# CONTRIBUTED SOFTWARE

contrib: vmec ftracer

contrib_clean: vmec_clean ftracer_clean

vmec:
ifdef VMEC
	$(MAKE) -e -C contrib/vmec/LIBSTELL release
endif

vmec_clean:
ifdef VMEC
	$(MAKE) -e -C contrib/vmec/LIBSTELL/Release -f makelibstell clean
endif

ftracer:
	$(MAKE) -e -C contrib/field_tracer

ftracer_clean:
	$(MAKE) -e -C contrib/field_tracer distclean


# CLEAN ALL

distclean: contrib_clean
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir distclean;  done
