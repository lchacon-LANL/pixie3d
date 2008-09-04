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
HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5
#MPI_HOME  =/usr/local/mpich2-1.0.7/absoft_

PREPROC = -D

HDF5 = t

LIBS = -llapack -lblas

# System-dependent variables

HOST?=$(HOSTNAME)

ifeq ($(HOST),quevedo.ornl.gov)
  ifdef BOPT
    PETSC_DIR =$(HOME)/lib/petsc-2.3.3

#    PETSC_ARCH = linux64_absoft_
#    HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/absoft_/parallel

     FC=gfortran
     PETSC_ARCH = linux64_openmpi
     HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/gfortran/parallel
  else
    LIBS := -L/usr/lib64 -llapack -lblas
    HDF5_HOME =$(HOME)/lib/hdf5-1.6.7/absoft_/serial
  endif
  HDF5_LIBS = -L$(HDF5_HOME)/lib -lhdf5_fortran -lhdf5 -lz -lm
endif

ifeq ($(HOST),nip.lanl.gov)
  LIBS := -llapack -lblas -lg2c 
endif

ifeq ($(HOST),ra4)
  PETSC_DIR =/ricercatori/ft/petsc-2.2.0
  HDF5_HOME =
  CPPFLAGS += $(PREPROC)RFX
endif

ifeq ($(HOST),ra22)
  PETSC_DIR =/ricercatori/ft/petsc-2.3.3
#  MPI_HOME  =/ricercatori/ft/mpich2-1.0.5/gcc-pgf90
  HDF5_HOME =
  CPPFLAGS += $(PREPROC)RFX
  FC        = pgf95
endif

ifeq ($(HOST),cayenne1)
#  MPI_HOME  =/packages/mpi/mpich2-1.0.5p4-gcc-4.1-pgi-7.0-debug
  OPT       = g
  FC        = pgf95
  PETSC_C   = t
  BOPT      = O
#  FFLAGS    += -i8
endif

ifeq ($(HOST),gongora.lanl.gov)
  PETSC_DIR =/usr/local/petsc-2.2.0
  HDF5_HOME =/usr/local/hdf5/parallel/f90_
#  MPI_HOME  =/usr/local/mpich-1.2.5.2/f90_
endif

ifeq ($(HOST),bassi)
  OPTIMIZATION = -O3 -qstrict -qarch=pwr3 -qtune=pwr3
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

  BOPT = t
endif

ifeq ($(HOST),franklin)
  FC = ftn
#  MPI_HOME  = $(MPICH_DIR)
  HDF5_HOME = $(HDF5_PAR_DIR)
  HDF5_LIBS:= $(HDF5)
  override HDF5 = t

  OPTIMIZATION = -fastsse -Mipa=fast
  DEBUG        = -g -Mbounds -Mchkptr -Ktrap=fp
  PROFILE      = -pg
  STATIC       =
  MODFLAG      = -module 
  ADDMODFLAG   = -module 
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)pgf90

  BOPT = t
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

  FFLAGS      += -w -YEXT_NAMES=LCS -YEXT_SFX=_ -YCFRL=1 -m64 -lU77
  LDFLAGS      = -lU77 
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  OPTIMIZATION = --nap --nchk --fast --npca --nsav --ntrace
  DEBUG        = -g --chk ase --warn --f95 --trap
  PROFILE      =
  STATIC       = 
  MODFLAG      = -M
  ADDMODFLAG   = -I
  VERBOSE      = --verbose
  CPPFLAGS    += $(PREPROC)lf95 

  ifdef BOPT
    FFLAGS    += --ml cdecl
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
endif

# Flags for g95
ifeq ($(FC),g95)
  OPTIMIZATION = -O2
  DEBUG        = -g -fbounds-check -ftrace=full
  PROFILE      = -pg
  STATIC       = -fstatic
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)g95
  FFLAGS      += -fstatic
endif

# Flags for gfortran
ifeq ($(FC),gfortran)
  OPTIMIZATION = -O3
  DEBUG        = -g -fbounds-check -fbacktrace
  PROFILE      = -pg
  STATIC       = -fno-automatic
  MODFLAG      = -I
  ADDMODFLAG   = -I
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)gfortran
  FFLAGS      += -fno-automatic
endif

# Flags for Portland Group f90
ifeq ($(FC),pgf95)
  OPTIMIZATION = -fastsse -Mipa=fast
  DEBUG        = -g -Mbounds -Mchkptr -Ktrap=fp
  PROFILE      = -pg
  STATIC       = 
  MODFLAG      = -module 
  ADDMODFLAG   = -module 
  VERBOSE      = -v
  CPPFLAGS    += $(PREPROC)pgf90
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
  CONTRIBLIBS += -L../contrib/vmec/lib -lstell
  CPPFLAGS    += $(PREPROC)vmec
  MODPATH     += $(ADDMODFLAG)../contrib/vmec/lib
endif

# ARPACK setup

ifdef ARPACK
  CONTRIBLIBS += -L../contrib/arpack/ -larpack_$(FC)
  CPPFLAGS    += $(PREPROC)arpack
endif

# Petsc setup

ifdef BOPT
  ifdef VECPOT
    TARGET = petsc_a
  else
    TARGET = petsc
  endif

#  CPPFLAGS += $(PREPROC)petsc $(PREPROC)NVAR=8 -I$(PETSC_DIR)/include -I${PETSC_DIR}/bmake/$(PETSC_ARCH) -I$(MPI_HOME)/include
  CPPFLAGS += $(PREPROC)petsc $(PREPROC)NVAR=8 -I$(PETSC_DIR)/include -I${PETSC_DIR}/bmake/$(PETSC_ARCH)

  ifdef PETSC_C
    CPPFLAGS += $(PREPROC)petsc_c
    SNES_OPT = -snes_mf
  endif
endif

#Export required variables

export FC FFLAGS CPPFLAGS MODFLAG ADDMODFLAG MODPATH LIBS LDFLAGS HDF5_HOME \
       H5LIBS MPI_HOME BOPT PETSC_DIR PETSC_ARCH VECPOT VMEC ARPACK SNES_OPT \
       BINDIR CONTRIBLIBS

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

contrib: vmec ftracer arpack

contrib_clean: vmec_clean ftracer_clean arpack_clean

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
