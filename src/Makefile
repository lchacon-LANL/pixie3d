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
# A call example that uses lf95 in debugging mode is:
#
#        make FC='lf95' OPT='g'

PROD =	3dmhd.x
BETA =	3dmhd.beta.x
PETSCCODE = 3dmhd.petsc.x

SRC  := $(wildcard *.[F,c,f])

MODS := $(wildcard *_mod.[F,f])
MODS := mhd_mod.f imposeBC_mod.F precond_mod.f

OBJ  := $(filter %.o, $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))\
                      $(patsubst %.c,%.o,$(filter-out $(MODS),$(SRC)))\
                      $(patsubst %.F,%.o,$(filter-out $(MODS),$(SRC))))

OBJMOD := $(filter %.o, $(patsubst %.f,%.o,$(MODS))\
                        $(patsubst %.F,%.o,$(MODS)))

COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.[f,F,c])))
COMMON_OBJS = $(filter %.o, $(patsubst %.f,%.o,$(COMMON_SRC))\
                            $(patsubst %.c,%.o,$(COMMON_SRC))\
                            $(patsubst %.F,%.o,$(COMMON_SRC)))

prefix = .

# Petsc include

DRIVERDIR = driver
ifneq ($(strip $(BOPT)),)
  include ${PETSC_DIR}/bmake/common/base
  DRIVERDIR = driver-petsc
  PETSCTRUE = true
  FFLAGS    = -Dpetsc
endif

#Define compiler flags

FC = f90

# Flags for Absoft f90
ifeq ($(FC),f90)
  OPTIMIZATION = -O2 -cpu:host
  DEBUG        = -g -et -Rb -Rp -Rc
#  DEBUG        = -g
  PROFILE      = -P
  STATIC       = -s
  MODPATH      = -p
  FFLAGS       += -w
  LIBS         = -llapack_f90 -lblas_f90

  ifdef PETSCTRUE
    CPPFLAGS += -Dabsoft -DNVAR=8
    LIBS    = 
    FFLAGS += -DNVAR=8
  endif

  MODDIRS = $(patsubst $(COMMONDIR)%,$(MODPATH)$(COMMONDIR)%,$(SUBDIRS))
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  OPTIMIZATION = -O
#  DEBUG = -g --chkglobal
  DEBUG = -g --chk ase
  PROFILE =
  STATIC = 
  MODPATH = -M
  FFLAGS +=

  LIBS = -llapackmt -lblasmt
#  LIBS = -llapack -lblas -lf2c

  ifdef PETSCTRUE
    CPPFLAGS += -Dlahey -DNVAR=8
    FFLAGS += -DNVAR=8 --ml cdecl
    PETSC_ARCH = linux_lahey
  endif

  MODDIRS = $(MODPATH) . $(patsubst $(COMMONDIR)%,-I $(COMMONDIR)%,$(SUBDIRS))
endif

# Flags for Intel ifort
ifeq ($(FC),ifort)
  OPTIMIZATION = -O2 -mp -axW
#  DEBUG = -g
  DEBUG = -g -check -traceback
  PROFILE = -gp
  STATIC =
  MODPATH =
  FFLAGS += -vec_report0

  LIBS = -llapack -lblas -lf2c
#  LIBS = -llapack_intel -lblas_intel

  ifdef PETSCTRUE
    CPPFLAGS += -Dlahey -DNVAR=8
    FFLAGS += -DNVAR=8
    PETSC_ARCH = linux_intel
  endif

  MODDIRS = $(patsubst $(COMMONDIR)%,-I $(COMMONDIR)%,$(SUBDIRS))
endif

#Setup compiler options

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

export FC FFLAGS CPPFLAGS

#Define linker flags

LDFLAGS = 

#Define subdirectories for common code

COMMONDIR = ./commondirs

SUBDIRS = $(COMMONDIR)/$(DRIVERDIR)\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/interpolation/slatec/1d_spline\
	  $(COMMONDIR)/interpolation/slatec/2d_spline\
	  $(COMMONDIR)/interpolation/slatec/3d_spline\
	  $(COMMONDIR)/interpolation/slatec/spline_common\
	  $(COMMONDIR)/interpolation/slatec/xerror\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG

#Define targets

beta: common $(BETA)

$(BETA):  $(OBJMOD) $(OBJ)  $(COMMON_SRC)
	@echo 'Linking BETA 3DMHD code'
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(LIBS)

production: common $(PROD)

$(PROD):  $(OBJMOD) $(OBJ)  $(COMMON_SRC)
	@echo ''
	@echo 'Linking PRODUCTION 3DMHD code'
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(LIBS)

petsc: common $(PETSCCODE)

$(PETSCCODE): $(OBJMOD) $(OBJ)  $(COMMON_SRC)
	@echo ''
	@echo 'Linking PETSC code'
	-${FLINKER} -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(LIBS)  ${PETSC_SNES_LIB}

common: ;
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -e -C $$subdir target; done

clean:
	-rm -f *.o *.mod

recclean: clean
	for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Define dependencies

$(OBJ):  $(MODS) $(COMMOM_MODS)

#Define patterns

%.o : %.f
	$(FC) -c $(MODDIRS) $(FFLAGS) $<

%.o : %.F
	$(FC) -c $(MODDIRS) $(FFLAGS) $<