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

SRC  := $(wildcard *.f) 

MODS := $(wildcard *_mod.f)

OBJS := $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))

OBJMOD := $(MODS:.f=.o)

COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.f)))
COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.f))
COMMON_OBJS = $(COMMON_SRC:.f=.o)

prefix = .

#Define compiler flags

FC = f90

# Flags for Absoft f90
ifeq ($(FC),f90)
  OPTIMIZATION = -O2 -cpu:host
  DEBUG = -g -et -Rb
  #DEBUG = -g
  PROFILE = -P
  STATIC = -s
  MODPATH = -p
  FFLAGS = -w

  LIBS = -llapack_f90 -lblas_f90

  MODDIRS = $(patsubst $(COMMONDIR)%,$(MODPATH)$(COMMONDIR)%,$(SUBDIRS))
endif

# Flags for Lahey lf95
ifeq ($(FC),lf95)
  OPTIMIZATION = -O
#  DEBUG = -g --chkglobal
  DEBUG = -g --chk s
  PROFILE =
  STATIC = 
  MODPATH = -M
  FFLAGS =

#  LIBS = -llapackmt -lblasmt
  LIBS = -llapack -lblas -lf2c

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
  FFLAGS = -vec_report0

  LIBS = -llapack -lblas -lf2c
#  LIBS = -llapack_intel -lblas_intel

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

export FC FFLAGS

#Define linker flags

LDFLAGS = 

#Define subdirectories for common code

COMMONDIR = ./commondirs

SUBDIRS = $(COMMONDIR)/driver\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/interpolation/1d_spline\
	  $(COMMONDIR)/interpolation/2d_spline\
	  $(COMMONDIR)/interpolation/3d_spline\
	  $(COMMONDIR)/interpolation/spline_common\
	  $(COMMONDIR)/interpolation/xerror\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG

#Define targets

beta: common $(BETA)

$(BETA):  $(OBJMOD) $(OBJS)  $(COMMON_SRC)
	@echo 'Linking BETA 3DMHD code'
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(COMMON_OBJS) $(OBJS)  $(OBJMOD) $(LIBS)

production: common $(PROD)

$(PROD):  $(OBJMOD) $(OBJS)  $(COMMON_SRC)
	@echo ''
	@echo 'Linking PRODUCTION 3DMHD code'
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(COMMON_OBJS) $(OBJS)  $(OBJMOD) $(LIBS)

common: ;
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -e -C $$subdir ; done

clean:
	-rm -f *.o *.mod

recclean: clean
	for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Define dependencies

$(OBJS):  $(MODS) $(COMMOM_MODS)

#Define patterns

%.o : %.f
	$(FC) -c $(MODDIRS) $(FFLAGS) $<
