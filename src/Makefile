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

SRC  := modules.f \
	explicit.f \
	readInput.f \
	allocateSpecificVariables.f \
	initializeTimePlots.f \
	initialize3Dplots.f \
	evaluateDiagnostics.f \
	output.f \
	setupNonlinearFunction.f \
	nonlinearRHS.f \
	applyBoundaryConditions.f \
	setupPreconditioner.f \
	applyPreconditioner.f \
	killPreconditioner.f \
	setEquilibrium.f

OBJS := $(SRC:.f=.o)

COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.f)))
#COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*.f))
#COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.f))
#COMMON_OBJS = $(patsubst %.f,%.o,$(filter-out $(COMMON_MODS),$(COMMON_SRC)))
COMMON_OBJS = $(COMMON_SRC:.f=.o)
#COMMON_OBJS = $(patsubst %.f,%.o,$(filter-out test.f,$(COMMON_SRC)))

LIBS = 

prefix = .

#Define compiler flags

FC = f90

ifeq ($(FC),f90)

# Flags for Absoft f90
OPTIMIZATION = -O -B101 -cpu:host
#DEBUG = -g -et -trap=ALL
DEBUG = -g
PROFILE = -P
STATIC = -s

else

# Flags for Lahey lf95
OPTIMIZATION = -O --tpp
DEBUG = -g --trace --trap
#DEBUG = -g 
PROFILE =
STATIC =

endif

#Setup compiler options

OPT = O s

FFLAGS =

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

DRIVERDIR = ../../3d_driver

SUBDIRS = $(DRIVERDIR) \
	  $(DRIVERDIR)/newtonGMRES \

#Define targets

beta: common $(BETA)

$(BETA):  $(OBJS)  $(COMMON_SRC)
	@echo 'Linking BETA 3DMHD code'
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(COMMON_OBJS) $(OBJS) $(LIBS)

production: common $(PROD)

$(PROD):  $(OBJS)  $(COMMON_SRC)
	@echo ''
	@echo 'Linking PRODUCTION 3DMHD code'
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(COMMON_OBJS) $(OBJS) $(LIBS)

common: ;
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -e -C $$subdir ; \
		cp $$subdir/*.mod $(prefix) >& /dev/null ; done

clean:
	-rm -f *.o *.mod

recclean: clean
	for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Define dependencies

$(OBJS):  modules.f $(DRIVERDIR)/3d_driver_mod.f $(DRIVERDIR)/graphics_mod.f $(DRIVERDIR)/grid_mod.f

#Define patterns

%.o : %.f
	$(FC) -c $(FFLAGS) $<
