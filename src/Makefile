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
	createNonlinearFunction.f \
	initializePlots.f \
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
COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.f))
COMMON_OBJS = $(COMMON_SRC:.f=.o)

prefix = .

#Define compiler flags

FC = f90

ifeq ($(FC),f90)

# Flags for Absoft f90
OPTIMIZATION = -O -B101 -cpu:host
#DEBUG = -g -et
#DEBUG = -g -et -trap=ALL
DEBUG = -g
PROFILE = -P
STATIC = -s
MODPATH = -p

LIBS = -llapack_f90 -lblas_f90

else

# Flags for Lahey lf95
OPTIMIZATION = -O --tpp
DEBUG = -g --trace --trap
#DEBUG = -g 
PROFILE =
STATIC =
MODPATH = -M

LIBS = -llapack_lf95 -lblas_lf95

endif

#Setup compiler options

OPT = O 

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

COMMONDIR = ./commondirs

SUBDIRS = $(COMMONDIR)/driver \
	  $(COMMONDIR)/grid \
	  $(COMMONDIR)/newtonGMRES \
	  $(COMMONDIR)/interpolation/1d_spline \
	  $(COMMONDIR)/interpolation/xerror

MODDIRS = $(patsubst $(COMMONDIR)%,$(MODPATH) $(COMMONDIR)%,$(SUBDIRS))

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
		$(MAKE) -e -C $$subdir ; done

clean:
	-rm -f *.o *.mod

recclean: clean
	for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Define dependencies

$(OBJS):  modules.f $(COMMOM_MODS)

#Define patterns

%.o : %.f
	$(FC) -c $(MODDIRS) $(FFLAGS) $<
