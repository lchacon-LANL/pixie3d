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

PROD =	3dplot.x

SRC  := $(wildcard *.[F,c,f])

MODS := $(wildcard *_mod.[F,f])
MODS := mhd_mod.f imposeBC_mod.F

OBJS := $(filter %.o, $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))\
                      $(patsubst %.c,%.o,$(filter-out $(MODS),$(SRC)))\
                      $(patsubst %.F,%.o,$(filter-out $(MODS),$(SRC))))

OBJMOD := $(filter %.o, $(patsubst %.f,%.o,$(MODS))\
                        $(patsubst %.F,%.o,$(MODS)))

COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.[f,F]))
COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.[f,F,c])))
COMMON_OBJS = $(filter %.o, $(patsubst %.f,%.o,$(COMMON_SRC))\
                            $(patsubst %.c,%.o,$(COMMON_SRC))\
                            $(patsubst %.F,%.o,$(COMMON_SRC)))

prefix = .

#Define compiler flags

FC = f90

ifeq ($(FC),f90)

# Flags for Absoft f90
OPTIMIZATION = -O1
DEBUG = -g -et -Rp
#DEBUG = -g -et --trace=ALL -Rb
#DEBUG = -g
PROFILE = -P
STATIC = -s
MODPATH = -p
FFLAGS = -w

LIBS = -lU77 -llapack_f90 -lblas_f90

MODDIRS = $(patsubst $(COMMONDIR)%,$(MODPATH)$(COMMONDIR)%,$(SUBDIRS))

else

# Flags for Lahey lf95
OPTIMIZATION = -O --tpp
DEBUG = -g --chkglobal
#DEBUG = -g --trace --trap
PROFILE =
STATIC =
MODPATH = -M
FFLAGS =

LIBS = 

MODDIRS = $(MODPATH) . $(patsubst $(COMMONDIR)%,-I $(COMMONDIR)%,$(SUBDIRS))

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

SUBDIRS = $(COMMONDIR)/grid\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/interpolation/slatec/1d_spline\
	  $(COMMONDIR)/interpolation/slatec/2d_spline\
	  $(COMMONDIR)/interpolation/slatec/3d_spline\
	  $(COMMONDIR)/interpolation/slatec/spline_common\
	  $(COMMONDIR)/interpolation/slatec/xerror\
	  $(COMMONDIR)/3dplot

#Define targets

.PHONY: clean recclean common petsc $(SUBDIRS)

production: common $(PROD)

$(PROD):  $(OBJMOD) $(OBJS) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking 3DPLOT code'
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(COMMON_OBJS) $(OBJS)  $(OBJMOD) $(LIBS)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -e -C $@ target

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

%.o : %.F
	$(FC) -c $(MODDIRS) $(FFLAGS) $<
