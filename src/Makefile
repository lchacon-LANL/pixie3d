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

# Source files

PROD      = pixie3d.x
BETA      = pixie3d.beta.x
PETSCCODE = pixie3d.petsc.x

SRC  := $(wildcard *.[F,c,f])

MODS := $(wildcard *_mod.[F,f])
MODS := mhd_mod.f imposeBC_mod.F precond_mod.f

prefix = .

# Petsc include

DRIVERDIR = driver
ifdef BOPT
  include ${PETSC_DIR}/bmake/common/base
  DRIVERDIR = driver-petsc
  SRC      := $(filter-out explicit.f,$(wildcard *.[F,c,f]))
endif

#Objects

OBJ    := $(filter %.o, $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.c,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.F,%.o,$(filter-out $(MODS),$(SRC))))

OBJMOD := $(filter %.o, $(patsubst %.f,%.o,$(MODS))\
                        $(patsubst %.F,%.o,$(MODS)))

COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.[f,F]))
COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(dir)/test.f,$(wildcard $(dir)/*.[f,F,c])))
COMMON_OBJS = $(filter %.o, $(patsubst %.f,%.o,$(COMMON_SRC))\
                            $(patsubst %.c,%.o,$(COMMON_SRC))\
                            $(patsubst %.F,%.o,$(COMMON_SRC)))


#export FC FFLAGS CPPFLAGS MODPATH ADDMODPATH GMODPATH PETSC_DIR PETSC_ARCH

#Define linker flags

#LDFLAGS = 

#Define relevant directories

BINDIR = ../bin

COMMONDIR = ../common

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

#Module search path

MODDIRS = $(MODPATH) $(patsubst $(COMMONDIR)%,$(ADDMODFLAG)$(COMMONDIR)%,$(SUBDIRS))

#Define targets

.PHONY: clean distclean common beta production petsc $(SUBDIRS)

beta: common $(BETA)

$(BETA):  $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	@$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(LIBS)
	@mv $@ $(BINDIR)

production: common $(PROD)

$(PROD):  $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	@$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(LIBS)
	@mv $@ $(BINDIR)

petsc: common $(PETSCCODE)

$(PETSCCODE): $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	${FLINKER} -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(LIBS) ${PETSC_FORTRAN_LIB} ${PETSC_SNES_LIB}
	@mv $@ $(BINDIR)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -e -C $@ target

clean:
	-rm -f *.o *.mod

distclean: clean
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Define dependencies

$(OBJ):  $(MODS) $(COMMON_MODS) $(COMMONDIR)/BC/imposeBCinterface_mod.F
imposeBC_mod.o: $(COMMONDIR)/BC/imposeBCinterface_mod.F

#Define patterns

%.o : %.f
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $<

%.o : %.F
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $(CPPFLAGS) $<
