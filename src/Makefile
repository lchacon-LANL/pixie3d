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

-include make.inc

# Exclude files

ifeq ($(VMEC),t)
  EXCLUDE += %test.f
else
  EXCLUDE += %test.f %vmec.f90
endif

# Source files

SRC ?= $(filter-out $(EXCLUDE),$(wildcard *.[F,c,f] *.?90))

MODS ?= $(filter-out $(EXCLUDE),$(wildcard *_mod.[F,f] *_mod.?90))

prefix = .

# Petsc include

DRIVERDIR = driver
ifdef BOPT
  SRC2      = $(SRC)
  ifdef PETSC_C
    EXCLUDE  += %explicit.F
  endif
  SRC      := $(filter-out $(EXCLUDE),$(SRC2))
endif

ifdef SAMR
    # Include all dependencies from SAMRSOLVERS (this will include SAMRUTILS and SAMRAI)
    include ${SAMRSOLVERS_HOME}/Makefile.config
    SRC2      = $(SRC)
    EXCLUDE  += %output.F
    ifdef SAMR_IMP
      EXCLUDE += %SAMRAIDriver.C
    else
      EXCLUDE += %SAMRAIDriverImplicit.C %ImplicitPixie3dApplication.C
    endif
    SRC      := $(filter-out $(EXCLUDE),$(SRC2))
endif

# Objects

OBJ    := $(filter %.o, $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.c,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.C,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.F,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.f90,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.F90,%.o,$(filter-out $(MODS),$(SRC))))


OBJMOD := $(patsubst %.f,%.o,$(patsubst %.F,%.f,$(patsubst %.f90,%.f,$(patsubst %.F90,%.f,$(MODS)))))

COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.[f,F] $(dir)/*_mod.f90 $(dir)/*_mod.F90))
COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(EXCLUDE),$(wildcard $(dir)/*.[f,F,c,C] $(dir)/*.f90 $(dir)/*.F90)))
COMMON_OBJS = $(filter %.o, $(patsubst %.f,%.o,$(COMMON_SRC))\
                            $(patsubst %.c,%.o,$(COMMON_SRC))\
                            $(patsubst %.C,%.o,$(COMMON_SRC))\
                            $(patsubst %.F,%.o,$(COMMON_SRC))\
                            $(patsubst %.f90,%.o,$(COMMON_SRC))\
                            $(patsubst %.F90,%.o,$(COMMON_SRC)))

#Define linker flags

#LDFLAGS = 

#Module search path

MODDIRS = $(MODPATH) $(patsubst $(COMMONDIR)%,$(ADDMODFLAG)$(COMMONDIR)%,$(SUBDIRS))

#Define targets

.PHONY: clean distclean common code code_a petsc petsc_a $(SUBDIRS)

ifdef SAMR
    code: $(COMMON) $(SAMRAICODE) 
else
    code: $(COMMON) $(CODE) 
endif 

$(CODE):  $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	$(FC) $(FFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(CONTRIBLIBS) $(LIBS) $(LDFLAGS) 
	@mv $@ $(BINDIR)
	@cd $(BINDIR); ln -s -f $@ $(BASE).x

code_a: $(COMMON) $(A_CODE) 

$(A_CODE):  $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	$(FC) $(FFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(CONTRIBLIBS) $(LIBS) $(LDFLAGS)
	@mv $@ $(BINDIR)
	@cd $(BINDIR); ln -s -f $@ $(BASE)_a.x

petsc: $(COMMON) $(PETSCCODE) 

$(PETSCCODE): $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	${FLINKER} $(FFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(CONTRIBLIBS) ${PETSC_SNES_LIB} $(LIBS)
	@mv $@ $(BINDIR)
	@cd $(BINDIR); ln -s -f $@ $(BASE).petsc.x

petsc_a: $(COMMON) $(PETSC_A_CODE) 

$(PETSC_A_CODE): $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	${FLINKER} $(FFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(CONTRIBLIBS) ${PETSC_SNES_LIB} $(LIBS)
	@mv $@ $(BINDIR)
	@cd $(BINDIR); ln -s -f $@ $(BASE)_a.petsc.x

samrai: $(COMMON) $(SAMRAICODE)

$(SAMRAICODE): $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	${CXX} ${CXXFLAGS} $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(CONTRIBLIBS) $(VMECLIBS) ${PETSC_SNES_LIB} ${H5LIBS} ${LDFLAGS} ${LIBSAMRAI3D} ${LIBSAMRAI} ${LDLIBS} -o $@
	@mv $@ $(BINDIR)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) --no-print-directory -e -C $@ target

clean:
	-@rm -rf *.o *.mod *.dSYM

distclean: clean
	-@for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

setup: ;
	-@for file in $(LNK_FILES) ; do \
		ln -s -f $$file 2>/dev/null ; done

#Define dependencies

$(OBJMOD): $(COMMON_MODS)
$(OBJ):    $(MODS) $(COMMON_MODS)

#Define patterns

%.o : %.f
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $<

%.o : %.F
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $(CPPFLAGS) $<

%.o : %.f90
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $<

%.o : %.F90
	@echo 'Compiling' $@
	$(FC) -c $(MODDIRS) $(FFLAGS) $(CPPFLAGS) $<
