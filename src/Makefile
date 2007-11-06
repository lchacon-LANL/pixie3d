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

ifdef VMEC
  EXCLUDE = %test.f
else
  EXCLUDE = %test.f %vmec.f90
endif

# Source files

ifndef SRC
  SRC := $(filter-out $(EXCLUDE),$(wildcard *.[F,c,f] *.f90))
endif

ifndef MODS
  MODS := $(filter-out $(EXCLUDE),$(wildcard *_mod.[F,f] *_mod.f90))
endif

prefix = .

# Petsc include

DRIVERDIR = driver
ifdef BOPT
  include ${PETSC_DIR}/bmake/common/base
  SRC2      = $(SRC)
  ifdef PETSC_C
    EXCLUDE  += %explicit.f
  endif
  SRC      := $(filter-out $(EXCLUDE),$(SRC2))
endif

# Objects

OBJ    := $(filter %.o, $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.c,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.F,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.f90,%.o,$(filter-out $(MODS),$(SRC))))


OBJMOD := $(patsubst %.f,%.o,$(patsubst %.F,%.f,$(patsubst %.f90,%.f,$(MODS))))

COMMON_MODS = $(foreach dir,$(SUBDIRS),$(wildcard $(dir)/*_mod.[f,F] $(dir)/*_mod.f90))
COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(EXCLUDE),$(wildcard $(dir)/*.[f,F,c] $(dir)/*.f90)))
COMMON_OBJS = $(filter %.o, $(patsubst %.f,%.o,$(COMMON_SRC))\
                            $(patsubst %.c,%.o,$(COMMON_SRC))\
                            $(patsubst %.F,%.o,$(COMMON_SRC))\
                            $(patsubst %.f90,%.o,$(COMMON_SRC)))

#Define linker flags

#LDFLAGS = 

#Module search path

MODDIRS = $(MODPATH) $(patsubst $(COMMONDIR)%,$(ADDMODFLAG)$(COMMONDIR)%,$(SUBDIRS))

#Define targets

.PHONY: clean distclean common code code_a petsc petsc_a $(SUBDIRS)

code: common $(CODE)

$(CODE):  $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	$(FC) $(FFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(LIBS) $(VMECLIBS) $(H5LIBS) $(LDFLAGS)
	@mv $@ $(BINDIR)

code_a: common $(A_CODE)

$(A_CODE):  $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	$(FC) $(FFLAGS) -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(LIBS) $(VMECLIBS) $(H5LIBS) $(LDFLAGS)
	@mv $@ $(BINDIR)

petsc: common $(PETSCCODE)

$(PETSCCODE): $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	${FLINKER} -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(VMECLIBS) ${PETSC_SNES_LIB} $(H5LIBS)
	@mv $@ $(BINDIR)

petsc_a: common $(PETSC_A_CODE)

$(PETSC_A_CODE): $(OBJMOD) $(OBJ) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	${FLINKER} -o $@ $(COMMON_OBJS) $(OBJ)  $(OBJMOD) $(VMECLIBS) ${PETSC_SNES_LIB} $(H5LIBS)
	@mv $@ $(BINDIR)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -e -C $@ target

clean:
	-rm -f *.o *.mod

distclean: clean
	-for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

setup: ;
	-for file in $(LNK_FILES) ; do \
		ln -s $$file 2>/dev/null ; done

#Define dependencies

$(OBJMOD): $(COMMON_MODS)
$(OBJ):    $(MODS) $(COMMON_MODS)
#$(OBJ):  $(MODS) $(COMMON_MODS) $(COMMONDIR)/BC/*_mod.[f,F]
#imposeBC_mod.o: $(COMMONDIR)/BC/*_mod.[f,F] 

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
