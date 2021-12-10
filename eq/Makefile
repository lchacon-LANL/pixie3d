# Read application specific variables

-include make.inc

# GENERAL PURPOSE MAKEFILE

SRC  ?= $(filter-out $(EXCLUDE),$(wildcard *.[F,c,f] *.f90 *.F90))

MODS ?= $(filter-out $(EXCLUDE),$(wildcard *_mod.[F,f] *_mod.f90 *_mod.F90))

prefix = .

# Objects

OBJS   := $(filter %.o, $(patsubst %.f,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.c,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.C,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.F,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.f90,%.o,$(filter-out $(MODS),$(SRC)))\
                        $(patsubst %.F90,%.o,$(filter-out $(MODS),$(SRC))))

OBJMOD := $(patsubst %.f,%.o,$(patsubst %.F,%.f,$(patsubst %.f90,%.f,$(patsubst %.F90,%.f,$(MODS)))))

COMMON_MODS = $(foreach dir,$(SUBDIRS),$(filter-out $(EXCLUDE),$(wildcard $(dir)/*_mod.[f,F] $(dir)/*_mod.f90 $(dir)/*_mod.F90)))
COMMON_SRC  = $(foreach dir,$(SUBDIRS),$(filter-out $(EXCLUDE) %test.f,$(wildcard $(dir)/*.[f,F,c,C] $(dir)/*.f90 $(dir)/*.F90)))
COMMON_OBJS = $(filter %.o, $(patsubst %.f,%.o,$(COMMON_SRC))\
                            $(patsubst %.c,%.o,$(COMMON_SRC))\
                            $(patsubst %.C,%.o,$(COMMON_SRC))\
                            $(patsubst %.F,%.o,$(COMMON_SRC))\
                            $(patsubst %.f90,%.o,$(COMMON_SRC))\
                            $(patsubst %.F90,%.o,$(COMMON_SRC)))

#Module search path

MODDIRS = $(MODPATH) $(patsubst $(COMMONDIR)%,$(ADDMODFLAG)$(COMMONDIR)%,$(SUBDIRS))

#Define targets

.PHONY: clean distclean common code code_a petsc petsc_a $(SUBDIRS)

code code_a: common $(CODE)

petsc petsc_a: common $(PETSCCODE) 

$(CODE): $(OBJMOD) $(OBJS) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	$(FC) $(FFLAGS) -o $@ $(OBJMOD) $(OBJS) $(COMMON_OBJS) $(CONTRIBLIBS) $(LIBS) $(LDFLAGS)
	@mv $@ $(BINDIR)
	@cd $(BINDIR); ln -s -f $@ $(BASE).x

$(PETSCCODE): $(OBJMOD) $(OBJS) $(COMMON_OBJS)
	@echo ''
	@echo 'Linking ' $@ ' code'
	@echo ''
	${FLINKER} $(FFLAGS) -o $@ $(OBJMOD) $(OBJS) $(COMMON_OBJS) $(CONTRIBLIBS) ${PETSC_SNES_LIB} $(LIBS) $(LDFLAGS)
	@mv $@ $(BINDIR)
	@cd $(BINDIR); ln -s -f $@ $(BASE).petsc.x

target: common message $(OBJMOD) $(OBJS) $(COMMON_OBJS)

message: ;
	@echo ''
#	@echo $(MESSAGE)

common: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -e -C $@ target

#Cleaning targets

clean: ;
	-@rm -f *.o *.*mod *.a *.x

distclean: clean
	-@for subdir in $(SUBDIRS) ; do \
		$(MAKE) -C $$subdir clean;  done

#Main setup targets

setup: ;
	-@for file in $(LNK_FILES) ; do \
		ln -s $$file 2>/dev/null ; done

#Define dependencies

$(OBJMOD): $(COMMON_MODS)
$(OBJS) : $(MODS) $(COMMOM_MODS)

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

%.o : %.c
	@echo 'Compiling' $@
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $<

%.o : %.C
	@echo 'Compiling' $@
	$(CXX) -c $(CXXFLAGS) $(CPPFLAGS) $<


