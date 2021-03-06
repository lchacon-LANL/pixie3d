BASE         = pixplot
EXEC         = $(BASE)-$(MACHINE)-$(REL1)-$(LABEL)_$(OPT)
EXECA        = $(BASE)_a-$(MACHINE)-$(REL1)-$(LABEL)_$(OPT)
CODE         = $(EXEC).x
A_CODE       = $(EXECA).x
ifdef ADIOS_VER
  PETSCCODE    = $(EXEC)-adios-$(ADIOS_VER).petsc.x
  PETSC_A_CODE = $(EXECA)-adios-$(ADIOS_VER).petsc.x
else
  PETSCCODE    = $(EXEC)-posix.petsc.x
  PETSC_A_CODE = $(EXECA)-posix.petsc.x
endif

MODS := io_mod.F90 mhd_mod.F BC_mod.F rw_mod.F PC_mod.F

EXCLUDE += %ts_if_mod.F %ts_mod.F %proj_int_mod.F

override CPPFLAGS += -Dplotting

# Subdirectories

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/io\
	  $(COMMONDIR)/math\
	  $(COMMONDIR)/fft\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
	  $(COMMONDIR)/OMP\
	  $(COMMONDIR)/orbit

ifeq ($(MK),t)
  SUBDIRS += $(COMMONDIR)/mk
  EXCLUDE += %mk_pc_mod.F
endif

SUBDIRS += $(COMMONDIR)/3dplot

# Set link files

LNK_FILES = ../src/applyBC.F\
	    ../src/BC_mod.F\
	    ../src/mhd_mod.F\
	    ../src/PC_mod.F\
	    ../src/readInput.F\
	    ../src/vmec.f90\
	    ../src/efit.F90\
	    ../src/io_mod.F90\
	    ../src/rw_mod.F

# Dependencies

postProcessSolution.o: rfx_module
