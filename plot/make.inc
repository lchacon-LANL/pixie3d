CODE         = pixplot.x
A_CODE       = pixplot_a.x
PETSCCODE    = pixplot.petsc.x
PETSC_A_CODE = pixplot_a.petsc.x

MODS := mhd_mod.F imposeBC_mod.F

# Subdirectories

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/io\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
	  $(COMMONDIR)/orbit

ifeq ($(MK),t)
  SUBDIRS += $(COMMONDIR)/mk
  EXCLUDE = %mk_pc_mod.F
endif

SUBDIRS += $(COMMONDIR)/3dplot

# Set link files

LNK_FILES = ../src/applyBoundaryConditions.F\
	    ../src/imposeBC_mod.F\
	    ../src/mhd_mod.F\
	    ../src/readInput.F

# Dependencies

#postProcessSolution.o: diagnostics.o
