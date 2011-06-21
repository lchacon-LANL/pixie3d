CODE         = pixplot.x
A_CODE       = pixplot_a.x
PETSCCODE    = pixplot.petsc.x
PETSC_A_CODE = pixplot_a.petsc.x

MODS := mhd_mod.F imposeBC_mod.F

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/io\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
	  $(COMMONDIR)/orbit\
	  $(COMMONDIR)/mk\
	  $(COMMONDIR)/3dplot

LNK_FILES = ../src/applyBoundaryConditions.F\
	    ../src/imposeBC_mod.F\
	    ../src/mhd_mod.F\
	    ../src/readInput.F

EXCLUDE = %mk_pc_mod.F