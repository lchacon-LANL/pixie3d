CODE         = pixplot.x
A_CODE       = pixplot_a.x
PETSCCODE    = pixplot.petsc.x
PETSC_A_CODE = pixplot_a.petsc.x

MODS := mhd_mod.F imposeBC_mod.F

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/interpolation/slatec/1d_spline\
	  $(COMMONDIR)/interpolation/slatec/2d_spline\
	  $(COMMONDIR)/interpolation/slatec/3d_spline\
	  $(COMMONDIR)/interpolation/slatec/spline_common\
	  $(COMMONDIR)/interpolation/slatec/xerror\
	  $(COMMONDIR)/3dplot

LNK_FILES = applyBoundaryConditions.F\
	    imposeBC_mod.F\
	    mhd_mod.F\
	    readInput.F             
