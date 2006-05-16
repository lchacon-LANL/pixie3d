PROD      = pixplot.x
BETA      = pixplot.beta.x
PETSCCODE = pixplot.petsc.x

MODS := mhd_mod.f imposeBC_mod.F

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
	    mhd_mod.f\
	    readInput.f              
