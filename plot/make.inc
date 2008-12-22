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
	  $(COMMONDIR)/interpolation/slatec/1d_spline\
	  $(COMMONDIR)/interpolation/slatec/2d_spline\
	  $(COMMONDIR)/interpolation/slatec/3d_spline\
	  $(COMMONDIR)/interpolation/slatec/spline_common\
	  $(COMMONDIR)/interpolation/slatec/xerror\
	  $(COMMONDIR)/3dplot

LNK_FILES = ../src/applyBoundaryConditions.F\
	    ../src/imposeBC_mod.F\
	    ../src/mhd_mod.F\
	    ../src/readInput.F\
            ../src/vmec.f90 
