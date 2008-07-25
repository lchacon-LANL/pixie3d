CODE         = pixie3d.x
A_CODE       = pixie3d_a.x
PETSCCODE    = pixie3d.petsc.x
PETSC_A_CODE = pixie3d_a.petsc.x

MODS := mhd_mod.F imposeBC_mod.F precond_mod.F

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/driver\
	  $(COMMONDIR)/math\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/eigen\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/interpolation/slatec/1d_spline\
	  $(COMMONDIR)/interpolation/slatec/2d_spline\
	  $(COMMONDIR)/interpolation/slatec/3d_spline\
	  $(COMMONDIR)/interpolation/slatec/spline_common\
	  $(COMMONDIR)/interpolation/slatec/xerror\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG

ifdef BOPT 
ifdef PETSC_C
  SUBDIRS = $(COMMONDIR)/driver-petsc\
	    $(COMMONDIR)/math\
	    $(COMMONDIR)/grid\
	    $(COMMONDIR)/interpolation/1d_int\
	    $(COMMONDIR)/interpolation/slatec/1d_spline\
	    $(COMMONDIR)/interpolation/slatec/2d_spline\
	    $(COMMONDIR)/interpolation/slatec/3d_spline\
	    $(COMMONDIR)/interpolation/slatec/spline_common\
	    $(COMMONDIR)/interpolation/slatec/xerror\
	    $(COMMONDIR)/matrixLightSolver\
	    $(COMMONDIR)/MG
endif
endif
