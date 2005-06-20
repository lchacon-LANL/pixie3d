PROD      = pixie3d_a.x
BETA      = pixie3d_a.beta.x
PETSCCODE = pixie3d_a.petsc.x

MODS := mhd_mod.f imposeBC_mod.F precond_mod.F

COMMONDIR = ../common

ifdef BOPT
  SUBDIRS = $(COMMONDIR)/driver-petsc\
	    $(COMMONDIR)/grid\
	    $(COMMONDIR)/interpolation/1d_int\
	    $(COMMONDIR)/interpolation/slatec/1d_spline\
	    $(COMMONDIR)/interpolation/slatec/2d_spline\
	    $(COMMONDIR)/interpolation/slatec/3d_spline\
	    $(COMMONDIR)/interpolation/slatec/spline_common\
	    $(COMMONDIR)/interpolation/slatec/xerror\
	    $(COMMONDIR)/matrixLightSolver\
	    $(COMMONDIR)/MG
else
  SUBDIRS = $(COMMONDIR)/driver\
	    $(COMMONDIR)/grid\
	    $(COMMONDIR)/newtonGMRES\
	    $(COMMONDIR)/interpolation/1d_int\
	    $(COMMONDIR)/interpolation/slatec/1d_spline\
	    $(COMMONDIR)/interpolation/slatec/2d_spline\
	    $(COMMONDIR)/interpolation/slatec/3d_spline\
	    $(COMMONDIR)/interpolation/slatec/spline_common\
	    $(COMMONDIR)/interpolation/slatec/xerror\
	    $(COMMONDIR)/matrixLightSolver\
	    $(COMMONDIR)/MG
endif