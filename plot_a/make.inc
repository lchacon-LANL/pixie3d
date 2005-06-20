PROD      = pixplot_a.x
BETA      = pixplot_a.beta.x
PETSCCODE = pixplot_a.petsc.x

MODS := mhd_mod.f imposeBC_mod.F

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/grid\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/interpolation/slatec/1d_spline\
	  $(COMMONDIR)/interpolation/slatec/2d_spline\
	  $(COMMONDIR)/interpolation/slatec/3d_spline\
	  $(COMMONDIR)/interpolation/slatec/spline_common\
	  $(COMMONDIR)/interpolation/slatec/xerror\
	  $(COMMONDIR)/3dplot

