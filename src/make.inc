CODE         = pixie3d.x
A_CODE       = pixie3d_a.x
PETSCCODE    = pixie3d.petsc.x
PETSC_A_CODE = pixie3d_a.petsc.x
SAMRAICODE   = pixie3d.samrai.x

MODS := mhd_mod.F precond_mod.F imposeBC_mod.F 

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/io\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/driver\
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/eigen\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/interpolation/slatec/1d_spline\
	  $(COMMONDIR)/interpolation/slatec/2d_spline\
	  $(COMMONDIR)/interpolation/slatec/3d_spline\
	  $(COMMONDIR)/interpolation/slatec/spline_common\
	  $(COMMONDIR)/interpolation/slatec/xerror\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
	  $(COMMONDIR)/orbit\
	  $(COMMONDIR)/lsode

ifdef SAMR
  SUBDIRS = $(COMMONDIR)/math\
  	    $(COMMONDIR)/io\
	    $(COMMONDIR)/grid\
	    $(COMMONDIR)/var_setup\
	    $(COMMONDIR)/driver_setup\
	    $(COMMONDIR)/driver-samrai\
	    $(COMMONDIR)/interpolation/1d_int\
	    $(COMMONDIR)/interpolation/slatec/1d_spline\
	    $(COMMONDIR)/interpolation/slatec/2d_spline\
	    $(COMMONDIR)/interpolation/slatec/3d_spline\
	    $(COMMONDIR)/interpolation/slatec/spline_common\
	    $(COMMONDIR)/interpolation/slatec/xerror\
	    $(COMMONDIR)/matrixLightSolver\
	    $(COMMONDIR)/MG\
	    $(COMMONDIR)/orbit\
	    $(COMMONDIR)/lsode
endif

ifdef BOPT 
ifdef PETSC_C
  SUBDIRS = $(COMMONDIR)/math\
	    $(COMMONDIR)/io\
	    $(COMMONDIR)/grid\
	    $(COMMONDIR)/var_setup\
	    $(COMMONDIR)/driver_setup\
	    $(COMMONDIR)/driver-petsc\
	    $(COMMONDIR)/interpolation/1d_int\
	    $(COMMONDIR)/interpolation/slatec/1d_spline\
	    $(COMMONDIR)/interpolation/slatec/2d_spline\
	    $(COMMONDIR)/interpolation/slatec/3d_spline\
	    $(COMMONDIR)/interpolation/slatec/spline_common\
	    $(COMMONDIR)/interpolation/slatec/xerror\
	    $(COMMONDIR)/matrixLightSolver\
	    $(COMMONDIR)/MG\
	    $(COMMONDIR)/orbit
	    $(COMMONDIR)/lsode
endif
endif