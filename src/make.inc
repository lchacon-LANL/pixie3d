CODE         = pixie3d.x
A_CODE       = pixie3d_a.x
PETSCCODE    = pixie3d.petsc.x
PETSC_A_CODE = pixie3d_a.petsc.x
SAMRAICODE   = pixie3d.samrai.x

MODS := mhd_mod.F imposeBC_mod.F precond_mod.F B_tools_mod.F par_int_mod.F

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/io\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/eigen\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
	  $(COMMONDIR)/orbit

ifeq ($(MK),t)
  SUBDIRS += $(COMMONDIR)/mk
endif

SUBDIRS += $(COMMONDIR)/driver

ifdef SAMR
  SUBDIRS = $(COMMONDIR)/math\
  	    $(COMMONDIR)/io\
	    $(COMMONDIR)/grid\
	    $(COMMONDIR)/var_setup\
	    $(COMMONDIR)/driver_setup\
	    $(COMMONDIR)/interpolation/1d_int\
	    $(COMMONDIR)/matrixLightSolver\
	    $(COMMONDIR)/MG\
	    $(COMMONDIR)/orbit\
	    $(COMMONDIR)/driver-samrai
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
	    $(COMMONDIR)/matrixLightSolver\
	    $(COMMONDIR)/MG\
	    $(COMMONDIR)/orbit

ifeq ($(MK),t)
  SUBDIRS += $(COMMONDIR)/mk
endif

SUBDIRS += $(COMMONDIR)/driver-petsc
endif
endif
