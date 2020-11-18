BASE         = pixie3d
EXEC         = $(BASE)-$(MACHINE)-$(REL1)-$(LABEL)_$(OPT)
EXECA        = $(BASE)_a-$(MACHINE)-$(REL1)-$(LABEL)_$(OPT)
CODE         = $(EXEC).x
A_CODE       = $(EXECA).x
ifdef ADIOS_VER
  PETSCCODE    = $(EXEC)-adios-$(ADIOS_VER).petsc.x
  PETSC_A_CODE = $(EXECA)-adios-$(ADIOS_VER).petsc.x
else
  PETSCCODE    = $(EXEC)-posix.petsc.x
  PETSC_A_CODE = $(EXECA)-posix.petsc.x
endif
SAMRAICODE   = $(EXEC).samrai.x

MODS := mhd_mod.F BC_mod.F PC_mod.F B_tools_mod.F

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/io\
	  $(COMMONDIR)/math\
	  $(COMMONDIR)/fft\
	  $(COMMONDIR)/OMP\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/eigen\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
	  $(COMMONDIR)/orbit

ifeq ($(MK),t)
  SUBDIRS += $(COMMONDIR)/mk
endif

SUBDIRS += $(COMMONDIR)/driver

ifdef SAMR
  SUBDIRS = $(COMMONDIR)/io\
	    $(COMMONDIR)/math\
  	    $(COMMONDIR)/OMP\
	    $(COMMONDIR)/interpolation/1d_int\
	    $(COMMONDIR)/grid\
	    $(COMMONDIR)/var_setup\
	    $(COMMONDIR)/driver_setup\
	    $(COMMONDIR)/matrixLightSolver\
	    $(COMMONDIR)/MG\
	    $(COMMONDIR)/orbit\
	    $(COMMONDIR)/driver-samrai
endif

ifdef BOPT 
ifdef PETSC_C
  SUBDIRS = $(COMMONDIR)/io\
  	    $(COMMONDIR)/math\
	    $(COMMONDIR)/interpolation/1d_int\
	    $(COMMONDIR)/OMP\
	    $(COMMONDIR)/grid\
	    $(COMMONDIR)/var_setup\
	    $(COMMONDIR)/driver_setup\
	    $(COMMONDIR)/driver-petsc\
	    $(COMMONDIR)/matrixLightSolver\
	    $(COMMONDIR)/MG\
	    $(COMMONDIR)/orbit

ifeq ($(MK),t)
  SUBDIRS += $(COMMONDIR)/mk
endif

SUBDIRS += $(COMMONDIR)/driver-petsc
endif
endif
