BASE         = pixeq_xfer
EXEC         = $(BASE)-$(MACHINE)-$(REL1)-$(LABEL)_$(OPT)
CODE         = $(EXEC)_xfer.x
ifdef ADIOS_VER
  PETSCCODE    = $(EXEC)-adios-$(ADIOS_VER).petsc.x
else
  PETSCCODE    = $(EXEC)-posix.petsc.x
endif
SAMRAICODE   = $(EXEC).samrai.x

#Define message

MESSAGE = "Compiling EQUILIBRIUM  modules"

override CPPFLAGS += -D_equ

#Modules

MODS := io_mod.F90 mhd_mod.F BC_mod.F PC_mod.F rw_bim_mod.F90 rw_mod.F

#Define subdirectories for common code

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/io\
	  $(COMMONDIR)/math\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/fft1d\
	  $(COMMONDIR)/OMP\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/newtonGMRES\
	  $(COMMONDIR)/eigen\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/orbit
	  
# Set link files

LNK_FILES = ../src/applyBC.F\
	    ../src/BC_mod.F\
	    ../src/mhd_mod.F\
	    ../src/PC_mod.F\
	    ../src/rw_mod.F\
	    ../src/io_mod.F90\
	    ../src/readInput.F

#Exclusions

EXCLUDE += %ts_io_mod.F %ts_if_mod.F %ts_base_mod.F %ts_mod.F %proj_int_mod.F %nk_mod.F %eigen.F
