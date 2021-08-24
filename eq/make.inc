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

MODS := mhd_mod.F BC_mod.F PC_mod.F

#Define subdirectories for common code

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/math\
          $(COMMONDIR)/io\
          $(COMMONDIR)/fft\
          $(COMMONDIR)/interpolation/1d_int\
          $(COMMONDIR)/grid\
	  $(COMMONDIR)/driver_setup\
          $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
	  $(COMMONDIR)/OMP\
	  $(COMMONDIR)/orbit

# Set link files

LNK_FILES = ../src/applyBC.F\
	    ../src/BC_mod.F\
	    ../src/mhd_mod.F\
	    ../src/PC_mod.F\
	    ../src/readInput.F

#Exclusions

EXCLUDE += %ts_io_mod.F %ts_if_mod.F %ts_base_mod.F %ts_mod.F %proj_int_mod.F
