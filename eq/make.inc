#Define message

MESSAGE = "Compiling EQUILIBRIUM  modules"

override CPPFLAGS += -D_equ

#Modules

MODS := mhd_mod.F BC_mod.F

#Define subdirectories for common code

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/math\
          $(COMMONDIR)/io\
          $(COMMONDIR)/interpolation/1d_int\
          $(COMMONDIR)/grid\
	  $(COMMONDIR)/driver_setup\
          $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG

# Set link files

LNK_FILES = ../src/applyBC.F\
	    ../src/BC_mod.F\
	    ../src/mhd_mod.F\
	    ../src/readInput.F\
	    ../src/applyBC.F

#Exclusions

EXCLUDE += %ts_io_mod.F %ts_if_mod.F %ts_base_mod.F %ts_mod.F %proj_int_mod.F
