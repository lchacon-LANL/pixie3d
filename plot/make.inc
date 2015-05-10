CODE         = pixplot.x
A_CODE       = pixplot_a.x
PETSCCODE    = pixplot.petsc.x
PETSC_A_CODE = pixplot_a.petsc.x

MODS := mhd_mod.F BC_mod.F

EXCLUDE += %ts_io_mod.F %ts_if_mod.F %ts_base_mod.F %ts_mod.F %proj_int_mod.F

override CPPFLAGS += -Dplotting

# Subdirectories

COMMONDIR = ../common

SUBDIRS = $(COMMONDIR)/math\
	  $(COMMONDIR)/io\
	  $(COMMONDIR)/grid\
	  $(COMMONDIR)/var_setup\
	  $(COMMONDIR)/driver_setup\
	  $(COMMONDIR)/interpolation/1d_int\
	  $(COMMONDIR)/matrixLightSolver\
	  $(COMMONDIR)/MG\
	  $(COMMONDIR)/orbit

ifeq ($(MK),t)
  SUBDIRS += $(COMMONDIR)/mk
  EXCLUDE = %mk_pc_mod.F
endif

SUBDIRS += $(COMMONDIR)/3dplot

# Set link files

LNK_FILES = ../src/applyBC.F\
	    ../src/BC_mod.F\
	    ../src/mhd_mod.F\
	    ../src/readInput.F

# Dependencies

#postProcessSolution.o: diagnostics.o
