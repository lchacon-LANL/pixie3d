c allocateApplicationVariables
c######################################################################
      subroutine allocateApplicationVariables

c----------------------------------------------------------------------
c     Allocates graphics-related arrays
c----------------------------------------------------------------------

      use app_iosetup

      use auxiliaryVariables

      use auxPlotVariables

      use equilibrium

      use local_BCS_variables, ONLY: allocBCVariables

      implicit none

c Call variables

c Local variables

      INTERFACE
        subroutine vmec_map(equ_file,g_def)
          use grid
          type(grid_mg_def),pointer :: g_def
          character(*) :: equ_file
        end subroutine vmec_map

        subroutine efit_map(equ_file,g_def,coord)
          use grid
          type(grid_mg_def),pointer :: g_def
          character(*) :: equ_file,coord
        end subroutine efit_map
      END INTERFACE

c Begin program

      call allocAuxPlotVar

c Allocate auxiliary variables

      call defineAuxVariables(gv%aux)

c Read coordinate map
      
#if defined(vmec)
      if (trim(equil)=='vmec') call vmec_map(equ_file,gv%gparams)
#endif

      if (trim(equil)=='efit') call efit_map(equ_file,gv%gparams
     .                                      ,efit_coords)

c Allocate BC variables
      
      call allocBCVariables
 
c End programs

      end subroutine allocateApplicationVariables

c deallocateApplicationVariables
c######################################################################
      subroutine deallocateApplicationVariables

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use auxiliaryVariables

      use auxPlotVariables

      use local_BCS_variables, ONLY: bc_vars_allcd,deallocBCVariables

      implicit none

c Call variables

c Local variables

c Begin program

      call deallocLocalAuxVars

      call deallocAuxPlotVar

      if (bc_vars_allcd) call deallocBCVariables

c End programs

      end subroutine deallocateApplicationVariables
