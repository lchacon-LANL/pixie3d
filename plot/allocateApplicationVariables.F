c allocateApplicationVariables
c######################################################################
      subroutine allocateApplicationVariables

c----------------------------------------------------------------------
c     Allocates graphics-related arrays
c----------------------------------------------------------------------

      use app_iosetup

      use auxiliaryVariables

      use auxPlotVariables

      implicit none

c Call variables

c Local variables

      INTERFACE
         subroutine vmec_map(equ_file,g_def)
         use grid
         type(grid_mg_def),pointer :: g_def
         character(*) :: equ_file
         end subroutine vmec_map

         subroutine efit_map(equ_file,g_def)
           use grid
           type(grid_mg_def),pointer :: g_def
           character(*) :: equ_file
         end subroutine efit_map

      END INTERFACE

c Begin program

      call allocAuxPlotVar

c Allocate auxiliary variables

      call defineAuxVariables(gv%aux)

c Read coordinate map
      
#if defined(vmec)
      if (trim(equil) == 'vmec') call vmec_map(equ_file,gv%gparams)
#endif

      if (trim(equil) == 'efit') call efit_map(equ_file,gv%gparams)

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

      implicit none

c Call variables

c Local variables

c Begin program

      call deallocLocalAuxVars

      call deallocAuxPlotVar

c End programs

      end subroutine deallocateApplicationVariables
