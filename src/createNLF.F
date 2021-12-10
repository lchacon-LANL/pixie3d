c createNLF
c######################################################################
      subroutine createNLF

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use app_iosetup

      use problem_def

      use equilibrium

      use rw_bc, ONLY: rw_allocBCVariables

      implicit none

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

c Call variables

c Local variables

c Begin program

#if defined(vmec)
      if (trim(equil)=='vmec') call vmec_map(equ_file,gv%gparams)
#endif

      if (trim(equil)=='efit') call efit_map(equ_file,gv%gparams
     .                                      ,efit_coords)

      call rw_allocBCVariables
      
c     End programs

      end subroutine createNLF

c destroyNLF
c######################################################################
      subroutine destroyNLF

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use PC_var

      use rw_bc, ONLY: rw_deallocBCVariables

      implicit none

c Call variables

c Local variables

c Begin program

#if !defined(samrai)
      if (pc_vars_allcd) call deallocPCVariables
#endif

      call rw_deallocBCVariables
      
c End programs

      end subroutine destroyNLF
