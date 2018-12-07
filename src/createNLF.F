c createNLF
c######################################################################
      subroutine createNLF

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use app_iosetup

      use problem_def

      implicit none

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

c Call variables

c Local variables

c Begin program

#if defined(vmec)
      if (trim(equil) == 'vmec') call vmec_map(equ_file,gv%gparams)
#endif

      if (trim(equil) == 'efit') call efit_map(equ_file,gv%gparams)

c     End programs

      end subroutine createNLF

c destroyNLF
c######################################################################
      subroutine destroyNLF

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

cc      use nlfunction_setup

      use PC_var

      implicit none

c Call variables

c Local variables

c Begin program

#if !defined(samrai)
      if (pc_vars_allcd) call deallocPCVariables
#endif

c End programs

      end subroutine destroyNLF
