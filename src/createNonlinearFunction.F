c createNonlinearFunction
c######################################################################
      subroutine createNonlinearFunction

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use nlfunction_setup

      use precond_variables

      implicit none

c Call variables

c Local variables

c Begin program

c Auxiliary variables

#if !defined(samrai)
      call allocAuxVariables
#endif

c Preconditioner

      call allocPrecVariables

c Grid aliases

      if (.not.associated(xx)) then
        xx => grid_params%xx
        yy => grid_params%yy
        zz => grid_params%zz

        dxh => grid_params%dxh
        dyh => grid_params%dyh
        dzh => grid_params%dzh

        dx => grid_params%dx
        dy => grid_params%dy
        dz => grid_params%dz
      endif

c End programs

      end subroutine createNonlinearFunction

c destroyNonlinearFunction
c######################################################################
      subroutine destroyNonlinearFunction

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use nlfunction_setup

      use precond_variables

      implicit none

c Call variables

c Local variables

c Begin program

#if !defined(samrai)
      call deallocAuxVariables
#endif

      call deallocPrecVariables

c End programs

      end subroutine destroyNonlinearFunction
