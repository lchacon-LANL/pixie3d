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
