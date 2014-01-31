c createNLF
c######################################################################
      subroutine createNLF

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use nlfunction_setup

      use PC_var

      implicit none

c Call variables

c Local variables

c Begin program

c Auxiliary variables

#if !defined(samrai)
cc      call defineAuxVariables(gv%aux)

c Preconditioner

      call allocPCVariables
#endif

c     End programs

      end subroutine createNLF

c destroyNLF
c######################################################################
      subroutine destroyNLF

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use nlfunction_setup

      use PC_var

      implicit none

c Call variables

c Local variables

c Begin program

#if !defined(samrai)
cc      call deallocLocalAuxVars

      call deallocPCVariables
#endif

c End programs

      end subroutine destroyNLF
