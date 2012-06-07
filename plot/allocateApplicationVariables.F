c allocateApplicationVariables
c######################################################################
      subroutine allocateApplicationVariables

c----------------------------------------------------------------------
c     Allocates graphics-related arrays
c----------------------------------------------------------------------

      use parameters

      use auxiliaryVariables

      use auxPlotVariables

      implicit none

c Call variables

c Local variables

c Begin program

      call allocAuxPlotVar

c Allocate auxiliary variables

      call defineAuxVariables(gv%aux)

c End programs

      end subroutine allocateApplicationVariables

c deallocateApplicationVariables
c######################################################################
      subroutine deallocateApplicationVariables

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use parameters

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
