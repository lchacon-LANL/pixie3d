c allocateApplicationVariables
c######################################################################
      subroutine allocateApplicationVariables

c----------------------------------------------------------------------
c     Allocates graphics-related arrays
c----------------------------------------------------------------------

      use parameters

      use grid_aliases

      use auxiliaryVariables

      use auxPlotVariables

      implicit none

c Call variables

c Local variables

c Begin program

      call allocAuxPlotVar

c Allocate auxiliary variables

      call allocAuxVariables(gv%aux)

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

#if !defined(samrai)
      call deallocAuxVariables(gv%aux)
#endif

      call deallocAuxPlotVar

c End programs

      end subroutine deallocateApplicationVariables
