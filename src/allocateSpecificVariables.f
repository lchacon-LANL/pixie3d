c allocateSpecificVariables
c######################################################################
      subroutine allocateSpecificVariables

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use parameters

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

c Begin program

        allocate (bx_cov(0:nxdp,0:nydp,0:nzdp)
     .           ,by_cov(0:nxdp,0:nydp,0:nzdp)
     .           ,bz_cov(0:nxdp,0:nydp,0:nzdp)
     .           ,divrgB(0:nxdp,0:nydp,0:nzdp))

        allocate (jx(0:nxdp,0:nydp,0:nzdp)
     .           ,jy(0:nxdp,0:nydp,0:nzdp)
     .           ,jz(0:nxdp,0:nydp,0:nzdp)
     .           ,jx_cov(0:nxdp,0:nydp,0:nzdp)
     .           ,jy_cov(0:nxdp,0:nydp,0:nzdp)
     .           ,jz_cov(0:nxdp,0:nydp,0:nzdp))

        allocate (vx(0:nxdp,0:nydp,0:nzdp)
     .           ,vy(0:nxdp,0:nydp,0:nzdp)
     .           ,vz(0:nxdp,0:nydp,0:nzdp))

        allocate (eeta(0:nxdp,0:nydp,0:nzdp)
     .           ,nuu (0:nxdp,0:nydp,0:nzdp))

c End program

      end subroutine allocateSpecificVariables
