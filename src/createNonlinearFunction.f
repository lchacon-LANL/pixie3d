c createNonlinearFunction
c######################################################################
      subroutine createNonlinearFunction

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      implicit none

c Call variables

c Local variables

c Begin program

      call allocateSpecificVariables

      call defineGridAliases

c End programs

      end subroutine createNonlinearFunction

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
     .           ,bz_cov(0:nxdp,0:nydp,0:nzdp))

        allocate (jx    (0:nxdp,0:nydp,0:nzdp)
     .           ,jy    (0:nxdp,0:nydp,0:nzdp)
     .           ,jz    (0:nxdp,0:nydp,0:nzdp)
     .           ,jx_cov(0:nxdp,0:nydp,0:nzdp)
     .           ,jy_cov(0:nxdp,0:nydp,0:nzdp)
     .           ,jz_cov(0:nxdp,0:nydp,0:nzdp))

        allocate (vx    (0:nxdp,0:nydp,0:nzdp)
     .           ,vy    (0:nxdp,0:nydp,0:nzdp)
     .           ,vz    (0:nxdp,0:nydp,0:nzdp)
     .           ,vx_cov(0:nxdp,0:nydp,0:nzdp)
     .           ,vy_cov(0:nxdp,0:nydp,0:nzdp)
     .           ,vz_cov(0:nxdp,0:nydp,0:nzdp))

        allocate (eeta  (0:nxdp,0:nydp,0:nzdp)
     .           ,nuu   (0:nxdp,0:nydp,0:nzdp))

c End program

      end subroutine allocateSpecificVariables

c defineGridAliases
c####################################################################
      subroutine defineGridAliases
c--------------------------------------------------------------------
c     Sets aliases for grid quantities
c--------------------------------------------------------------------

      use grid

      use grid_aliases

      implicit none

c Call variables

c Local variables

c Begin program

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

c End

      end subroutine defineGridAliases
