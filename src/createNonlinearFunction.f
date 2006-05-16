c createNonlinearFunction
c######################################################################
      subroutine createNonlinearFunction

c----------------------------------------------------------------------
c     Allocates application-related arrays
c----------------------------------------------------------------------

      use grid_aliases

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

      use precond_variables

      implicit none

c Call variables

c Local variables

c Begin program

        allocate (bx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,by_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,bz_cov(ilom:ihip,jlom:jhip,klom:khip)
     $           ,bcnv  (ilom:ihip,jlom:jhip,klom:khip,3)
     $           ,vcnv  (ilom:ihip,jlom:jhip,klom:khip,3)
     $           ,vecnv (ilom:ihip,jlom:jhip,klom:khip,3))
cc     $           ,vcov  (ilom:ihip,jlom:jhip,klom:khip,3))

        allocate (jx    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,jy    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,jz    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,jx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,jy_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,jz_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,ejx   (ilom:ihip,jlom:jhip,klom:khip)
     .           ,ejy   (ilom:ihip,jlom:jhip,klom:khip)
     .           ,ejz   (ilom:ihip,jlom:jhip,klom:khip))

        allocate (vx    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,vy    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,vz    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,vx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,vy_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,vz_cov(ilom:ihip,jlom:jhip,klom:khip))

        allocate (eeta  (ilom:ihip,jlom:jhip,klom:khip)
     .           ,nuu   (ilom:ihip,jlom:jhip,klom:khip))

c SI operator
        allocate(b_n(ilom:ihip,jlom:jhip,klom:khip,3)
     .          ,p_n(ilom:ihip,jlom:jhip,klom:khip,1)
     .          ,dv_dt (ilom:ihip,jlom:jhip,klom:khip,3))
c SI operator

        call allocPrecVariables

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
