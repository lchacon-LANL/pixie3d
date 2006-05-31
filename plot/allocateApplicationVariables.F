c allocateApplicationVariables
c######################################################################
      subroutine allocateApplicationVariables

c----------------------------------------------------------------------
c     Allocates graphics-related arrays
c----------------------------------------------------------------------

      use parameters

      use grid_aliases

      use auxiliaryVariables

      implicit none

c Call variables

c Local variables

c Begin program

#if !defined(vec_pot)
      allocate (bx_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,by_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz_cov (ilom:ihip,jlom:jhip,klom:khip))
#else
      allocate (bx     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,by     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,bx_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,by_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,ax_cnv (ilom:ihip,jlom:jhip,klom:khip)
     .         ,ay_cnv (ilom:ihip,jlom:jhip,klom:khip)
     .         ,az_cnv (ilom:ihip,jlom:jhip,klom:khip)
     .         ,acnv   (ilom:ihip,jlom:jhip,klom:khip,3))
#endif

      allocate (jx     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jy     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jz     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jx_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jy_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jz_cov (ilom:ihip,jlom:jhip,klom:khip))

#if !defined(vec_pot)
      allocate (vx     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vx_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vex    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vey    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vez    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vcnv   (ilom:ihip,jlom:jhip,klom:khip,3))
#else
      allocate (vx     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vx_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,ax_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,ay_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,az_car (ilom:ihip,jlom:jhip,klom:khip))
#endif

      allocate (bx_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,by_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jx_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jy_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jz_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vx_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,divrgB (ilom:ihip,jlom:jhip,klom:khip)
     .         ,divrgV (ilom:ihip,jlom:jhip,klom:khip)
     .         ,divrgJ (ilom:ihip,jlom:jhip,klom:khip)
     .         ,Pflux  (ilom:ihip,jlom:jhip,klom:khip)
     .         ,p_tot  (ilom:ihip,jlom:jhip,klom:khip)
     .         ,qfactor(ilom:ihip,jlom:jhip,klom:khip)
     .         ,lambda (ilom:ihip,jlom:jhip,klom:khip)
     .         ,qfactr2(ilom:ihip,jlom:jhip,klom:khip)
     .         ,lambda2(ilom:ihip,jlom:jhip,klom:khip))

        allocate (eeta (ilom:ihip,jlom:jhip,klom:khip)
     .           ,nuu  (ilom:ihip,jlom:jhip,klom:khip))

      call defineGridAliases

c End programs

      end subroutine allocateApplicationVariables

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
