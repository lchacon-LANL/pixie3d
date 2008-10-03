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

      allocate (vcnv  (ilom:ihip,jlom:jhip,klom:khip,3)
     .         ,vcov  (ilom:ihip,jlom:jhip,klom:khip,3)
     .         ,vecnv (ilom:ihip,jlom:jhip,klom:khip,3)
     .         ,jcov_0(ilom:ihip,jlom:jhip,klom:khip,3))

#if defined(vec_pot)
      allocate (acnv   (ilom:ihip,jlom:jhip,klom:khip,3)
     .         ,acov   (ilom:ihip,jlom:jhip,klom:khip,3))

      allocate (bx     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,by     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,ax_cnv (ilom:ihip,jlom:jhip,klom:khip)
     .         ,ay_cnv (ilom:ihip,jlom:jhip,klom:khip)
     .         ,az_cnv (ilom:ihip,jlom:jhip,klom:khip))

      allocate (ax_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,ay_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,az_car (ilom:ihip,jlom:jhip,klom:khip))

#endif

      allocate (bcnv  (ilom:ihip,jlom:jhip,klom:khip,3)
     .         ,bcov  (ilom:ihip,jlom:jhip,klom:khip,3))

      allocate (bx_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,by_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz_cov (ilom:ihip,jlom:jhip,klom:khip))

      allocate (jx     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jy     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jz     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jx_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jy_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jz_cov (ilom:ihip,jlom:jhip,klom:khip))

      allocate (vx     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz     (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vx_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz_cov (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vex    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vey    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vez    (ilom:ihip,jlom:jhip,klom:khip))

      allocate (bx_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,by_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jx_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jy_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jz_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vx_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz_car (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vpar   (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vpsi   (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vperp  (ilom:ihip,jlom:jhip,klom:khip)
     .         ,divrgB (ilom:ihip,jlom:jhip,klom:khip)
     .         ,divrgV (ilom:ihip,jlom:jhip,klom:khip)
     .         ,divrgJ (ilom:ihip,jlom:jhip,klom:khip)
     .         ,Pflux  (ilom:ihip,jlom:jhip,klom:khip)
     .         ,p_tot  (ilom:ihip,jlom:jhip,klom:khip))

      allocate (qfactor(ilom:ihip,jlom:jhip,klom:khip)
     .         ,lambda (ilom:ihip,jlom:jhip,klom:khip)
     .         ,qfactr2(ilom:ihip,jlom:jhip,klom:khip)
     .         ,lambda2(ilom:ihip,jlom:jhip,klom:khip))

      allocate (eeta (ilom:ihip,jlom:jhip,klom:khip)
     .         ,nuu  (ilom:ihip,jlom:jhip,klom:khip))

c Setup grid aliases

      call defineGridAliases

c Allocate auxiliary variables

      call allocAuxVariables

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

      deallocate(vcnv,vcov,vecnv,jcov_0)

#if defined(vec_pot)
      deallocate (acnv)

      deallocate (bx,by,bz,ax_cnv,ay_cnv,az_cnv)

      deallocate (ax_car,ay_car,az_car)
#endif

      deallocate (bcnv,bcov)

      deallocate (bx_cov,by_cov,bz_cov)

      deallocate (jx,jy,jz,jx_cov,jy_cov,jz_cov)

      deallocate (vx,vy,vz,vx_cov,vy_cov,vz_cov,vex,vey,vez)

      deallocate (bx_car,by_car,bz_car,jx_car,jy_car,jz_car
     .           ,vx_car,vy_car,vz_car,vpar,vpsi,vperp 
     .           ,divrgB,divrgV,divrgJ,Pflux,p_tot )

      deallocate (qfactor,lambda,qfactr2,lambda2)

      deallocate (eeta,nuu )

      call deallocAuxVariables

c End programs

      end subroutine deallocateApplicationVariables
