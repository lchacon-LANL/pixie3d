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

      call findHeta

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

        allocate (bx    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,by    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,bz    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,bx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,by_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,bz_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,acnv(ilom:ihip,jlom:jhip,klom:khip,3)
     $           ,acov(ilom:ihip,jlom:jhip,klom:khip,3)
     $           ,vcnv(ilom:ihip,jlom:jhip,klom:khip,3))

        allocate (jx    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,jy    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,jz    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,jx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,jy_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,jz_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,jcov(ilom:ihip,jlom:jhip,klom:khip,3))

        allocate (vx    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,vy    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,vz    (ilom:ihip,jlom:jhip,klom:khip)
     .           ,vx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,vy_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,vz_cov(ilom:ihip,jlom:jhip,klom:khip)
     .           ,vex   (ilom:ihip,jlom:jhip,klom:khip)
     .           ,vey   (ilom:ihip,jlom:jhip,klom:khip)
     .           ,vez   (ilom:ihip,jlom:jhip,klom:khip))

        allocate (eeta  (ilom:ihip,jlom:jhip,klom:khip)
     .           ,nuu   (ilom:ihip,jlom:jhip,klom:khip))

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

c findHeta
c####################################################################
      subroutine findHeta
c--------------------------------------------------------------------
c     Sets aliases for grid quantities
c--------------------------------------------------------------------

      use grid

      use grid_aliases

      use transport_params

      use parameters

      implicit none

c Call variables

c Local variables

      real(8) :: kk,norm,idx,idy,idz
      integer(4) :: ig,jg,kg,i,j,k

c Begin program

      igx = 1
      igy = 1
      igz = 1

      kk = 0d0

      do k=1,nzd
        do j=1,nyd
          do i=1,nxd

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            !Maximum kk
            idx  = 1./dx(ig)
            if (nxd == 1) idx = 1d-2
            idy  = 1./dy(jg)
            if (nyd == 1) idy = 1d-2
            idz  = 1./dz(kg)
            if (nzd == 1) idz = 1d-2

            norm = vectorNorm(i,j,k,igx,igy,igz,idx,idy,idz,.true.)
            kk   = max(kk,norm)

          enddo
        enddo
      enddo

      heta = .2*di/kk

c End

      end subroutine findHeta
