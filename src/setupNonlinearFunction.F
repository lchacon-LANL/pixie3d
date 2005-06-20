c setupNonlinearFunction
c#################################################################
      subroutine setupNonlinearFunction(varray)
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use parameters

      use variables

      use nlfunction_setup

      use equilibrium

      use constants

      use timeStepping

cc      use debug

      use error

      implicit none

c Call variables

      type (var_array),target :: varray

c Local variables

      integer(4) :: i,j,k,ig,jg,kg

      integer(4) :: dim,loc,ibc,ieq,bctype

      real(8)    :: dh,vmax

c Begin program

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Impose boundary conditions and find auxiliary quantities

      call imposeBoundaryConditions(varray,igx,igy,igz)

c Set aliases

      rho => varray%array_var(IRHO)%array
      rvx => varray%array_var(IVX )%array
      rvy => varray%array_var(IVY )%array
      rvz => varray%array_var(IVZ )%array
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
      tmp => varray%array_var(ITMP)%array

c Safeguards

      if (minval(rho) < 0d0) then
        call pstop('setupNonlinearFunction',
     .             'Warning: negative densities are occurring')
      endif

      if (minval(tmp) < 0d0) then
        call pstop('setupNonlinearFunction',
     .             'Warning: negative temperatures are occurring')
      endif

c Transport parameters

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            nuu  (i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
            eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

      ejx = jx_cov*eeta
      ejy = jy_cov*eeta
      ejz = jz_cov*eeta

      bcnv(:,:,:,1) = bx
      bcnv(:,:,:,2) = by
      bcnv(:,:,:,3) = bz

      vcnv(:,:,:,1) = vx
      vcnv(:,:,:,2) = vy
      vcnv(:,:,:,3) = vz

c End program

      end subroutine setupNonlinearFunction

c killNonlinearFunction
c#################################################################
      subroutine killNonlinearFunction
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

c Begin program

c End program

      end subroutine killNonlinearFunction
