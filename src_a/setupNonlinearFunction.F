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
      ax  => varray%array_var(IAX )%array
      ay  => varray%array_var(IAY )%array
      az  => varray%array_var(IAZ )%array
      tmp => varray%array_var(ITMP)%array

      acov(:,:,:,1) = ax
      acov(:,:,:,2) = ay
      acov(:,:,:,3) = az

      vcnv(:,:,:,1) = vx
      vcnv(:,:,:,2) = vy
      vcnv(:,:,:,3) = vz

c Safeguards

      if (minval(rho) < 0d0) then
        call pstop('setupNonlinearFunction',
     .             'Warning: negative densities are occurring')
      endif

      if (minval(tmp) < 0d0) then
        call pstop('setupNonlinearFunction',
     .             'Warning: negative temperatures are occurring')
      endif

ccc Current
cc
cc      do k=0,nz+1
cc        do j=0,ny+1
cc          do i=0,nx+1
cc            jx(i,j,k)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                     ,bx_cov,by_cov,bz_cov,1)
cc            jy(i,j,k)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                     ,bx_cov,by_cov,bz_cov,2)
cc            jz(i,j,k)=curl2(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                     ,bx_cov,by_cov,bz_cov,3)
cc
cc            call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .            ,jx_cov(i,j,k),jy_cov(i,j,k),jz_cov(i,j,k)
cc     .            ,jx    (i,j,k),jy    (i,j,k),jz    (i,j,k),.false.)
cc          enddo
cc        enddo
cc      enddo

c Electron velocity

      vex = vx -di*jx/rho
      vey = vy -di*jy/rho
      vez = vz -di*jz/rho

c Transport parameters

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            nuu  (i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
            eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

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
