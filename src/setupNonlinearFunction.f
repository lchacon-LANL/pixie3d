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

      implicit none

c Call variables

      type (var_array),target :: varray

c Local variables

      integer*4    i,j,k,ig,jg,kg

c Externals

      real*8       vis,res
      external     vis,res

c Begin program

c Aliases

      rho => varray%array_var(IRHO)%array
      rvx => varray%array_var(IVX )%array
      rvy => varray%array_var(IVY )%array
      rvz => varray%array_var(IVZ )%array
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
      tmp => varray%array_var(ITMP)%array

      where (rho /= 0d0)
        vx = rvx/rho
        vy = rvy/rho
        vz = rvz/rho
      end where

      xx => grid_params%xx
      yy => grid_params%yy
      zz => grid_params%zz

      dxh => grid_params%dxh
      dyh => grid_params%dyh
      dzh => grid_params%dzh

      dx => grid_params%dx
      dy => grid_params%dy
      dz => grid_params%dz

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Safeguards

      if (minval(rho) < 0d0) then
        write (*,*) 'Warning: negative densities are occurring'
      endif

      if (minval(tmp) < 0d0) then
        write (*,*) 'Warning: negative temperatures are occurring'
      endif

c Calculate auxiliary quantities

      !Covariant magnetic field
      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
            gsub = G_sub(xx(ig),yy(jg),zz(kg))
            cnv = (/ bx(i,j,k),by(i,j,k),bz(i,j,k) /)
            cov = matmul(gsub,cnv)
            bx_cov(i,j,k) = cov(1)
            by_cov(i,j,k) = cov(2)
            bz_cov(i,j,k) = cov(3)
          enddo
        enddo
      enddo

      !Contravariant current
      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            jx(i,j,k) = (bz_cov(i,j+1,k)-bz_cov(i,j-1,k))/2./dyh(jg)
     .                 -(by_cov(i,j,k+1)-by_cov(i,j,k-1))/2./dzh(kg)
            jy(i,j,k) = (bx_cov(i,j,k+1)-bx_cov(i,j,k-1))/2./dzh(kg)
     .                 -(bz_cov(i+1,j,k)-bz_cov(i-1,j,k))/2./dxh(ig)
            jz(i,j,k) = (by_cov(i+1,j,k)-by_cov(i-1,j,k))/2./dxh(ig)
     .                 -(bx_cov(i,j+1,k)-bx_cov(i,j-1,k))/2./dyh(jg)
          enddo
        enddo
      enddo

      !Transport parameters
      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            nuu  (i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
            eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

c End program

      return
      end
