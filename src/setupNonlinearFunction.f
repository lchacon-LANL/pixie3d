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

      implicit none

c Call variables

      type (var_array),target :: varray

c Local variables

      integer*4    i,j,k,ig,jg,kg

      integer(4) :: dim,loc,ibc,ieq,bctype

      real(8)    :: dh,vmax

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
        write (*,*) rho
        stop
      endif

      if (minval(tmp) < 0d0) then
        write (*,*) 'Warning: negative temperatures are occurring'
      endif

c Calculate auxiliary quantities

cc      !Covariant magnetic field
cccc      do k = 1,nz
cccc        do j = 1,ny
cccc          do i = 1,nx
cc      do k = 0,nz+1
cc        do j = 0,ny+1
cc          do i = 0,nx+1
cc            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc            gsub = G_sub(xx(ig),yy(jg),zz(kg))
cc            cnv = (/ bx(i,j,k),by(i,j,k),bz(i,j,k) /)
cc            cov = matmul(gsub,cnv)
cc            bx_cov(i,j,k) = cov(1)
cc            by_cov(i,j,k) = cov(2)
cc            bz_cov(i,j,k) = cov(3)
cc          enddo
cc        enddo
cc      enddo
cc
cc      !Contravariant current
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc            jx(i,j,k) = (bz_cov(i,j+1,k)-bz_cov(i,j-1,k))/2./dyh(jg)
cc     .                 -(by_cov(i,j,k+1)-by_cov(i,j,k-1))/2./dzh(kg)
cc            jy(i,j,k) = (bx_cov(i,j,k+1)-bx_cov(i,j,k-1))/2./dzh(kg)
cc     .                 -(bz_cov(i+1,j,k)-bz_cov(i-1,j,k))/2./dxh(ig)
cc            jz(i,j,k) = (by_cov(i+1,j,k)-by_cov(i-1,j,k))/2./dxh(ig)
cc     .                 -(bx_cov(i,j+1,k)-bx_cov(i,j-1,k))/2./dyh(jg)
cc          enddo
cc        enddo
cc      enddo
cc
cc      do dim=1,3
cc        do loc=0,1
cc          ibc = (1+loc)+2*(dim-1)
cc
cc          bctype=PER
cc
cc          ieq = IBX
cc          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
cc            call FillGhostNodes(ieq,dim,loc,bctype,jx,zeros)
cc          endif
cc
cc          ieq = IBY
cc          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
cc            call FillGhostNodes(ieq,dim,loc,bctype,jy,zeros)
cc          endif
cc
cc          ieq = IBZ
cc          if (varray%array_var(ieq)%bconds(ibc) == bctype) then
cc            call FillGhostNodes(ieq,dim,loc,bctype,jz,zeros)
cc          endif
cc
cc        enddo
cc      enddo
cc
cc      !Covariant current
cc      do k = 0,nz+1
cc        do j = 0,ny+1
cc          do i = 0,nx+1
cc            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc            gsub = G_sub(xx(ig),yy(jg),zz(kg))
cc            cnv = (/ jx(i,j,k),jy(i,j,k),jz(i,j,k) /)
cc            cov = matmul(gsub,cnv)
cccc            jx_cov(i,j,k) = cov(1)
cccc            jy_cov(i,j,k) = cov(2)
cccc            jz_cov(i,j,k) = cov(3)
cc            jx_cov(i,j,k) = fj_ip(i,j,k,1)
cc            jy_cov(i,j,k) = fj_ip(i,j,k,2)
cc            jz_cov(i,j,k) = fj_ip(i,j,k,3)
cc          enddo
cc        enddo
cc      enddo

c Transport parameters

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            nuu  (i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
            eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

cc      dh = 1./sqrt((1./minval(dx(1:nx)*2)
cc     .             +1./minval(dy(1:ny)*2)
cc     .             +1./minval(dz(1:nz)*2)))
cc
cc      vmax = sqrt(maxval(vx**2+vy**2+vz**2))
cc
cc      kdiv = min(vmax*dh,maxval(eeta))

cc      kdiv = max(dh**2/dt,kdiv)

cc      write (*,*) dh,kdiv

      kdiv = 0e-3

cc      kdiv = sum(eeta)/size(eeta)
cc      write (*,*) kdiv

      eeta = eeta - kdiv

c End program

      return
      end
