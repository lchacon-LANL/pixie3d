c evaluateDiagnostics
c####################################################################
      subroutine evaluateDiagnostics

c--------------------------------------------------------------------
c     Calculates diagnostics.
c--------------------------------------------------------------------

      use variables

      use timeStepping

      use graphics

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

      integer(4) :: i,j,k,ieq
      real(8)    :: dpert(neqd),mag(neqd),dmag1,dmag2
      real(8)    :: array(0:nxdp,0:nydp,0:nzdp)
     .             ,vx0  (0:nxdp,0:nydp,0:nzdp)
     .             ,vy0  (0:nxdp,0:nydp,0:nzdp)
     .             ,vz0  (0:nxdp,0:nydp,0:nzdp)
cc      real(8)    :: cs(1,nxd-1),sn(1,nxd-1)

      real(8)    :: Em,Ek,Et,energy,diverB,Npar,px,py,pz

      real(8),pointer,dimension(:,:,:) :: rho0,rvx0,rvy0,rvz0
     .                                   ,bx0,by0,bz0,tmp0

c Externals

      real(8)    :: integral,vectorNorm
      external   :: integral,vectorNorm

c Begin program

      rho0 => u_0%array_var(IRHO)%array
      rvx0 => u_0%array_var(IVX )%array
      rvy0 => u_0%array_var(IVY )%array
      rvz0 => u_0%array_var(IVZ )%array
      bx0  => u_0%array_var(IBX )%array
      by0  => u_0%array_var(IBY )%array
      bz0  => u_0%array_var(IBZ )%array
      tmp0 => u_0%array_var(ITMP)%array

      array = 0d0

      where (rho0 /= 0d0)
        vx0 = rvx0/rho0
        vy0 = rvy0/rho0
        vz0 = rvz0/rho0
      end where

c Particle diagnostics

      Npar = integral(nxd,nyd,nzd,rho ,igx,igy,igz)
     .      -integral(nxd,nyd,nzd,rho0,igx,igy,igz)

c Momentum diagnostics

      px = integral(nxd,nyd,nzd,rvx ,igx,igy,igz)
     .    -integral(nxd,nyd,nzd,rvx0,igx,igy,igz)

      py = integral(nxd,nyd,nzd,rvy ,igx,igy,igz)
     .    -integral(nxd,nyd,nzd,rvy0,igx,igy,igz)

      pz = integral(nxd,nyd,nzd,rvz ,igx,igy,igz)
     .    -integral(nxd,nyd,nzd,rvz0,igx,igy,igz)

c Energy diagnostics

      !Magnetic energy

      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            array(i,j,k) = vectorNorm(i,j,k,igx,igy,igz
     .                               ,bx (i,j,k),by (i,j,k),bz (i,j,k))
     .                    -vectorNorm(i,j,k,igx,igy,igz
     .                               ,bx0(i,j,k),by0(i,j,k),bz0(i,j,k))
          enddo
        enddo
      enddo

      Em = 0.5*integral(nxd,nyd,nzd,array,igx,igy,igz)

      !Ion kinetic energy

      do k = 0,nz+1
        do j = 0,ny+1
          do i = 0,nx+1
            array(i,j,k) = rho (i,j,k)*vectorNorm(i,j,k,igx,igy,igz
     .                               ,vx (i,j,k),vy (i,j,k),vz (i,j,k))
     .                    -rho0(i,j,k)*vectorNorm(i,j,k,igx,igy,igz
     .                               ,vx0(i,j,k),vy0(i,j,k),vz0(i,j,k))
          enddo
        enddo
      enddo

      Ek = 0.5*integral(nxd,nyd,nzd,array,igx,igy,igz)

      !Thermal energy

      if (gamma /= 1d0) then
        array = 2*rho *tmp /(gamma-1.)
     .         -2*rho0*tmp0/(gamma-1.)
      else
        array = 0d0
      endif

      Et = integral(nxd,nyd,nzd,array,igx,igy,igz)

      !Total energy

      energy = Ek + Em + Et

c div(B) diagnostics

      !Do NOT change to 0:nx+1, etc.
      do k = 1,nz-1
        do j = 1,ny-1
          do i = 1,nx-1
            divrgB(i,j,k) = divB(i,j,k)
          enddo
        enddo
      enddo

      divrgB = divrgB**2

      diverB = integral(nxd,nyd,nzd,divrgB,igx,igy,igz)

      diverB = sqrt(diverB)

cc      do k = 0,nz+1
cc        do j = 0,ny+1
cc          do i = 0,nx+1
cc            array(i,j,k) = divB(i,j,k)**2
cc          enddo
cc        enddo
cc      enddo
cc
cc      diverB = integral(nxd,nyd,nzd,array,igx,igy,igz)
cc
cc      diverB = sqrt(diverB)

c Growth rate diagnostics

      do ieq=1,neqd

        array = (u_np%array_var(ieq)%array
     .          -u_0 %array_var(ieq)%array )**2

        dpert(ieq) = integral(nxd,nyd,nzd,array,igx,igy,igz)

        if (dpert(ieq).gt.0d0) dpert(ieq) = log(sqrt(dpert(ieq)))

      enddo

c B x rot(jz) diagnostic
cc
cc      do i=1,nxd
cc        do j=1,nyd
cc          array(i,j) = sqrt( ( (curr(i+1,j)-curr(i-1,j))/2./dxx )**2
cc     .                      +( (curr(i,j+1)-curr(i,j-1))/2./dyy )**2 )
cc        enddo
cc      enddo
cc
cc      array = eeta*array
cc
cc      brotj = integral(nxd,nyd,dxx,dyy,array)

ccc Cosine and sine transform of flux
cc
cc      cs(1,1:nxd-1) = cos(2*pi*xx(1:nxd-1)/xlength)
cc      sn(1,1:nxd-1) = sin(2*pi*xx(1:nxd-1)/xlength)
cc
cc      psic = sum(matmul(matmul(cs,t(1:nxd-1,:)),sin(pi*yy)))
cc     .       *dxx*dyy/xlength
cc      psis = sum(matmul(matmul(sn,t(1:nxd-1,:)),sin(pi*yy)))
cc     .       *dxx*dyy/xlength
cc
cccc      psic = log(abs(psic))
cccc      psis = log(abs(psis))

ccc RMS flux diagnostic
cc
cc      psi_rms = sum( (sum(t(1:nxd-1,:)**2,1)*dxx/xlength 
cc     .             - (sum(t(1:nxd-1,:)   ,1)*dxx/xlength)**2) )*dyy
cc
cccc      write (*,*) psi_rms
cccc      psi_rms = log(sqrt(psi_rms))
cc      psi_rms = sqrt(psi_rms)
cc
ccc Flux penetration diagnostic
cc
cc      fpen = 0.5*(maxval(t(1:nxd-1,nyd+1))-minval(t(1:nxd-1,nyd+1)))

c Calculation of local growth rate for CN

      do ieq=1,neqd

        array = (u_n%array_var(ieq)%array
     .          -u_0%array_var(ieq)%array )**2

        dmag1 = integral(nxd,nyd,nzd,array,igx,igy,igz)

        array = (u_np%array_var(ieq)%array
     .          +u_n %array_var(ieq)%array
     .       -2.*u_0 %array_var(ieq)%array )**2

        dmag2 = integral(nxd,nyd,nzd,array,igx,igy,igz)

        if (dpert(ieq).ne.0d0) then
          mag(ieq) = .5*dt*sqrt(dmag2)/(exp(dpert(ieq))-sqrt(dmag1))
        else
          mag(ieq) = 1e30
        endif

      enddo

      gammat = 1./minval(abs(mag))

c Diagnostic assignment

      diagnostics(1:neqd) = dpert(:)
      diagnostics(neqd+1) = Em
      diagnostics(neqd+2) = Ek
      diagnostics(neqd+3) = Et
      diagnostics(neqd+4) = energy
      diagnostics(neqd+5) = dt
      diagnostics(neqd+6) = gammat
      diagnostics(neqd+7) = diverB
      diagnostics(neqd+8) = Npar
      diagnostics(neqd+9) = px
      diagnostics(neqd+10)= py
      diagnostics(neqd+11)= pz

      diagnostics(neqd+12:20) = 0d0

cc      do i=1,20
cc        write (*,*) diag_desc(i),diagnostics(i)
cc      enddo

c End 

      end subroutine evaluateDiagnostics

c vectorNorm
c ####################################################################
      real(8) function vectorNorm(i,j,k,igx,igy,igz,ax,ay,az)

c -------------------------------------------------------------------
c     Finds norm of vector A given its contravariant components.
c -------------------------------------------------------------------

      use grid

      use icond

      implicit none

c Call variables

      integer(4) :: i,j,k,igx,igy,igz
      real(8)    :: ax,ay,az

c Local variables

      integer(4) :: ig,jg,kg
      real(8)    :: x1,y1,z1,gsub(3,3),cnv(3),cov(3)

c Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      x1 = grid_params%xx(ig)
      y1 = grid_params%yy(jg)
      z1 = grid_params%zz(kg)

      gsub = G_sub(x1,y1,z1)

      cnv = (/ ax,ay,az /)
      cov = matmul(gsub,cnv)

      vectorNorm = dot_product(cov,cnv)/jacobian(x1,y1,z1)

c End 

      end function vectorNorm

c integral
c ####################################################################
      real(8) function integral(nx,ny,nz,array,igx,igy,igz)

c -------------------------------------------------------------------
c     Integrates array(i,j) on domain (nx)x(ny) and returns result.
c -------------------------------------------------------------------

      use grid

      implicit none

c Call variables

      integer(4) :: igx,igy,igz,nx,ny,nz
      real(8)    :: array(0:nx+1,0:ny+1,0:nz+1)

c Local variables

      integer(4) :: i,j,k
      integer(4) :: imin,imax,jmin,jmax,kmin,kmax

c Begin program

c Find integral limits

      imin = 1
      imax = nx

      jmin = 1
      jmax = ny

      kmin = 1
      kmax = nz

cc      if (bcond(1) == PER) then
cc        imin = 1
cc        imax = nx
cc      elseif (bcond(1) == SP) then
cc        imin = 1
cc        imax = nx+1
cc      else
cc        imin = 0
cc        imax = nx+1
cc      endif
cc
cc      if (bcond(3) == PER) then
cc        jmin = 1
cc        jmax = ny
cc      elseif (bcond(3) == SP) then
cc        jmin = 1
cc        jmax = ny+1
cc      else
cccc        jmin = 0
cccc        jmax = ny+1
cc        jmin = 1
cc        jmax = ny
cc      endif
cc
cc      if (bcond(5) == PER) then
cc        kmin = 1
cc        kmax = nz
cc      elseif (bcond(5) == SP) then
cc        kmin = 1
cc        kmax = nz+1
cc      else
cc        kmin = 0
cc        kmax = nz+1
cc      endif

c Integrate

      integral = 0d0

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax
            integral = integral
     .               + array(i,j,k)*volume(i,j,k,igx,igy,igz)
          enddo
        enddo
      enddo

c End 

      end function integral

