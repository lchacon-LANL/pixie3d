c nonlinearRHS
c#################################################################
      subroutine nonlinearRHS(i,j,k,varray,ff)

c-----------------------------------------------------------------
c     This function computes the nonlinear rhs for the 
c     2-D Hall MHD problem. Nonlinear functions Fi should be
c     coded with the following convention:
c        dt Ui + Fi(Uj) = 0
c     where Ui are the dependent variables, given in varray.
c     The result of Fi(Uj) at node (i,j) for all equations neqd
c     is returned in the vector ff of dimension neqd.
c-----------------------------------------------------------------

      use variables

      use timeStepping

      use nlfunction_setup

      use equilibrium

      implicit none

c Call variables

      real(8) :: ff(neqd)

      integer :: i,j,k
      
      type (var_array),target :: varray

c Local variables

      integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km

      real(8)    :: dvol,dS1,dS2,dS3,dxx,dyy,dzz

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,dummy

      ! rho equation
      real(8)    :: advec,diffus

      ! Faraday's law
      real(8)    :: Ez_jp,Ez_jm,Ey_kp,Ey_km,Ex_kp,Ex_km
     .             ,Ez_ip,Ez_im,Ey_ip,Ey_im,Ex_jp,Ex_jm

      ! EOM

      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o

      real(8)    :: hess(3,3),msource

      ! tmp equation
      real(8)    :: heat_flx,heat_src,joule,viscous

c Begin program

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

c     Grid parameters

cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      x0 = xx(ig)
cc      y0 = yy(jg)
cc      z0 = zz(kg)
cc
cc      cartesian = .false.

      call getCoordinates(im,j,k,igx,igy,igz,ig,jg,kg,xim,yim,zim
     .                   ,cartesian)
      call getCoordinates(ip,j,k,igx,igy,igz,ig,jg,kg,xip,yip,zip
     .                   ,cartesian)
      call getCoordinates(i,jm,k,igx,igy,igz,ig,jg,kg,xjm,yjm,zjm
     .                   ,cartesian)
      call getCoordinates(i,jp,k,igx,igy,igz,ig,jg,kg,xjp,yjp,zjp
     .                   ,cartesian)
      call getCoordinates(i,j,km,igx,igy,igz,ig,jg,kg,xkm,ykm,zkm
     .                   ,cartesian)
      call getCoordinates(i,j,kp,igx,igy,igz,ig,jg,kg,xkp,ykp,zkp
     .                   ,cartesian)

      call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x0,y0,z0,cartesian)

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      dS1 = dyy*dzz
      dS2 = dxx*dzz
      dS3 = dxx*dyy

      dvol = dxx*dyy*dzz

      sing_point = .false.
      if (i == 1 .and. bcond(1) == SP) sing_point = .true.

      gsub   = g_sub   (x0,y0,z0,cartesian)
      gsuper = g_super (x0,y0,z0,cartesian)
      jac    = jacobian(x0,y0,z0,cartesian)

      cnv1   = contravariantVector(1,x0,y0,z0,cartesian)
      cnv2   = contravariantVector(2,x0,y0,z0,cartesian)
      cnv3   = contravariantVector(3,x0,y0,z0,cartesian)

c     Rho

      flxip = 0.5*(vx(ip,j,k)*rho(i,j,k) + vx(i,j,k)*rho(ip,j,k))
      flxim = 0.5*(vx(im,j,k)*rho(i,j,k) + vx(i,j,k)*rho(im,j,k))

      flxjp = 0.5*(vy(i,jp,k)*rho(i,j,k) + vy(i,j,k)*rho(i,jp,k))
      flxjm = 0.5*(vy(i,jm,k)*rho(i,j,k) + vy(i,j,k)*rho(i,jm,k))

      flxkp = 0.5*(vz(i,j,kp)*rho(i,j,k) + vz(i,j,k)*rho(i,j,kp))
      flxkm = 0.5*(vz(i,j,km)*rho(i,j,k) + vz(i,j,k)*rho(i,j,km))

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .                      ,varray%array_var(IRHO)%bconds)

      advec = dS1*(flxip - flxim)
     .      + dS2*(flxjp - flxjm)
     .      + dS3*(flxkp - flxkm)

cc      diffus = dd*laplacian(i,j,k,rho)

      ff(IRHO) = advec - diffus

c     Bx

      flxip = 0d0
      flxim = 0d0

      !Resistive contribution
      Ez_jp = 0.5*(eeta(i,jp,k)*jz_cov(i,jp,k)
     .           + eeta(i,j ,k)*jz_cov(i,j ,k))
      Ez_jm = 0.5*(eeta(i,jm,k)*jz_cov(i,jm,k)
     .           + eeta(i,j ,k)*jz_cov(i,j ,k))

      Ey_kp = 0.5*(eeta(i,j,kp)*jy_cov(i,j,kp)
     .           + eeta(i,j,k )*jy_cov(i,j,k ))
      Ey_km = 0.5*(eeta(i,j,km)*jy_cov(i,j,km)
     .           + eeta(i,j,k )*jy_cov(i,j,k ))

      !Ideal contribution
      flxjp = Ez_jp
     .       +0.5*(vy(i,jp,k)*bx(i,jp,k) + vy(i,j,k)*bx(i,j,k))
     .       -0.5*(vx(i,jp,k)*by(i,jp,k) + vx(i,j,k)*by(i,j,k))
      flxjm = Ez_jm
     .       +0.5*(vy(i,jm,k)*bx(i,jm,k) + vy(i,j,k)*bx(i,j,k))
     .       -0.5*(vx(i,jm,k)*by(i,jm,k) + vx(i,j,k)*by(i,j,k))

      flxkp =-Ey_kp
     .       +0.5*(vz(i,j,kp)*bx(i,j,kp) + vz(i,j,k)*bx(i,j,k))
     .       -0.5*(vx(i,j,kp)*bz(i,j,kp) + vx(i,j,k)*bz(i,j,k))
      flxkm =-Ey_km
     .       +0.5*(vz(i,j,km)*bx(i,j,km) + vz(i,j,k)*bx(i,j,k))
     .       -0.5*(vx(i,j,km)*bz(i,j,km) + vx(i,j,k)*bz(i,j,k))

cc      flxjp = Ez_jp
cc     .       +0.5*(vy(i,jp,k)*bx(i,j,k) + vy(i,j,k)*bx(i,jp,k))
cc     .       -0.5*(vx(i,jp,k)*by(i,j,k) + vx(i,j,k)*by(i,jp,k))
cc      flxjm = Ez_jm
cc     .       +0.5*(vy(i,jm,k)*bx(i,j,k) + vy(i,j,k)*bx(i,jm,k))
cc     .       -0.5*(vx(i,jm,k)*by(i,j,k) + vx(i,j,k)*by(i,jm,k))
cc
cc      flxkp =-Ey_kp
cc     .       +0.5*(vz(i,j,kp)*bx(i,j,k) + vz(i,j,k)*bx(i,j,kp))
cc     .       -0.5*(vx(i,j,kp)*bz(i,j,k) + vx(i,j,k)*bz(i,j,kp))
cc      flxkm =-Ey_km
cc     .       +0.5*(vz(i,j,km)*bx(i,j,k) + vz(i,j,k)*bx(i,j,km))
cc     .       -0.5*(vx(i,j,km)*bz(i,j,k) + vx(i,j,k)*bz(i,j,km))

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IBX)%bconds)

      ff(IBX) = jac*( dS1*(flxip-flxim)
     .               +dS2*(flxjp-flxjm)
     .               +dS3*(flxkp-flxkm) )

c     By

      flxjp = 0d0
      flxjm = 0d0

      !Resistive contribution
      Ez_ip = 0.5*(eeta(ip,j,k)*jz_cov(ip,j,k)
     .           + eeta(i ,j,k)*jz_cov(i ,j,k))
      Ez_im = 0.5*(eeta(im,j,k)*jz_cov(im,j,k)
     .           + eeta(i ,j,k)*jz_cov(i ,j,k))

      Ex_kp = 0.5*(eeta(i,j,kp)*jx_cov(i,j,kp)
     .           + eeta(i,j,k )*jx_cov(i,j,k ))
      Ex_km = 0.5*(eeta(i,j,km)*jx_cov(i,j,km)
     .           + eeta(i,j,k )*jx_cov(i,j,k ))

      !Ideal contribution
      flxip =-Ez_ip
     .       +0.5*(vx(ip,j,k)*by(ip,j,k) + vx(i,j,k)*by(i,j,k))
     .       -0.5*(vy(ip,j,k)*bx(ip,j,k) + vy(i,j,k)*bx(i,j,k))
      flxim =-Ez_im
     .       +0.5*(vx(im,j,k)*by(im,j,k) + vx(i,j,k)*by(i,j,k))
     .       -0.5*(vy(im,j,k)*bx(im,j,k) + vy(i,j,k)*bx(i,j,k))

      flxkp = Ex_kp
     .       +0.5*(vz(i,j,kp)*by(i,j,kp) + vz(i,j,k)*by(i,j,k))
     .       -0.5*(vy(i,j,kp)*bz(i,j,kp) + vy(i,j,k)*bz(i,j,k))
      flxkm = Ex_km
     .       +0.5*(vz(i,j,km)*by(i,j,km) + vz(i,j,k)*by(i,j,k))
     .       -0.5*(vy(i,j,km)*bz(i,j,km) + vy(i,j,k)*bz(i,j,k))

cc      flxip =-Ez_ip
cc     .       +0.5*(vx(ip,j,k)*by(i,j,k) + vx(i,j,k)*by(ip,j,k))
cc     .       -0.5*(vy(ip,j,k)*bx(i,j,k) + vy(i,j,k)*bx(ip,j,k))
cc      flxim =-Ez_im
cc     .       +0.5*(vx(im,j,k)*by(i,j,k) + vx(i,j,k)*by(im,j,k))
cc     .       -0.5*(vy(im,j,k)*bx(i,j,k) + vy(i,j,k)*bx(im,j,k))
cc
cc      flxkp = Ex_kp
cc     .       +0.5*(vz(i,j,kp)*by(i,j,k) + vz(i,j,k)*by(i,j,kp))
cc     .       -0.5*(vy(i,j,kp)*bz(i,j,k) + vy(i,j,k)*bz(i,j,kp))
cc      flxkm = Ex_km
cc     .       +0.5*(vz(i,j,km)*by(i,j,k) + vz(i,j,k)*by(i,j,km))
cc     .       -0.5*(vy(i,j,km)*bz(i,j,k) + vy(i,j,k)*bz(i,j,km))

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IBY)%bconds)

      ff(IBY) = jac*( dS1*(flxip-flxim)
     .               +dS2*(flxjp-flxjm)
     .               +dS3*(flxkp-flxkm) )

c     Bz

      flxkp = 0d0
      flxkm = 0d0

      !Resistive contribution
      Ey_ip = 0.5*(eeta(ip,j,k)*jy_cov(ip,j,k)
     .           + eeta(i ,j,k)*jy_cov(i ,j,k))
      Ey_im = 0.5*(eeta(im,j,k)*jy_cov(im,j,k)
     .           + eeta(i ,j,k)*jy_cov(i ,j,k))

      Ex_jp = 0.5*(eeta(i,jp,k)*jx_cov(i,jp,k)
     .           + eeta(i,j ,k)*jx_cov(i,j ,k))
      Ex_jm = 0.5*(eeta(i,jm,k)*jx_cov(i,jm,k)
     .           + eeta(i,j ,k)*jx_cov(i,j ,k))

      !Ideal contribution
      flxip = Ey_ip
     .       +0.5*(vx(ip,j,k)*bz(ip,j,k) + vx(i,j,k)*bz(i,j,k))
     .       -0.5*(vz(ip,j,k)*bx(ip,j,k) + vz(i,j,k)*bx(i,j,k))
      flxim = Ey_im
     .       +0.5*(vx(im,j,k)*bz(im,j,k) + vx(i,j,k)*bz(i,j,k))
     .       -0.5*(vz(im,j,k)*bx(im,j,k) + vz(i,j,k)*bx(i,j,k))

      flxjp =-Ex_jp
     .       +0.5*(vy(i,jp,k)*bz(i,jp,k) + vy(i,j,k)*bz(i,j,k))
     .       -0.5*(vz(i,jp,k)*by(i,jp,k) + vz(i,j,k)*by(i,j,k))
      flxjm =-Ex_jm
     .       +0.5*(vy(i,jm,k)*bz(i,jm,k) + vy(i,j,k)*bz(i,j,k))
     .       -0.5*(vz(i,jm,k)*by(i,jm,k) + vz(i,j,k)*by(i,j,k))

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IBZ)%bconds)

      ff(IBZ) = jac*( dS1*(flxip-flxim)
     .               +dS2*(flxjp-flxjm)
     .               +dS3*(flxkp-flxkm) )

c     Temperature

      flxip = (vx(ip,j,k)*tmp(i,j,k) + vx(i,j,k)*tmp(ip,j,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vx(ip,j,k)+vx(i,j,k))/2.
      flxim = (vx(im,j,k)*tmp(i,j,k) + vx(i,j,k)*tmp(im,j,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vx(im,j,k)+vx(i,j,k))/2.

      flxjp = (vy(i,jp,k)*tmp(i,j,k) + vy(i,j,k)*tmp(i,jp,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vy(i,jp,k)+vy(i,j,k))/2.
      flxjm = (vy(i,jm,k)*tmp(i,j,k) + vy(i,j,k)*tmp(i,jm,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vy(i,jm,k)+vy(i,j,k))/2.

      flxkp = (vz(i,j,kp)*tmp(i,j,k) + vz(i,j,k)*tmp(i,j,kp))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vz(i,j,kp)+vz(i,j,k))/2.
      flxkm = (vz(i,j,km)*tmp(i,j,k) + vz(i,j,k)*tmp(i,j,km))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vz(i,j,km)+vz(i,j,k))/2.

      !Heat flux
cc      heat_flx = -chi*laplacian(i,j,k,tmp)

      !Joule heating
      joule = dxh(ig)*dyh(jg)*dzh(kg)*eeta(i,j,k)
     .                               *( jx(i,j,k)*jx_cov(i,j,k)
     .                                 +jy(i,j,k)*jy_cov(i,j,k)
     .                                 +jz(i,j,k)*jz_cov(i,j,k) )
cc      joule = jouleHeating(i,j,k)

      !Viscous heating
      if (sing_point) then
        nabla_v = fnabla_v(i,j,k,0.5*(xip+x0),0.5*(yip+y0)
     .                    ,0.5*(zip+z0),1)
      else
        nabla_v = fnabla_v(i,j,k,x0,y0,z0,0)
      endif
      cov_tnsr =           matmul(gsub  ,nabla_v)
      cnv_tnsr = transpose(matmul(gsuper,nabla_v))
      viscous = dxh(ig)*dyh(jg)*dzh(kg)
     .         *nuu(i,j,k)/jac*sum(cov_tnsr*cnv_tnsr)

      !Heat source
      heat_src = joule/rho(i,j,k) + viscous

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(ITMP)%bconds)

      ff(ITMP) = dS1*(flxip-flxim)
     .          +dS2*(flxjp-flxjm)
     .          +dS3*(flxkp-flxkm)
cc     .          +(gamma-1.)*heat_flx
cc     .          -(gamma-1.)*heat_src

c     Vx

      call vtensor_x(i  ,j,k,t11p,t12p,t13p, 1)
      call vtensor_x(i-1,j,k,t11m,t12m,t13m,-1)
      if (coords /= 'car') call vtensor_x(i  ,j,k,t11o,t12o,t13o, 0)

      call vtensor_y(i,j  ,k,t21p,t22p,t23p, 1)
      call vtensor_y(i,j-1,k,t21m,t22m,t23m,-1)
      if (coords /= 'car') call vtensor_y(i,j  ,k,t21o,t22o,t23o, 0)

      call vtensor_z(i,j,k  ,t31p,t32p,t33p, 1)
      call vtensor_z(i,j,k-1,t31m,t32m,t33m,-1)
      if (coords /= 'car') call vtensor_z(i,j,k  ,t31o,t32o,t33o, 0)

      cov   = covariantVector(1,x0,y0,z0,cartesian)

      flxip = vflx_x(t11p,t12p,t13p,cov)
      flxim = vflx_x(t11m,t12m,t13m,cov)

      flxjp = vflx_y(t21p,t22p,t23p,cov)
      flxjm = vflx_y(t21m,t22m,t23m,cov)

      flxkp = vflx_z(t31p,t32p,t33p,cov)
      flxkm = vflx_z(t31m,t32m,t33m,cov)

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IVX)%bconds)

      if (coords == 'car') then
        msource = 0d0
      else
        hess = hessian(1,x0,y0,z0,cartesian)
        msource = dvol*( t11o*hess(1,1)+t12o*hess(1,2)+t13o*hess(1,3)
     .                  +t21o*hess(2,1)+t22o*hess(2,2)+t23o*hess(2,3)
     .                  +t31o*hess(3,1)+t32o*hess(3,2)+t33o*hess(3,3))
      endif

      ff(IVX) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm)
     .              - msource )

c     Vy

      cov   = covariantVector(2,x0,y0,z0,cartesian)

      flxip = vflx_x(t11p,t12p,t13p,cov)
      flxim = vflx_x(t11m,t12m,t13m,cov)

      flxjp = vflx_y(t21p,t22p,t23p,cov)
      flxjm = vflx_y(t21m,t22m,t23m,cov)

      flxkp = vflx_z(t31p,t32p,t33p,cov)
      flxkm = vflx_z(t31m,t32m,t33m,cov)

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IVY)%bconds)

      if (coords == 'car') then
        msource = 0d0
      else
        hess = hessian(2,x0,y0,z0,cartesian)
        msource =dvol*( t11o*hess(1,1)+t12o*hess(1,2)+t13o*hess(1,3)
     .                 +t21o*hess(2,1)+t22o*hess(2,2)+t23o*hess(2,3)
     .                 +t31o*hess(3,1)+t32o*hess(3,2)+t33o*hess(3,3))
      endif

      ff(IVY) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm)
     .              - msource )

c     Vz

      cov   = covariantVector(3,x0,y0,z0,cartesian)

      flxip = vflx_x(t11p,t12p,t13p,cov)
      flxim = vflx_x(t11m,t12m,t13m,cov)

      flxjp = vflx_y(t21p,t22p,t23p,cov)
      flxjm = vflx_y(t21m,t22m,t23m,cov)

      flxkp = vflx_z(t31p,t32p,t33p,cov)
      flxkm = vflx_z(t31m,t32m,t33m,cov)

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IVZ)%bconds)

      if (coords == 'car') then
        msource = 0d0
      else
        hess = hessian(3,x0,y0,z0,cartesian)
        msource =dvol*( t11o*hess(1,1)+t12o*hess(1,2)+t13o*hess(1,3)
     .                 +t21o*hess(2,1)+t22o*hess(2,2)+t23o*hess(2,3)
     .                 +t31o*hess(3,1)+t32o*hess(3,2)+t33o*hess(3,3))
      endif

      ff(IVZ) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm)
     .              - msource )

c End

      end subroutine

c defineTSParameters
c####################################################################
      subroutine defineTSParameters(cnf,one_over_dt)
c--------------------------------------------------------------------
c     Calculates nonlinear residuals. 
c--------------------------------------------------------------------

      use parameters

      use equilibrium

      use timeStepping

      implicit none

c Call variables

      real(8) :: cnf(neqd),one_over_dt(neqd)

c Local variables

c Begin program

      cnf = cnfactor
      one_over_dt = 1d0/dt

c End program

      end subroutine

c res
c##################################################################
      real*8 function res(i,j,k,nx,ny,nz,igx,igy,igz)
c------------------------------------------------------------------
c     This function computes the resistivity on the grid level
c     igrid.
c------------------------------------------------------------------

      use grid

      use transport_params

      implicit none

c Call variables

      integer*4    i,j,k,nx,ny,nz,igx,igy,igz

c Local variables

c Begin program

      res = eta

c End program

      end function res

c vis
c##################################################################
      real*8 function vis(i,j,k,nx,ny,nz,igx,igy,igz)
c------------------------------------------------------------------
c     This function computes the viscosity on the grid level
c     igrid.
c------------------------------------------------------------------

      use grid

      use transport_params

      implicit none

c Call variables

      integer*4    i,j,k,nx,ny,nz,igx,igy,igz

c Local variables

c Begin program

      vis = nu

c End program

      end function vis
