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

      type (var_array) :: varray

c Local variables

      integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km

      real(8)    :: dvol,ivol,dS1,dS2,dS3,dxx,dyy,dzz,xhm,yhm,zhm

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

      real(8)    :: hess(3,3,3),msource

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

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      dS1 = dyy*dzz
      dS2 = dxx*dzz
      dS3 = dxx*dyy

      dvol = dxx*dyy*dzz

      sing_point = isSP(i,j,k,igx,igy,igz)

      gsub   = gmetric%grid(igx)%gsub(i,j,k,:,:)
      gsuper = gmetric%grid(igx)%gsup(i,j,k,:,:)
      jac    = gmetric%grid(igx)%jac (i,j,k)

      jacip  = gmetric%grid(igx)%jac(ip,j,k)
      jacim  = gmetric%grid(igx)%jac(im,j,k)
      jacjp  = gmetric%grid(igx)%jac(i,jp,k)
      jacjm  = gmetric%grid(igx)%jac(i,jm,k)
      jackp  = gmetric%grid(igx)%jac(i,j,kp)
      jackm  = gmetric%grid(igx)%jac(i,j,km)

      ivol = 1d0/jac/dvol

      if (coords /= 'car') then
        hess = gmetric%grid(igx)%Gamma(i,j,k,:,:,:)
      endif

c     Rho

      if (bcond(1) == SP) then
cc        if (i == 1) then      !Upwind around singular point
cc          jach = 0.5*(jac+jacip)
cc          flxip = 0.25*jach
cc     .        *( (    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)
cc     .            +abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *rho(i ,j,k)
cc     .          +(    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)          
cc     .            -abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *rho(ip,j,k) )
cc          flxim = 0d0
cc        elseif (i < 5) then  !Upwind around singular point
cc          jach = 0.5*(jac+jacip)
cc          flxip = 0.25*jach
cc     .        *( (    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)
cc     .            +abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *rho(i ,j,k)
cc     .          +(    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)          
cc     .            -abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *rho(ip,j,k) )
cc
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.25*jach
cc     .         *( (    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
cc     .             +abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *rho(im,j,k)
cc     .           +(    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
cc     .             -abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *rho(i ,j,k) )
cc        elseif (i == 5) then  !Upwind around singular point
cc          jach = 0.5*(jac+jacip)
cc          flxip = 0.5*(vx(ip,j,k)*rho(i ,j,k)/jacip
cc     .               + vx(i ,j,k)*rho(ip,j,k)/jac  )*jach
cc
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.25*jach
cc     .         *( (    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
cc     .             +abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *rho(im,j,k)
cc     .           +(    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
cc     .             -abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *rho(i ,j,k) )
cc        elseif (i < nx) then
        if (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
          jach = 0.5*(jac+jacip)
          flxip = 0.5*(vx(ip,j,k)*rho(i ,j,k)/jacip
     .               + vx(i ,j,k)*rho(ip,j,k)/jac  )*jach
          jach = 0.5*(jac+jacim)
          flxim = 0.5*(vx(im,j,k)*rho(i ,j,k)/jacim
     .               + vx(i ,j,k)*rho(im,j,k)/jac  )*jach
        elseif (i+grid_params%ilo(igx)-1 == grid_params%nxgl(igx)) then
          flxip = 0.5*(vx(ip,j,k)*rho(i,j,k) + vx(i,j,k)*rho(ip,j,k))

          jach = 0.5*(jac+jacim)
          flxim = 0.5*(vx(im,j,k)*rho(i ,j,k)/jacim
     .               + vx(i ,j,k)*rho(im,j,k)/jac  )*jach
        else
          flxip = 0.5*(vx(ip,j,k)*rho(i,j,k) + vx(i,j,k)*rho(ip,j,k))
          flxim = 0.5*(vx(im,j,k)*rho(i,j,k) + vx(i,j,k)*rho(im,j,k))
        endif
      else
        flxip = 0.5*(vx(ip,j,k)*rho(i,j,k) + vx(i,j,k)*rho(ip,j,k))
        flxim = 0.5*(vx(im,j,k)*rho(i,j,k) + vx(i,j,k)*rho(im,j,k))
      endif

cc      if (sing_point) then  !Upwind around singular point
cc        flxjp = 0.25*( (    (vy(i,j,k)+vy(i,jp,k))
cc     .                  +abs(vy(i,j,k)+vy(i,jp,k)) ) *rho(i,j ,k)
cc     .                +(    (vy(i,j,k)+vy(i,jp,k))
cc     .                  -abs(vy(i,j,k)+vy(i,jp,k)) ) *rho(i,jp,k) )
cc        flxjm = 0.25*( (    (vy(i,j,k)+vy(i,jm,k))
cc     .                  +abs(vy(i,j,k)+vy(i,jm,k)) ) *rho(i,jm,k)
cc     .                +(    (vy(i,j,k)+vy(i,jm,k))
cc     .                  -abs(vy(i,j,k)+vy(i,jm,k)) ) *rho(i,j ,k) )
cc      else
        flxjp = 0.5*(vy(i,jp,k)*rho(i,j,k) + vy(i,j,k)*rho(i,jp,k))
        flxjm = 0.5*(vy(i,jm,k)*rho(i,j,k) + vy(i,j,k)*rho(i,jm,k))
cc      endif

      flxkp = 0.5*(vz(i,j,kp)*rho(i,j,k) + vz(i,j,k)*rho(i,j,kp))
      flxkm = 0.5*(vz(i,j,km)*rho(i,j,k) + vz(i,j,k)*rho(i,j,km))

      if (sing_point) flxim = 0d0

      advec = dS1*(flxip - flxim)
     .      + dS2*(flxjp - flxjm)
     .      + dS3*(flxkp - flxkm)

      if (dd /= 0d0) then
        diffus = dd*laplacian(i,j,k,nx,ny,nz,rho)
      else
        diffus = 0d0
      endif

      ff(IRHO) = advec - diffus

c     Bx

      flxip = 0d0
      flxim = 0d0

      flxjp = 0.5*(vy(i,jp,k)*bx(i,jp,k)/jacjp
     .           + vy(i,j ,k)*bx(i,j ,k)/jac  )
     .       -0.5*(vx(i,jp,k)*by(i,jp,k)/jacjp
     .           + vx(i,j ,k)*by(i,j ,k)/jac  )
      flxjm = 0.5*(vy(i,jm,k)*bx(i,jm,k)/jacjm
     .           + vy(i,j ,k)*bx(i,j ,k)/jac  )
     .       -0.5*(vx(i,jm,k)*by(i,jm,k)/jacjm
     .           + vx(i,j ,k)*by(i,j ,k)/jac  )
      Ez_jp = 0.5*(eeta(i,jp,k)*jz_cov(i,jp,k)
     .           + eeta(i,j ,k)*jz_cov(i,j ,k))
      Ez_jm = 0.5*(eeta(i,jm,k)*jz_cov(i,jm,k)
     .           + eeta(i,j ,k)*jz_cov(i,j ,k))

      flxkp = 0.5*(vz(i,j,kp)*bx(i,j,kp)/jackp
     .           + vz(i,j,k )*bx(i,j,k )/jac  )
     .       -0.5*(vx(i,j,kp)*bz(i,j,kp)/jackp
     .           + vx(i,j,k )*bz(i,j,k )/jac  )
      flxkm = 0.5*(vz(i,j,km)*bx(i,j,km)/jackm
     .           + vz(i,j,k )*bx(i,j,k )/jac  )
     .       -0.5*(vx(i,j,km)*bz(i,j,km)/jackm
     .           + vx(i,j,k )*bz(i,j,k )/jac  )
      Ey_kp = 0.5*(eeta(i,j,kp)*jy_cov(i,j,kp)
     .           + eeta(i,j,k )*jy_cov(i,j,k ))
      Ey_km = 0.5*(eeta(i,j,km)*jy_cov(i,j,km)
     .           + eeta(i,j,k )*jy_cov(i,j,k ))

      !Add resistive contribution
      flxjp = flxjp + Ez_jp
      flxjm = flxjm + Ez_jm
      flxkp = flxkp - Ey_kp
      flxkm = flxkm - Ey_km

      !Group fluxes
      ff(IBX) = jac*( dS1*(flxip-flxim)
     .               +dS2*(flxjp-flxjm)
     .               +dS3*(flxkp-flxkm) )

c     By

      flxjp = 0d0
      flxjm = 0d0

      !Ideal contribution
      flxip = 0.5*(vx(ip,j,k)*by(ip,j,k)/jacip
     .            +vx(i ,j,k)*by(i ,j,k)/jac  )
     .       -0.5*(vy(ip,j,k)*bx(ip,j,k)/jacip
     .            +vy(i ,j,k)*bx(i ,j,k)/jac  )
      Ez_ip = 0.5*(eeta(ip,j,k)*jz_cov(ip,j,k)
     .           + eeta(i ,j,k)*jz_cov(i ,j,k))
      if (.not.sing_point) then
        flxim = 0.5*(vx(im,j,k)*by(im,j,k)/jacim
     .             + vx(i ,j,k)*by(i ,j,k)/jac  )
     .         -0.5*(vy(im,j,k)*bx(im,j,k)/jacim
     .             + vy(i ,j,k)*bx(i ,j,k)/jac  )
        Ez_im = 0.5*(eeta(im,j,k)*jz_cov(im,j,k)
     .             + eeta(i ,j,k)*jz_cov(i ,j,k))
      else
        flxim = vx(im,j,k)*by(im,j,k)/jacim
     .         -vy(im,j,k)*bx(im,j,k)/jacim
        Ez_im = eeta(im,j,k)*jz_cov(im,j,k)
      endif


      flxkp = 0.5*(vz(i,j,kp)*by(i,j,kp)/jackp
     .           + vz(i,j,k )*by(i,j,k )/jac  )
     .       -0.5*(vy(i,j,kp)*bz(i,j,kp)/jackp
     .           + vy(i,j,k )*bz(i,j,k )/jac  )
      flxkm = 0.5*(vz(i,j,km)*by(i,j,km)/jackm
     .           + vz(i,j,k )*by(i,j,k )/jac  )
     .       -0.5*(vy(i,j,km)*bz(i,j,km)/jackm
     .           + vy(i,j,k )*bz(i,j,k )/jac)
      Ex_kp = 0.5*(eeta(i,j,kp)*jx_cov(i,j,kp)
     .           + eeta(i,j,k )*jx_cov(i,j,k ))
      Ex_km = 0.5*(eeta(i,j,km)*jx_cov(i,j,km)
     .           + eeta(i,j,k )*jx_cov(i,j,k ))

      !Add resistive contribution
      flxip = flxip - Ez_ip
      flxim = flxim - Ez_im

      flxkp = flxkp + Ex_kp
      flxkm = flxkm + Ex_km

      !Group fluxes
      ff(IBY) = jac*( dS1*(flxip-flxim)
     .               +dS2*(flxjp-flxjm)
     .               +dS3*(flxkp-flxkm) )

c     Bz

      flxkp = 0d0
      flxkm = 0d0

      if (sing_point) then
        jach = 0.5*(jac+jacip)

        flxip = 0.5*(vx(ip,j,k)*bz(ip,j,k)/jacip**2
     .             + vx(i ,j,k)*bz(i ,j,k)/jac**2  )*jach
     .         -0.5*(vz(ip,j,k)*bx(ip,j,k)/jacip**2
     .             + vz(i ,j,k)*bx(i ,j,k)/jac**2  )*jach
cc        flxip = 0.5*(vx(ip,j,k)*bz(ip,j,k)/jacip
cc     .             + vx(i ,j,k)*bz(i ,j,k)/jac  )
cc     .         -0.5*(vz(ip,j,k)*bx(ip,j,k)/jacip
cc     .             + vz(i ,j,k)*bx(i ,j,k)/jac  )
        Ey_ip = 0.5*(eeta(ip,j,k)*jy_cov(ip,j,k)/jacip
     .             + eeta(i ,j,k)*jy_cov(i ,j,k)/jac  )*jach

cc        flxim = 0d0 !SP BC on ideal flux (since fxl_x ~ r)
        flxim = vx(im,j,k)*bz(im,j,k)/jacim
     .         -vz(im,j,k)*bx(im,j,k)/jacim
        Ey_im = eeta(im,j,k)*jy_cov(im,j,k) !Current at SP
      else
        flxip = 0.5*(vx(ip,j,k)*bz(ip,j,k)/jacip
     .             + vx(i ,j,k)*bz(i ,j,k)/jac  )
     .         -0.5*(vz(ip,j,k)*bx(ip,j,k)/jacip
     .             + vz(i ,j,k)*bx(i ,j,k)/jac  )
        flxim = 0.5*(vx(im,j,k)*bz(im,j,k)/jacim
     .             + vx(i ,j,k)*bz(i ,j,k)/jac  )
     .         -0.5*(vz(im,j,k)*bx(im,j,k)/jacim
     .             + vz(i ,j,k)*bx(i ,j,k)/jac  )

        Ey_ip = 0.5*(eeta(ip,j,k)*jy_cov(ip,j,k)
     .             + eeta(i ,j,k)*jy_cov(i ,j,k))
        Ey_im = 0.5*(eeta(im,j,k)*jy_cov(im,j,k)
     .             + eeta(i ,j,k)*jy_cov(i ,j,k))
      endif


      flxjp = 0.5*(vy(i,jp,k)*bz(i,jp,k)/jacjp
     .           + vy(i,j ,k)*bz(i,j ,k)/jac  )
     .       -0.5*(vz(i,jp,k)*by(i,jp,k)/jacjp
     .           + vz(i,j ,k)*by(i,j ,k)/jac  )
      flxjm = 0.5*(vy(i,jm,k)*bz(i,jm,k)/jacjm
     .           + vy(i,j ,k)*bz(i,j ,k)/jac  )
     .       -0.5*(vz(i,jm,k)*by(i,jm,k)/jacjm
     .           + vz(i,j ,k)*by(i,j ,k)/jac  )

      Ex_jp = 0.5*(eeta(i,jp,k)*jx_cov(i,jp,k)
     .           + eeta(i,j ,k)*jx_cov(i,j ,k))
      Ex_jm = 0.5*(eeta(i,jm,k)*jx_cov(i,jm,k)
     .           + eeta(i,j ,k)*jx_cov(i,j ,k))

      !Add resistive contribution
      flxip = flxip + Ey_ip
      flxim = flxim + Ey_im

      flxjp = flxjp - Ex_jp
      flxjm = flxjm - Ex_jm

      !Group fluxes
      ff(IBZ) = jac*( dS1*(flxip-flxim)
     .               +dS2*(flxjp-flxjm)
     .               +dS3*(flxkp-flxkm) )

c     Temperature

      if (bcond(1) == SP) then
        if (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
          jach = 0.5*(jac+jacip)
          flxip = (vx(ip,j,k)*tmp(i ,j,k)/jacip
     .           + vx(i ,j,k)*tmp(ip,j,k)/jac  )*jach/2.
     .           +(gamma-2.)*tmp(i,j,k)*(vx(ip,j,k)/jacip
     .                                  +vx(i ,j,k)/jac  )*jach/2.
          jach = 0.5*(jac+jacim)
          flxim = (vx(im,j,k)*tmp(i ,j,k)/jacim
     .           + vx(i ,j,k)*tmp(im,j,k)/jac  )*jach/2.
     .           +(gamma-2.)*tmp(i,j,k)*(vx(im,j,k)/jacim
     .                                  +vx(i ,j,k)/jac  )*jach/2.
        elseif (i+grid_params%ilo(igx)-1 == grid_params%nxgl(igx)) then
          flxip = (vx(ip,j,k)*tmp(i,j,k) + vx(i,j,k)*tmp(ip,j,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vx(ip,j,k)+vx(i,j,k))/2.

          jach = 0.5*(jac+jacim)
          flxim = (vx(im,j,k)*tmp(i ,j,k)/jacim
     .           + vx(i ,j,k)*tmp(im,j,k)/jac  )*jach/2.
     .           +(gamma-2.)*tmp(i,j,k)*(vx(im,j,k)/jacim
     .                                  +vx(i ,j,k)/jac  )*jach/2.
        endif
      else
        flxip = (vx(ip,j,k)*tmp(i,j,k) + vx(i,j,k)*tmp(ip,j,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vx(ip,j,k)+vx(i,j,k))/2.
        flxim = (vx(im,j,k)*tmp(i,j,k) + vx(i,j,k)*tmp(im,j,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vx(im,j,k)+vx(i,j,k))/2.
      endif

      flxjp = (vy(i,jp,k)*tmp(i,j,k) + vy(i,j,k)*tmp(i,jp,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vy(i,jp,k)+vy(i,j,k))/2.
      flxjm = (vy(i,jm,k)*tmp(i,j,k) + vy(i,j,k)*tmp(i,jm,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vy(i,jm,k)+vy(i,j,k))/2.

      flxkp = (vz(i,j,kp)*tmp(i,j,k) + vz(i,j,k)*tmp(i,j,kp))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vz(i,j,kp)+vz(i,j,k))/2.
      flxkm = (vz(i,j,km)*tmp(i,j,k) + vz(i,j,k)*tmp(i,j,km))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vz(i,j,km)+vz(i,j,k))/2.

      !Heat flux
cc      heat_flx = -chi*laplacian(i,j,k,nx,ny,nz,tmp)

      !Joule heating
cc      joule = dxh(ig)*dyh(jg)*dzh(kg)*eeta(i,j,k)
cc     .                               *( jx(i,j,k)*jx_cov(i,j,k)
cc     .                                 +jy(i,j,k)*jy_cov(i,j,k)
cc     .                                 +jz(i,j,k)*jz_cov(i,j,k) )
cccc      joule = jouleHeating(i,j,k)
cc
cc      !Viscous heating
cc      if (sing_point) then
cc        nabla_v = fnabla_v(i,j,k,0.5*(xip+x0)
cc     .                          ,0.5*(yip+y0)
cc     .                          ,0.5*(zip+z0),1)
cc      else
cc        nabla_v = fnabla_v(i,j,k,x0,y0,z0,0)
cc      endif
cc      cov_tnsr =           matmul(gsub  ,nabla_v)
cc      cnv_tnsr = transpose(matmul(gsuper,nabla_v))
cc      viscous = dxh(ig)*dyh(jg)*dzh(kg)
cc     .         *nuu(i,j,k)/jac*sum(cov_tnsr*cnv_tnsr)
cc
cc      !Heat source
cc      heat_src = joule/rho(i,j,k) + viscous

      if (sing_point) flxim = 0d0

      ff(ITMP) = dS1*(flxip-flxim)
     .          +dS2*(flxjp-flxjm)
     .          +dS3*(flxkp-flxkm)
cc     .          +(gamma-1.)*heat_flx
cc     .          -(gamma-1.)*heat_src

c     Vx

      call vtensor_x(i ,j,k,t11p,t12p,t13p, 1)
      call vtensor_x(im,j,k,t11m,t12m,t13m,-1)
      if (coords /= 'car') call vtensor_x(i,j,k,t11o,t12o,t13o, 0)

      call vtensor_y(i,j ,k,t21p,t22p,t23p, 1)
      call vtensor_y(i,jm,k,t21m,t22m,t23m,-1)
      if (coords /= 'car') call vtensor_y(i,j,k,t21o,t22o,t23o, 0)

      call vtensor_z(i,j,k ,t31p,t32p,t33p, 1)
      call vtensor_z(i,j,km,t31m,t32m,t33m,-1)
      if (coords /= 'car') call vtensor_z(i,j,k,t31o,t32o,t33o, 0)

      flxip = t11p
      flxim = t11m

      flxjp = t21p
      flxjm = t21m

      flxkp = t31p
      flxkm = t31m

      if (sing_point) flxim = 0d0

      if (coords == 'car') then
        msource = 0d0
      else
        msource =dvol
     .          *(t11o*hess(1,1,1)+t12o*hess(1,1,2)+t13o*hess(1,1,3)
     .           +t21o*hess(1,2,1)+t22o*hess(1,2,2)+t23o*hess(1,2,3)
     .           +t31o*hess(1,3,1)+t32o*hess(1,3,2)+t33o*hess(1,3,3))
      endif

      ff(IVX) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm) ) + msource

c     Vy

      flxip = t12p
      flxim = t12m

      flxjp = t22p
      flxjm = t22m

      flxkp = t32p
      flxkm = t32m

      if (coords == 'car') then
        ff(IVY) =  dS1*(flxip - flxim)
     .           + dS2*(flxjp - flxjm)
     .           + dS3*(flxkp - flxkm)
      else
        if (alt_eom) then

          if (sing_point) flxim = 0d0

          msource=dvol
     .           *(t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
     .            +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
     .            +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3)
     .            -t12o*(hess(1,1,1)+hess(2,1,2)+hess(3,1,3))
     .            -t22o*(hess(1,2,1)+hess(2,2,2)+hess(3,2,3))
     .            -t32o*(hess(1,3,1)+hess(2,3,2)+hess(3,3,3)))

          ff(IVY) =   dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm) + msource
        else
          msource=dvol
     .            *(t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
     .             +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
     .             +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3))

          ff(IVY) = jac*( dS1*(flxip - flxim)
     .                  + dS2*(flxjp - flxjm)
     .                  + dS3*(flxkp - flxkm) ) + msource
        endif
      endif

c     Vz

      flxip = t13p
      flxim = t13m

      flxjp = t23p
      flxjm = t23m

      flxkp = t33p
      flxkm = t33m

      if (sing_point) flxim = 0d0

      if (coords == 'car') then
        msource = 0d0
      else
        msource =dvol
     .          *(t11o*hess(3,1,1)+t12o*hess(3,1,2)+t13o*hess(3,1,3)
     .           +t21o*hess(3,2,1)+t22o*hess(3,2,2)+t23o*hess(3,2,3)
     .           +t31o*hess(3,3,1)+t32o*hess(3,3,2)+t33o*hess(3,3,3))
      endif

      ff(IVZ) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm) ) + msource

c     Divide by cell volume factor

      ff = ff*ivol

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

      use equilibrium

      implicit none

c Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz

c Local variables

      integer(4) :: nn,ig,jg,kg
      real(8)    :: x1,y1,z1,aa
      logical    :: cartsn

c Begin program

c Resistivity profile eta*(1 + aa*x^nn)
c     Coefficient aa is set so that res = 20*eta at wall
c     and nn so that res=??*eta at sing. surf. xs ~ 0.33

      select case (equil)
      case ('rfp1')

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
        nn = 4
        aa = 19.
        res = eta*(1. + aa*grid_params%xx(ig)**nn)

      case default

        res = eta

      end select

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
