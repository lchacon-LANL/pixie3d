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

      logical    :: alt_eom_b

      ! EOM
      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o

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
        diffus = dd*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,rho)
      else
        diffus = 0d0
      endif

      ff(IRHO) = advec - diffus

c     Faraday's law

      alt_eom_b = .false.

      ff(IBX:IBZ)= div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom_b
     .                       ,btensor_x,btensor_y,btensor_z)
     .           + jac*dvol*curl(i,j,k,nx,ny,nz,igx,igy,igz,ejx,ejy,ejz)
cc     $           - eeta(i,j,k)*veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz
cc     $                                     ,bcnv,alt_eom_b)

c     Bx

cc      flxip = 0d0
cc      flxim = 0d0
cc
cc      flxjp = 0.5*(vy(i,jp,k)*bx(i,jp,k)/jacjp
cc     .           + vy(i,j ,k)*bx(i,j ,k)/jac  )
cc     .       -0.5*(vx(i,jp,k)*by(i,jp,k)/jacjp
cc     .           + vx(i,j ,k)*by(i,j ,k)/jac  )
cc      flxjm = 0.5*(vy(i,jm,k)*bx(i,jm,k)/jacjm
cc     .           + vy(i,j ,k)*bx(i,j ,k)/jac  )
cc     .       -0.5*(vx(i,jm,k)*by(i,jm,k)/jacjm
cc     .           + vx(i,j ,k)*by(i,j ,k)/jac  )
cc      Ez_jp = 0.5*(eeta(i,jp,k)*jz_cov(i,jp,k)
cc     .           + eeta(i,j ,k)*jz_cov(i,j ,k))
cc      Ez_jm = 0.5*(eeta(i,jm,k)*jz_cov(i,jm,k)
cc     .           + eeta(i,j ,k)*jz_cov(i,j ,k))
cc
cc      flxkp = 0.5*(vz(i,j,kp)*bx(i,j,kp)/jackp
cc     .           + vz(i,j,k )*bx(i,j,k )/jac  )
cc     .       -0.5*(vx(i,j,kp)*bz(i,j,kp)/jackp
cc     .           + vx(i,j,k )*bz(i,j,k )/jac  )
cc      flxkm = 0.5*(vz(i,j,km)*bx(i,j,km)/jackm
cc     .           + vz(i,j,k )*bx(i,j,k )/jac  )
cc     .       -0.5*(vx(i,j,km)*bz(i,j,km)/jackm
cc     .           + vx(i,j,k )*bz(i,j,k )/jac  )
cc      Ey_kp = 0.5*(eeta(i,j,kp)*jy_cov(i,j,kp)
cc     .           + eeta(i,j,k )*jy_cov(i,j,k ))
cc      Ey_km = 0.5*(eeta(i,j,km)*jy_cov(i,j,km)
cc     .           + eeta(i,j,k )*jy_cov(i,j,k ))
cc
cc      !Add resistive contribution
cc      flxjp = flxjp + Ez_jp
cc      flxjm = flxjm + Ez_jm
cc      flxkp = flxkp - Ey_kp
cc      flxkm = flxkm - Ey_km
cc
cc      !Group fluxes
cc      ff(IBX) = jac*( dS1*(flxip-flxim)
cc     .               +dS2*(flxjp-flxjm)
cc     .               +dS3*(flxkp-flxkm) )
cc
ccc     By
cc
cc      flxjp = 0d0
cc      flxjm = 0d0
cc
cc      !Ideal contribution
cc      flxip = 0.5*(vx(ip,j,k)*by(ip,j,k)/jacip
cc     .            +vx(i ,j,k)*by(i ,j,k)/jac  )
cc     .       -0.5*(vy(ip,j,k)*bx(ip,j,k)/jacip
cc     .            +vy(i ,j,k)*bx(i ,j,k)/jac  )
cc      Ez_ip = 0.5*(eeta(ip,j,k)*jz_cov(ip,j,k)
cc     .           + eeta(i ,j,k)*jz_cov(i ,j,k))
cc      if (.not.sing_point) then
cc        flxim = 0.5*(vx(im,j,k)*by(im,j,k)/jacim
cc     .             + vx(i ,j,k)*by(i ,j,k)/jac  )
cc     .         -0.5*(vy(im,j,k)*bx(im,j,k)/jacim
cc     .             + vy(i ,j,k)*bx(i ,j,k)/jac  )
cc        Ez_im = 0.5*(eeta(im,j,k)*jz_cov(im,j,k)
cc     .             + eeta(i ,j,k)*jz_cov(i ,j,k))
cc      else
cc        flxim = vx(im,j,k)*by(im,j,k)/jacim
cc     .         -vy(im,j,k)*bx(im,j,k)/jacim
cc        Ez_im = eeta(im,j,k)*jz_cov(im,j,k)
cc      endif
cc
cc
cc      flxkp = 0.5*(vz(i,j,kp)*by(i,j,kp)/jackp
cc     .           + vz(i,j,k )*by(i,j,k )/jac  )
cc     .       -0.5*(vy(i,j,kp)*bz(i,j,kp)/jackp
cc     .           + vy(i,j,k )*bz(i,j,k )/jac  )
cc      flxkm = 0.5*(vz(i,j,km)*by(i,j,km)/jackm
cc     .           + vz(i,j,k )*by(i,j,k )/jac  )
cc     .       -0.5*(vy(i,j,km)*bz(i,j,km)/jackm
cc     .           + vy(i,j,k )*bz(i,j,k )/jac)
cc      Ex_kp = 0.5*(eeta(i,j,kp)*jx_cov(i,j,kp)
cc     .           + eeta(i,j,k )*jx_cov(i,j,k ))
cc      Ex_km = 0.5*(eeta(i,j,km)*jx_cov(i,j,km)
cc     .           + eeta(i,j,k )*jx_cov(i,j,k ))
cc
cc      !Add resistive contribution
cc      flxip = flxip - Ez_ip
cc      flxim = flxim - Ez_im
cc
cc      flxkp = flxkp + Ex_kp
cc      flxkm = flxkm + Ex_km
cc
cc      !Group fluxes
cc      ff(IBY) = jac*( dS1*(flxip-flxim)
cc     .               +dS2*(flxjp-flxjm)
cc     .               +dS3*(flxkp-flxkm) )
cc
ccc     Bz
cc
cc      flxkp = 0d0
cc      flxkm = 0d0
cc
cc      if (sing_point) then
cc        jach = 0.5*(jac+jacip)
cc
cc        flxip = 0.5*(vx(ip,j,k)*bz(ip,j,k)/jacip**2
cc     .             + vx(i ,j,k)*bz(i ,j,k)/jac**2  )*jach
cc     .         -0.5*(vz(ip,j,k)*bx(ip,j,k)/jacip**2
cc     .             + vz(i ,j,k)*bx(i ,j,k)/jac**2  )*jach
cccc        flxip = 0.5*(vx(ip,j,k)*bz(ip,j,k)/jacip
cccc     .             + vx(i ,j,k)*bz(i ,j,k)/jac  )
cccc     .         -0.5*(vz(ip,j,k)*bx(ip,j,k)/jacip
cccc     .             + vz(i ,j,k)*bx(i ,j,k)/jac  )
cc        Ey_ip = 0.5*(eeta(ip,j,k)*jy_cov(ip,j,k)/jacip
cc     .             + eeta(i ,j,k)*jy_cov(i ,j,k)/jac  )*jach
cc
cccc        flxim = 0d0 !SP BC on ideal flux (since flx_x ~ r)
cc        flxim = vx(im,j,k)*bz(im,j,k)/jacim
cc     .         -vz(im,j,k)*bx(im,j,k)/jacim
cc        Ey_im = eeta(im,j,k)*jy_cov(im,j,k) !Current at SP
cc      else
cc        flxip = 0.5*(vx(ip,j,k)*bz(ip,j,k)/jacip
cc     .             + vx(i ,j,k)*bz(i ,j,k)/jac  )
cc     .         -0.5*(vz(ip,j,k)*bx(ip,j,k)/jacip
cc     .             + vz(i ,j,k)*bx(i ,j,k)/jac  )
cc        flxim = 0.5*(vx(im,j,k)*bz(im,j,k)/jacim
cc     .             + vx(i ,j,k)*bz(i ,j,k)/jac  )
cc     .         -0.5*(vz(im,j,k)*bx(im,j,k)/jacim
cc     .             + vz(i ,j,k)*bx(i ,j,k)/jac  )
cc
cc        Ey_ip = 0.5*(eeta(ip,j,k)*jy_cov(ip,j,k)
cc     .             + eeta(i ,j,k)*jy_cov(i ,j,k))
cc        Ey_im = 0.5*(eeta(im,j,k)*jy_cov(im,j,k)
cc     .             + eeta(i ,j,k)*jy_cov(i ,j,k))
cc      endif
cc
cc
cc      flxjp = 0.5*(vy(i,jp,k)*bz(i,jp,k)/jacjp
cc     .           + vy(i,j ,k)*bz(i,j ,k)/jac  )
cc     .       -0.5*(vz(i,jp,k)*by(i,jp,k)/jacjp
cc     .           + vz(i,j ,k)*by(i,j ,k)/jac  )
cc      flxjm = 0.5*(vy(i,jm,k)*bz(i,jm,k)/jacjm
cc     .           + vy(i,j ,k)*bz(i,j ,k)/jac  )
cc     .       -0.5*(vz(i,jm,k)*by(i,jm,k)/jacjm
cc     .           + vz(i,j ,k)*by(i,j ,k)/jac  )
cc
cc      Ex_jp = 0.5*(eeta(i,jp,k)*jx_cov(i,jp,k)
cc     .           + eeta(i,j ,k)*jx_cov(i,j ,k))
cc      Ex_jm = 0.5*(eeta(i,jm,k)*jx_cov(i,jm,k)
cc     .           + eeta(i,j ,k)*jx_cov(i,j ,k))
cc
cc      !Add resistive contribution
cc      flxip = flxip + Ey_ip
cc      flxim = flxim + Ey_im
cc
cc      flxjp = flxjp - Ex_jp
cc      flxjm = flxjm - Ex_jm
cc
cc      !Group fluxes
cc      ff(IBZ) = jac*( dS1*(flxip-flxim)
cc     .               +dS2*(flxjp-flxjm)
cc     .               +dS3*(flxkp-flxkm) )

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
cc      heat_flx = -chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,tmp)

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

c     EOM

      ff(IVX:IVZ) = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                        ,vtensor_x,vtensor_y,vtensor_z)

ccc     Vx
cc
cc      call vtensor_x(i ,j,k,igx,igy,igz,t11p,t12p,t13p, 1)
cc      call vtensor_x(im,j,k,igx,igy,igz,t11m,t12m,t13m,-1)
cc      if (coords /= 'car')
cc     .     call vtensor_x(i,j,k,igx,igy,igz,t11o,t12o,t13o, 0)
cc
cc      call vtensor_y(i,j ,k,igx,igy,igz,t21p,t22p,t23p, 1)
cc      call vtensor_y(i,jm,k,igx,igy,igz,t21m,t22m,t23m,-1)
cc      if (coords /= 'car')
cc     .     call vtensor_y(i,j,k,igx,igy,igz,t21o,t22o,t23o, 0)
cc
cc      call vtensor_z(i,j,k ,igx,igy,igz,t31p,t32p,t33p, 1)
cc      call vtensor_z(i,j,km,igx,igy,igz,t31m,t32m,t33m,-1)
cc      if (coords /= 'car')
cc     .     call vtensor_z(i,j,k,igx,igy,igz,t31o,t32o,t33o, 0)
cc
cc      flxip = t11p
cc      flxim = t11m
cc
cc      flxjp = t21p
cc      flxjm = t21m
cc
cc      flxkp = t31p
cc      flxkm = t31m
cc
cc      if (sing_point) flxim = 0d0
cc
cc      if (coords == 'car') then
cc        msource = 0d0
cc      else
cc        msource =dvol
cc     .          *(t11o*hess(1,1,1)+t12o*hess(1,1,2)+t13o*hess(1,1,3)
cc     .           +t21o*hess(1,2,1)+t22o*hess(1,2,2)+t23o*hess(1,2,3)
cc     .           +t31o*hess(1,3,1)+t32o*hess(1,3,2)+t33o*hess(1,3,3))
cc      endif
cc
cc      ff(IVX) = jac*( dS1*(flxip - flxim)
cc     .              + dS2*(flxjp - flxjm)
cc     .              + dS3*(flxkp - flxkm) ) + msource
cc
ccc     Vy
cc
cc      flxip = t12p
cc      flxim = t12m
cc
cc      flxjp = t22p
cc      flxjm = t22m
cc
cc      flxkp = t32p
cc      flxkm = t32m
cc
cc      if (coords == 'car') then
cc        ff(IVY) =  dS1*(flxip - flxim)
cc     .           + dS2*(flxjp - flxjm)
cc     .           + dS3*(flxkp - flxkm)
cc      else
cc        if (alt_eom) then
cc
cc          if (sing_point) flxim = 0d0
cc
cc          msource=dvol
cc     .           *(t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
cc     .            +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
cc     .            +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3)
cc     .            -t12o*(hess(1,1,1)+hess(2,1,2)+hess(3,1,3))
cc     .            -t22o*(hess(1,2,1)+hess(2,2,2)+hess(3,2,3))
cc     .            -t32o*(hess(1,3,1)+hess(2,3,2)+hess(3,3,3)))
cc
cc          ff(IVY) =   dS1*(flxip - flxim)
cc     .              + dS2*(flxjp - flxjm)
cc     .              + dS3*(flxkp - flxkm) + msource
cc        else
cc          msource=dvol
cc     .            *(t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
cc     .             +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
cc     .             +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3))
cc
cc          ff(IVY) = jac*( dS1*(flxip - flxim)
cc     .                  + dS2*(flxjp - flxjm)
cc     .                  + dS3*(flxkp - flxkm) ) + msource
cc        endif
cc      endif
cc
ccc     Vz
cc
cc      flxip = t13p
cc      flxim = t13m
cc
cc      flxjp = t23p
cc      flxjm = t23m
cc
cc      flxkp = t33p
cc      flxkm = t33m
cc
cc      if (sing_point) flxim = 0d0
cc
cc      if (coords == 'car') then
cc        msource = 0d0
cc      else
cc        msource =dvol
cc     .          *(t11o*hess(3,1,1)+t12o*hess(3,1,2)+t13o*hess(3,1,3)
cc     .           +t21o*hess(3,2,1)+t22o*hess(3,2,2)+t23o*hess(3,2,3)
cc     .           +t31o*hess(3,3,1)+t32o*hess(3,3,2)+t33o*hess(3,3,3))
cc      endif
cc
cc      ff(IVZ) = jac*( dS1*(flxip - flxim)
cc     .              + dS2*(flxjp - flxjm)
cc     .              + dS3*(flxkp - flxkm) ) + msource

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

      end subroutine defineTSParameters
