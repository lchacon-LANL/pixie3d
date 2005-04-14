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

      ! tmp equation
      real(8)    :: heat_flx,heat_src,joule,viscous

      integer(4) :: ieq
      real(8)    :: cnv(3),kappa

c Begin program

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      alt_eom = alt__eom()
cc      alt_eom = .false.

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

      ivol = 1d0/volume(i,j,k,igx,igy,igz)

c     Rho

      !X flux
csp      if (bcond(1) == SP) then
      if (bcSP()) then
        if (i+grid_params%ilo(igx)-1 == 1) then      !Upwind around singular point
          jach = 0.5*(jac+jacip)
          flxip = 0.25*jach
     .        *( (    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)
     .            +abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *rho(i ,j,k)
     .          +(    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)          
     .            -abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *rho(ip,j,k) )
          flxim = 0d0
        elseif (i+grid_params%ilo(igx)-1 < 5) then  !Upwind around singular point
          jach = 0.5*(jac+jacip)
          flxip = 0.25*jach
     .        *( (    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)
     .            +abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *rho(i ,j,k)
     .          +(    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)          
     .            -abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *rho(ip,j,k) )

          jach = 0.5*(jac+jacim)
          flxim = 0.25*jach
     .         *( (    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
     .             +abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *rho(im,j,k)
     .           +(    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
     .             -abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *rho(i ,j,k) )
        elseif (i+grid_params%ilo(igx)-1 == 5) then  !Upwind around singular point
          jach = 0.5*(jac+jacip)
          flxip = 0.5*(vx(ip,j,k)*rho(i ,j,k)/jacip
     .               + vx(i ,j,k)*rho(ip,j,k)/jac  )*jach

          jach = 0.5*(jac+jacim)
          flxim = 0.25*jach
     .         *( (    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
     .             +abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *rho(im,j,k)
     .           +(    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
     .             -abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *rho(i ,j,k) )
        elseif (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
cc        if (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
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

      !Y flux
      if (bcSP() .and. (i+grid_params%ilo(igx)-1 < 5) ) then !Upwind around singular point
cc      if (sing_point) then  !Upwind around singular point
        flxjp = 0.25*( (    (vy(i,j,k)+vy(i,jp,k))
     .                  +abs(vy(i,j,k)+vy(i,jp,k)) ) *rho(i,j ,k)
     .                +(    (vy(i,j,k)+vy(i,jp,k))
     .                  -abs(vy(i,j,k)+vy(i,jp,k)) ) *rho(i,jp,k) )
        flxjm = 0.25*( (    (vy(i,j,k)+vy(i,jm,k))
     .                  +abs(vy(i,j,k)+vy(i,jm,k)) ) *rho(i,jm,k)
     .                +(    (vy(i,j,k)+vy(i,jm,k))
     .                  -abs(vy(i,j,k)+vy(i,jm,k)) ) *rho(i,j ,k) )
      else
        flxjp = 0.5*(vy(i,jp,k)*rho(i,j,k) + vy(i,j,k)*rho(i,jp,k))
        flxjm = 0.5*(vy(i,jm,k)*rho(i,j,k) + vy(i,j,k)*rho(i,jm,k))
      endif

      !Z flux
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

      vec1 => vcnv
      vec2 => bcnv

      if (solenoidal) then
        ff(IBX:IBZ)= div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .                         ,btensor_x,btensor_y,btensor_z)
cc     .           + jac*dvol*curl(i,j,k,nx,ny,nz,igx,igy,igz,ejx,ejy,ejz)
     .             - eeta(i,j,k)*veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                                       ,bcnv,.false.)
      else
        ff(IBX:IBZ)= div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .                         ,btensor_x,btensor_y,btensor_z)
     .             - eeta(i,j,k)*veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                                       ,bcnv,.false.)

        !Marder divergence cleaning

cc        flxip = div(i ,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=1)
cc        if (.not.sing_point) then
cc          flxim = div(im,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=1)
cc        else
cccc          flxim = div(im,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
cc          flxim = 0d0
cc        endif
cc
cc        flxjp = div(i,j ,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=2)
cc        flxjm = div(i,jm,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=2)
cc
cc        flxkp = div(i,j,k ,nx,ny,nz,igx,igy,igz,bx,by,bz,he=3) 
cc        flxkm = div(i,j,km,nx,ny,nz,igx,igy,igz,bx,by,bz,he=3) 

        flxip = 0.5*(div(ip,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
     .              +div(i ,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
        if (.not.sing_point) then
          flxim = 0.5*(div(im,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
     .                +div(i ,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
        else
cc          flxim = 0.5*div(im,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
          flxim = 0d0
        endif

        flxjp = 0.5*(div(i,jp,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
     .              +div(i,j ,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
        flxjm = 0.5*(div(i,jm,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
     .              +div(i,j ,k,nx,ny,nz,igx,igy,igz,bx,by,bz))

        flxkp = 0.5*(div(i,j,kp,nx,ny,nz,igx,igy,igz,bx,by,bz) 
     .              +div(i,j,k ,nx,ny,nz,igx,igy,igz,bx,by,bz))
        flxkm = 0.5*(div(i,j,km,nx,ny,nz,igx,igy,igz,bx,by,bz) 
     .              +div(i,j,k ,nx,ny,nz,igx,igy,igz,bx,by,bz))

        !Take homogeneous dirichlet BCs at boundaries
        if (i == nx) flxip = 0d0
        if (i == 1 ) flxim = 0d0
        if (j == ny) flxjp = 0d0
        if (j == 1 ) flxjm = 0d0
        if (k == nz) flxkp = 0d0
        if (k == 1 ) flxkm = 0d0

        kappa = 10*eeta(i,j,k)

        ff(IBX:IBZ) = ff(IBX:IBZ)
     .              - kappa*gsuper(:,1)*jac*dS1*(flxip-flxim)
     .              - kappa*gsuper(:,2)*jac*dS2*(flxjp-flxjm)
     .              - kappa*gsuper(:,3)*jac*dS3*(flxkp-flxkm)
      endif

      nullify(vec1,vec2)

c     Temperature

csp      if (bcond(1) == SP) then
      if (bcSP()) then
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
      if (chi /= 0d0) then
        heat_flx =-chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,tmp)
      else
        heat_flx = 0d0
      endif

      !Heat sources

      !!!Joule heating
      joule = dvol*eeta(i,j,k)*( jx(i,j,k)*jx_cov(i,j,k)
     .                          +jy(i,j,k)*jy_cov(i,j,k)
     .                          +jz(i,j,k)*jz_cov(i,j,k) )

      !!!Viscous heating
      nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)

      cov_tnsr = matmul(nabla_v,gsub   )
      cnv_tnsr = matmul(gsuper ,nabla_v)
cc      cnv_tnsr = transpose(matmul(gsuper,nabla_v))

      viscous = dvol*rho(i,j,k)*nuu(i,j,k)/jac*sum(cov_tnsr*cnv_tnsr)

      heat_src = joule + viscous

cc      heat_src = 0d0

      if (heat_src < 0d0) then
        write (*,*) 'Heat source is negative'
        write (*,*) 'Aborting...'
        stop
      endif

      if (sing_point) flxim = 0d0

      ff(ITMP) = dS1*(flxip-flxim)
     .          +dS2*(flxjp-flxjm)
     .          +dS3*(flxkp-flxkm)
     .          +0.5*(gamma-1.)*(heat_flx - heat_src)/rho(i,j,k)

c     EOM

      ff(IVX:IVZ) = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                        ,vtensor_x,vtensor_y,vtensor_z)

c diag **** Non-conservative advection
cc      nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)
cc
cc      do ieq=1,3
cc        cnv(ieq) =( vx(i,j,k)*nabla_v(1,ieq)
cc     .             +vy(i,j,k)*nabla_v(2,ieq)
cc     .             +vz(i,j,k)*nabla_v(3,ieq))/jac
cc      enddo
cc
cc      ff(IVX:IVZ) = ff(IVX:IVZ) + rho(i,j,k)*cnv/ivol
c diag ****

      if (nc_eom) then
        cnv(1) = jy_cov(i,j,k)*bz_cov(i,j,k)
     .         - jz_cov(i,j,k)*by_cov(i,j,k)

        cnv(2) = jz_cov(i,j,k)*bx_cov(i,j,k)
     .         - jx_cov(i,j,k)*bz_cov(i,j,k)

        cnv(3) = jx_cov(i,j,k)*by_cov(i,j,k)
     .         - jy_cov(i,j,k)*bx_cov(i,j,k)

        ff(IVX:IVZ) = ff(IVX:IVZ) - cnv/ivol
      endif

c     Divide by cell volume factor

      ff = ff*ivol

c End

      end subroutine

c defineTSParameters
c####################################################################
      subroutine defineTSParameters
c--------------------------------------------------------------------
c     Calculates nonlinear residuals. 
c--------------------------------------------------------------------

      use parameters

      use equilibrium

      use timeStepping

      implicit none

c Call variables

c Local variables

c Begin program

      cnf = cnfactor
      one_over_dt = 1d0/dt

c End program

      end subroutine defineTSParameters
