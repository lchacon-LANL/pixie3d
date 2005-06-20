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
     .             ,vlap(3),vlap2(3)

      integer(4) :: ieq
      real(8)    :: cnv(3),cov(3),kappa

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

cc      advec = c_advec(i,j,k,nx,ny,nz,igx,igy,igz,rho,sp_upwind=5)
      advec = c_advec(i,j,k,nx,ny,nz,igx,igy,igz,rho)

      if (dd /= 0d0) then
        diffus = dd*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,rho)
      else
        diffus = 0d0
      endif

      ff(IRHO) = advec - diffus

c     Vector potential

      cov(1) = vey(i,j,k)*bz(i,j,k)
     .       - vez(i,j,k)*by(i,j,k)

      cov(2) = vez(i,j,k)*bx(i,j,k)
     .       - vex(i,j,k)*bz(i,j,k)

      cov(3) = vex(i,j,k)*by(i,j,k)
     .       - vey(i,j,k)*bx(i,j,k)

      vlap  = veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,acov,vol=.false.)
cc      cnv   = veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,acnv
cc     .                     ,.false.,vol=.false.)
cc
cc      call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,vlap(1),vlap(2),vlap(3)
cc     .                        ,cnv (1),cnv (2),cnv (3)
cc     .                        ,.false.)

cc      vlap2 = veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,jcov,vol=.false.)
      vlap2 = 0d0

      cov = cov/jac + eeta(i,j,k)*vlap - heta*vlap2

      ff(IAX:IAZ) = -cov/ivol

c     Temperature

      !Advective part
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

      advec = dS1*(flxip-flxim)
     .       +dS2*(flxjp-flxjm)
     .       +dS3*(flxkp-flxkm)

cc      advec = c_advec(i,j,k,nx,ny,nz,igx,igy,igz,tmp)

cc      advec = advec
cc     .       +(gamma-2.)/ivol*tmp(i,j,k)
cc     .                  *div(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz)

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

      viscous = dvol*rho(i,j,k)*nuu(i,j,k)/jac*sum(cov_tnsr*cnv_tnsr)

      heat_src = joule + viscous
      heat_src = 0d0

      if (heat_src < 0d0) then
        write (*,*) 'Heat source is negative'
        write (*,*) 'Aborting...'
        stop
      endif

      ff(ITMP) = advec
     .          +0.5*(gamma-1.)*(heat_flx - heat_src)/rho(i,j,k)

c     EOM

      if (nc_eom_v.and.nc_eom_b) then
        ff(IVX:IVZ) = 0d0
      else 
        ff(IVX:IVZ) = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                          ,vtensor_x,vtensor_y,vtensor_z)
      endif

cc      if (.not.(nc_eom_v.or.nc_eom_b)) then
cc        ff(IVX:IVZ) = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                          ,vtensor_x,vtensor_y,vtensor_z)
cc      else 
cc        ff(IVX:IVZ) = 0d0
cc      endif

      if (nc_eom_v) then
        nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)

        do ieq=1,3
          cnv(ieq) =( vx(i,j,k)*nabla_v(1,ieq)
     .               +vy(i,j,k)*nabla_v(2,ieq)
     .               +vz(i,j,k)*nabla_v(3,ieq))/jac
        enddo

        cnv = cnv - nuu(i,j,k)
     .          *veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,vcnv
     .                       ,alt_eom,vol=.false.)

        ff(IVX:IVZ) = ff(IVX:IVZ) + rho(i,j,k)*cnv/ivol
      endif

      if (nc_eom_b) then
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

      end subroutine nonlinearRHS

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
