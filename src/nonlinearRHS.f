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

      integer(4) :: ig,jg,kg

      real(8)    :: dS1,dS2,dS3,jac

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      ! rho equation
      real(8)    :: advec,diffus

      ! Faraday's law
      real(8)    :: Ez_ipjp,Ez_imjp,Ez_ipjm,Ez_imjm
     .             ,Ey_ipkp,Ey_imkp,Ey_ipkm,Ey_imkm
     .             ,Ex_jpkp,Ex_jmkp,Ex_jpkm,Ex_jmkm
      real(8)    :: Ez_jp,Ez_jm,Ey_kp,Ey_km,Ex_kp,Ex_km
     .             ,Ez_ip,Ez_im,Ey_ip,Ey_im,Ex_jp,Ex_jm

      ! EOM

      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m
     .             ,t21p,t22p,t23p,t21m,t22m,t23m
     .             ,t31p,t32p,t33p,t31m,t32m,t33m

      ! tmp equation
      real(8)    :: heat_flx,heat_src,joule,viscous

c Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dS1 = dyh(jg)*dzh(kg)
      dS2 = dxh(ig)*dzh(kg)
      dS3 = dxh(ig)*dyh(jg)

      jac = jacobian(xx(ig),yy(jg),zz(kg))

      sing_point = .false.

      if (jac == 0d0) then
        jac = jacobian(dxh(ig)/2.,yy(jg),zz(kg)) !To get volume right
                                                 !at singular point
        sing_point = .true.
      endif

c     Rho

cc      flxip = 0.5*(vx(i+1,j,k)*rho(i,j,k) + vx(i,j,k)*rho(i+1,j,k))
cc      flxim = 0.5*(vx(i-1,j,k)*rho(i,j,k) + vx(i,j,k)*rho(i-1,j,k))
cc
cc      flxjp = 0.5*(vy(i,j+1,k)*rho(i,j,k) + vy(i,j,k)*rho(i,j+1,k))
cc      flxjm = 0.5*(vy(i,j-1,k)*rho(i,j,k) + vy(i,j,k)*rho(i,j-1,k))
cc
cc      flxkp = 0.5*(vz(i,j,k+1)*rho(i,j,k) + vz(i,j,k)*rho(i,j,k+1))
cc      flxkm = 0.5*(vz(i,j,k-1)*rho(i,j,k) + vz(i,j,k)*rho(i,j,k-1))

      flxip = 0.5*(rvx(i+1,j,k) + rvx(i,j,k))
      flxim = 0.5*(rvx(i-1,j,k) + rvx(i,j,k))

      flxjp = 0.5*(rvy(i,j+1,k) + rvy(i,j,k))
      flxjm = 0.5*(rvy(i,j-1,k) + rvy(i,j,k))

      flxkp = 0.5*(rvz(i,j,k+1) + rvz(i,j,k))
      flxkm = 0.5*(rvz(i,j,k-1) + rvz(i,j,k))

cc      if (.not.sing_point) then
cc        advec = dS1*(flxip - flxim)
cc     .        + dS2*(flxjp - flxjm)
cc     .        + dS3*(flxkp - flxkm)
cc      else
cc        advec = dS1*(flxip)
cc     .        + dS3*(flxkp - flxkm)
cc      endif

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IRHO)%bconds)

      advec = dS1*(flxip - flxim)
     .      + dS2*(flxjp - flxjm)
     .      + dS3*(flxkp - flxkm)

      diffus = dd*laplacian(i,j,k,rho)

      ff(IRHO) = advec - diffus

c     Bx

      if (.not.sing_point) then

        Ez_ipjp = fEz_ipjp(i  ,j  ,k)
        Ez_imjp = fEz_ipjp(i-1,j  ,k)
        Ez_ipjm = fEz_ipjp(i  ,j-1,k)
        Ez_imjm = fEz_ipjp(i-1,j-1,k)

        Ez_jp = 0.5*(Ez_ipjp + Ez_imjp)
        Ez_jm = 0.5*(Ez_ipjm + Ez_imjm)

        Ey_ipkp = fEy_ipkp(i  ,j,k  )
        Ey_imkp = fEy_ipkp(i-1,j,k  )
        Ey_ipkm = fEy_ipkp(i  ,j,k-1)
        Ey_imkm = fEy_ipkp(i-1,j,k-1)

        Ey_kp = 0.5*(Ey_ipkp + Ey_imkp)
        Ey_km = 0.5*(Ey_ipkm + Ey_imkm)

cc        ff(IBX) = jac*( dS2*( Ez_jp - Ez_jm ) - dS3*( Ey_kp - Ey_km ) )

        Ez_ip = 0d0

      else

        Ez_ipjp = fEz_ipjp(i,1,k)
        Ez_ipjm = fEz_ipjp(i,0,k)

        Ez_ip = 0.5*(Ez_ipjp + Ez_ipjm)*sin(yy(jg))

        Ey_ipkp = fEy_ipkp(i,1,k  )
        Ey_ipkm = fEy_ipkp(i,1,k-1)

        Ey_kp = Ey_ipkp
        Ey_km = Ey_ipkm

cc        ff(IBX) = jac*( dS1*( Ez_ip ) )
cc     .                - dS3*( Ey_kp - Ey_km ) )
      endif

      call imposeBConfluxes (i,j,k,Ez_ip,flxim,Ez_jp,Ez_jm
     .                      ,Ey_kp,Ey_km,varray%array_var(IBX)%bconds)

      ff(IBX) = jac*(  dS1*( Ez_ip )
     .               + dS2*( Ez_jp - Ez_jm ) - dS3*( Ey_kp - Ey_km ) )

c     By

      if (.not.sing_point) then
        Ex_jpkp = fEx_jpkp(i,j  ,k  )
        Ex_jmkp = fEx_jpkp(i,j-1,k  )
        Ex_jpkm = fEx_jpkp(i,j  ,k-1)
        Ex_jmkm = fEx_jpkp(i,j-1,k-1)

        Ex_kp = 0.5*(Ex_jpkp + Ex_jmkp)
        Ex_km = 0.5*(Ex_jpkm + Ex_jmkm)

        Ez_ip = 0.5*(Ez_ipjp + Ez_ipjm)
        Ez_im = 0.5*(Ez_imjp + Ez_imjm)

cc        ff(IBY) = jac*( dS3*( Ex_kp - Ex_km ) - dS1*( Ez_ip - Ez_im ) )
      else
cc        Ex_jpkp = fEx_jpkp(i,1,k  )
cc        Ex_jpkm = fEx_jpkp(i,1,k-1)
cc
cc        Ex_kp = Ex_jpkp
cc        Ex_km = Ex_jpkm

        Ex_jpkp = fEx_jpkp(i,1,k  )
        Ex_jmkp = fEx_jpkp(i,0,k  )
        Ex_jpkm = fEx_jpkp(i,1,k-1)
        Ex_jmkm = fEx_jpkp(i,0,k-1)

        Ex_kp = 0.5*(Ex_jpkp + Ex_jmkp)
        Ex_km = 0.5*(Ex_jpkm + Ex_jmkm)

        Ez_ip = 0.5*(Ez_ipjp + Ez_ipjm)*cos(yy(jg))

cc        ff(IBY) = jac*( dS3*( Ex_kp - Ex_km )
cc     .                - dS1*( Ez_ip*cos(yy(jg))) )
      endif

      call imposeBConfluxes (i,j,k,Ez_ip,Ez_im,flxjp,flxjm
     .                      ,Ex_kp,Ex_km,varray%array_var(IBY)%bconds)

      ff(IBY) = jac*( dS3*( Ex_kp - Ex_km ) - dS1*( Ez_ip - Ez_im ) )

c     Bz

      if (.not.sing_point) then
        Ey_ip = 0.5*(Ey_ipkp + Ey_ipkm)
        Ey_im = 0.5*(Ey_imkp + Ey_imkm)

        Ex_jp = 0.5*(Ex_jpkp + Ex_jpkm)
        Ex_jm = 0.5*(Ex_jmkp + Ex_jmkm)

cc        ff(IBZ) = jac*( dS1*( Ey_ip - Ey_im ) - dS2*( Ex_jp - Ex_jm ) )
      else
        Ey_ip = 0.5*(Ey_ipkp + Ey_ipkm)

cc        ff(IBZ) = jac*( dS1*( Ey_ip )  )
      endif

      call imposeBConfluxes (i,j,k,Ey_ip,Ey_im,Ex_jp,Ex_jm
     .                      ,flxkp,flxkm,varray%array_var(IBZ)%bconds)

      ff(IBZ) = jac*( dS1*( Ey_ip - Ey_im ) - dS2*( Ex_jp - Ex_jm ) )

c     Temperature

      heat_flx = -chi*laplacian(i,j,k,tmp)

      gsub   = g_sub  (xx(ig),yy(jg),zz(kg))
      gsuper = g_super(xx(ig),yy(jg),zz(kg))

      !Joule heating
      cnv = (/ jx(i,j,k), jy(i,j,k), jz(i,j,k) /)
      cov = matmul(gsub,cnv)
      joule = dxh(ig)*dyh(jg)*dzh(kg)*eeta(i,j,k)*sum(cov*cnv)
cc      joule = jouleHeating(i,j,k)


      flxip = vx(i,j,k)*(tmp(i+1,j,k)+tmp(i,j,k))/2.
     .              +(gamma-1.)*tmp(i,j,k)*(vx(i+1,j,k)+vx(i,j,k))/2.
      flxim = vx(i,j,k)*(tmp(i-1,j,k)+tmp(i,j,k))/2.
     .              +(gamma-1.)*tmp(i,j,k)*(vx(i-1,j,k)+vx(i,j,k))/2.

      flxjp = vy(i,j,k)*(tmp(i,j+1,k)+tmp(i,j,k))/2.
     .              +(gamma-1.)*tmp(i,j,k)*(vy(i,j+1,k)+vy(i,j,k))/2.
      flxjm = vy(i,j,k)*(tmp(i,j-1,k)+tmp(i,j,k))/2.
     .              +(gamma-1.)*tmp(i,j,k)*(vy(i,j-1,k)+vy(i,j,k))/2.

      flxkp = vz(i,j,k)*(tmp(i,j,k+1)+tmp(i,j,k))/2.
     .              +(gamma-1.)*tmp(i,j,k)*(vz(i,j,k+1)+vz(i,j,k))/2.
      flxkm = vz(i,j,k)*(tmp(i,j,k-1)+tmp(i,j,k))/2.
     .              +(gamma-1.)*tmp(i,j,k)*(vz(i,j,k-1)+vz(i,j,k))/2.

      if (.not.sing_point) then

        !Viscous heating
        nabla_v = fnabla_v(i,j,k,0)
        cov_tnsr =           matmul(gsub  ,nabla_v)
        cnv_tnsr = transpose(matmul(gsuper,nabla_v))
        viscous = dxh(ig)*dyh(jg)*dzh(kg)
     .           *nuu(i,j,k)/jac*sum(cov_tnsr*cnv_tnsr)

        heat_src = joule/rho(i,j,k) + viscous

cc        ff(ITMP) = dS1*(flxip-flxim)
cc     .            +dS2*(flxjp-flxjm)
cc     .            +dS3*(flxkp-flxkm)
cc     .            +(gamma-1.)*heat_flx
cc     .            -(gamma-1.)*heat_src

      else

        !Viscous heating
        nabla_v = fnabla_v(i,j,k,1)
        cov_tnsr =           matmul(gsub  ,nabla_v)
        cnv_tnsr = transpose(matmul(gsuper,nabla_v))
        viscous = dxh(ig)*dyh(jg)*dzh(kg)
     .           *nuu(i,j,k)/jac*sum(cov_tnsr*cnv_tnsr)

        heat_src = joule/rho(i,j,k) + viscous

cc        ff(ITMP) = dS1*(flxip)
cc     .            +dS3*(flxkp-flxkm)
cc     .            +(gamma-1.)*heat_flx
cc     .            -(gamma-1.)*heat_src

      endif

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(ITMP)%bconds)

      ff(ITMP) = dS1*(flxip-flxim)
     .          +dS2*(flxjp-flxjm)
     .          +dS3*(flxkp-flxkm)
     .          +(gamma-1.)*heat_flx
     .          -(gamma-1.)*heat_src

c     Vx

      call vtensor_x(i  ,j,k,t11p,t12p,t13p)
      call vtensor_x(i-1,j,k,t11m,t12m,t13m)

      call vtensor_y(i,j  ,k,t21p,t22p,t23p)
      call vtensor_y(i,j-1,k,t21m,t22m,t23m)

      call vtensor_z(i,j,k  ,t31p,t32p,t33p)
      call vtensor_z(i,j,k-1,t31m,t32m,t33m)

      cov   = covariantVector(1,xx(ig),yy(ig),zz(ig))

      if (.not.sing_point) then
        flxip = vflx_x(i  ,j,k,t11p,t12p,t13p,cov)
        flxim = vflx_x(i-1,j,k,t11m,t12m,t13m,cov)

        flxjp = vflx_y(i,j  ,k,t21p,t22p,t23p,cov)
        flxjm = vflx_y(i,j-1,k,t21m,t22m,t23m,cov)

        flxkp = vflx_z(i,j,k  ,t31p,t32p,t33p,cov)
        flxkm = vflx_z(i,j,k-1,t31m,t32m,t33m,cov)
      else
        flxip = vflx_x(i,1,k,t11p,t12p,t13p,cov)
        flxim = 0d0

        flxjp = 0d0
        flxjm = 0d0

        flxkp = vflx_z(i,1,k  ,t31p,t32p,t33p,cov)
        flxkm = vflx_z(i,1,k-1,t31m,t32m,t33m,cov)
      endif
      
      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IVX)%bconds)

      ff(IVX) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm) )

c     Vy

      cov   = covariantVector(2,xx(ig),yy(ig),zz(ig))

      if (.not.sing_point) then
        flxip = vflx_x(i  ,j,k,t11p,t12p,t13p,cov)
        flxim = vflx_x(i-1,j,k,t11m,t12m,t13m,cov)

        flxjp = vflx_y(i,j  ,k,t21p,t22p,t23p,cov)
        flxjm = vflx_y(i,j-1,k,t21m,t22m,t23m,cov)

        flxkp = vflx_z(i,j,k  ,t31p,t32p,t33p,cov)
        flxkm = vflx_z(i,j,k-1,t31m,t32m,t33m,cov)
      else
        flxip = vflx_x(i,1,k,t11p,t12p,t13p,cov)
        flxim = 0d0

        flxjp = 0d0
        flxjm = 0d0

        flxkp = vflx_z(i,1,k  ,t31p,t32p,t33p,cov)
        flxkm = vflx_z(i,1,k-1,t31m,t32m,t33m,cov)
      endif

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IVY)%bconds)

      ff(IVY) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm) )

c     Vz

      cov   = covariantVector(3,xx(ig),yy(ig),zz(ig))

      if (.not.sing_point) then
        flxip = vflx_x(i  ,j,k,t11p,t12p,t13p,cov)
        flxim = vflx_x(i-1,j,k,t11m,t12m,t13m,cov)

        flxjp = vflx_y(i,j  ,k,t21p,t22p,t23p,cov)
        flxjm = vflx_y(i,j-1,k,t21m,t22m,t23m,cov)

        flxkp = vflx_z(i,j,k  ,t31p,t32p,t33p,cov)
        flxkm = vflx_z(i,j,k-1,t31m,t32m,t33m,cov)
      else
        flxip = vflx_x(i,1,k,t11p,t12p,t13p,cov)
        flxim = 0d0

        flxjp = 0d0
        flxjm = 0d0

        flxkp = vflx_z(i,1,k  ,t31p,t32p,t33p,cov)
        flxkm = vflx_z(i,1,k-1,t31m,t32m,t33m,cov)
      endif


      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IVZ)%bconds)

      ff(IVZ) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm) )

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
