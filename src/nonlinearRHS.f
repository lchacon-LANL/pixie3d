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

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,dummy

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

c     Grid parameters

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dS1 = dyh(jg)*dzh(kg)
      dS2 = dxh(ig)*dzh(kg)
      dS3 = dxh(ig)*dyh(jg)

      gsub   = g_sub  (xx(ig),yy(jg),zz(kg))
      gsuper = g_super(xx(ig),yy(jg),zz(kg))

      jac = jacobian(xx(ig),yy(jg),zz(kg))

      sing_point = .false.

      if (jac == 0d0) then
        jac = jacobian(dxh(ig)/2.,yy(jg),zz(kg)) !To get volume right
                                                 !at singular point
        sing_point = .true.
      endif

c     Rho

      flxip = 0.5*(vx(i+1,j,k)*rho(i,j,k) + vx(i,j,k)*rho(i+1,j,k))
      flxim = 0.5*(vx(i-1,j,k)*rho(i,j,k) + vx(i,j,k)*rho(i-1,j,k))

      flxjp = 0.5*(vy(i,j+1,k)*rho(i,j,k) + vy(i,j,k)*rho(i,j+1,k))
      flxjm = 0.5*(vy(i,j-1,k)*rho(i,j,k) + vy(i,j,k)*rho(i,j-1,k))

      flxkp = 0.5*(vz(i,j,k+1)*rho(i,j,k) + vz(i,j,k)*rho(i,j,k+1))
      flxkm = 0.5*(vz(i,j,k-1)*rho(i,j,k) + vz(i,j,k)*rho(i,j,k-1))

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

cc      call imposeBConfluxes (i,j,k,Ez_ip,flxim,Ez_jp,Ez_jm
cc     .                      ,Ey_kp,Ey_km,varray%array_var(IBX)%bconds)

cc      ff(IBX) = jac*(  dS1*( Ez_ip )
cc     .               + dS2*( Ez_jp - Ez_jm ) - dS3*( Ey_kp - Ey_km ) )

      flxip = 0d0
      flxim = 0d0

      !Resistive contribution
      Ez_jp = 0.5*(eeta(i,j+1,k)*jz_cov(i,j+1,k)
     .           + eeta(i,j  ,k)*jz_cov(i,j  ,k))
      Ez_jm = 0.5*(eeta(i,j-1,k)*jz_cov(i,j-1,k)
     .           + eeta(i,j  ,k)*jz_cov(i,j  ,k))

      Ey_kp = 0.5*(eeta(i,j,k+1)*jy_cov(i,j,k+1)
     .           + eeta(i,j,k  )*jy_cov(i,j,k  ))
      Ey_km = 0.5*(eeta(i,j,k-1)*jy_cov(i,j,k-1)
     .           + eeta(i,j,k  )*jy_cov(i,j,k  ))

cc      Ez_jp = 2./(1./eeta(i,j+1,k)+1./eeta(i,j,k))*fj_jp(i,j  ,k,3)
cc      Ez_jm = 2./(1./eeta(i,j-1,k)+1./eeta(i,j,k))*fj_jp(i,j-1,k,3)
cc
cc      Ey_kp = 2./(1./eeta(i,j,k+1)+1./eeta(i,j,k))*fj_kp(i,j,k  ,2)
cc      Ey_km = 2./(1./eeta(i,j,k-1)+1./eeta(i,j,k))*fj_kp(i,j,k-1,2)
cc
cc      Ez_jp = 0d0
cc      Ez_jm = 0d0
cc
cc      Ey_kp = 0d0
cc      Ey_km = 0d0

      !Ideal contribution
      flxjp = Ez_jp
     .       +0.5*(vy(i,j+1,k)*bx(i,j+1,k) + vy(i,j,k)*bx(i,j,k))
     .       -0.5*(vx(i,j+1,k)*by(i,j+1,k) + vx(i,j,k)*by(i,j,k))
      flxjm = Ez_jm
     .       +0.5*(vy(i,j-1,k)*bx(i,j-1,k) + vy(i,j,k)*bx(i,j,k))
     .       -0.5*(vx(i,j-1,k)*by(i,j-1,k) + vx(i,j,k)*by(i,j,k))

      flxkp =-Ey_kp
     .       +0.5*(vz(i,j,k+1)*bx(i,j,k+1) + vz(i,j,k)*bx(i,j,k))
     .       -0.5*(vx(i,j,k+1)*bz(i,j,k+1) + vx(i,j,k)*bz(i,j,k))
      flxkm =-Ey_km
     .       +0.5*(vz(i,j,k-1)*bx(i,j,k-1) + vz(i,j,k)*bx(i,j,k))
     .       -0.5*(vx(i,j,k-1)*bz(i,j,k-1) + vx(i,j,k)*bz(i,j,k))

cc      flxjp = Ez_jp
cc     .       +0.5*(vy(i,j+1,k)*bx(i,j,k) + vy(i,j,k)*bx(i,j+1,k))
cc     .       -0.5*(vx(i,j+1,k)*by(i,j,k) + vx(i,j,k)*by(i,j+1,k))
cc      flxjm = Ez_jm
cc     .       +0.5*(vy(i,j-1,k)*bx(i,j,k) + vy(i,j,k)*bx(i,j-1,k))
cc     .       -0.5*(vx(i,j-1,k)*by(i,j,k) + vx(i,j,k)*by(i,j-1,k))
cc
cc      flxkp =-Ey_kp
cc     .       +0.5*(vz(i,j,k+1)*bx(i,j,k) + vz(i,j,k)*bx(i,j,k+1))
cc     .       -0.5*(vx(i,j,k+1)*bz(i,j,k) + vx(i,j,k)*bz(i,j,k+1))
cc      flxkm =-Ey_km
cc     .       +0.5*(vz(i,j,k-1)*bx(i,j,k) + vz(i,j,k)*bx(i,j,k-1))
cc     .       -0.5*(vx(i,j,k-1)*bz(i,j,k) + vx(i,j,k)*bz(i,j,k-1))

cc      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
cc     .                      ,flxkp,flxkm,varray%array_var(IBX)%bconds)

      ff(IBX) = jac*( dS1*(flxip-flxim)
     .               +dS2*(flxjp-flxjm)
     .               +dS3*(flxkp-flxkm) )
cc     .         -jac*gradivB(i,j,k,1)
cc     .         +jac*curlcurl2(i,j,k,bx_cov,by_cov,bz_cov,1)
cc     .         +jac*curlcurl(i,j,k,bx_cov,by_cov,bz_cov,1)
     .         -kdiv*laplacian(i,j,k,bx)
cc     .         -eta*laplacian(i,j,k,bx)

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

cc      call imposeBConfluxes (i,j,k,Ez_ip,Ez_im,flxjp,flxjm
cc     .                      ,Ex_kp,Ex_km,varray%array_var(IBY)%bconds)

cc      ff(IBY) = jac*( dS3*( Ex_kp - Ex_km ) - dS1*( Ez_ip - Ez_im ) )

      flxjp = 0d0
      flxjm = 0d0

      !Resistive contribution
      Ez_ip = 0.5*(eeta(i+1,j,k)*jz_cov(i+1,j,k)
     .           + eeta(i  ,j,k)*jz_cov(i  ,j,k))
      Ez_im = 0.5*(eeta(i-1,j,k)*jz_cov(i-1,j,k)
     .           + eeta(i  ,j,k)*jz_cov(i  ,j,k))

      Ex_kp = 0.5*(eeta(i,j,k+1)*jx_cov(i,j,k+1)
     .           + eeta(i,j,k  )*jx_cov(i,j,k  ))
      Ex_km = 0.5*(eeta(i,j,k-1)*jx_cov(i,j,k-1)
     .           + eeta(i,j,k  )*jx_cov(i,j,k  ))

cc      Ez_ip = 2./(1./eeta(i+1,j,k)+1./eeta(i,j,k))*fj_ip(i  ,j,k,3)
cc      Ez_im = 2./(1./eeta(i-1,j,k)+1./eeta(i,j,k))*fj_ip(i-1,j,k,3)
cc
cc      Ex_kp = 2./(1./eeta(i,j,k+1)+1./eeta(i,j,k))*fj_kp(i,j,k  ,1)
cc      Ex_km = 2./(1./eeta(i,j,k-1)+1./eeta(i,j,k))*fj_kp(i,j,k-1,1)
cc
cc      Ez_ip = 0d0
cc      Ez_im = 0d0
cc
cc      Ex_kp = 0d0
cc      Ex_km = 0d0

      !Ideal contribution
      flxip =-Ez_ip
     .       +0.5*(vx(i+1,j,k)*by(i+1,j,k) + vx(i,j,k)*by(i,j,k))
     .       -0.5*(vy(i+1,j,k)*bx(i+1,j,k) + vy(i,j,k)*bx(i,j,k))
      flxim =-Ez_im
     .       +0.5*(vx(i-1,j,k)*by(i-1,j,k) + vx(i,j,k)*by(i,j,k))
     .       -0.5*(vy(i-1,j,k)*bx(i-1,j,k) + vy(i,j,k)*bx(i,j,k))

      flxkp = Ex_kp
     .       +0.5*(vz(i,j,k+1)*by(i,j,k+1) + vz(i,j,k)*by(i,j,k))
     .       -0.5*(vy(i,j,k+1)*bz(i,j,k+1) + vy(i,j,k)*bz(i,j,k))
      flxkm = Ex_km
     .       +0.5*(vz(i,j,k-1)*by(i,j,k-1) + vz(i,j,k)*by(i,j,k))
     .       -0.5*(vy(i,j,k-1)*bz(i,j,k-1) + vy(i,j,k)*bz(i,j,k))

cc      flxip = -Ez_ip
cc     .       +0.5*(vx(i+1,j,k)*by(i,j,k) + vx(i,j,k)*by(i+1,j,k))
cc     .       -0.5*(vy(i+1,j,k)*bx(i,j,k) + vy(i,j,k)*bx(i+1,j,k))
cc      flxim = -Ez_im
cc     .       +0.5*(vx(i-1,j,k)*by(i,j,k) + vx(i,j,k)*by(i-1,j,k))
cc     .       -0.5*(vy(i-1,j,k)*bx(i,j,k) + vy(i,j,k)*bx(i-1,j,k))
cc
cc      flxkp = Ex_kp
cc     .       +0.5*(vz(i,j,k+1)*by(i,j,k) + vz(i,j,k)*by(i,j,k+1))
cc     .       -0.5*(vy(i,j,k+1)*bz(i,j,k) + vy(i,j,k)*bz(i,j,k+1))
cc      flxkm = Ex_km
cc     .       +0.5*(vz(i,j,k-1)*by(i,j,k) + vz(i,j,k)*by(i,j,k-1))
cc     .       -0.5*(vy(i,j,k-1)*bz(i,j,k) + vy(i,j,k)*bz(i,j,k-1))

cc      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
cc     .                      ,flxkp,flxkm,varray%array_var(IBY)%bconds)

      ff(IBY) = jac*( dS1*(flxip-flxim)
     .               +dS2*(flxjp-flxjm)
     .               +dS3*(flxkp-flxkm) )
cc     .         -jac*gradivB(i,j,k,2)
cc     .         +jac*curlcurl2(i,j,k,bx_cov,by_cov,bz_cov,2)
cc     .         +jac*curlcurl(i,j,k,bx_cov,by_cov,bz_cov,2)
     .         -kdiv*laplacian(i,j,k,by)
cc     .         -eta*laplacian(i,j,k,by)

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

cc      call imposeBConfluxes (i,j,k,Ey_ip,Ey_im,Ex_jp,Ex_jm
cc     .                      ,flxkp,flxkm,varray%array_var(IBZ)%bconds)

cc      ff(IBZ) = jac*( dS1*( Ey_ip - Ey_im ) - dS2*( Ex_jp - Ex_jm ) )

      flxkp = 0d0
      flxkm = 0d0

      !Resistive contribution
      Ey_ip = 0.5*(eeta(i+1,j,k)*jy_cov(i+1,j,k)
     .           + eeta(i  ,j,k)*jy_cov(i  ,j,k))
      Ey_im = 0.5*(eeta(i-1,j,k)*jy_cov(i-1,j,k)
     .           + eeta(i  ,j,k)*jy_cov(i  ,j,k))

      Ex_jp = 0.5*(eeta(i,j+1,k)*jx_cov(i,j+1,k)
     .           + eeta(i,j  ,k)*jx_cov(i,j  ,k))
      Ex_jm = 0.5*(eeta(i,j-1,k)*jx_cov(i,j-1,k)
     .           + eeta(i,j  ,k)*jx_cov(i,j  ,k))

cc      Ey_ip = 2./(1./eeta(i+1,j,k)+1./eeta(i,j,k))*fj_ip(i  ,j,k,2)
cc      Ey_im = 2./(1./eeta(i-1,j,k)+1./eeta(i,j,k))*fj_ip(i-1,j,k,2)
cc
cc      Ex_jp = 2./(1./eeta(i,j+1,k)+1./eeta(i,j,k))*fj_jp(i,j  ,k,1)
cc      Ex_jm = 2./(1./eeta(i,j-1,k)+1./eeta(i,j,k))*fj_jp(i,j-1,k,1)
cc
cc      Ey_ip = 0d0
cc      Ey_im = 0d0
cc
cc      Ex_jp = 0d0
cc      Ex_jm = 0d0

      !Ideal contribution
      flxip = Ey_ip
     .       +0.5*(vx(i+1,j,k)*bz(i+1,j,k) + vx(i,j,k)*bz(i,j,k))
     .       -0.5*(vz(i+1,j,k)*bx(i+1,j,k) + vz(i,j,k)*bx(i,j,k))
      flxim = Ey_im
     .       +0.5*(vx(i-1,j,k)*bz(i-1,j,k) + vx(i,j,k)*bz(i,j,k))
     .       -0.5*(vz(i-1,j,k)*bx(i-1,j,k) + vz(i,j,k)*bx(i,j,k))

      flxjp =-Ex_jp
     .       +0.5*(vy(i,j+1,k)*bz(i,j+1,k) + vy(i,j,k)*bz(i,j,k))
     .       -0.5*(vz(i,j+1,k)*by(i,j+1,k) + vz(i,j,k)*by(i,j,k))
      flxjm =-Ex_jm
     .       +0.5*(vy(i,j-1,k)*bz(i,j-1,k) + vy(i,j,k)*bz(i,j,k))
     .       -0.5*(vz(i,j-1,k)*by(i,j-1,k) + vz(i,j,k)*by(i,j,k))

cc      flxip = Ey_ip
cc     .       +0.5*(vx(i+1,j,k)*bz(i,j,k) + vx(i,j,k)*bz(i+1,j,k))
cc     .       -0.5*(vz(i+1,j,k)*bx(i,j,k) + vz(i,j,k)*bx(i+1,j,k))
cc      flxim = Ey_im
cc     .       +0.5*(vx(i-1,j,k)*bz(i,j,k) + vx(i,j,k)*bz(i-1,j,k))
cc     .       -0.5*(vz(i-1,j,k)*bx(i,j,k) + vz(i,j,k)*bx(i-1,j,k))
cc
cc      flxjp =-Ex_jp
cc     .       +0.5*(vy(i,j+1,k)*bz(i,j,k) + vy(i,j,k)*bz(i,j+1,k))
cc     .       -0.5*(vz(i,j+1,k)*by(i,j,k) + vz(i,j,k)*by(i,j+1,k))
cc      flxjm =-Ex_jm
cc     .       +0.5*(vy(i,j-1,k)*bz(i,j,k) + vy(i,j,k)*bz(i,j-1,k))
cc     .       -0.5*(vz(i,j-1,k)*by(i,j,k) + vz(i,j,k)*by(i,j-1,k))

cc      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
cc     .                      ,flxkp,flxkm,varray%array_var(IBZ)%bconds)

      ff(IBZ) = jac*( dS1*(flxip-flxim)
     .               +dS2*(flxjp-flxjm)
     .               +dS3*(flxkp-flxkm) )
cc     .         -jac*gradivB(i,j,k,3)
cc     .         +jac*curlcurl2(i,j,k,bx_cov,by_cov,bz_cov,3)
cc     .         +jac*curlcurl(i,j,k,bx_cov,by_cov,bz_cov,3)
     .         -kdiv*laplacian(i,j,k,bz)
cc     .        -eta*laplacian(i,j,k,bz)

c     Temperature

      flxip = (vx(i+1,j,k)*tmp(i,j,k) + vx(i,j,k)*tmp(i+1,j,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vx(i+1,j,k)+vx(i,j,k))/2.
      flxim = (vx(i-1,j,k)*tmp(i,j,k) + vx(i,j,k)*tmp(i-1,j,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vx(i-1,j,k)+vx(i,j,k))/2.

      flxjp = (vy(i,j+1,k)*tmp(i,j,k) + vy(i,j,k)*tmp(i,j+1,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vy(i,j+1,k)+vy(i,j,k))/2.
      flxjm = (vy(i,j-1,k)*tmp(i,j,k) + vy(i,j,k)*tmp(i,j-1,k))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vy(i,j-1,k)+vy(i,j,k))/2.

      flxkp = (vz(i,j,k+1)*tmp(i,j,k) + vz(i,j,k)*tmp(i,j,k+1))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vz(i,j,k+1)+vz(i,j,k))/2.
      flxkm = (vz(i,j,k-1)*tmp(i,j,k) + vz(i,j,k)*tmp(i,j,k-1))/2.
     .              +(gamma-2.)*tmp(i,j,k)*(vz(i,j,k-1)+vz(i,j,k))/2.

      !Heat flux
      heat_flx = -chi*laplacian(i,j,k,tmp)

      !Joule heating
      joule = dxh(ig)*dyh(jg)*dzh(kg)*eeta(i,j,k)
     .                               *( jx(i,j,k)*jx_cov(i,j,k)
     .                                 +jy(i,j,k)*jy_cov(i,j,k)
     .                                 +jz(i,j,k)*jz_cov(i,j,k) )
cc      joule = jouleHeating(i,j,k)

      !Viscous heating
      if (sing_point) then
        nabla_v = fnabla_v(i,j,k,1)
      else
        nabla_v = fnabla_v(i,j,k,0)
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

      call vtensor_x(i  ,j,k,t11p,t12p,t13p)
      call vtensor_x(i-1,j,k,t11m,t12m,t13m)

      call vtensor_y(i,j  ,k,t21p,t22p,t23p)
      call vtensor_y(i,j-1,k,t21m,t22m,t23m)

      call vtensor_z(i,j,k  ,t31p,t32p,t33p)
      call vtensor_z(i,j,k-1,t31m,t32m,t33m)

      cov   = covariantVector(1,xx(ig),yy(ig),zz(ig))

      flxip = vflx_x(i  ,j,k,t11p,t12p,t13p,cov)
      flxim = vflx_x(i-1,j,k,t11m,t12m,t13m,cov)

      flxjp = vflx_y(i,j  ,k,t21p,t22p,t23p,cov)
      flxjm = vflx_y(i,j-1,k,t21m,t22m,t23m,cov)

      flxkp = vflx_z(i,j,k  ,t31p,t32p,t33p,cov)
      flxkm = vflx_z(i,j,k-1,t31m,t32m,t33m,cov)

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IVX)%bconds)

      ff(IVX) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm) )

c     Vy

      cov   = covariantVector(2,xx(ig),yy(ig),zz(ig))

      flxip = vflx_x(i  ,j,k,t11p,t12p,t13p,cov)
      flxim = vflx_x(i-1,j,k,t11m,t12m,t13m,cov)

      flxjp = vflx_y(i,j  ,k,t21p,t22p,t23p,cov)
      flxjm = vflx_y(i,j-1,k,t21m,t22m,t23m,cov)

      flxkp = vflx_z(i,j,k  ,t31p,t32p,t33p,cov)
      flxkm = vflx_z(i,j,k-1,t31m,t32m,t33m,cov)

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                      ,flxkp,flxkm,varray%array_var(IVY)%bconds)

      ff(IVY) = jac*( dS1*(flxip - flxim)
     .              + dS2*(flxjp - flxjm)
     .              + dS3*(flxkp - flxkm) )

c     Vz

      cov   = covariantVector(3,xx(ig),yy(ig),zz(ig))

      flxip = vflx_x(i  ,j,k,t11p,t12p,t13p,cov)
      flxim = vflx_x(i-1,j,k,t11m,t12m,t13m,cov)

      flxjp = vflx_y(i,j  ,k,t21p,t22p,t23p,cov)
      flxjm = vflx_y(i,j-1,k,t21m,t22m,t23m,cov)

      flxkp = vflx_z(i,j,k  ,t31p,t32p,t33p,cov)
      flxkm = vflx_z(i,j,k-1,t31m,t32m,t33m,cov)

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
