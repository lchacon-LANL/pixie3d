c module equilibrium
c ######################################################################
      module equilibrium

        integer(4)       :: IRHO,IVX,IVY,IVZ,IBX,IBY,IBZ,ITMP
        parameter (IRHO=1,IVX=2,IVY=3,IVZ=4,IBX=5,IBY=6,IBZ=7,ITMP=8)

        real(8) :: dlambda,rshear,vparflow,vperflow

        character*(5)    :: equil

      end module equilibrium

c module transport_params
c ######################################################################
      module transport_params

        real(8) :: nu,eta,dd,chi,gamma

      end module transport_params

c module nlfunction_setup
c ######################################################################
      module nlfunction_setup

        use grid

        use parameters

        use transport_params

        use equilibrium

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_cov,by_cov,bz_cov,eeta,nuu,jx,jy,jz,divrgB,vx,vy,vz

        real(8),pointer,dimension(:,:,:):: rho,rvx,rvy,rvz,bx,by,bz,tmp
        real(8),pointer,dimension(:)    :: xx,yy,zz,dxh,dyh,dzh,dx,dy,dz

        integer(4) :: igx,igy,igz,nx,ny,nz

        real(8)    :: gsub(3,3),gsuper(3,3),cnv(3),cov(3)

        real(8)    :: nabla_v(3,3),hessian1(3,3)
     .               ,hessian2(3,3),hessian3(3,3)
     .               ,cov_tnsr(3,3),cnv_tnsr(3,3)

        logical    :: sing_point

      contains

c     fjz_ipjp
c     #############################################################
      real*8 function fjz_ipjp(i,j,k)
c     -------------------------------------------------------------
c     Calculates z-component of covariant current at (i+1/2,j+1/2,k)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: jx,jy,jz,x,y,z,gsub(3,3)

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x = 0.5*(xx(ig)+xx(ig+1))
        y = 0.5*(yy(jg)+yy(jg+1))
        z = zz(kg)

        gsub = g_sub(x,y,z)

        jx = (bz_cov(i  ,j+1,k) - bz_cov(i  ,j,k)
     .       +bz_cov(i+1,j+1,k) - bz_cov(i+1,j,k) )/2./dy(jg)
     .      -(by_cov(i,j+1,k+1) + by_cov(i  ,j  ,k+1)
     .       +by_cov(i+1,j,k+1) + by_cov(i+1,j+1,k+1)
     .       -by_cov(i,j+1,k-1) - by_cov(i  ,j  ,k-1)
     .       -by_cov(i+1,j,k-1) - by_cov(i+1,j+1,k-1))/8./dzh(kg)

        jy = (bx_cov(i,j+1,k+1) + bx_cov(i  ,j  ,k+1)
     .       +bx_cov(i+1,j,k+1) + bx_cov(i+1,j+1,k+1)
     .       -bx_cov(i,j+1,k-1) - bx_cov(i  ,j  ,k-1)
     .       -bx_cov(i+1,j,k-1) - bx_cov(i+1,j+1,k-1))/8./dzh(kg)
     .      -(bz_cov(i+1,j  ,k) - bz_cov(i,j  ,k)
     .       +bz_cov(i+1,j+1,k) - bz_cov(i,j+1,k) )/2./dx(ig)

        jz = (by_cov(i+1,j  ,k) - by_cov(i,j  ,k)
     .       +by_cov(i+1,j+1,k) - by_cov(i,j+1,k) )/2./dx(ig)
     .      -(bx_cov(i  ,j+1,k) - bx_cov(i  ,j,k)
     .       +bx_cov(i+1,j+1,k) - bx_cov(i+1,j,k) )/2./dy(jg)

        fjz_ipjp = gsub(3,1)*jx + gsub(3,2)*jy + gsub(3,3)*jz

c     End program

      end function fjz_ipjp

c     fjy_ipkp
c     #############################################################
      real*8 function fjy_ipkp(i,j,k)
c     -------------------------------------------------------------
c     Calculates y-component of covariant current at (i+1/2,j,k+1/2)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: jx,jy,jz,x,y,z,gsub(3,3)

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x = 0.5*(xx(ig)+xx(ig+1))
        y =      yy(jg)
        z = 0.5*(zz(kg)+zz(kg+1))

        gsub = g_sub(x,y,z)

        jx = (bz_cov(i  ,j+1,k+1) + bz_cov(i  ,j+1,k  )
     .       +bz_cov(i+1,j+1,k  ) + bz_cov(i+1,j+1,k+1)
     .       -bz_cov(i  ,j-1,k+1) - bz_cov(i  ,j-1,k  )
     .       -bz_cov(i+1,j-1,k  ) - bz_cov(i+1,j-1,k+1))/8./dyh(jg)
     .      -(by_cov(i  ,j,k+1) - by_cov(i  ,j,k)
     .       +by_cov(i+1,j,k+1) - by_cov(i+1,j,k) )/2./dz(kg)

        jy = (bx_cov(i  ,j,k+1) - bx_cov(i  ,j,k)
     .       +bx_cov(i+1,j,k+1) - bx_cov(i+1,j,k) )/2./dz(kg)
     .      -(bz_cov(i+1,j,k  ) - bz_cov(i,j,k  )
     .       +bz_cov(i+1,j,k+1) - bz_cov(i,j,k+1) )/2./dx(ig)

        jz = (by_cov(i+1,j,k  ) - by_cov(i  ,j,k  )
     .       +by_cov(i+1,j,k+1) - by_cov(i  ,j,k+1) )/2./dx(ig)
     .      -(bx_cov(i  ,j+1,k+1) + bx_cov(i  ,j+1,k  )
     .       +bx_cov(i+1,j+1,k  ) + bx_cov(i+1,j+1,k+1)
     .       -bx_cov(i  ,j-1,k+1) - bx_cov(i  ,j-1,k  )
     .       -bx_cov(i+1,j-1,k  ) - bx_cov(i+1,j-1,k+1))/8./dyh(jg)

        fjy_ipkp = gsub(2,1)*jx + gsub(2,2)*jy + gsub(2,3)*jz

c     End program

      end function fjy_ipkp

c     fjx_jpkp
c     #############################################################
      real*8 function fjx_jpkp(i,j,k)
c     -------------------------------------------------------------
c     Calculates x-component of covariant current at (i,j+1/2,k+1/2)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: jx,jy,jz,x,y,z,gsub(3,3)

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x = 0.5*(xx(ig)+xx(ig+1))
        y =      yy(jg)
        z = 0.5*(zz(kg)+zz(kg+1))

        gsub = g_sub(x,y,z)

        jx = (bz_cov(i,j+1,k  ) - bz_cov(i,j,k  )
     .       +bz_cov(i,j+1,k+1) - bz_cov(i,j,k+1) )/2./dy(jg)
     .      -(by_cov(i,j  ,k+1) - by_cov(i,j  ,k)
     .       +by_cov(i,j+1,k+1) - by_cov(i,j+1,k) )/2./dz(kg)

        jy = (bx_cov(i,j  ,k+1) - bx_cov(i,j  ,k)
     .       +bx_cov(i,j+1,k+1) - bx_cov(i,j+1,k) )/2./dz(kg)
     .      -(bz_cov(i+1,j  ,k+1) + bz_cov(i+1,j+1,k  )
     .       +bz_cov(i+1,j+1,k+1) + bz_cov(i+1,j  ,k  )
     .       -bz_cov(i-1,j  ,k+1) + bz_cov(i-1,j+1,k  )
     .       +bz_cov(i-1,j+1,k+1) + bz_cov(i-1,j  ,k  ))/8./dxh(ig)

        jz = (by_cov(i+1,j  ,k+1) + by_cov(i+1,j+1,k  )
     .       +by_cov(i+1,j+1,k+1) + by_cov(i+1,j  ,k  )
     .       -by_cov(i-1,j  ,k+1) + by_cov(i-1,j+1,k  )
     .       +by_cov(i-1,j+1,k+1) + by_cov(i-1,j  ,k  ))/8./dxh(ig)
     .      -(bx_cov(i,j+1,k  ) - bx_cov(i,j,k  )
     .       +bx_cov(i,j+1,k+1) - bx_cov(i,j,k+1) )/2./dy(jg)

        fjx_jpkp = gsub(1,1)*jx + gsub(1,2)*jy + gsub(1,3)*jz

c     End program

      end function fjx_jpkp

c     fEz_ipjp
c     #############################################################
      real*8 function fEz_ipjp(i,j,k)
c     -------------------------------------------------------------
c     Calculates z-component of covariant MHD electric field 
c     at (i+1/2,j+1/2,k).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k
        real(8)    :: eta_ipjp,jz_ipjp

c     Local variables

c     Begin program

        eta_ipjp = 4./(1./eeta(i+1,j+1,k) + 1./eeta(i,j+1,k)
     .                +1./eeta(i+1,j  ,k) + 1./eeta(i,j  ,k))

        jz_ipjp = fjz_ipjp(i,j,k)

        fEz_ipjp = eta_ipjp*jz_ipjp
     .     -(vx(i+1,j+1,k)*by(i  ,j,k) + vx(i  ,j,k)*by(i+1,j+1,k)
     .      +vx(i  ,j+1,k)*by(i+1,j,k) + vx(i+1,j,k)*by(i  ,j+1,k))/4.
     .     +(vy(i+1,j+1,k)*bx(i  ,j,k) + vy(i  ,j,k)*bx(i+1,j+1,k)
     .      +vy(i  ,j+1,k)*bx(i+1,j,k) + vy(i+1,j,k)*bx(i  ,j+1,k))/4.

c     End program

      end function fEz_ipjp

c     fEy_ipkp
c     #############################################################
      real*8 function fEy_ipkp(i,j,k)
c     -------------------------------------------------------------
c     Calculates y-component of covariant ideal MHD electric field 
c     at (i+1/2,j,k+1/2).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k
        real(8)    :: eta_ipkp,jy_ipkp

c     Local variables

c     Begin program

        eta_ipkp = 4./(1./eeta(i+1,j,k+1) + 1./eeta(i,j,k+1)
     .                +1./eeta(i+1,j,k  ) + 1./eeta(i,j,k  ))

        jy_ipkp = fjy_ipkp(i,j,k)

        fEy_ipkp = eta_ipkp*jy_ipkp
     .     -(vz(i+1,j,k+1)*bx(i  ,j,k) + vz(i  ,j,k)*bx(i+1,j,k+1)
     .      +vz(i  ,j,k+1)*bx(i+1,j,k) + vz(i+1,j,k)*bx(i  ,j,k+1))/4.
     .     +(vx(i+1,j,k+1)*bz(i  ,j,k) + vx(i  ,j,k)*bz(i+1,j,k+1)
     .      +vx(i  ,j,k+1)*bz(i+1,j,k) + vx(i+1,j,k)*bz(i  ,j,k+1))/4.

c     End program

      end function fEy_ipkp

c     fEx_jpkp
c     #############################################################
      real*8 function fEx_jpkp(i,j,k)
c     -------------------------------------------------------------
c     Calculates y-component of covariant ideal MHD electric field 
c     at (i+1/2,j,k+1/2).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k
        real(8)    :: eta_jpkp,jx_jpkp

c     Local variables

c     Begin program

        eta_jpkp = 4./(1./eeta(i,j+1,k+1) + 1./eeta(i,j,k+1)
     .                +1./eeta(i,j+1,k  ) + 1./eeta(i,j,k  ))

        jx_jpkp = fjx_jpkp(i,j,k)

        fEx_jpkp = eta_jpkp*jx_jpkp
     .     -(vy(i,j+1,k+1)*bz(i,j,k  ) + vy(i,j,k  )*bz(i,j+1,k+1)
     .      +vy(i,j+1,k  )*bz(i,j,k+1) + vy(i,j,k+1)*bz(i,j+1,k  ))/4.
     .     +(vz(i,j+1,k+1)*by(i,j,k  ) + vz(i,j,k  )*by(i,j+1,k+1)
     .      +vz(i,j+1,k  )*by(i,j,k+1) + vz(i,j,k+1)*by(i,j+1,k  ))/4.

c     End program

      end function fEx_jpkp

c     divB
c     ###############################################################
      real(8) function divB(i,j,k)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of magnetic field at grid vertices.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: i,j,k

c     Local variables

      integer(4) :: ig,jg,kg

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      divB =
     .       (bx(i+1,j  ,k  )-bx(i,j  ,k  ))*dyh(jg  )*dzh(kg  )
     .      +(bx(i+1,j+1,k  )-bx(i,j+1,k  ))*dyh(jg+1)*dzh(kg  )
     .      +(bx(i+1,j  ,k+1)-bx(i,j  ,k+1))*dyh(jg  )*dzh(kg+1)
     .      +(bx(i+1,j+1,k+1)-bx(i,j+1,k+1))*dyh(jg+1)*dzh(kg+1)
     .      +(by(i  ,j+1,k  )-by(i  ,j,k  ))*dxh(ig  )*dzh(kg  )
     .      +(by(i+1,j+1,k  )-by(i+1,j,k  ))*dxh(ig+1)*dzh(kg  )
     .      +(by(i  ,j+1,k+1)-by(i  ,j,k+1))*dxh(ig  )*dzh(kg+1)
     .      +(by(i+1,j+1,k+1)-by(i+1,j,k+1))*dxh(ig+1)*dzh(kg+1)
     .      +(bz(i  ,j  ,k+1)-bz(i  ,j  ,k))*dyh(jg  )*dxh(ig  )
     .      +(bz(i  ,j+1,k+1)-bz(i  ,j+1,k))*dyh(jg+1)*dxh(ig  )
     .      +(bz(i+1,j  ,k+1)-bz(i+1,j  ,k))*dyh(jg  )*dxh(ig+1)
     .      +(bz(i+1,j+1,k+1)-bz(i+1,j+1,k))*dyh(jg+1)*dxh(ig+1)
      
c End 

      end function divB

c     jouleHeating
c     #############################################################
      real*8 function jouleHeating(i,j,k) result(joule)
c     -------------------------------------------------------------
c     Calculates joule heating at cell volume (i,j,k) based on
c     edge currents. Result includes cell volume.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: joule1,joule2,joule3

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        joule1 = (dx(ig  )*dy(jg  )*dzh(kg)*joule_ipjp(i  ,j  ,k)
     .           +dx(ig-1)*dy(jg  )*dzh(kg)*joule_ipjp(i-1,j  ,k)
     .           +dx(ig  )*dy(jg-1)*dzh(kg)*joule_ipjp(i  ,j-1,k)
     .           +dx(ig-1)*dy(jg-1)*dzh(kg)*joule_ipjp(i-1,j-1,k))/4.

        joule2 = (dx(ig  )*dyh(jg)*dz(kg  )*joule_ipkp(i  ,j,k  )
     .           +dx(ig-1)*dyh(jg)*dz(kg  )*joule_ipkp(i-1,j,k  )
     .           +dx(ig  )*dyh(jg)*dz(kg-1)*joule_ipkp(i  ,j,k-1)
     .           +dx(ig-1)*dyh(jg)*dz(kg-1)*joule_ipkp(i-1,j,k-1))/4.

        joule3 = (dxh(ig)*dy(jg  )*dz(kg  )*joule_jpkp(i,j  ,k  )
     .           +dxh(ig)*dy(jg-1)*dz(kg  )*joule_jpkp(i,j-1,k  )
     .           +dxh(ig)*dy(jg  )*dz(kg-1)*joule_jpkp(i,j  ,k-1)
     .           +dxh(ig)*dy(jg-1)*dz(kg-1)*joule_jpkp(i,j-1,k-1))/4.

        joule = (joule1 + joule2 + joule3)/3.

c     End program

      contains

c     joule_ipjp
c     #############################################################
      real(8) function joule_ipjp(i,j,k)

c     -------------------------------------------------------------
c     Calculates joule heating at cell volume (i+1/2,j+1/2,k) based
c     on edge currents. Result includes cell volume.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: jj(3),x,y,z,gsub(3,3),eta_ipjp

c     Begin program

        eta_ipjp = 4./(1./eeta(i+1,j+1,k) + 1./eeta(i,j+1,k)
     .                +1./eeta(i+1,j  ,k) + 1./eeta(i,j  ,k))

cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x = 0.5*(xx(ig)+xx(ig+1))
        y = 0.5*(yy(jg)+yy(jg+1))
        z = zz(kg)

        gsub = g_sub(x,y,z)

        !Contravarian current components at (ip,jp,k)

        jj(1) = (bz_cov(i  ,j+1,k) - bz_cov(i  ,j,k)
     .          +bz_cov(i+1,j+1,k) - bz_cov(i+1,j,k) )/2./dy(jg)
     .         -(by_cov(i,j+1,k+1) + by_cov(i  ,j  ,k+1)
     .          +by_cov(i+1,j,k+1) + by_cov(i+1,j+1,k+1)
     .          -by_cov(i,j+1,k-1) - by_cov(i  ,j  ,k-1)
     .          -by_cov(i+1,j,k-1) - by_cov(i+1,j+1,k-1))/8./dzh(kg)

        jj(2) = (bx_cov(i,j+1,k+1) + bx_cov(i  ,j  ,k+1)
     .          +bx_cov(i+1,j,k+1) + bx_cov(i+1,j+1,k+1)
     .          -bx_cov(i,j+1,k-1) - bx_cov(i  ,j  ,k-1)
     .          -bx_cov(i+1,j,k-1) - bx_cov(i+1,j+1,k-1))/8./dzh(kg)
     .         -(bz_cov(i+1,j  ,k) - bz_cov(i,j  ,k)
     .          +bz_cov(i+1,j+1,k) - bz_cov(i,j+1,k) )/2./dx(ig)

        jj(3) = (by_cov(i+1,j  ,k) - by_cov(i,j  ,k)
     .          +by_cov(i+1,j+1,k) - by_cov(i,j+1,k) )/2./dx(ig)
     .         -(bx_cov(i  ,j+1,k) - bx_cov(i  ,j,k)
     .          +bx_cov(i+1,j+1,k) - bx_cov(i+1,j,k) )/2./dy(jg)

        joule_ipjp = eta_ipjp*dot_product(jj,matmul(gsub,jj))

      end function joule_ipjp

c     joule_ipkp
c     #############################################################
      real(8) function joule_ipkp(i,j,k)

c     -------------------------------------------------------------
c     Calculates joule heating at cell volume (i+1/2,j,k+1/2) based
c     on edge currents.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: jj(3),x,y,z,gsub(3,3),eta_ipkp

c     Begin program

        eta_ipkp = 4./(1./eeta(i+1,j,k+1) + 1./eeta(i,j,k+1)
     .                +1./eeta(i+1,j,k  ) + 1./eeta(i,j,k  ))

cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x = 0.5*(xx(ig)+xx(ig+1))
        y =      yy(jg)
        z = 0.5*(zz(kg)+zz(kg+1))

        gsub = g_sub(x,y,z)

        jj(1) = (bz_cov(i  ,j+1,k+1) + bz_cov(i  ,j+1,k  )
     .          +bz_cov(i+1,j+1,k  ) + bz_cov(i+1,j+1,k+1)
     .          -bz_cov(i  ,j-1,k+1) - bz_cov(i  ,j-1,k  )
     .          -bz_cov(i+1,j-1,k  ) - bz_cov(i+1,j-1,k+1))/8./dyh(jg)
     .         -(by_cov(i  ,j,k+1) - by_cov(i  ,j,k)
     .          +by_cov(i+1,j,k+1) - by_cov(i+1,j,k) )/2./dz(kg)

        jj(2) = (bx_cov(i  ,j,k+1) - bx_cov(i  ,j,k)
     .          +bx_cov(i+1,j,k+1) - bx_cov(i+1,j,k) )/2./dz(kg)
     .         -(bz_cov(i+1,j,k  ) - bz_cov(i,j,k  )
     .          +bz_cov(i+1,j,k+1) - bz_cov(i,j,k+1) )/2./dx(ig)

        jj(3) = (by_cov(i+1,j,k  ) - by_cov(i  ,j,k  )
     .          +by_cov(i+1,j,k+1) - by_cov(i  ,j,k+1) )/2./dx(ig)
     .         -(bx_cov(i  ,j+1,k+1) + bx_cov(i  ,j+1,k  )
     .          +bx_cov(i+1,j+1,k  ) + bx_cov(i+1,j+1,k+1)
     .          -bx_cov(i  ,j-1,k+1) - bx_cov(i  ,j-1,k  )
     .          -bx_cov(i+1,j-1,k  ) - bx_cov(i+1,j-1,k+1))/8./dyh(jg)

        !Contravarian current components at (ip,j,kp)

        joule_ipkp = eta_ipkp*dot_product(jj,matmul(gsub,jj))

      end function joule_ipkp

c     joule_jpkp
c     #############################################################
      real(8) function joule_jpkp(i,j,k)

c     -------------------------------------------------------------
c     Calculates joule heating at cell volume (i,j+1/2,k+1/2) based
c     on edge currents.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: jj(3),x,y,z,gsub(3,3),eta_jpkp

c     Begin program

        eta_jpkp = 4./(1./eeta(i,j+1,k+1) + 1./eeta(i,j,k+1)
     .                +1./eeta(i,j+1,k  ) + 1./eeta(i,j,k  ))

cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x = 0.5*(xx(ig)+xx(ig+1))
        y =      yy(jg)
        z = 0.5*(zz(kg)+zz(kg+1))

        gsub = g_sub(x,y,z)

        jj(1) = (bz_cov(i,j+1,k  ) - bz_cov(i,j,k  )
     .          +bz_cov(i,j+1,k+1) - bz_cov(i,j,k+1) )/2./dy(jg)
     .         -(by_cov(i,j  ,k+1) - by_cov(i,j  ,k)
     .          +by_cov(i,j+1,k+1) - by_cov(i,j+1,k) )/2./dz(kg)

        jj(2) = (bx_cov(i,j  ,k+1) - bx_cov(i,j  ,k)
     .          +bx_cov(i,j+1,k+1) - bx_cov(i,j+1,k) )/2./dz(kg)
     .         -(bz_cov(i+1,j  ,k+1) + bz_cov(i+1,j+1,k  )
     .          +bz_cov(i+1,j+1,k+1) + bz_cov(i+1,j  ,k  )
     .          -bz_cov(i-1,j  ,k+1) + bz_cov(i-1,j+1,k  )
     .          +bz_cov(i-1,j+1,k+1) + bz_cov(i-1,j  ,k  ))/8./dxh(ig)

        jj(3) = (by_cov(i+1,j  ,k+1) + by_cov(i+1,j+1,k  )
     .          +by_cov(i+1,j+1,k+1) + by_cov(i+1,j  ,k  )
     .          -by_cov(i-1,j  ,k+1) + by_cov(i-1,j+1,k  )
     .          +by_cov(i-1,j+1,k+1) + by_cov(i-1,j  ,k  ))/8./dxh(ig)
     .         -(bx_cov(i,j+1,k  ) - bx_cov(i,j,k  )
     .          +bx_cov(i,j+1,k+1) - bx_cov(i,j,k+1) )/2./dy(jg)

        !Contravarian current components at (ip,j,kp)

        joule_jpkp = eta_jpkp*dot_product(jj,matmul(gsub,jj))

      end function joule_jpkp

      end function jouleHeating

c     fnabla_v
c     #############################################################
      function fnabla_v(i,j,k,half_elem) result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor nabla(vec v) at the following positions:
c       + half_elem =/ 1,2,3 => i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,half_elem
        real(8)    :: tensor(3,3)

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: jx,jy,jz,x,y,z
        real(8)    :: dhx,dhy,dhz,vxx,vyy,vzz

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        !Defaults

        x = xx(ig)
        y = yy(jg)
        z = zz(kg)

        dhx = 2.*dxh(ig)
        dhy = 2.*dyh(jg)
        dhz = 2.*dzh(kg)

        vxx = vx(i,j,k)
        vyy = vy(i,j,k)
        vzz = vz(i,j,k)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        !Exceptions

        select case(half_elem)
        case (1)
          dhx = dx(ig)
          x   = (xx(ig) + xx(ig+1))/2.
          im  = i
          vxx = 0.5*(vx(i,j,k)+vx(ip,j,k))
          vyy = 0.5*(vy(i,j,k)+vy(ip,j,k))
          vzz = 0.5*(vz(i,j,k)+vz(ip,j,k))
        case (2)
          dhy = dy(jg)
          y   = (yy(jg) + yy(jg+1))/2.
          jm  = j
          vxx = 0.5*(vx(i,j,k)+vx(i,jp,k))
          vyy = 0.5*(vy(i,j,k)+vy(i,jp,k))
          vzz = 0.5*(vz(i,j,k)+vz(i,jp,k))
        case (3)
          dhz = dz(kg)
          z   = (zz(kg) + zz(kg+1))/2.
          km  = k
          vxx = 0.5*(vx(i,j,k)+vx(i,j,kp))
          vyy = 0.5*(vy(i,j,k)+vy(i,j,kp))
          vzz = 0.5*(vz(i,j,k)+vz(i,j,kp))
        end select

c     !Calculate nabla_v tensor

        hessian1 = hessian(1,x,y,z)
        hessian2 = hessian(2,x,y,z)
        hessian3 = hessian(3,x,y,z)

      ! l = 1, m = 1
        tensor(1,1) = (vx(ip,j,k)-vx(im,j,k))/dhx
     .               + vxx*(hessian1(1,1)
     .                    + hessian2(2,1) + hessian3(3,1))
     .               - vxx*hessian1(1,1)
     .               - vyy*hessian1(2,1)
     .               - vzz*hessian1(3,1)

      ! l = 1, m = 2
        tensor(1,2) = (vy(ip,j,k)-vy(im,j,k))/dhx
     .               + vyy*(hessian1(1,1)
     .                    + hessian2(2,1) + hessian3(3,1))
     .               - vxx*hessian2(1,1)
     .               - vyy*hessian2(2,1)
     .               - vzz*hessian2(3,1)

      ! l = 1, m = 3
        tensor(1,3) = (vz(ip,j,k)-vz(im,j,k))/dhx
     .               + vzz*(hessian1(1,1)
     .                    + hessian2(2,1) + hessian3(3,1))
     .               - vxx*hessian3(1,1)
     .               - vyy*hessian3(2,1)
     .               - vzz*hessian3(3,1)

      ! l = 2, m = 1
        tensor(2,1) = (vx(i,jp,k)-vx(i,jm,k))/dhy
     .               + vxx*(hessian1(1,2)
     .                    + hessian2(2,2) + hessian3(3,2))
     .               - vxx*hessian1(1,2)
     .               - vyy*hessian1(2,2)
     .               - vzz*hessian1(3,2)

      ! l = 2, m = 2
        tensor(2,2) = (vy(i,jp,k)-vy(i,jm,k))/dhy
     .               + vyy*(hessian1(1,2)
     .                    + hessian2(2,2) + hessian3(3,2))
     .               - vxx*hessian2(1,2)
     .               - vyy*hessian2(2,2)
     .               - vzz*hessian2(3,2)

      ! l = 2, m = 3
        tensor(2,3) = (vz(i,jp,k)-vz(i,jm,k))/dhy
     .               + vzz*(hessian1(1,2)
     .                    + hessian2(2,2) + hessian3(3,2))
     .               - vxx*hessian3(1,2)
     .               - vyy*hessian3(2,2)
     .               - vzz*hessian3(3,2)

      ! l = 3, m = 1
        tensor(3,1) = (vx(i,j,kp)-vx(i,j,km))/dhz
     .               + vxx*(hessian1(1,3)
     .                    + hessian2(2,3) + hessian3(3,3))
     .               - vxx*hessian1(1,3)
     .               - vyy*hessian1(2,3)
     .               - vzz*hessian1(3,3)

      ! l = 3, m = 2
        tensor(3,2) = (vy(i,j,kp)-vy(i,j,km))/dhz
     .               + vyy*(hessian1(1,3)
     .                    + hessian2(2,3) + hessian3(3,3))
     .               - vxx*hessian2(1,3)
     .               - vyy*hessian2(2,3)
     .               - vzz*hessian2(3,3)

      ! l = 3, m = 3
        tensor(3,3) = (vz(i,j,kp)-vz(i,j,km))/dhz
     .               + vzz*(hessian1(1,3)
     .                    + hessian2(2,3) + hessian3(3,3))
     .               - vxx*hessian3(1,3)
     .               - vyy*hessian3(2,3)
     .               - vzz*hessian3(3,3)

c     End program

      end function fnabla_v

c     fnabla_v_bc
c     #############################################################
      function fnabla_v_bc(i,j,k,half_elem,coeff) result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor nabla(vec v) at the following positions:
c       + half_elem =/ 1,2,3 => i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,half_elem
        real(8)    :: tensor(3,3),coeff

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: jx,jy,jz,x,y,z
        real(8)    :: dhx,dhy,dhz,vxx,vyy,vzz
        real(8)    :: coeffx,coeffy,coeffz

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        !Defaults

        x = xx(ig)
        y = yy(jg)
        z = zz(kg)

        dhx = 2.*dxh(ig)
        dhy = 2.*dyh(jg)
        dhz = 2.*dzh(kg)

        vxx = rvx(i,j,k)/rho(i,j,k)
        vyy = rvy(i,j,k)/rho(i,j,k)
        vzz = rvz(i,j,k)/rho(i,j,k)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        coeffx = 1d0
        coeffy = 1d0
        coeffz = 1d0

        !Exceptions

        select case(half_elem)
        case (1)

cc          dhx = dx(ig)
cc          x   = (xx(ig) + xx(ig+1))/2.
cc          im  = i

cc          vxx = 0.5*(rvx(i ,j,k)/rho(i ,j,k)
cc     .              +rvx(ip,j,k)/rho(ip,j,k))
cc          vyy = 0.5*(rvy(i ,j,k)/rho(i ,j,k)
cc     .              +rvy(ip,j,k)/rho(ip,j,k))
cc          vzz = 0.5*(rvz(i ,j,k)/rho(i ,j,k)
cc     .              +rvz(ip,j,k)/rho(ip,j,k))

          jp = min(j+1,ny)
          jm = max(j-1,1)
          kp = min(k+1,nz)
          km = max(k-1,1)

          if (j == ny) dhy = dy(jg-1)
          if (j == 1 ) dhy = dy(jg)
          if (k == nz) dhz = dz(kg-1)
          if (k == 1 ) dhz = dz(kg)

          coeffx = 0d0

        case (2)
cc
cc          dhy = dy(jg)
cc          y   = (yy(jg) + yy(jg+1))/2.
cc          jm  = j

cc          vxx = 0.5*(rvx(i,j ,k)/rho(i,j ,k)
cc     .              +rvx(i,jp,k)/rho(i,jp,k))
cc          vyy = 0.5*(rvy(i,j ,k)/rho(i,j ,k)
cc     .              +rvy(i,jp,k)/rho(i,jp,k))
cc          vzz = 0.5*(rvz(i,j ,k)/rho(i,j ,k)
cc     .              +rvz(i,jp,k)/rho(i,jp,k))

          ip = min(i+1,nx)
          im = max(i-1,1)
          kp = min(k+1,nz)
          km = max(k-1,1)

          if (i == nx) dhx = dx(ig-1)
          if (i == 1 ) dhx = dx(ig)
          if (k == nz) dhz = dz(kg-1)
          if (k == 1 ) dhz = dz(kg)

          coeffy = 0d0

        case (3)
cc
cc          dhz = dz(kg)
cc          z   = (zz(kg) + zz(kg+1))/2.
cc          km  = k

cc          vxx = 0.5*(rvx(i,j,k )/rho(i,j,k )
cc     .              +rvx(i,j,kp)/rho(i,j,kp))
cc          vyy = 0.5*(rvy(i,j,k )/rho(i,j,k )
cc     .              +rvy(i,j,kp)/rho(i,j,kp))
cc          vzz = 0.5*(rvz(i,j,k )/rho(i,j,k )
cc     .              +rvz(i,j,kp)/rho(i,j,kp))

          ip = min(i+1,nx)
          im = max(i-1,1)
          jp = min(j+1,ny)
          jm = max(j-1,1)

          if (i == nx) dhx = dx(ig-1)
          if (i == 1 ) dhx = dx(ig)
          if (j == ny) dhy = dy(jg-1)
          if (j == 1 ) dhy = dy(jg)

          coeffz = 0d0

        end select

c     !Calculate nabla_v tensor

        hessian1 = hessian(1,x,y,z)
        hessian2 = hessian(2,x,y,z)
        hessian3 = hessian(3,x,y,z)

      ! l = 1, m = 1
        tensor(1,1) = coeffx*(rvx(ip,j,k)/rho(ip,j,k)
     .                       -rvx(im,j,k)/rho(im,j,k))/dhx
     .               + vxx*(hessian1(1,1)
     .                    + hessian2(2,1) + hessian3(3,1))
     .               - vxx*hessian1(1,1)
     .               - vyy*hessian1(2,1)
     .               - vzz*hessian1(3,1)

      ! l = 1, m = 2
        tensor(1,2) = coeffx*(rvy(ip,j,k)/rho(ip,j,k)
     .                       -rvy(im,j,k)/rho(im,j,k))/dhx
     .               + vyy*(hessian1(1,1)
     .                    + hessian2(2,1) + hessian3(3,1))
     .               - vxx*hessian2(1,1)
     .               - vyy*hessian2(2,1)
     .               - vzz*hessian2(3,1)

      ! l = 1, m = 3
        tensor(1,3) = coeffx*(rvz(ip,j,k)/rho(ip,j,k)
     .                       -rvz(im,j,k)/rho(im,j,k))/dhx
     .               + vzz*(hessian1(1,1)
     .                    + hessian2(2,1) + hessian3(3,1))
     .               - vxx*hessian3(1,1)
     .               - vyy*hessian3(2,1)
     .               - vzz*hessian3(3,1)

      ! l = 2, m = 1
        tensor(2,1) = coeffy*(rvx(i,jp,k)/rho(i,jp,k)
     .                       -rvx(i,jm,k)/rho(i,jm,k))/dhy
     .               + vxx*(hessian1(1,2)
     .                    + hessian2(2,2) + hessian3(3,2))
     .               - vxx*hessian1(1,2)
     .               - vyy*hessian1(2,2)
     .               - vzz*hessian1(3,2)

      ! l = 2, m = 2
        tensor(2,2) = coeffy*(rvy(i,jp,k)/rho(i,jp,k)
     .                       -rvy(i,jm,k)/rho(i,jm,k))/dhy
     .               + vyy*(hessian1(1,2)
     .                    + hessian2(2,2) + hessian3(3,2))
     .               - vxx*hessian2(1,2)
     .               - vyy*hessian2(2,2)
     .               - vzz*hessian2(3,2)

      ! l = 2, m = 3
        tensor(2,3) = coeffy*(rvz(i,jp,k)/rho(i,jp,k)
     .                       -rvz(i,jm,k)/rho(i,jm,k))/dhy
     .               + vzz*(hessian1(1,2)
     .                    + hessian2(2,2) + hessian3(3,2))
     .               - vxx*hessian3(1,2)
     .               - vyy*hessian3(2,2)
     .               - vzz*hessian3(3,2)

      ! l = 3, m = 1
        tensor(3,1) = coeffz*(rvx(i,j,kp)/rho(i,j,kp)
     .                       -rvx(i,j,km)/rho(i,j,km))/dhz
     .               + vxx*(hessian1(1,3)
     .                    + hessian2(2,3) + hessian3(3,3))
     .               - vxx*hessian1(1,3)
     .               - vyy*hessian1(2,3)
     .               - vzz*hessian1(3,3)

      ! l = 3, m = 2
        tensor(3,2) = coeffz*(rvy(i,j,kp)/rho(i,j,kp)
     .                       -rvy(i,j,km)/rho(i,j,km))/dhz
     .               + vyy*(hessian1(1,3)
     .                    + hessian2(2,3) + hessian3(3,3))
     .               - vxx*hessian2(1,3)
     .               - vyy*hessian2(2,3)
     .               - vzz*hessian2(3,3)

      ! l = 3, m = 3
        tensor(3,3) = coeffz*(rvz(i,j,kp)/rho(i,j,kp)
     .                       -rvz(i,j,km)/rho(i,j,km))/dhz
     .               + vzz*(hessian1(1,3)
     .                    + hessian2(2,3) + hessian3(3,3))
     .               - vxx*hessian3(1,3)
     .               - vyy*hessian3(2,3)
     .               - vzz*hessian3(3,3)

c     End program

      end function fnabla_v_bc

c     vtensor_x
c     #############################################################
      subroutine vtensor_x(i,j,k,t11,t12,t13)
c     -------------------------------------------------------------
c     Calculates fluxes in X-direction for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k
        real(8)    :: t11,t12,t13

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x,y,z
        real(8)    :: jac,ptot,vis

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x = (xx(ig) + xx(ig+1))/2.
        y =  yy(jg)
        z =  zz(kg)

        gsuper = g_super(x,y,z)

        jac = jacobian(x,y,z)

        nabla_v = fnabla_v(i,j,k,1)

        !Recall p=2nT
        ptot = jac*(rho(i+1,j,k)*tmp(i  ,j,k)
     .             +rho(i  ,j,k)*tmp(i+1,j,k))
     .       +(bx(i+1,j,k)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i+1,j,k)
     .        +by(i+1,j,k)*by_cov(i,j,k)+by(i,j,k)*by_cov(i+1,j,k)
     .        +bz(i+1,j,k)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i+1,j,k))/2.

        vis = 2./(1./rho(i+1,j,k)/nuu(i+1,j,k)
     .          + 1./rho(i  ,j,k)/nuu(i  ,j,k))

        t11 =
cc     .     0.25*( rho(i+1,j,k)*vx(i+1,j,k)*vx(i  ,j,k)
cc     .           +rho(i+1,j,k)*vx(i  ,j,k)*vx(i+1,j,k)
cc     .           +rho(i  ,j,k)*vx(i+1,j,k)*vx(i  ,j,k)
cc     .           +rho(i  ,j,k)*vx(i  ,j,k)*vx(i+1,j,k) )
     .     0.25*( rvx(i+1,j,k)*vx (i  ,j,k) +rvx(i  ,j,k)*vx (i+1,j,k)
     .           + vx(i  ,j,k)*rvx(i+1,j,k) +vx (i+1,j,k)*rvx(i  ,j,k))
     .     -0.5*( bx(i+1,j,k)*bx(i  ,j,k)
     .           +bx(i  ,j,k)*bx(i+1,j,k) )
     .     +gsuper(1,1)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,1)
     .           +gsuper(1,2)*nabla_v(2,1)
     .           +gsuper(1,3)*nabla_v(3,1) )

        t12 =
cc     .     0.25*( rho(i+1,j,k)*vx(i+1,j,k)*vy(i  ,j,k)
cc     .           +rho(i+1,j,k)*vx(i  ,j,k)*vy(i+1,j,k)
cc     .           +rho(i  ,j,k)*vx(i+1,j,k)*vy(i  ,j,k)
cc     .           +rho(i  ,j,k)*vx(i  ,j,k)*vy(i+1,j,k) )
     .     0.25*( rvx(i+1,j,k)*vy (i  ,j,k) +rvx(i  ,j,k)*vy (i+1,j,k)
     .           + vx(i  ,j,k)*rvy(i+1,j,k) +vx (i+1,j,k)*rvy(i  ,j,k))
     .     -0.5*( bx(i+1,j,k)*by(i  ,j,k)
     .           +bx(i  ,j,k)*by(i+1,j,k) )
     .     +gsuper(1,2)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,2)
     .           +gsuper(1,2)*nabla_v(2,2)
     .           +gsuper(1,3)*nabla_v(3,2) )

        t13 =
cc     .     0.25*( rho(i+1,j,k)*vx(i+1,j,k)*vz(i  ,j,k)
cc     .           +rho(i+1,j,k)*vx(i  ,j,k)*vz(i+1,j,k)
cc     .           +rho(i  ,j,k)*vx(i+1,j,k)*vz(i  ,j,k)
cc     .           +rho(i  ,j,k)*vx(i  ,j,k)*vz(i+1,j,k) )
     .     0.25*( rvx(i+1,j,k)*vz (i  ,j,k) +rvx(i  ,j,k)*vz (i+1,j,k)
     .           + vx(i  ,j,k)*rvz(i+1,j,k) +vx (i+1,j,k)*rvz(i  ,j,k))
     .     -0.5*( bx(i+1,j,k)*bz(i  ,j,k)
     .           +bx(i  ,j,k)*bz(i+1,j,k) )
     .     +gsuper(1,3)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,3)
     .           +gsuper(1,2)*nabla_v(2,3)
     .           +gsuper(1,3)*nabla_v(3,3) )

c     End program

      end subroutine vtensor_x

c     vtensor_y
c     #############################################################
      subroutine vtensor_y(i,j,k,t21,t22,t23)
c     -------------------------------------------------------------
c     Calculates fluxes in Y-direction for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k
        real(8)    :: t21,t22,t23

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x,y,z
        real(8)    :: jac,ptot,vis

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x =  xx(ig)
        y = (yy(jg) + yy(jg+1))/2.
        z =  zz(kg)

        gsuper = g_super(x,y,z)

        jac = jacobian(x,y,z)

        nabla_v = fnabla_v(i,j,k,2)

        !Recall p=2nT
        ptot = jac*(rho(i,j+1,k)*tmp(i,j  ,k)
     .             +rho(i,j  ,k)*tmp(i,j+1,k))
     .       +(bx(i,j+1,k)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,j+1,k)
     .        +by(i,j+1,k)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,j+1,k)
     .        +bz(i,j+1,k)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,j+1,k))/2.

        vis = 2./(1./rho(i,j+1,k)/nuu(i,j+1,k)
     .          + 1./rho(i,j  ,k)/nuu(i,j  ,k))

        t21 =
cc     .     0.25*( rho(i,j+1,k)*vy(i,j+1,k)*vx(i  ,j,k)
cc     .           +rho(i,j+1,k)*vy(i  ,j,k)*vx(i,j+1,k)
cc     .           +rho(i  ,j,k)*vy(i,j+1,k)*vx(i  ,j,k)
cc     .           +rho(i  ,j,k)*vy(i  ,j,k)*vx(i,j+1,k) )
     .     0.25*( rvy(i,j+1,k)*vx(i,j,k) + rvy(i,j,k)*vx(i,j+1,k)
     .           +rvx(i,j+1,k)*vy(i,j,k) + rvx(i,j,k)*vy(i,j+1,k))
     .     -0.5*( by(i,j+1,k)*bx(i,j  ,k)
     .           +by(i,j  ,k)*bx(i,j+1,k) )
     .     +gsuper(2,1)*ptot
     .     -vis*( gsuper(2,1)*nabla_v(1,1)
     .           +gsuper(2,2)*nabla_v(2,1)
     .           +gsuper(2,3)*nabla_v(3,1) )

        t22 =
cc     .     0.25*( rho(i,j+1,k)*vy(i,j+1,k)*vy(i  ,j,k)
cc     .           +rho(i,j+1,k)*vy(i  ,j,k)*vy(i,j+1,k)
cc     .           +rho(i  ,j,k)*vy(i,j+1,k)*vy(i  ,j,k)
cc     .           +rho(i  ,j,k)*vy(i  ,j,k)*vy(i,j+1,k) )
     .     0.25*( rvy(i,j+1,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,j+1,k)
     .           +rvy(i,j+1,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,j+1,k))
     .     -0.5*( by(i,j+1,k)*by(i,j  ,k)
     .           +by(i,j  ,k)*by(i,j+1,k) )
     .     +gsuper(2,2)*ptot
     .     -vis*( gsuper(2,1)*nabla_v(1,2)
     .           +gsuper(2,2)*nabla_v(2,2)
     .           +gsuper(2,3)*nabla_v(3,2) )

        t23 =
cc     .     0.25*( rho(i,j+1,k)*vy(i,j+1,k)*vz(i  ,j,k)
cc     .           +rho(i,j+1,k)*vy(i  ,j,k)*vz(i,j+1,k)
cc     .           +rho(i  ,j,k)*vy(i,j+1,k)*vz(i  ,j,k)
cc     .           +rho(i  ,j,k)*vy(i  ,j,k)*vz(i,j+1,k) )
     .     0.25*( rvy(i,j+1,k)*vz(i,j,k) + rvy(i,j,k)*vz(i,j+1,k)
     .           +rvz(i,j+1,k)*vy(i,j,k) + rvz(i,j,k)*vy(i,j+1,k))
     .     -0.5*( by(i,j+1,k)*bz(i,j  ,k)
     .           +by(i,j  ,k)*bz(i,j+1,k) )
     .     +gsuper(2,3)*ptot
     .     -vis*( gsuper(2,1)*nabla_v(1,3)
     .           +gsuper(2,2)*nabla_v(2,3)
     .           +gsuper(2,3)*nabla_v(3,3) )

c     End program

      end subroutine vtensor_y

c     vtensor_z
c     #############################################################
      subroutine vtensor_z(i,j,k,t31,t32,t33)
c     -------------------------------------------------------------
c     Calculates fluxes in Z-direction for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k
        real(8)    :: t31,t32,t33

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x,y,z
        real(8)    :: jac,ptot,vis

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x =  xx(ig)
        y =  yy(jg)
        z = (zz(kg) + zz(kg+1))/2.

        gsuper = g_super(x,y,z)

        jac = jacobian(x,y,z)

        nabla_v = fnabla_v(i,j,k,3)

        !Recall p=2nT
        ptot = jac*(rho(i,j,k+1)*tmp(i,j,k  )
     .             +rho(i,j,k  )*tmp(i,j,k+1))
     .       +(bx(i,j,k+1)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,j,k+1)
     .        +by(i,j,k+1)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,j,k+1)
     .        +bz(i,j,k+1)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,j,k+1))/2.

        vis = 2./(1./rho(i,j,k+1)/nuu(i,j,k+1)
     .          + 1./rho(i,j,k  )/nuu(i,j,k  ))

        t31 =
cc     .     0.25*( rho(i,j,k+1)*vz(i,j,k+1)*vx(i,j,k  )
cc     .           +rho(i,j,k+1)*vz(i,j,k  )*vx(i,j,k+1)
cc     .           +rho(i,j,k  )*vz(i,j,k+1)*vx(i,j,k  )
cc     .           +rho(i,j,k  )*vz(i,j,k  )*vx(i,j,k+1) )
     .     0.25*( rvz(i,j,k+1)*vx(i,j,k) + rvz(i,j,k)*vx(i,j,k+1)
     .           +rvx(i,j,k+1)*vz(i,j,k) + rvx(i,j,k)*vz(i,j,k+1) )
     .     -0.5*( bz(i,j,k+1)*bx(i,j,k  )
     .           +bz(i,j,k  )*bx(i,j,k+1) )
     .     +gsuper(3,1)*ptot
     .     -vis*( gsuper(3,1)*nabla_v(1,1)
     .           +gsuper(3,2)*nabla_v(2,1)
     .           +gsuper(3,3)*nabla_v(3,1) )

        t32 =
cc     .     0.25*( rho(i,j,k+1)*vz(i,j,k+1)*vy(i,j,k  )
cc     .           +rho(i,j,k+1)*vz(i,j,k  )*vy(i,j,k+1)
cc     .           +rho(i,j,k  )*vz(i,j,k+1)*vy(i,j,k  )
cc     .           +rho(i,j,k  )*vz(i,j,k  )*vy(i,j,k+1) )
     .     0.25*( rvz(i,j,k+1)*vy(i,j,k) + rvz(i,j,k)*vy(i,j,k+1)
     .           +rvy(i,j,k+1)*vz(i,j,k) + rvy(i,j,k)*vz(i,j,k+1) )
     .     -0.5*( bz(i,j,k+1)*by(i,j,k  )
     .           +bz(i,j,k  )*by(i,j,k+1) )
     .     +gsuper(3,2)*ptot
     .     -vis*( gsuper(3,1)*nabla_v(1,2)
     .           +gsuper(3,2)*nabla_v(2,2)
     .           +gsuper(3,3)*nabla_v(3,2) )

        t33 =
cc     .     0.25*( rho(i,j,k+1)*vz(i,j,k+1)*vz(i,j,k  )
cc     .           +rho(i,j,k+1)*vz(i,j,k  )*vz(i,j,k+1)
cc     .           +rho(i,j,k  )*vz(i,j,k+1)*vz(i,j,k  )
cc     .           +rho(i,j,k  )*vz(i,j,k  )*vz(i,j,k+1) )
     .     0.25*( rvz(i,j,k+1)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,k+1)
     .           +rvz(i,j,k+1)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,k+1) )
     .     -0.5*( bz(i,j,k+1)*bz(i  ,j,k)
     .           +bz(i  ,j,k)*bz(i,j,k+1) )
     .     +gsuper(3,3)*ptot
     .     -vis*( gsuper(3,1)*nabla_v(1,3)
     .           +gsuper(3,2)*nabla_v(2,3)
     .           +gsuper(3,3)*nabla_v(3,3) )

c     End program

      end subroutine vtensor_z

c     vflx_x
c     #############################################################
      function vflx_x(i,j,k,t11,t12,t13,cov) result(flx_x)
c     -------------------------------------------------------------
c     Calculates fluxes in X-direction for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k
        real(8)    :: flx_x,cov(3),t11,t12,t13

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x,y,z
        real(8)    :: cnv1(3),cnv2(3),cnv3(3)

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x = (xx(ig) + xx(ig+1))/2.
        y =  yy(jg)
        z =  zz(kg)

        cnv1   = contravariantVector(1,x,y,z)
        cnv2   = contravariantVector(2,x,y,z)
        cnv3   = contravariantVector(3,x,y,z)

        flx_x =  t11*dot_product(cnv1,cov)
     .          +t12*dot_product(cnv2,cov)
     .          +t13*dot_product(cnv3,cov)

c     End program

      end function vflx_x

c     vflx_y
c     #############################################################
      function vflx_y(i,j,k,t21,t22,t23,cov) result(flx_y)
c     -------------------------------------------------------------
c     Calculates fluxes in Y-direction for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k
        real(8)    :: flx_y,cov(3),t21,t22,t23

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x,y,z
        real(8)    :: cnv1(3),cnv2(3),cnv3(3)

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x =  xx(ig)
        y = (yy(jg) + yy(jg+1))/2.
        z =  zz(kg)

        cnv1   = contravariantVector(1,x,y,z)
        cnv2   = contravariantVector(2,x,y,z)
        cnv3   = contravariantVector(3,x,y,z)

        flx_y =  t21*dot_product(cnv1,cov)
     .          +t22*dot_product(cnv2,cov)
     .          +t23*dot_product(cnv3,cov)

c     End program

      end function vflx_y

c     vflx_z
c     #############################################################
      function vflx_z(i,j,k,t31,t32,t33,cov) result(flx_z)
c     -------------------------------------------------------------
c     Calculates fluxes in Z-direction for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k
        real(8)    :: flx_z,cov(3),t31,t32,t33

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: x,y,z
        real(8)    :: cnv1(3),cnv2(3),cnv3(3)

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        x =  xx(ig)
        y =  yy(jg)
        z = (zz(kg) + zz(kg+1))/2.

        cnv1   = contravariantVector(1,x,y,z)
        cnv2   = contravariantVector(2,x,y,z)
        cnv3   = contravariantVector(3,x,y,z)

        flx_z =  t31*dot_product(cnv1,cov)
     .          +t32*dot_product(cnv2,cov)
     .          +t33*dot_product(cnv3,cov)

c     End program

      end function vflx_z

c     laplacian
c     ###############################################################
      real*8 function laplacian(i,j,k,arr)

c     ---------------------------------------------------------------
c     Calculates lap(arr) in general non-orthogonal coordinates.
c     ---------------------------------------------------------------

      use grid

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1
      
      if (i == 0 .or. j == 0 .or. k == 0 ) then
        write (*,*) 'Error in laplace; i,j,k=0'
      elseif (i == nx+1 .or. j == ny+1 .or. k == nz+1) then
        write (*,*) 'Error in laplace; i,j,k=nmax+1'
      elseif (.not.sing_point) then
        d_xx_ip = g_super_elem(1,1,0.5*(xx(ig+1)+xx(ig)),yy(jg),zz(kg))
        d_xx_im = g_super_elem(1,1,0.5*(xx(ig-1)+xx(ig)),yy(jg),zz(kg))
        d_yy_jp = g_super_elem(2,2,xx(ig),0.5*(yy(jg+1)+yy(jg)),zz(kg))
        d_yy_jm = g_super_elem(2,2,xx(ig),0.5*(yy(jg-1)+yy(jg)),zz(kg))
        d_zz_kp = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg+1)+zz(kg)))
        d_zz_km = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg-1)+zz(kg)))

        d_xy_ipjp = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
     .                              ,0.5*(yy(jg+1)+yy(jg)),zz(kg))
        d_xy_ipjm = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
     .                              ,0.5*(yy(jg-1)+yy(jg)),zz(kg))
        d_xy_imjp = g_super_elem(1,2,0.5*(xx(ig-1)+xx(ig))
     .                              ,0.5*(yy(jg+1)+yy(jg)),zz(kg))
        d_xy_imjm = g_super_elem(1,2,0.5*(xx(ig-1)+xx(ig))
     .                              ,0.5*(yy(jg-1)+yy(jg)),zz(kg))
                 
        d_xz_ipkp = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg+1)+zz(kg)))
        d_xz_ipkm = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg-1)+zz(kg)))
        d_xz_imkp = g_super_elem(1,3,0.5*(xx(ig-1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg+1)+zz(kg)))
        d_xz_imkm = g_super_elem(1,3,0.5*(xx(ig-1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg-1)+zz(kg)))
                 
        d_yz_jpkp = g_super_elem(2,3,xx(ig),0.5*(yy(jg+1)+yy(jg))
     .                                     ,0.5*(zz(kg+1)+zz(kg)))
        d_yz_jpkm = g_super_elem(2,3,xx(ig),0.5*(yy(jg+1)+yy(jg))
     .                                     ,0.5*(zz(kg-1)+zz(kg)))
        d_yz_jmkp = g_super_elem(2,3,xx(ig),0.5*(yy(jg-1)+yy(jg))
     .                                     ,0.5*(zz(kg+1)+zz(kg)))
        d_yz_jmkm = g_super_elem(2,3,xx(ig),0.5*(yy(jg-1)+yy(jg))
     .                                     ,0.5*(zz(kg-1)+zz(kg)))

      else

        d_xx_ip = g_super_elem(1,1,0.5*(xx(ig+1)+xx(ig)),yy(jg),zz(kg))
        d_xx_im = 0d0
        d_yy_jp = 0d0
        d_yy_jm = 0d0
        d_zz_kp = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg+1)+zz(kg)))
        d_zz_km = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg-1)+zz(kg)))

        d_xy_ipjp = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
     .                              ,0.5*(yy(jg+1)+yy(jg)),zz(kg))
        d_xy_ipjm = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
     .                              ,0.5*(yy(jg-1)+yy(jg)),zz(kg))
        d_xy_imjp = 0d0
        d_xy_imjm = 0d0
                 
        d_xz_ipkp = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg+1)+zz(kg)))
        d_xz_ipkm = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg-1)+zz(kg)))
        d_xz_imkp = 0d0
        d_xz_imkm = 0d0
                 
        d_yz_jpkp = g_super_elem(2,3,xx(ig),0.5*(yy(jg+1)+yy(jg))
     .                                     ,0.5*(zz(kg+1)+zz(kg)))
        d_yz_jpkm = g_super_elem(2,3,xx(ig),0.5*(yy(jg+1)+yy(jg))
     .                                     ,0.5*(zz(kg-1)+zz(kg)))
        d_yz_jmkp = g_super_elem(2,3,xx(ig),0.5*(yy(jg-1)+yy(jg))
     .                                     ,0.5*(zz(kg+1)+zz(kg)))
        d_yz_jmkm = g_super_elem(2,3,xx(ig),0.5*(yy(jg-1)+yy(jg))
     .                                     ,0.5*(zz(kg-1)+zz(kg)))

      endif

      laplacian =
     .     dyh(jg)*dzh(kg)*( d_xx_ip*(arr(ip,j,k)-arr(i,j,k))/dx(ig)
     .                      +d_xx_im*(arr(im,j,k)-arr(i,j,k))/dx(ig-1) )
     .    +dxh(ig)*dzh(kg)*( d_yy_jp*(arr(i,jp,k)-arr(i,j,k))/dy(jg)
     .                      +d_yy_jm*(arr(i,jm,k)-arr(i,j,k))/dy(jg-1) )
     .    +dxh(ig)*dyh(jg)*( d_zz_kp*(arr(i,j,kp)-arr(i,j,k))/dz(kg)
     .                      +d_zz_km*(arr(i,j,km)-arr(i,j,k))/dz(kg-1) )
     .    +0.5*
     .      ( dzh(kg)*( d_xy_ipjp*(arr(ip,jp,k)-arr(i,j,k))
     .                 +d_xy_imjm*(arr(im,jm,k)-arr(i,j,k))
     .                 -d_xy_ipjm*(arr(ip,jm,k)-arr(i,j,k))
     .                 -d_xy_imjp*(arr(im,jp,k)-arr(i,j,k)) )
     .       +dyh(jg)*( d_xz_ipkp*(arr(ip,j,kp)-arr(i,j,k))
     .                 +d_xz_imkm*(arr(im,j,km)-arr(i,j,k))
     .                 -d_xz_ipkm*(arr(ip,j,km)-arr(i,j,k))
     .                 -d_xz_imkp*(arr(im,j,kp)-arr(i,j,k)) )
     .       +dxh(jg)*( d_yz_jpkp*(arr(i,jp,kp)-arr(i,j,k))
     .                 +d_yz_jmkm*(arr(i,jm,km)-arr(i,j,k))
     .                 -d_yz_jpkm*(arr(i,jp,km)-arr(i,j,k))
     .                 -d_yz_jmkp*(arr(i,jm,kp)-arr(i,j,k)) ) )

c     End program

      end function laplacian

      end module nlfunction_setup

c module precond_variables
c ######################################################################
      module precond_variables

        use timeStepping

        use precond_setup

        use transport_params

cc        use vectorToArrayXfer

        use equilibrium

        use parameters

        use constants

cc        real(8), allocatable, dimension(:)
cc     .                       :: bxx,byy,vxx,vyy,vxe,vye,diag_mu
cc
cc        real(8), allocatable, dimension(:,:,:) :: pp
cc
cc        integer, allocatable, dimension(:,:) :: idiagp,bcs
cc
cc        real(8),pointer,    dimension(:,:) :: uu,vv,tt,ww,nn
cc        real(8),allocatable,dimension(:,:) :: ue
cc
cc        real(8),allocatable,dimension(:,:):: d_sc,d_bc,d_pc,d_bh,d_alf
cc
cc        logical :: formdiag
cc
cc        integer :: iiout
cc
cc        integer :: bc_sc(4,1),bc_cpld(4,2)

      end module precond_variables
