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
     .          bx_cov,by_cov,bz_cov,bx_car,by_car,bz_car
     .         ,jx,jy,jz,jx_cov,jy_cov,jz_cov,jx_car,jy_car,jz_car
     .         ,vx,vy,vz,vx_cov,vy_cov,vz_cov,vx_car,vy_car,vz_car
     .         ,eeta,nuu,divrgB

        real(8),pointer,dimension(:,:,:):: rho,rvx,rvy,rvz,bx,by,bz,tmp
        real(8),pointer,dimension(:)    :: xx,yy,zz,dxh,dyh,dzh,dx,dy,dz

        integer(4) :: igx,igy,igz,nx,ny,nz

        real(8)    :: x0,y0,z0,xim,yim,zim,xip,yip,zip
     .                        ,xjm,yjm,zjm,xjp,yjp,zjp
     .                        ,xkm,ykm,zkm,xkp,ykp,zkp

        real(8)    :: gsub(3,3),gsuper(3,3)
     .               ,cnv1(3),cnv2(3),cnv3(3),cov(3),jac

        real(8)    :: nabla_v(3,3),hessian1(3,3)
     .               ,hessian2(3,3),hessian3(3,3)
     .               ,cov_tnsr(3,3),cnv_tnsr(3,3)

        logical    :: sing_point,cartesian

      contains

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

      if (i == 1 .and. bcond(1) == SP) then
        divB =
     .       (bx(i+1,j  ,k  )+bx(i  ,j  ,k  ))/2./dxh(ig)
     .      +(by(i  ,j+1,k  )-by(i  ,j-1,k  ))/2./dyh(jg)
     .      +(bz(i  ,j  ,k+1)-bz(i  ,j  ,k-1))/2./dzh(kg)
      else
        divB =
     .       (bx(i+1,j  ,k  )-bx(i-1,j  ,k  ))/2./dxh(ig)
     .      +(by(i  ,j+1,k  )-by(i  ,j-1,k  ))/2./dyh(jg)
     .      +(bz(i  ,j  ,k+1)-bz(i  ,j  ,k-1))/2./dzh(kg)
      endif

c     End 

      end function divB

c     divV
c     ###############################################################
      real(8) function divV(i,j,k)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of velocity field at grid vertices.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: i,j,k

c     Local variables

      integer(4) :: ig,jg,kg

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      if (i == 1 .and. bcond(1) == SP) then
        divV =
     .       (rvx(i+1,j  ,k  )/rho(i+1,j  ,k  )
     .       +rvx(i  ,j  ,k  )/rho(i  ,j  ,k  ))/2./dxh(ig)
     .      +(rvy(i  ,j+1,k  )/rho(i  ,j+1,k  )
     .       -rvy(i  ,j-1,k  )/rho(i  ,j-1,k  ))/2./dyh(jg)
     .      +(rvz(i  ,j  ,k+1)/rho(i  ,j  ,k+1)
     .       -rvz(i  ,j  ,k-1)/rho(i  ,j  ,k-1))/2./dzh(kg)
      else
        divV =
     .       (rvx(i+1,j  ,k  )/rho(i+1,j  ,k  )
     .       -rvx(i-1,j  ,k  )/rho(i-1,j  ,k  ))/2./dxh(ig)
     .      +(rvy(i  ,j+1,k  )/rho(i  ,j+1,k  )
     .       -rvy(i  ,j-1,k  )/rho(i  ,j-1,k  ))/2./dyh(jg)
     .      +(rvz(i  ,j  ,k+1)/rho(i  ,j  ,k+1)
     .       -rvz(i  ,j  ,k-1)/rho(i  ,j  ,k-1))/2./dzh(kg)
      endif

c     End 

      end function divV

ccc     jouleHeating
ccc     #############################################################
cc      real*8 function jouleHeating(i,j,k) result(joule)
ccc     -------------------------------------------------------------
ccc     Calculates joule heating at cell volume (i,j,k) based on
ccc     edge currents. Result includes cell volume.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: joule1,joule2,joule3
cc
ccc     Begin program
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        joule1 = (dx(ig  )*dy(jg  )*dzh(kg)*joule_ipjp(i  ,j  ,k)
cc     .           +dx(ig-1)*dy(jg  )*dzh(kg)*joule_ipjp(i-1,j  ,k)
cc     .           +dx(ig  )*dy(jg-1)*dzh(kg)*joule_ipjp(i  ,j-1,k)
cc     .           +dx(ig-1)*dy(jg-1)*dzh(kg)*joule_ipjp(i-1,j-1,k))/4.
cc
cc        joule2 = (dx(ig  )*dyh(jg)*dz(kg  )*joule_ipkp(i  ,j,k  )
cc     .           +dx(ig-1)*dyh(jg)*dz(kg  )*joule_ipkp(i-1,j,k  )
cc     .           +dx(ig  )*dyh(jg)*dz(kg-1)*joule_ipkp(i  ,j,k-1)
cc     .           +dx(ig-1)*dyh(jg)*dz(kg-1)*joule_ipkp(i-1,j,k-1))/4.
cc
cc        joule3 = (dxh(ig)*dy(jg  )*dz(kg  )*joule_jpkp(i,j  ,k  )
cc     .           +dxh(ig)*dy(jg-1)*dz(kg  )*joule_jpkp(i,j-1,k  )
cc     .           +dxh(ig)*dy(jg  )*dz(kg-1)*joule_jpkp(i,j  ,k-1)
cc     .           +dxh(ig)*dy(jg-1)*dz(kg-1)*joule_jpkp(i,j-1,k-1))/4.
cc
cc        joule = (joule1 + joule2 + joule3)/3.
cc
ccc     End program
cc
cc      contains
cc
ccc     joule_ipjp
ccc     #############################################################
cc      real(8) function joule_ipjp(i,j,k)
cc
ccc     -------------------------------------------------------------
ccc     Calculates joule heating at cell volume (i+1/2,j+1/2,k) based
ccc     on edge currents. Result includes cell volume.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: jj(3),x,y,z,gsub(3,3),eta_ipjp
cc
ccc     Begin program
cc
cc        eta_ipjp = 4./(1./eeta(i+1,j+1,k) + 1./eeta(i,j+1,k)
cc     .                +1./eeta(i+1,j  ,k) + 1./eeta(i,j  ,k))
cc
cccc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        x = 0.5*(xx(ig)+xx(ig+1))
cc        y = 0.5*(yy(jg)+yy(jg+1))
cc        z = zz(kg)
cc
cc        gsub = g_sub(x,y,z)
cc
cc        !Contravarian current components at (ip,jp,k)
cc
cc        jj(1) = (bz_cov(i  ,j+1,k) - bz_cov(i  ,j,k)
cc     .          +bz_cov(i+1,j+1,k) - bz_cov(i+1,j,k) )/2./dy(jg)
cc     .         -(by_cov(i,j+1,k+1) + by_cov(i  ,j  ,k+1)
cc     .          +by_cov(i+1,j,k+1) + by_cov(i+1,j+1,k+1)
cc     .          -by_cov(i,j+1,k-1) - by_cov(i  ,j  ,k-1)
cc     .          -by_cov(i+1,j,k-1) - by_cov(i+1,j+1,k-1))/8./dzh(kg)
cc
cc        jj(2) = (bx_cov(i,j+1,k+1) + bx_cov(i  ,j  ,k+1)
cc     .          +bx_cov(i+1,j,k+1) + bx_cov(i+1,j+1,k+1)
cc     .          -bx_cov(i,j+1,k-1) - bx_cov(i  ,j  ,k-1)
cc     .          -bx_cov(i+1,j,k-1) - bx_cov(i+1,j+1,k-1))/8./dzh(kg)
cc     .         -(bz_cov(i+1,j  ,k) - bz_cov(i,j  ,k)
cc     .          +bz_cov(i+1,j+1,k) - bz_cov(i,j+1,k) )/2./dx(ig)
cc
cc        jj(3) = (by_cov(i+1,j  ,k) - by_cov(i,j  ,k)
cc     .          +by_cov(i+1,j+1,k) - by_cov(i,j+1,k) )/2./dx(ig)
cc     .         -(bx_cov(i  ,j+1,k) - bx_cov(i  ,j,k)
cc     .          +bx_cov(i+1,j+1,k) - bx_cov(i+1,j,k) )/2./dy(jg)
cc
cc        joule_ipjp = eta_ipjp*dot_product(jj,matmul(gsub,jj))
cc
cc      end function joule_ipjp
cc
ccc     joule_ipkp
ccc     #############################################################
cc      real(8) function joule_ipkp(i,j,k)
cc
ccc     -------------------------------------------------------------
ccc     Calculates joule heating at cell volume (i+1/2,j,k+1/2) based
ccc     on edge currents.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: jj(3),x,y,z,gsub(3,3),eta_ipkp
cc
ccc     Begin program
cc
cc        eta_ipkp = 4./(1./eeta(i+1,j,k+1) + 1./eeta(i,j,k+1)
cc     .                +1./eeta(i+1,j,k  ) + 1./eeta(i,j,k  ))
cc
cccc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        x = 0.5*(xx(ig)+xx(ig+1))
cc        y =      yy(jg)
cc        z = 0.5*(zz(kg)+zz(kg+1))
cc
cc        gsub = g_sub(x,y,z)
cc
cc        jj(1) = (bz_cov(i  ,j+1,k+1) + bz_cov(i  ,j+1,k  )
cc     .          +bz_cov(i+1,j+1,k  ) + bz_cov(i+1,j+1,k+1)
cc     .          -bz_cov(i  ,j-1,k+1) - bz_cov(i  ,j-1,k  )
cc     .          -bz_cov(i+1,j-1,k  ) - bz_cov(i+1,j-1,k+1))/8./dyh(jg)
cc     .         -(by_cov(i  ,j,k+1) - by_cov(i  ,j,k)
cc     .          +by_cov(i+1,j,k+1) - by_cov(i+1,j,k) )/2./dz(kg)
cc
cc        jj(2) = (bx_cov(i  ,j,k+1) - bx_cov(i  ,j,k)
cc     .          +bx_cov(i+1,j,k+1) - bx_cov(i+1,j,k) )/2./dz(kg)
cc     .         -(bz_cov(i+1,j,k  ) - bz_cov(i,j,k  )
cc     .          +bz_cov(i+1,j,k+1) - bz_cov(i,j,k+1) )/2./dx(ig)
cc
cc        jj(3) = (by_cov(i+1,j,k  ) - by_cov(i  ,j,k  )
cc     .          +by_cov(i+1,j,k+1) - by_cov(i  ,j,k+1) )/2./dx(ig)
cc     .         -(bx_cov(i  ,j+1,k+1) + bx_cov(i  ,j+1,k  )
cc     .          +bx_cov(i+1,j+1,k  ) + bx_cov(i+1,j+1,k+1)
cc     .          -bx_cov(i  ,j-1,k+1) - bx_cov(i  ,j-1,k  )
cc     .          -bx_cov(i+1,j-1,k  ) - bx_cov(i+1,j-1,k+1))/8./dyh(jg)
cc
cc        !Contravarian current components at (ip,j,kp)
cc
cc        joule_ipkp = eta_ipkp*dot_product(jj,matmul(gsub,jj))
cc
cc      end function joule_ipkp
cc
ccc     joule_jpkp
ccc     #############################################################
cc      real(8) function joule_jpkp(i,j,k)
cc
ccc     -------------------------------------------------------------
ccc     Calculates joule heating at cell volume (i,j+1/2,k+1/2) based
ccc     on edge currents.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: jj(3),x,y,z,gsub(3,3),eta_jpkp
cc
ccc     Begin program
cc
cc        eta_jpkp = 4./(1./eeta(i,j+1,k+1) + 1./eeta(i,j,k+1)
cc     .                +1./eeta(i,j+1,k  ) + 1./eeta(i,j,k  ))
cc
cccc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        x = 0.5*(xx(ig)+xx(ig+1))
cc        y =      yy(jg)
cc        z = 0.5*(zz(kg)+zz(kg+1))
cc
cc        gsub = g_sub(x,y,z)
cc
cc        jj(1) = (bz_cov(i,j+1,k  ) - bz_cov(i,j,k  )
cc     .          +bz_cov(i,j+1,k+1) - bz_cov(i,j,k+1) )/2./dy(jg)
cc     .         -(by_cov(i,j  ,k+1) - by_cov(i,j  ,k)
cc     .          +by_cov(i,j+1,k+1) - by_cov(i,j+1,k) )/2./dz(kg)
cc
cc        jj(2) = (bx_cov(i,j  ,k+1) - bx_cov(i,j  ,k)
cc     .          +bx_cov(i,j+1,k+1) - bx_cov(i,j+1,k) )/2./dz(kg)
cc     .         -(bz_cov(i+1,j  ,k+1) + bz_cov(i+1,j+1,k  )
cc     .          +bz_cov(i+1,j+1,k+1) + bz_cov(i+1,j  ,k  )
cc     .          -bz_cov(i-1,j  ,k+1) + bz_cov(i-1,j+1,k  )
cc     .          +bz_cov(i-1,j+1,k+1) + bz_cov(i-1,j  ,k  ))/8./dxh(ig)
cc
cc        jj(3) = (by_cov(i+1,j  ,k+1) + by_cov(i+1,j+1,k  )
cc     .          +by_cov(i+1,j+1,k+1) + by_cov(i+1,j  ,k  )
cc     .          -by_cov(i-1,j  ,k+1) + by_cov(i-1,j+1,k  )
cc     .          +by_cov(i-1,j+1,k+1) + by_cov(i-1,j  ,k  ))/8./dxh(ig)
cc     .         -(bx_cov(i,j+1,k  ) - bx_cov(i,j,k  )
cc     .          +bx_cov(i,j+1,k+1) - bx_cov(i,j,k+1) )/2./dy(jg)
cc
cc        !Contravarian current components at (ip,j,kp)
cc
cc        joule_jpkp = eta_jpkp*dot_product(jj,matmul(gsub,jj))
cc
cc      end function joule_jpkp
cc
cc      end function jouleHeating

c     fnabla_v
c     #############################################################
      function fnabla_v(i,j,k,x,y,z,half_elem) result(tensor)
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
        real(8)    :: tensor(3,3),x,y,z

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: jx,jy,jz
        real(8)    :: dhx,dhy,dhz,vxx,vyy,vzz
cc        logical    :: cartesian

c     Begin program

        !Defaults

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

cc        x = xx(ig)
cc        y = yy(jg)
cc        z = zz(kg)
cc
cc        cartesian = .false.

cc        call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x,y,z,cartesian)

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
cc          x   = (xx(ig) + xx(ig+1))/2.
          im  = i
          vxx = 0.5*(vx(i,j,k)+vx(ip,j,k))
          vyy = 0.5*(vy(i,j,k)+vy(ip,j,k))
          vzz = 0.5*(vz(i,j,k)+vz(ip,j,k))
        case (2)
          dhy = dy(jg)
cc          y   = (yy(jg) + yy(jg+1))/2.
          jm  = j
          vxx = 0.5*(vx(i,j,k)+vx(i,jp,k))
          vyy = 0.5*(vy(i,j,k)+vy(i,jp,k))
          vzz = 0.5*(vz(i,j,k)+vz(i,jp,k))
        case (3)
          dhz = dz(kg)
cc          z   = (zz(kg) + zz(kg+1))/2.
          km  = k
          vxx = 0.5*(vx(i,j,k)+vx(i,j,kp))
          vyy = 0.5*(vy(i,j,k)+vy(i,j,kp))
          vzz = 0.5*(vz(i,j,k)+vz(i,j,kp))
        end select

c     !Calculate nabla_v tensor

        hessian1 = hessian(1,x,y,z,cartesian)
        hessian2 = hessian(2,x,y,z,cartesian)
        hessian3 = hessian(3,x,y,z,cartesian)

      ! l = 1, m = 1
        tensor(1,1) = (vx(ip,j,k)-vx(im,j,k))/dhx
     .               + vxx*(hessian1(1,1)
     .                    + hessian2(2,1)
     .                    + hessian3(3,1))
     .               - vxx*hessian1(1,1)
     .               - vyy*hessian1(2,1)
     .               - vzz*hessian1(3,1)

      ! l = 1, m = 2
        tensor(1,2) = (vy(ip,j,k)-vy(im,j,k))/dhx
     .               + vyy*(hessian1(1,1)
     .                    + hessian2(2,1)
     .                    + hessian3(3,1))
     .               - vxx*hessian2(1,1)
     .               - vyy*hessian2(2,1)
     .               - vzz*hessian2(3,1)

      ! l = 1, m = 3
        tensor(1,3) = (vz(ip,j,k)-vz(im,j,k))/dhx
     .               + vzz*(hessian1(1,1)
     .                    + hessian2(2,1)
     .                    + hessian3(3,1))
     .               - vxx*hessian3(1,1)
     .               - vyy*hessian3(2,1)
     .               - vzz*hessian3(3,1)

      ! l = 2, m = 1
        tensor(2,1) = (vx(i,jp,k)-vx(i,jm,k))/dhy
     .               + vxx*(hessian1(1,2)
     .                    + hessian2(2,2)
     .                    + hessian3(3,2))
     .               - vxx*hessian1(1,2)
     .               - vyy*hessian1(2,2)
     .               - vzz*hessian1(3,2)

      ! l = 2, m = 2
        tensor(2,2) = (vy(i,jp,k)-vy(i,jm,k))/dhy
     .               + vyy*(hessian1(1,2)
     .                    + hessian2(2,2)
     .                    + hessian3(3,2))
     .               - vxx*hessian2(1,2)
     .               - vyy*hessian2(2,2)
     .               - vzz*hessian2(3,2)

      ! l = 2, m = 3
        tensor(2,3) = (vz(i,jp,k)-vz(i,jm,k))/dhy
     .               + vzz*(hessian1(1,2)
     .                    + hessian2(2,2)
     .                    + hessian3(3,2))
     .               - vxx*hessian3(1,2)
     .               - vyy*hessian3(2,2)
     .               - vzz*hessian3(3,2)

      ! l = 3, m = 1
        tensor(3,1) = (vx(i,j,kp)-vx(i,j,km))/dhz
     .               + vxx*(hessian1(1,3)
     .                    + hessian2(2,3)
     .                    + hessian3(3,3))
     .               - vxx*hessian1(1,3)
     .               - vyy*hessian1(2,3)
     .               - vzz*hessian1(3,3)

      ! l = 3, m = 2
        tensor(3,2) = (vy(i,j,kp)-vy(i,j,km))/dhz
     .               + vyy*(hessian1(1,3)
     .                    + hessian2(2,3)
     .                    + hessian3(3,3))
     .               - vxx*hessian2(1,3)
     .               - vyy*hessian2(2,3)
     .               - vzz*hessian2(3,3)

      ! l = 3, m = 3
        tensor(3,3) = (vz(i,j,kp)-vz(i,j,km))/dhz
     .               + vzz*(hessian1(1,3)
     .                    + hessian2(2,3)
     .                    + hessian3(3,3))
     .               - vxx*hessian3(1,3)
     .               - vyy*hessian3(2,3)
     .               - vzz*hessian3(3,3)

c     End program

      end function fnabla_v

c     fnabla_v_bc
c     #############################################################
      function fnabla_v_bc(i,j,k,x,y,z,cartesian,half_elem)
     .         result(tensor)
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
        real(8)    :: tensor(3,3),x,y,z
        logical    :: cartesian

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: dhx,dhy,dhz
        real(8)    :: vxx,vyy,vzz,vxip,vyip,vzip,vxim,vyim,vzim
     .                           ,vxjp,vyjp,vzjp,vxjm,vyjm,vzjm
     .                           ,vxkp,vykp,vzkp,vxkm,vykm,vzkm
        real(8)    :: coeffx,coeffy,coeffz

c     Begin program

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc        cartesian = .false.

        !Defaults

cc        x = xx(ig)
cc        y = yy(jg)
cc        z = zz(kg)

        dhx = 2.*dxh(ig)
        dhy = 2.*dyh(jg)
        dhz = 2.*dzh(kg)

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

cc          dhy = dy(jg)
cc          y   = (yy(jg) + yy(jg+1))/2.
cc          jm  = j

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

cc          dhz = dz(kg)
cc          z   = (zz(kg) + zz(kg+1))/2.
cc          km  = k

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

c    Calculate grid quantities

cc        gsuper   = G_super(x,y,z,cartesian)
        hessian1 = hessian(1,x,y,z,cartesian)
        hessian2 = hessian(2,x,y,z,cartesian)
        hessian3 = hessian(3,x,y,z,cartesian)

c    Calculate velocities

        vxx = rvx(i,j,k)/rho(i,j,k)
        vyy = rvy(i,j,k)/rho(i,j,k)
        vzz = rvz(i,j,k)/rho(i,j,k)

        vxip = rvx(ip,j,k)/rho(ip,j,k)
        vyip = rvy(ip,j,k)/rho(ip,j,k)
        vzip = rvz(ip,j,k)/rho(ip,j,k)

        vxim = rvx(im,j,k)/rho(im,j,k)
        vyim = rvy(im,j,k)/rho(im,j,k)
        vzim = rvz(im,j,k)/rho(im,j,k)

        vxjp = rvx(i,jp,k)/rho(i,jp,k)
        vyjp = rvy(i,jp,k)/rho(i,jp,k)
        vzjp = rvz(i,jp,k)/rho(i,jp,k)
                                    
        vxjm = rvx(i,jm,k)/rho(i,jm,k)
        vyjm = rvy(i,jm,k)/rho(i,jm,k)
        vzjm = rvz(i,jm,k)/rho(i,jm,k)

        vxkp = rvx(i,j,kp)/rho(i,j,kp)
        vykp = rvy(i,j,kp)/rho(i,j,kp)
        vzkp = rvz(i,j,kp)/rho(i,j,kp)
                                     
        vxkm = rvx(i,j,km)/rho(i,j,km)
        vykm = rvy(i,j,km)/rho(i,j,km)
        vzkm = rvz(i,j,km)/rho(i,j,km)

c     Calculate nabla_v tensor

      ! l = 1, m = 1
        tensor(1,1) = coeffx*(vxip-vxim)/dhx
     .               + vxx*(hessian1(1,1)
     .                    + hessian2(2,1)
     .                    + hessian3(3,1))
     .               - vxx*hessian1(1,1)
     .               - vyy*hessian1(2,1)
     .               - vzz*hessian1(3,1)

      ! l = 1, m = 2
        tensor(1,2) = coeffx*(vyip-vyim)/dhx
     .               + vyy*(hessian1(1,1)
     .                    + hessian2(2,1)
     .                    + hessian3(3,1))
     .               - vxx*hessian2(1,1)
     .               - vyy*hessian2(2,1)
     .               - vzz*hessian2(3,1)

      ! l = 1, m = 3
        tensor(1,3) = coeffx*(vzip-vzim)/dhx
     .               + vzz*(hessian1(1,1)
     .                    + hessian2(2,1)
     .                    + hessian3(3,1))
     .               - vxx*hessian3(1,1)
     .               - vyy*hessian3(2,1)
     .               - vzz*hessian3(3,1)

      ! l = 2, m = 1
        tensor(2,1) = coeffy*(vxjp-vxjm)/dhy
     .               + vxx*(hessian1(1,2)
     .                    + hessian2(2,2)
     .                    + hessian3(3,2))
     .               - vxx*hessian1(1,2)
     .               - vyy*hessian1(2,2)
     .               - vzz*hessian1(3,2)

      ! l = 2, m = 2
        tensor(2,2) = coeffy*(vyjp-vyjm)/dhy
     .               + vyy*(hessian1(1,2)
     .                    + hessian2(2,2)
     .                    + hessian3(3,2))
     .               - vxx*hessian2(1,2)
     .               - vyy*hessian2(2,2)
     .               - vzz*hessian2(3,2)

      ! l = 2, m = 3
        tensor(2,3) = coeffy*(vzjp-vzjm)/dhy
     .               + vzz*(hessian1(1,2)
     .                    + hessian2(2,2)
     .                    + hessian3(3,2))
     .               - vxx*hessian3(1,2)
     .               - vyy*hessian3(2,2)
     .               - vzz*hessian3(3,2)

      ! l = 3, m = 1
        tensor(3,1) = coeffz*(vxkp-vxkm)/dhz
     .               + vxx*(hessian1(1,3)
     .                    + hessian2(2,3)
     .                    + hessian3(3,3))
     .               - vxx*hessian1(1,3)
     .               - vyy*hessian1(2,3)
     .               - vzz*hessian1(3,3)

      ! l = 3, m = 2
        tensor(3,2) = coeffz*(vykp-vykm)/dhz
     .               + vyy*(hessian1(1,3)
     .                    + hessian2(2,3)
     .                    + hessian3(3,3))
     .               - vxx*hessian2(1,3)
     .               - vyy*hessian2(2,3)
     .               - vzz*hessian2(3,3)

      ! l = 3, m = 3
        tensor(3,3) = coeffz*(vzkp-vzkm)/dhz
     .               + vzz*(hessian1(1,3)
     .                    + hessian2(2,3)
     .                    + hessian3(3,3))
     .               - vxx*hessian3(1,3)
     .               - vyy*hessian3(2,3)
     .               - vzz*hessian3(3,3)

c     End program

      end function fnabla_v_bc

c     vtensor_x
c     #############################################################
      subroutine vtensor_x(i,j,k,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag
        real(8)    :: t11,t12,t13

c     Local variables

        integer(4) :: ig,jg,kg,ip
        real(8)    :: x,y,z
        real(8)    :: jac,ptot,vis
cc        logical    :: cartesian

c     Begin program

        ip = i+1

cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        x = (xx(ig) + xx(ig+1))/2.
cc        y =  yy(jg)
cc        z =  zz(kg)
cc
cc        if (flag == 0) then
cc          ip = i
cc          x = xx(ig)
cc        endif
cc
cc        cartesian = .false.

        if (flag == 0) then
          ip = i
          x = x0
          y = y0
          z = z0
        elseif (flag == 1) then
          x = (xip+x0)/2.
          y = (yip+y0)/2.
          z = (zip+z0)/2.
        else
          x = (xim+x0)/2.
          y = (yim+y0)/2.
          z = (zim+z0)/2.
        endif

        gsuper = g_super (x,y,z,cartesian)

        jac    = jacobian(x,y,z,cartesian)

        if (flag /= 0 .or. sing_point) then
          nabla_v = fnabla_v(i,j,k,x,y,z,1)
        else
          nabla_v = fnabla_v(i,j,k,x,y,z,0)
        endif

        !Recall p=2nT
        ptot = jac*(rho(ip,j,k)*tmp(ip,j,k)
     .             +rho(i ,j,k)*tmp(i ,j,k))
cc        ptot = jac*(rho(ip,j,k)*tmp(i ,j,k)
cc     .             +rho(i ,j,k)*tmp(ip,j,k))
     .       +(bx(ip,j,k)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(ip,j,k)
     .        +by(ip,j,k)*by_cov(i,j,k)+by(i,j,k)*by_cov(ip,j,k)
     .        +bz(ip,j,k)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(ip,j,k))/4.

        vis = 2./(1./rho(ip,j,k)/nuu(ip,j,k)
     .          + 1./rho(i ,j,k)/nuu(i ,j,k))

        t11 =
     .     0.25*( rvx(ip,j,k)*vx (i ,j,k) +rvx(i ,j,k)*vx (ip,j,k)
     .           + vx(i ,j,k)*rvx(ip,j,k) +vx (ip,j,k)*rvx(i ,j,k))
     .     -0.5*( bx(ip,j,k)*bx(i ,j,k)
     .           +bx(i ,j,k)*bx(ip,j,k) )
     .     +gsuper(1,1)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,1)
     .           +gsuper(1,2)*nabla_v(2,1)
     .           +gsuper(1,3)*nabla_v(3,1) )

        t12 =
     .     0.25*( rvx(ip,j,k)*vy (i ,j,k) +rvx(i ,j,k)*vy (ip,j,k)
     .           + vx(i ,j,k)*rvy(ip,j,k) +vx (ip,j,k)*rvy(i ,j,k))
     .     -0.5*( bx(ip,j,k)*by(i ,j,k)
     .           +bx(i ,j,k)*by(ip,j,k) )
     .     +gsuper(1,2)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,2)
     .           +gsuper(1,2)*nabla_v(2,2)
     .           +gsuper(1,3)*nabla_v(3,2) )

cc        t12 = 0d0

        t13 =
     .     0.25*( rvx(ip,j,k)*vz (i ,j,k) +rvx(i ,j,k)*vz (ip,j,k)
     .           + vx(i ,j,k)*rvz(ip,j,k) +vx (ip,j,k)*rvz(i ,j,k))
     .     -0.5*( bx(ip,j,k)*bz(i ,j,k)
     .           +bx(i ,j,k)*bz(ip,j,k) )
     .     +gsuper(1,3)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,3)
     .           +gsuper(1,2)*nabla_v(2,3)
     .           +gsuper(1,3)*nabla_v(3,3) )

        t11 = t11/jac
        t12 = t12/jac
        t13 = t13/jac

c     End program

      end subroutine vtensor_x

c     vtensor_y
c     #############################################################
      subroutine vtensor_y(i,j,k,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag
        real(8)    :: t21,t22,t23

c     Local variables

        integer(4) :: ig,jg,kg,jp
        real(8)    :: x,y,z
        real(8)    :: jac,ptot,vis
cc        logical    :: cartesian

c     Begin program

        jp = j+1

cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        x =  xx(ig)
cc        y = (yy(jg) + yy(jg+1))/2.
cc        z =  zz(kg)
cc
cc        cartesian = .false.
cc
cc        if (flag == 0) then
cc          jp = j
cc          y = yy(jg)
cc        endif

        if (flag == 0) then
          jp = j
          x = x0
          y = y0
          z = z0
        elseif (flag == 1) then
          x = (xjp+x0)/2.
          y = (yjp+y0)/2.
          z = (zjp+z0)/2.
        else
          x = (xjm+x0)/2.
          y = (yjm+y0)/2.
          z = (zjm+z0)/2.
        endif

        gsuper = g_super (x,y,z,cartesian)

        jac    = jacobian(x,y,z,cartesian)

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,x,y,z,2)
        else
          nabla_v = fnabla_v(i,j,k,x,y,z,0)
        endif

        !Recall p=2nT
        ptot = jac*(rho(i,jp,k)*tmp(i,jp,k)
     .             +rho(i,j ,k)*tmp(i,j ,k))
     .       +(bx(i,jp,k)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,jp,k)
     .        +by(i,jp,k)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,jp,k)
     .        +bz(i,jp,k)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,jp,k))/4.

        vis = 2./(1./rho(i,jp,k)/nuu(i,jp,k)
     .          + 1./rho(i,j ,k)/nuu(i,j ,k))

        t21 =
     .     0.25*( rvy(i,jp,k)*vx(i,j,k) + rvy(i,j,k)*vx(i,jp,k)
     .           +rvx(i,jp,k)*vy(i,j,k) + rvx(i,j,k)*vy(i,jp,k))
     .     -0.5*( by(i,jp,k)*bx(i,j ,k)
     .           +by(i,j ,k)*bx(i,jp,k) )
     .     +gsuper(2,1)*ptot
     .     -vis*( gsuper(2,1)*nabla_v(1,1)
     .           +gsuper(2,2)*nabla_v(2,1)
     .           +gsuper(2,3)*nabla_v(3,1) )

        t22 =
     .     0.25*( rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k)
     .           +rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k))
     .     -0.5*( by(i,jp,k)*by(i,j ,k)
     .           +by(i,j ,k)*by(i,jp,k) )
     .     +gsuper(2,2)*ptot
     .     -vis*( gsuper(2,1)*nabla_v(1,2)
     .           +gsuper(2,2)*nabla_v(2,2)
     .           +gsuper(2,3)*nabla_v(3,2) )

        t23 =
     .     0.25*( rvy(i,jp,k)*vz(i,j,k) + rvy(i,j,k)*vz(i,jp,k)
     .           +rvz(i,jp,k)*vy(i,j,k) + rvz(i,j,k)*vy(i,jp,k))
     .     -0.5*( by(i,jp,k)*bz(i,j ,k)
     .           +by(i,j ,k)*bz(i,jp,k) )
     .     +gsuper(2,3)*ptot
     .     -vis*( gsuper(2,1)*nabla_v(1,3)
     .           +gsuper(2,2)*nabla_v(2,3)
     .           +gsuper(2,3)*nabla_v(3,3) )

        t21 = t21/jac
        t22 = t22/jac
        t23 = t23/jac

c     End program

      end subroutine vtensor_y

c     vtensor_z
c     #############################################################
      subroutine vtensor_z(i,j,k,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag
        real(8)    :: t31,t32,t33

c     Local variables

        integer(4) :: ig,jg,kg,kp
        real(8)    :: x,y,z
        real(8)    :: jac,ptot,vis
cc        logical    :: cartesian

c     Begin program

        kp=k+1

cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        cartesian = .false.
cc
cc        x =  xx(ig)
cc        y =  yy(jg)
cc        z = (zz(kg) + zz(kg+1))/2.
cc
cc        if (flag == 0) then
cc          kp = k
cc          z = zz(kg)
cc        endif

        if (flag == 0) then
          kp = k
          x = x0
          y = y0
          z = z0
        elseif (flag == 1) then
          x = (xkp+x0)/2.
          y = (ykp+y0)/2.
          z = (zkp+z0)/2.
        else
          x = (xkm+x0)/2.
          y = (ykm+y0)/2.
          z = (zkm+z0)/2.
        endif

        gsuper = g_super (x,y,z,cartesian)

        jac    = jacobian(x,y,z,cartesian)

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,x,y,z,3)
        else
          nabla_v = fnabla_v(i,j,k,x,y,z,0)
        endif

        !Recall p=2nT
        ptot = jac*(rho(i,j,kp)*tmp(i,j,kp)
     .             +rho(i,j,k )*tmp(i,j,k ))
     .       +(bx(i,j,kp)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,j,kp)
     .        +by(i,j,kp)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,j,kp)
     .        +bz(i,j,kp)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,j,kp))/4.

        vis = 2./(1./rho(i,j,kp)/nuu(i,j,kp)
     .          + 1./rho(i,j,k )/nuu(i,j,k ))

        t31 =
     .     0.25*( rvz(i,j,kp)*vx(i,j,k) + rvz(i,j,k)*vx(i,j,kp)
     .           +rvx(i,j,kp)*vz(i,j,k) + rvx(i,j,k)*vz(i,j,kp) )
     .     -0.5*( bz(i,j,kp)*bx(i,j,k )
     .           +bz(i,j,k )*bx(i,j,kp) )
     .     +gsuper(3,1)*ptot
     .     -vis*( gsuper(3,1)*nabla_v(1,1)
     .           +gsuper(3,2)*nabla_v(2,1)
     .           +gsuper(3,3)*nabla_v(3,1) )

        t32 =
     .     0.25*( rvz(i,j,kp)*vy(i,j,k) + rvz(i,j,k)*vy(i,j,kp)
     .           +rvy(i,j,kp)*vz(i,j,k) + rvy(i,j,k)*vz(i,j,kp) )
     .     -0.5*( bz(i,j,kp)*by(i,j,k )
     .           +bz(i,j,k )*by(i,j,kp) )
     .     +gsuper(3,2)*ptot
     .     -vis*( gsuper(3,1)*nabla_v(1,2)
     .           +gsuper(3,2)*nabla_v(2,2)
     .           +gsuper(3,3)*nabla_v(3,2) )

        t33 =
     .     0.25*( rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp)
     .           +rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp) )
     .     -0.5*( bz(i,j,kp)*bz(i  ,j,k)
     .           +bz(i  ,j,k)*bz(i,j,kp) )
     .     +gsuper(3,3)*ptot
     .     -vis*( gsuper(3,1)*nabla_v(1,3)
     .           +gsuper(3,2)*nabla_v(2,3)
     .           +gsuper(3,3)*nabla_v(3,3) )

        t31 = t31/jac
        t32 = t32/jac
        t33 = t33/jac

c     End program

      end subroutine vtensor_z

ccc     vflx_x2
ccc     #############################################################
cc      function vflx_x2(i,j,k,t11,t12,t13,cov,flag) result(flx_x)
ccc     -------------------------------------------------------------
ccc     Calculates fluxes in X-direction for EOM
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,flag
cc        real(8)    :: flx_x,cov(3),t11,t12,t13
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: x,y,z,jac,jacp,ptot
cc        real(8)    :: cnv1(3),cnv2(3),cnv3(3),cov1(3)
cc        logical    :: cartesian
cc
ccc     Begin program
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        cartesian = .false.
cc
cc        x = (xx(ig) + xx(ig+1))/2.
cc        y =  yy(jg)
cc        z =  zz(kg)
cc
cc        jac = 1d0
cc        if (.not.conserv) then
cc          if (flag == 1) then
cc            x = xx(ig)
cc          else
cc            x = xx(ig+1)
cc          endif
cc          jac = jacobian(x,y,z,cartesian)
cc        endif
cc
cc        cnv1   = contravariantVector(1,x,y,z,cartesian)
cc        cnv2   = contravariantVector(2,x,y,z,cartesian)
cc        cnv3   = contravariantVector(3,x,y,z,cartesian)
cc
cc        flx_x =  t11*dot_product(cnv1,cov)*jac
cc     .          +t12*dot_product(cnv2,cov)*jac
cc     .          +t13*dot_product(cnv3,cov)*jac
cc
ccc     End program
cc
cc      end function vflx_x2

c     vflx_x
c     #############################################################
      function vflx_x(t11,t12,t13,cov) result(flx_x)
c     -------------------------------------------------------------
c     Calculates fluxes in X-direction for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag
        real(8)    :: x,y,z,flx_x,cov(3),t11,t12,t13

c     Local variables

c     Begin program

        flx_x =  t11*dot_product(cnv1,cov)*jac
     .          +t12*dot_product(cnv2,cov)*jac
     .          +t13*dot_product(cnv3,cov)*jac

c     End program

      end function vflx_x

ccc     vflx_y2
ccc     #############################################################
cc      function vflx_y2(i,j,k,t21,t22,t23,cov,flag) result(flx_y)
ccc     -------------------------------------------------------------
ccc     Calculates fluxes in Y-direction for EOM
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,flag
cc        real(8)    :: flx_y,cov(3),t21,t22,t23
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: x,y,z,jac,jacp,ptot
cc        real(8)    :: cnv1(3),cnv2(3),cnv3(3),cov1(3)
cc        logical    :: cartesian
cc
ccc     Begin program
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        cartesian = .false.
cc
cc        x =  xx(ig)
cc        y = (yy(jg) + yy(jg+1))/2.
cc        z =  zz(kg)
cc
cc        jac = 1d0
cc        if (.not.conserv) then
cc          if (flag == 1) then
cc            y = yy(jg)
cc          else
cc            y = yy(jg+1)
cc          endif
cc          jac = jacobian(x,y,z,cartesian)
cc        endif
cc
cc        cnv1   = contravariantVector(1,x,y,z,cartesian)
cc        cnv2   = contravariantVector(2,x,y,z,cartesian)
cc        cnv3   = contravariantVector(3,x,y,z,cartesian)
cc
cc        flx_y =  t21*dot_product(cnv1,cov)*jac
cc     .          +t22*dot_product(cnv2,cov)*jac
cc     .          +t23*dot_product(cnv3,cov)*jac
cc
ccc     End program
cc
cc      end function vflx_y2

c     vflx_y
c     #############################################################
      function vflx_y(t21,t22,t23,cov) result(flx_y)
c     -------------------------------------------------------------
c     Calculates fluxes in Y-direction for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        real(8)    :: flx_y,cov(3),t21,t22,t23

c     Local variables

c     Begin program

        flx_y =  t21*dot_product(cnv1,cov)*jac
     .          +t22*dot_product(cnv2,cov)*jac
     .          +t23*dot_product(cnv3,cov)*jac

c     End program

      end function vflx_y

ccc     vflx_z2
ccc     #############################################################
cc      function vflx_z2(i,j,k,t31,t32,t33,cov,flag) result(flx_z)
ccc     -------------------------------------------------------------
ccc     Calculates fluxes in Z-direction for EOM
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,flag
cc        real(8)    :: flx_z,cov(3),t31,t32,t33
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg
cc        real(8)    :: x,y,z,ptot,jac,jacp
cc        real(8)    :: cnv1(3),cnv2(3),cnv3(3),cov1(3)
cc        logical    :: cartesian
cc
ccc     Begin program
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        cartesian = .false.
cc
cc        x =  xx(ig)
cc        y =  yy(jg)
cc        z = (zz(kg) + zz(kg+1))/2.
cc
cc        jac = 1d0
cc        if (.not.conserv) then
cc          if (flag == 1) then
cc            z = zz(kg)
cc          else
cc            z = zz(kg+1)
cc          endif
cc          jac = jacobian(x,y,z,cartesian)
cc        endif
cc
cc        cnv1 = contravariantVector(1,x,y,z,cartesian)
cc        cnv2 = contravariantVector(2,x,y,z,cartesian)
cc        cnv3 = contravariantVector(3,x,y,z,cartesian)
cc
cc        flx_z =  t31*dot_product(cnv1,cov)*jac
cc     .          +t32*dot_product(cnv2,cov)*jac
cc     .          +t33*dot_product(cnv3,cov)*jac
cc
ccc     End program
cc
cc      end function vflx_z2

c     vflx_z
c     #############################################################
      function vflx_z(t31,t32,t33,cov) result(flx_z)
c     -------------------------------------------------------------
c     Calculates fluxes in Z-direction for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        real(8)    :: flx_z,cov(3),t31,t32,t33

c     Local variables

c     Begin program

        flx_z =  t31*dot_product(cnv1,cov)*jac
     .          +t32*dot_product(cnv2,cov)*jac
     .          +t33*dot_product(cnv3,cov)*jac

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
        d_xx_ip = g_super_elem(1,1,0.5*(xx(ig+1)+xx(ig)),yy(jg),zz(kg)
     .                        ,.false.)
        d_xx_im = g_super_elem(1,1,0.5*(xx(ig-1)+xx(ig)),yy(jg),zz(kg)
     .                        ,.false.)
        d_yy_jp = g_super_elem(2,2,xx(ig),0.5*(yy(jg+1)+yy(jg)),zz(kg)
     .                        ,.false.)
        d_yy_jm = g_super_elem(2,2,xx(ig),0.5*(yy(jg-1)+yy(jg)),zz(kg)
     .                        ,.false.)
        d_zz_kp = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg+1)+zz(kg))
     .                        ,.false.)
        d_zz_km = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg-1)+zz(kg))
     .                        ,.false.)

        d_xy_ipjp = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
     .                              ,0.5*(yy(jg+1)+yy(jg))
     .                              ,zz(kg),.false.)
        d_xy_ipjm = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
     .                              ,0.5*(yy(jg-1)+yy(jg))
     .                              ,zz(kg),.false.)
        d_xy_imjp = g_super_elem(1,2,0.5*(xx(ig-1)+xx(ig))
     .                              ,0.5*(yy(jg+1)+yy(jg))
     .                              ,zz(kg),.false.)
        d_xy_imjm = g_super_elem(1,2,0.5*(xx(ig-1)+xx(ig))
     .                              ,0.5*(yy(jg-1)+yy(jg))
     .                              ,zz(kg),.false.)
                 
        d_xz_ipkp = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg+1)+zz(kg)),.false.)
        d_xz_ipkm = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg-1)+zz(kg)),.false.)
        d_xz_imkp = g_super_elem(1,3,0.5*(xx(ig-1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg+1)+zz(kg)),.false.)
        d_xz_imkm = g_super_elem(1,3,0.5*(xx(ig-1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg-1)+zz(kg)),.false.)
                 
        d_yz_jpkp = g_super_elem(2,3,xx(ig)
     .                              ,0.5*(yy(jg+1)+yy(jg))
     .                              ,0.5*(zz(kg+1)+zz(kg)),.false.)
        d_yz_jpkm = g_super_elem(2,3,xx(ig)
     .                              ,0.5*(yy(jg+1)+yy(jg))
     .                              ,0.5*(zz(kg-1)+zz(kg)),.false.)
        d_yz_jmkp = g_super_elem(2,3,xx(ig)
     .                              ,0.5*(yy(jg-1)+yy(jg))
     .                              ,0.5*(zz(kg+1)+zz(kg)),.false.)
        d_yz_jmkm = g_super_elem(2,3,xx(ig)
     .                              ,0.5*(yy(jg-1)+yy(jg))
     .                              ,0.5*(zz(kg-1)+zz(kg)),.false.)

      else

        d_xx_ip = g_super_elem(1,1,0.5*(xx(ig+1)+xx(ig)),yy(jg),zz(kg)
     .                         ,.false.)
        d_xx_im = 0d0
        d_yy_jp = 0d0
        d_yy_jm = 0d0
        d_zz_kp = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg+1)+zz(kg))
     .                        ,.false.)
        d_zz_km = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg-1)+zz(kg))
     .                        ,.false.)

        d_xy_ipjp = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
     .                              ,0.5*(yy(jg+1)+yy(jg))
     .                              ,zz(kg),.false.)
        d_xy_ipjm = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
     .                              ,0.5*(yy(jg-1)+yy(jg))
     .                              ,zz(kg),.false.)
        d_xy_imjp = 0d0
        d_xy_imjm = 0d0
                 
        d_xz_ipkp = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg+1)+zz(kg)),.false.)
        d_xz_ipkm = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
     .                              ,0.5*(zz(kg-1)+zz(kg)),.false.)
        d_xz_imkp = 0d0
        d_xz_imkm = 0d0
                 
        d_yz_jpkp = g_super_elem(2,3,xx(ig)
     .                              ,0.5*(yy(jg+1)+yy(jg))
     .                              ,0.5*(zz(kg+1)+zz(kg)),.false.)
        d_yz_jpkm = g_super_elem(2,3,xx(ig)
     .                              ,0.5*(yy(jg+1)+yy(jg))
     .                              ,0.5*(zz(kg-1)+zz(kg)),.false.)
        d_yz_jmkp = g_super_elem(2,3,xx(ig)
     .                              ,0.5*(yy(jg-1)+yy(jg))
     .                              ,0.5*(zz(kg+1)+zz(kg)),.false.)
        d_yz_jmkm = g_super_elem(2,3,xx(ig)
     .                              ,0.5*(yy(jg-1)+yy(jg))
     .                              ,0.5*(zz(kg-1)+zz(kg)),.false.)

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
     .       +dxh(ig)*( d_yz_jpkp*(arr(i,jp,kp)-arr(i,j,k))
     .                 +d_yz_jmkm*(arr(i,jm,km)-arr(i,j,k))
     .                 -d_yz_jpkm*(arr(i,jp,km)-arr(i,j,k))
     .                 -d_yz_jmkp*(arr(i,jm,kp)-arr(i,j,k)) ) )

cc      laplacian =
cc     .     dyh(jg)*dzh(kg)*( (arr(ip,j,k)-arr(i,j,k))/dx(ig)
cc     .                      +(arr(im,j,k)-arr(i,j,k))/dx(ig-1) )
cc     .    +dxh(ig)*dzh(kg)*( (arr(i,jp,k)-arr(i,j,k))/dy(jg)
cc     .                      +(arr(i,jm,k)-arr(i,j,k))/dy(jg-1) )
cc     .    +dxh(ig)*dyh(jg)*( (arr(i,j,kp)-arr(i,j,k))/dz(kg)
cc     .                      +(arr(i,j,km)-arr(i,j,k))/dz(kg-1) )

c     End program

      end function laplacian

c     curl
c     ###############################################################
      real*8 function curl(i,j,k,ax,ay,az,comp)

c     ---------------------------------------------------------------
c     Calculates curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      use grid

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,comp

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg
      real(8)    :: dx,dy,dz,flxip,flxim,flxjp,flxjm,flxkp,flxkm

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dx = grid_params%dxh(ig)
      dy = grid_params%dyh(jg)
      dz = grid_params%dzh(kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      select case(comp)
      case(1)

        flxip = 0d0
        flxim = 0d0

        flxjp = 0.5*(az(i,jp,k)+az(i,j,k))
        flxjm = 0.5*(az(i,jm,k)+az(i,j,k))

        flxkp =-0.5*(ay(i,j,kp)+ay(i,j,k))
        flxkm =-0.5*(ay(i,j,km)+ay(i,j,k))

      case(2)

        flxjp = 0d0
        flxjm = 0d0

        flxip =-0.5*(az(ip,j,k)+az(i,j,k))
        flxim =-0.5*(az(im,j,k)+az(i,j,k))

        flxkp = 0.5*(ax(i,j,kp)+ax(i,j,k))
        flxkm = 0.5*(ax(i,j,km)+ax(i,j,k))

      case(3)

        flxkp = 0d0
        flxkm = 0d0

        flxip = 0.5*(ay(ip,j,k)+ay(i,j,k))
        flxim = 0.5*(ay(im,j,k)+ay(i,j,k))

        flxjp =-0.5*(ax(i,jp,k)+ax(i,j,k))
        flxjm =-0.5*(ax(i,jm,k)+ax(i,j,k))

      case default

        write (*,*) 'Error in component in curl'
        write (*,*) 'Aborting...'
        stop

      end select

      call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .                      ,bcond)

      curl = (flxip-flxim)/dx
     .      +(flxjp-flxjm)/dy
     .      +(flxkp-flxkm)/dz

c     End program

      end function curl

c     curl2
c     ###############################################################
      real*8 function curl2(i,j,k,ax,ay,az,comp)

c     ---------------------------------------------------------------
c     Calculates curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      use grid

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,comp

      real(8)       :: ax(0:nx+1,0:ny+1,0:nz+1)
     .                ,ay(0:nx+1,0:ny+1,0:nz+1)
     .                ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg
      real(8)    :: dx,dy,dz,dx1,dx2

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dx = grid_params%dxh(ig)
      dy = grid_params%dyh(jg)
      dz = grid_params%dzh(kg)

      if (i == 1 .and. bcond(1) == SP) then
        dx1=grid_params%dx(ig)
        dx2=grid_params%dx(ig+1)
      endif

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      select case(comp)
      case(1)

        curl2 = (az(i,jp,k)-az(i,jm,k))/2./dy
     .         -(ay(i,j,kp)-ay(i,j,km))/2./dz

      case(2)

        if (i == 1 .and. bcond(1) == SP) then
          curl2 = (ax(i,j,kp)-ax(i,j,km))/2./dz
     .          -(-az(i+2,j,k)*dx1/dx2/(dx1+dx2)
     .            +az(i+1,j,k)*(dx1+dx2)/dx1/dx2
     .            -az(i  ,j,k)*(1./dx1+1./(dx1+dx2)))
        else
          curl2 = (ax(i,j,kp)-ax(i,j,km))/2./dz
     .           -(az(ip,j,k)-az(im,j,k))/2./dx
        endif

      case(3)

        if (i == 1 .and. bcond(1) == SP) then
          curl2 = (-ay(i+2,j,k)*dx1/dx2/(dx1+dx2)
     .             +ay(i+1,j,k)*(dx1+dx2)/dx1/dx2
     .             -ay(i  ,j,k)*(1./dx1+1/(dx1+dx2)))
     .            -(ax(i,jp,k)-ax(i,jm,k))/2./dy
        else
          curl2 = (ay(ip,j,k)-ay(im,j,k))/2./dx
     .           -(ax(i,jp,k)-ax(i,jm,k))/2./dy
        endif

      case default

        write (*,*) 'Error in component in curl'
        write (*,*) 'Aborting...'
        stop

      end select

c     End program

      end function curl2

c     curlcurl
c     ###############################################################
      real*8 function curlcurl(i,j,k,ax,ay,az,comp)

c     ---------------------------------------------------------------
c     Calculates curl(eta curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      use grid

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,comp

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: x0,xp,xm,y0,yp,ym,z0,zp,zm
      real(8)    :: dwa,dwb,dwc,dwd,dwe,dwf,dwg,dwh,dwi,dwj,dwk

      !Component one
      real(8)    :: d_yy_kp,d_yy_km,d_zz_jp,d_zz_jm
     .             ,d_xy_kp,d_xy_km,d_yz_kp,d_yz_km
     .             ,d_xz_jp,d_xz_jm,d_yz_jp,d_yz_jm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm

      !Component two
      real(8)    :: d_zz_ip,d_zz_im,d_xx_kp,d_xx_km
     .             ,d_xz_kp,d_xz_km
     .             ,d_xz_ip,d_xz_im,d_yz_ip,d_yz_im
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm

      !Component three
      real(8)    :: d_yy_ip,d_yy_im,d_xx_jp,d_xx_jm
     .             ,d_xy_jp,d_xy_jm,d_xy_ip,d_xy_im
cc     .             ,d_xz_jp,d_xz_jm,d_yz_ip,d_yz_im
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      x0 = xx(ig)
      xp = 0.5*(xx(ig+1)+xx(ig))
      xm = 0.5*(xx(ig-1)+xx(ig))

      y0 = yy(jg)
      yp = 0.5*(yy(jg+1)+yy(jg))
      ym = 0.5*(yy(jg-1)+yy(jg))

      z0 = zz(kg)
      zp = 0.5*(zz(kg+1)+zz(kg))
      zm = 0.5*(zz(kg-1)+zz(kg))

      select case(comp)
      case(1)

        d_yy_kp = 2./(1./eeta(i,j,kp)+1./eeta(i,j,k))
     .           *g_sub_elem(2,2,x0,y0,zp,.false.)
        d_yy_km = 2./(1./eeta(i,j,km)+1./eeta(i,j,k))
     .           *g_sub_elem(2,2,x0,y0,zm,.false.)

        d_zz_jp = 2./(1./eeta(i,jp,k)+1./eeta(i,j,k))
     .           *g_sub_elem(3,3,x0,yp,z0,.false.)
        d_zz_jm = 2./(1./eeta(i,jm,k)+1./eeta(i,j,k))
     .           *g_sub_elem(3,3,x0,ym,z0,.false.)

        d_yz_jpkp = 4./(1./eeta(i,jp,kp)+1./eeta(i,jp,k )
     .                 +1./eeta(i,j ,kp)+1./eeta(i,j ,k ))
     .             *g_sub_elem(2,3,x0,yp,zp,.false.)
        d_yz_jpkm = 4./(1./eeta(i,jp,km)+1./eeta(i,jp,k )
     .                 +1./eeta(i,j ,km)+1./eeta(i,j ,k ))
     .             *g_sub_elem(2,3,x0,yp,zm,.false.)
        d_yz_jmkp = 4./(1./eeta(i,jm,kp)+1./eeta(i,jm,k )
     .                 +1./eeta(i,j ,kp)+1./eeta(i,j ,k ))
     .             *g_sub_elem(2,3,x0,ym,zp,.false.)
        d_yz_jmkm = 4./(1./eeta(i,jm,km)+1./eeta(i,jm,k )
     .                 +1./eeta(i,j ,km)+1./eeta(i,j ,k ))
     .             *g_sub_elem(2,3,x0,ym,zm,.false.)

        d_xy_kp = 2./(1./eeta(i,j,kp)+1./eeta(i,j,k))
     .           *g_sub_elem(1,2,x0,y0,zp,.false.)
        d_xy_km = 2./(1./eeta(i,j,km)+1./eeta(i,j,k))
     .           *g_sub_elem(1,2,x0,y0,zm,.false.)

        d_yz_kp = 2./(1./eeta(i,j,kp)+1./eeta(i,j,k))
     .           *g_sub_elem(2,3,x0,y0,zp,.false.)
        d_yz_km = 2./(1./eeta(i,j,km)+1./eeta(i,j,k))
     .           *g_sub_elem(2,3,x0,y0,zm,.false.)

        d_xz_jp = 2./(1./eeta(i,jp,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,3,x0,yp,z0,.false.)
        d_xz_jm = 2./(1./eeta(i,jm,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,3,x0,ym,z0,.false.)

        d_yz_jp = 2./(1./eeta(i,jp,k)+1./eeta(i,j,k))
     .           *g_sub_elem(2,3,x0,yp,z0,.false.)
        d_yz_jm = 2./(1./eeta(i,jm,k)+1./eeta(i,j,k))
     .           *g_sub_elem(2,3,x0,ym,z0,.false.)

        dwa = -dxh(ig)*dyh(jg)
     .              *(d_yy_kp*(ax(i,j,kp)-ax(i,j,k ))/dz(kg  )
     .               -d_yy_km*(ax(i,j,k )-ax(i,j,km))/dz(kg-1) )
        dwb = -dxh(ig)*dzh(kg)
     .              *(d_zz_jp*(ax(i,jp,k)-ax(i,j ,k))/dy(jg  )
     .               -d_zz_jm*(ax(i,j ,k)-ax(i,jm,k))/dy(jg-1) )
        dwc = -dxh(ig)/2.
     .              *(-d_yz_jpkp*(ax(i,jp,kp)-ax(i,j ,k ))
     .                +d_yz_jmkm*(ax(i,j ,k )-ax(i,jm,km))
     .                +d_yz_jpkm*(ax(i,jp,km)-ax(i,j ,k ))
     .                -d_yz_jmkp*(ax(i,j ,k )-ax(i,jm,kp)))
        dwd =  dyh(jg)/4.
     .              *(d_yy_kp*(az(ip,j,kp)-az(im,j,kp)
     .                        +az(ip,j,k )-az(im,j,k ))
     .               -d_yy_km*(az(ip,j,km)-az(im,j,km)
     .                        +az(ip,j,k )-az(im,j,k )))
        dwe = -dxh(ig)/4.
     .              *(d_xy_kp*(az(i,jp,kp)-az(i,jm,kp)
     .                        +az(i,jp,k )-az(i,jm,k ))
     .               -d_xy_km*(az(i,jp,km)-az(i,jm,km)
     .                        +az(i,jp,k )-az(i,jm,k )))
        dwf =  dxh(ig)*dyh(jg)
     .              *(d_xy_kp*(ay(i,j,kp)-ay(i,j,k ))/dz(kg)
     .               -d_xy_km*(ay(i,j,k )-ay(i,j,km))/dz(kg-1))
        dwg = -dyh(jg)/4.
     .              *(d_yz_kp*(ay(ip,j,kp)-ay(im,j,kp)
     .                        +ay(ip,j,k )-ay(im,j,k ))
     .               -d_yz_km*(ay(ip,j,km)-ay(im,j,km)
     .                        +ay(ip,j,k )-ay(im,j,k )))
        dwh =  dzh(kg)/4.
     .              *(d_zz_jp*(ay(ip,jp,k)-ay(im,jp,k)
     .                        +ay(ip,j ,k)-ay(im,j ,k))
     .               -d_zz_jm*(ay(ip,jm,k)-ay(im,jm,k)
     .                        +ay(ip,j ,k)-ay(im,j ,k)))
        dwi =  dxh(ig)*dzh(kg)
     .              *(d_xz_jp*(az(i,jp,k)-az(i,j ,k))/dy(jg  )
     .               -d_xz_jm*(az(i,j ,k)-az(i,jm,k))/dy(jg-1))
        dwj = -dxh(ig)/4.
     .              *(d_xz_jp*(ay(i,jp,kp)-ay(i,jp,km)
     .                        +ay(i,j ,kp)-ay(i,j ,km))
     .               -d_xz_jm*(ay(i,jm,kp)-ay(i,jm,km)
     .                        +ay(i,j ,kp)-ay(i,j ,km)))
        dwk = -dzh(kg)/4.
     .              *(d_yz_jp*(az(ip,jp,k)-az(im,jp,k)
     .                        +az(ip,j ,k)-az(im,j ,k))
     .               -d_yz_jm*(az(ip,jm,k)-az(im,jm,k)
     .                        +az(ip,j ,k)-az(im,j ,k)))

      case(2)

        d_xx_kp = 2./(1./eeta(i,j,kp)+1./eeta(i,j,k))
     .           *g_sub_elem(1,1,x0,y0,zp,.false.)
        d_xx_km = 2./(1./eeta(i,j,km)+1./eeta(i,j,k))
     .           *g_sub_elem(1,1,x0,y0,zm,.false.)

        d_zz_ip = 2./(1./eeta(ip,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(3,3,xp,y0,z0,.false.)
        d_zz_im = 2./(1./eeta(im,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(3,3,xm,y0,z0,.false.)

        d_xz_ipkp = 4./(1./eeta(ip,j,kp)+1./eeta(ip,j,k )
     .                 +1./eeta(i ,j,kp)+1./eeta(i ,j,k ))
     .             *g_sub_elem(1,3,xp,y0,zp,.false.)
        d_xz_ipkm = 4./(1./eeta(ip,j,km)+1./eeta(ip,j,k )
     .                 +1./eeta(i ,j,km)+1./eeta(i ,j,k ))
     .             *g_sub_elem(1,3,xp,y0,zm,.false.)
        d_xz_imkp = 4./(1./eeta(im,j,kp)+1./eeta(im,j,k )
     .                 +1./eeta(i ,j,kp)+1./eeta(i ,j,k ))
     .             *g_sub_elem(1,3,xm,y0,zp,.false.)
        d_xz_imkm = 4./(1./eeta(im,j,km)+1./eeta(im,j,k )
     .                 +1./eeta(i ,j,km)+1./eeta(i ,j,k ))
     .             *g_sub_elem(1,3,xm,y0,zm,.false.)

        d_xy_kp = 2./(1./eeta(i,j,kp)+1./eeta(i,j,k))
     .           *g_sub_elem(1,2,x0,y0,zp,.false.)
        d_xy_km = 2./(1./eeta(i,j,km)+1./eeta(i,j,k))
     .           *g_sub_elem(1,2,x0,y0,zm,.false.)

        d_yz_ip = 2./(1./eeta(ip,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(2,3,xp,y0,z0,.false.)
        d_yz_im = 2./(1./eeta(im,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(2,3,xm,y0,z0,.false.)

        d_xz_kp = 2./(1./eeta(i,j,kp)+1./eeta(i,j,k))
     .           *g_sub_elem(1,3,x0,y0,zp,.false.)
        d_xz_km = 2./(1./eeta(i,j,km)+1./eeta(i,j,k))
     .           *g_sub_elem(1,3,x0,y0,zm,.false.)

        d_xz_ip = 2./(1./eeta(ip,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,3,xp,y0,z0,.false.)
        d_xz_im = 2./(1./eeta(im,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,3,xm,y0,z0,.false.)

        dwa = -dzh(kg)*dyh(jg)
     .              *(d_zz_ip*(ay(ip,j,k)-ay(i ,j,k))/dx(ig  )
     .               -d_zz_im*(ay(i ,j,k)-ay(im,j,k))/dx(ig-1) )
        dwb = -dxh(ig)*dyh(jg)
     .              *(d_xx_kp*(ay(i,j,kp)-ay(i,j,k ))/dz(kg  )
     .               -d_xx_km*(ay(i,j,k )-ay(i,j,km))/dz(kg-1) )
        dwc = -dyh(jg)/2.
     .              *(-d_xz_ipkp*(ay(ip,j,kp)-ay(i ,j,k ))
     .                +d_xz_imkm*(ay(i ,j,k )-ay(im,j,km))
     .                +d_xz_ipkm*(ay(ip,j,km)-ay(i ,j,k ))
     .                -d_xz_imkp*(ay(i ,j,k )-ay(im,j,kp)))
        dwd =  dxh(ig)/4.
     .              *(d_xx_kp*(az(i,jp,kp)-az(i,jm,kp)
     .                        +az(i,jp,k )-az(i,jm,k ))
     .               -d_xx_km*(az(i,jp,km)-az(i,jm,km)
     .                        +az(i,jp,k )-az(i,jm,k )))
        dwe = -dyh(jg)/4.
     .              *(d_xy_kp*(az(ip,j,kp)-az(im,j,kp)
     .                        +az(ip,j,k )-az(im,j,k ))
     .               -d_xy_km*(az(ip,j,km)-az(im,j,km)
     .                        +az(ip,j,k )-az(im,j,k )))
        dwf =  dxh(ig)*dyh(jg)
     .              *(d_xy_kp*(ax(i,j,kp)-ax(i,j,k ))/dz(kg)
     .               -d_xy_km*(ax(i,j,k )-ax(i,j,km))/dz(kg-1))
        dwg = -dxh(ig)/4.
     .              *(d_xz_kp*(ax(i,jp,kp)-ax(i,jm,kp)
     .                        +ax(i,jp,k )-ax(i,jm,k ))
     .               -d_xz_km*(ax(i,jp,km)-ax(i,jm,km)
     .                        +ax(i,jp,k )-ax(i,jm,k )))
        dwh =  dzh(kg)/4.
     .              *(d_zz_ip*(ax(ip,jp,k)-ax(ip,jm,k)
     .                        +ax(i ,jp,k)-ax(i ,jm,k))
     .               -d_zz_im*(ax(im,jp,k)-ax(im,jm,k)
     .                        +ax(i ,jp,k)-ax(i ,jm,k)))
        dwi =  dyh(jg)*dzh(kg)
     .              *(d_yz_ip*(az(ip,j,k)-az(i ,j,k))/dx(ig  )
     .               -d_yz_im*(az(i ,j,k)-az(im,j,k))/dx(ig-1))
        dwj = -dzh(kg)/4.
     .              *(d_xz_ip*(az(ip,jp,k)-az(ip,jm,k)
     .                        +az(i ,jp,k)-az(i ,jm,k))
     .               -d_xz_im*(az(im,jp,k)-az(im,jm,k)
     .                        +az(i ,jp,k)-az(i ,jm,k)))
        dwk = -dyh(jg)/4.
     .              *(d_yz_ip*(ax(ip,j,kp)-ax(ip,j,km)
     .                        +ax(i ,j,kp)-ax(i ,j,km))
     .               -d_yz_im*(ax(im,j,kp)-ax(im,j,km)
     .                        +ax(i ,j,kp)-ax(i ,j,km)))

      case(3)

        d_xx_jp = 2./(1./eeta(i,jp,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,1,x0,yp,z0,.false.)
        d_xx_jm = 2./(1./eeta(i,jm,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,1,x0,ym,z0,.false.)

        d_yy_ip = 2./(1./eeta(ip,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(2,2,xp,y0,z0,.false.)
        d_yy_im = 2./(1./eeta(im,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(2,2,xm,y0,z0,.false.)

        d_xy_ipjp = 4./(1./eeta(ip,jp,k)+1./eeta(ip,j,k )
     .                 +1./eeta(i ,jp,k)+1./eeta(i ,j,k ))
     .             *g_sub_elem(1,2,xp,yp,z0,.false.)
        d_xy_ipjm = 4./(1./eeta(ip,jm,k)+1./eeta(ip,j,k )
     .                 +1./eeta(i ,jm,k)+1./eeta(i ,j,k ))
     .             *g_sub_elem(1,2,xp,ym,z0,.false.)
        d_xy_imjp = 4./(1./eeta(im,jp,k)+1./eeta(im,j,k )
     .                 +1./eeta(i ,jp,k)+1./eeta(i ,j,k ))
     .             *g_sub_elem(1,2,xm,yp,z0,.false.)
        d_xy_imjm = 4./(1./eeta(im,jm,k)+1./eeta(im,j,k )
     .                 +1./eeta(i ,jm,k)+1./eeta(i ,j,k ))
     .             *g_sub_elem(1,2,xm,ym,z0,.false.)

        d_xy_jp = 2./(1./eeta(i,jp,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,2,x0,yp,z0,.false.)
        d_xy_jm = 2./(1./eeta(i,jm,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,2,x0,ym,z0,.false.)

        d_yz_ip = 2./(1./eeta(ip,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(2,3,xp,y0,z0,.false.)
        d_yz_im = 2./(1./eeta(im,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(2,3,xm,y0,z0,.false.)

        d_xy_ip = 2./(1./eeta(ip,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,2,xp,y0,z0,.false.)
        d_xy_im = 2./(1./eeta(im,j,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,2,xm,y0,z0,.false.)

        d_xz_jp = 2./(1./eeta(i,jp,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,3,x0,yp,z0,.false.)
        d_xz_jm = 2./(1./eeta(i,jm,k)+1./eeta(i,j,k))
     .           *g_sub_elem(1,3,x0,ym,z0,.false.)

        dwa = -dzh(kg)*dyh(jg)
     .              *(d_yy_ip*(az(ip,j,k)-az(i ,j,k))/dx(ig  )
     .               -d_yy_im*(az(i ,j,k)-az(im,j,k))/dx(ig-1) )
        dwb = -dxh(ig)*dzh(kg)
     .              *(d_xx_jp*(az(i,jp,k)-az(i,j ,k))/dy(jg  )
     .               -d_xx_jm*(az(i,j ,k)-az(i,jm,k))/dy(jg-1) )
        dwc = -dzh(kg)/2.
     .              *(-d_xy_ipjp*(az(ip,jp,k)-az(i ,j ,k))
     .                +d_xy_imjm*(az(i ,j ,k)-az(im,jm,k))
     .                +d_xy_ipjm*(az(ip,jm,k)-az(i ,j ,k))
     .                -d_xy_imjp*(az(i ,j ,k)-az(im,jp,k)))
        dwd =  dxh(ig)/4.
     .              *(d_xx_jp*(ay(i,jp,kp)-ay(i,jp,km)
     .                        +ay(i,j ,kp)-ay(i,j ,km))
     .               -d_xx_jm*(ay(i,jm,kp)-ay(i,jm,km)
     .                        +ay(i,j ,kp)-ay(i,j ,km)))
        dwe=  -dxh(ig)/4.
     .              *(d_xy_jp*(ax(i,jp,kp)-ax(i,jp,km)
     .                        +ax(i,j ,kp)-ax(i,j ,km))
     .               -d_xy_jm*(ax(i,jm,kp)-ax(i,jm,km)
     .                        +ax(i,j ,kp)-ax(i,j ,km)))
        dwf =  dxh(ig)*dzh(kg)
     .              *(d_xz_jp*(ax(i,jp,k)-ax(i,j ,k))/dy(jg)
     .               -d_xz_jm*(ax(i,j ,k)-ax(i,jm,k))/dy(jg-1))
        dwg = -dzh(kg)/4.
     .              *(d_xz_jp*(ay(ip,jp,k)-ay(im,jp,k)
     .                        +ay(ip,j ,k)-ay(im,j ,k))
     .               -d_xz_jm*(ay(ip,jm,k)-ay(im,jm,k)
     .                        +ay(ip,j ,k)-ay(im,j ,k)))
        dwh =  dyh(jg)/4.
     .              *(d_yy_ip*(ax(ip,j,kp)-ax(ip,j,km)
     .                        +ax(i ,j,kp)-ax(i ,j,km))
     .               -d_yy_im*(ax(im,j,kp)-ax(im,j,km)
     .                        +ax(i ,j,kp)-ax(i ,j,km)))
        dwi =  dyh(jg)*dzh(kg)
     .              *(d_yz_ip*(ay(ip,j,k)-ay(i ,j,k))/dx(ig  )
     .               -d_yz_im*(ay(i ,j,k)-ay(im,j,k))/dx(ig-1))
        dwj = -dzh(kg)/4.
     .              *(d_yz_ip*(ax(ip,jp,k)-ax(ip,jm,k)
     .                        +ax(i ,jp,k)-ax(i ,jm,k))
     .               -d_yz_im*(ax(im,jp,k)-ax(im,jm,k)
     .                        +ax(i ,jp,k)-ax(i ,jm,k)))
        dwk = -dyh(jg)/4.
     .              *(d_xy_ip*(ay(ip,j,kp)-ay(ip,j,km)
     .                        +ay(i ,j,kp)-ay(i ,j,km))
     .               -d_xy_im*(ay(im,j,kp)-ay(im,j,km)
     .                        +ay(i ,j,kp)-ay(i ,j,km)))
      case default

        write (*,*) 'Error in component in curlcurl'
        write (*,*) 'Aborting...'
        stop

      end select

      curlcurl = dwa+dwb+dwc+dwd+dwe+dwf+dwg+dwh+dwi+dwj+dwk

cc      if (i == 0 .or. j == 0 .or. k == 0 ) then
cc        write (*,*) 'Error in laplace; i,j,k=0'
cc      elseif (i == nx+1 .or. j == ny+1 .or. k == nz+1) then
cc        write (*,*) 'Error in laplace; i,j,k=nmax+1'
cc      elseif (.not.sing_point) then
cc        d_xx_ip = g_super_elem(1,1,0.5*(xx(ig+1)+xx(ig)),yy(jg),zz(kg))
cc        d_xx_im = g_super_elem(1,1,0.5*(xx(ig-1)+xx(ig)),yy(jg),zz(kg))
cc        d_yy_jp = g_super_elem(2,2,xx(ig),0.5*(yy(jg+1)+yy(jg)),zz(kg))
cc        d_yy_jm = g_super_elem(2,2,xx(ig),0.5*(yy(jg-1)+yy(jg)),zz(kg))
cc        d_zz_kp = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg+1)+zz(kg)))
cc        d_zz_km = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg-1)+zz(kg)))
cc
cc        d_xy_ipjp = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
cc     .                              ,0.5*(yy(jg+1)+yy(jg)),zz(kg))
cc        d_xy_ipjm = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
cc     .                              ,0.5*(yy(jg-1)+yy(jg)),zz(kg))
cc        d_xy_imjp = g_super_elem(1,2,0.5*(xx(ig-1)+xx(ig))
cc     .                              ,0.5*(yy(jg+1)+yy(jg)),zz(kg))
cc        d_xy_imjm = g_super_elem(1,2,0.5*(xx(ig-1)+xx(ig))
cc     .                              ,0.5*(yy(jg-1)+yy(jg)),zz(kg))
cc                 
cc        d_xz_ipkp = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
cc     .                              ,0.5*(zz(kg+1)+zz(kg)))
cc        d_xz_ipkm = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
cc     .                              ,0.5*(zz(kg-1)+zz(kg)))
cc        d_xz_imkp = g_super_elem(1,3,0.5*(xx(ig-1)+xx(ig)),yy(jg)
cc     .                              ,0.5*(zz(kg+1)+zz(kg)))
cc        d_xz_imkm = g_super_elem(1,3,0.5*(xx(ig-1)+xx(ig)),yy(jg)
cc     .                              ,0.5*(zz(kg-1)+zz(kg)))
cc                 
cc        d_yz_jpkp = g_super_elem(2,3,xx(ig),0.5*(yy(jg+1)+yy(jg))
cc     .                                     ,0.5*(zz(kg+1)+zz(kg)))
cc        d_yz_jpkm = g_super_elem(2,3,xx(ig),0.5*(yy(jg+1)+yy(jg))
cc     .                                     ,0.5*(zz(kg-1)+zz(kg)))
cc        d_yz_jmkp = g_super_elem(2,3,xx(ig),0.5*(yy(jg-1)+yy(jg))
cc     .                                     ,0.5*(zz(kg+1)+zz(kg)))
cc        d_yz_jmkm = g_super_elem(2,3,xx(ig),0.5*(yy(jg-1)+yy(jg))
cc     .                                     ,0.5*(zz(kg-1)+zz(kg)))
cc
cc      else
cc
cc        d_xx_ip = g_super_elem(1,1,0.5*(xx(ig+1)+xx(ig)),yy(jg),zz(kg))
cc        d_xx_im = 0d0
cc        d_yy_jp = 0d0
cc        d_yy_jm = 0d0
cc        d_zz_kp = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg+1)+zz(kg)))
cc        d_zz_km = g_super_elem(3,3,xx(ig),yy(jg),0.5*(zz(kg-1)+zz(kg)))
cc
cc        d_xy_ipjp = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
cc     .                              ,0.5*(yy(jg+1)+yy(jg)),zz(kg))
cc        d_xy_ipjm = g_super_elem(1,2,0.5*(xx(ig+1)+xx(ig))
cc     .                              ,0.5*(yy(jg-1)+yy(jg)),zz(kg))
cc        d_xy_imjp = 0d0
cc        d_xy_imjm = 0d0
cc                 
cc        d_xz_ipkp = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
cc     .                              ,0.5*(zz(kg+1)+zz(kg)))
cc        d_xz_ipkm = g_super_elem(1,3,0.5*(xx(ig+1)+xx(ig)),yy(jg)
cc     .                              ,0.5*(zz(kg-1)+zz(kg)))
cc        d_xz_imkp = 0d0
cc        d_xz_imkm = 0d0
cc                 
cc        d_yz_jpkp = g_super_elem(2,3,xx(ig),0.5*(yy(jg+1)+yy(jg))
cc     .                                     ,0.5*(zz(kg+1)+zz(kg)))
cc        d_yz_jpkm = g_super_elem(2,3,xx(ig),0.5*(yy(jg+1)+yy(jg))
cc     .                                     ,0.5*(zz(kg-1)+zz(kg)))
cc        d_yz_jmkp = g_super_elem(2,3,xx(ig),0.5*(yy(jg-1)+yy(jg))
cc     .                                     ,0.5*(zz(kg+1)+zz(kg)))
cc        d_yz_jmkm = g_super_elem(2,3,xx(ig),0.5*(yy(jg-1)+yy(jg))
cc     .                                     ,0.5*(zz(kg-1)+zz(kg)))
cc
cc      endif
cc
cc      curlcurl =
cc     .     dyh(jg)*dzh(kg)*( d_xx_ip*(arr(ip,j,k)-arr(i,j,k))/dx(ig)
cc     .                      +d_xx_im*(arr(im,j,k)-arr(i,j,k))/dx(ig-1) )
cc     .    +dxh(ig)*dzh(kg)*( d_yy_jp*(arr(i,jp,k)-arr(i,j,k))/dy(jg)
cc     .                      +d_yy_jm*(arr(i,jm,k)-arr(i,j,k))/dy(jg-1) )
cc     .    +dxh(ig)*dyh(jg)*( d_zz_kp*(arr(i,j,kp)-arr(i,j,k))/dz(kg)
cc     .                      +d_zz_km*(arr(i,j,km)-arr(i,j,k))/dz(kg-1) )
cc     .    +0.5*
cc     .      ( dzh(kg)*( d_xy_ipjp*(arr(ip,jp,k)-arr(i,j,k))
cc     .                 +d_xy_imjm*(arr(im,jm,k)-arr(i,j,k))
cc     .                 -d_xy_ipjm*(arr(ip,jm,k)-arr(i,j,k))
cc     .                 -d_xy_imjp*(arr(im,jp,k)-arr(i,j,k)) )
cc     .       +dyh(jg)*( d_xz_ipkp*(arr(ip,j,kp)-arr(i,j,k))
cc     .                 +d_xz_imkm*(arr(im,j,km)-arr(i,j,k))
cc     .                 -d_xz_ipkm*(arr(ip,j,km)-arr(i,j,k))
cc     .                 -d_xz_imkp*(arr(im,j,kp)-arr(i,j,k)) )
cc     .       +dxh(jg)*( d_yz_jpkp*(arr(i,jp,kp)-arr(i,j,k))
cc     .                 +d_yz_jmkm*(arr(i,jm,km)-arr(i,j,k))
cc     .                 -d_yz_jpkm*(arr(i,jp,km)-arr(i,j,k))
cc     .                 -d_yz_jmkp*(arr(i,jm,kp)-arr(i,j,k)) ) )

c     End program

      end function curlcurl

c     curlcurl2
c     ###############################################################
      real*8 function curlcurl2(i,j,k,ax,ay,az,comp)

c     ---------------------------------------------------------------
c     Calculates curl(eta curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      use grid

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,comp

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: x0,xp,xm,y0,yp,ym,z0,zp,zm

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      select case(comp)
      case(1)

        curlcurl2 = dzh(kg)*(ay(ip,jp,k)-ay(im,jp,k)
     .                      -ay(ip,jm,k)+ay(im,jm,k))/4.
     .             -dxh(ig)*dzh(kg)
     .                     *(ax(i,jp,k)-2*ax(i,j,k)+ax(i,jm,k))/dyh(jg)

      case(2)

        curlcurl2 = dzh(kg)*(ax(ip,jp,k)-ax(im,jp,k)
     .                      -ax(ip,jm,k)+ax(im,jm,k))/4.
     .             -dyh(jg)*dzh(kg)
     .                     *(ay(ip,j,k)-2*ay(i,j,k)+ay(im,j,k))/dxh(ig)

      case(3)

        curlcurl2 =-dyh(jg)*dzh(kg)
     .                     *(az(ip,j,k)-2*az(i,j,k)+az(im,j,k))/dxh(ig)
     .             -dxh(ig)*dzh(kg)
     .                     *(az(i,jp,k)-2*az(i,j,k)+az(i,jm,k))/dyh(jg)

      case default

        write (*,*) 'Error in component in curlcurl'
        write (*,*) 'Aborting...'
        stop

      end select

      curlcurl2 = eeta(i,j,k)*curlcurl2

c     End program

      end function curlcurl2

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
