c module equilibrium
c ######################################################################
      module equilibrium

        integer(4) :: IRHO,IVX,IVY,IVZ,IBX,IBY,IBZ,ITMP,IJX,IJY,IJZ
        parameter    (IRHO=1,IVX=2,IVY=3,IVZ=4,IBX=5,IBY=6,IBZ=7,ITMP=8
     .               ,IJX=9,IJY=10,IJZ=11)

        real(8)    :: dlambda,rshear,vparflow,vperflow

        character*(5) :: equil

      end module equilibrium

c module transport_params
c ######################################################################
      module transport_params

        real(8) :: nu,eta,dd,chi,gamma

      end module transport_params

c module grid_aliases
c ######################################################################
      module grid_aliases

        real(8),pointer,dimension(:) :: xx,yy,zz,dxh,dyh,dzh,dx,dy,dz

        integer(4) :: igx,igy,igz,nx,ny,nz

      end module grid_aliases

c module auxiliaryVariables
c ######################################################################
      module auxiliaryVariables

cc        use grid_aliases

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_cov,by_cov,bz_cov
     .         ,jx,jy,jz,jx_cov,jy_cov,jz_cov
     .         ,vx,vy,vz,vx_cov,vy_cov,vz_cov
     .         ,eeta,nuu

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_car,by_car,bz_car
     .         ,jx_car,jy_car,jz_car
     .         ,vx_car,vy_car,vz_car
     .         ,divrgB,divrgV,divrgJ,Pflux,qfactor,p_tot

        real(8),pointer,dimension(:,:,:):: rho,rvx,rvy,rvz,bx,by,bz,tmp

      end module auxiliaryVariables

c module operators
c ######################################################################
      module operators

        use grid

        use grid_aliases

        use generalOperators

      contains

c     div
c     ###############################################################
      real(8) function div(i,j,k,ax,ay,az)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: i,j,k
      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ig,jg,kg
      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacp,jac0,xh,yh,zh,jach
     .             ,xip,yip,zip
      logical    :: cartesian

c     Begin program

      call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x0,y0,z0,cartesian)

      jac0 = jacobian(x0 ,y0 ,z0 ,cartesian)

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      if (i == 1 .and. bcond(1) == SP) then
        call getCoordinates(i+1,j,k,igx,igy,igz,ig,jg,kg,xip,yip,zip
     .                     ,cartesian)

        jacp = jacobian(xip,yip,zip,cartesian)

        xh = (xip+x0)/2.
        yh = (yip+y0)/2.
        zh = (zip+z0)/2.
        jach = jacobian(xh ,yh ,zh ,cartesian)
cc        jach = 0.5*(jacp+jac0)   !Only good for cylindrical-like geom.

        div =  ((ax(i+1,j  ,k  )/jacp
     .          +ax(i  ,j  ,k  )/jac0)*jach)    /2./dxx
     .        +(ay(i  ,j+1,k  )-ay(i  ,j-1,k  ))/2./dyy
     .        +(az(i  ,j  ,k+1)-az(i  ,j  ,k-1))/2./dzz
      else
        div =  (ax(i+1,j  ,k  )-ax(i-1,j  ,k  ))/2./dxx
     .        +(ay(i  ,j+1,k  )-ay(i  ,j-1,k  ))/2./dyy
     .        +(az(i  ,j  ,k+1)-az(i  ,j  ,k-1))/2./dzz
      endif

      div = div/jac0
      
c     End 

      end function div

c     laplacian
c     ###############################################################
      real*8 function laplacian(i,j,k,arr)

c     ---------------------------------------------------------------
c     Calculates dvol*lap(arr) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

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

      logical    :: sing_point

c     Begin program

      sing_point = .false.
      if (i == 1 .and. bcond(1) == SP) sing_point = .true.

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
     .                        ,.false.)
        d_xx_im = 0d0

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
     .    +0.5*( dzh(kg)*( d_xy_ipjp*(arr(ip,jp,k)-arr(i,j,k))
     .                    +d_xy_imjm*(arr(im,jm,k)-arr(i,j,k))
     .                    -d_xy_ipjm*(arr(ip,jm,k)-arr(i,j,k))
     .                    -d_xy_imjp*(arr(im,jp,k)-arr(i,j,k)) )
     .          +dyh(jg)*( d_xz_ipkp*(arr(ip,j,kp)-arr(i,j,k))
     .                    +d_xz_imkm*(arr(im,j,km)-arr(i,j,k))
     .                    -d_xz_ipkm*(arr(ip,j,km)-arr(i,j,k))
     .                    -d_xz_imkp*(arr(im,j,kp)-arr(i,j,k)) )
     .          +dxh(ig)*( d_yz_jpkp*(arr(i,jp,kp)-arr(i,j,k))
     .                    +d_yz_jmkm*(arr(i,jm,km)-arr(i,j,k))
     .                    -d_yz_jpkm*(arr(i,jp,km)-arr(i,j,k))
     .                    -d_yz_jmkp*(arr(i,jm,kp)-arr(i,j,k)) ) )

c     End program

      end function laplacian

c     veclaplacian
c     ###############################################################
      function veclaplacian(i,j,k,ax,ay,az,diff,alteom,icomp)
     .         result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*lap(vector) at cell centers in general non-orthog.
c     coordinates.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,icomp

      real(8)    :: ax  (0:nx+1,0:ny+1,0:nz+1)
     .             ,ay  (0:nx+1,0:ny+1,0:nz+1)
     .             ,az  (0:nx+1,0:ny+1,0:nz+1)
     .             ,diff(0:nx+1,0:ny+1,0:nz+1)

      real(8)    :: vlap

      logical    :: alteom

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: xim,yim,zim,xip,yip,zip
     .             ,xjm,yjm,zjm,xjp,yjp,zjp
     .             ,xkm,ykm,zkm,xkp,ykp,zkp

      real(8)    :: dvol,dS1,dS2,dS3,dxx,dyy,dzz,x0,y0,z0,jac

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o

      real(8)    :: nabla_v(3,3),gsuper(3,3)
     .             ,hess1(3,3),hess2(3,3),hess3(3,3),msource

      logical    :: sing_point,cartesian

c     Begin program

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      sing_point = .false.
      if (i == 1 .and. bcond(1) == SP) sing_point = .true.

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

      jac  = jacobian(x0,y0,z0,cartesian)

      if (coords /= 'car') then
        hess1 = hessian(1,x0,y0,z0,cartesian)
        hess2 = hessian(2,x0,y0,z0,cartesian)
        hess3 = hessian(3,x0,y0,z0,cartesian)
      endif

      call nabtensor_x(i ,j,k,t11p,t12p,t13p, 1)
      call nabtensor_x(im,j,k,t11m,t12m,t13m,-1)
      if (coords /= 'car') call nabtensor_x(i,j,k,t11o,t12o,t13o, 0)

      call nabtensor_y(i,j ,k,t21p,t22p,t23p, 1)
      call nabtensor_y(i,jm,k,t21m,t22m,t23m,-1)
      if (coords /= 'car') call nabtensor_y(i,j,k,t21o,t22o,t23o, 0)

      call nabtensor_z(i,j,k ,t31p,t32p,t33p, 1)
      call nabtensor_z(i,j,km,t31m,t32m,t33m,-1)
      if (coords /= 'car') call nabtensor_z(i,j,k,t31o,t32o,t33o, 0)

      msource = 0d0

      select case (icomp)
      case(1)

        flxip = t11p
        flxim = t11m

        flxjp = t21p
        flxjm = t21m

        flxkp = t31p
        flxkm = t31m

        call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                              ,flxkp,flxkm,bcond)

        if (coords /= 'car') then
          msource =dvol
     .            *(t11o*hess1(1,1)+t12o*hess1(1,2)+t13o*hess1(1,3)
     .             +t21o*hess1(2,1)+t22o*hess1(2,2)+t23o*hess1(2,3)
     .             +t31o*hess1(3,1)+t32o*hess1(3,2)+t33o*hess1(3,3))
        endif

      case(2)

        flxip = t12p
        flxim = t12m

        flxjp = t22p
        flxjm = t22m

        flxkp = t32p
        flxkm = t32m

        if (coords /= 'car') then
          if (alteom) then

            call imposeBConfluxes (i,j,k,flxip,flxim
     .                                  ,flxjp,flxjm
     .                                  ,flxkp,flxkm,bcond)

            msource=dvol
     .             *(t11o*hess2(1,1)+t12o*hess2(1,2)+t13o*hess2(1,3)
     .              +t21o*hess2(2,1)+t22o*hess2(2,2)+t23o*hess2(2,3)
     .              +t31o*hess2(3,1)+t32o*hess2(3,2)+t33o*hess2(3,3)
     .              -t12o*(hess1(1,1)+hess2(1,2)+hess3(1,3))
     .              -t22o*(hess1(2,1)+hess2(2,2)+hess3(2,3))
     .              -t32o*(hess1(3,1)+hess2(3,2)+hess3(3,3)))

            flxip = flxip/jac
            flxim = flxim/jac
                         
            flxjp = flxjp/jac
            flxjm = flxjm/jac
                         
            flxkp = flxkp/jac
            flxkm = flxkm/jac

          else
            msource=dvol
     .             *(t11o*hess2(1,1)+t12o*hess2(1,2)+t13o*hess2(1,3)
     .              +t21o*hess2(2,1)+t22o*hess2(2,2)+t23o*hess2(2,3)
     .              +t31o*hess2(3,1)+t32o*hess2(3,2)+t33o*hess2(3,3))
          endif
        endif

      case(3)

        flxip = t13p
        flxim = t13m

        flxjp = t23p
        flxjm = t23m

        flxkp = t33p
        flxkm = t33m

        call imposeBConfluxes (i,j,k,flxip,flxim,flxjp,flxjm
     .                              ,flxkp,flxkm,bcond)

        if (coords /= 'car') then
          msource =dvol
     .           *(t11o*hess3(1,1)+t12o*hess3(1,2)+t13o*hess3(1,3)
     .            +t21o*hess3(2,1)+t22o*hess3(2,2)+t23o*hess3(2,3)
     .            +t31o*hess3(3,1)+t32o*hess3(3,2)+t33o*hess3(3,3))
        endif

      end select

      vlap = jac*( dS1*(flxip - flxim)
     .           + dS2*(flxjp - flxjm)
     .           + dS3*(flxkp - flxkm) ) -msource

c     End program

      contains

c     nabtensor_x
c     #############################################################
      subroutine nabtensor_x(i,j,k,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag
        real(8)    :: t11,t12,t13

c     Local variables

        integer(4) :: ig,jg,kg,ip
        real(8)    :: x,y,z,jac,jac0,jacp,vis

c     Begin program

        ip = i+1

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
          if (.not.alteom .and. sing_point) then
            ip = i
            x = xim
            y = yim
            z = zim
          else
            x = (xim+x0)/2.
            y = (yim+y0)/2.
            z = (zim+z0)/2.
          endif
        endif

        jac    = jacobian(x,y,z,cartesian)

        gsuper = g_super (x,y,z,cartesian)

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,x,y,z,ax,ay,az,cartesian,1)
        else
          nabla_v = fnabla_v(i,j,k,x,y,z,ax,ay,az,cartesian,0)
        endif

        vis = 2./( 1./diff(ip,j,k) + 1./diff(i,j,k) )

        t11 =vis*( gsuper(1,1)*nabla_v(1,1)
     .            +gsuper(1,2)*nabla_v(2,1)
     .            +gsuper(1,3)*nabla_v(3,1) )

        t12 =vis*( gsuper(1,1)*nabla_v(1,2)
     .            +gsuper(1,2)*nabla_v(2,2)
     .            +gsuper(1,3)*nabla_v(3,2) )

        t13 =vis*( gsuper(1,1)*nabla_v(1,3)
     .            +gsuper(1,2)*nabla_v(2,3)
     .            +gsuper(1,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t11 = t11/jac
          if (.not.alteom) t12 = t12/jac
          t13 = t13/jac
        endif

cc        t11 = 0d0
cc        t12 = 0d0
cc        t13 = 0d0
c     End program

      end subroutine nabtensor_x

c     nabtensor_y
c     #############################################################
      subroutine nabtensor_y(i,j,k,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag
        real(8)    :: t21,t22,t23

c     Local variables

        integer(4) :: ig,jg,kg,jp
        real(8)    :: x,y,z,jac,vis

c     Begin program

        jp = j+1

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
          nabla_v = fnabla_v(i,j,k,x,y,z,ax,ay,az,cartesian,2)
        else
          nabla_v = fnabla_v(i,j,k,x,y,z,ax,ay,az,cartesian,0)
        endif

        vis = 2./( 1./diff(i,jp,k) + 1./diff(i,j,k) )

        t21 =vis*( gsuper(2,1)*nabla_v(1,1)
     .            +gsuper(2,2)*nabla_v(2,1)
     .            +gsuper(2,3)*nabla_v(3,1) )

        t22 =vis*( gsuper(2,1)*nabla_v(1,2)
     .            +gsuper(2,2)*nabla_v(2,2)
     .            +gsuper(2,3)*nabla_v(3,2) )

        t23 =vis*( gsuper(2,1)*nabla_v(1,3)
     .            +gsuper(2,2)*nabla_v(2,3)
     .            +gsuper(2,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t21 = t21/jac
          if (.not.alteom) t22 = t22/jac
          t23 = t23/jac
        endif

cc        t21 = 0d0
cc        t22 = 0d0
cc        t23 = 0d0
c     End program

      end subroutine nabtensor_y

c     nabtensor_z
c     #############################################################
      subroutine nabtensor_z(i,j,k,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag
        real(8)    :: t31,t32,t33

c     Local variables

        integer(4) :: ig,jg,kg,kp
        real(8)    :: x,y,z,jac,vis

c     Begin program

        kp=k+1

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
          nabla_v = fnabla_v(i,j,k,x,y,z,ax,ay,az,cartesian,3)
        else
          nabla_v = fnabla_v(i,j,k,x,y,z,ax,ay,az,cartesian,0)
        endif

        vis = 2./( 1./diff(i,j,kp) + 1./diff(i,j,k) )

        t31 =vis*( gsuper(3,1)*nabla_v(1,1)
     .            +gsuper(3,2)*nabla_v(2,1)
     .            +gsuper(3,3)*nabla_v(3,1) )

        t32 =vis*( gsuper(3,1)*nabla_v(1,2)
     .            +gsuper(3,2)*nabla_v(2,2)
     .            +gsuper(3,3)*nabla_v(3,2) )

        t33 =vis*( gsuper(3,1)*nabla_v(1,3)
     .            +gsuper(3,2)*nabla_v(2,3)
     .            +gsuper(3,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alteom) t32 = t32/jac
          t33 = t33/jac
        endif

cc        t31 = 0d0
cc        t32 = 0d0
cc        t33 = 0d0

c     End program

      end subroutine nabtensor_z

      end function veclaplacian

c     curl
c     ###############################################################
      real*8 function curl(i,j,k,ax,ay,az,comp)

c     ---------------------------------------------------------------
c     Calculates curl(A) at cell centers in general non-orthogonal
c     coordinates. The vector A is covariant, and returns the 
c     contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,comp

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg
      real(8)    :: dx,dy,dz,flxip,flxim,flxjp,flxjm,flxkp,flxkm
      real(8)    :: dS1,dS2,dS3,x0,y0,z0,jac
     .             ,xip,yip,zip,jacp,xh,yh,zh,jach
      logical    :: cartesian

c     Begin program

      call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x0,y0,z0,cartesian)

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
        if (i == 1 .and. bcond(1) == SP) then
          flxim =-az(im,j,k)
        else
          flxim =-0.5*(az(im,j,k)+az(i,j,k))
        endif

        flxkp = 0.5*(ax(i,j,kp)+ax(i,j,k))
        flxkm = 0.5*(ax(i,j,km)+ax(i,j,k))

      case(3)

        flxkp = 0d0
        flxkm = 0d0

        if (i == 1 .and. bcond(1) == SP) then
          call getCoordinates(i+1,j,k,igx,igy,igz,ig,jg,kg,xip,yip,zip
     .                       ,cartesian)
          jacp = jacobian(xip,yip,zip,cartesian)
          jac  = jacobian(x0 ,y0 ,z0 ,cartesian)

          xh = (xip+x0)/2.
          yh = (yip+y0)/2.
          zh = (zip+z0)/2.
          jach = jacobian(xh ,yh ,zh ,cartesian)
cc          jach = 0.5*(jacp+jac) !Only good for cylindrical-like geom.

          flxip = 0.5*(ay(ip,j,k)/jacp+ay(i,j,k)/jac)*jach
          flxim = ay(im,j,k)

cc        elseif (i == 2 .and. bcond(1) == SP) then
cc          call getCoordinates(i-1,j,k,igx,igy,igz,ig,jg,kg,xim,yim,zim
cc     .                       ,cartesian)
cc          xh = (xim+x0)/2.
cc          yh = (yim+y0)/2.
cc          zh = (zim+z0)/2.
cc          jacm = jacobian(xim,yim,zim,cartesian)
cc          jach = jacobian(xh ,yh ,zh ,cartesian)
cc
cc          flxip = 0.5*(ay(ip,j,k)+ay(i,j,k))
cc          flxim = 0.5*(ay(im,j,k)/jacm+ay(i,j,k)/jac)*jach
cccc          flxim = 0.5*(ay(im,j,k)+ay(i,j,k))
        else
          flxip = 0.5*(ay(ip,j,k)+ay(i,j,k))
          flxim = 0.5*(ay(im,j,k)+ay(i,j,k))
        endif

        flxjp =-0.5*(ax(i,jp,k)+ax(i,j,k))
        flxjm =-0.5*(ax(i,jm,k)+ax(i,j,k))

      case default

        write (*,*) 'Error in component in curl'
        write (*,*) 'Aborting...'
        stop

      end select

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

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,comp

      real(8)       :: ax(0:nx+1,0:ny+1,0:nz+1)
     .                ,ay(0:nx+1,0:ny+1,0:nz+1)
     .                ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg
      real(8)    :: dx,dy,dz,dx1,dx2,dy1,dy2,dz1,dz2
      real(8)    :: daxdz,dazdx,daydx,daxdy,dazdy,daydz
      real(8)    :: jac0,jacp,jacm,x0,y0,z0,xp,yp,zp,xm,ym,zm
      logical    :: cartesian

c     Begin program

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      if (i == 1 .and. bcond(1) == SP) then

        call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x0,y0,z0
     .                     ,cartesian)
        call getCoordinates(ip,j,k,igx,igy,igz,ig,jg,kg,xp,yp,zp
     .                     ,cartesian)
        xm = 0.5*(x0+xp)
        ym = 0.5*(y0+yp)
        zm = 0.5*(z0+zp)

        jacp = jacobian(xp,yp,zp,cartesian)
        jac0 = jacobian(x0,y0,z0,cartesian)
        jacm = jacobian(xm,ym,zm,cartesian)

      else
        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
      endif

      select case(comp)
      case(1)

        if (j==0) then
          dy1=grid_params%dy(jg)
          dy2=grid_params%dy(jg+1)
          dazdy = (-az(i,j+2,k)*dy1/dy2/(dy1+dy2)
     .             +az(i,j+1,k)*(dy1+dy2)/dy1/dy2
     .             -az(i,j  ,k)*(1./dy1+1./(dy1+dy2)))
c1st          dazdy = (az(i,jp,k)-az(i,j,k))/dy1
        elseif (j==ny+1) then
          dy1=grid_params%dy(jg-1)
          dy2=grid_params%dy(jg-2)
          dazdy = -(-az(i,j-2,k)*dy1/dy2/(dy1+dy2)
     .              +az(i,j-1,k)*(dy1+dy2)/dy1/dy2
     .              -az(i,j  ,k)*(1./dy1+1./(dy1+dy2)))
c1st          dazdy = (az(i,j,k)-az(i,jm,k))/dy1
        else
          dy = grid_params%dyh(jg)
          dazdy = (az(i,jp,k)-az(i,jm,k))/2./dy
        endif

        if (k==0) then
          dz1=grid_params%dz(kg)
          dz2=grid_params%dz(kg+1)
          daydz = (-ay(i,j,k+2)*dz1/dz2/(dz1+dz2)
     .             +ay(i,j,k+1)*(dz1+dz2)/dz1/dz2
     .             -ay(i,j,k  )*(1./dz1+1./(dz1+dz2)))
c1st          daydz = (ay(i,j,kp)-ay(i,j,k))/dz1
        elseif (k==nz+1) then
          dz1=grid_params%dz(kg-1)
          dz2=grid_params%dz(kg-2)
          daydz = -(-ay(i,j,k-2)*dz1/dz2/(dz1+dz2)
     .              +ay(i,j,k-1)*(dz1+dz2)/dz1/dz2
     .              -ay(i,j,k  )*(1./dz1+1./(dz1+dz2)))
c1st          daydz = (ay(i,j,k)-ay(i,j,km))/dz1
        else
          dz = grid_params%dzh(kg)
          daydz = (ay(i,j,kp)-ay(i,j,km))/2./dz
        endif

        curl2 = dazdy - daydz

      case(2)

        if (i==0) then
          dx1=grid_params%dx(ig)
          dx2=grid_params%dx(ig+1)
          dazdx = (-az(i+2,j,k)*dx1/dx2/(dx1+dx2)
     .             +az(i+1,j,k)*(dx1+dx2)/dx1/dx2
     .             -az(i  ,j,k)*(1./dx1+1./(dx1+dx2)))
c1st          dazdx = (az(ip,j,k)-az(i,j,k))/dx1
        elseif (i==nx+1) then
          dx1=grid_params%dx(ig-1)
          dx2=grid_params%dx(ig-2)
          dazdx = -(-az(i-2,j,k)*dx1/dx2/(dx1+dx2)
     .              +az(i-1,j,k)*(dx1+dx2)/dx1/dx2
     .              -az(i  ,j,k)*(1./dx1+1./(dx1+dx2)))
c1st          dazdx = (az(i,j,k)-az(im,j,k))/dx1
        elseif (i == 1 .and. bcond(1) == SP) then
          dx = grid_params%dxh(ig)
          dazdx = ((az(ip,j,k)+az(i,j,k))/2.-az(im,j,k))/dx
        else
          dx = grid_params%dxh(ig)
          dazdx = (az(ip,j,k)-az(im,j,k))/2./dx
        endif

        if (k==0) then
          dz1=grid_params%dz(kg)
          dz2=grid_params%dz(kg+1)
          daxdz = (-ax(i,j,k+2)*dz1/dz2/(dz1+dz2)
     .             +ax(i,j,k+1)*(dz1+dz2)/dz1/dz2
     .             -ax(i,j,k  )*(1./dz1+1./(dz1+dz2)))
c1st          daxdz = (ax(i,j,kp)-ax(i,j,k))/dz1
        elseif (k==nz+1) then
          dz1=grid_params%dz(kg-1)
          dz2=grid_params%dz(kg-2)
          daxdz = -(-ax(i,j,k-2)*dz1/dz2/(dz1+dz2)
     .              +ax(i,j,k-1)*(dz1+dz2)/dz1/dz2
     .              -ax(i,j,k  )*(1./dz1+1./(dz1+dz2)))
c1st          daxdz = (ax(i,j,k)-ax(i,j,km))/dz1
        else
          dz = grid_params%dzh(kg)
          daxdz = (ax(i,j,kp)-ax(i,j,km))/2./dz
        endif

        curl2 = daxdz - dazdx

      case(3)

        if (i==0) then
          dx1=grid_params%dx(ig)
          dx2=grid_params%dx(ig+1)
          daydx = (-ay(i+2,j,k)*dx1/dx2/(dx1+dx2)
     .             +ay(i+1,j,k)*(dx1+dx2)/dx1/dx2
     .             -ay(i  ,j,k)*(1./dx1+1./(dx1+dx2)))
c1st          daydx = (ay(ip,j,k)-ay(i,j,k))/dx1
        elseif (i==nx+1) then
          dx1=grid_params%dx(ig-1)
          dx2=grid_params%dx(ig-2)
          daydx =-(-ay(i-2,j,k)*dx1/dx2/(dx1+dx2)
     .             +ay(i-1,j,k)*(dx1+dx2)/dx1/dx2
     .             -ay(i  ,j,k)*(1./dx1+1/(dx1+dx2)))
c1st          daydx = (ay(i,j,k)-ay(im,j,k))/dx1
        elseif (i == 1 .and. bcond(1) == SP) then
          dx = grid_params%dxh(ig)
          daydx = ((ay(ip,j,k)/jacp+ay(i,j,k)/jac0)/2.*jacm
     .            - ay(im,j,k))/dx
        else
          dx = grid_params%dxh(ig)
          daydx = (ay(ip,j,k)-ay(im,j,k))/2./dx
        endif

        if (j==0) then
          dy1=grid_params%dy(jg)
          dy2=grid_params%dy(jg+1)
          daxdy = (-ax(i,j+2,k)*dy1/dy2/(dy1+dy2)
     .             +ax(i,j+1,k)*(dy1+dy2)/dy1/dy2
     .             -ax(i,j  ,k)*(1./dy1+1./(dy1+dy2)))
c1st          daxdy = (ax(i,jp,k)-ax(i,j,k))/dy1
        elseif (j==ny+1) then
          dy1=grid_params%dy(jg-1)
          dy2=grid_params%dy(jg-2)
          daxdy = -(-ax(i,j-2,k)*dy1/dy2/(dy1+dy2)
     .              +ax(i,j-1,k)*(dy1+dy2)/dy1/dy2
     .              -ax(i,j  ,k)*(1./dy1+1./(dy1+dy2)))
c1st          daxdy = (ax(i,j,k)-ax(i,jm,k))/dy1
        else
          dy = grid_params%dyh(jg)
          daxdy = (ax(i,jp,k)-ax(i,jm,k))/2./dy
        endif

        curl2 = daydx - daxdy

      case default

        write (*,*) 'Error in component in curl'
        write (*,*) 'Aborting...'
        stop

      end select

c     End program

      end function curl2

c     curlcurl
c     ###############################################################
      real*8 function curlcurl(i,j,k,ax,ay,az,diff,comp)

c     ---------------------------------------------------------------
c     Calculates curl(eta curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      use grid

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,comp

      real(8)    :: ax  (0:nx+1,0:ny+1,0:nz+1)
     .             ,ay  (0:nx+1,0:ny+1,0:nz+1)
     .             ,az  (0:nx+1,0:ny+1,0:nz+1)
     .             ,diff(0:nx+1,0:ny+1,0:nz+1)

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

        d_yy_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
     .           *g_sub_elem(2,2,x0,y0,zp,.false.)
        d_yy_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
     .           *g_sub_elem(2,2,x0,y0,zm,.false.)

        d_zz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
     .           *g_sub_elem(3,3,x0,yp,z0,.false.)
        d_zz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
     .           *g_sub_elem(3,3,x0,ym,z0,.false.)

        d_yz_jpkp = 4./(1./diff(i,jp,kp)+1./diff(i,jp,k )
     .                 +1./diff(i,j ,kp)+1./diff(i,j ,k ))
     .             *g_sub_elem(2,3,x0,yp,zp,.false.)
        d_yz_jpkm = 4./(1./diff(i,jp,km)+1./diff(i,jp,k )
     .                 +1./diff(i,j ,km)+1./diff(i,j ,k ))
     .             *g_sub_elem(2,3,x0,yp,zm,.false.)
        d_yz_jmkp = 4./(1./diff(i,jm,kp)+1./diff(i,jm,k )
     .                 +1./diff(i,j ,kp)+1./diff(i,j ,k ))
     .             *g_sub_elem(2,3,x0,ym,zp,.false.)
        d_yz_jmkm = 4./(1./diff(i,jm,km)+1./diff(i,jm,k )
     .                 +1./diff(i,j ,km)+1./diff(i,j ,k ))
     .             *g_sub_elem(2,3,x0,ym,zm,.false.)

        d_xy_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
     .           *g_sub_elem(1,2,x0,y0,zp,.false.)
        d_xy_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
     .           *g_sub_elem(1,2,x0,y0,zm,.false.)

        d_yz_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
     .           *g_sub_elem(2,3,x0,y0,zp,.false.)
        d_yz_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
     .           *g_sub_elem(2,3,x0,y0,zm,.false.)

        d_xz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,3,x0,yp,z0,.false.)
        d_xz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,3,x0,ym,z0,.false.)

        d_yz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
     .           *g_sub_elem(2,3,x0,yp,z0,.false.)
        d_yz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
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

        d_xx_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
     .           *g_sub_elem(1,1,x0,y0,zp,.false.)
        d_xx_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
     .           *g_sub_elem(1,1,x0,y0,zm,.false.)

        d_zz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(3,3,xp,y0,z0,.false.)
        d_zz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(3,3,xm,y0,z0,.false.)

        d_xz_ipkp = 4./(1./diff(ip,j,kp)+1./diff(ip,j,k )
     .                 +1./diff(i ,j,kp)+1./diff(i ,j,k ))
     .             *g_sub_elem(1,3,xp,y0,zp,.false.)
        d_xz_ipkm = 4./(1./diff(ip,j,km)+1./diff(ip,j,k )
     .                 +1./diff(i ,j,km)+1./diff(i ,j,k ))
     .             *g_sub_elem(1,3,xp,y0,zm,.false.)
        d_xz_imkp = 4./(1./diff(im,j,kp)+1./diff(im,j,k )
     .                 +1./diff(i ,j,kp)+1./diff(i ,j,k ))
     .             *g_sub_elem(1,3,xm,y0,zp,.false.)
        d_xz_imkm = 4./(1./diff(im,j,km)+1./diff(im,j,k )
     .                 +1./diff(i ,j,km)+1./diff(i ,j,k ))
     .             *g_sub_elem(1,3,xm,y0,zm,.false.)

        d_xy_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
     .           *g_sub_elem(1,2,x0,y0,zp,.false.)
        d_xy_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
     .           *g_sub_elem(1,2,x0,y0,zm,.false.)

        d_yz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(2,3,xp,y0,z0,.false.)
        d_yz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(2,3,xm,y0,z0,.false.)

        d_xz_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
     .           *g_sub_elem(1,3,x0,y0,zp,.false.)
        d_xz_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
     .           *g_sub_elem(1,3,x0,y0,zm,.false.)

        d_xz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,3,xp,y0,z0,.false.)
        d_xz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
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

        d_xx_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,1,x0,yp,z0,.false.)
        d_xx_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,1,x0,ym,z0,.false.)

        d_yy_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(2,2,xp,y0,z0,.false.)
        d_yy_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(2,2,xm,y0,z0,.false.)

        d_xy_ipjp = 4./(1./diff(ip,jp,k)+1./diff(ip,j,k )
     .                 +1./diff(i ,jp,k)+1./diff(i ,j,k ))
     .             *g_sub_elem(1,2,xp,yp,z0,.false.)
        d_xy_ipjm = 4./(1./diff(ip,jm,k)+1./diff(ip,j,k )
     .                 +1./diff(i ,jm,k)+1./diff(i ,j,k ))
     .             *g_sub_elem(1,2,xp,ym,z0,.false.)
        d_xy_imjp = 4./(1./diff(im,jp,k)+1./diff(im,j,k )
     .                 +1./diff(i ,jp,k)+1./diff(i ,j,k ))
     .             *g_sub_elem(1,2,xm,yp,z0,.false.)
        d_xy_imjm = 4./(1./diff(im,jm,k)+1./diff(im,j,k )
     .                 +1./diff(i ,jm,k)+1./diff(i ,j,k ))
     .             *g_sub_elem(1,2,xm,ym,z0,.false.)

        d_xy_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,2,x0,yp,z0,.false.)
        d_xy_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,2,x0,ym,z0,.false.)

        d_yz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(2,3,xp,y0,z0,.false.)
        d_yz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(2,3,xm,y0,z0,.false.)

        d_xy_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,2,xp,y0,z0,.false.)
        d_xy_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,2,xm,y0,z0,.false.)

        d_xz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
     .           *g_sub_elem(1,3,x0,yp,z0,.false.)
        d_xz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
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

c     End program

      end function curlcurl

c     fnabla_v
c     #############################################################
      function fnabla_v(i,j,k,x,y,z,ax,ay,az,cartesian,half_elem)
     .         result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor nabla(vec v) at the following positions:
c       + half_elem /= 1,2,3 => i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,half_elem
        real(8)    :: tensor(3,3),x,y,z
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)
        logical    :: cartesian

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: dhx,dhy,dhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: hessian1(3,3),hessian2(3,3),hessian3(3,3)
        logical    :: sing_point

c     Begin program

        sing_point = .false.
        if (i == 1 .and. bcond(1) == SP) sing_point = .true.

        !Defaults

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        dhx = 2.*dxh(ig)
        dhy = 2.*dyh(jg)
        dhz = 2.*dzh(kg)

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

          vxx = 0.5*(ax(i,j,k)+ax(ip,j,k))
          vyy = 0.5*(ay(i,j,k)+ay(ip,j,k))
          vzz = 0.5*(az(i,j,k)+az(ip,j,k))

          vxip = ax(ip,j,k)
          vxim = ax(i ,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(i ,j,k)
          vzip = az(ip,j,k)
          vzim = az(i ,j,k)

          vxjp = (ax(ip,jp,k)+ax(i,jp,k))/2.
          vxjm = (ax(ip,jm,k)+ax(i,jm,k))/2.
          vyjp = (ay(ip,jp,k)+ay(i,jp,k))/2.
          vyjm = (ay(ip,jm,k)+ay(i,jm,k))/2.
          vzjp = (az(ip,jp,k)+az(i,jp,k))/2.
          vzjm = (az(ip,jm,k)+az(i,jm,k))/2.

          vxkp = (ax(ip,j,kp)+ax(i,j,kp))/2.
          vxkm = (ax(ip,j,km)+ax(i,j,km))/2.
          vykp = (ay(ip,j,kp)+ay(i,j,kp))/2.
          vykm = (ay(ip,j,km)+ay(i,j,km))/2.
          vzkp = (az(ip,j,kp)+az(i,j,kp))/2.
          vzkm = (az(ip,j,km)+az(i,j,km))/2.
        case (2)
          dhy = dy(jg)

          vxx = 0.5*(ax(i,j,k)+ax(i,jp,k))
          vyy = 0.5*(ay(i,j,k)+ay(i,jp,k))
          vzz = 0.5*(az(i,j,k)+az(i,jp,k))

          if (sing_point) then
            vxip = (ax(ip,j,k)+ax(ip,jp,k))/2.
     .            +(ax(i ,j,k)+ax(i ,jp,k))/2.
            vxim = (ax(im,j,k)+ax(im,jp,k))
            vyip = (ay(ip,j,k)+ay(ip,jp,k))/2.
     .            +(ay(i ,j,k)+ay(i ,jp,k))/2.
            vyim = (ay(im,j,k)+ay(im,jp,k))
            vzip = (az(ip,j,k)+az(ip,jp,k))/2
     .            +(az(i ,j,k)+az(i ,jp,k))/2.
            vzim = (az(im,j,k)+az(im,jp,k))
          else
            vxip = (ax(ip,j,k)+ax(ip,jp,k))/2.
            vxim = (ax(im,j,k)+ax(im,jp,k))/2.
            vyip = (ay(ip,j,k)+ay(ip,jp,k))/2.
            vyim = (ay(im,j,k)+ay(im,jp,k))/2.
            vzip = (az(ip,j,k)+az(ip,jp,k))/2.
            vzim = (az(im,j,k)+az(im,jp,k))/2.
          endif

	  vxjp = ax(i,jp,k)
	  vxjm = ax(i,j ,k)
	  vyjp = ay(i,jp,k)
	  vyjm = ay(i,j ,k)
	  vzjp = az(i,jp,k)
	  vzjm = az(i,j ,k)

          vxkp = (ax(i,j,kp)+ax(i,jp,kp))/2.
          vxkm = (ax(i,j,km)+ax(i,jp,km))/2.
          vykp = (ay(i,j,kp)+ay(i,jp,kp))/2.
          vykm = (ay(i,j,km)+ay(i,jp,km))/2.
          vzkp = (az(i,j,kp)+az(i,jp,kp))/2.
          vzkm = (az(i,j,km)+az(i,jp,km))/2.
        case (3)
          dhz = dz(kg)

          vxx = 0.5*(ax(i,j,k)+ax(i,j,kp))
          vyy = 0.5*(ay(i,j,k)+ay(i,j,kp))
          vzz = 0.5*(az(i,j,k)+az(i,j,kp))

          if (sing_point) then
            vxip = (ax(ip,j,k)+ax(ip,j,kp))/2.
     .            +(ax(i ,j,k)+ax(i ,j,kp))/2.
            vxim = (ax(im,j,k)+ax(im,j,kp))
            vyip = (ay(ip,j,k)+ay(ip,j,kp))/2.
     .            +(ay(i ,j,k)+ay(i ,j,kp))/2.
            vyim = (ay(im,j,k)+ay(im,j,kp))
            vzip = (az(ip,j,k)+az(ip,j,kp))/2.
     .            +(az(i ,j,k)+az(i ,j,kp))/2.
            vzim = (az(im,j,k)+az(im,j,kp))
          else
            vxip = (ax(ip,j,k)+ax(ip,j,kp))/2.
            vxim = (ax(im,j,k)+ax(im,j,kp))/2.
            vyip = (ay(ip,j,k)+ay(ip,j,kp))/2.
            vyim = (ay(im,j,k)+ay(im,j,kp))/2.
            vzip = (az(ip,j,k)+az(ip,j,kp))/2.
            vzim = (az(im,j,k)+az(im,j,kp))/2.
          endif

          vxjp = (ax(i,jp,k)+ax(i,jp,kp))/2.
          vxjm = (ax(i,jm,k)+ax(i,jm,kp))/2.
          vyjp = (ay(i,jp,k)+ay(i,jp,kp))/2.
          vyjm = (ay(i,jm,k)+ay(i,jm,kp))/2.
          vzjp = (az(i,jp,k)+az(i,jp,kp))/2.
          vzjm = (az(i,jm,k)+az(i,jm,kp))/2.

	  vxkp = ax(i,j,kp)
	  vxkm = ax(i,j,k )
	  vykp = ay(i,j,kp)
	  vykm = ay(i,j,k )
	  vzkp = az(i,j,kp)
	  vzkm = az(i,j,k )

        case default

          vxx = ax(i,j,k)
          vyy = ay(i,j,k)
          vzz = az(i,j,k)

          if (sing_point) then
            vxip = ax(ip,j,k)+ax(i,j,k)
            vxim = 2.*ax(i ,j,k)
            vyip = ay(ip,j,k)+ay(i,j,k)
            vyim = 2.*ay(im,j,k)
            vzip = az(ip,j,k)+az(i,j,k)
            vzim = 2.*az(im,j,k)
          else
            vxip = ax(ip,j,k)
            vxim = ax(im,j,k)
            vyip = ay(ip,j,k)
            vyim = ay(im,j,k)
            vzip = az(ip,j,k)
            vzim = az(im,j,k)
          endif

	  vxjp = ax(i,jp,k)
	  vxjm = ax(i,jm,k)
	  vyjp = ay(i,jp,k)
	  vyjm = ay(i,jm,k)
	  vzjp = az(i,jp,k)
	  vzjm = az(i,jm,k)

	  vxkp = ax(i,j,kp)
	  vxkm = ax(i,j,km)
	  vykp = ay(i,j,kp)
	  vykm = ay(i,j,km)
	  vzkp = az(i,j,kp)
	  vzkm = az(i,j,km)

        end select

c     !Calculate nabla_v tensor

        hessian1 = hessian(1,x,y,z,cartesian)
        hessian2 = hessian(2,x,y,z,cartesian)
        hessian3 = hessian(3,x,y,z,cartesian)

      ! l = 1, m = 1
        tensor(1,1) = (vxip-vxim)/dhx
     .               + vxx*(hessian1(1,1)
     .                    + hessian2(2,1)
     .                    + hessian3(3,1))
     .               - vxx*hessian1(1,1)
     .               - vyy*hessian1(2,1)
     .               - vzz*hessian1(3,1)

      ! l = 1, m = 2
        tensor(1,2) = (vyip-vyim)/dhx
     .               + vyy*(hessian1(1,1)
     .                    + hessian2(2,1)
     .                    + hessian3(3,1))
     .               - vxx*hessian2(1,1)
     .               - vyy*hessian2(2,1)
     .               - vzz*hessian2(3,1)

      ! l = 1, m = 3
        tensor(1,3) = (vzip-vzim)/dhx
     .               + vzz*(hessian1(1,1)
     .                    + hessian2(2,1)
     .                    + hessian3(3,1))
     .               - vxx*hessian3(1,1)
     .               - vyy*hessian3(2,1)
     .               - vzz*hessian3(3,1)

      ! l = 2, m = 1
        tensor(2,1) = (vxjp-vxjm)/dhy
     .               + vxx*(hessian1(1,2)
     .                    + hessian2(2,2)
     .                    + hessian3(3,2))
     .               - vxx*hessian1(1,2)
     .               - vyy*hessian1(2,2)
     .               - vzz*hessian1(3,2)

      ! l = 2, m = 2
        tensor(2,2) = (vyjp-vyjm)/dhy
     .               + vyy*(hessian1(1,2)
     .                    + hessian2(2,2)
     .                    + hessian3(3,2))
     .               - vxx*hessian2(1,2)
     .               - vyy*hessian2(2,2)
     .               - vzz*hessian2(3,2)

      ! l = 2, m = 3
        tensor(2,3) = (vzjp-vzjm)/dhy
     .               + vzz*(hessian1(1,2)
     .                    + hessian2(2,2)
     .                    + hessian3(3,2))
     .               - vxx*hessian3(1,2)
     .               - vyy*hessian3(2,2)
     .               - vzz*hessian3(3,2)

      ! l = 3, m = 1
        tensor(3,1) = (vxkp-vxkm)/dhz
     .               + vxx*(hessian1(1,3)
     .                    + hessian2(2,3)
     .                    + hessian3(3,3))
     .               - vxx*hessian1(1,3)
     .               - vyy*hessian1(2,3)
     .               - vzz*hessian1(3,3)

      ! l = 3, m = 2
        tensor(3,2) = (vykp-vykm)/dhz
     .               + vyy*(hessian1(1,3)
     .                    + hessian2(2,3)
     .                    + hessian3(3,3))
     .               - vxx*hessian2(1,3)
     .               - vyy*hessian2(2,3)
     .               - vzz*hessian2(3,3)

      ! l = 3, m = 3
        tensor(3,3) = (vzkp-vzkm)/dhz
     .               + vzz*(hessian1(1,3)
     .                    + hessian2(2,3)
     .                    + hessian3(3,3))
     .               - vxx*hessian3(1,3)
     .               - vyy*hessian3(2,3)
     .               - vzz*hessian3(3,3)

c     End program

      end function fnabla_v

      end module operators

c module nlfunction_setup
c ######################################################################
      module nlfunction_setup

        use grid

        use grid_aliases

        use parameters

        use transport_params

        use equilibrium

        use auxiliaryVariables

        use operators

        real(8)    :: xim,yim,zim,xip,yip,zip
     .               ,xjm,yjm,zjm,xjp,yjp,zjp
     .               ,xkm,ykm,zkm,xkp,ykp,zkp

        real(8)    :: x0,y0,z0,xh,yh,zh

        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .               ,jacp,jacm,jach,jac0

        real(8)    :: gsub(3,3),gsuper(3,3)
     .               ,cnv1(3),cnv2(3),cnv3(3),jac

        real(8)    :: nabla_v(3,3),hessian1(3,3)
     .               ,hessian2(3,3),hessian3(3,3)
     .               ,cov_tnsr(3,3),cnv_tnsr(3,3)

        logical    :: sing_point,cartesian,alt_eom

      contains

ccc     divB
ccc     ###############################################################
cc      real(8) function divB(i,j,k)
cc      implicit none
ccc     ---------------------------------------------------------------
ccc     Calculates divergence of magnetic field at grid vertices.
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer(4) :: i,j,k
cc
ccc     Local variables
cc
cc      integer(4) :: ig,jg,kg
cc      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacp,jac0,xh,yh,zh,jach
cc     .             ,xip,yip,zip
cc      logical    :: cartesian
cc
ccc     Begin program
cc
cc      call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x0,y0,z0,cartesian)
cc
cc      dxx = dxh(ig)
cc      dyy = dyh(jg)
cc      dzz = dzh(kg)
cc
cc      if (i == 1 .and. bcond(1) == SP) then
cc        call getCoordinates(i+1,j,k,igx,igy,igz,ig,jg,kg,xip,yip,zip
cc     .                     ,cartesian)
cc
cc        jacp = jacobian(xip,yip,zip,cartesian)
cc        jac0 = jacobian(x0 ,y0 ,z0 ,cartesian)
cc
cc        xh = (xip+x0)/2.
cc        yh = (yip+y0)/2.
cc        zh = (zip+z0)/2.
cc        jach = jacobian(xh ,yh ,zh ,cartesian)
cccc        jach = 0.5*(jacp+jac0)   !Only good for cylindrical-like geom.
cc
cc        divB =  ((bx(i+1,j  ,k  )/jacp
cc     .          +bx(i  ,j  ,k  )/jac0)*jach)    /2./dxx
cc     .        +(by(i  ,j+1,k  )-by(i  ,j-1,k  ))/2./dyy
cc     .        +(bz(i  ,j  ,k+1)-bz(i  ,j  ,k-1))/2./dzz
cc      else
cc        divB =  (bx(i+1,j  ,k  )-bx(i-1,j  ,k  ))/2./dxx
cc     .        +(by(i  ,j+1,k  )-by(i  ,j-1,k  ))/2./dyy
cc     .        +(bz(i  ,j  ,k+1)-bz(i  ,j  ,k-1))/2./dzz
cc      endif
cc
ccc     End 
cc
cc      end function divB
cc
ccc     divJ
ccc     ###############################################################
cc      real(8) function divJ(i,j,k)
cc      implicit none
ccc     ---------------------------------------------------------------
ccc     Calculates divergence of magnetic field at grid vertices.
ccc     ---------------------------------------------------------------
cc
ccc     Call variables
cc
cc      integer(4) :: i,j,k
cc
ccc     Local variables
cc
cc      integer(4) :: ig,jg,kg
cc      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacp,jac0,xh,yh,zh,jach
cc     .             ,xip,yip,zip
cc      logical    :: cartesian
cc
ccc     Begin program
cc
cc      call getCoordinates(i,j,k,igx,igy,igz,ig,jg,kg,x0,y0,z0,cartesian)
cc
cc      dxx = dxh(ig)
cc      dyy = dyh(jg)
cc      dzz = dzh(kg)
cc
cc      if (i == 1 .and. bcond(1) == SP) then
cc        call getCoordinates(i+1,j,k,igx,igy,igz,ig,jg,kg,xip,yip,zip
cc     .                     ,cartesian)
cc
cc        jacp = jacobian(xip,yip,zip,cartesian)
cc        jac0 = jacobian(x0 ,y0 ,z0 ,cartesian)
cc
cc        xh = (xip+x0)/2.
cc        yh = (yip+y0)/2.
cc        zh = (zip+z0)/2.
cc        jach = jacobian(xh ,yh ,zh ,cartesian)
cccc        jach = 0.5*(jacp+jac0)   !Only good for cylindrical-like geom.
cc
cc        divJ =  ((jx(i+1,j  ,k  )/jacp
cc     .          +jx(i  ,j  ,k  )/jac0)*jach)    /2./dxx
cc     .        +(jy(i  ,j+1,k  )-jy(i  ,j-1,k  ))/2./dyy
cc     .        +(jz(i  ,j  ,k+1)-jz(i  ,j  ,k-1))/2./dzz
cc      else
cc        divJ =  (jx(i+1,j  ,k  )-jx(i-1,j  ,k  ))/2./dxx
cc     .        +(jy(i  ,j+1,k  )-jy(i  ,j-1,k  ))/2./dyy
cc     .        +(jz(i  ,j  ,k+1)-jz(i  ,j  ,k-1))/2./dzz
cc      endif
cc
ccc     End 
cc
cc      end function divJ

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

ccc     fnabla_v
ccc     #############################################################
cc      function fnabla_v(i,j,k,x,y,z,half_elem) result(tensor)
ccc     -------------------------------------------------------------
ccc     Calculates the tensor nabla(vec v) at the following positions:
ccc       + half_elem /= 1,2,3 => i,j,k
ccc       + half_elem=1 --> i+1/2,j,k
ccc       + half_elem=2 --> i,j+1/2,k
ccc       + half_elem=3 --> i,j,k+1/2
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,half_elem
cc        real(8)    :: tensor(3,3),x,y,z
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
cc        real(8)    :: jx,jy,jz
cc        real(8)    :: dhx,dhy,dhz,vxx,vyy,vzz
cc
ccc     Begin program
cc
cc        !Defaults
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        dhx = 2.*dxh(ig)
cc        dhy = 2.*dyh(jg)
cc        dhz = 2.*dzh(kg)
cc
cc        vxx = vx(i,j,k)
cc        vyy = vy(i,j,k)
cc        vzz = vz(i,j,k)
cc
cc        ip = i+1
cc        im = i-1
cc        jp = j+1
cc        jm = j-1
cc        kp = k+1
cc        km = k-1
cc
cc        !Exceptions
cc
cc        select case(half_elem)
cc        case (1)
cc          dhx = dx(ig)
cc          im  = i
cc          vxx = 0.5*(vx(i,j,k)+vx(ip,j,k))
cc          vyy = 0.5*(vy(i,j,k)+vy(ip,j,k))
cc          vzz = 0.5*(vz(i,j,k)+vz(ip,j,k))
cc        case (2)
cc          dhy = dy(jg)
cc          jm  = j
cc          vxx = 0.5*(vx(i,j,k)+vx(i,jp,k))
cc          vyy = 0.5*(vy(i,j,k)+vy(i,jp,k))
cc          vzz = 0.5*(vz(i,j,k)+vz(i,jp,k))
cc        case (3)
cc          dhz = dz(kg)
cc          km  = k
cc          vxx = 0.5*(vx(i,j,k)+vx(i,j,kp))
cc          vyy = 0.5*(vy(i,j,k)+vy(i,j,kp))
cc          vzz = 0.5*(vz(i,j,k)+vz(i,j,kp))
cc        end select
cc
ccc     !Calculate nabla_v tensor
cc
cc        hessian1 = hessian(1,x,y,z,cartesian)
cc        hessian2 = hessian(2,x,y,z,cartesian)
cc        hessian3 = hessian(3,x,y,z,cartesian)
cc
cc      ! l = 1, m = 1
cc        if (sing_point .and. half_elem /= 1) then
cc          tensor(1,1) = (vx(ip,j,k)+vx(i,j,k)-2.*vx(im,j,k))/dhx
cc     .                 + vxx*(hessian1(1,1)
cc     .                      + hessian2(2,1)
cc     .                      + hessian3(3,1))
cc     .                 - vxx*hessian1(1,1)
cc     .                 - vyy*hessian1(2,1)
cc     .                 - vzz*hessian1(3,1)
cc        else
cc          tensor(1,1) = (vx(ip,j,k)-vx(im,j,k))/dhx
cc     .                 + vxx*(hessian1(1,1)
cc     .                      + hessian2(2,1)
cc     .                      + hessian3(3,1))
cc     .                 - vxx*hessian1(1,1)
cc     .                 - vyy*hessian1(2,1)
cc     .                 - vzz*hessian1(3,1)
cc        endif
cc
cc      ! l = 1, m = 2
cc        if (sing_point .and. half_elem /= 1) then
cc          tensor(1,2) = (vy(ip,j,k)+vy(i,j,k)-2.*vy(im,j,k))/dhx
cc     .                 + vyy*(hessian1(1,1)
cc     .                      + hessian2(2,1)
cc     .                      + hessian3(3,1))
cc     .                 - vxx*hessian2(1,1)
cc     .                 - vyy*hessian2(2,1)
cc     .                 - vzz*hessian2(3,1)
cc
cc        else
cc          tensor(1,2) = (vy(ip,j,k)-vy(im,j,k))/dhx
cc     .                 + vyy*(hessian1(1,1)
cc     .                      + hessian2(2,1)
cc     .                      + hessian3(3,1))
cc     .                 - vxx*hessian2(1,1)
cc     .                 - vyy*hessian2(2,1)
cc     .                 - vzz*hessian2(3,1)
cc        endif
cc
cc      ! l = 1, m = 3
cc        if (sing_point .and. half_elem /= 1) then
cc          tensor(1,3) = (vz(ip,j,k)+vz(i,j,k)-2.*vz(im,j,k))/dhx
cc     .                 + vzz*(hessian1(1,1)
cc     .                      + hessian2(2,1)
cc     .                      + hessian3(3,1))
cc     .                 - vxx*hessian3(1,1)
cc     .                 - vyy*hessian3(2,1)
cc     .                 - vzz*hessian3(3,1)
cc        else
cc          tensor(1,3) = (vz(ip,j,k)-vz(im,j,k))/dhx
cc     .                 + vzz*(hessian1(1,1)
cc     .                      + hessian2(2,1)
cc     .                      + hessian3(3,1))
cc     .                 - vxx*hessian3(1,1)
cc     .                 - vyy*hessian3(2,1)
cc     .                 - vzz*hessian3(3,1)
cc        endif
cc
cc      ! l = 2, m = 1
cc        tensor(2,1) = (vx(i,jp,k)-vx(i,jm,k))/dhy
cc     .               + vxx*(hessian1(1,2)
cc     .                    + hessian2(2,2)
cc     .                    + hessian3(3,2))
cc     .               - vxx*hessian1(1,2)
cc     .               - vyy*hessian1(2,2)
cc     .               - vzz*hessian1(3,2)
cc
cc      ! l = 2, m = 2
cc        tensor(2,2) = (vy(i,jp,k)-vy(i,jm,k))/dhy
cc     .               + vyy*(hessian1(1,2)
cc     .                    + hessian2(2,2)
cc     .                    + hessian3(3,2))
cc     .               - vxx*hessian2(1,2)
cc     .               - vyy*hessian2(2,2)
cc     .               - vzz*hessian2(3,2)
cc
cc      ! l = 2, m = 3
cc        tensor(2,3) = (vz(i,jp,k)-vz(i,jm,k))/dhy
cc     .               + vzz*(hessian1(1,2)
cc     .                    + hessian2(2,2)
cc     .                    + hessian3(3,2))
cc     .               - vxx*hessian3(1,2)
cc     .               - vyy*hessian3(2,2)
cc     .               - vzz*hessian3(3,2)
cc
cc      ! l = 3, m = 1
cc        tensor(3,1) = (vx(i,j,kp)-vx(i,j,km))/dhz
cc     .               + vxx*(hessian1(1,3)
cc     .                    + hessian2(2,3)
cc     .                    + hessian3(3,3))
cc     .               - vxx*hessian1(1,3)
cc     .               - vyy*hessian1(2,3)
cc     .               - vzz*hessian1(3,3)
cc
cc      ! l = 3, m = 2
cc        tensor(3,2) = (vy(i,j,kp)-vy(i,j,km))/dhz
cc     .               + vyy*(hessian1(1,3)
cc     .                    + hessian2(2,3)
cc     .                    + hessian3(3,3))
cc     .               - vxx*hessian2(1,3)
cc     .               - vyy*hessian2(2,3)
cc     .               - vzz*hessian2(3,3)
cc
cc      ! l = 3, m = 3
cc        tensor(3,3) = (vz(i,j,kp)-vz(i,j,km))/dhz
cc     .               + vzz*(hessian1(1,3)
cc     .                    + hessian2(2,3)
cc     .                    + hessian3(3,3))
cc     .               - vxx*hessian3(1,3)
cc     .               - vyy*hessian3(2,3)
cc     .               - vzz*hessian3(3,3)
cc
ccc     End program
cc
cc      end function fnabla_v

ccc     fnabla_v
ccc     #############################################################
cc      function fnabla_v(i,j,k,x,y,z,half_elem) result(tensor)
ccc     -------------------------------------------------------------
ccc     Calculates the tensor nabla(vec v) at the following positions:
ccc       + half_elem /= 1,2,3 => i,j,k
ccc       + half_elem=1 --> i+1/2,j,k
ccc       + half_elem=2 --> i,j+1/2,k
ccc       + half_elem=3 --> i,j,k+1/2
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,half_elem
cc        real(8)    :: tensor(3,3),x,y,z
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
cc        real(8)    :: dhx,dhy,dhz
cc        real(8)    :: vxx,vyy,vzz
cc     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
cc     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
cc     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
cc
ccc     Begin program
cc
cc        !Defaults
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        dhx = 2.*dxh(ig)
cc        dhy = 2.*dyh(jg)
cc        dhz = 2.*dzh(kg)
cc
cc        ip = i+1
cc        im = i-1
cc        jp = j+1
cc        jm = j-1
cc        kp = k+1
cc        km = k-1
cc
cc        !Exceptions
cc
cc        select case(half_elem)
cc        case (1)
cc          dhx = dx(ig)
cc
cc          vxx = 0.5*(vx(i,j,k)+vx(ip,j,k))
cc          vyy = 0.5*(vy(i,j,k)+vy(ip,j,k))
cc          vzz = 0.5*(vz(i,j,k)+vz(ip,j,k))
cc
cc          vxip = vx(ip,j,k)
cc          vxim = vx(i ,j,k)
cc          vyip = vy(ip,j,k)
cc          vyim = vy(i ,j,k)
cc          vzip = vz(ip,j,k)
cc          vzim = vz(i ,j,k)
cc
cc          vxjp = (vx(ip,jp,k)+vx(i,jp,k))/2.
cc          vxjm = (vx(ip,jm,k)+vx(i,jm,k))/2.
cc          vyjp = (vy(ip,jp,k)+vy(i,jp,k))/2.
cc          vyjm = (vy(ip,jm,k)+vy(i,jm,k))/2.
cc          vzjp = (vz(ip,jp,k)+vz(i,jp,k))/2.
cc          vzjm = (vz(ip,jm,k)+vz(i,jm,k))/2.
cc
cc          vxkp = (vx(ip,j,kp)+vx(i,j,kp))/2.
cc          vxkm = (vx(ip,j,km)+vx(i,j,km))/2.
cc          vykp = (vy(ip,j,kp)+vy(i,j,kp))/2.
cc          vykm = (vy(ip,j,km)+vy(i,j,km))/2.
cc          vzkp = (vz(ip,j,kp)+vz(i,j,kp))/2.
cc          vzkm = (vz(ip,j,km)+vz(i,j,km))/2.
cc        case (2)
cc          dhy = dy(jg)
cc
cc          vxx = 0.5*(vx(i,j,k)+vx(i,jp,k))
cc          vyy = 0.5*(vy(i,j,k)+vy(i,jp,k))
cc          vzz = 0.5*(vz(i,j,k)+vz(i,jp,k))
cc
cc          if (sing_point) then
cc            vxip = (vx(ip,j,k)+vx(ip,jp,k))/2.
cc     .            +(vx(i ,j,k)+vx(i ,jp,k))/2.
cc            vxim = (vx(im,j,k)+vx(im,jp,k))
cc            vyip = (vy(ip,j,k)+vy(ip,jp,k))/2.
cc     .            +(vy(i ,j,k)+vy(i ,jp,k))/2.
cc            vyim = (vy(im,j,k)+vy(im,jp,k))
cc            vzip = (vz(ip,j,k)+vz(ip,jp,k))/2
cc     .            +(vz(i ,j,k)+vz(i ,jp,k))/2.
cc            vzim = (vz(im,j,k)+vz(im,jp,k))
cc          else
cc            vxip = (vx(ip,j,k)+vx(ip,jp,k))/2.
cc            vxim = (vx(im,j,k)+vx(im,jp,k))/2.
cc            vyip = (vy(ip,j,k)+vy(ip,jp,k))/2.
cc            vyim = (vy(im,j,k)+vy(im,jp,k))/2.
cc            vzip = (vz(ip,j,k)+vz(ip,jp,k))/2.
cc            vzim = (vz(im,j,k)+vz(im,jp,k))/2.
cc          endif
cc
cc	  vxjp = vx(i,jp,k)
cc	  vxjm = vx(i,j ,k)
cc	  vyjp = vy(i,jp,k)
cc	  vyjm = vy(i,j ,k)
cc	  vzjp = vz(i,jp,k)
cc	  vzjm = vz(i,j ,k)
cc
cc          vxkp = (vx(i,j,kp)+vx(i,jp,kp))/2.
cc          vxkm = (vx(i,j,km)+vx(i,jp,km))/2.
cc          vykp = (vy(i,j,kp)+vy(i,jp,kp))/2.
cc          vykm = (vy(i,j,km)+vy(i,jp,km))/2.
cc          vzkp = (vz(i,j,kp)+vz(i,jp,kp))/2.
cc          vzkm = (vz(i,j,km)+vz(i,jp,km))/2.
cc        case (3)
cc          dhz = dz(kg)
cc
cc          vxx = 0.5*(vx(i,j,k)+vx(i,j,kp))
cc          vyy = 0.5*(vy(i,j,k)+vy(i,j,kp))
cc          vzz = 0.5*(vz(i,j,k)+vz(i,j,kp))
cc
cc          if (sing_point) then
cc            vxip = (vx(ip,j,k)+vx(ip,j,kp))/2.
cc     .            +(vx(i ,j,k)+vx(i ,j,kp))/2.
cc            vxim = (vx(im,j,k)+vx(im,j,kp))
cc            vyip = (vy(ip,j,k)+vy(ip,j,kp))/2.
cc     .            +(vy(i ,j,k)+vy(i ,j,kp))/2.
cc            vyim = (vy(im,j,k)+vy(im,j,kp))
cc            vzip = (vz(ip,j,k)+vz(ip,j,kp))/2.
cc     .            +(vz(i ,j,k)+vz(i ,j,kp))/2.
cc            vzim = (vz(im,j,k)+vz(im,j,kp))
cc          else
cc            vxip = (vx(ip,j,k)+vx(ip,j,kp))/2.
cc            vxim = (vx(im,j,k)+vx(im,j,kp))/2.
cc            vyip = (vy(ip,j,k)+vy(ip,j,kp))/2.
cc            vyim = (vy(im,j,k)+vy(im,j,kp))/2.
cc            vzip = (vz(ip,j,k)+vz(ip,j,kp))/2.
cc            vzim = (vz(im,j,k)+vz(im,j,kp))/2.
cc          endif
cc
cc          vxjp = (vx(i,jp,k)+vx(i,jp,kp))/2.
cc          vxjm = (vx(i,jm,k)+vx(i,jm,kp))/2.
cc          vyjp = (vy(i,jp,k)+vy(i,jp,kp))/2.
cc          vyjm = (vy(i,jm,k)+vy(i,jm,kp))/2.
cc          vzjp = (vz(i,jp,k)+vz(i,jp,kp))/2.
cc          vzjm = (vz(i,jm,k)+vz(i,jm,kp))/2.
cc
cc	  vxkp = vx(i,j,kp)
cc	  vxkm = vx(i,j,k )
cc	  vykp = vy(i,j,kp)
cc	  vykm = vy(i,j,k )
cc	  vzkp = vz(i,j,kp)
cc	  vzkm = vz(i,j,k )
cc
cc        case default
cc
cc          vxx = vx(i,j,k)
cc          vyy = vy(i,j,k)
cc          vzz = vz(i,j,k)
cc
cc          if (sing_point) then
cc            vxip = vx(ip,j,k)+vx(i,j,k)
cc            vxim = 2.*vx(i ,j,k)
cc            vyip = vy(ip,j,k)+vy(i,j,k)
cc            vyim = 2.*vy(im,j,k)
cc            vzip = vz(ip,j,k)+vz(i,j,k)
cc            vzim = 2.*vz(im,j,k)
cc          else
cc            vxip = vx(ip,j,k)
cc            vxim = vx(im,j,k)
cc            vyip = vy(ip,j,k)
cc            vyim = vy(im,j,k)
cc            vzip = vz(ip,j,k)
cc            vzim = vz(im,j,k)
cc          endif
cc
cc	  vxjp = vx(i,jp,k)
cc	  vxjm = vx(i,jm,k)
cc	  vyjp = vy(i,jp,k)
cc	  vyjm = vy(i,jm,k)
cc	  vzjp = vz(i,jp,k)
cc	  vzjm = vz(i,jm,k)
cc
cc	  vxkp = vx(i,j,kp)
cc	  vxkm = vx(i,j,km)
cc	  vykp = vy(i,j,kp)
cc	  vykm = vy(i,j,km)
cc	  vzkp = vz(i,j,kp)
cc	  vzkm = vz(i,j,km)
cc
cc        end select
cc
ccc     !Calculate nabla_v tensor
cc
cc        hessian1 = hessian(1,x,y,z,cartesian)
cc        hessian2 = hessian(2,x,y,z,cartesian)
cc        hessian3 = hessian(3,x,y,z,cartesian)
cc
cc      ! l = 1, m = 1
cc        tensor(1,1) = (vxip-vxim)/dhx
cc     .               + vxx*(hessian1(1,1)
cc     .                    + hessian2(2,1)
cc     .                    + hessian3(3,1))
cc     .               - vxx*hessian1(1,1)
cc     .               - vyy*hessian1(2,1)
cc     .               - vzz*hessian1(3,1)
cc
cc      ! l = 1, m = 2
cc        tensor(1,2) = (vyip-vyim)/dhx
cc     .               + vyy*(hessian1(1,1)
cc     .                    + hessian2(2,1)
cc     .                    + hessian3(3,1))
cc     .               - vxx*hessian2(1,1)
cc     .               - vyy*hessian2(2,1)
cc     .               - vzz*hessian2(3,1)
cc
cc      ! l = 1, m = 3
cc        tensor(1,3) = (vzip-vzim)/dhx
cc     .               + vzz*(hessian1(1,1)
cc     .                    + hessian2(2,1)
cc     .                    + hessian3(3,1))
cc     .               - vxx*hessian3(1,1)
cc     .               - vyy*hessian3(2,1)
cc     .               - vzz*hessian3(3,1)
cc
cc      ! l = 2, m = 1
cc        tensor(2,1) = (vxjp-vxjm)/dhy
cc     .               + vxx*(hessian1(1,2)
cc     .                    + hessian2(2,2)
cc     .                    + hessian3(3,2))
cc     .               - vxx*hessian1(1,2)
cc     .               - vyy*hessian1(2,2)
cc     .               - vzz*hessian1(3,2)
cc
cc      ! l = 2, m = 2
cc        tensor(2,2) = (vyjp-vyjm)/dhy
cc     .               + vyy*(hessian1(1,2)
cc     .                    + hessian2(2,2)
cc     .                    + hessian3(3,2))
cc     .               - vxx*hessian2(1,2)
cc     .               - vyy*hessian2(2,2)
cc     .               - vzz*hessian2(3,2)
cc
cc      ! l = 2, m = 3
cc        tensor(2,3) = (vzjp-vzjm)/dhy
cc     .               + vzz*(hessian1(1,2)
cc     .                    + hessian2(2,2)
cc     .                    + hessian3(3,2))
cc     .               - vxx*hessian3(1,2)
cc     .               - vyy*hessian3(2,2)
cc     .               - vzz*hessian3(3,2)
cc
cc      ! l = 3, m = 1
cc        tensor(3,1) = (vxkp-vxkm)/dhz
cc     .               + vxx*(hessian1(1,3)
cc     .                    + hessian2(2,3)
cc     .                    + hessian3(3,3))
cc     .               - vxx*hessian1(1,3)
cc     .               - vyy*hessian1(2,3)
cc     .               - vzz*hessian1(3,3)
cc
cc      ! l = 3, m = 2
cc        tensor(3,2) = (vykp-vykm)/dhz
cc     .               + vyy*(hessian1(1,3)
cc     .                    + hessian2(2,3)
cc     .                    + hessian3(3,3))
cc     .               - vxx*hessian2(1,3)
cc     .               - vyy*hessian2(2,3)
cc     .               - vzz*hessian2(3,3)
cc
cc      ! l = 3, m = 3
cc        tensor(3,3) = (vzkp-vzkm)/dhz
cc     .               + vzz*(hessian1(1,3)
cc     .                    + hessian2(2,3)
cc     .                    + hessian3(3,3))
cc     .               - vxx*hessian3(1,3)
cc     .               - vyy*hessian3(2,3)
cc     .               - vzz*hessian3(3,3)
cc
ccc     End program
cc
cc      end function fnabla_v

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
        real(8)    :: jac,jac0,jacp,ptot,vis

c     Begin program

        ip = i+1

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
          if (.not.alt_eom .and. sing_point) then
            ip = i
            x = xim
            y = yim
            z = zim
          else
            x = (xim+x0)/2.
            y = (yim+y0)/2.
            z = (zim+z0)/2.
          endif
        endif

        jac    = jacobian(x,y,z,cartesian)

        gsuper = g_super (x,y,z,cartesian)

cc        if (i == 1 .and. bcond(1) == SP .and. flag == 1) then
cc          jacp = jacobian(xip,yip,zip,cartesian)
cc          jac0 = jacobian(x0 ,y0 ,z0 ,cartesian)
cc        elseif (i == 2 .and. bcond(1) == SP .and. flag == -1) then
cc          jacp = jacobian(x0 ,y0 ,z0 ,cartesian)
cc          jac0 = jacobian(xim,yim,zim,cartesian)
cc        else
cc          jacp = jac
cc          jac0 = jac
cc        endif
        if (i < nx .and. bcond(1) == SP .and. flag == 1) then
          jacp = jacobian(xip,yip,zip,cartesian)
          jac0 = jacobian(x0 ,y0 ,z0 ,cartesian)
        elseif (i < nx .and. bcond(1) == SP .and. flag == -1) then
          jacp = jacobian(x0 ,y0 ,z0 ,cartesian)
          jac0 = jacobian(xim,yim,zim,cartesian)
        else
          jacp = jac
          jac0 = jac
        endif

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,x,y,z,vx,vy,vz,cartesian,1)
        else
          nabla_v = fnabla_v(i,j,k,x,y,z,vx,vy,vz,cartesian,0)
        endif

        !Recall p=2nT
        ptot = jac*(rho(ip,j,k)*tmp(ip,j,k)
     .             +rho(i ,j,k)*tmp(i ,j,k))
     .       +jac*(bx(ip,j,k)*bx_cov(i ,j,k)/jacp
     .            +bx(i ,j,k)*bx_cov(ip,j,k)/jac0
     .            +by(ip,j,k)*by_cov(i ,j,k)/jac0
     .            +by(i ,j,k)*by_cov(ip,j,k)/jacp
     .            +bz(ip,j,k)*bz_cov(i ,j,k)/jacp
     .            +bz(i ,j,k)*bz_cov(ip,j,k)/jac0)/4.
cc        ptot = jac*(rho(ip,j,k)*tmp(i ,j,k)
cc     .             +rho(i ,j,k)*tmp(ip,j,k))
cc     .       +(bx(ip,j,k)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(ip,j,k)
cc     .        +by(ip,j,k)*by_cov(i,j,k)+by(i,j,k)*by_cov(ip,j,k)
cc     .        +bz(ip,j,k)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(ip,j,k))/4.

        vis = 2./(1./rho(ip,j,k)/nuu(ip,j,k)
     .          + 1./rho(i ,j,k)/nuu(i ,j,k))

        t11 =
     .      0.5*( rvx(ip,j,k)*vx (i ,j,k)/jacp/jac0
     .           +rvx(i ,j,k)*vx (ip,j,k)/jacp/jac0)*jac**2
     .     -    ( bx(ip,j,k)*bx(i ,j,k)/jacp/jac0 )*jac**2
     .     +gsuper(1,1)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,1)
     .           +gsuper(1,2)*nabla_v(2,1)
     .           +gsuper(1,3)*nabla_v(3,1) )

        t12 =
     .     0.25*( rvx(ip,j,k)*vy (i ,j,k)/jacp
     .           +rvx(i ,j,k)*vy (ip,j,k)/jac0
     .           + vx(i ,j,k)*rvy(ip,j,k)/jac0
     .           + vx(ip,j,k)*rvy(i ,j,k)/jacp)*jac
     .     -0.5*( bx(ip,j,k)*by(i ,j,k)/jacp
     .           +bx(i ,j,k)*by(ip,j,k)/jac0 )*jac
     .     +gsuper(1,2)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,2)
     .           +gsuper(1,2)*nabla_v(2,2)
     .           +gsuper(1,3)*nabla_v(3,2) )

        t13 =
     .      0.25*( rvx(ip,j,k)*vz (i ,j,k)/jacp/jac0
     .            +rvx(i ,j,k)*vz (ip,j,k)/jacp/jac0
     .            + vx(i ,j,k)*rvz(ip,j,k)/jacp/jac0
     .            + vx(ip,j,k)*rvz(i ,j,k)/jacp/jac0)*jac**2
     .     -0.5*( bx(ip,j,k)*bz(i ,j,k)/jacp/jac0
     .           +bx(i ,j,k)*bz(ip,j,k)/jacp/jac0)*jac**2
     .     +gsuper(1,3)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,3)
     .           +gsuper(1,2)*nabla_v(2,3)
     .           +gsuper(1,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t11 = t11/jac
          if (.not.alt_eom) t12 = t12/jac
          t13 = t13/jac
        endif

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

c     Begin program

        jp = j+1

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
          nabla_v = fnabla_v(i,j,k,x,y,z,vx,vy,vz,cartesian,2)
        else
          nabla_v = fnabla_v(i,j,k,x,y,z,vx,vy,vz,cartesian,0)
        endif

        !Recall p=2nT
        ptot = jac*(rho(i,jp,k)*tmp(i,jp,k)
     .             +rho(i,j ,k)*tmp(i,j ,k))
cc        ptot = jac*(rho(i,jp,k)*tmp(i,j ,k)
cc     .             +rho(i,j ,k)*tmp(i,jp,k))
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

        if (flag /= 0) then
          t21 = t21/jac
          if (.not.alt_eom) t22 = t22/jac
          t23 = t23/jac
        endif

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

c     Begin program

        kp=k+1

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
          nabla_v = fnabla_v(i,j,k,x,y,z,vx,vy,vz,cartesian,3)
        else
          nabla_v = fnabla_v(i,j,k,x,y,z,vx,vy,vz,cartesian,0)
        endif

        !Recall p=2nT
        ptot = jac*(rho(i,j,kp)*tmp(i,j,kp)
     .             +rho(i,j,k )*tmp(i,j,k ))
cc        ptot = jac*(rho(i,j,kp)*tmp(i,j,k )
cc     .             +rho(i,j,k )*tmp(i,j,kp))
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
     .     -0.5*( bz(i,j,kp)*bz(i,j,k )
     .           +bz(i,j,k )*bz(i,j,kp) )
     .     +gsuper(3,3)*ptot
     .     -vis*( gsuper(3,1)*nabla_v(1,3)
     .           +gsuper(3,2)*nabla_v(2,3)
     .           +gsuper(3,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alt_eom) t32 = t32/jac
          t33 = t33/jac
        endif

c     End program

      end subroutine vtensor_z

      end module nlfunction_setup

c module precond_setup
c ######################################################################
      module precond_setup

        integer(4)    :: nsweep,maxvcyc
        real(8)       :: mgtol

        character*(10):: precon

      end module precond_setup
