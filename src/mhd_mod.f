c map
c #####################################################################
      subroutine map(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1)

c ---------------------------------------------------------------------
c     Give Cartesian coordinates corresponding to node (i,j,k) at grid
c     level (igx,igy,igz).
c ---------------------------------------------------------------------

      implicit none

c Input variables

      integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg

c Local variables

      real(8)    :: x1,y1,z1

c Begin program

      write (*,*) 'Error: subroutine map called!'
      write (*,*) 'Aborting...'
      stop

      end subroutine map

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

        use grid

        real(8),pointer,dimension(:) :: xx,yy,zz,dxh,dyh,dzh,dx,dy,dz

        integer(4) :: igx,igy,igz,nx,ny,nz

        real(8)    :: xim,yim,zim,xip,yip,zip
     .               ,xjm,yjm,zjm,xjp,yjp,zjp
     .               ,xkm,ykm,zkm,xkp,ykp,zkp

        real(8)    :: x0,y0,z0,xh,yh,zh

        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .               ,jacp,jacm,jach,jac0

        real(8)    :: gsub(3,3),gsuper(3,3),jac

        real(8)    :: nabla_v(3,3),hessian(3,3,3)
     .               ,cov_tnsr(3,3),cnv_tnsr(3,3)

        logical    :: sing_point,cartesian

      end module grid_aliases

c module auxiliaryVariables
c ######################################################################
      module auxiliaryVariables

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
      real(8) function div(i,j,k,nx,ny,nz,ax,ay,az)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz
      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ig,jg,kg
      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacp,jac0,jach

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      jac0 = gmetric%grid(igx)%jac(i,j,k)

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      if (i == 1 .and. bcond(1) == SP) then
        jacp = gmetric%grid(igx)%jac(i+1,j,k)
        jach = 0.5*(jacp+jac0)   !Only good for cylindrical-like geom.

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
      real*8 function laplacian(i,j,k,nx,ny,nz,arr)

c     ---------------------------------------------------------------
c     Calculates dvol*lap(arr) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz

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
        d_xx_ip = 0.5*(gmetric%grid(igx)%gsup(i ,j,k,1,1)
     .                +gmetric%grid(igx)%gsup(ip,j,k,1,1))
        d_xx_im = 0.5*(gmetric%grid(igx)%gsup(i ,j,k,1,1)
     .                +gmetric%grid(igx)%gsup(im,j,k,1,1))
        d_yy_jp = 0.5*(gmetric%grid(igx)%gsup(i,j ,k,2,2)
     .                +gmetric%grid(igx)%gsup(i,jp,k,2,2))
        d_yy_jm = 0.5*(gmetric%grid(igx)%gsup(i,j ,k,2,2)
     .                +gmetric%grid(igx)%gsup(i,jm,k,2,2))
        d_zz_kp = 0.5*(gmetric%grid(igx)%gsup(i,j,k ,3,3)
     .                +gmetric%grid(igx)%gsup(i,j,kp,3,3))
        d_zz_km = 0.5*(gmetric%grid(igx)%gsup(i,j,k ,3,3)
     .                +gmetric%grid(igx)%gsup(i,j,km,3,3))

        d_xy_ipjp = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(ip,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,jp,k,1,2))
        d_xy_ipjm = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(ip,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,jm,k,1,2))
        d_xy_imjp = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(im,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,jp,k,1,2))
        d_xy_imjm = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(im,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,jm,k,1,2))
                 
        d_xz_ipkp = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(ip,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,kp,1,3))
        d_xz_ipkm = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(ip,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,km,1,3))
        d_xz_imkp = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(ip,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,km,1,3))
        d_xz_imkm = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(im,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,km,1,3))
                 
      else

        d_xx_ip = 0.5*(gmetric%grid(igx)%gsup(i ,j,k,1,1)
     .                +gmetric%grid(igx)%gsup(ip,j,k,1,1))
        d_xx_im = 0d0
        d_yy_jp = 0.5*(gmetric%grid(igx)%gsup(i,j ,k,2,2)
     .                +gmetric%grid(igx)%gsup(i,jp,k,2,2))
        d_yy_jm = 0.5*(gmetric%grid(igx)%gsup(i,j ,k,2,2)
     .                +gmetric%grid(igx)%gsup(i,jm,k,2,2))
        d_zz_kp = 0.5*(gmetric%grid(igx)%gsup(i,j,k ,3,3)
     .                +gmetric%grid(igx)%gsup(i,j,kp,3,3))
        d_zz_km = 0.5*(gmetric%grid(igx)%gsup(i,j,k ,3,3)
     .                +gmetric%grid(igx)%gsup(i,j,km,3,3))

        d_xy_ipjp = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(ip,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,jp,k,1,2))
        d_xy_ipjm = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(ip,j,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igx)%gsup(i,jm,k,1,2))
        d_xy_imjp = 0d0
        d_xy_imjm = 0d0

        d_xz_ipkp = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(ip,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,kp,1,3))
        d_xz_ipkm = 0.25*(gmetric%grid(igx)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(ip,j,k,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igx)%gsup(i,j,km,1,3))
        d_xz_imkp = 0d0
        d_xz_imkm = 0d0

      endif

      d_yz_jpkp = 0.25*(gmetric%grid(igx)%gsup(i,j ,k,2,3)
     .                 +gmetric%grid(igx)%gsup(i,jp,k,2,3)
     .                 +gmetric%grid(igx)%gsup(i,j,k ,2,3)
     .                 +gmetric%grid(igx)%gsup(i,j,kp,2,3))
      d_yz_jpkm = 0.25*(gmetric%grid(igx)%gsup(i,j ,k,2,3)
     .                 +gmetric%grid(igx)%gsup(i,jp,k,2,3)
     .                 +gmetric%grid(igx)%gsup(i,j,k ,2,3)
     .                 +gmetric%grid(igx)%gsup(i,j,km,2,3))
      d_yz_jmkp = 0.25*(gmetric%grid(igx)%gsup(i,j ,k,2,3)
     .                 +gmetric%grid(igx)%gsup(i,jm,k,2,3)
     .                 +gmetric%grid(igx)%gsup(i,j,k ,2,3)
     .                 +gmetric%grid(igx)%gsup(i,j,kp,2,3))
      d_yz_jmkm = 0.25*(gmetric%grid(igx)%gsup(i,j ,k,2,3)
     .                 +gmetric%grid(igx)%gsup(i,jm,k,2,3)
     .                 +gmetric%grid(igx)%gsup(i,j,k ,2,3)
     .                 +gmetric%grid(igx)%gsup(i,j,km,2,3))

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
      function veclaplacian(i,j,k,nx,ny,nz,ax,ay,az,diff,alteom,icomp)
     .         result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*lap(vector) at cell centers in general non-orthog.
c     coordinates.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,icomp,nx,ny,nz

      real(8)    :: ax  (0:nx+1,0:ny+1,0:nz+1)
     .             ,ay  (0:nx+1,0:ny+1,0:nz+1)
     .             ,az  (0:nx+1,0:ny+1,0:nz+1)
     .             ,diff(0:nx+1,0:ny+1,0:nz+1)

      real(8)    :: vlap

      logical    :: alteom

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: dvol,dS1,dS2,dS3,dxx,dyy,dzz,jac

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o

      real(8)    :: nabla_v(3,3),gsuper(3,3),hess(3,3,3),msource

      logical    :: sing_point

c     Begin program

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      sing_point = .false.
      if (i == 1 .and. bcond(1) == SP) sing_point = .true.

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      dS1 = dyy*dzz
      dS2 = dxx*dzz
      dS3 = dxx*dyy

      dvol = dxx*dyy*dzz

      jac  = gmetric%grid(igx)%jac(i,j,k)

      if (coords /= 'car') then
        hess = gmetric%grid(igx)%Gamma(i,j,k,:,:,:)
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
     .            *(t11o*hess(1,1,1)+t12o*hess(1,1,2)+t13o*hess(1,1,3)
     .             +t21o*hess(1,2,1)+t22o*hess(1,2,2)+t23o*hess(1,2,3)
     .             +t31o*hess(1,3,1)+t32o*hess(1,3,2)+t33o*hess(1,3,3))
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
     .             *(t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3)
     .              -t12o*(hess(1,1,1)+hess(2,1,2)+hess(3,1,3))
     .              -t22o*(hess(1,2,1)+hess(2,2,2)+hess(3,2,3))
     .              -t32o*(hess(1,3,1)+hess(2,3,2)+hess(3,3,3)))

            flxip = flxip/jac
            flxim = flxim/jac
                         
            flxjp = flxjp/jac
            flxjm = flxjm/jac
                         
            flxkp = flxkp/jac
            flxkm = flxkm/jac

          else
            msource=dvol
     .             *(t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3))
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
     .           *(t11o*hess(3,1,1)+t12o*hess(3,1,2)+t13o*hess(3,1,3)
     .            +t21o*hess(3,2,1)+t22o*hess(3,2,2)+t23o*hess(3,2,3)
     .            +t31o*hess(3,3,1)+t32o*hess(3,3,2)+t33o*hess(3,3,3))
        endif

      end select

      vlap = jac*( dS1*(flxip - flxim)
     .           + dS2*(flxjp - flxjm)
     .           + dS3*(flxkp - flxkm) ) + msource

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
        if (flag == 0 .or. (.not.alteom .and. sing_point) ) ip = i

        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
     .               +gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))

        if (i < nx .and. bcond(1) == SP .and. flag /= 0) then
          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac0 = gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,ax,ay,az,1)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,ax,ay,az,0)
        endif

        vis = 2d0/( 1d0/diff(ip,j,k) + 1d0/diff(i,j,k) )

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
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,ax,ay,az,2)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,ax,ay,az,0)
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
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,ax,ay,az,3)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,ax,ay,az,0)
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

c     End program

      end subroutine nabtensor_z

      end function veclaplacian

c     curl
c     ###############################################################
      real*8 function curl(i,j,k,nx,ny,nz,ax,ay,az,comp)

c     ---------------------------------------------------------------
c     Calculates curl(A) at cell centers in general non-orthogonal
c     coordinates. The vector A is covariant, and returns the 
c     contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,comp,nx,ny,nz

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg
      real(8)    :: dx,dy,dz,flxip,flxim,flxjp,flxjm,flxkp,flxkm
      real(8)    :: dS1,dS2,dS3,x0,y0,z0,jac
     .             ,xip,yip,zip,jacp,xh,yh,zh,jach

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
          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac  = gmetric%grid(igx)%jac(i ,j,k)
          jach = 0.5*(jacp+jac)

          flxip = 0.5*(ay(ip,j,k)/jacp+ay(i,j,k)/jac)*jach
          flxim = ay(im,j,k)

cc        elseif (i == 2 .and. bcond(1) == SP) then
cc          jacm = gmetric%grid(igx)%jac(im,j,k)
cc          jac  = gmetric%grid(igx)%jac(i ,j,k)
cc          jach = 0.5*(jacm+jac)
cc
cc          flxip = 0.5*(ay(ip,j,k)+ay(i,j,k))
cc          flxim = 0.5*(ay(im,j,k)/jacm+ay(i,j,k)/jac)*jach
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
      real*8 function curl2(i,j,k,nx,ny,nz,ax,ay,az,comp)

c     ---------------------------------------------------------------
c     Calculates curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,comp,nx,ny,nz

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg
      real(8)    :: dx,dy,dz,dx1,dx2,dy1,dy2,dz1,dz2
      real(8)    :: daxdz,dazdx,daydx,daxdy,dazdy,daydz
      real(8)    :: jac0,jacp,jacm

c     Begin program

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

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

          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac0 = gmetric%grid(igx)%jac(i ,j,k)
          jacm = 0.5*(jacp+jac0)

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

ccc     curlcurl
ccc     ###############################################################
cc      real*8 function curlcurl(i,j,k,ax,ay,az,diff,comp)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates curl(eta curl(A)) in general non-orthogonal
ccc     coordinates, preserving the SPD property. The vector A is
ccc     covariant, and returns the contravariant component "comp".
ccc     ---------------------------------------------------------------
cc
cc      use grid
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer(4) :: i,j,k,comp
cc
cc      real(8)    :: ax  (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,ay  (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,az  (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,diff(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg
cc
cc      real(8)    :: x0,xp,xm,y0,yp,ym,z0,zp,zm
cc      real(8)    :: dwa,dwb,dwc,dwd,dwe,dwf,dwg,dwh,dwi,dwj,dwk
cc
cc      !Component one
cc      real(8)    :: d_yy_kp,d_yy_km,d_zz_jp,d_zz_jm
cc     .             ,d_xy_kp,d_xy_km,d_yz_kp,d_yz_km
cc     .             ,d_xz_jp,d_xz_jm,d_yz_jp,d_yz_jm
cc     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm
cc
cc      !Component two
cc      real(8)    :: d_zz_ip,d_zz_im,d_xx_kp,d_xx_km
cc     .             ,d_xz_kp,d_xz_km
cc     .             ,d_xz_ip,d_xz_im,d_yz_ip,d_yz_im
cc     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
cc
cc      !Component three
cc      real(8)    :: d_yy_ip,d_yy_im,d_xx_jp,d_xx_jm
cc     .             ,d_xy_jp,d_xy_jm,d_xy_ip,d_xy_im
cccc     .             ,d_xz_jp,d_xz_jm,d_yz_ip,d_yz_im
cc     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
cc
ccc     Begin program
cc
cc      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc      ip = i+1
cc      im = i-1
cc      jp = j+1
cc      jm = j-1
cc      kp = k+1
cc      km = k-1
cc
cc      x0 = xx(ig)
cc      xp = 0.5*(xx(ig+1)+xx(ig))
cc      xm = 0.5*(xx(ig-1)+xx(ig))
cc
cc      y0 = yy(jg)
cc      yp = 0.5*(yy(jg+1)+yy(jg))
cc      ym = 0.5*(yy(jg-1)+yy(jg))
cc
cc      z0 = zz(kg)
cc      zp = 0.5*(zz(kg+1)+zz(kg))
cc      zm = 0.5*(zz(kg-1)+zz(kg))
cc
cc      select case(comp)
cc      case(1)
cc
cc        d_yy_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,2,x0,y0,zp,.false.)
cc        d_yy_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,2,x0,y0,zm,.false.)
cc
cc        d_zz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(3,3,x0,yp,z0,.false.)
cc        d_zz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(3,3,x0,ym,z0,.false.)
cc
cc        d_yz_jpkp = 4./(1./diff(i,jp,kp)+1./diff(i,jp,k )
cc     .                 +1./diff(i,j ,kp)+1./diff(i,j ,k ))
cc     .             *g_sub_elem(2,3,x0,yp,zp,.false.)
cc        d_yz_jpkm = 4./(1./diff(i,jp,km)+1./diff(i,jp,k )
cc     .                 +1./diff(i,j ,km)+1./diff(i,j ,k ))
cc     .             *g_sub_elem(2,3,x0,yp,zm,.false.)
cc        d_yz_jmkp = 4./(1./diff(i,jm,kp)+1./diff(i,jm,k )
cc     .                 +1./diff(i,j ,kp)+1./diff(i,j ,k ))
cc     .             *g_sub_elem(2,3,x0,ym,zp,.false.)
cc        d_yz_jmkm = 4./(1./diff(i,jm,km)+1./diff(i,jm,k )
cc     .                 +1./diff(i,j ,km)+1./diff(i,j ,k ))
cc     .             *g_sub_elem(2,3,x0,ym,zm,.false.)
cc
cc        d_xy_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,2,x0,y0,zp,.false.)
cc        d_xy_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,2,x0,y0,zm,.false.)
cc
cc        d_yz_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,3,x0,y0,zp,.false.)
cc        d_yz_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,3,x0,y0,zm,.false.)
cc
cc        d_xz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,3,x0,yp,z0,.false.)
cc        d_xz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,3,x0,ym,z0,.false.)
cc
cc        d_yz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,3,x0,yp,z0,.false.)
cc        d_yz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,3,x0,ym,z0,.false.)
cc
cc        dwa = -dxh(ig)*dyh(jg)
cc     .              *(d_yy_kp*(ax(i,j,kp)-ax(i,j,k ))/dz(kg  )
cc     .               -d_yy_km*(ax(i,j,k )-ax(i,j,km))/dz(kg-1) )
cc        dwb = -dxh(ig)*dzh(kg)
cc     .              *(d_zz_jp*(ax(i,jp,k)-ax(i,j ,k))/dy(jg  )
cc     .               -d_zz_jm*(ax(i,j ,k)-ax(i,jm,k))/dy(jg-1) )
cc        dwc = -dxh(ig)/2.
cc     .              *(-d_yz_jpkp*(ax(i,jp,kp)-ax(i,j ,k ))
cc     .                +d_yz_jmkm*(ax(i,j ,k )-ax(i,jm,km))
cc     .                +d_yz_jpkm*(ax(i,jp,km)-ax(i,j ,k ))
cc     .                -d_yz_jmkp*(ax(i,j ,k )-ax(i,jm,kp)))
cc        dwd =  dyh(jg)/4.
cc     .              *(d_yy_kp*(az(ip,j,kp)-az(im,j,kp)
cc     .                        +az(ip,j,k )-az(im,j,k ))
cc     .               -d_yy_km*(az(ip,j,km)-az(im,j,km)
cc     .                        +az(ip,j,k )-az(im,j,k )))
cc        dwe = -dxh(ig)/4.
cc     .              *(d_xy_kp*(az(i,jp,kp)-az(i,jm,kp)
cc     .                        +az(i,jp,k )-az(i,jm,k ))
cc     .               -d_xy_km*(az(i,jp,km)-az(i,jm,km)
cc     .                        +az(i,jp,k )-az(i,jm,k )))
cc        dwf =  dxh(ig)*dyh(jg)
cc     .              *(d_xy_kp*(ay(i,j,kp)-ay(i,j,k ))/dz(kg)
cc     .               -d_xy_km*(ay(i,j,k )-ay(i,j,km))/dz(kg-1))
cc        dwg = -dyh(jg)/4.
cc     .              *(d_yz_kp*(ay(ip,j,kp)-ay(im,j,kp)
cc     .                        +ay(ip,j,k )-ay(im,j,k ))
cc     .               -d_yz_km*(ay(ip,j,km)-ay(im,j,km)
cc     .                        +ay(ip,j,k )-ay(im,j,k )))
cc        dwh =  dzh(kg)/4.
cc     .              *(d_zz_jp*(ay(ip,jp,k)-ay(im,jp,k)
cc     .                        +ay(ip,j ,k)-ay(im,j ,k))
cc     .               -d_zz_jm*(ay(ip,jm,k)-ay(im,jm,k)
cc     .                        +ay(ip,j ,k)-ay(im,j ,k)))
cc        dwi =  dxh(ig)*dzh(kg)
cc     .              *(d_xz_jp*(az(i,jp,k)-az(i,j ,k))/dy(jg  )
cc     .               -d_xz_jm*(az(i,j ,k)-az(i,jm,k))/dy(jg-1))
cc        dwj = -dxh(ig)/4.
cc     .              *(d_xz_jp*(ay(i,jp,kp)-ay(i,jp,km)
cc     .                        +ay(i,j ,kp)-ay(i,j ,km))
cc     .               -d_xz_jm*(ay(i,jm,kp)-ay(i,jm,km)
cc     .                        +ay(i,j ,kp)-ay(i,j ,km)))
cc        dwk = -dzh(kg)/4.
cc     .              *(d_yz_jp*(az(ip,jp,k)-az(im,jp,k)
cc     .                        +az(ip,j ,k)-az(im,j ,k))
cc     .               -d_yz_jm*(az(ip,jm,k)-az(im,jm,k)
cc     .                        +az(ip,j ,k)-az(im,j ,k)))
cc
cc      case(2)
cc
cc        d_xx_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,1,x0,y0,zp,.false.)
cc        d_xx_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,1,x0,y0,zm,.false.)
cc
cc        d_zz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(3,3,xp,y0,z0,.false.)
cc        d_zz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(3,3,xm,y0,z0,.false.)
cc
cc        d_xz_ipkp = 4./(1./diff(ip,j,kp)+1./diff(ip,j,k )
cc     .                 +1./diff(i ,j,kp)+1./diff(i ,j,k ))
cc     .             *g_sub_elem(1,3,xp,y0,zp,.false.)
cc        d_xz_ipkm = 4./(1./diff(ip,j,km)+1./diff(ip,j,k )
cc     .                 +1./diff(i ,j,km)+1./diff(i ,j,k ))
cc     .             *g_sub_elem(1,3,xp,y0,zm,.false.)
cc        d_xz_imkp = 4./(1./diff(im,j,kp)+1./diff(im,j,k )
cc     .                 +1./diff(i ,j,kp)+1./diff(i ,j,k ))
cc     .             *g_sub_elem(1,3,xm,y0,zp,.false.)
cc        d_xz_imkm = 4./(1./diff(im,j,km)+1./diff(im,j,k )
cc     .                 +1./diff(i ,j,km)+1./diff(i ,j,k ))
cc     .             *g_sub_elem(1,3,xm,y0,zm,.false.)
cc
cc        d_xy_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,2,x0,y0,zp,.false.)
cc        d_xy_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,2,x0,y0,zm,.false.)
cc
cc        d_yz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,3,xp,y0,z0,.false.)
cc        d_yz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,3,xm,y0,z0,.false.)
cc
cc        d_xz_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,3,x0,y0,zp,.false.)
cc        d_xz_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,3,x0,y0,zm,.false.)
cc
cc        d_xz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,3,xp,y0,z0,.false.)
cc        d_xz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,3,xm,y0,z0,.false.)
cc
cc        dwa = -dzh(kg)*dyh(jg)
cc     .              *(d_zz_ip*(ay(ip,j,k)-ay(i ,j,k))/dx(ig  )
cc     .               -d_zz_im*(ay(i ,j,k)-ay(im,j,k))/dx(ig-1) )
cc        dwb = -dxh(ig)*dyh(jg)
cc     .              *(d_xx_kp*(ay(i,j,kp)-ay(i,j,k ))/dz(kg  )
cc     .               -d_xx_km*(ay(i,j,k )-ay(i,j,km))/dz(kg-1) )
cc        dwc = -dyh(jg)/2.
cc     .              *(-d_xz_ipkp*(ay(ip,j,kp)-ay(i ,j,k ))
cc     .                +d_xz_imkm*(ay(i ,j,k )-ay(im,j,km))
cc     .                +d_xz_ipkm*(ay(ip,j,km)-ay(i ,j,k ))
cc     .                -d_xz_imkp*(ay(i ,j,k )-ay(im,j,kp)))
cc        dwd =  dxh(ig)/4.
cc     .              *(d_xx_kp*(az(i,jp,kp)-az(i,jm,kp)
cc     .                        +az(i,jp,k )-az(i,jm,k ))
cc     .               -d_xx_km*(az(i,jp,km)-az(i,jm,km)
cc     .                        +az(i,jp,k )-az(i,jm,k )))
cc        dwe = -dyh(jg)/4.
cc     .              *(d_xy_kp*(az(ip,j,kp)-az(im,j,kp)
cc     .                        +az(ip,j,k )-az(im,j,k ))
cc     .               -d_xy_km*(az(ip,j,km)-az(im,j,km)
cc     .                        +az(ip,j,k )-az(im,j,k )))
cc        dwf =  dxh(ig)*dyh(jg)
cc     .              *(d_xy_kp*(ax(i,j,kp)-ax(i,j,k ))/dz(kg)
cc     .               -d_xy_km*(ax(i,j,k )-ax(i,j,km))/dz(kg-1))
cc        dwg = -dxh(ig)/4.
cc     .              *(d_xz_kp*(ax(i,jp,kp)-ax(i,jm,kp)
cc     .                        +ax(i,jp,k )-ax(i,jm,k ))
cc     .               -d_xz_km*(ax(i,jp,km)-ax(i,jm,km)
cc     .                        +ax(i,jp,k )-ax(i,jm,k )))
cc        dwh =  dzh(kg)/4.
cc     .              *(d_zz_ip*(ax(ip,jp,k)-ax(ip,jm,k)
cc     .                        +ax(i ,jp,k)-ax(i ,jm,k))
cc     .               -d_zz_im*(ax(im,jp,k)-ax(im,jm,k)
cc     .                        +ax(i ,jp,k)-ax(i ,jm,k)))
cc        dwi =  dyh(jg)*dzh(kg)
cc     .              *(d_yz_ip*(az(ip,j,k)-az(i ,j,k))/dx(ig  )
cc     .               -d_yz_im*(az(i ,j,k)-az(im,j,k))/dx(ig-1))
cc        dwj = -dzh(kg)/4.
cc     .              *(d_xz_ip*(az(ip,jp,k)-az(ip,jm,k)
cc     .                        +az(i ,jp,k)-az(i ,jm,k))
cc     .               -d_xz_im*(az(im,jp,k)-az(im,jm,k)
cc     .                        +az(i ,jp,k)-az(i ,jm,k)))
cc        dwk = -dyh(jg)/4.
cc     .              *(d_yz_ip*(ax(ip,j,kp)-ax(ip,j,km)
cc     .                        +ax(i ,j,kp)-ax(i ,j,km))
cc     .               -d_yz_im*(ax(im,j,kp)-ax(im,j,km)
cc     .                        +ax(i ,j,kp)-ax(i ,j,km)))
cc
cc      case(3)
cc
cc        d_xx_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,1,x0,yp,z0,.false.)
cc        d_xx_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,1,x0,ym,z0,.false.)
cc
cc        d_yy_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,2,xp,y0,z0,.false.)
cc        d_yy_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,2,xm,y0,z0,.false.)
cc
cc        d_xy_ipjp = 4./(1./diff(ip,jp,k)+1./diff(ip,j,k )
cc     .                 +1./diff(i ,jp,k)+1./diff(i ,j,k ))
cc     .             *g_sub_elem(1,2,xp,yp,z0,.false.)
cc        d_xy_ipjm = 4./(1./diff(ip,jm,k)+1./diff(ip,j,k )
cc     .                 +1./diff(i ,jm,k)+1./diff(i ,j,k ))
cc     .             *g_sub_elem(1,2,xp,ym,z0,.false.)
cc        d_xy_imjp = 4./(1./diff(im,jp,k)+1./diff(im,j,k )
cc     .                 +1./diff(i ,jp,k)+1./diff(i ,j,k ))
cc     .             *g_sub_elem(1,2,xm,yp,z0,.false.)
cc        d_xy_imjm = 4./(1./diff(im,jm,k)+1./diff(im,j,k )
cc     .                 +1./diff(i ,jm,k)+1./diff(i ,j,k ))
cc     .             *g_sub_elem(1,2,xm,ym,z0,.false.)
cc
cc        d_xy_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,2,x0,yp,z0,.false.)
cc        d_xy_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,2,x0,ym,z0,.false.)
cc
cc        d_yz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,3,xp,y0,z0,.false.)
cc        d_yz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(2,3,xm,y0,z0,.false.)
cc
cc        d_xy_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,2,xp,y0,z0,.false.)
cc        d_xy_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,2,xm,y0,z0,.false.)
cc
cc        d_xz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,3,x0,yp,z0,.false.)
cc        d_xz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .           *g_sub_elem(1,3,x0,ym,z0,.false.)
cc
cc        dwa = -dzh(kg)*dyh(jg)
cc     .              *(d_yy_ip*(az(ip,j,k)-az(i ,j,k))/dx(ig  )
cc     .               -d_yy_im*(az(i ,j,k)-az(im,j,k))/dx(ig-1) )
cc        dwb = -dxh(ig)*dzh(kg)
cc     .              *(d_xx_jp*(az(i,jp,k)-az(i,j ,k))/dy(jg  )
cc     .               -d_xx_jm*(az(i,j ,k)-az(i,jm,k))/dy(jg-1) )
cc        dwc = -dzh(kg)/2.
cc     .              *(-d_xy_ipjp*(az(ip,jp,k)-az(i ,j ,k))
cc     .                +d_xy_imjm*(az(i ,j ,k)-az(im,jm,k))
cc     .                +d_xy_ipjm*(az(ip,jm,k)-az(i ,j ,k))
cc     .                -d_xy_imjp*(az(i ,j ,k)-az(im,jp,k)))
cc        dwd =  dxh(ig)/4.
cc     .              *(d_xx_jp*(ay(i,jp,kp)-ay(i,jp,km)
cc     .                        +ay(i,j ,kp)-ay(i,j ,km))
cc     .               -d_xx_jm*(ay(i,jm,kp)-ay(i,jm,km)
cc     .                        +ay(i,j ,kp)-ay(i,j ,km)))
cc        dwe=  -dxh(ig)/4.
cc     .              *(d_xy_jp*(ax(i,jp,kp)-ax(i,jp,km)
cc     .                        +ax(i,j ,kp)-ax(i,j ,km))
cc     .               -d_xy_jm*(ax(i,jm,kp)-ax(i,jm,km)
cc     .                        +ax(i,j ,kp)-ax(i,j ,km)))
cc        dwf =  dxh(ig)*dzh(kg)
cc     .              *(d_xz_jp*(ax(i,jp,k)-ax(i,j ,k))/dy(jg)
cc     .               -d_xz_jm*(ax(i,j ,k)-ax(i,jm,k))/dy(jg-1))
cc        dwg = -dzh(kg)/4.
cc     .              *(d_xz_jp*(ay(ip,jp,k)-ay(im,jp,k)
cc     .                        +ay(ip,j ,k)-ay(im,j ,k))
cc     .               -d_xz_jm*(ay(ip,jm,k)-ay(im,jm,k)
cc     .                        +ay(ip,j ,k)-ay(im,j ,k)))
cc        dwh =  dyh(jg)/4.
cc     .              *(d_yy_ip*(ax(ip,j,kp)-ax(ip,j,km)
cc     .                        +ax(i ,j,kp)-ax(i ,j,km))
cc     .               -d_yy_im*(ax(im,j,kp)-ax(im,j,km)
cc     .                        +ax(i ,j,kp)-ax(i ,j,km)))
cc        dwi =  dyh(jg)*dzh(kg)
cc     .              *(d_yz_ip*(ay(ip,j,k)-ay(i ,j,k))/dx(ig  )
cc     .               -d_yz_im*(ay(i ,j,k)-ay(im,j,k))/dx(ig-1))
cc        dwj = -dzh(kg)/4.
cc     .              *(d_yz_ip*(ax(ip,jp,k)-ax(ip,jm,k)
cc     .                        +ax(i ,jp,k)-ax(i ,jm,k))
cc     .               -d_yz_im*(ax(im,jp,k)-ax(im,jm,k)
cc     .                        +ax(i ,jp,k)-ax(i ,jm,k)))
cc        dwk = -dyh(jg)/4.
cc     .              *(d_xy_ip*(ay(ip,j,kp)-ay(ip,j,km)
cc     .                        +ay(i ,j,kp)-ay(i ,j,km))
cc     .               -d_xy_im*(ay(im,j,kp)-ay(im,j,km)
cc     .                        +ay(i ,j,kp)-ay(i ,j,km)))
cc      case default
cc
cc        write (*,*) 'Error in component in curlcurl'
cc        write (*,*) 'Aborting...'
cc        stop
cc
cc      end select
cc
cc      curlcurl = dwa+dwb+dwc+dwd+dwe+dwf+dwg+dwh+dwi+dwj+dwk
cc
ccc     End program
cc
cc      end function curlcurl

c     fnabla_v
c     #############################################################
      function fnabla_v(i,j,k,nx,ny,nz,ax,ay,az,half_elem)
     .         result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor nabla(vec v) at the following positions:
c       + half_elem=0 --> i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,half_elem,nx,ny,nz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: tsrc(3,3)
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

          tsrc = 0.5*(nabla_v_src(ip,j,k)
     .               +nabla_v_src(i ,j,k))

        case (2)
          dhy = dy(jg)

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

          tsrc = 0.5*(nabla_v_src(i,jp,k)
     .               +nabla_v_src(i,j ,k))

        case (3)
          dhz = dz(kg)

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

          tsrc = 0.5*(nabla_v_src(i,j,kp)
     .               +nabla_v_src(i,j,k ))

        case default

          if (sing_point) then
            vxip = ax(ip,j,k)+ax(i,j,k)
            vxim = 2.*ax(im,j,k)
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

          tsrc = nabla_v_src(i,j,k)

        end select

        idhx = 1./dhx
        idhy = 1./dhy
        idhz = 1./dhz

      ! l = 1, m = 1
        tensor(1,1) = (vxip-vxim)*idhx

      ! l = 1, m = 2
        tensor(1,2) = (vyip-vyim)*idhx

      ! l = 1, m = 3
        tensor(1,3) = (vzip-vzim)*idhx

      ! l = 2, m = 1
        tensor(2,1) = (vxjp-vxjm)*idhy

      ! l = 2, m = 2
        tensor(2,2) = (vyjp-vyjm)*idhy

      ! l = 2, m = 3
        tensor(2,3) = (vzjp-vzjm)*idhy

      ! l = 3, m = 1
        tensor(3,1) = (vxkp-vxkm)*idhz

      ! l = 3, m = 2
        tensor(3,2) = (vykp-vykm)*idhz

      ! l = 3, m = 3
        tensor(3,3) = (vzkp-vzkm)*idhz

      ! Add geometric source

        tensor = tensor - tsrc

c     End program

      contains

c     nabla_v_src
c     #############################################################
      function nabla_v_src(i,j,k) result(tensor)

c     -------------------------------------------------------------
c     Finds geometric source of tensor nabla(v) at cell (i,j,k)
c     -------------------------------------------------------------

        implicit none

        integer(4) :: i,j,k
        real(8)    :: tensor(3,3)

        real(8)    :: hessian(3,3,3)
     .               ,vxx,vyy,vzz

c     Begin program

        vxx = ax(i,j,k)
        vyy = ay(i,j,k)
        vzz = az(i,j,k)

        hessian = gmetric%grid(igx)%Gamma(i,j,k,:,:,:)

cc      ! l = 1, m = 1
cc        tensor(1,1) =  vxx*(hessian(1,1,1)
cc     .                    + hessian(2,2,1)
cc     .                    + hessian(3,3,1))
cc     .               - vxx *hessian(1,1,1)
cc     .               - vyy *hessian(1,2,1)
cc     .               - vzz *hessian(1,3,1)
cc
cc      ! l = 1, m = 2
cc        tensor(1,2) =  vyy*(hessian(1,1,1)
cc     .                    + hessian(2,2,1)
cc     .                    + hessian(3,3,1))
cc     .               - vxx* hessian(2,1,1)
cc     .               - vyy* hessian(2,2,1)
cc     .               - vzz* hessian(2,3,1)
cc
cc      ! l = 1, m = 3
cc        tensor(1,3) =  vzz*(hessian(1,1,1)
cc     .                    + hessian(2,2,1)
cc     .                    + hessian(3,3,1))
cc     .               - vxx* hessian(3,1,1)
cc     .               - vyy* hessian(3,2,1)
cc     .               - vzz* hessian(3,3,1)
cc
cc      ! l = 2, m = 1
cc        tensor(2,1) =  vxx*(hessian(1,1,2)
cc     .                    + hessian(2,2,2)
cc     .                    + hessian(3,3,2))
cc     .               - vxx* hessian(1,1,2)
cc     .               - vyy* hessian(1,2,2)
cc     .               - vzz* hessian(1,3,2)
cc
cc      ! l = 2, m = 2
cc        tensor(2,2) =  vyy*(hessian(1,1,2)
cc     .                    + hessian(2,2,2)
cc     .                    + hessian(3,3,2))
cc     .               - vxx* hessian(2,1,2)
cc     .               - vyy* hessian(2,2,2)
cc     .               - vzz* hessian(2,3,2)
cc
cc      ! l = 2, m = 3
cc        tensor(2,3) =  vzz*(hessian(1,1,2)
cc     .                    + hessian(2,2,2)
cc     .                    + hessian(3,3,2))
cc     .               - vxx* hessian(3,1,2)
cc     .               - vyy* hessian(3,2,2)
cc     .               - vzz* hessian(3,3,2)
cc
cc      ! l = 3, m = 1
cc        tensor(3,1) =  vxx*(hessian(1,1,3)
cc     .                    + hessian(2,2,3)
cc     .                    + hessian(3,3,3))
cc     .               - vxx* hessian(1,1,3)
cc     .               - vyy* hessian(1,2,3)
cc     .               - vzz* hessian(1,3,3)
cc
cc      ! l = 3, m = 2
cc        tensor(3,2) =  vyy*(hessian(1,1,3)
cc     .                    + hessian(2,2,3)
cc     .                    + hessian(3,3,3))
cc     .               - vxx* hessian(2,1,3)
cc     .               - vyy* hessian(2,2,3)
cc     .               - vzz* hessian(2,3,3)
cc
cc      ! l = 3, m = 3
cc        tensor(3,3) =  vzz*(hessian(1,1,3)
cc     .                    + hessian(2,2,3)
cc     .                    + hessian(3,3,3))
cc     .               - vxx* hessian(3,1,3)
cc     .               - vyy* hessian(3,2,3)
cc     .               - vzz* hessian(3,3,3)

      ! l = 1, m = 1
        tensor(1,1) =  vxx*(hessian(2,2,1)
     .                    + hessian(3,3,1))
     .               - vyy *hessian(1,2,1)
     .               - vzz *hessian(1,3,1)

      ! l = 1, m = 2
        tensor(1,2) =  vyy*(hessian(1,1,1)
     .                    + hessian(3,3,1))
     .               - vxx* hessian(2,1,1)
     .               - vzz* hessian(2,3,1)

      ! l = 1, m = 3
        tensor(1,3) =  vzz*(hessian(1,1,1)
     .                    + hessian(2,2,1))
     .               - vxx* hessian(3,1,1)
     .               - vyy* hessian(3,2,1)

      ! l = 2, m = 1
        tensor(2,1) =  vxx*(hessian(2,2,2)
     .                    + hessian(3,3,2))
     .               - vyy* hessian(1,2,2)
     .               - vzz* hessian(1,3,2)

      ! l = 2, m = 2
        tensor(2,2) =  vyy*(hessian(1,1,2)
     .                    + hessian(3,3,2))
     .               - vxx* hessian(2,1,2)
     .               - vzz* hessian(2,3,2)

      ! l = 2, m = 3
        tensor(2,3) =  vzz*(hessian(1,1,2)
     .                    + hessian(2,2,2))
     .               - vxx* hessian(3,1,2)
     .               - vyy* hessian(3,2,2)

      ! l = 3, m = 1
        tensor(3,1) =  vxx*(hessian(2,2,3)
     .                    + hessian(3,3,3))
     .               - vyy* hessian(1,2,3)
     .               - vzz* hessian(1,3,3)

      ! l = 3, m = 2
        tensor(3,2) =  vyy*(hessian(1,1,3)
     .                    + hessian(3,3,3))
     .               - vxx* hessian(2,1,3)
     .               - vzz* hessian(2,3,3)

      ! l = 3, m = 3
        tensor(3,3) =  vzz*(hessian(1,1,3)
     .                    + hessian(2,2,3))
     .               - vxx* hessian(3,1,3)
     .               - vyy* hessian(3,2,3)

c     End program

      end function nabla_v_src

      end function fnabla_v

c     fnabla_v_upwd
c     #############################################################
      function fnabla_v_upwd(i,j,k,nx,ny,nz,ax,ay,az,hex,hey,hez)
     $         result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor nabla(vec v) at the following positions:
c       + hex,hey,hez = 0 => i,j,k
c       + hex=+-1 --> i+-1/2
c       + hey=+-1 --> j+-1/2
c       + hez=+-1 --> k+-1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,hex,hey,hez,nx,ny,nz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: dhx,dhy,dhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: hessian(3,3,3)
        logical    :: sing_point,cartsn

c     Begin program

        sing_point = .false.
        if (i == 1 .and. bcond(1) == SP) sing_point = .true.

c     Defaults

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        dhx = 2.*dxh(ig)
        dhy = 2.*dyh(jg)
        dhz = 2.*dzh(kg)

        hessian = -gmetric%grid(igx)%Gamma(i,j,k,:,:,:)

c     Exceptions

        if (hex == 1) then
          im = i
          dhx = dx(ig)
        elseif (hex == -1) then
          ip = i
          dhx = dx(ig-1)
        endif

        if (hey == 1) then
          jm = j
          dhy = dy(jg)
        elseif (hey == -1) then
          jp = j
          dhy = dy(jg-1)
        endif

        if (hez == 1) then
          km = k
          dhz = dz(kg)
        elseif (hez == -1) then
          kp = k
          dhz = dz(kg-1)
        endif

c     Vectors

        vxx = ax(i,j,k)
        vyy = ay(i,j,k)
        vzz = az(i,j,k)

        if (sing_point) then
          vxip = ax(ip,j,k)+ax(i,j,k)
          vxim = 2.*ax(im,j,k)
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

c     Calculate nabla_v tensor

      ! l = 1, m = 1
        tensor(1,1) = (vxip-vxim)/dhx
     .               + vxx*(hessian(1,1,1)
     .                    + hessian(2,2,1)
     .                    + hessian(3,3,1))
     .               - vxx* hessian(1,1,1)
     .               - vyy* hessian(1,2,1)
     .               - vzz* hessian(1,3,1)

      ! l = 1, m = 2
        tensor(1,2) = (vyip-vyim)/dhx
     .               + vyy*(hessian(1,1,1)
     .                    + hessian(2,2,1)
     .                    + hessian(3,3,1))
     .               - vxx* hessian(2,1,1)
     .               - vyy* hessian(2,2,1)
     .               - vzz* hessian(2,3,1)

      ! l = 1, m = 3
        tensor(1,3) = (vzip-vzim)/dhx
     .               + vzz*(hessian(1,1,1)
     .                    + hessian(2,2,1)
     .                    + hessian(3,3,1))
     .               - vxx* hessian(3,1,1)
     .               - vyy* hessian(3,2,1)
     .               - vzz* hessian(3,3,1)

      ! l = 2, m = 1
        tensor(2,1) = (vxjp-vxjm)/dhy
     .               + vxx*(hessian(1,1,2)
     .                    + hessian(2,2,2)
     .                    + hessian(3,3,2))
     .               - vxx* hessian(1,1,2)
     .               - vyy* hessian(1,2,2)
     .               - vzz* hessian(1,3,2)

      ! l = 2, m = 2
        tensor(2,2) = (vyjp-vyjm)/dhy
     .               + vyy*(hessian(1,1,2)
     .                    + hessian(2,2,2)
     .                    + hessian(3,3,2))
     .               - vxx* hessian(2,1,2)
     .               - vyy* hessian(2,2,2)
     .               - vzz* hessian(2,3,2)

      ! l = 2, m = 3
        tensor(2,3) = (vzjp-vzjm)/dhy
     .               + vzz*(hessian(1,1,2)
     .                    + hessian(2,2,2)
     .                    + hessian(3,3,2))
     .               - vxx* hessian(3,1,2)
     .               - vyy* hessian(3,2,2)
     .               - vzz* hessian(3,3,2)

      ! l = 3, m = 1
        tensor(3,1) = (vxkp-vxkm)/dhz
     .               + vxx*(hessian(1,1,3)
     .                    + hessian(2,2,3)
     .                    + hessian(3,3,3))
     .               - vxx* hessian(1,1,3)
     .               - vyy* hessian(1,2,3)
     .               - vzz* hessian(1,3,3)

      ! l = 3, m = 2
        tensor(3,2) = (vykp-vykm)/dhz
     .               + vyy*(hessian(1,1,3)
     .                    + hessian(2,2,3)
     .                    + hessian(3,3,3))
     .               - vxx* hessian(2,1,3)
     .               - vyy* hessian(2,2,3)
     .               - vzz* hessian(2,3,3)

      ! l = 3, m = 3
        tensor(3,3) = (vzkp-vzkm)/dhz
     .               + vzz*(hessian(1,1,3)
     .                    + hessian(2,2,3)
     .                    + hessian(3,3,3))
     .               - vxx* hessian(3,1,3)
     .               - vyy* hessian(3,2,3)
     .               - vzz* hessian(3,3,3)

c     End program

      end function fnabla_v_upwd

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

        logical :: alt_eom

      contains

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
        if (flag == 0 .or. (.not.alt_eom .and. sing_point) ) ip = i

        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
     .               +gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))

        if (i < nx .and. bcond(1) == SP .and. flag /= 0) then
          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac0 = gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,vx,vy,vz,1)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,vx,vy,vz,0)
        endif

        !Recall p=2nT
        ptot = jac*(rho(ip,j,k)*tmp(ip,j,k)
     .             +rho(i ,j,k)*tmp(i ,j,k))
cc        ptot = jac*(rho(ip,j,k)*tmp(i ,j,k)
cc     .             +rho(i ,j,k)*tmp(ip,j,k))
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
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,vx,vy,vz,2)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,vx,vy,vz,0)
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
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,vx,vy,vz,3)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,vx,vy,vz,0)
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

        integer(4)    :: nsweep,maxvcyc,precpass
        real(8)       :: mgtol

        character*(10):: precon

      end module precond_setup
