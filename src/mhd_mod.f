c map
c #####################################################################
      subroutine map(i,j,k,igx,igy,igz,ig,jg,kg,x1,y1,z1)

c ---------------------------------------------------------------------
c     Give Cartesian coordinates corresponding to node (i,j,k) at grid
c     level (igx,igy,igz).
c ---------------------------------------------------------------------

      use error

      implicit none

c Input variables

      integer(4) :: i,j,k,igx,igy,igz,ig,jg,kg

c Local variables

      real(8)    :: x1,y1,z1

c Begin program

      messg = 'This subroutine should not be called here'
      call pstop('map',messg)

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

      contains

c     res
c     #############################################################
      function res(i,j,k,nx,ny,nz,igx,igy,igz)
c     -------------------------------------------------------------
c     This function computes the resistivity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: res
      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz

c     Local variables

      integer(4) :: nn,ig,jg,kg
      real(8)    :: x1,y1,z1,aa
      logical    :: cartsn

c     Begin program

c     Resistivity profile eta*(1 + aa*x^nn)
c       Coefficient aa is set so that res = 20*eta at wall
c       and nn so that res=??*eta at sing. surf. xs ~ 0.33

cc      select case (equil)
cc      case ('rfp1')
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc        nn = 4
cc        aa = 19.
cc        res = eta*(1. + aa*grid_params%xx(ig)**nn)
cc
cc      case default

        res = eta

cc      end select

c     End program

      end function res

c     vis
c     #############################################################
      function vis(i,j,k,nx,ny,nz,igx,igy,igz)
c     -------------------------------------------------------------
c     This function computes the viscosity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: vis
      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz

c     Local variables

c     Begin program

      vis = nu

c     End program

      end function vis

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

        logical    :: cartesian

      end module grid_aliases

c module auxiliaryVariables
c ######################################################################
      module auxiliaryVariables

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_cov,by_cov,bz_cov
     .         ,jx,jy,jz,jx_cov,jy_cov,jz_cov
     .         ,vx,vy,vz,vx_cov,vy_cov,vz_cov
     .         ,eeta,nuu,ejx,ejy,ejz

        real(8),target,allocatable,dimension(:,:,:,:) :: bcnv,vcnv

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_car,by_car,bz_car
     .         ,jx_car,jy_car,jz_car
     .         ,vx_car,vy_car,vz_car
     .         ,divrgB,divrgV,divrgJ,Pflux,qfactor,lambda,p_tot

        logical :: alt_eom

      end module auxiliaryVariables

c module operators
c ######################################################################
      module operators

        use grid

        use grid_aliases

        use error

        real(8),pointer,dimension(:,:,:,:) :: vec,vec1,vec2

        logical    :: solenoidal=.true.

      contains

c     div
c     ###############################################################
      function div(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az,he,sp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1),div
      integer(4),optional :: he
      logical,optional    :: sp

c     Local variables

      integer(4) :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km
      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacp,jacm,jac0,jach
      logical    :: spoint

c     Begin program

      if (PRESENT(he)) then
        half_elem = he
      else
        half_elem = 0
      endif

      if (PRESENT(sp)) then
        spoint = sp
      else
        spoint = .false.
      endif

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      dxx = 2*dxh(ig)
      dyy = 2*dyh(jg)
      dzz = 2*dzh(kg)

      jac0 = gmetric%grid(igrid)%jac(i,j,k)

      select case(half_elem)
      case(1)
        im = i
        dxx = dx(ig)
        jac0 = 0.5*(gmetric%grid(igrid)%jac(ip,j,k)
     .             +gmetric%grid(igrid)%jac(i ,j,k))
      case(2)
        jm = j
        dyy = dy(jg)
        jac0 = 0.5*(gmetric%grid(igrid)%jac(i,jp,k)
     .             +gmetric%grid(igrid)%jac(i,j ,k))
      case(3)
        km = k
        dzz = dz(kg)
        jac0 = 0.5*(gmetric%grid(igrid)%jac(i,j,kp)
     .             +gmetric%grid(igrid)%jac(i,j,k ))
      end select

      if (isSP(i,j,k,igx,igy,igz).and.half_elem == 0.and.spoint) then
          jacp = gmetric%grid(igrid)%jac(i+1,j,k)
          jach = 0.5*(jacp+jac0)   !Only good for cylindrical-like geom.

          div = ((ax(i+1,j  ,k  )/jacp
     .           +ax(i  ,j  ,k  )/jac0)*jach     )/dxx
     .          +(ay(i  ,j+1,k  )-ay(i  ,j-1,k  ))/dyy
     .          +(az(i  ,j  ,k+1)-az(i  ,j  ,k-1))/dzz
      else
        div =  (ax(ip,j ,k )-ax(im,j ,k ))/dxx
     .        +(ay(i ,jp,k )-ay(i ,jm,k ))/dyy
     .        +(az(i ,j ,kp)-az(i ,j ,km))/dzz
      endif

      div = div/jac0
      
c     End 

      end function div

c     laplacian
c     ###############################################################
      function laplacian(i,j,k,nx,ny,nz,igx,igy,igz,arr,vol)

c     ---------------------------------------------------------------
c     Calculates lap(arr) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian

      logical,optional,intent(IN) :: vol

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

      igrid = igx

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
cc      elseif (.not.sing_point) then
      else
        d_xx_ip = 0.5*(gmetric%grid(igrid)%gsup(i ,j,k,1,1)
     .                +gmetric%grid(igrid)%gsup(ip,j,k,1,1))
        d_xx_im = 0.5*(gmetric%grid(igrid)%gsup(i ,j,k,1,1)
     .                +gmetric%grid(igrid)%gsup(im,j,k,1,1))
        d_yy_jp = 0.5*(gmetric%grid(igrid)%gsup(i,j ,k,2,2)
     .                +gmetric%grid(igrid)%gsup(i,jp,k,2,2))
        d_yy_jm = 0.5*(gmetric%grid(igrid)%gsup(i,j ,k,2,2)
     .                +gmetric%grid(igrid)%gsup(i,jm,k,2,2))
        d_zz_kp = 0.5*(gmetric%grid(igrid)%gsup(i,j,k ,3,3)
     .                +gmetric%grid(igrid)%gsup(i,j,kp,3,3))
        d_zz_km = 0.5*(gmetric%grid(igrid)%gsup(i,j,k ,3,3)
     .                +gmetric%grid(igrid)%gsup(i,j,km,3,3))

        d_xy_ipjp = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(ip,j,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(i,jp,k,1,2))
        d_xy_ipjm = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(ip,j,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(i,jm,k,1,2))
        d_xy_imjp = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(im,j,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(i,jp,k,1,2))
        d_xy_imjm = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(im,j,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(i,j ,k,1,2)
     .                   +gmetric%grid(igrid)%gsup(i,jm,k,1,2))
                 
        d_xz_ipkp = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igrid)%gsup(ip,j,k,1,3)
     .                   +gmetric%grid(igrid)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igrid)%gsup(i,j,kp,1,3))
        d_xz_ipkm = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igrid)%gsup(ip,j,k,1,3)
     .                   +gmetric%grid(igrid)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igrid)%gsup(i,j,km,1,3))
        d_xz_imkp = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igrid)%gsup(ip,j,k,1,3)
     .                   +gmetric%grid(igrid)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igrid)%gsup(i,j,km,1,3))
        d_xz_imkm = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,3)
     .                   +gmetric%grid(igrid)%gsup(im,j,k,1,3)
     .                   +gmetric%grid(igrid)%gsup(i,j,k ,1,3)
     .                   +gmetric%grid(igrid)%gsup(i,j,km,1,3))
                 
cc      else
cc
cc        d_xx_ip = 0.5*(gmetric%grid(igrid)%gsup(i ,j,k,1,1)
cc     .                +gmetric%grid(igrid)%gsup(ip,j,k,1,1))
cc        d_xx_im = 0d0
cc        d_yy_jp = 0.5*(gmetric%grid(igrid)%gsup(i,j ,k,2,2)
cc     .                +gmetric%grid(igrid)%gsup(i,jp,k,2,2))
cc        d_yy_jm = 0.5*(gmetric%grid(igrid)%gsup(i,j ,k,2,2)
cc     .                +gmetric%grid(igrid)%gsup(i,jm,k,2,2))
cc        d_zz_kp = 0.5*(gmetric%grid(igrid)%gsup(i,j,k ,3,3)
cc     .                +gmetric%grid(igrid)%gsup(i,j,kp,3,3))
cc        d_zz_km = 0.5*(gmetric%grid(igrid)%gsup(i,j,k ,3,3)
cc     .                +gmetric%grid(igrid)%gsup(i,j,km,3,3))
cc
cc        d_xy_ipjp = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,2)
cc     .                   +gmetric%grid(igrid)%gsup(ip,j,k,1,2)
cc     .                   +gmetric%grid(igrid)%gsup(i,j ,k,1,2)
cc     .                   +gmetric%grid(igrid)%gsup(i,jp,k,1,2))
cc        d_xy_ipjm = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,2)
cc     .                   +gmetric%grid(igrid)%gsup(ip,j,k,1,2)
cc     .                   +gmetric%grid(igrid)%gsup(i,j ,k,1,2)
cc     .                   +gmetric%grid(igrid)%gsup(i,jm,k,1,2))
cc        d_xy_imjp = 0d0
cc        d_xy_imjm = 0d0
cc
cc        d_xz_ipkp = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,3)
cc     .                   +gmetric%grid(igrid)%gsup(ip,j,k,1,3)
cc     .                   +gmetric%grid(igrid)%gsup(i,j,k ,1,3)
cc     .                   +gmetric%grid(igrid)%gsup(i,j,kp,1,3))
cc        d_xz_ipkm = 0.25*(gmetric%grid(igrid)%gsup(i ,j,k,1,3)
cc     .                   +gmetric%grid(igrid)%gsup(ip,j,k,1,3)
cc     .                   +gmetric%grid(igrid)%gsup(i,j,k ,1,3)
cc     .                   +gmetric%grid(igrid)%gsup(i,j,km,1,3))
cc        d_xz_imkp = 0d0
cc        d_xz_imkm = 0d0

      endif

      d_yz_jpkp = 0.25*(gmetric%grid(igrid)%gsup(i,j ,k,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,jp,k,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,j,k ,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,j,kp,2,3))
      d_yz_jpkm = 0.25*(gmetric%grid(igrid)%gsup(i,j ,k,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,jp,k,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,j,k ,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,j,km,2,3))
      d_yz_jmkp = 0.25*(gmetric%grid(igrid)%gsup(i,j ,k,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,jm,k,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,j,k ,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,j,kp,2,3))
      d_yz_jmkm = 0.25*(gmetric%grid(igrid)%gsup(i,j ,k,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,jm,k,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,j,k ,2,3)
     .                 +gmetric%grid(igrid)%gsup(i,j,km,2,3))

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

      if (.not.vol_wgt) laplacian=laplacian/volume(i,j,k,igx,igy,igz)

c     End program

      end function laplacian

c     div_tensor
c     ###############################################################
      function div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                   ,tsrx,tsry,tsrz,vol) result (divt)

c     ---------------------------------------------------------------
c     Calculates div(tensor) at cell centers in general non-orthogonal
c     coordinates.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8)    :: divt(3)

      integer(4) :: i,j,k,igx,igy,igz,nx,ny,nz

      logical    :: alt_eom

      external   tsrx,tsry,tsrz

      logical,optional,intent(IN) :: vol

c     Local variables

      integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid

      real(8)    :: jac,dvol,dS1,dS2,dS3,dxx,dyy,dzz

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o

      real(8)    :: hess(3,3,3),msource,dx1,dx2,ll

      logical    :: vol_wgt

      real(8)    :: dum1,dum2

c     Begin program

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

      igrid = igx

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      jac = gmetric%grid(igrid)%jac(i,j,k)

cc      hess = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      dS1 = dyy*dzz
      dS2 = dxx*dzz
      dS3 = dxx*dyy

      dvol = dxx*dyy*dzz

      if (coords /= 'car') then
        hess = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
      endif

      call tsrx(i ,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11p,t12p,t13p, 1)
      call tsrx(im,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11m,t12m,t13m,-1)
cc      if (sing_point) then
cc        t11m = 0d0
cc        if (alt_eom) t12m = 0d0
cc        t13m = 0d0
cc      endif
      if (coords /= 'car')
     .    call tsrx(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11o,t12o,t13o,0)

      call tsry(i,j ,k,nx,ny,nz,igx,igy,igz,alt_eom,t21p,t22p,t23p, 1)
      call tsry(i,jm,k,nx,ny,nz,igx,igy,igz,alt_eom,t21m,t22m,t23m,-1)
      if (coords /= 'car')
     .    call tsry(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t21o,t22o,t23o,0)

      call tsrz(i,j,k ,nx,ny,nz,igx,igy,igz,alt_eom,t31p,t32p,t33p, 1)
      call tsrz(i,j,km,nx,ny,nz,igx,igy,igz,alt_eom,t31m,t32m,t33m,-1)
      if (coords /= 'car')
     .    call tsrz(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t31o,t32o,t33o,0)

      msource = 0d0

      !Component 1
      flxip = t11p
      flxim = t11m

      flxjp = t21p
      flxjm = t21m

      flxkp = t31p
      flxkm = t31m

      if (coords /= 'car') then
        msource =  (t11o*hess(1,1,1)+t12o*hess(1,1,2)+t13o*hess(1,1,3)
     .             +t21o*hess(1,2,1)+t22o*hess(1,2,2)+t23o*hess(1,2,3)
     .             +t31o*hess(1,3,1)+t32o*hess(1,3,2)+t33o*hess(1,3,3))
      endif

      divt(1) =  (flxip - flxim)/dxx
     .          +(flxjp - flxjm)/dyy
     .          +(flxkp - flxkm)/dzz + msource/jac

      !Component 2
      flxip = t12p
      flxim = t12m

      flxjp = t22p
      flxjm = t22m

      flxkp = t32p
      flxkm = t32m

      if (coords /= 'car') then
        if (alt_eom) then

          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
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
          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3))
        endif
      endif

      divt(2) =  (flxip - flxim)/dxx
     .          +(flxjp - flxjm)/dyy
     .          +(flxkp - flxkm)/dzz  + msource/jac

      !Component 3
      flxip = t13p
      flxim = t13m

      flxjp = t23p
      flxjm = t23m

      flxkp = t33p
      flxkm = t33m

      if (coords /= 'car') then
        msource = (t11o*hess(3,1,1)+t12o*hess(3,1,2)+t13o*hess(3,1,3)
     .            +t21o*hess(3,2,1)+t22o*hess(3,2,2)+t23o*hess(3,2,3)
     .            +t31o*hess(3,3,1)+t32o*hess(3,3,2)+t33o*hess(3,3,3))
      endif

      divt(3) =  (flxip - flxim)/dxx
     .          +(flxjp - flxjm)/dyy
     .          +(flxkp - flxkm)/dzz + msource/jac

      !Volume factor
      if (vol_wgt) divt=divt*volume(i,j,k,igx,igy,igz)

c     End program

      end function div_tensor

c     veclaplacian
c     ###############################################################
      function veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                     ,alteom,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*lap(vector) at cell centers in general non-orthog.
c     coordinates. Vector is assumed in contravariant representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)

      real(8)    :: vlap(3)

      logical    :: alteom

      logical,optional,intent(IN) :: vol

c     Local variables

      integer(4) :: icomp

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

      if (coords == 'car') then
        do icomp=1,3
          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,vfield(:,:,:,icomp),vol=vol_wgt)
        enddo
        return
      else
        vec => vfield !Pointer passed to nabtensor routines
        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
     $                   ,vol=vol_wgt)
      endif

c     End program

      end function veclaplacian

c     nabtensor_x
c     #############################################################
      subroutine nabtensor_x(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for nabla(vec)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alteom

c     Local variables

        integer(4) :: ig,jg,kg,ip,igrid
        real(8)    :: x,y,z,jac,jac0,jacp
        real(8)    :: nabla_v(3,3),gsuper(3,3)

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

cc        if (bcond(1) == SP .and. i == 0) jac = 1d-10
        if (isSP(i+1,j,k,igx,igy,igz)) jac = 1d-10

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),1)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),0)
        endif

        t11 =( gsuper(1,1)*nabla_v(1,1)
     .        +gsuper(1,2)*nabla_v(2,1)
     .        +gsuper(1,3)*nabla_v(3,1) )

        t12 =( gsuper(1,1)*nabla_v(1,2)
     .        +gsuper(1,2)*nabla_v(2,2)
     .        +gsuper(1,3)*nabla_v(3,2) )

        t13 =( gsuper(1,1)*nabla_v(1,3)
     .        +gsuper(1,2)*nabla_v(2,3)
     .        +gsuper(1,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t11 = t11/jac
          if (.not.alteom) t12 = t12/jac
          t13 = t13/jac
        endif

c     End program

      end subroutine nabtensor_x

c     nabtensor_y
c     #############################################################
      subroutine nabtensor_y(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alteom

c     Local variables

        integer(4) :: ig,jg,kg,jp,igrid
        real(8)    :: x,y,z,jac
        real(8)    :: nabla_v(3,3),gsuper(3,3)

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),2)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),0)
        endif

        t21 =( gsuper(2,1)*nabla_v(1,1)
     .        +gsuper(2,2)*nabla_v(2,1)
     .        +gsuper(2,3)*nabla_v(3,1) )

        t22 =( gsuper(2,1)*nabla_v(1,2)
     .        +gsuper(2,2)*nabla_v(2,2)
     .        +gsuper(2,3)*nabla_v(3,2) )

        t23 =( gsuper(2,1)*nabla_v(1,3)
     .        +gsuper(2,2)*nabla_v(2,3)
     .        +gsuper(2,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t21 = t21/jac
          if (.not.alteom) t22 = t22/jac
          t23 = t23/jac
        endif

c     End program

      end subroutine nabtensor_y

c     nabtensor_z
c     #############################################################
      subroutine nabtensor_z(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alteom

c     Local variables

        integer(4) :: ig,jg,kg,kp,igrid
        real(8)    :: x,y,z,jac
        real(8)    :: nabla_v(3,3),gsuper(3,3)

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),3)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),0)
        endif

        t31 =( gsuper(3,1)*nabla_v(1,1)
     .        +gsuper(3,2)*nabla_v(2,1)
     .        +gsuper(3,3)*nabla_v(3,1) )

        t32 =( gsuper(3,1)*nabla_v(1,2)
     .        +gsuper(3,2)*nabla_v(2,2)
     .        +gsuper(3,3)*nabla_v(3,2) )

        t33 =( gsuper(3,1)*nabla_v(1,3)
     .        +gsuper(3,2)*nabla_v(2,3)
     .        +gsuper(3,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alteom) t32 = t32/jac
          t33 = t33/jac
        endif

c     End program

      end subroutine nabtensor_z

c     curl
c     ###############################################################
      function curl(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az) result(crl)

c     ---------------------------------------------------------------
c     Calculates curl(A) at cell centers in general non-orthogonal
c     coordinates. The vector components (ax,ay,az) are covariant.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8)    :: crl(3)

      integer(4) :: i,j,k,comp,nx,ny,nz,igx,igy,igz

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid
      real(8)    :: dx,dy,dz,flxip,flxim,flxjp,flxjm,flxkp,flxkm
      real(8)    :: dS1,dS2,dS3,x0,y0,z0,jac,dx1,dx2,ll
     .             ,xip,yip,zip,jacp,xh,yh,zh,jach

c     Begin program

      igrid = igx

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

      !X comp

      flxjp = 0.5*(az(i,jp,k)+az(i,j,k))
      flxjm = 0.5*(az(i,jm,k)+az(i,j,k))

      flxkp =-0.5*(ay(i,j,kp)+ay(i,j,k))
      flxkm =-0.5*(ay(i,j,km)+ay(i,j,k))

      crl(1) = (flxjp-flxjm)/dy
     .        +(flxkp-flxkm)/dz

      !Y comp

cc      if (sing_point) then
cc        !Second order formula
cccc        dx1=grid_params%dx(ig-1)
cccc        dx2=grid_params%dx(ig  )
cccc        dx = 0.5*(dx1+dx2)
cccc        ll = dx1/dx2
cccc        flxip =-0.5*(az(ip,j,k)*ll + az(i,j,k)*(1.+1./ll))
cccc        flxim =-0.5*(az(im,j,k)/ll + az(i,j,k)*(1.+   ll))
cc        flxip = -0.5*(az(ip,j,k)+az(i,j,k))
cc        flxim = -az(im,j,k)
cc      else
        flxip =-0.5*(az(ip,j,k)+az(i,j,k))
        flxim =-0.5*(az(im,j,k)+az(i,j,k))
cc      endif

      flxkp = 0.5*(ax(i,j,kp)+ax(i,j,k))
      flxkm = 0.5*(ax(i,j,km)+ax(i,j,k))

      crl(2) = (flxip-flxim)/dx
     .        +(flxkp-flxkm)/dz

      !Z comp

cc      if (sing_point) then
cccc        dx1=grid_params%dx(ig-1)
cccc        dx2=grid_params%dx(ig  )
cccc        dx = 0.5*(dx1+dx2)
cccc        ll = dx1/dx2
cccc        flxip = 0.5*(ay(ip,j,k)*ll + ay(i,j,k)*(1.+1./ll))
cccc        flxim = 0.5*(ay(im,j,k)/ll + ay(i,j,k)*(1.+   ll))
cc        jacp = gmetric%grid(igrid)%jac(ip,j,k)
cc        jac  = gmetric%grid(igrid)%jac(i ,j,k)
cc        jach = 0.5*(jacp+jac)
cc
cc        flxip = 0.5*(ay(ip,j,k)/jacp+ay(i,j,k)/jac)*jach
cc        flxim = ay(im,j,k)
cc      else
        flxip = 0.5*(ay(ip,j,k)+ay(i,j,k))
        flxim = 0.5*(ay(im,j,k)+ay(i,j,k))
cc      endif

      flxjp =-0.5*(ax(i,jp,k)+ax(i,j,k))
      flxjm =-0.5*(ax(i,jm,k)+ax(i,j,k))

      crl(3) = (flxip-flxim)/dx
     .        +(flxjp-flxjm)/dy

c     End program

      end function curl

c     curl2
c     ###############################################################
      function curl2(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az,comp)

c     ---------------------------------------------------------------
c     Calculates curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8)    :: curl2

      integer(4) :: i,j,k,comp,nx,ny,nz,igx,igy,igz

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid
      real(8)    :: dx,dy,dz,dx1,dx2,dy1,dy2,dz1,dz2
      real(8)    :: daxdz,dazdx,daydx,daxdy,dazdy,daydz
      real(8)    :: jac0,jacp,jacm,ll,al,ar

c     Begin program

      igrid = igx

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
          dazdy = 0.5*(az(i,jp,k)-az(i,jm,k))/dy
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
          daydz = 0.5*(ay(i,j,kp)-ay(i,j,km))/dz
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
cc        elseif (sing_point) then
cc          !Second order formula
cc          dx1=grid_params%dx(ig-1)
cc          dx2=grid_params%dx(ig  )
cc          ll = dx1/dx2
cc          ar = 0.5*(az(ip,j,k)*ll + az(i,j,k)*(1.+1./ll))
cc          al = 0.5*(az(im,j,k)/ll + az(i,j,k)*(1.+   ll))
cc          dazdx = 2.*(ar-al)/(dx1+dx2)
cccc          dx = grid_params%dxh(ig)
cccc          dazdx = (0.5*(az(ip,j,k)+az(i,j,k))-az(im,j,k))/dx
        else
          dx = grid_params%dxh(ig)
          dazdx = 0.5*(az(ip,j,k)-az(im,j,k))/dx
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
          daxdz = 0.5*(ax(i,j,kp)-ax(i,j,km))/dz
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
cc        elseif (sing_point) then
cc          dx = grid_params%dxh(ig)
cc
cc          jacp = gmetric%grid(igrid)%jac(ip,j,k)
cc          jac0 = gmetric%grid(igrid)%jac(i ,j,k)
cc          jacm = 0.5*(jacp+jac0)
cc
cc          daydx = ((ay(ip,j,k)/jacp+ay(i,j,k)/jac0)/2.*jacm
cc     .            - ay(im,j,k))/dx
cc
cccc          dx1=grid_params%dx(ig-1)
cccc          dx2=grid_params%dx(ig  )
cccc          ll = dx1/dx2
cccc          ar = 0.5*(ay(ip,j,k)*ll + ay(i,j,k)*(1.+1./ll))
cccc          al = 0.5*(ay(im,j,k)/ll + ay(i,j,k)*(1.+   ll))
cccc          daydx = 2.*(ar-al)/(dx1+dx2)
        else
          dx = grid_params%dxh(ig)
          daydx = 0.5*(ay(ip,j,k)-ay(im,j,k))/dx
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
          daxdy = 0.5*(ax(i,jp,k)-ax(i,jm,k))/dy
        endif

        curl2 = daydx - daxdy

      case default

        messg = 'Error in component in curl'
        call pstop('curl2',messg)

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

ccc     fnabla_v
ccc     #############################################################
cc      function fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az,half_elem)
cc     .         result(tensor)
ccc     -------------------------------------------------------------
ccc     Calculates the tensor nabla(vec v) at the following positions:
ccc       + half_elem=0 --> i,j,k
ccc       + half_elem=1 --> i+1/2,j,k
ccc       + half_elem=2 --> i,j+1/2,k
ccc       + half_elem=3 --> i,j,k+1/2
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
cc        real(8)    :: tensor(3,3)
cc        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
cc     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
cc     .               ,az(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
cc        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
cc        real(8)    :: vxx,vyy,vzz
cc     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
cc     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
cc     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
cc        real(8)    :: tsrc(3,3)
cc        logical    :: sing_point
cc
ccc     Begin program
cc
cc        igrid = igx
cc
cc        sing_point = isSP(i,j,k,igx,igy,igz)
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
cc          vxip = ax(ip,j,k)
cc          vxim = ax(i ,j,k)
cc          vyip = ay(ip,j,k)
cc          vyim = ay(i ,j,k)
cc          vzip = az(ip,j,k)
cc          vzim = az(i ,j,k)
cc
cc          vxjp = (ax(ip,jp,k)+ax(i,jp,k))/2.
cc          vxjm = (ax(ip,jm,k)+ax(i,jm,k))/2.
cc          vyjp = (ay(ip,jp,k)+ay(i,jp,k))/2.
cc          vyjm = (ay(ip,jm,k)+ay(i,jm,k))/2.
cc          vzjp = (az(ip,jp,k)+az(i,jp,k))/2.
cc          vzjm = (az(ip,jm,k)+az(i,jm,k))/2.
cc
cc          vxkp = (ax(ip,j,kp)+ax(i,j,kp))/2.
cc          vxkm = (ax(ip,j,km)+ax(i,j,km))/2.
cc          vykp = (ay(ip,j,kp)+ay(i,j,kp))/2.
cc          vykm = (ay(ip,j,km)+ay(i,j,km))/2.
cc          vzkp = (az(ip,j,kp)+az(i,j,kp))/2.
cc          vzkm = (az(ip,j,km)+az(i,j,km))/2.
cc
cc          tsrc = 0.5*(nabla_v_src(ip,j,k)
cc     .               +nabla_v_src(i ,j,k))
cc
cc        case (2)
cc          dhy = dy(jg)
cc
cc          if (sing_point) then
cc            vxip = (ax(ip,j,k)+ax(ip,jp,k))/2.
cc     .            +(ax(i ,j,k)+ax(i ,jp,k))/2.
cc            vxim = (ax(im,j,k)+ax(im,jp,k))
cc            vyip = (ay(ip,j,k)+ay(ip,jp,k))/2.
cc     .            +(ay(i ,j,k)+ay(i ,jp,k))/2.
cc            vyim = (ay(im,j,k)+ay(im,jp,k))
cc            vzip = (az(ip,j,k)+az(ip,jp,k))/2
cc     .            +(az(i ,j,k)+az(i ,jp,k))/2.
cc            vzim = (az(im,j,k)+az(im,jp,k))
cc          else
cc            vxip = (ax(ip,j,k)+ax(ip,jp,k))/2.
cc            vxim = (ax(im,j,k)+ax(im,jp,k))/2.
cc            vyip = (ay(ip,j,k)+ay(ip,jp,k))/2.
cc            vyim = (ay(im,j,k)+ay(im,jp,k))/2.
cc            vzip = (az(ip,j,k)+az(ip,jp,k))/2.
cc            vzim = (az(im,j,k)+az(im,jp,k))/2.
cc          endif
cc
cc          vxjp = ax(i,jp,k)
cc          vxjm = ax(i,j ,k)
cc          vyjp = ay(i,jp,k)
cc          vyjm = ay(i,j ,k)
cc          vzjp = az(i,jp,k)
cc          vzjm = az(i,j ,k)
cc
cc          vxkp = (ax(i,j,kp)+ax(i,jp,kp))/2.
cc          vxkm = (ax(i,j,km)+ax(i,jp,km))/2.
cc          vykp = (ay(i,j,kp)+ay(i,jp,kp))/2.
cc          vykm = (ay(i,j,km)+ay(i,jp,km))/2.
cc          vzkp = (az(i,j,kp)+az(i,jp,kp))/2.
cc          vzkm = (az(i,j,km)+az(i,jp,km))/2.
cc
cc          tsrc = 0.5*(nabla_v_src(i,jp,k)
cc     .               +nabla_v_src(i,j ,k))
cc
cc        case (3)
cc          dhz = dz(kg)
cc
cc          if (sing_point) then
cc            vxip = (ax(ip,j,k)+ax(ip,j,kp))/2.
cc     .            +(ax(i ,j,k)+ax(i ,j,kp))/2.
cc            vxim = (ax(im,j,k)+ax(im,j,kp))
cc            vyip = (ay(ip,j,k)+ay(ip,j,kp))/2.
cc     .            +(ay(i ,j,k)+ay(i ,j,kp))/2.
cc            vyim = (ay(im,j,k)+ay(im,j,kp))
cc            vzip = (az(ip,j,k)+az(ip,j,kp))/2.
cc     .            +(az(i ,j,k)+az(i ,j,kp))/2.
cc            vzim = (az(im,j,k)+az(im,j,kp))
cc          else
cc            vxip = (ax(ip,j,k)+ax(ip,j,kp))/2.
cc            vxim = (ax(im,j,k)+ax(im,j,kp))/2.
cc            vyip = (ay(ip,j,k)+ay(ip,j,kp))/2.
cc            vyim = (ay(im,j,k)+ay(im,j,kp))/2.
cc            vzip = (az(ip,j,k)+az(ip,j,kp))/2.
cc            vzim = (az(im,j,k)+az(im,j,kp))/2.
cc          endif
cc
cc          vxjp = (ax(i,jp,k)+ax(i,jp,kp))/2.
cc          vxjm = (ax(i,jm,k)+ax(i,jm,kp))/2.
cc          vyjp = (ay(i,jp,k)+ay(i,jp,kp))/2.
cc          vyjm = (ay(i,jm,k)+ay(i,jm,kp))/2.
cc          vzjp = (az(i,jp,k)+az(i,jp,kp))/2.
cc          vzjm = (az(i,jm,k)+az(i,jm,kp))/2.
cc
cc          vxkp = ax(i,j,kp)
cc          vxkm = ax(i,j,k )
cc          vykp = ay(i,j,kp)
cc          vykm = ay(i,j,k )
cc          vzkp = az(i,j,kp)
cc          vzkm = az(i,j,k )
cc
cc          tsrc = 0.5*(nabla_v_src(i,j,kp)
cc     .               +nabla_v_src(i,j,k ))
cc
cc        case default
cc
cc          if (sing_point) then
cc            vxip = ax(ip,j,k)+ax(i,j,k)
cc            vxim = 2.*ax(im,j,k)
cc            vyip = ay(ip,j,k)+ay(i,j,k)
cc            vyim = 2.*ay(im,j,k)
cc            vzip = az(ip,j,k)+az(i,j,k)
cc            vzim = 2.*az(im,j,k)
cc          else
cc            vxip = ax(ip,j,k)
cc            vxim = ax(im,j,k)
cc            vyip = ay(ip,j,k)
cc            vyim = ay(im,j,k)
cc            vzip = az(ip,j,k)
cc            vzim = az(im,j,k)
cc          endif
cc
cc          vxjp = ax(i,jp,k)
cc          vxjm = ax(i,jm,k)
cc          vyjp = ay(i,jp,k)
cc          vyjm = ay(i,jm,k)
cc          vzjp = az(i,jp,k)
cc          vzjm = az(i,jm,k)
cc
cc          vxkp = ax(i,j,kp)
cc          vxkm = ax(i,j,km)
cc          vykp = ay(i,j,kp)
cc          vykm = ay(i,j,km)
cc          vzkp = az(i,j,kp)
cc          vzkm = az(i,j,km)
cc
cc          tsrc = nabla_v_src(i,j,k)
cc
cc        end select
cc
cc        idhx = 1./dhx
cc        idhy = 1./dhy
cc        idhz = 1./dhz
cc
cc      ! l = 1, m = 1
cc        tensor(1,1) = (vxip-vxim)*idhx
cc
cc      ! l = 1, m = 2
cc        tensor(1,2) = (vyip-vyim)*idhx
cc
cc      ! l = 1, m = 3
cc        tensor(1,3) = (vzip-vzim)*idhx
cc
cc      ! l = 2, m = 1
cc        tensor(2,1) = (vxjp-vxjm)*idhy
cc
cc      ! l = 2, m = 2
cc        tensor(2,2) = (vyjp-vyjm)*idhy
cc
cc      ! l = 2, m = 3
cc        tensor(2,3) = (vzjp-vzjm)*idhy
cc
cc      ! l = 3, m = 1
cc        tensor(3,1) = (vxkp-vxkm)*idhz
cc
cc      ! l = 3, m = 2
cc        tensor(3,2) = (vykp-vykm)*idhz
cc
cc      ! l = 3, m = 3
cc        tensor(3,3) = (vzkp-vzkm)*idhz
cc
cc      ! Add geometric source
cc
cc        tensor = tensor - tsrc
cc
ccc     End program
cc
cc      contains
cc
ccc     nabla_v_src
ccc     #############################################################
cc      function nabla_v_src(i,j,k) result(tensor)
cc
ccc     -------------------------------------------------------------
ccc     Finds geometric source of tensor nabla(v) at cell (i,j,k)
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
cc        integer(4) :: i,j,k
cc        real(8)    :: tensor(3,3)
cc
cc        real(8)    :: vxx,vyy,vzz,hessian(3,3,3)
cc
ccc     Begin program
cc
cccc        vxx = ax(i,j,k)
cccc        vyy = ay(i,j,k)
cccc        vzz = az(i,j,k)
cc
cc        hessian = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc
cccc      ! l = 1, m = 1
cccc        tensor(1,1) =  vxx*(hessian(1,1,1)
cccc     .                    + hessian(2,2,1)
cccc     .                    + hessian(3,3,1))
cccc     .               - vxx *hessian(1,1,1)
cccc     .               - vyy *hessian(1,2,1)
cccc     .               - vzz *hessian(1,3,1)
cccc
cccc      ! l = 1, m = 2
cccc        tensor(1,2) =  vyy*(hessian(1,1,1)
cccc     .                    + hessian(2,2,1)
cccc     .                    + hessian(3,3,1))
cccc     .               - vxx* hessian(2,1,1)
cccc     .               - vyy* hessian(2,2,1)
cccc     .               - vzz* hessian(2,3,1)
cccc
cccc      ! l = 1, m = 3
cccc        tensor(1,3) =  vzz*(hessian(1,1,1)
cccc     .                    + hessian(2,2,1)
cccc     .                    + hessian(3,3,1))
cccc     .               - vxx* hessian(3,1,1)
cccc     .               - vyy* hessian(3,2,1)
cccc     .               - vzz* hessian(3,3,1)
cccc
cccc      ! l = 2, m = 1
cccc        tensor(2,1) =  vxx*(hessian(1,1,2)
cccc     .                    + hessian(2,2,2)
cccc     .                    + hessian(3,3,2))
cccc     .               - vxx* hessian(1,1,2)
cccc     .               - vyy* hessian(1,2,2)
cccc     .               - vzz* hessian(1,3,2)
cccc
cccc      ! l = 2, m = 2
cccc        tensor(2,2) =  vyy*(hessian(1,1,2)
cccc     .                    + hessian(2,2,2)
cccc     .                    + hessian(3,3,2))
cccc     .               - vxx* hessian(2,1,2)
cccc     .               - vyy* hessian(2,2,2)
cccc     .               - vzz* hessian(2,3,2)
cccc
cccc      ! l = 2, m = 3
cccc        tensor(2,3) =  vzz*(hessian(1,1,2)
cccc     .                    + hessian(2,2,2)
cccc     .                    + hessian(3,3,2))
cccc     .               - vxx* hessian(3,1,2)
cccc     .               - vyy* hessian(3,2,2)
cccc     .               - vzz* hessian(3,3,2)
cccc
cccc      ! l = 3, m = 1
cccc        tensor(3,1) =  vxx*(hessian(1,1,3)
cccc     .                    + hessian(2,2,3)
cccc     .                    + hessian(3,3,3))
cccc     .               - vxx* hessian(1,1,3)
cccc     .               - vyy* hessian(1,2,3)
cccc     .               - vzz* hessian(1,3,3)
cccc
cccc      ! l = 3, m = 2
cccc        tensor(3,2) =  vyy*(hessian(1,1,3)
cccc     .                    + hessian(2,2,3)
cccc     .                    + hessian(3,3,3))
cccc     .               - vxx* hessian(2,1,3)
cccc     .               - vyy* hessian(2,2,3)
cccc     .               - vzz* hessian(2,3,3)
cccc
cccc      ! l = 3, m = 3
cccc        tensor(3,3) =  vzz*(hessian(1,1,3)
cccc     .                    + hessian(2,2,3)
cccc     .                    + hessian(3,3,3))
cccc     .               - vxx* hessian(3,1,3)
cccc     .               - vyy* hessian(3,2,3)
cccc     .               - vzz* hessian(3,3,3)
cc
cc      ! l = 1, m = 1
cc        tensor(1,1) =  ax(i,j,k)*(hessian(2,2,1)
cc     .                          + hessian(3,3,1))
cc     .               - ay(i,j,k) *hessian(1,2,1)
cc     .               - az(i,j,k) *hessian(1,3,1)
cc
cc      ! l = 1, m = 2
cc        tensor(1,2) =  ay(i,j,k)*(hessian(1,1,1)
cc     .                          + hessian(3,3,1))
cc     .               - ax(i,j,k)* hessian(2,1,1)
cc     .               - az(i,j,k)* hessian(2,3,1)
cc
cc      ! l = 1, m = 3
cc        tensor(1,3) =  az(i,j,k)*(hessian(1,1,1)
cc     .                          + hessian(2,2,1))
cc     .               - ax(i,j,k)* hessian(3,1,1)
cc     .               - ay(i,j,k)* hessian(3,2,1)
cc
cc      ! l = 2, m = 1
cc        tensor(2,1) =  ax(i,j,k)*(hessian(2,2,2)
cc     .                          + hessian(3,3,2))
cc     .               - ay(i,j,k)* hessian(1,2,2)
cc     .               - az(i,j,k)* hessian(1,3,2)
cc
cc      ! l = 2, m = 2
cc        tensor(2,2) =  ay(i,j,k)*(hessian(1,1,2)
cc     .                          + hessian(3,3,2))
cc     .               - ax(i,j,k)* hessian(2,1,2)
cc     .               - az(i,j,k)* hessian(2,3,2)
cc
cc      ! l = 2, m = 3
cc        tensor(2,3) =  az(i,j,k)*(hessian(1,1,2)
cc     .                          + hessian(2,2,2))
cc     .               - ax(i,j,k)* hessian(3,1,2)
cc     .               - ay(i,j,k)* hessian(3,2,2)
cc
cc      ! l = 3, m = 1
cc        tensor(3,1) =  ax(i,j,k)*(hessian(2,2,3)
cc     .                          + hessian(3,3,3))
cc     .               - ay(i,j,k)* hessian(1,2,3)
cc     .               - az(i,j,k)* hessian(1,3,3)
cc
cc      ! l = 3, m = 2
cc        tensor(3,2) =  ay(i,j,k)*(hessian(1,1,3)
cc     .                          + hessian(3,3,3))
cc     .               - ax(i,j,k)* hessian(2,1,3)
cc     .               - az(i,j,k)* hessian(2,3,3)
cc
cc      ! l = 3, m = 3
cc        tensor(3,3) =  az(i,j,k)*(hessian(1,1,3)
cc     .                          + hessian(2,2,3))
cc     .               - ax(i,j,k)* hessian(3,1,3)
cc     .               - ay(i,j,k)* hessian(3,2,3)
cc
ccc     End program
cc
cc      end function nabla_v_src
cc
cc      end function fnabla_v

c     fnabla_v
c     #############################################################
      function fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az,half_elem)
     .         result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor components T_l^m of nabla(vec v) and
c     fills tensor(l,m) at the following positions:
c       + half_elem=0 --> i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: tsrc(3,3),jacip,jacjp,jacjm,jackp,jackm
     .               ,jacipkp,jacipkm,jacipjp,jacipjm,jac

c     Begin program

        igrid = igx

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

cc          if (isSP(ip,j,k,igrid,igrid,igrid)) ip = i

cc          if (sing_point) then
cc            jacip  = gmetric%grid(igrid)%jac(ip,j,k )
cc            jacipjp= gmetric%grid(igrid)%jac(ip,jp,k)
cc            jacipjm= gmetric%grid(igrid)%jac(ip,jm,k)
cc            jacjp  = gmetric%grid(igrid)%jac(i ,jp,k)
cc            jacjm  = gmetric%grid(igrid)%jac(i ,jm,k)
cc            jacipkp= gmetric%grid(igrid)%jac(ip,j,kp)
cc            jacipkm= gmetric%grid(igrid)%jac(ip,j,km)
cc            jackp  = gmetric%grid(igrid)%jac(i ,j,kp)
cc            jackm  = gmetric%grid(igrid)%jac(i ,j,km)
cc            jac    = gmetric%grid(igrid)%jac(i ,j,k )
cc
cc            vxjp = 0.5*(ax(ip,jp,k)/jacipjp+ax(i,jp,k)/jacjp)
cc     .            *0.5*(jacipjp+jacjp)
cc            vxjm = 0.5*(ax(ip,jm,k)/jacipjm+ax(i,jm,k)/jacjm)
cc     .            *0.5*(jacipjm+jacjm)
cc            vyjp = 0.5*(ay(ip,jp,k)+ay(i,jp,k))
cc            vyjm = 0.5*(ay(ip,jm,k)+ay(i,jm,k))
cc            vzjp = 0.5*(az(ip,jp,k)/jacipjp+az(i,jp,k)/jacjp)
cc     .            *0.5*(jacipjp+jacjp)
cc            vzjm = 0.5*(az(ip,jm,k)/jacipjm+az(i,jm,k)/jacjm)
cc     .            *0.5*(jacipjm+jacjm)
cc
cc            vxkp = 0.5*(ax(ip,j,kp)/jacipkp+ax(i,j,kp)/jackp)
cc     .            *0.5*(jacipkp+jackp)
cc            vxkm = 0.5*(ax(ip,j,km)/jacipkm+ax(i,j,km)/jackm)
cc     .            *0.5*(jacipkm+jackm)
cc            vykp = 0.5*(ay(ip,j,kp)+ay(i,j,kp))
cc            vykm = 0.5*(ay(ip,j,km)+ay(i,j,km))
cc            vzkp = 0.5*(az(ip,j,kp)/jacipkp+az(i,j,kp)/jackp)
cc     .            *0.5*(jacipkp+jackp)
cc            vzkm = 0.5*(az(ip,j,km)/jacipkm+az(i,j,km)/jackm)
cc     .            *0.5*(jacipkm+jackm)
cc          else
            vxjp = 0.5*(ax(ip,jp,k)+ax(i,jp,k))
            vxjm = 0.5*(ax(ip,jm,k)+ax(i,jm,k))
            vyjp = 0.5*(ay(ip,jp,k)+ay(i,jp,k))
            vyjm = 0.5*(ay(ip,jm,k)+ay(i,jm,k))
            vzjp = 0.5*(az(ip,jp,k)+az(i,jp,k))
            vzjm = 0.5*(az(ip,jm,k)+az(i,jm,k))

            vxkp = 0.5*(ax(ip,j,kp)+ax(i,j,kp))
            vxkm = 0.5*(ax(ip,j,km)+ax(i,j,km))
            vykp = 0.5*(ay(ip,j,kp)+ay(i,j,kp))
            vykm = 0.5*(ay(ip,j,km)+ay(i,j,km))
            vzkp = 0.5*(az(ip,j,kp)+az(i,j,kp))
            vzkm = 0.5*(az(ip,j,km)+az(i,j,km))
cc          endif

          tsrc = 0.5*(nabla_v_src(ip,j,k)
     .               +nabla_v_src(i ,j,k))

        case (2)
          dhy = dy(jg)

cc          if (sing_point) then
cc            dhx  = dxh(ig)
cc
cc            jacip  = gmetric%grid(igrid)%jac(ip,j,k )
cc            jacipjp= gmetric%grid(igrid)%jac(ip,jp,k)
cc            jacjp  = gmetric%grid(igrid)%jac(i ,jp,k)
cc            jac    = gmetric%grid(igrid)%jac(i,j,k )
cc
cc            vxip = 0.25*(ax(ip,j,k)/jacip+ax(ip,jp,k)/jacipjp
cc     .                  +ax(i ,j,k)/jac  +ax(i ,jp,k)/jacjp  )
cc     .            *0.25*(jacip+jacipjp+jac+jacjp)
cc            vxim = 0.5*(ax(im,j,k)+ax(im,jp,k))
cc            vyip = 0.25*(ay(ip,j,k)+ay(ip,jp,k)
cc     .                  +ay(i ,j,k)+ay(i ,jp,k))
cc            vyim = 0.5*(ay(im,j,k)+ay(im,jp,k))
cc            vzip = 0.25*(az(ip,j,k)/jacip+az(ip,jp,k)/jacipjp
cc     .                  +az(i ,j,k)/jac  +az(i ,jp,k)/jacjp  )
cc     .            *0.25*(jacip+jacipjp+jac+jacjp)
cc            vzim = 0.5*(az(im,j,k)+az(im,jp,k))
cc          else
            vxip = 0.5*(ax(ip,j,k)+ax(ip,jp,k))
            vxim = 0.5*(ax(im,j,k)+ax(im,jp,k))
            vyip = 0.5*(ay(ip,j,k)+ay(ip,jp,k))
            vyim = 0.5*(ay(im,j,k)+ay(im,jp,k))
            vzip = 0.5*(az(ip,j,k)+az(ip,jp,k))
            vzim = 0.5*(az(im,j,k)+az(im,jp,k))
cc          endif

          vxjp = ax(i,jp,k)
          vxjm = ax(i,j ,k)
          vyjp = ay(i,jp,k)
          vyjm = ay(i,j ,k)
          vzjp = az(i,jp,k)
          vzjm = az(i,j ,k)

          vxkp = 0.5*(ax(i,j,kp)+ax(i,jp,kp))
          vxkm = 0.5*(ax(i,j,km)+ax(i,jp,km))
          vykp = 0.5*(ay(i,j,kp)+ay(i,jp,kp))
          vykm = 0.5*(ay(i,j,km)+ay(i,jp,km))
          vzkp = 0.5*(az(i,j,kp)+az(i,jp,kp))
          vzkm = 0.5*(az(i,j,km)+az(i,jp,km))

          tsrc = 0.5*(nabla_v_src(i,jp,k)
     .               +nabla_v_src(i,j ,k))

        case (3)
          dhz = dz(kg)

cc          if (sing_point) then
cc            dhx  = dxh(ig)
cc
cc            jacip  = gmetric%grid(igrid)%jac(ip,j,k )
cc            jacipkp= gmetric%grid(igrid)%jac(ip,j,kp)
cc            jackp  = gmetric%grid(igrid)%jac(i ,j,kp)
cc            jac    = gmetric%grid(igrid)%jac(i,j,k )
cc
cc            vxip = 0.25*(ax(ip,j,k)/jacip+ax(ip,j,kp)/jacipkp
cc     .                 + ax(i ,j,k)/jac  +ax(i ,j,kp)/jackp  )
cc     .            *0.25*(jacip+jacipkp+jackp+jac)
cc            vxim = 0.5*(ax(im,j,k)+ax(im,j,kp))
cc            vyip = 0.25*(ay(ip,j,k)+ay(ip,j,kp)
cc     .                  +ay(i ,j,k)+ay(i ,j,kp))
cc            vyim = 0.5*(ay(im,j,k)+ay(im,j,kp))
cc            vzip = 0.25*(az(ip,j,k)/jacip+az(ip,j,kp)/jacipkp
cc     .                  +az(i ,j,k)/jac  +az(i ,j,kp)/jackp  )
cc     .            *0.25*(jacip+jacipkp+jackp+jac)
cc            vzim = 0.5*(az(im,j,k)+az(im,j,kp))
cc          else
            vxip = 0.5*(ax(ip,j,k)+ax(ip,j,kp))
            vxim = 0.5*(ax(im,j,k)+ax(im,j,kp))
            vyip = 0.5*(ay(ip,j,k)+ay(ip,j,kp))
            vyim = 0.5*(ay(im,j,k)+ay(im,j,kp))
            vzip = 0.5*(az(ip,j,k)+az(ip,j,kp))
            vzim = 0.5*(az(im,j,k)+az(im,j,kp))
cc          endif

          vxjp = 0.5*(ax(i,jp,k)+ax(i,jp,kp))
          vxjm = 0.5*(ax(i,jm,k)+ax(i,jm,kp))
          vyjp = 0.5*(ay(i,jp,k)+ay(i,jp,kp))
          vyjm = 0.5*(ay(i,jm,k)+ay(i,jm,kp))
          vzjp = 0.5*(az(i,jp,k)+az(i,jp,kp))
          vzjm = 0.5*(az(i,jm,k)+az(i,jm,kp))

          vxkp = ax(i,j,kp)
          vxkm = ax(i,j,k )
          vykp = ay(i,j,kp)
          vykm = ay(i,j,k )
          vzkp = az(i,j,kp)
          vzkm = az(i,j,k )

          tsrc = 0.5*(nabla_v_src(i,j,kp)
     .               +nabla_v_src(i,j,k ))

        case default

cc          if (sing_point) then
cc            dhx  = dxh(ig)
cc
cc            jacip= gmetric%grid(igrid)%jac(ip,j,k)
cc            jac  = gmetric%grid(igrid)%jac(i,j,k )
cc
cc            vxip = 0.5*(ax(ip,j,k)/jacip+ax(i,j,k)/jac)
cc     .            *0.5*(jacip+jac)
cc            vxim = ax(im,j,k)
cc            vyip = 0.5*(ay(ip,j,k)+ay(i,j,k))
cc            vyim = ay(im,j,k)
cc            vzip = 0.5*(az(ip,j,k)/jacip+az(i,j,k)/jac)
cc     .            *0.5*(jacip+jac)
cc            vzim = az(im,j,k)
cc          else
            vxip = ax(ip,j,k)
            vxim = ax(im,j,k)
            vyip = ay(ip,j,k)
            vyim = ay(im,j,k)
            vzip = az(ip,j,k)
            vzim = az(im,j,k)
cc          endif

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

c     Begin program

        hessian = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

      ! l = 1, m = 1
        tensor(1,1) =  ax(i,j,k)*(hessian(2,2,1)
     .                          + hessian(3,3,1))
     .               - ay(i,j,k)* hessian(1,2,1)
     .               - az(i,j,k)* hessian(1,3,1)

      ! l = 1, m = 2
        tensor(1,2) =  ay(i,j,k)*(hessian(1,1,1)
     .                          + hessian(3,3,1))
     .               - ax(i,j,k)* hessian(2,1,1)
     .               - az(i,j,k)* hessian(2,3,1)

      ! l = 1, m = 3
        tensor(1,3) =  az(i,j,k)*(hessian(1,1,1)
     .                          + hessian(2,2,1))
     .               - ax(i,j,k)* hessian(3,1,1)
     .               - ay(i,j,k)* hessian(3,2,1)

      ! l = 2, m = 1
        tensor(2,1) =  ax(i,j,k)*(hessian(2,2,2)
     .                          + hessian(3,3,2))
     .               - ay(i,j,k)* hessian(1,2,2)
     .               - az(i,j,k)* hessian(1,3,2)

      ! l = 2, m = 2
        tensor(2,2) =  ay(i,j,k)*(hessian(1,1,2)
     .                          + hessian(3,3,2))
     .               - ax(i,j,k)* hessian(2,1,2)
     .               - az(i,j,k)* hessian(2,3,2)

      ! l = 2, m = 3
        tensor(2,3) =  az(i,j,k)*(hessian(1,1,2)
     .                          + hessian(2,2,2))
     .               - ax(i,j,k)* hessian(3,1,2)
     .               - ay(i,j,k)* hessian(3,2,2)

      ! l = 3, m = 1
        tensor(3,1) =  ax(i,j,k)*(hessian(2,2,3)
     .                          + hessian(3,3,3))
     .               - ay(i,j,k)* hessian(1,2,3)
     .               - az(i,j,k)* hessian(1,3,3)

      ! l = 3, m = 2
        tensor(3,2) =  ay(i,j,k)*(hessian(1,1,3)
     .                          + hessian(3,3,3))
     .               - ax(i,j,k)* hessian(2,1,3)
     .               - az(i,j,k)* hessian(2,3,3)

      ! l = 3, m = 3
        tensor(3,3) =  az(i,j,k)*(hessian(1,1,3)
     .                          + hessian(2,2,3))
     .               - ax(i,j,k)* hessian(3,1,3)
     .               - ay(i,j,k)* hessian(3,2,3)

c     End program

      end function nabla_v_src

      end function fnabla_v

c     fnabla_v_upwd
c     #############################################################
      function fnabla_v_upwd(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
     .                      ,hex,hey,hez) result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor nabla(vec v) at the following positions:
c       + hex,hey,hez = 0 => i,j,k
c       + hex=+-1 --> i+-1/2
c       + hey=+-1 --> j+-1/2
c       + hez=+-1 --> k+-1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,hex,hey,hez,nx,ny,nz,igx,igy,igz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
        real(8)    :: dhx,dhy,dhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: hessian(3,3,3)
        logical    :: cartsn

c     Begin program

        igrid = igx

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

        hessian = -gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

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

cc        if (sing_point) then
cc          vxip = ax(ip,j,k)+ax(i,j,k)
cc          vxim = 2.*ax(im,j,k)
cc          vyip = ay(ip,j,k)+ay(i,j,k)
cc          vyim = 2.*ay(im,j,k)
cc          vzip = az(ip,j,k)+az(i,j,k)
cc          vzim = 2.*az(im,j,k)
cc        else
          vxip = ax(ip,j,k)
          vxim = ax(im,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(im,j,k)
          vzip = az(ip,j,k)
          vzim = az(im,j,k)
cc        endif

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

c     btensor_x
c     #############################################################
      subroutine btensor_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer(4) :: ig,jg,kg,ip
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp,ptot,vis

c     Begin program

        ip = i+1
        if (flag == 0) ip = i

        jacp = gmetric%grid(igx)%jac(ip,j,k)
        jac0 = gmetric%grid(igx)%jac(i ,j,k)
        jac  = 0.5*(jacp+jac0)

cc        if (bcond(1) == SP .and. i == 0) jac = 1d-10
        if (isSP(i+1,j,k,igx,igy,igz)) jac = 1d-10

        t11 = 0d0

        if (solenoidal) then

          t12 = 0.5*( vec1(i ,j,k,1)*vec2(i ,j,k,2)/jac0
     .               +vec1(ip,j,k,1)*vec2(ip,j,k,2)/jacp)*jac
     .         -0.5*( vec1(i ,j,k,2)*vec2(i ,j,k,1)/jac0
     .               +vec1(ip,j,k,2)*vec2(ip,j,k,1)/jacp)*jac

cc          if (sing_point .and. flag == 1) then
cc            t13 = 0.5*( vec1(i ,j,k,1)*vec2(i ,j,k,3)/jac0**2
cc     .                 +vec1(ip,j,k,1)*vec2(ip,j,k,3)/jacp**2)*jac**2
cc     .           -0.5*( vec1(i ,j,k,3)*vec2(i ,j,k,1)/jac0**2
cc     .                 +vec1(ip,j,k,3)*vec2(ip,j,k,1)/jacp**2)*jac**2
cc          else
            t13 = 0.5*( vec1(i ,j,k,1)*vec2(i ,j,k,3)/jac0
     .                 +vec1(ip,j,k,1)*vec2(ip,j,k,3)/jacp)*jac
     .           -0.5*( vec1(i ,j,k,3)*vec2(i ,j,k,1)/jac0
     .                 +vec1(ip,j,k,3)*vec2(ip,j,k,1)/jacp)*jac
cc          endif
        else
          t12 = 0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,2)/jacp
     .               +vec1(i ,j,k,1)*vec2(ip,j,k,2)/jac0)*jac
     .         -0.5*( vec1(ip,j,k,2)*vec2(i ,j,k,1)/jac0
     .               +vec1(i ,j,k,2)*vec2(ip,j,k,1)/jacp)*jac

cc          if (sing_point .and. flag == 1) then
cc            t13 = 0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,3)/jac0/jacp
cc     .                 +vec1(i ,j,k,1)*vec2(ip,j,k,3)/jacp/jac0)*jac**2
cc     .           -0.5*( vec1(ip,j,k,3)*vec2(i ,j,k,1)/jac0/jacp
cc     .                 +vec1(i ,j,k,3)*vec2(ip,j,k,1)/jacp/jac0)*jac**2
cc          else
            t13 = 0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,3)/jac0
     .                 +vec1(i ,j,k,1)*vec2(ip,j,k,3)/jacp)*jac
     .           -0.5*( vec1(ip,j,k,3)*vec2(i ,j,k,1)/jac0
     .                 +vec1(i ,j,k,3)*vec2(ip,j,k,1)/jacp)*jac

cc          endif
        endif

        if (flag /= 0) then
          t11 = t11/jac
          if (.not.alt_eom) t12 = t12/jac
          t13 = t13/jac
        endif

c     End program

      end subroutine btensor_x

c     btensor_y
c     #############################################################
      subroutine btensor_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer(4) :: ig,jg,kg,jp
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp

c     Begin program

        jp = j+1
        if (flag == 0) jp = j

        jacp = gmetric%grid(igx)%jac(i,jp,k)
        jac0 = gmetric%grid(igx)%jac(i,j ,k)
        jac  = 0.5*(jacp+jac0)

        if (solenoidal) then
          t21 = 0.5*( vec1(i,j ,k,2)*vec2(i,j ,k,1)/jac0
     .               +vec1(i,jp,k,2)*vec2(i,jp,k,1)/jacp)*jac
     .         -0.5*( vec1(i,j ,k,1)*vec2(i,j ,k,2)/jac0
     .               +vec1(i,jp,k,1)*vec2(i,jp,k,2)/jacp)*jac

          t22 = 0d0

          t23 = 0.5*( vec1(i,j ,k,2)*vec2(i,j ,k,3)/jac0
     .               +vec1(i,jp,k,2)*vec2(i,jp,k,3)/jacp)*jac
     .         -0.5*( vec1(i,j ,k,3)*vec2(i,j ,k,2)/jac0
     .               +vec1(i,jp,k,3)*vec2(i,jp,k,2)/jacp)*jac
        else
          t21 = 0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,1)/jac0
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,1)/jacp)*jac
     .         -0.5*( vec1(i,jp,k,1)*vec2(i,j ,k,2)/jacp
     .               +vec1(i,j ,k,1)*vec2(i,jp,k,2)/jac0)*jac

          t22 = 0d0

          t23 = 0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,3)/jac0
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,3)/jacp)*jac
     .         -0.5*( vec1(i,jp,k,3)*vec2(i,j ,k,2)/jac0
     .               +vec1(i,j ,k,3)*vec2(i,jp,k,2)/jacp)*jac
        endif

        if (flag /= 0) then
          t21 = t21/jac
          if (.not.alt_eom) t22 = t22/jac
          t23 = t23/jac
        endif

c     End program

      end subroutine btensor_y

c     btensor_z
c     #############################################################
      subroutine btensor_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer(4) :: ig,jg,kg,kp
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp

c     Begin program

        kp = k+1
        if (flag == 0) kp = k

        jacp = gmetric%grid(igx)%jac(i,j,kp)
        jac0 = gmetric%grid(igx)%jac(i,j,k )
        jac  = 0.5*(jacp+jac0)

        if (solenoidal) then
          t31 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,k ,1)/jac0
     .               +vec1(i,j,kp,3)*vec2(i,j,kp,1)/jacp)*jac
     .         -0.5*( vec1(i,j,k ,1)*vec2(i,j,k ,3)/jac0
     .               +vec1(i,j,kp,1)*vec2(i,j,kp,3)/jacp)*jac

          t32 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,k ,2)/jac0
     .               +vec1(i,j,kp,3)*vec2(i,j,kp,2)/jacp)*jac
     .         -0.5*( vec1(i,j,k ,2)*vec2(i,j,k ,3)/jac0
     .               +vec1(i,j,kp,2)*vec2(i,j,kp,3)/jacp)*jac
        else
          t31 = 0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,1)/jac0
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,1)/jacp)*jac
     .         -0.5*( vec1(i,j,kp,1)*vec2(i,j,k ,3)/jac0
     .               +vec1(i,j,k ,1)*vec2(i,j,kp,3)/jacp)*jac

          t32 = 0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,2)/jac0
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,2)/jacp)*jac
     .         -0.5*( vec1(i,j,kp,2)*vec2(i,j,k ,3)/jac0
     .               +vec1(i,j,k ,2)*vec2(i,j,kp,3)/jacp)*jac
        endif

        t33 = 0d0

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alt_eom) t32 = t32/jac
          t33 = t33/jac
        endif

c     End program

      end subroutine btensor_z

c     lf_x
c     #############################################################
      subroutine lf_x(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .               ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for linearized Lorentz
c     force term in conservative form.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer(4) :: ip,igrid
        real(8)    :: jac,jac0,jacp,gsuper(3,3),scalar_prod
     .               ,acnv(3),acnvp(3),bcov(3),bcovp(3),bcnv(3),bcnvp(3)

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

cc        if (bcond(1) == SP .and. i == 0) jac = 1d-10
        if (isSP(i+1,j,k,igx,igy,igz)) jac = 1d-10

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igrid)%jac(ip,j,k)
          jac0 = gmetric%grid(igrid)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

cc        acnv(1) = 0.5*(vec1(ip,j,k,1)/jacp+vec1(i,j,k,1)/jac0)*jac
cc        acnv(2) = 0.5*(vec1(ip,j,k,2)     +vec1(i,j,k,2))
cc        acnv(3) = 0.5*(vec1(ip,j,k,3)/jacp+vec1(i,j,k,3)/jac0)*jac
cc
cc        bcnv(1) = 0.5*(vec2(ip,j,k,1)/jacp+vec2(i,j,k,1)/jac0)*jac
cc        bcnv(2) = 0.5*(vec2(ip,j,k,2)     +vec2(i,j,k,2))
cc        bcnv(3) = 0.5*(vec2(ip,j,k,3)/jacp+vec2(i,j,k,3)/jac0)*jac
cc
cc        if (flag /= 0) then
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=1)
cc        else
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=0)
cc        endif
cc
cc        scalar_prod = dot_product(acnv,bcov)
cc
cc        t11 =( acnv(1)*bcnv(1)
cc     .        +acnv(1)*bcnv(1)
cc     .        -gsuper(1,1)*scalar_prod )
cc
cc        t12 =( acnv(1)*bcnv(2)
cc     .        +acnv(2)*bcnv(1)
cc     .        -gsuper(1,2)*scalar_prod )
cc
cc        t13 =( acnv(1)*bcnv(3)
cc     .        +acnv(3)*bcnv(1)
cc     .        -gsuper(1,3)*scalar_prod )

        acnv(1) = vec1(i,j,k,1)/jac0*jac
        acnv(2) = vec1(i,j,k,2)
        acnv(3) = vec1(i,j,k,3)/jac0*jac

        acnvp(1) = vec1(ip,j,k,1)/jacp*jac
        acnvp(2) = vec1(ip,j,k,2)
        acnvp(3) = vec1(ip,j,k,3)/jacp*jac

cc        bcnv(1) = vec2(i,j,k,1)/jac0*jac
cc        bcnv(2) = vec2(i,j,k,2)
cc        bcnv(3) = vec2(i,j,k,3)/jac0*jac
cc
cc        bcnvp(1) = vec2(ip,j,k,1)/jacp*jac
cc        bcnvp(2) = vec2(ip,j,k,2)
cc        bcnvp(3) = vec2(ip,j,k,3)/jacp*jac

        bcnv(1) = vec2(i,j,k,1)
        bcnv(2) = vec2(i,j,k,2)
        bcnv(3) = vec2(i,j,k,3)

        bcnvp(1) = vec2(ip,j,k,1)
        bcnvp(2) = vec2(ip,j,k,2)
        bcnvp(3) = vec2(ip,j,k,3)

        call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,bcov(1),bcov(2),bcov(3)
     .                        ,bcnv(1),bcnv(2),bcnv(3)
     .                        ,.false.,half_elem=0)

        call transformFromCurvToCurv(ip,j,k,igx,igy,igz
     .                        ,bcovp(1),bcovp(2),bcovp(3)
     .                        ,bcnvp(1),bcnvp(2),bcnvp(3)
     .                        ,.false.,half_elem=0)

        bcnv(1) = bcnv(1)/jac0*jac
        bcnv(3) = bcnv(3)/jac0*jac

        bcnvp(1) = bcnvp(1)/jacp*jac
        bcnvp(3) = bcnvp(3)/jacp*jac

        bcov (2) = bcov (2)/jac0*jac
        bcovp(2) = bcovp(2)/jacp*jac

        scalar_prod = dot_product(acnvp,bcov)
     .               +dot_product(acnv,bcovp)

        t11 =0.5*( 2*acnvp(1)*bcnv(1)
     .            +2*acnv(1)*bcnvp(1)
     .            -gsuper(1,1)*scalar_prod )

        t12 =0.5*( acnvp(1)*bcnv(2) + acnv(1)*bcnvp(2)
     .            +acnvp(2)*bcnv(1) + acnv(2)*bcnvp(1)
     .            -gsuper(1,2)*scalar_prod )

        t13 =0.5*( acnvp(1)*bcnv(3) + acnv(1)*bcnvp(3)
     .            +acnvp(3)*bcnv(1) + acnv(3)*bcnvp(1)
     .            -gsuper(1,3)*scalar_prod )

        if (flag /= 0) then
          t11 = t11/jac
          if (.not.alt_eom) t12 = t12/jac
          t13 = t13/jac
        endif

c     End program

      end subroutine lf_x

c     lf_y
c     #############################################################
      subroutine lf_y(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .               ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for linearized Lorentz
c     force term in conservative form.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer(4) :: jp,igrid
        real(8)    :: jac,gsuper(3,3),scalar_prod
     .               ,acnv(3),acnvp(3),bcov(3),bcovp(3),bcnv(3),bcnvp(3)

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

cc        acnv(1) = 0.5*(vec1(i,jp,k,1)+vec1(i,j,k,1))
cc        acnv(2) = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cc        acnv(3) = 0.5*(vec1(i,jp,k,3)+vec1(i,j,k,3))
cc                                
cc        bcnv(1) = 0.5*(vec2(i,jp,k,1)+vec2(i,j,k,1))
cc        bcnv(2) = 0.5*(vec2(i,jp,k,2)+vec2(i,j,k,2))
cc        bcnv(3) = 0.5*(vec2(i,jp,k,3)+vec2(i,j,k,3))
cc
cc        if (flag /= 0) then
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=2)
cc        else
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=0)
cc        endif
cc
cc
cc        scalar_prod = dot_product(acnv,bcov)
cc
cc        t21 =( acnv(2)*bcnv(1)
cc     .        +acnv(1)*bcnv(2)
cc     .        -gsuper(2,1)*scalar_prod )
cc
cc        t22 =( acnv(2)*bcnv(2)
cc     .        +acnv(2)*bcnv(2)
cc     .        -gsuper(2,2)*scalar_prod )
cc
cc        t23 =( acnv(2)*bcnv(3)
cc     .        +acnv(3)*bcnv(2)
cc     .        -gsuper(2,3)*scalar_prod )

        acnv(1) = vec1(i,j,k,1)
        acnv(2) = vec1(i,j,k,2)
        acnv(3) = vec1(i,j,k,3)

        acnvp(1) = vec1(i,jp,k,1)
        acnvp(2) = vec1(i,jp,k,2)
        acnvp(3) = vec1(i,jp,k,3)

        bcnv(1) = vec2(i,j,k,1)
        bcnv(2) = vec2(i,j,k,2)
        bcnv(3) = vec2(i,j,k,3)

        bcnvp(1) = vec2(i,jp,k,1)
        bcnvp(2) = vec2(i,jp,k,2)
        bcnvp(3) = vec2(i,jp,k,3)

        call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,bcov(1),bcov(2),bcov(3)
     .                        ,bcnv(1),bcnv(2),bcnv(3)
     .                        ,.false.,half_elem=0)

        call transformFromCurvToCurv(i,jp,k,igx,igy,igz
     .                        ,bcovp(1),bcovp(2),bcovp(3)
     .                        ,bcnvp(1),bcnvp(2),bcnvp(3)
     .                        ,.false.,half_elem=0)

        scalar_prod = dot_product(acnvp,bcov)
     .               +dot_product(acnv,bcovp)

        t21 =0.5*( acnvp(2)*bcnv(1) + acnv(2)*bcnvp(1)
     .            +acnvp(1)*bcnv(2) + acnv(1)*bcnvp(2)
     .            -gsuper(2,1)*scalar_prod )

        t22 =0.5*( 2*acnvp(2)*bcnv(2) + 2*acnv(2)*bcnvp(2)
     .            -gsuper(2,2)*scalar_prod )

        t23 =0.5*( acnvp(2)*bcnv(3) + acnv(2)*bcnvp(3)
     .            +acnvp(3)*bcnv(2) + acnv(3)*bcnvp(2)
     .            -gsuper(2,3)*scalar_prod )

        if (flag /= 0) then
          t21 = t21/jac
          if (.not.alt_eom) t22 = t22/jac
          t23 = t23/jac
        endif


c     End program

      end subroutine lf_y

c     lf_z
c     #############################################################
      subroutine lf_z(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for linearized Lorentz
c     force term in conservative form.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer(4) :: kp,igrid
        real(8)    :: jac,gsuper(3,3),scalar_prod
     .               ,acnv(3),acnvp(3),bcov(3),bcovp(3),bcnv(3),bcnvp(3)

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

cc        acnv(1) = 0.5*(vec1(i,j,kp,1)+vec1(i,j,k,1))
cc        acnv(2) = 0.5*(vec1(i,j,kp,2)+vec1(i,j,k,2))
cc        acnv(3) = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cc                                  
cc        bcnv(1) = 0.5*(vec2(i,j,kp,1)+vec2(i,j,k,1))
cc        bcnv(2) = 0.5*(vec2(i,j,kp,2)+vec2(i,j,k,2))
cc        bcnv(3) = 0.5*(vec2(i,j,kp,3)+vec2(i,j,k,3))
cc
cc        if (flag /= 0) then
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=3)
cc        else
cc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=0)
cc        endif
cc
cc        scalar_prod = dot_product(acnv,bcov)
cc
cc        t31 =( acnv(3)*bcnv(1)
cc     .        +acnv(1)*bcnv(3)
cc     .        -gsuper(3,1)*scalar_prod )
cc
cc        t32 =( acnv(3)*bcnv(2)
cc     .        +acnv(2)*bcnv(3)
cc     .        -gsuper(3,2)*scalar_prod )
cc
cc        t33 =( acnv(3)*bcnv(3)
cc     .        +acnv(3)*bcnv(3)
cc     .        -gsuper(3,3)*scalar_prod )

        acnv(1) = vec1(i,j,k,1)
        acnv(2) = vec1(i,j,k,2)
        acnv(3) = vec1(i,j,k,3)

        acnvp(1) = vec1(i,j,kp,1)
        acnvp(2) = vec1(i,j,kp,2)
        acnvp(3) = vec1(i,j,kp,3)

        bcnv(1) = vec2(i,j,k,1)
        bcnv(2) = vec2(i,j,k,2)
        bcnv(3) = vec2(i,j,k,3)

        bcnvp(1) = vec2(i,j,kp,1)
        bcnvp(2) = vec2(i,j,kp,2)
        bcnvp(3) = vec2(i,j,kp,3)

        call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,bcov(1),bcov(2),bcov(3)
     .                        ,bcnv(1),bcnv(2),bcnv(3)
     .                        ,.false.,half_elem=0)

        call transformFromCurvToCurv(i,j,kp,igx,igy,igz
     .                        ,bcovp(1),bcovp(2),bcovp(3)
     .                        ,bcnvp(1),bcnvp(2),bcnvp(3)
     .                        ,.false.,half_elem=0)

        scalar_prod = dot_product(acnvp,bcov)
     .               +dot_product(acnv,bcovp)

        t31 =0.5*( acnvp(3)*bcnv(1) + acnv(3)*bcnvp(1)
     .            +acnvp(1)*bcnv(3) + acnv(1)*bcnvp(3)
     .            -gsuper(3,1)*scalar_prod )

        t32 =0.5*( acnvp(3)*bcnv(2) + acnv(3)*bcnvp(2)
     .            +acnvp(2)*bcnv(3) + acnv(2)*bcnvp(3)
     .            -gsuper(3,2)*scalar_prod )

        t33 =0.5*( 2*acnvp(3)*bcnv(3)
     .           + 2*acnv(3)*bcnvp(3)
     .            -gsuper(3,3)*scalar_prod )

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alt_eom) t32 = t32/jac
          t33 = t33/jac
        endif

c     End program

      end subroutine lf_z

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

        real(8),pointer,dimension(:,:,:):: rho,rvx,rvy,rvz,bx,by,bz,tmp

        logical :: nc_eom=.false.

      contains

c     c_advec
c     ###############################################################
      function c_advec(i,j,k,nx,ny,nz,igx,igy,igz,phi,sp_upwind,vol_wgt)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.phi) at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1),c_advec

      integer(4),optional :: sp_upwind
      logical,optional :: vol_wgt

c     Local variables

      integer(4) :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km,su
      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacip,jacim,jac,jach
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical    :: vw

c     Begin program

      if (PRESENT(vol_wgt)) then
        vw = vol_wgt
      else
        vw = .true.
      endif

      if (PRESENT(sp_upwind)) then
        su = sp_upwind
      else
        su = 0
      endif

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      jac   = gmetric%grid(igrid)%jac(i ,j,k)
      jacip = gmetric%grid(igrid)%jac(ip,j,k)
      jacim = gmetric%grid(igrid)%jac(im,j,k)

      if (bcSP()) then
        if (i+grid_params%ilo(igx)-1 < su) then  !Upwind around singular point
          jach = 0.5*(jac+jacip)
          flxip = 0.25*jach
     .        *( (    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)
     .            +abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *phi(i ,j,k)
     .          +(    (vx(i,j,k)/jac+vx(ip,j,k)/jacip)          
     .            -abs(vx(i,j,k)/jac+vx(ip,j,k)/jacip) ) *phi(ip,j,k) )

          jach = 0.5*(jac+jacim)
          flxim = 0.25*jach
     .         *( (    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
     .             +abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *phi(im,j,k)
     .           +(    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
     .             -abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *phi(i ,j,k) )
        elseif (i+grid_params%ilo(igx)-1 == su) then  !Upwind around singular point
          jach = 0.5*(jac+jacip)
          flxip = 0.5*(vx(ip,j,k)*phi(i ,j,k)/jacip
     .               + vx(i ,j,k)*phi(ip,j,k)/jac  )*jach

          jach = 0.5*(jac+jacim)
          flxim = 0.25*jach
     .         *( (    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
     .             +abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *phi(im,j,k)
     .           +(    (vx(i,j,k)/jac+vx(im,j,k)/jacim)          
     .             -abs(vx(i,j,k)/jac+vx(im,j,k)/jacim) ) *phi(i ,j,k) )
        elseif (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
cc        if (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
          jach = 0.5*(jac+jacip)
          flxip = 0.5*(vx(ip,j,k)*phi(i ,j,k)/jacip
     .               + vx(i ,j,k)*phi(ip,j,k)/jac  )*jach
          jach = 0.5*(jac+jacim)
          flxim = 0.5*(vx(im,j,k)*phi(i ,j,k)/jacim
     .               + vx(i ,j,k)*phi(im,j,k)/jac  )*jach
        elseif (i+grid_params%ilo(igx)-1 == grid_params%nxgl(igx)) then
          flxip = 0.5*(vx(ip,j,k)*phi(i,j,k) + vx(i,j,k)*phi(ip,j,k))

          jach = 0.5*(jac+jacim)
          flxim = 0.5*(vx(im,j,k)*phi(i ,j,k)/jacim
     .               + vx(i ,j,k)*phi(im,j,k)/jac  )*jach
        else
          flxip = 0.5*(vx(ip,j,k)*phi(i,j,k) + vx(i,j,k)*phi(ip,j,k))
          flxim = 0.5*(vx(im,j,k)*phi(i,j,k) + vx(i,j,k)*phi(im,j,k))
        endif
      else
        flxip = 0.5*(vx(ip,j,k)*phi(i,j,k) + vx(i,j,k)*phi(ip,j,k))
        flxim = 0.5*(vx(im,j,k)*phi(i,j,k) + vx(i,j,k)*phi(im,j,k))
      endif

      !Y flux
      if (bcSP() .and. (i+grid_params%ilo(igx)-1 < su) ) then !Upwind around singular point
        flxjp = 0.25*( (    (vy(i,j,k)+vy(i,jp,k))
     .                  +abs(vy(i,j,k)+vy(i,jp,k)) ) *phi(i,j ,k)
     .                +(    (vy(i,j,k)+vy(i,jp,k))
     .                  -abs(vy(i,j,k)+vy(i,jp,k)) ) *phi(i,jp,k) )
        flxjm = 0.25*( (    (vy(i,j,k)+vy(i,jm,k))
     .                  +abs(vy(i,j,k)+vy(i,jm,k)) ) *phi(i,jm,k)
     .                +(    (vy(i,j,k)+vy(i,jm,k))
     .                  -abs(vy(i,j,k)+vy(i,jm,k)) ) *phi(i,j ,k) )
      else
        flxjp = 0.5*(vy(i,jp,k)*phi(i,j,k) + vy(i,j,k)*phi(i,jp,k))
        flxjm = 0.5*(vy(i,jm,k)*phi(i,j,k) + vy(i,j,k)*phi(i,jm,k))
      endif

      !Z flux
      flxkp = 0.5*(vz(i,j,kp)*phi(i,j,k) + vz(i,j,k)*phi(i,j,kp))
      flxkm = 0.5*(vz(i,j,km)*phi(i,j,k) + vz(i,j,k)*phi(i,j,km))

      c_advec =( (flxip - flxim)/dxx
     .         + (flxjp - flxjm)/dyy
     .         + (flxkp - flxkm)/dzz )/jac
      
      if (vw) c_advec = c_advec*volume(i,j,k,igx,igy,igz)

c     End 

      end function c_advec

c     vtensor_x
c     #############################################################
      subroutine vtensor_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t11,t12,t13: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer(4) :: ig,jg,kg,ip
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp,ptot,vis

c     Begin program

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
     .               +gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))

cc        if (bcond(1) == SP .and. i == 0) jac = 1d-10
        if (isSP(i+1,j,k,igx,igy,igz)) jac = 1d-10

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac0 = gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,1)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)
        endif

        !Harmonic average for calculation of viscosity coeff. at faces
        vis = 2./(1./rho(ip,j,k)/nuu(ip,j,k)
     .          + 1./rho(i ,j,k)/nuu(i ,j,k))

        if (nc_eom) then
         !Recall p=2nT
          ptot = jac*(rho(ip,j,k)*tmp(i ,j,k)
     .               +rho(i ,j,k)*tmp(ip,j,k))
          t11 =
     .      0.5*( rvx(ip,j,k)*vx (i ,j,k)/jacp/jac0
     .           +rvx(i ,j,k)*vx (ip,j,k)/jacp/jac0)*jac**2
     .     +gsuper(1,1)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,1)
     .           +gsuper(1,2)*nabla_v(2,1)
     .           +gsuper(1,3)*nabla_v(3,1) )

          t12 =
     .     0.25*( rvx(ip,j,k)*vy (i ,j,k)/jacp
     .           +rvx(i ,j,k)*vy (ip,j,k)/jac0
     .           + vx(i ,j,k)*rvy(ip,j,k)/jac0
     .           + vx(ip,j,k)*rvy(i ,j,k)/jacp)*jac
     .     +gsuper(1,2)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,2)
     .           +gsuper(1,2)*nabla_v(2,2)
     .           +gsuper(1,3)*nabla_v(3,2) )

          t13 =
     .      0.25*( rvx(ip,j,k)*vz (i ,j,k)/jacp/jac0
     .            +rvx(i ,j,k)*vz (ip,j,k)/jacp/jac0
     .            + vx(i ,j,k)*rvz(ip,j,k)/jacp/jac0
     .            + vx(ip,j,k)*rvz(i ,j,k)/jacp/jac0)*jac**2
     .     +gsuper(1,3)*ptot
     .     -vis*( gsuper(1,1)*nabla_v(1,3)
     .           +gsuper(1,2)*nabla_v(2,3)
     .           +gsuper(1,3)*nabla_v(3,3) )
        else
         !Recall p=2nT
          ptot = jac*(rho(ip,j,k)*tmp(i ,j,k)
     .               +rho(i ,j,k)*tmp(ip,j,k))
cc          ptot = jac*(rho(ip,j,k)*tmp(ip,j,k)
cc     .               +rho(i ,j,k)*tmp(i ,j,k))
     .       +jac*(bx(ip,j,k)*bx_cov(i ,j,k)/jacp
     .            +bx(i ,j,k)*bx_cov(ip,j,k)/jac0
     .            +by(ip,j,k)*by_cov(i ,j,k)/jac0
     .            +by(i ,j,k)*by_cov(ip,j,k)/jacp
     .            +bz(ip,j,k)*bz_cov(i ,j,k)/jacp
     .            +bz(i ,j,k)*bz_cov(ip,j,k)/jac0)/4.

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
        endif

        if (flag /= 0) then
          t11 = t11/jac
          if (.not.alt_eom) t12 = t12/jac
          t13 = t13/jac
        endif

c     End program

      end subroutine vtensor_x

c     vtensor_y
c     #############################################################
      subroutine vtensor_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EOM. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t21,t22,t23: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

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
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,2)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)
        endif

        vis = 2./(1./rho(i,jp,k)/nuu(i,jp,k)
     .          + 1./rho(i,j ,k)/nuu(i,j ,k))

        if (nc_eom) then
          ptot = jac*(rho(i,jp,k)*tmp(i,j ,k)
     .               +rho(i,j ,k)*tmp(i,jp,k))

          t21 =
     .       0.25*( rvy(i,jp,k)*vx(i,j,k) + rvy(i,j,k)*vx(i,jp,k)
     .             +rvx(i,jp,k)*vy(i,j,k) + rvx(i,j,k)*vy(i,jp,k))
     .       +gsuper(2,1)*ptot
     .       -vis*( gsuper(2,1)*nabla_v(1,1)
     .             +gsuper(2,2)*nabla_v(2,1)
     .             +gsuper(2,3)*nabla_v(3,1) )

          t22 =
     .       0.25*( rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k)
     .             +rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k))
     .       +gsuper(2,2)*ptot
     .       -vis*( gsuper(2,1)*nabla_v(1,2)
     .             +gsuper(2,2)*nabla_v(2,2)
     .             +gsuper(2,3)*nabla_v(3,2) )

          t23 =
     .       0.25*( rvy(i,jp,k)*vz(i,j,k) + rvy(i,j,k)*vz(i,jp,k)
     .             +rvz(i,jp,k)*vy(i,j,k) + rvz(i,j,k)*vy(i,jp,k))
     .       +gsuper(2,3)*ptot
     .       -vis*( gsuper(2,1)*nabla_v(1,3)
     .             +gsuper(2,2)*nabla_v(2,3)
     .             +gsuper(2,3)*nabla_v(3,3) )
        else
          !Recall p=2nT
cc          ptot = jac*(rho(i,jp,k)*tmp(i,jp,k)
cc     .               +rho(i,j ,k)*tmp(i,j ,k))
          ptot = jac*(rho(i,jp,k)*tmp(i,j ,k)
     .               +rho(i,j ,k)*tmp(i,jp,k))
     .         +(bx(i,jp,k)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,jp,k)
     .          +by(i,jp,k)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,jp,k)
     .          +bz(i,jp,k)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,jp,k))*0.25

          t21 =
     .       0.25*( rvy(i,jp,k)*vx(i,j,k) + rvy(i,j,k)*vx(i,jp,k)
     .             +rvx(i,jp,k)*vy(i,j,k) + rvx(i,j,k)*vy(i,jp,k))
     .       -0.5*( by(i,jp,k)*bx(i,j ,k)
     .             +by(i,j ,k)*bx(i,jp,k) )
     .       +gsuper(2,1)*ptot
     .       -vis*( gsuper(2,1)*nabla_v(1,1)
     .             +gsuper(2,2)*nabla_v(2,1)
     .             +gsuper(2,3)*nabla_v(3,1) )

          t22 =
     .       0.25*( rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k)
     .             +rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k))
     .       -0.5*( by(i,jp,k)*by(i,j ,k)
     .             +by(i,j ,k)*by(i,jp,k) )
     .       +gsuper(2,2)*ptot
     .       -vis*( gsuper(2,1)*nabla_v(1,2)
     .             +gsuper(2,2)*nabla_v(2,2)
     .             +gsuper(2,3)*nabla_v(3,2) )

          t23 =
     .       0.25*( rvy(i,jp,k)*vz(i,j,k) + rvy(i,j,k)*vz(i,jp,k)
     .             +rvz(i,jp,k)*vy(i,j,k) + rvz(i,j,k)*vy(i,jp,k))
     .       -0.5*( by(i,jp,k)*bz(i,j ,k)
     .             +by(i,j ,k)*bz(i,jp,k) )
     .       +gsuper(2,3)*ptot
     .       -vis*( gsuper(2,1)*nabla_v(1,3)
     .             +gsuper(2,2)*nabla_v(2,3)
     .             +gsuper(2,3)*nabla_v(3,3) )
        endif

        if (flag /= 0) then
          t21 = t21/jac
          if (.not.alt_eom) t22 = t22/jac
          t23 = t23/jac
        endif

c     End program

      end subroutine vtensor_y

c     vtensor_z
c     #############################################################
      subroutine vtensor_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EOM. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t31,t32,t33: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer(4) :: ig,jg,kg,kp
        real(8)    :: x,y,z
        real(8)    :: jac,ptot,vis

c     Begin program

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,3)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)
        endif

        vis = 2./(1./rho(i,j,kp)/nuu(i,j,kp)
     .          + 1./rho(i,j,k )/nuu(i,j,k ))

        if (nc_eom) then
          !Recall p=2nT
cc          ptot = jac*(rho(i,j,kp)*tmp(i,j,kp)
cc     .               +rho(i,j,k )*tmp(i,j,k ))
          ptot = jac*(rho(i,j,kp)*tmp(i,j,k )
     .               +rho(i,j,k )*tmp(i,j,kp))

          t31 =
     .       0.25*( rvz(i,j,kp)*vx(i,j,k) + rvz(i,j,k)*vx(i,j,kp)
     .             +rvx(i,j,kp)*vz(i,j,k) + rvx(i,j,k)*vz(i,j,kp) )
     .       +gsuper(3,1)*ptot
     .       -vis*( gsuper(3,1)*nabla_v(1,1)
     .             +gsuper(3,2)*nabla_v(2,1)
     .             +gsuper(3,3)*nabla_v(3,1) )

          t32 =
     .       0.25*( rvz(i,j,kp)*vy(i,j,k) + rvz(i,j,k)*vy(i,j,kp)
     .             +rvy(i,j,kp)*vz(i,j,k) + rvy(i,j,k)*vz(i,j,kp) )
     .       +gsuper(3,2)*ptot
     .       -vis*( gsuper(3,1)*nabla_v(1,2)
     .             +gsuper(3,2)*nabla_v(2,2)
     .             +gsuper(3,3)*nabla_v(3,2) )

          t33 =
     .       0.25*( rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp)
     .             +rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp) )
     .       +gsuper(3,3)*ptot
     .       -vis*( gsuper(3,1)*nabla_v(1,3)
     .             +gsuper(3,2)*nabla_v(2,3)
     .             +gsuper(3,3)*nabla_v(3,3) )
        else
          !Recall p=2nT
cc          ptot = jac*(rho(i,j,kp)*tmp(i,j,kp)
cc     .               +rho(i,j,k )*tmp(i,j,k ))
          ptot = jac*(rho(i,j,kp)*tmp(i,j,k )
     .               +rho(i,j,k )*tmp(i,j,kp))
     .         +(bx(i,j,kp)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,j,kp)
     .          +by(i,j,kp)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,j,kp)
     .          +bz(i,j,kp)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,j,kp))*0.25

          t31 =
     .       0.25*( rvz(i,j,kp)*vx(i,j,k) + rvz(i,j,k)*vx(i,j,kp)
     .             +rvx(i,j,kp)*vz(i,j,k) + rvx(i,j,k)*vz(i,j,kp) )
     .       -0.5*( bz(i,j,kp)*bx(i,j,k )
     .             +bz(i,j,k )*bx(i,j,kp) )
     .       +gsuper(3,1)*ptot
     .       -vis*( gsuper(3,1)*nabla_v(1,1)
     .             +gsuper(3,2)*nabla_v(2,1)
     .             +gsuper(3,3)*nabla_v(3,1) )

          t32 =
     .       0.25*( rvz(i,j,kp)*vy(i,j,k) + rvz(i,j,k)*vy(i,j,kp)
     .             +rvy(i,j,kp)*vz(i,j,k) + rvy(i,j,k)*vz(i,j,kp) )
     .       -0.5*( bz(i,j,kp)*by(i,j,k )
     .             +bz(i,j,k )*by(i,j,kp) )
     .       +gsuper(3,2)*ptot
     .       -vis*( gsuper(3,1)*nabla_v(1,2)
     .             +gsuper(3,2)*nabla_v(2,2)
     .             +gsuper(3,3)*nabla_v(3,2) )

          t33 =
     .       0.25*( rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp)
     .             +rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp) )
     .       -0.5*( bz(i,j,kp)*bz(i,j,k )
     .             +bz(i,j,k )*bz(i,j,kp) )
     .       +gsuper(3,3)*ptot
     .       -vis*( gsuper(3,1)*nabla_v(1,3)
     .             +gsuper(3,2)*nabla_v(2,3)
     .             +gsuper(3,3)*nabla_v(3,3) )
        endif

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alt_eom) t32 = t32/jac
          t33 = t33/jac
        endif

c     End program

      end subroutine vtensor_z

      end module nlfunction_setup
