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

        real(8) :: nu,eta,dd,chi,gamma,E0(3)=0d0,di,a_p

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
     .         ,vex,vey,vez
     .         ,eeta,nuu,ejx,ejy,ejz,pe,gpex,gpey,gpez

        real(8),target,allocatable,dimension(:,:,:,:) :: bcnv,vcnv,vecnv

c SI operator
        real(8),target,allocatable,dimension(:,:,:,:) :: b_n,p_n,dv_dt
c SI operator

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_car,by_car,bz_car
     .         ,jx_car,jy_car,jz_car
     .         ,vx_car,vy_car,vz_car
     .         ,divrgB,divrgV,divrgJ,Pflux
     .         ,qfactor,lambda,qfactr2,lambda2,p_tot

        logical :: alt_eom

        logical :: nc_eom_f=.false.,nc_eom_v=.false.

        real(8) :: k_si=0d0

      end module auxiliaryVariables

c module operators
c ######################################################################
      module operators

        use grid

        use grid_aliases

        use error

        real(8),pointer,dimension(:,:,:,:) :: vec,vec1,vec2
        real(8),pointer,dimension(:,:,:)   :: coef

        logical    :: solenoidal=.true.

        INTERFACE veclaplacian
          module procedure veclap_diff,veclap_ndiff
        end INTERFACE

        INTERFACE laplacian
          module procedure lap_diff,lap_ndiff
        end INTERFACE

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

c     advec
c     ###############################################################
      function advec(i,j,k,nx,ny,nz,igx,igy,igz,v1,v2,v3,phi,vol
     .              ,upwind)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.phi) at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
     .             ,v1 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v2 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v3 (0:nx+1,0:ny+1,0:nz+1),advec

      logical,optional :: vol,upwind

c     Local variables

      integer(4) :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km,su
      real(8)    :: dxx,dyy,dzz,x0,y0,z0
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical    :: vw,upw

c     Begin program

      if (PRESENT(vol)) then
        vw = vol
      else
        vw = .true.
      endif

      if (PRESENT(upwind)) then
        upw = upwind
      else
        upw = .false.
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
cc      jacip = gmetric%grid(igrid)%jac(ip,j,k)
cc      jacim = gmetric%grid(igrid)%jac(im,j,k)

      if (upw) then
        flxip = 0.25*(
     .         (    (v1(i,j,k)+v1(ip,j,k))
     .          +abs(v1(i,j,k)+v1(ip,j,k)) )*phi(i ,j,k)
     .        +(    (v1(i,j,k)+v1(ip,j,k))          
     .          -abs(v1(i,j,k)+v1(ip,j,k)) )*phi(ip,j,k)) 

        if (isSP(i,j,k,igx,igy,igz)) then
          flxim = 0d0
        else
          flxim = 0.25*(
     .         (    (v1(i,j,k)+v1(im,j,k))
     .          +abs(v1(i,j,k)+v1(im,j,k)) )*phi(im,j,k)
     .        +(    (v1(i,j,k)+v1(im,j,k))          
     .          -abs(v1(i,j,k)+v1(im,j,k)) )*phi(i ,j,k))
        endif

        flxjp = 0.25*(
     .         (    (v2(i,j,k)+v2(i,jp,k))
     .          +abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,j ,k)
     .        +(    (v2(i,j,k)+v2(i,jp,k))          
     .          -abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,jp,k))
        flxjm = 0.25*(
     .         (    (v2(i,j,k)+v2(i,jm,k))
     .          +abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,jm,k)
     .        +(    (v2(i,j,k)+v2(i,jm,k))          
     .          -abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,j ,k))

        flxkp = 0.25*(
     .         (    (v3(i,j,k)+v3(i,j,kp))
     .          +abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,k )
     .        +(    (v3(i,j,k)+v3(i,j,kp))             
     .          -abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,kp))
        flxkm = 0.25*(
     .         (    (v3(i,j,k)+v3(i,j,km))
     .          +abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,km)
     .        +(    (v3(i,j,k)+v3(i,j,km))             
     .          -abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,k ))

      else
cc        flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc        if (isSP(i,j,k,igx,igy,igz)) then
cc          flxim = 0d0
cc        else
cc          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cc        endif

cc        if (isSP(i,j,k,igx,igy,igz)) then
cc          flxip = 0.5*(v1(ip,j,k)*phi(ip,j,k) + v1(i,j,k)*phi(i,j,k))
cc          flxim = 0d0
cc        elseif (isSP(i-1,j,k,igx,igy,igz)) then
cc          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k) + v1(i,j,k)*phi(ip,j,k))
cc          flxim = 0.5*(v1(im,j,k)*phi(im,j,k) + v1(i,j,k)*phi(i ,j,k))
cc        else
          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cc        endif

        flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
        flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))

        flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
        flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))
      endif

      advec =( (flxip - flxim)/dxx
     .       + (flxjp - flxjm)/dyy
     .       + (flxkp - flxkm)/dzz )/jac
      
      if (vw) advec = advec*gmetric%grid(igrid)%dvol(i,j,k)

c     End 

      end function advec

c     c_advec
c     ###############################################################
      function c_advec(i,j,k,nx,ny,nz,igx,igy,igz,v1,v2,v3,phi
     .                ,sp,upwind,vol)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.phi) at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
     .             ,v1 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v2 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v3 (0:nx+1,0:ny+1,0:nz+1),c_advec

      logical,optional :: sp,vol,upwind

c     Local variables

      integer(4) :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km,su
      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacip,jacim,jac,jach
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical    :: vw,spoint,upw

c     Begin program

      if (PRESENT(vol)) then
        vw = vol
      else
        vw = .true.
      endif

      if (PRESENT(upwind)) then
        upw = upwind
      else
        upw = .false.
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

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      jac   = gmetric%grid(igrid)%jac(i,j,k)

      if (spoint) then
        jacip = gmetric%grid(igrid)%jac(ip,j,k)
        jacim = gmetric%grid(igrid)%jac(im,j,k)
      else
        jacip = jac
        jacim = jac
      endif

      if (upw) then
        !X flux
        if (     i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .or. bcond(2) == PER) then
          jach = 0.5*(jac+jacip)
          flxip = 0.25*jach
     .      *( (    (v1(i,j,k)/jac+v1(ip,j,k)/jacip)
     .          +abs(v1(i,j,k)/jac+v1(ip,j,k)/jacip) ) *phi(i ,j,k)
     .        +(    (v1(i,j,k)/jac+v1(ip,j,k)/jacip)          
     .          -abs(v1(i,j,k)/jac+v1(ip,j,k)/jacip) ) *phi(ip,j,k) )
        else
          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
        endif

        if (     i+grid_params%ilo(igx)-1 > 1
     .      .or. spoint
     .      .or. bcond(1) == PER) then
          jach = 0.5*(jac+jacim)
          flxim = 0.25*jach
     .       *( (    (v1(i,j,k)/jac+v1(im,j,k)/jacim)          
     .           +abs(v1(i,j,k)/jac+v1(im,j,k)/jacim) ) *phi(im,j,k)
     .         +(    (v1(i,j,k)/jac+v1(im,j,k)/jacim)          
     .           -abs(v1(i,j,k)/jac+v1(im,j,k)/jacim) ) *phi(i ,j,k) )
        else
          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
        endif

        if (     j+grid_params%jlo(igy)-1 < grid_params%nygl(igy)
     .      .or. bcond(4) == PER) then
          flxjp = 0.25*(
     .         (    (v2(i,j,k)+v2(i,jp,k))
     .          +abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,j ,k)
     .        +(    (v2(i,j,k)+v2(i,jp,k))          
     .          -abs(v2(i,j,k)+v2(i,jp,k)) )*phi(i,jp,k))
        else
          flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
        endif

        if (     j+grid_params%jlo(igy)-1 > 1
     .      .or. bcond(3) == PER) then
          flxjm = 0.25*(
     .         (    (v2(i,j,k)+v2(i,jm,k))
     .          +abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,jm,k)
     .        +(    (v2(i,j,k)+v2(i,jm,k))          
     .          -abs(v2(i,j,k)+v2(i,jm,k)) )*phi(i,j ,k))
        else
          flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))
        endif

        if (     k+grid_params%klo(igz)-1 < grid_params%nzgl(igz)
     .      .or. bcond(6) == PER) then
          flxkp = 0.25*(
     .         (    (v3(i,j,k)+v3(i,j,kp))
     .          +abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,k )
     .        +(    (v3(i,j,k)+v3(i,j,kp))             
     .          -abs(v3(i,j,k)+v3(i,j,kp)) )*phi(i,j,kp))
        else
          flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
        endif

        if (     k+grid_params%klo(igz)-1 > 1
     .      .or. bcond(5) == PER) then
          flxkm = 0.25*(
     .         (    (v3(i,j,k)+v3(i,j,km))
     .          +abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,km)
     .        +(    (v3(i,j,k)+v3(i,j,km))             
     .          -abs(v3(i,j,k)+v3(i,j,km)) )*phi(i,j,k ))
        else
          flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))
        endif

      else
        !X flux
        if (i+grid_params%ilo(igx)-1 < grid_params%nxgl(igx)) then
          jach = 0.5*(jac+jacip)
          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k)/jacip
     .               + v1(i ,j,k)*phi(ip,j,k)/jac  )*jach
          jach = 0.5*(jac+jacim)
          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
        elseif (i+grid_params%ilo(igx)-1 == grid_params%nxgl(igx)) then
          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))

          jach = 0.5*(jac+jacim)
          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
        endif
  
        !Y flux
        flxjp = 0.5*(v2(i,jp,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jp,k))
        flxjm = 0.5*(v2(i,jm,k)*phi(i,j,k) + v2(i,j,k)*phi(i,jm,k))

        !Z flux
        flxkp = 0.5*(v3(i,j,kp)*phi(i,j,k) + v3(i,j,k)*phi(i,j,kp))
        flxkm = 0.5*(v3(i,j,km)*phi(i,j,k) + v3(i,j,k)*phi(i,j,km))

      endif

      c_advec =( (flxip - flxim)/dxx
     .         + (flxjp - flxjm)/dyy
     .         + (flxkp - flxkm)/dzz )/jac
      
      if (vw) c_advec = c_advec*volume(i,j,k,igx,igy,igz)

c     End 

      end function c_advec

c     lap_diff
c     ###############################################################
      function lap_diff(i,j,k,nx,ny,nz,igx,igy,igz,arr,dff,vol)
     .         result (laplacian)
c     ---------------------------------------------------------------
c     Calculates lap(arr) at cell centers in general non-orthog.
c     coordinates, preserving the SPD property.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian

      real(8)    :: dff(0:nx+1,0:ny+1,0:nz+1)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup

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
      endif

      gsup => gmetric%grid(igrid)%gsup

      d_xx_ip = 0.5*(dff(i ,j,k)*gsup(i ,j,k,1,1)
     .              +dff(ip,j,k)*gsup(ip,j,k,1,1))
      d_xx_im = 0.5*(dff(i ,j,k)*gsup(i ,j,k,1,1)
     .              +dff(im,j,k)*gsup(im,j,k,1,1))
      d_yy_jp = 0.5*(dff(i,j ,k)*gsup(i,j ,k,2,2)
     .              +dff(i,jp,k)*gsup(i,jp,k,2,2))
      d_yy_jm = 0.5*(dff(i,j ,k)*gsup(i,j ,k,2,2)
     .              +dff(i,jm,k)*gsup(i,jm,k,2,2))
      d_zz_kp = 0.5*(dff(i,j,k )*gsup(i,j,k ,3,3)
     .              +dff(i,j,kp)*gsup(i,j,kp,3,3))
      d_zz_km = 0.5*(dff(i,j,k )*gsup(i,j,k ,3,3)
     .              +dff(i,j,km)*gsup(i,j,km,3,3))

      d_xy_ipjp = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,2)
     .                 +dff(ip,j,k)*gsup(ip,j,k,1,2)
     .                 +dff(i,j ,k)*gsup(i,j ,k,1,2)
     .                 +dff(i,jp,k)*gsup(i,jp,k,1,2))
      d_xy_ipjm = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,2)
     .                 +dff(ip,j,k)*gsup(ip,j,k,1,2)
     .                 +dff(i,j ,k)*gsup(i,j ,k,1,2)
     .                 +dff(i,jm,k)*gsup(i,jm,k,1,2))
      d_xy_imjp = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,2)
     .                 +dff(im,j,k)*gsup(im,j,k,1,2)
     .                 +dff(i,j ,k)*gsup(i,j ,k,1,2)
     .                 +dff(i,jp,k)*gsup(i,jp,k,1,2))
      d_xy_imjm = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,2)
     .                 +dff(im,j,k)*gsup(im,j,k,1,2)
     .                 +dff(i,j ,k)*gsup(i,j ,k,1,2)
     .                 +dff(i,jm,k)*gsup(i,jm,k,1,2))

      d_xz_ipkp = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,3)
     .                 +dff(ip,j,k)*gsup(ip,j,k,1,3)
     .                 +dff(i,j,k )*gsup(i,j,k ,1,3)
     .                 +dff(i,j,kp)*gsup(i,j,kp,1,3))
      d_xz_ipkm = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,3)
     .                 +dff(ip,j,k)*gsup(ip,j,k,1,3)
     .                 +dff(i,j,k )*gsup(i,j,k ,1,3)
     .                 +dff(i,j,km)*gsup(i,j,km,1,3))
      d_xz_imkp = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,3)
     .                 +dff(ip,j,k)*gsup(ip,j,k,1,3)
     .                 +dff(i,j,k )*gsup(i,j,k ,1,3)
     .                 +dff(i,j,km)*gsup(i,j,km,1,3))
      d_xz_imkm = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,3)
     .                 +dff(im,j,k)*gsup(im,j,k,1,3)
     .                 +dff(i,j,k )*gsup(i,j,k ,1,3)
     .                 +dff(i,j,km)*gsup(i,j,km,1,3))

      d_yz_jpkp = 0.25*(dff(i,j ,k)*gsup(i,j ,k,2,3)
     .                 +dff(i,jp,k)*gsup(i,jp,k,2,3)
     .                 +dff(i,j,k )*gsup(i,j,k ,2,3)
     .                 +dff(i,j,kp)*gsup(i,j,kp,2,3))
      d_yz_jpkm = 0.25*(dff(i,j ,k)*gsup(i,j ,k,2,3)
     .                 +dff(i,jp,k)*gsup(i,jp,k,2,3)
     .                 +dff(i,j,k )*gsup(i,j,k ,2,3)
     .                 +dff(i,j,km)*gsup(i,j,km,2,3))
      d_yz_jmkp = 0.25*(dff(i,j ,k)*gsup(i,j ,k,2,3)
     .                 +dff(i,jm,k)*gsup(i,jm,k,2,3)
     .                 +dff(i,j,k )*gsup(i,j,k ,2,3)
     .                 +dff(i,j,kp)*gsup(i,j,kp,2,3))
      d_yz_jmkm = 0.25*(dff(i,j ,k)*gsup(i,j ,k,2,3)
     .                 +dff(i,jm,k)*gsup(i,jm,k,2,3)
     .                 +dff(i,j,k )*gsup(i,j,k ,2,3)
     .                 +dff(i,j,km)*gsup(i,j,km,2,3))

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

      end function lap_diff

c     lap_ndiff
c     ###############################################################
      function lap_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,arr,vol)
     .         result (laplacian)
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

      real(8),dimension(:,:,:,:,:),pointer :: gsup

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
      endif

      gsup => gmetric%grid(igrid)%gsup

      d_xx_ip = 0.5*(gsup(i ,j,k,1,1)
     .              +gsup(ip,j,k,1,1))
      d_xx_im = 0.5*(gsup(i ,j,k,1,1)
     .              +gsup(im,j,k,1,1))
      d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
     .              +gsup(i,jp,k,2,2))
      d_yy_jm = 0.5*(gsup(i,j ,k,2,2)
     .              +gsup(i,jm,k,2,2))
      d_zz_kp = 0.5*(gsup(i,j,k ,3,3)
     .              +gsup(i,j,kp,3,3))
      d_zz_km = 0.5*(gsup(i,j,k ,3,3)
     .              +gsup(i,j,km,3,3))

      d_xy_ipjp = 0.25*(gsup(i ,j,k,1,2)
     .                 +gsup(ip,j,k,1,2)
     .                 +gsup(i,j ,k,1,2)
     .                 +gsup(i,jp,k,1,2))
      d_xy_ipjm = 0.25*(gsup(i ,j,k,1,2)
     .                 +gsup(ip,j,k,1,2)
     .                 +gsup(i,j ,k,1,2)
     .                 +gsup(i,jm,k,1,2))
      d_xy_imjp = 0.25*(gsup(i ,j,k,1,2)
     .                 +gsup(im,j,k,1,2)
     .                 +gsup(i,j ,k,1,2)
     .                 +gsup(i,jp,k,1,2))
      d_xy_imjm = 0.25*(gsup(i ,j,k,1,2)
     .                 +gsup(im,j,k,1,2)
     .                 +gsup(i,j ,k,1,2)
     .                 +gsup(i,jm,k,1,2))

      d_xz_ipkp = 0.25*(gsup(i ,j,k,1,3)
     .                 +gsup(ip,j,k,1,3)
     .                 +gsup(i,j,k ,1,3)
     .                 +gsup(i,j,kp,1,3))
      d_xz_ipkm = 0.25*(gsup(i ,j,k,1,3)
     .                 +gsup(ip,j,k,1,3)
     .                 +gsup(i,j,k ,1,3)
     .                 +gsup(i,j,km,1,3))
      d_xz_imkp = 0.25*(gsup(i ,j,k,1,3)
     .                 +gsup(ip,j,k,1,3)
     .                 +gsup(i,j,k ,1,3)
     .                 +gsup(i,j,km,1,3))
      d_xz_imkm = 0.25*(gsup(i ,j,k,1,3)
     .                 +gsup(im,j,k,1,3)
     .                 +gsup(i,j,k ,1,3)
     .                 +gsup(i,j,km,1,3))

      d_yz_jpkp = 0.25*(gsup(i,j ,k,2,3)
     .                 +gsup(i,jp,k,2,3)
     .                 +gsup(i,j,k ,2,3)
     .                 +gsup(i,j,kp,2,3))
      d_yz_jpkm = 0.25*(gsup(i,j ,k,2,3)
     .                 +gsup(i,jp,k,2,3)
     .                 +gsup(i,j,k ,2,3)
     .                 +gsup(i,j,km,2,3))
      d_yz_jmkp = 0.25*(gsup(i,j ,k,2,3)
     .                 +gsup(i,jm,k,2,3)
     .                 +gsup(i,j,k ,2,3)
     .                 +gsup(i,j,kp,2,3))
      d_yz_jmkm = 0.25*(gsup(i,j ,k,2,3)
     .                 +gsup(i,jm,k,2,3)
     .                 +gsup(i,j,k ,2,3)
     .                 +gsup(i,j,km,2,3))

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

      end function lap_ndiff

ccc     laplacian
ccc     ###############################################################
cc      function laplacian(i,j,k,nx,ny,nz,igx,igy,igz,arr,vol,diff)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates lap(arr) at cell centers in general non-orthog.
ccc     coordinates, preserving the SPD property.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian
cc
cc      logical,optional,intent(IN) :: vol
cc
cc      real(8),optional :: diff(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid
cc
cc      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
cc     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
cc     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
cc     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm
cc
cc      real(8),target :: dff(0:nx+1,0:ny+1,0:nz+1)
cc
cc      real(8),dimension(:,:,:,:,:),pointer :: gsup
cc
cc      logical    :: vol_wgt
cc
ccc     Begin program
cc
cc      vol_wgt = .true.
cc      if (PRESENT(vol)) vol_wgt = vol
cc
cc      igrid = igx
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
cc      if (i == 0 .or. j == 0 .or. k == 0 ) then
cc        write (*,*) 'Error in laplace; i,j,k=0'
cc      elseif (i == nx+1 .or. j == ny+1 .or. k == nz+1) then
cc        write (*,*) 'Error in laplace; i,j,k=nmax+1'
cc      endif
cc
cc      gsup => gmetric%grid(igrid)%gsup
cc
cc      if (PRESENT(diff)) then
cc        dff = diff
cc
cc        d_xx_ip = 0.5*(dff(i ,j,k)*gsup(i ,j,k,1,1)
cc     .                +dff(ip,j,k)*gsup(ip,j,k,1,1))
cc        d_xx_im = 0.5*(dff(i ,j,k)*gsup(i ,j,k,1,1)
cc     .                +dff(im,j,k)*gsup(im,j,k,1,1))
cc        d_yy_jp = 0.5*(dff(i,j ,k)*gsup(i,j ,k,2,2)
cc     .                +dff(i,jp,k)*gsup(i,jp,k,2,2))
cc        d_yy_jm = 0.5*(dff(i,j ,k)*gsup(i,j ,k,2,2)
cc     .                +dff(i,jm,k)*gsup(i,jm,k,2,2))
cc        d_zz_kp = 0.5*(dff(i,j,k )*gsup(i,j,k ,3,3)
cc     .                +dff(i,j,kp)*gsup(i,j,kp,3,3))
cc        d_zz_km = 0.5*(dff(i,j,k )*gsup(i,j,k ,3,3)
cc     .                +dff(i,j,km)*gsup(i,j,km,3,3))
cc
cc        d_xy_ipjp = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,2)
cc     .                   +dff(ip,j,k)*gsup(ip,j,k,1,2)
cc     .                   +dff(i,j ,k)*gsup(i,j ,k,1,2)
cc     .                   +dff(i,jp,k)*gsup(i,jp,k,1,2))
cc        d_xy_ipjm = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,2)
cc     .                   +dff(ip,j,k)*gsup(ip,j,k,1,2)
cc     .                   +dff(i,j ,k)*gsup(i,j ,k,1,2)
cc     .                   +dff(i,jm,k)*gsup(i,jm,k,1,2))
cc        d_xy_imjp = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,2)
cc     .                   +dff(im,j,k)*gsup(im,j,k,1,2)
cc     .                   +dff(i,j ,k)*gsup(i,j ,k,1,2)
cc     .                   +dff(i,jp,k)*gsup(i,jp,k,1,2))
cc        d_xy_imjm = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,2)
cc     .                   +dff(im,j,k)*gsup(im,j,k,1,2)
cc     .                   +dff(i,j ,k)*gsup(i,j ,k,1,2)
cc     .                   +dff(i,jm,k)*gsup(i,jm,k,1,2))
cc
cc        d_xz_ipkp = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,3)
cc     .                   +dff(ip,j,k)*gsup(ip,j,k,1,3)
cc     .                   +dff(i,j,k )*gsup(i,j,k ,1,3)
cc     .                   +dff(i,j,kp)*gsup(i,j,kp,1,3))
cc        d_xz_ipkm = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,3)
cc     .                   +dff(ip,j,k)*gsup(ip,j,k,1,3)
cc     .                   +dff(i,j,k )*gsup(i,j,k ,1,3)
cc     .                   +dff(i,j,km)*gsup(i,j,km,1,3))
cc        d_xz_imkp = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,3)
cc     .                   +dff(ip,j,k)*gsup(ip,j,k,1,3)
cc     .                   +dff(i,j,k )*gsup(i,j,k ,1,3)
cc     .                   +dff(i,j,km)*gsup(i,j,km,1,3))
cc        d_xz_imkm = 0.25*(dff(i ,j,k)*gsup(i ,j,k,1,3)
cc     .                   +dff(im,j,k)*gsup(im,j,k,1,3)
cc     .                   +dff(i,j,k )*gsup(i,j,k ,1,3)
cc     .                   +dff(i,j,km)*gsup(i,j,km,1,3))
cc
cc        d_yz_jpkp = 0.25*(dff(i,j ,k)*gsup(i,j ,k,2,3)
cc     .                   +dff(i,jp,k)*gsup(i,jp,k,2,3)
cc     .                   +dff(i,j,k )*gsup(i,j,k ,2,3)
cc     .                   +dff(i,j,kp)*gsup(i,j,kp,2,3))
cc        d_yz_jpkm = 0.25*(dff(i,j ,k)*gsup(i,j ,k,2,3)
cc     .                   +dff(i,jp,k)*gsup(i,jp,k,2,3)
cc     .                   +dff(i,j,k )*gsup(i,j,k ,2,3)
cc     .                   +dff(i,j,km)*gsup(i,j,km,2,3))
cc        d_yz_jmkp = 0.25*(dff(i,j ,k)*gsup(i,j ,k,2,3)
cc     .                   +dff(i,jm,k)*gsup(i,jm,k,2,3)
cc     .                   +dff(i,j,k )*gsup(i,j,k ,2,3)
cc     .                   +dff(i,j,kp)*gsup(i,j,kp,2,3))
cc        d_yz_jmkm = 0.25*(dff(i,j ,k)*gsup(i,j ,k,2,3)
cc     .                   +dff(i,jm,k)*gsup(i,jm,k,2,3)
cc     .                   +dff(i,j,k )*gsup(i,j,k ,2,3)
cc     .                   +dff(i,j,km)*gsup(i,j,km,2,3))
cc      else
cc        d_xx_ip = 0.5*(gsup(i ,j,k,1,1)
cc     .                +gsup(ip,j,k,1,1))
cc        d_xx_im = 0.5*(gsup(i ,j,k,1,1)
cc     .                +gsup(im,j,k,1,1))
cc        d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
cc     .                +gsup(i,jp,k,2,2))
cc        d_yy_jm = 0.5*(gsup(i,j ,k,2,2)
cc     .                +gsup(i,jm,k,2,2))
cc        d_zz_kp = 0.5*(gsup(i,j,k ,3,3)
cc     .                +gsup(i,j,kp,3,3))
cc        d_zz_km = 0.5*(gsup(i,j,k ,3,3)
cc     .                +gsup(i,j,km,3,3))
cc
cc        d_xy_ipjp = 0.25*(gsup(i ,j,k,1,2)
cc     .                   +gsup(ip,j,k,1,2)
cc     .                   +gsup(i,j ,k,1,2)
cc     .                   +gsup(i,jp,k,1,2))
cc        d_xy_ipjm = 0.25*(gsup(i ,j,k,1,2)
cc     .                   +gsup(ip,j,k,1,2)
cc     .                   +gsup(i,j ,k,1,2)
cc     .                   +gsup(i,jm,k,1,2))
cc        d_xy_imjp = 0.25*(gsup(i ,j,k,1,2)
cc     .                   +gsup(im,j,k,1,2)
cc     .                   +gsup(i,j ,k,1,2)
cc     .                   +gsup(i,jp,k,1,2))
cc        d_xy_imjm = 0.25*(gsup(i ,j,k,1,2)
cc     .                   +gsup(im,j,k,1,2)
cc     .                   +gsup(i,j ,k,1,2)
cc     .                   +gsup(i,jm,k,1,2))
cc
cc        d_xz_ipkp = 0.25*(gsup(i ,j,k,1,3)
cc     .                   +gsup(ip,j,k,1,3)
cc     .                   +gsup(i,j,k ,1,3)
cc     .                   +gsup(i,j,kp,1,3))
cc        d_xz_ipkm = 0.25*(gsup(i ,j,k,1,3)
cc     .                   +gsup(ip,j,k,1,3)
cc     .                   +gsup(i,j,k ,1,3)
cc     .                   +gsup(i,j,km,1,3))
cc        d_xz_imkp = 0.25*(gsup(i ,j,k,1,3)
cc     .                   +gsup(ip,j,k,1,3)
cc     .                   +gsup(i,j,k ,1,3)
cc     .                   +gsup(i,j,km,1,3))
cc        d_xz_imkm = 0.25*(gsup(i ,j,k,1,3)
cc     .                   +gsup(im,j,k,1,3)
cc     .                   +gsup(i,j,k ,1,3)
cc     .                   +gsup(i,j,km,1,3))
cc
cc        d_yz_jpkp = 0.25*(gsup(i,j ,k,2,3)
cc     .                   +gsup(i,jp,k,2,3)
cc     .                   +gsup(i,j,k ,2,3)
cc     .                   +gsup(i,j,kp,2,3))
cc        d_yz_jpkm = 0.25*(gsup(i,j ,k,2,3)
cc     .                   +gsup(i,jp,k,2,3)
cc     .                   +gsup(i,j,k ,2,3)
cc     .                   +gsup(i,j,km,2,3))
cc        d_yz_jmkp = 0.25*(gsup(i,j ,k,2,3)
cc     .                   +gsup(i,jm,k,2,3)
cc     .                   +gsup(i,j,k ,2,3)
cc     .                   +gsup(i,j,kp,2,3))
cc        d_yz_jmkm = 0.25*(gsup(i,j ,k,2,3)
cc     .                   +gsup(i,jm,k,2,3)
cc     .                   +gsup(i,j,k ,2,3)
cc     .                   +gsup(i,j,km,2,3))
cc      endif
cc
cc      laplacian = 
cc     .     dyh(jg)*dzh(kg)*( d_xx_ip*(arr(ip,j,k)-arr(i,j,k))/dx(ig)
cc     .                      +d_xx_im*(arr(im,j,k)-arr(i,j,k))/dx(ig-1) )
cc     .    +dxh(ig)*dzh(kg)*( d_yy_jp*(arr(i,jp,k)-arr(i,j,k))/dy(jg)
cc     .                      +d_yy_jm*(arr(i,jm,k)-arr(i,j,k))/dy(jg-1) )
cc     .    +dxh(ig)*dyh(jg)*( d_zz_kp*(arr(i,j,kp)-arr(i,j,k))/dz(kg)
cc     .                      +d_zz_km*(arr(i,j,km)-arr(i,j,k))/dz(kg-1) )
cc     .    +0.5*( dzh(kg)*( d_xy_ipjp*(arr(ip,jp,k)-arr(i,j,k))
cc     .                    +d_xy_imjm*(arr(im,jm,k)-arr(i,j,k))
cc     .                    -d_xy_ipjm*(arr(ip,jm,k)-arr(i,j,k))
cc     .                    -d_xy_imjp*(arr(im,jp,k)-arr(i,j,k)) )
cc     .          +dyh(jg)*( d_xz_ipkp*(arr(ip,j,kp)-arr(i,j,k))
cc     .                    +d_xz_imkm*(arr(im,j,km)-arr(i,j,k))
cc     .                    -d_xz_ipkm*(arr(ip,j,km)-arr(i,j,k))
cc     .                    -d_xz_imkp*(arr(im,j,kp)-arr(i,j,k)) )
cc     .          +dxh(ig)*( d_yz_jpkp*(arr(i,jp,kp)-arr(i,j,k))
cc     .                    +d_yz_jmkm*(arr(i,jm,km)-arr(i,j,k))
cc     .                    -d_yz_jpkm*(arr(i,jp,km)-arr(i,j,k))
cc     .                    -d_yz_jmkp*(arr(i,jm,kp)-arr(i,j,k)) ) )
cc
cc      if (.not.vol_wgt) laplacian=laplacian/volume(i,j,k,igx,igy,igz)
cc
ccc     End program
cc
cc      end function laplacian

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

      real(8)    :: jac,ijac,dvol,dS1,dS2,dS3,dxx,dyy,dzz

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o

      real(8)    :: hess(3,3,3),msource,dx1,dx2,ll

      logical    :: vol_wgt

      real(8)    :: dum1,dum2,coeff

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
      ijac = 1./jac

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      if (coords /= 'car') then
        hess = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
      endif

      call tsrx(i ,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11p,t12p,t13p, 1)
      call tsrx(im,j,k,nx,ny,nz,igx,igy,igz,alt_eom,t11m,t12m,t13m,-1)
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
      if (coords /= 'car') then
        msource =  (t11o*hess(1,1,1)+t12o*hess(1,1,2)+t13o*hess(1,1,3)
     .             +t21o*hess(1,2,1)+t22o*hess(1,2,2)+t23o*hess(1,2,3)
     .             +t31o*hess(1,3,1)+t32o*hess(1,3,2)+t33o*hess(1,3,3))
     .             *ijac
      endif

      divt(1) =  (t11p - t11m)/dxx
     .          +(t21p - t21m)/dyy
     .          +(t31p - t31m)/dzz + msource

      !Component 2
      coeff = 1d0

      if (coords /= 'car') then
        if (alt_eom) then

          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3)
     .              -t12o*(hess(1,1,1)+hess(2,1,2)+hess(3,1,3))
     .              -t22o*(hess(1,2,1)+hess(2,2,2)+hess(3,2,3))
     .              -t32o*(hess(1,3,1)+hess(2,3,2)+hess(3,3,3)))*ijac

          coeff = ijac

        else
          msource=  (t11o*hess(2,1,1)+t12o*hess(2,1,2)+t13o*hess(2,1,3)
     .              +t21o*hess(2,2,1)+t22o*hess(2,2,2)+t23o*hess(2,2,3)
     .              +t31o*hess(2,3,1)+t32o*hess(2,3,2)+t33o*hess(2,3,3))
     .              *ijac
        endif

      endif

      divt(2) = ( (t12p - t12m)/dxx
     .           +(t22p - t22m)/dyy
     .           +(t32p - t32m)/dzz)*coeff  + msource

      !Component 3
      if (coords /= 'car') then
        msource = (t11o*hess(3,1,1)+t12o*hess(3,1,2)+t13o*hess(3,1,3)
     .            +t21o*hess(3,2,1)+t22o*hess(3,2,2)+t23o*hess(3,2,3)
     .            +t31o*hess(3,3,1)+t32o*hess(3,3,2)+t33o*hess(3,3,3))
     .            *ijac
      endif

      divt(3) =  (t13p - t13m)/dxx
     .          +(t23p - t23m)/dyy
     .          +(t33p - t33m)/dzz + msource

      !Volume factor
cc      if (vol_wgt) divt=divt*volume(i,j,k,igx,igy,igz)
      if (vol_wgt) divt=divt*gmetric%grid(igrid)%dvol(i,j,k)

c     End program

      end function div_tensor

c     veclap_diff
c     ###############################################################
      function veclap_diff(i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                     ,alteom,diff,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*lap(vector) at cell centers in general non-orthog.
c     coordinates. Vector is assumed in contravariant representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,diff   (0:nx+1,0:ny+1,0:nz+1)

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
     .                           ,vfield(:,:,:,icomp),diff
     .                           ,vol=vol_wgt)
        enddo
        return
      else
        vec  => vfield !Pointer passed to nabtensor routines
        coef => diff   !Pointer passed to nabtensor routines
        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
     $                   ,vol=vol_wgt)
      endif

c     End program

      end function veclap_diff

c     veclap_ndiff
c     ###############################################################
      function veclap_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,vfield
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
     .                           ,vfield(:,:,:,icomp)
     .                           ,vol=vol_wgt)
        enddo
        return
      else
        vec  => vfield !Pointer passed to nabtensor routines
        nullify(coef)
        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
     $                   ,vol=vol_wgt)
      endif

c     End program

      end function veclap_ndiff

ccc     veclaplacian
ccc     ###############################################################
cc      function veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,vfield
cc     .                     ,alteom,vol,diff) result (vlap)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates dvol*lap(vector) at cell centers in general non-orthog.
ccc     coordinates. Vector is assumed in contravariant representation.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
cc
cc      real(8)    :: vlap(3)
cc
cc      logical    :: alteom
cc
cc      logical,optional,intent(IN) :: vol
cc
cc      real(8),optional :: diff(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
cc      integer(4) :: icomp
cc
cc      real(8),target :: dff(0:nx+1,0:ny+1,0:nz+1)
cc
cc      logical    :: vol_wgt
cc
ccc     Begin program
cc
cc      vol_wgt = .true.
cc      if (PRESENT(vol)) vol_wgt = vol
cc
cc      if (coords == 'car') then
cc        if (PRESENT(diff)) then
cc          dff = diff
cc          do icomp=1,3
cc            vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                           ,vfield(:,:,:,icomp)
cc     .                           ,vol=vol_wgt,diff=dff)
cc          enddo
cc        else
cc          do icomp=1,3
cc            vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                           ,vfield(:,:,:,icomp)
cc     .                           ,vol=vol_wgt)
cc          enddo
cc        endif
cc        return
cc      else
cc        if (PRESENT(diff)) then
cc          dff = diff
cc        else
cc          dff = 1d0
cc        endif
cc        vec  => vfield !Pointer passed to nabtensor routines
cc        coef => dff    !Pointer passed to nabtensor routines
cc        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alteom
cc     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
cc     $                   ,vol=vol_wgt)
cc      endif
cc
ccc     End program
cc
cc      end function veclaplacian

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
        real(8)    :: x,y,z,jac,jac0,jacp,ijac
        real(8)    :: nabla_v(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
        ijac = 1d0/jac

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),1)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),0)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(ip,j,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t11 =dd*( gsuper(1,1)*nabla_v(1,1)
     .           +gsuper(1,2)*nabla_v(2,1)
     .           +gsuper(1,3)*nabla_v(3,1) )

        t12 =dd*( gsuper(1,1)*nabla_v(1,2)
     .           +gsuper(1,2)*nabla_v(2,2)
     .           +gsuper(1,3)*nabla_v(3,2) )

        t13 =dd*( gsuper(1,1)*nabla_v(1,3)
     .           +gsuper(1,2)*nabla_v(2,3)
     .           +gsuper(1,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t11 = t11*ijac
          if (.not.alteom) t12 = t12*ijac
          t13 = t13*ijac
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
        real(8)    :: x,y,z,jac,ijac
        real(8)    :: nabla_v(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        ijac = 1d0/jac

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),2)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),0)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,jp,k) + 1./coef(i,j ,k))
        else
          dd = 1d0
        endif

        t21 =dd*( gsuper(2,1)*nabla_v(1,1)
     .           +gsuper(2,2)*nabla_v(2,1)
     .           +gsuper(2,3)*nabla_v(3,1) )

        t22 =dd*( gsuper(2,1)*nabla_v(1,2)
     .           +gsuper(2,2)*nabla_v(2,2)
     .           +gsuper(2,3)*nabla_v(3,2) )

        t23 =dd*( gsuper(2,1)*nabla_v(1,3)
     .           +gsuper(2,2)*nabla_v(2,3)
     .           +gsuper(2,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t21 = t21*ijac
          if (.not.alteom) t22 = t22*ijac
          t23 = t23*ijac
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
        real(8)    :: x,y,z,jac,ijac
        real(8)    :: nabla_v(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        ijac = 1d0/jac

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),3)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vec(:,:,:,1)
     .                                                 ,vec(:,:,:,2)
     .                                                 ,vec(:,:,:,3),0)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,j,kp) + 1./coef(i,j,k ))
        else
          dd = 1d0
        endif

        t31 =dd*( gsuper(3,1)*nabla_v(1,1)
     .           +gsuper(3,2)*nabla_v(2,1)
     .           +gsuper(3,3)*nabla_v(3,1) )

        t32 =dd*( gsuper(3,1)*nabla_v(1,2)
     .           +gsuper(3,2)*nabla_v(2,2)
     .           +gsuper(3,3)*nabla_v(3,2) )

        t33 =dd*( gsuper(3,1)*nabla_v(1,3)
     .           +gsuper(3,2)*nabla_v(2,3)
     .           +gsuper(3,3)*nabla_v(3,3) )

        if (flag /= 0) then
          t31 = t31*ijac
          if (.not.alteom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine nabtensor_z

c     grad
c     ###############################################################
      subroutine grad(i,j,k,nx,ny,nz,igx,igy,igz,arr,gx,gy,gz)

c     ---------------------------------------------------------------
c     Calculates grad(A)) in general non-orthogonal coordinates,
c     preserving the SPD property. The vector grad(A) is covariant.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer(4) :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),gx,gy,gz

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

      !X comp
      if (i==0) then
        dx1=grid_params%dx(ig)
        dx2=grid_params%dx(ig+1)
        gx = (-arr(i+2,j,k)*dx1/dx2/(dx1+dx2)
     .        +arr(i+1,j,k)*(dx1+dx2)/dx1/dx2
     .        -arr(i  ,j,k)*(1./dx1+1./(dx1+dx2)))
c1st        gx = (arr(ip,j,k)-arr(i,j,k))/dx1
      elseif (i==nx+1) then
        dx1=grid_params%dx(ig-1)
        dx2=grid_params%dx(ig-2)
        gx = -(-arr(i-2,j,k)*dx1/dx2/(dx1+dx2)
     .         +arr(i-1,j,k)*(dx1+dx2)/dx1/dx2
     .         -arr(i  ,j,k)*(1./dx1+1./(dx1+dx2)))
c1st        gx = (arr(i,j,k)-arr(im,j,k))/dx1
cc      elseif (sing_point) then
cc        !Second order formula
cc        dx1=grid_params%dx(ig-1)
cc        dx2=grid_params%dx(ig  )
cc        ll = dx1/dx2
cc        ar = 0.5*(arr(ip,j,k)*ll + arr(i,j,k)*(1.+1./ll))
cc        al = 0.5*(arr(im,j,k)/ll + arr(i,j,k)*(1.+   ll))
cc        gx = 2.*(ar-al)/(dx1+dx2)
cccc        dx = grid_params%dxh(ig)
cccc        gx = (0.5*(arr(ip,j,k)+arr(i,j,k))-arr(im,j,k))/dx
      else
        dx = grid_params%dxh(ig)
        gx = 0.5*(arr(ip,j,k)-arr(im,j,k))/dx
      endif

      !Y comp
      if (j==0) then
        dy1=grid_params%dy(jg)
        dy2=grid_params%dy(jg+1)
        gy = (-arr(i,j+2,k)*dy1/dy2/(dy1+dy2)
     .        +arr(i,j+1,k)*(dy1+dy2)/dy1/dy2
     .        -arr(i,j  ,k)*(1./dy1+1./(dy1+dy2)))
c1st        gy = (arr(i,jp,k)-arr(i,j,k))/dy1
      elseif (j==ny+1) then
        dy1=grid_params%dy(jg-1)
        dy2=grid_params%dy(jg-2)
        gy = -(-arr(i,j-2,k)*dy1/dy2/(dy1+dy2)
     .         +arr(i,j-1,k)*(dy1+dy2)/dy1/dy2
     .         -arr(i,j  ,k)*(1./dy1+1./(dy1+dy2)))
c1st        gy = (arr(i,j,k)-arr(i,jm,k))/dy1
      else
        dy = grid_params%dyh(jg)
        gy = 0.5*(arr(i,jp,k)-arr(i,jm,k))/dy
      endif

      !Z comp
      if (k==0) then
        dz1=grid_params%dz(kg)
        dz2=grid_params%dz(kg+1)
        gz = (-arr(i,j,k+2)*dz1/dz2/(dz1+dz2)
     .        +arr(i,j,k+1)*(dz1+dz2)/dz1/dz2
     .        -arr(i,j,k  )*(1./dz1+1./(dz1+dz2)))
c1st        gz = (arr(i,j,kp)-arr(i,j,k))/dz1
      elseif (k==nz+1) then
        dz1=grid_params%dz(kg-1)
        dz2=grid_params%dz(kg-2)
        gz = -(-arr(i,j,k-2)*dz1/dz2/(dz1+dz2)
     .         +arr(i,j,k-1)*(dz1+dz2)/dz1/dz2
     .         -arr(i,j,k  )*(1./dz1+1./(dz1+dz2)))
c1st        gz = (arr(i,j,k)-arr(i,j,km))/dz1
      else
        dz = grid_params%dzh(kg)
        gz = 0.5*(arr(i,j,kp)-arr(i,j,km))/dz
      endif

c     End program

      end subroutine grad

c     curl
c     ###############################################################
      function curl(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az,he) result(crl)

c     ---------------------------------------------------------------
c     Calculates curl(A) at cell centers in general non-orthogonal
c     coordinates. The vector components (ax,ay,az) are covariant.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8)    :: crl(3)

      integer(4) :: i,j,k,comp,nx,ny,nz,igx,igy,igz
      integer(4), optional :: he

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid,half_elem
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
     .             ,idhx,idhy,idhz
      real(8)    :: dS1,dS2,dS3,x0,y0,z0,jac,dx1,dx2,ll
     .             ,xip,yip,zip,jacp,xh,yh,zh,jach

c     Begin program

      if (PRESENT(he)) then
        half_elem = he
      else
        half_elem = 0
      endif

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      idhx = 1d0/grid_params%dxh(ig)
      idhy = 1d0/grid_params%dyh(jg)
      idhz = 1d0/grid_params%dzh(kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      select case(half_elem)
      case (1)
        idhx = 1./grid_params%dx(ig)

        !X comp

        flxjp = 0.25*(az(i ,jp,k)+az(i ,j,k)
     .               +az(ip,jp,k)+az(ip,j,k))
        flxjm = 0.25*(az(i ,jm,k)+az(i ,j,k)
     .               +az(ip,jm,k)+az(ip,j,k))

        flxkp =-0.25*(ay(i ,j,kp)+ay(i ,j,k)
     .               +ay(ip,j,kp)+ay(ip,j,k))
        flxkm =-0.25*(ay(i ,j,km)+ay(i ,j,k)
     .               +ay(ip,j,km)+ay(ip,j,k))

        crl(1) = (flxjp-flxjm)*idhy
     .          +(flxkp-flxkm)*idhz

        !Y comp

        flxip =-az(ip,j,k)
        flxim =-az(i ,j,k)

        flxkp = 0.25*(ax(i ,j,kp)+ax(i ,j,k)
     .               +ax(ip,j,kp)+ax(ip,j,k))
        flxkm = 0.25*(ax(i ,j,km)+ax(i ,j,k)
     .               +ax(ip,j,km)+ax(ip,j,k))

        crl(2) = (flxip-flxim)*idhx
     .          +(flxkp-flxkm)*idhz

        !Z comp

        flxip = ay(ip,j,k)
        flxim = ay(i ,j,k)

        flxjp =-0.25*(ax(i ,jp,k)+ax(i ,j,k)
     .               +ax(ip,jp,k)+ax(ip,j,k))
        flxjm =-0.25*(ax(i ,jm,k)+ax(i ,j,k)
     .               +ax(ip,jm,k)+ax(ip,j,k))

        crl(3) = (flxip-flxim)*idhx
     .          +(flxjp-flxjm)*idhy

      case (2)

        idhy = 1./grid_params%dy(jg)

        !X comp

        flxjp = az(i,jp,k)
        flxjm = az(i,j ,k)

        flxkp =-0.25*(ay(i,j ,kp)+ay(i,j ,k)
     .               +ay(i,jp,kp)+ay(i,jp,k))
        flxkm =-0.25*(ay(i,j ,km)+ay(i,j ,k)
     .               +ay(i,jp,km)+ay(i,jp,k))

        crl(1) = (flxjp-flxjm)*idhy
     .          +(flxkp-flxkm)*idhz

        !Y comp

        flxip =-0.25*(az(ip,j ,k)+az(i,j ,k)
     .               +az(ip,jp,k)+az(i,jp,k))
        flxim =-0.25*(az(im,j ,k)+az(i,j ,k)
     .               +az(im,jp,k)+az(i,jp,k))

        flxkp = 0.25*(ax(i,j ,kp)+ax(i,j ,k)
     .               +ax(i,jp,kp)+ax(i,jp,k))
        flxkm = 0.25*(ax(i,j ,km)+ax(i,j ,k)
     .               +ax(i,jp,km)+ax(i,jp,k))

        crl(2) = (flxip-flxim)*idhx
     .          +(flxkp-flxkm)*idhz

        !Z comp

        flxip = 0.25*(ay(ip,j ,k)+ay(i,j ,k)
     .               +ay(ip,jp,k)+ay(i,jp,k))
        flxim = 0.25*(ay(im,j ,k)+ay(i,j ,k)
     .               +ay(im,jp,k)+ay(i,jp,k))

        flxjp =-ax(i,jp,k)
        flxjm =-ax(i,j ,k)

        crl(3) = (flxip-flxim)*idhx
     .          +(flxjp-flxjm)*idhy

      case (3)

        idhz = 1./grid_params%dz(kg)

        !X comp

        flxjp = 0.25*(az(i,jp,k )+az(i,j,k )
     .               +az(i,jp,kp)+az(i,j,kp))
        flxjm = 0.25*(az(i,jm,k )+az(i,j,k )
     .               +az(i,jm,kp)+az(i,j,kp))

        flxkp =-ay(i,j,kp)
        flxkm =-ay(i,j,k )

        crl(1) = (flxjp-flxjm)*idhy
     .          +(flxkp-flxkm)*idhz

        !Y comp

        flxip =-0.25*(az(ip,j,k )+az(i,j,k )
     .               +az(ip,j,kp)+az(i,j,kp))
        flxim =-0.25*(az(im,j,k )+az(i,j,k )
     .               +az(im,j,kp)+az(i,j,kp))

        flxkp = ax(i,j,kp)
        flxkm = ax(i,j,k )

        crl(2) = (flxip-flxim)*idhx
     .          +(flxkp-flxkm)*idhz

        !Z comp

        flxip = 0.25*(ay(ip,j,k )+ay(i,j,k )
     .               +ay(ip,j,kp)+ay(i,j,kp))
        flxim = 0.25*(ay(im,j,k )+ay(i,j,k )
     .               +ay(im,j,kp)+ay(i,j,kp))

        flxjp =-0.25*(ax(i,jp,k )+ax(i,j,k )
     .               +ax(i,jp,kp)+ax(i,j,kp))
        flxjm =-0.25*(ax(i,jm,k )+ax(i,j,k )
     .               +ax(i,jm,kp)+ax(i,j,kp))

        crl(3) = (flxip-flxim)*idhx
     .          +(flxjp-flxjm)*idhy

      case default

        !X comp

        flxjp = 0.5*(az(i,jp,k)+az(i,j,k))
        flxjm = 0.5*(az(i,jm,k)+az(i,j,k))

        flxkp =-0.5*(ay(i,j,kp)+ay(i,j,k))
        flxkm =-0.5*(ay(i,j,km)+ay(i,j,k))

        crl(1) = (flxjp-flxjm)*idhy
     .          +(flxkp-flxkm)*idhz

        !Y comp

        flxip =-0.5*(az(ip,j,k)+az(i,j,k))
        flxim =-0.5*(az(im,j,k)+az(i,j,k))

        flxkp = 0.5*(ax(i,j,kp)+ax(i,j,k))
        flxkm = 0.5*(ax(i,j,km)+ax(i,j,k))

        crl(2) = (flxip-flxim)*idhx
     .          +(flxkp-flxkm)*idhz

        !Z comp

        flxip = 0.5*(ay(ip,j,k)+ay(i,j,k))
        flxim = 0.5*(ay(im,j,k)+ay(i,j,k))

        flxjp =-0.5*(ax(i,jp,k)+ax(i,j,k))
        flxjm =-0.5*(ax(i,jm,k)+ax(i,j,k))

        crl(3) = (flxip-flxim)*idhx
     .          +(flxjp-flxjm)*idhy

      end select

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

c     curl_bxv
c     ###################################################################
      function curl_bxv(i,j,k,nx,ny,nz,igx,igy,igz,vv,bb,half_elem)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of -curl(vv x bb) at the
c     grid node (i,j,k). One sided derivatives are employed when half_elem=1
c     (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nx,ny,nz,half_elem,igx,igy,igz
        real(8)    :: cnv(3)
        real(8)    :: vv(0:nx+1,0:ny+1,0:nz+1,3)
     $               ,bb(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

        integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq
        integer(4) :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

        real(8)    :: idhx,idhy,idhz,a(3)
        real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .               ,jacp,jacm,jacph,jacmh,jach,jac0
     .               ,ijacip,ijacim,ijacjp,ijacjm,ijackp,ijackm
     .               ,ijacp,ijacm,ijacph,ijacmh,ijach,ijac0

        real(8)    :: vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm

        real(8)    :: bxip,bxim,bxjp,bxjm,bxkp,bxkm
     .               ,byip,byim,byjp,byjm,bykp,bykm
     .               ,bzip,bzim,bzjp,bzjm,bzkp,bzkm

c     Begin program

c     Defaults

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        idhx = 0.5/grid_params%dxh(ig)
        idhy = 0.5/grid_params%dyh(jg)
        idhz = 0.5/grid_params%dzh(kg)

        jac  = gmetric%grid(igx)%jac(i,j,k)

cc        sing_point = isSP(i,j,k,igx,igy,igz)

c     Exceptions

        select case(half_elem)
        case (1)
          idhx = 1./dx(ig)
          im = i

          jacip  = gmetric%grid(igx)%jac(ip,j,k)
          jacim  = gmetric%grid(igx)%jac(i ,j,k)
          jacjp  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
     .                 +gmetric%grid(igx)%jac(i ,jp,k))
          jacjm  = 0.5*(gmetric%grid(igx)%jac(ip,jm,k)
     .                 +gmetric%grid(igx)%jac(i ,jm,k))
          jackp  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
     .                 +gmetric%grid(igx)%jac(i ,j,kp))
          jackm  = 0.5*(gmetric%grid(igx)%jac(ip,j,km)
     .                 +gmetric%grid(igx)%jac(i ,j,km))

          if (isSP(i+1,j,k,igx,igy,igz)) then
            jacjp = SP_flsv
            jacjm = SP_flsv
            jackp = SP_flsv
            jackm = SP_flsv
cc            jacjp = 1d5!SP_flsv
cc            jacjm = 1d5!SP_flsv
cc            jackp = 1d5!SP_flsv
cc            jackm = 1d5!SP_flsv
          endif

          vxip = vv(ip,j,k,1)
          vxim = vv(i ,j,k,1)
          vyip = vv(ip,j,k,2)
          vyim = vv(i ,j,k,2)
          vzip = vv(ip,j,k,3)
          vzim = vv(i ,j,k,3)

          vxjp = 0.5*(vv(ip,jp,k,1)+vv(i,jp,k,1))
          vxjm = 0.5*(vv(ip,jm,k,1)+vv(i,jm,k,1))
          vyjp = 0.5*(vv(ip,jp,k,2)+vv(i,jp,k,2))
          vyjm = 0.5*(vv(ip,jm,k,2)+vv(i,jm,k,2))
          vzjp = 0.5*(vv(ip,jp,k,3)+vv(i,jp,k,3))
          vzjm = 0.5*(vv(ip,jm,k,3)+vv(i,jm,k,3))

          vxkp = 0.5*(vv(ip,j,kp,1)+vv(i,j,kp,1))
          vxkm = 0.5*(vv(ip,j,km,1)+vv(i,j,km,1))
          vykp = 0.5*(vv(ip,j,kp,2)+vv(i,j,kp,2))
          vykm = 0.5*(vv(ip,j,km,2)+vv(i,j,km,2))
          vzkp = 0.5*(vv(ip,j,kp,3)+vv(i,j,kp,3))
          vzkm = 0.5*(vv(ip,j,km,3)+vv(i,j,km,3))

          bxip = bb(ip,j,k,1)
          bxim = bb(i ,j,k,1)
          byip = bb(ip,j,k,2)
          byim = bb(i ,j,k,2)
          bzip = bb(ip,j,k,3)
          bzim = bb(i ,j,k,3)

          bxjp = 0.5*(bb(ip,jp,k,1)+bb(i,jp,k,1))
          bxjm = 0.5*(bb(ip,jm,k,1)+bb(i,jm,k,1))
          byjp = 0.5*(bb(ip,jp,k,2)+bb(i,jp,k,2))
          byjm = 0.5*(bb(ip,jm,k,2)+bb(i,jm,k,2))
          bzjp = 0.5*(bb(ip,jp,k,3)+bb(i,jp,k,3))
          bzjm = 0.5*(bb(ip,jm,k,3)+bb(i,jm,k,3))

          bxkp = 0.5*(bb(ip,j,kp,1)+bb(i,j,kp,1))
          bxkm = 0.5*(bb(ip,j,km,1)+bb(i,j,km,1))
          bykp = 0.5*(bb(ip,j,kp,2)+bb(i,j,kp,2))
          bykm = 0.5*(bb(ip,j,km,2)+bb(i,j,km,2))
          bzkp = 0.5*(bb(ip,j,kp,3)+bb(i,j,kp,3))
          bzkm = 0.5*(bb(ip,j,km,3)+bb(i,j,km,3))

        case (2)

          idhy = 1./dy(jg)
          jm = j

          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
     .                 +gmetric%grid(igx)%jac(ip,j ,k))
          jacim  = 0.5*(gmetric%grid(igx)%jac(im,jp,k)
     .                 +gmetric%grid(igx)%jac(im,j ,k))
          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
          jacjm  = gmetric%grid(igx)%jac(i,j ,k)
          jackp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
     .                 +gmetric%grid(igx)%jac(i,j ,kp))
          jackm  = 0.5*(gmetric%grid(igx)%jac(i,jp,km)
     .                 +gmetric%grid(igx)%jac(i,j ,km))

          vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1))*0.5
          vxim = (vv(im,j,k,1)+vv(im,jp,k,1))*0.5
          vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2))*0.5
          vyim = (vv(im,j,k,2)+vv(im,jp,k,2))*0.5
          vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3))*0.5
          vzim = (vv(im,j,k,3)+vv(im,jp,k,3))*0.5

          vxjp = vv(i,jp,k,1)
          vxjm = vv(i,j ,k,1)
          vyjp = vv(i,jp,k,2)
          vyjm = vv(i,j ,k,2)
          vzjp = vv(i,jp,k,3)
          vzjm = vv(i,j ,k,3)

          vxkp = (vv(i,j,kp,1)+vv(i,jp,kp,1))*0.5
          vxkm = (vv(i,j,km,1)+vv(i,jp,km,1))*0.5
          vykp = (vv(i,j,kp,2)+vv(i,jp,kp,2))*0.5
          vykm = (vv(i,j,km,2)+vv(i,jp,km,2))*0.5
          vzkp = (vv(i,j,kp,3)+vv(i,jp,kp,3))*0.5
          vzkm = (vv(i,j,km,3)+vv(i,jp,km,3))*0.5

          bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1))*0.5
          bxim = (bb(im,j,k,1)+bb(im,jp,k,1))*0.5
          byip = (bb(ip,j,k,2)+bb(ip,jp,k,2))*0.5
          byim = (bb(im,j,k,2)+bb(im,jp,k,2))*0.5
          bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3))*0.5
          bzim = (bb(im,j,k,3)+bb(im,jp,k,3))*0.5

          bxjp = bb(i,jp,k,1)
          bxjm = bb(i,j ,k,1)
          byjp = bb(i,jp,k,2)
          byjm = bb(i,j ,k,2)
          bzjp = bb(i,jp,k,3)
          bzjm = bb(i,j ,k,3)

          bxkp = (bb(i,j,kp,1)+bb(i,jp,kp,1))*0.5
          bxkm = (bb(i,j,km,1)+bb(i,jp,km,1))*0.5
          bykp = (bb(i,j,kp,2)+bb(i,jp,kp,2))*0.5
          bykm = (bb(i,j,km,2)+bb(i,jp,km,2))*0.5
          bzkp = (bb(i,j,kp,3)+bb(i,jp,kp,3))*0.5
          bzkm = (bb(i,j,km,3)+bb(i,jp,km,3))*0.5

        case (3)
          idhz = 1./dz(kg)
          km = k

          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
     .                 +gmetric%grid(igx)%jac(ip,j,k ))
          jacim  = 0.5*(gmetric%grid(igx)%jac(im,j,kp)
     .                 +gmetric%grid(igx)%jac(im,j,k ))
          jacjp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
     .                 +gmetric%grid(igx)%jac(i,jp,k ))
          jacjm  = 0.5*(gmetric%grid(igx)%jac(i,jm,kp)
     .                 +gmetric%grid(igx)%jac(i,jm,k ))
          jackp  = gmetric%grid(igx)%jac(i,j,kp)
          jackm  = gmetric%grid(igx)%jac(i,j,km)

          vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1))*0.5
          vxim = (vv(im,j,k,1)+vv(im,j,kp,1))*0.5
          vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2))*0.5
          vyim = (vv(im,j,k,2)+vv(im,j,kp,2))*0.5
          vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3))*0.5
          vzim = (vv(im,j,k,3)+vv(im,j,kp,3))*0.5

          vxjp = (vv(i,jp,k,1)+vv(i,jp,kp,1))*0.5
          vxjm = (vv(i,jm,k,1)+vv(i,jm,kp,1))*0.5
          vyjp = (vv(i,jp,k,2)+vv(i,jp,kp,2))*0.5
          vyjm = (vv(i,jm,k,2)+vv(i,jm,kp,2))*0.5
          vzjp = (vv(i,jp,k,3)+vv(i,jp,kp,3))*0.5
          vzjm = (vv(i,jm,k,3)+vv(i,jm,kp,3))*0.5

          vxkp = vv(i,j,kp,1)
          vxkm = vv(i,j,k ,1)
          vykp = vv(i,j,kp,2)
          vykm = vv(i,j,k ,2)
          vzkp = vv(i,j,kp,3)
          vzkm = vv(i,j,k ,3)

          bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1))*0.5
          bxim = (bb(im,j,k,1)+bb(im,j,kp,1))*0.5
          byip = (bb(ip,j,k,2)+bb(ip,j,kp,2))*0.5
          byim = (bb(im,j,k,2)+bb(im,j,kp,2))*0.5
          bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3))*0.5
          bzim = (bb(im,j,k,3)+bb(im,j,kp,3))*0.5

          bxjp = (bb(i,jp,k,1)+bb(i,jp,kp,1))*0.5
          bxjm = (bb(i,jm,k,1)+bb(i,jm,kp,1))*0.5
          byjp = (bb(i,jp,k,2)+bb(i,jp,kp,2))*0.5
          byjm = (bb(i,jm,k,2)+bb(i,jm,kp,2))*0.5
          bzjp = (bb(i,jp,k,3)+bb(i,jp,kp,3))*0.5
          bzjm = (bb(i,jm,k,3)+bb(i,jm,kp,3))*0.5

          bxkp = bb(i,j,kp,1)
          bxkm = bb(i,j,k ,1)
          bykp = bb(i,j,kp,2)
          bykm = bb(i,j,k ,2)
          bzkp = bb(i,j,kp,3)
          bzkm = bb(i,j,k ,3)

        case default
          
          jacip  = gmetric%grid(igx)%jac(ip,j,k)
          jacim  = gmetric%grid(igx)%jac(im,j,k)
          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
          jacjm  = gmetric%grid(igx)%jac(i,jm,k)
          jackp  = gmetric%grid(igx)%jac(i,j,kp)
          jackm  = gmetric%grid(igx)%jac(i,j,km)

          vxip = vv(ip,j,k,1)
          vxim = vv(im,j,k,1)
          vyip = vv(ip,j,k,2)
          vyim = vv(im,j,k,2)
          vzip = vv(ip,j,k,3)
          vzim = vv(im,j,k,3)

          vxjp = vv(i,jp,k,1)
          vxjm = vv(i,jm,k,1)
          vyjp = vv(i,jp,k,2)
          vyjm = vv(i,jm,k,2)
          vzjp = vv(i,jp,k,3)
          vzjm = vv(i,jm,k,3)

          vxkp = vv(i,j,kp,1)
          vxkm = vv(i,j,km,1)
          vykp = vv(i,j,kp,2)
          vykm = vv(i,j,km,2)
          vzkp = vv(i,j,kp,3)
          vzkm = vv(i,j,km,3)

          bxip = bb(ip,j,k,1)
          bxim = bb(im,j,k,1)
          byip = bb(ip,j,k,2)
          byim = bb(im,j,k,2)
          bzip = bb(ip,j,k,3)
          bzim = bb(im,j,k,3)

          bxjp = bb(i,jp,k,1)
          bxjm = bb(i,jm,k,1)
          byjp = bb(i,jp,k,2)
          byjm = bb(i,jm,k,2)
          bzjp = bb(i,jp,k,3)
          bzjm = bb(i,jm,k,3)

          bxkp = bb(i,j,kp,1)
          bxkm = bb(i,j,km,1)
          bykp = bb(i,j,kp,2)
          bykm = bb(i,j,km,2)
          bzkp = bb(i,j,kp,3)
          bzkm = bb(i,j,km,3)

        end select

c     Components

        ijacip  = 1d0/jacip
        ijacim  = 1d0/jacim
        ijacjp  = 1d0/jacjp
        ijacjm  = 1d0/jacjm
        ijackp  = 1d0/jackp
        ijackm  = 1d0/jackm

        !component 1

        flxjp = ( vyjp*bxjp-vxjp*byjp )*ijacjp
        flxjm = ( vyjm*bxjm-vxjm*byjm )*ijacjm

        flxkp = ( vzkp*bxkp-vxkp*bzkp )*ijackp
        flxkm = ( vzkm*bxkm-vxkm*bzkm )*ijackm

        cnv(1) =  (flxjp-flxjm)*idhy
     .           +(flxkp-flxkm)*idhz

        !component 2

        flxip = ( vxip*byip-vyip*bxip )*ijacip
        flxim = ( vxim*byim-vyim*bxim )*ijacim

        flxkp = ( vzkp*bykp-vykp*bzkp )*ijackp
        flxkm = ( vzkm*bykm-vykm*bzkm )*ijackm

        cnv(2) =  (flxip-flxim)*idhx
     .           +(flxkp-flxkm)*idhz

        !component 3

        flxip = ( vxip*bzip-vzip*bxip )*ijacip
        flxim = ( vxim*bzim-vzim*bxim )*ijacim

        flxjp = ( vyjp*bzjp-vzjp*byjp )*ijacjp
        flxjm = ( vyjm*bzjm-vzjm*byjm )*ijacjm

        cnv(3) =  (flxip-flxim)*idhx
     .           +(flxjp-flxjm)*idhy

      end function curl_bxv

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

          tsrc = 0.5*(nabla_v_src(ip,j,k)
     .               +nabla_v_src(i ,j,k))

        case (2)
          dhy = dy(jg)

          vxip = 0.5*(ax(ip,j,k)+ax(ip,jp,k))
          vxim = 0.5*(ax(im,j,k)+ax(im,jp,k))
          vyip = 0.5*(ay(ip,j,k)+ay(ip,jp,k))
          vyim = 0.5*(ay(im,j,k)+ay(im,jp,k))
          vzip = 0.5*(az(ip,j,k)+az(ip,jp,k))
          vzim = 0.5*(az(im,j,k)+az(im,jp,k))

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

          vxip = 0.5*(ax(ip,j,k)+ax(ip,j,kp))
          vxim = 0.5*(ax(im,j,k)+ax(im,j,kp))
          vyip = 0.5*(ay(ip,j,k)+ay(ip,j,kp))
          vyim = 0.5*(ay(im,j,k)+ay(im,j,kp))
          vzip = 0.5*(az(ip,j,k)+az(ip,j,kp))
          vzim = 0.5*(az(im,j,k)+az(im,j,kp))

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

          vxip = ax(ip,j,k)
          vxim = ax(im,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(im,j,k)
          vzip = az(ip,j,k)
          vzim = az(im,j,k)

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

        idhx = 1d0/dhx
        idhy = 1d0/dhy
        idhz = 1d0/dhz

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

c     Begin program

      ! l = 1, m = 1
        tensor(1,1) = ax(i,j,k)*(gmetric%grid(igrid)%Gamma(i,j,k,2,2,1)
     .                         + gmetric%grid(igrid)%Gamma(i,j,k,3,3,1))
     .               -ay(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,1,2,1)
     .               -az(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,1,3,1)

      ! l = 1, m = 2
        tensor(1,2) = ay(i,j,k)*(gmetric%grid(igrid)%Gamma(i,j,k,1,1,1)
     .                         + gmetric%grid(igrid)%Gamma(i,j,k,3,3,1))
     .               -ax(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,2,1,1)
     .               -az(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,2,3,1)

      ! l = 1, m = 3
        tensor(1,3) = az(i,j,k)*(gmetric%grid(igrid)%Gamma(i,j,k,1,1,1)
     .                         + gmetric%grid(igrid)%Gamma(i,j,k,2,2,1))
     .               -ax(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,3,1,1)
     .               -ay(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,3,2,1)

      ! l = 2, m = 1
        tensor(2,1) = ax(i,j,k)*(gmetric%grid(igrid)%Gamma(i,j,k,2,2,2)
     .                         + gmetric%grid(igrid)%Gamma(i,j,k,3,3,2))
     .               -ay(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,1,2,2)
     .               -az(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,1,3,2)

      ! l = 2, m = 2
        tensor(2,2) = ay(i,j,k)*(gmetric%grid(igrid)%Gamma(i,j,k,1,1,2)
     .                         + gmetric%grid(igrid)%Gamma(i,j,k,3,3,2))
     .               -ax(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,2,1,2)
     .               -az(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,2,3,2)

      ! l = 2, m = 3
        tensor(2,3) = az(i,j,k)*(gmetric%grid(igrid)%Gamma(i,j,k,1,1,2)
     .                         + gmetric%grid(igrid)%Gamma(i,j,k,2,2,2))
     .               -ax(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,3,1,2)
     .               -ay(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,3,2,2)

      ! l = 3, m = 1
        tensor(3,1) = ax(i,j,k)*(gmetric%grid(igrid)%Gamma(i,j,k,2,2,3)
     .                         + gmetric%grid(igrid)%Gamma(i,j,k,3,3,3))
     .               -ay(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,1,2,3)
     .               -az(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,1,3,3)

      ! l = 3, m = 2
        tensor(3,2) = ay(i,j,k)*(gmetric%grid(igrid)%Gamma(i,j,k,1,1,3)
     .                         + gmetric%grid(igrid)%Gamma(i,j,k,3,3,3))
     .               -ax(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,2,1,3)
     .               -az(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,2,3,3)

      ! l = 3, m = 3
        tensor(3,3) = az(i,j,k)*(gmetric%grid(igrid)%Gamma(i,j,k,1,1,3)
     .                         + gmetric%grid(igrid)%Gamma(i,j,k,2,2,3))
     .               -ax(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,3,1,3)
     .               -ay(i,j,k)* gmetric%grid(igrid)%Gamma(i,j,k,3,2,3)

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
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp

c     Begin program

        ip = i+1
        if (flag == 0) ip = i

        jacp = gmetric%grid(igx)%jac(ip,j,k)
        jac0 = gmetric%grid(igx)%jac(i ,j,k)
        jac  = 0.5*(jacp+jac0)

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        ijac  = 1d0/jac
        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

        t11 = 0d0

        if (solenoidal) then
          t12 = 0.5*( vec1(i ,j,k,1)*vec2(i ,j,k,2)*ijac0
     .               +vec1(ip,j,k,1)*vec2(ip,j,k,2)*ijacp)*jac
     .         -0.5*( vec1(i ,j,k,2)*vec2(i ,j,k,1)*ijac0
     .               +vec1(ip,j,k,2)*vec2(ip,j,k,1)*ijacp)*jac

          t13 = 0.5*( vec1(i ,j,k,1)*vec2(i ,j,k,3)*ijac0
     .               +vec1(ip,j,k,1)*vec2(ip,j,k,3)*ijacp)*jac
     .         -0.5*( vec1(i ,j,k,3)*vec2(i ,j,k,1)*ijac0
     .               +vec1(ip,j,k,3)*vec2(ip,j,k,1)*ijacp)*jac
        else
          t12 = 0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp
     .               +vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0)*jac
     .         -0.5*( vec1(ip,j,k,2)*vec2(i ,j,k,1)*ijac0
     .               +vec1(i ,j,k,2)*vec2(ip,j,k,1)*ijacp)*jac

          t13 = 0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,3)*ijac0
     .               +vec1(i ,j,k,1)*vec2(ip,j,k,3)*ijacp)*jac
     .         -0.5*( vec1(ip,j,k,3)*vec2(i ,j,k,1)*ijac0
     .               +vec1(i ,j,k,3)*vec2(ip,j,k,1)*ijacp)*jac
        endif

        if (flag /= 0) then
          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
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
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp

c     Begin program

        jp = j+1
        if (flag == 0) jp = j

        jacp = gmetric%grid(igx)%jac(i,jp,k)
        jac0 = gmetric%grid(igx)%jac(i,j ,k)
        jac  = 0.5*(jacp+jac0)

        ijac  = 1d0/jac
        ijacp = 1d0/jacp
        ijac0 = 1d0/jac0

        if (solenoidal) then
          t21 = 0.5*( vec1(i,j ,k,2)*vec2(i,j ,k,1)*ijac0
     .               +vec1(i,jp,k,2)*vec2(i,jp,k,1)*ijacp)*jac
     .         -0.5*( vec1(i,j ,k,1)*vec2(i,j ,k,2)*ijac0
     .               +vec1(i,jp,k,1)*vec2(i,jp,k,2)*ijacp)*jac

          t22 = 0d0

          t23 = 0.5*( vec1(i,j ,k,2)*vec2(i,j ,k,3)*ijac0
     .               +vec1(i,jp,k,2)*vec2(i,jp,k,3)*ijacp)*jac
     .         -0.5*( vec1(i,j ,k,3)*vec2(i,j ,k,2)*ijac0
     .               +vec1(i,jp,k,3)*vec2(i,jp,k,2)*ijacp)*jac
        else
          t21 = 0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,1)*ijac0
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,1)*ijacp)*jac
     .         -0.5*( vec1(i,jp,k,1)*vec2(i,j ,k,2)*ijacp
     .               +vec1(i,j ,k,1)*vec2(i,jp,k,2)*ijac0)*jac

          t22 = 0d0

          t23 = 0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,3)*ijac0
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,3)*ijacp)*jac
     .         -0.5*( vec1(i,jp,k,3)*vec2(i,j ,k,2)*ijac0
     .               +vec1(i,j ,k,3)*vec2(i,jp,k,2)*ijacp)*jac
        endif

        if (flag /= 0) then
          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
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
        real(8)    :: jac,jac0,jacp,ijac,ijacp,ijac0

c     Begin program

        kp = k+1
        if (flag == 0) kp = k

        jacp = gmetric%grid(igx)%jac(i,j,kp)
        jac0 = gmetric%grid(igx)%jac(i,j,k )
        jac  = 0.5*(jacp+jac0)

        ijac  = 1d0/jac
        ijacp = 1d0/jacp
        ijac0 = 1d0/jac0

        if (solenoidal) then
          t31 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,k ,1)*ijac0
     .               +vec1(i,j,kp,3)*vec2(i,j,kp,1)*ijacp)*jac
     .         -0.5*( vec1(i,j,k ,1)*vec2(i,j,k ,3)*ijac0
     .               +vec1(i,j,kp,1)*vec2(i,j,kp,3)*ijacp)*jac

          t32 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,k ,2)*ijac0
     .               +vec1(i,j,kp,3)*vec2(i,j,kp,2)*ijacp)*jac
     .         -0.5*( vec1(i,j,k ,2)*vec2(i,j,k ,3)*ijac0
     .               +vec1(i,j,kp,2)*vec2(i,j,kp,3)*ijacp)*jac
        else
          t31 = 0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,1)*ijac0
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,1)*ijacp)*jac
     .         -0.5*( vec1(i,j,kp,1)*vec2(i,j,k ,3)*ijac0
     .               +vec1(i,j,k ,1)*vec2(i,j,kp,3)*ijacp)*jac

          t32 = 0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,2)*ijac0
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,2)*ijacp)*jac
     .         -0.5*( vec1(i,j,kp,2)*vec2(i,j,k ,3)*ijac0
     .               +vec1(i,j,k ,2)*vec2(i,j,kp,3)*ijacp)*jac
        endif

        t33 = 0d0

        if (flag /= 0) then
          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine btensor_z

ccc     lf_x
ccc     #############################################################
cc      subroutine lf_x(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .               ,t11,t12,t13,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t11-t13 for linearized Lorentz
ccc     force term in conservative form.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
cc        real(8)    :: t11,t12,t13
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer(4) :: ip,igrid
cc        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,gsuper(3,3)
cc     .               ,scalar_prod,acnv(3),acnvp(3),bcov(3)
cc     .               ,bcovp(3),bcnv(3),bcnvp(3)
cc
ccc     Begin program
cc
cc        igrid = igx
cc
cc        ip = i+1
cc        if (flag == 0) ip = i
cc
cc        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
cc     .               +gmetric%grid(igrid)%jac (i ,j,k))
cc        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
cc     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))
cc
cc        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
cc     .      .and. bcSP()
cc     .      .and. flag /= 0           ) then
cc          jacp = gmetric%grid(igrid)%jac(ip,j,k)
cc          jac0 = gmetric%grid(igrid)%jac(i ,j,k)
cc        else
cc          jacp = jac
cc          jac0 = jac
cc        endif
cc
cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cc        ijac  = 1d0/jac
cc        ijacp = 1d0/jacp
cc        ijac0 = 1d0/jac0
cc
cccc        acnv(1) = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cccc        acnv(2) = 0.5*(vec1(ip,j,k,2)     +vec1(i,j,k,2))
cccc        acnv(3) = 0.5*(vec1(ip,j,k,3)*ijacp+vec1(i,j,k,3)*ijac0)*jac
cccc
cccc        bcnv(1) = 0.5*(vec2(ip,j,k,1)*ijacp+vec2(i,j,k,1)*ijac0)*jac
cccc        bcnv(2) = 0.5*(vec2(ip,j,k,2)     +vec2(i,j,k,2))
cccc        bcnv(3) = 0.5*(vec2(ip,j,k,3)*ijacp+vec2(i,j,k,3)*ijac0)*jac
cccc
cccc        if (flag /= 0) then
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=1)
cccc        else
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=0)
cccc        endif
cccc
cccc        scalar_prod = dot_product(acnv,bcov)
cccc
cccc        t11 =( acnv(1)*bcnv(1)
cccc     .        +acnv(1)*bcnv(1)
cccc     .        -gsuper(1,1)*scalar_prod )
cccc
cccc        t12 =( acnv(1)*bcnv(2)
cccc     .        +acnv(2)*bcnv(1)
cccc     .        -gsuper(1,2)*scalar_prod )
cccc
cccc        t13 =( acnv(1)*bcnv(3)
cccc     .        +acnv(3)*bcnv(1)
cccc     .        -gsuper(1,3)*scalar_prod )
cc
cc        acnv(1) = vec1(i,j,k,1)*ijac0*jac
cc        acnv(2) = vec1(i,j,k,2)
cc        acnv(3) = vec1(i,j,k,3)*ijac0*jac
cc
cc        acnvp(1) = vec1(ip,j,k,1)*ijacp*jac
cc        acnvp(2) = vec1(ip,j,k,2)
cc        acnvp(3) = vec1(ip,j,k,3)*ijacp*jac
cc
cccc        bcnv(1) = vec2(i,j,k,1)*ijac0*jac
cccc        bcnv(2) = vec2(i,j,k,2)
cccc        bcnv(3) = vec2(i,j,k,3)*ijac0*jac
cccc
cccc        bcnvp(1) = vec2(ip,j,k,1)*ijacp*jac
cccc        bcnvp(2) = vec2(ip,j,k,2)
cccc        bcnvp(3) = vec2(ip,j,k,3)*ijacp*jac
cc
cc        bcnv(1) = vec2(i,j,k,1)
cc        bcnv(2) = vec2(i,j,k,2)
cc        bcnv(3) = vec2(i,j,k,3)
cc
cc        bcnvp(1) = vec2(ip,j,k,1)
cc        bcnvp(2) = vec2(ip,j,k,2)
cc        bcnvp(3) = vec2(ip,j,k,3)
cc
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        call transformFromCurvToCurv(ip,j,k,igx,igy,igz
cc     .                        ,bcovp(1),bcovp(2),bcovp(3)
cc     .                        ,bcnvp(1),bcnvp(2),bcnvp(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        bcnv(1) = bcnv(1)*ijac0*jac
cc        bcnv(3) = bcnv(3)*ijac0*jac
cc
cc        bcnvp(1) = bcnvp(1)*ijacp*jac
cc        bcnvp(3) = bcnvp(3)*ijacp*jac
cc
cc        bcov (2) = bcov (2)*ijac0*jac
cc        bcovp(2) = bcovp(2)*ijacp*jac
cc
cc        scalar_prod = dot_product(acnvp,bcov)
cc     .               +dot_product(acnv,bcovp)
cc
cc        t11 =0.5*( 2*acnvp(1)*bcnv(1)
cc     .            +2*acnv(1)*bcnvp(1)
cc     .            -gsuper(1,1)*scalar_prod )
cc
cc        t12 =0.5*( acnvp(1)*bcnv(2) + acnv(1)*bcnvp(2)
cc     .            +acnvp(2)*bcnv(1) + acnv(2)*bcnvp(1)
cc     .            -gsuper(1,2)*scalar_prod )
cc
cc        t13 =0.5*( acnvp(1)*bcnv(3) + acnv(1)*bcnvp(3)
cc     .            +acnvp(3)*bcnv(1) + acnv(3)*bcnvp(1)
cc     .            -gsuper(1,3)*scalar_prod )
cc
cccc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cccc        ijac  = 1d0/jac
cc
cc        if (flag /= 0) then
cc          t11 = t11*ijac
cc          if (.not.alt_eom) t12 = t12*ijac
cc          t13 = t13*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine lf_x
cc
ccc     lf_y
ccc     #############################################################
cc      subroutine lf_y(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .               ,t21,t22,t23,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t21-t23 for linearized Lorentz
ccc     force term in conservative form.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
cc        real(8)    :: t21,t22,t23
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer(4) :: jp,igrid
cc        real(8)    :: jac,ijac,gsuper(3,3),scalar_prod
cc     .               ,acnv(3),acnvp(3),bcov(3),bcovp(3),bcnv(3),bcnvp(3)
cc
ccc     Begin program
cc
cc        igrid = igx
cc
cc        jp = j+1
cc        if (flag == 0) jp = j
cc
cc        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
cc     .               +gmetric%grid(igrid)%jac (i,j ,k))
cc        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
cc     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))
cc
cc        ijac = 1d0/jac
cc
cccc        acnv(1) = 0.5*(vec1(i,jp,k,1)+vec1(i,j,k,1))
cccc        acnv(2) = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cccc        acnv(3) = 0.5*(vec1(i,jp,k,3)+vec1(i,j,k,3))
cccc                                
cccc        bcnv(1) = 0.5*(vec2(i,jp,k,1)+vec2(i,j,k,1))
cccc        bcnv(2) = 0.5*(vec2(i,jp,k,2)+vec2(i,j,k,2))
cccc        bcnv(3) = 0.5*(vec2(i,jp,k,3)+vec2(i,j,k,3))
cccc
cccc        if (flag /= 0) then
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=2)
cccc        else
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=0)
cccc        endif
cccc
cccc
cccc        scalar_prod = dot_product(acnv,bcov)
cccc
cccc        t21 =( acnv(2)*bcnv(1)
cccc     .        +acnv(1)*bcnv(2)
cccc     .        -gsuper(2,1)*scalar_prod )
cccc
cccc        t22 =( acnv(2)*bcnv(2)
cccc     .        +acnv(2)*bcnv(2)
cccc     .        -gsuper(2,2)*scalar_prod )
cccc
cccc        t23 =( acnv(2)*bcnv(3)
cccc     .        +acnv(3)*bcnv(2)
cccc     .        -gsuper(2,3)*scalar_prod )
cc
cc        acnv(1) = vec1(i,j,k,1)
cc        acnv(2) = vec1(i,j,k,2)
cc        acnv(3) = vec1(i,j,k,3)
cc
cc        acnvp(1) = vec1(i,jp,k,1)
cc        acnvp(2) = vec1(i,jp,k,2)
cc        acnvp(3) = vec1(i,jp,k,3)
cc
cc        bcnv(1) = vec2(i,j,k,1)
cc        bcnv(2) = vec2(i,j,k,2)
cc        bcnv(3) = vec2(i,j,k,3)
cc
cc        bcnvp(1) = vec2(i,jp,k,1)
cc        bcnvp(2) = vec2(i,jp,k,2)
cc        bcnvp(3) = vec2(i,jp,k,3)
cc
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        call transformFromCurvToCurv(i,jp,k,igx,igy,igz
cc     .                        ,bcovp(1),bcovp(2),bcovp(3)
cc     .                        ,bcnvp(1),bcnvp(2),bcnvp(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        scalar_prod = dot_product(acnvp,bcov)
cc     .               +dot_product(acnv,bcovp)
cc
cc        t21 =0.5*( acnvp(2)*bcnv(1) + acnv(2)*bcnvp(1)
cc     .            +acnvp(1)*bcnv(2) + acnv(1)*bcnvp(2)
cc     .            -gsuper(2,1)*scalar_prod )
cc
cc        t22 =0.5*( 2*acnvp(2)*bcnv(2) + 2*acnv(2)*bcnvp(2)
cc     .            -gsuper(2,2)*scalar_prod )
cc
cc        t23 =0.5*( acnvp(2)*bcnv(3) + acnv(2)*bcnvp(3)
cc     .            +acnvp(3)*bcnv(2) + acnv(3)*bcnvp(2)
cc     .            -gsuper(2,3)*scalar_prod )
cc
cc        if (flag /= 0) then
cc          t21 = t21*ijac
cc          if (.not.alt_eom) t22 = t22*ijac
cc          t23 = t23*ijac
cc        endif
cc
cc
ccc     End program
cc
cc      end subroutine lf_y
cc
ccc     lf_z
ccc     #############################################################
cc      subroutine lf_z(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .                   ,t31,t32,t33,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t31-t33 for linearized Lorentz
ccc     force term in conservative form.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
cc        real(8)    :: t31,t32,t33
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer(4) :: kp,igrid
cc        real(8)    :: jac,ijac,gsuper(3,3),scalar_prod
cc     .               ,acnv(3),acnvp(3),bcov(3),bcovp(3),bcnv(3),bcnvp(3)
cc
ccc     Begin program
cc
cc        igrid = igx
cc
cc        kp=k+1
cc        if (flag == 0) kp = k
cc
cc        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
cc     .               +gmetric%grid(igrid)%jac (i,j,k ))
cc        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
cc     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))
cc
cc        ijac = 1d0/jac
cc
cccc        acnv(1) = 0.5*(vec1(i,j,kp,1)+vec1(i,j,k,1))
cccc        acnv(2) = 0.5*(vec1(i,j,kp,2)+vec1(i,j,k,2))
cccc        acnv(3) = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cccc                                  
cccc        bcnv(1) = 0.5*(vec2(i,j,kp,1)+vec2(i,j,k,1))
cccc        bcnv(2) = 0.5*(vec2(i,j,kp,2)+vec2(i,j,k,2))
cccc        bcnv(3) = 0.5*(vec2(i,j,kp,3)+vec2(i,j,k,3))
cccc
cccc        if (flag /= 0) then
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=3)
cccc        else
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=0)
cccc        endif
cccc
cccc        scalar_prod = dot_product(acnv,bcov)
cccc
cccc        t31 =( acnv(3)*bcnv(1)
cccc     .        +acnv(1)*bcnv(3)
cccc     .        -gsuper(3,1)*scalar_prod )
cccc
cccc        t32 =( acnv(3)*bcnv(2)
cccc     .        +acnv(2)*bcnv(3)
cccc     .        -gsuper(3,2)*scalar_prod )
cccc
cccc        t33 =( acnv(3)*bcnv(3)
cccc     .        +acnv(3)*bcnv(3)
cccc     .        -gsuper(3,3)*scalar_prod )
cc
cc        acnv(1) = vec1(i,j,k,1)
cc        acnv(2) = vec1(i,j,k,2)
cc        acnv(3) = vec1(i,j,k,3)
cc
cc        acnvp(1) = vec1(i,j,kp,1)
cc        acnvp(2) = vec1(i,j,kp,2)
cc        acnvp(3) = vec1(i,j,kp,3)
cc
cc        bcnv(1) = vec2(i,j,k,1)
cc        bcnv(2) = vec2(i,j,k,2)
cc        bcnv(3) = vec2(i,j,k,3)
cc
cc        bcnvp(1) = vec2(i,j,kp,1)
cc        bcnvp(2) = vec2(i,j,kp,2)
cc        bcnvp(3) = vec2(i,j,kp,3)
cc
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        call transformFromCurvToCurv(i,j,kp,igx,igy,igz
cc     .                        ,bcovp(1),bcovp(2),bcovp(3)
cc     .                        ,bcnvp(1),bcnvp(2),bcnvp(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        scalar_prod = dot_product(acnvp,bcov)
cc     .               +dot_product(acnv,bcovp)
cc
cc        t31 =0.5*( acnvp(3)*bcnv(1) + acnv(3)*bcnvp(1)
cc     .            +acnvp(1)*bcnv(3) + acnv(1)*bcnvp(3)
cc     .            -gsuper(3,1)*scalar_prod )
cc
cc        t32 =0.5*( acnvp(3)*bcnv(2) + acnv(3)*bcnvp(2)
cc     .            +acnvp(2)*bcnv(3) + acnv(2)*bcnvp(3)
cc     .            -gsuper(3,2)*scalar_prod )
cc
cc        t33 =0.5*( 2*acnvp(3)*bcnv(3)
cc     .           + 2*acnv(3)*bcnvp(3)
cc     .            -gsuper(3,3)*scalar_prod )
cc
cc        if (flag /= 0) then
cc          t31 = t31*ijac
cc          if (.not.alt_eom) t32 = t32*ijac
cc          t33 = t33*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine lf_z

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

        use timeStepping

        real(8),pointer,dimension(:,:,:):: rho,rvx,rvy,rvz,bx,by,bz,tmp

        real(8) :: max_dv_dt

      contains

ccc     vtensor_x
ccc     #############################################################
cc      subroutine vtensor_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t11,t12,t13,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t11-t13 for EOM. In the call
ccc     sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t11,t12,t13: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t11,t12,t13
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg,ip
cc        real(8)    :: x,y,z
cc        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,ptot,vis
cc
ccc     Begin program
cc
cc        ip = i+1
cc        if (flag == 0) ip = i
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
cc     .               +gmetric%grid(igx)%jac (i ,j,k))
cc        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
cc     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))
cc
cc        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
cc     .      .and. bcSP()
cc     .      .and. flag /= 0           ) then
cc          jacp = gmetric%grid(igx)%jac(ip,j,k)
cc          jac0 = gmetric%grid(igx)%jac(i ,j,k)
cc        else
cc          jacp = jac
cc          jac0 = jac
cc        endif
cc
cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cc
cc        ijac  = 1d0/jac
cc        ijac0 = 1d0/jac0
cc        ijacp = 1d0/jacp
cc
cc        if (flag /= 0) then
cc          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,1)
cc        else
cc          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)
cc        endif
cc
cc        !Harmonic average for calculation of viscosity coeff. at faces
cc        vis = 2./(1./nuu(ip,j,k) + 1./nuu(i ,j,k))
cc
cc        if (nc_eom_f.and.(.not.nc_eom_v)) then
cc         !Recall p=2nT
cccc          ptot = jac*(rho(ip,j,k)*tmp(i ,j,k)
cccc     .               +rho(i ,j,k)*tmp(ip,j,k))
cccc          ptot = jac*(rho(ip,j,k)*tmp(ip,j,k)
cccc     .               +rho(i ,j,k)*tmp(i ,j,k))
cc          ptot = 0d0
cc
cc          t11 =
cc     .      0.5*( rvx(ip,j,k)*vx (i ,j,k)*ijacp*ijac0
cc     .           +rvx(i ,j,k)*vx (ip,j,k)*ijacp*ijac0)*jac**2
cc     .     +gsuper(1,1)*ptot
cc     .     -vis*( gsuper(1,1)*nabla_v(1,1)
cc     .           +gsuper(1,2)*nabla_v(2,1)
cc     .           +gsuper(1,3)*nabla_v(3,1) )
cc
cc          t12 =
cc     .     0.25*( rvx(ip,j,k)*vy (i ,j,k)*ijacp
cc     .           +rvx(i ,j,k)*vy (ip,j,k)*ijac0
cc     .           + vx(i ,j,k)*rvy(ip,j,k)*ijac0
cc     .           + vx(ip,j,k)*rvy(i ,j,k)*ijacp)*jac
cc     .     +gsuper(1,2)*ptot
cc     .     -vis*( gsuper(1,1)*nabla_v(1,2)
cc     .           +gsuper(1,2)*nabla_v(2,2)
cc     .           +gsuper(1,3)*nabla_v(3,2) )
cc
cc          t13 =
cc     .      0.25*( rvx(ip,j,k)*vz (i ,j,k)*ijacp*ijac0
cc     .            +rvx(i ,j,k)*vz (ip,j,k)*ijacp*ijac0
cc     .            + vx(i ,j,k)*rvz(ip,j,k)*ijacp*ijac0
cc     .            + vx(ip,j,k)*rvz(i ,j,k)*ijacp*ijac0)*jac**2
cc     .     +gsuper(1,3)*ptot
cc     .     -vis*( gsuper(1,1)*nabla_v(1,3)
cc     .           +gsuper(1,2)*nabla_v(2,3)
cc     .           +gsuper(1,3)*nabla_v(3,3) )
cc        elseif (nc_eom_v.and.(.not.nc_eom_f)) then
cc         !Recall p=2nT
cc          ptot = jac*(bx(ip,j,k)*bx_cov(i ,j,k)*ijacp
cc     .               +bx(i ,j,k)*bx_cov(ip,j,k)*ijac0
cc     .               +by(ip,j,k)*by_cov(i ,j,k)*ijac0
cc     .               +by(i ,j,k)*by_cov(ip,j,k)*ijacp
cc     .               +bz(ip,j,k)*bz_cov(i ,j,k)*ijacp
cc     .               +bz(i ,j,k)*bz_cov(ip,j,k)*ijac0)/4.
cc
cc          t11 =
cc     .     -    ( bx(ip,j,k)*bx(i ,j,k)*ijacp*ijac0 )*jac**2
cc     .     +gsuper(1,1)*ptot
cc
cc          t12 =
cc     .     -0.5*( bx(ip,j,k)*by(i ,j,k)*ijacp
cc     .           +bx(i ,j,k)*by(ip,j,k)*ijac0 )*jac
cc     .     +gsuper(1,2)*ptot
cc
cc          t13 =
cc     .     -0.5*( bx(ip,j,k)*bz(i ,j,k)*ijacp*ijac0
cc     .           +bx(i ,j,k)*bz(ip,j,k)*ijacp*ijac0)*jac**2
cc     .     +gsuper(1,3)*ptot
cc        elseif (.not.(nc_eom_v.or.nc_eom_f)) then
cc         !Recall p=2nT
cc          ptot = jac*(rho(ip,j,k)*tmp(i ,j,k)
cc     .               +rho(i ,j,k)*tmp(ip,j,k))
cccc          ptot = jac*(rho(ip,j,k)*tmp(ip,j,k)
cccc     .               +rho(i ,j,k)*tmp(i ,j,k))
cc     .       +jac*(bx(ip,j,k)*bx_cov(i ,j,k)*ijacp
cc     .            +bx(i ,j,k)*bx_cov(ip,j,k)*ijac0
cc     .            +by(ip,j,k)*by_cov(i ,j,k)*ijac0
cc     .            +by(i ,j,k)*by_cov(ip,j,k)*ijacp
cc     .            +bz(ip,j,k)*bz_cov(i ,j,k)*ijacp
cc     .            +bz(i ,j,k)*bz_cov(ip,j,k)*ijac0)/4.
cc
cc          t11 =
cc     .      0.5*( rvx(ip,j,k)*vx (i ,j,k)*ijacp*ijac0
cc     .           +rvx(i ,j,k)*vx (ip,j,k)*ijacp*ijac0)*jac**2
cc     .     -    ( bx(ip,j,k)*bx(i ,j,k)*ijacp*ijac0 )*jac**2
cc     .     +gsuper(1,1)*ptot
cc     .     -vis*( gsuper(1,1)*nabla_v(1,1)
cc     .           +gsuper(1,2)*nabla_v(2,1)
cc     .           +gsuper(1,3)*nabla_v(3,1) )
cc
cc          t12 =
cc     .     0.25*( rvx(ip,j,k)*vy (i ,j,k)*ijacp
cc     .           +rvx(i ,j,k)*vy (ip,j,k)*ijac0
cc     .           + vx(i ,j,k)*rvy(ip,j,k)*ijac0
cc     .           + vx(ip,j,k)*rvy(i ,j,k)*ijacp)*jac
cc     .     -0.5*( bx(ip,j,k)*by(i ,j,k)*ijacp
cc     .           +bx(i ,j,k)*by(ip,j,k)*ijac0 )*jac
cc     .     +gsuper(1,2)*ptot
cc     .     -vis*( gsuper(1,1)*nabla_v(1,2)
cc     .           +gsuper(1,2)*nabla_v(2,2)
cc     .           +gsuper(1,3)*nabla_v(3,2) )
cc
cc          t13 =
cc     .      0.25*( rvx(ip,j,k)*vz (i ,j,k)*ijacp*ijac0
cc     .            +rvx(i ,j,k)*vz (ip,j,k)*ijacp*ijac0
cc     .            + vx(i ,j,k)*rvz(ip,j,k)*ijacp*ijac0
cc     .            + vx(ip,j,k)*rvz(i ,j,k)*ijacp*ijac0)*jac**2
cc     .     -0.5*( bx(ip,j,k)*bz(i ,j,k)*ijacp*ijac0
cc     .           +bx(i ,j,k)*bz(ip,j,k)*ijacp*ijac0)*jac**2
cc     .     +gsuper(1,3)*ptot
cc     .     -vis*( gsuper(1,1)*nabla_v(1,3)
cc     .           +gsuper(1,2)*nabla_v(2,3)
cc     .           +gsuper(1,3)*nabla_v(3,3) )
cc        else
cc          t11 = 0d0
cc          t12 = 0d0
cc          t13 = 0d0
cc        endif
cc
cccc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cccc        ijac = 1d0/jac
cc
cc        if (flag /= 0) then
cc          t11 = t11*ijac
cc          if (.not.alt_eom) t12 = t12*ijac
cc          t13 = t13*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine vtensor_x
cc
ccc     vtensor_y
ccc     #############################################################
cc      subroutine vtensor_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t21,t22,t23,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t21-t23 for EOM. In the call
ccc     sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t21,t22,t23: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t21,t22,t23
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg,jp
cc        real(8)    :: x,y,z
cc        real(8)    :: jac,ijac,ptot,vis
cc
ccc     Begin program
cc
cc        jp = j+1
cc        if (flag == 0) jp = j
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
cc     .               +gmetric%grid(igx)%jac (i,j ,k))
cc        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
cc     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))
cc
cc        if (flag /= 0) then
cc          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,2)
cc        else
cc          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)
cc        endif
cc
cc        vis = 2./(1./nuu(i,jp,k) + 1./nuu(i,j ,k))
cc
cc        if (nc_eom_f.and.(.not.nc_eom_v)) then
cccc          ptot = jac*(rho(i,jp,k)*tmp(i,j ,k)
cccc     .               +rho(i,j ,k)*tmp(i,jp,k))
cc          ptot = 0d0
cc
cc          t21 =
cc     .       0.25*( rvy(i,jp,k)*vx(i,j,k) + rvy(i,j,k)*vx(i,jp,k)
cc     .             +rvx(i,jp,k)*vy(i,j,k) + rvx(i,j,k)*vy(i,jp,k))
cc     .       +gsuper(2,1)*ptot
cc     .       -vis*( gsuper(2,1)*nabla_v(1,1)
cc     .             +gsuper(2,2)*nabla_v(2,1)
cc     .             +gsuper(2,3)*nabla_v(3,1) )
cc
cc          t22 =
cc     .       0.25*( rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k)
cc     .             +rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k))
cc     .       +gsuper(2,2)*ptot
cc     .       -vis*( gsuper(2,1)*nabla_v(1,2)
cc     .             +gsuper(2,2)*nabla_v(2,2)
cc     .             +gsuper(2,3)*nabla_v(3,2) )
cc
cc          t23 =
cc     .       0.25*( rvy(i,jp,k)*vz(i,j,k) + rvy(i,j,k)*vz(i,jp,k)
cc     .             +rvz(i,jp,k)*vy(i,j,k) + rvz(i,j,k)*vy(i,jp,k))
cc     .       +gsuper(2,3)*ptot
cc     .       -vis*( gsuper(2,1)*nabla_v(1,3)
cc     .             +gsuper(2,2)*nabla_v(2,3)
cc     .             +gsuper(2,3)*nabla_v(3,3) )
cc        elseif (nc_eom_v.and.(.not.nc_eom_f)) then
cc          !Recall p=2nT
cccc          ptot = jac*(rho(i,jp,k)*tmp(i,jp,k)
cccc     .               +rho(i,j ,k)*tmp(i,j ,k))
cc          ptot =(bx(i,jp,k)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,jp,k)
cc     .          +by(i,jp,k)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,jp,k)
cc     .          +bz(i,jp,k)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,jp,k))*0.25
cc
cc          t21 =
cc     .       -0.5*( by(i,jp,k)*bx(i,j ,k)
cc     .             +by(i,j ,k)*bx(i,jp,k) )
cc     .       +gsuper(2,1)*ptot
cc
cc          t22 =
cc     .       -0.5*( by(i,jp,k)*by(i,j ,k)
cc     .             +by(i,j ,k)*by(i,jp,k) )
cc     .       +gsuper(2,2)*ptot
cc
cc          t23 =
cc     .       -0.5*( by(i,jp,k)*bz(i,j ,k)
cc     .             +by(i,j ,k)*bz(i,jp,k) )
cc     .       +gsuper(2,3)*ptot
cc        elseif (.not.(nc_eom_v.or.nc_eom_f)) then
cc          !Recall p=2nT
cccc          ptot = jac*(rho(i,jp,k)*tmp(i,jp,k)
cccc     .               +rho(i,j ,k)*tmp(i,j ,k))
cc          ptot = jac*(rho(i,jp,k)*tmp(i,j ,k)
cc     .               +rho(i,j ,k)*tmp(i,jp,k))
cc     .         +(bx(i,jp,k)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,jp,k)
cc     .          +by(i,jp,k)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,jp,k)
cc     .          +bz(i,jp,k)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,jp,k))*0.25
cc
cc          t21 =
cc     .       0.25*( rvy(i,jp,k)*vx(i,j,k) + rvy(i,j,k)*vx(i,jp,k)
cc     .             +rvx(i,jp,k)*vy(i,j,k) + rvx(i,j,k)*vy(i,jp,k))
cc     .       -0.5*( by(i,jp,k)*bx(i,j ,k)
cc     .             +by(i,j ,k)*bx(i,jp,k) )
cc     .       +gsuper(2,1)*ptot
cc     .       -vis*( gsuper(2,1)*nabla_v(1,1)
cc     .             +gsuper(2,2)*nabla_v(2,1)
cc     .             +gsuper(2,3)*nabla_v(3,1) )
cc
cc          t22 =
cc     .       0.25*( rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k)
cc     .             +rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k))
cc     .       -0.5*( by(i,jp,k)*by(i,j ,k)
cc     .             +by(i,j ,k)*by(i,jp,k) )
cc     .       +gsuper(2,2)*ptot
cc     .       -vis*( gsuper(2,1)*nabla_v(1,2)
cc     .             +gsuper(2,2)*nabla_v(2,2)
cc     .             +gsuper(2,3)*nabla_v(3,2) )
cc
cc          t23 =
cc     .       0.25*( rvy(i,jp,k)*vz(i,j,k) + rvy(i,j,k)*vz(i,jp,k)
cc     .             +rvz(i,jp,k)*vy(i,j,k) + rvz(i,j,k)*vy(i,jp,k))
cc     .       -0.5*( by(i,jp,k)*bz(i,j ,k)
cc     .             +by(i,j ,k)*bz(i,jp,k) )
cc     .       +gsuper(2,3)*ptot
cc     .       -vis*( gsuper(2,1)*nabla_v(1,3)
cc     .             +gsuper(2,2)*nabla_v(2,3)
cc     .             +gsuper(2,3)*nabla_v(3,3) )
cc        else
cc          t21 = 0d0
cc          t22 = 0d0
cc          t23 = 0d0
cc        endif
cc
cc        ijac = 1d0/jac
cc
cc        if (flag /= 0) then
cc          t21 = t21*ijac
cc          if (.not.alt_eom) t22 = t22*ijac
cc          t23 = t23*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine vtensor_y
cc
ccc     vtensor_z
ccc     #############################################################
cc      subroutine vtensor_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t31,t32,t33,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t31-t33 for EOM. In the call
ccc     sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t31,t32,t33: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer(4) :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t31,t32,t33
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer(4) :: ig,jg,kg,kp
cc        real(8)    :: x,y,z
cc        real(8)    :: jac,ijac,ptot,vis
cc
ccc     Begin program
cc
cc        kp = k+1
cc        if (flag == 0) kp = k
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
cc     .               +gmetric%grid(igx)%jac (i,j,k ))
cc        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
cc     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))
cc
cc        if (flag /= 0) then
cc          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,3)
cc        else
cc          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)
cc        endif
cc
cc        vis = 2./(1./nuu(i,j,kp) + 1./nuu(i,j,k ))
cc
cc        if (nc_eom_f.and.(.not.nc_eom_v)) then
cc          !Recall p=2nT
cccc          ptot = jac*(rho(i,j,kp)*tmp(i,j,kp)
cccc     .               +rho(i,j,k )*tmp(i,j,k ))
cccc          ptot = jac*(rho(i,j,kp)*tmp(i,j,k )
cccc     .               +rho(i,j,k )*tmp(i,j,kp))
cc          ptot = 0d0
cc
cc          t31 =
cc     .       0.25*( rvz(i,j,kp)*vx(i,j,k) + rvz(i,j,k)*vx(i,j,kp)
cc     .             +rvx(i,j,kp)*vz(i,j,k) + rvx(i,j,k)*vz(i,j,kp) )
cc     .       +gsuper(3,1)*ptot
cc     .       -vis*( gsuper(3,1)*nabla_v(1,1)
cc     .             +gsuper(3,2)*nabla_v(2,1)
cc     .             +gsuper(3,3)*nabla_v(3,1) )
cc
cc          t32 =
cc     .       0.25*( rvz(i,j,kp)*vy(i,j,k) + rvz(i,j,k)*vy(i,j,kp)
cc     .             +rvy(i,j,kp)*vz(i,j,k) + rvy(i,j,k)*vz(i,j,kp) )
cc     .       +gsuper(3,2)*ptot
cc     .       -vis*( gsuper(3,1)*nabla_v(1,2)
cc     .             +gsuper(3,2)*nabla_v(2,2)
cc     .             +gsuper(3,3)*nabla_v(3,2) )
cc
cc          t33 =
cc     .       0.25*( rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp)
cc     .             +rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp) )
cc     .       +gsuper(3,3)*ptot
cc     .       -vis*( gsuper(3,1)*nabla_v(1,3)
cc     .             +gsuper(3,2)*nabla_v(2,3)
cc     .             +gsuper(3,3)*nabla_v(3,3) )
cc        elseif (nc_eom_v.and.(.not.nc_eom_f)) then
cc          !Recall p=2nT
cccc          ptot = jac*(rho(i,j,kp)*tmp(i,j,kp)
cccc     .               +rho(i,j,k )*tmp(i,j,k ))
cc          ptot =(bx(i,j,kp)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,j,kp)
cc     .          +by(i,j,kp)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,j,kp)
cc     .          +bz(i,j,kp)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,j,kp))*0.25
cc
cc          t31 =
cc     .       -0.5*( bz(i,j,kp)*bx(i,j,k )
cc     .             +bz(i,j,k )*bx(i,j,kp) )
cc     .       +gsuper(3,1)*ptot
cc
cc          t32 =
cc     .       -0.5*( bz(i,j,kp)*by(i,j,k )
cc     .             +bz(i,j,k )*by(i,j,kp) )
cc     .       +gsuper(3,2)*ptot
cc
cc          t33 =
cc     .       -0.5*( bz(i,j,kp)*bz(i,j,k )
cc     .             +bz(i,j,k )*bz(i,j,kp) )
cc     .       +gsuper(3,3)*ptot
cc        elseif (.not.(nc_eom_v.or.nc_eom_f)) then
cc          !Recall p=2nT
cccc          ptot = jac*(rho(i,j,kp)*tmp(i,j,kp)
cccc     .               +rho(i,j,k )*tmp(i,j,k ))
cc          ptot = jac*(rho(i,j,kp)*tmp(i,j,k )
cc     .               +rho(i,j,k )*tmp(i,j,kp))
cc     .         +(bx(i,j,kp)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,j,kp)
cc     .          +by(i,j,kp)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,j,kp)
cc     .          +bz(i,j,kp)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,j,kp))*0.25
cc
cc          t31 =
cc     .       0.25*( rvz(i,j,kp)*vx(i,j,k) + rvz(i,j,k)*vx(i,j,kp)
cc     .             +rvx(i,j,kp)*vz(i,j,k) + rvx(i,j,k)*vz(i,j,kp) )
cc     .       -0.5*( bz(i,j,kp)*bx(i,j,k )
cc     .             +bz(i,j,k )*bx(i,j,kp) )
cc     .       +gsuper(3,1)*ptot
cc     .       -vis*( gsuper(3,1)*nabla_v(1,1)
cc     .             +gsuper(3,2)*nabla_v(2,1)
cc     .             +gsuper(3,3)*nabla_v(3,1) )
cc
cc          t32 =
cc     .       0.25*( rvz(i,j,kp)*vy(i,j,k) + rvz(i,j,k)*vy(i,j,kp)
cc     .             +rvy(i,j,kp)*vz(i,j,k) + rvy(i,j,k)*vz(i,j,kp) )
cc     .       -0.5*( bz(i,j,kp)*by(i,j,k )
cc     .             +bz(i,j,k )*by(i,j,kp) )
cc     .       +gsuper(3,2)*ptot
cc     .       -vis*( gsuper(3,1)*nabla_v(1,2)
cc     .             +gsuper(3,2)*nabla_v(2,2)
cc     .             +gsuper(3,3)*nabla_v(3,2) )
cc
cc          t33 =
cc     .       0.25*( rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp)
cc     .             +rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp) )
cc     .       -0.5*( bz(i,j,kp)*bz(i,j,k )
cc     .             +bz(i,j,k )*bz(i,j,kp) )
cc     .       +gsuper(3,3)*ptot
cc     .       -vis*( gsuper(3,1)*nabla_v(1,3)
cc     .             +gsuper(3,2)*nabla_v(2,3)
cc     .             +gsuper(3,3)*nabla_v(3,3) )
cc        else
cc          t31 = 0d0
cc          t32 = 0d0
cc          t33 = 0d0
cc        endif
cc
cc        ijac = 1d0/jac
cc
cc        if (flag /= 0) then
cc          t31 = t31*ijac
cc          if (.not.alt_eom) t32 = t32*ijac
cc          t33 = t33*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine vtensor_z

c     eom_advc_x
c     #############################################################
      subroutine eom_advc_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,vis

c     Begin program

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
     .               +gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac0 = gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

        if (flag /= 0) then
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,1)
        else
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)
        endif

        !Harmonic average for calculation of viscosity coeff. at faces
        vis = 2./(1./nuu(ip,j,k) + 1./nuu(i ,j,k))

        t11 =
     .      0.5*( rvx(ip,j,k)*vx (i ,j,k)*ijacp*ijac0
     .           +rvx(i ,j,k)*vx (ip,j,k)*ijacp*ijac0)*jac**2
     .     -vis*( gsuper(1,1)*nabla_v(1,1)
     .           +gsuper(1,2)*nabla_v(2,1)
     .           +gsuper(1,3)*nabla_v(3,1) )

        t12 =
     .     0.25*( rvx(ip,j,k)*vy (i ,j,k)*ijacp
     .           +rvx(i ,j,k)*vy (ip,j,k)*ijac0
     .           + vx(i ,j,k)*rvy(ip,j,k)*ijac0
     .           + vx(ip,j,k)*rvy(i ,j,k)*ijacp)*jac
     .     -vis*( gsuper(1,1)*nabla_v(1,2)
     .           +gsuper(1,2)*nabla_v(2,2)
     .           +gsuper(1,3)*nabla_v(3,2) )

        t13 =
     .      0.25*( rvx(ip,j,k)*vz (i ,j,k)*ijacp*ijac0
     .            +rvx(i ,j,k)*vz (ip,j,k)*ijacp*ijac0
     .            + vx(i ,j,k)*rvz(ip,j,k)*ijacp*ijac0
     .            + vx(ip,j,k)*rvz(i ,j,k)*ijacp*ijac0)*jac**2
     .     -vis*( gsuper(1,1)*nabla_v(1,3)
     .           +gsuper(1,2)*nabla_v(2,3)
     .           +gsuper(1,3)*nabla_v(3,3) )

        if (flag /= 0) then
          ijac  = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine eom_advc_x

c     eom_advc_y
c     #############################################################
      subroutine eom_advc_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
        real(8)    :: jac,ijac,vis

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

        vis = 2./(1./nuu(i,jp,k) + 1./nuu(i,j ,k))

        t21 =
     .       0.25*( rvy(i,jp,k)*vx(i,j,k) + rvy(i,j,k)*vx(i,jp,k)
     .             +rvx(i,jp,k)*vy(i,j,k) + rvx(i,j,k)*vy(i,jp,k))
     .       -vis*( gsuper(2,1)*nabla_v(1,1)
     .             +gsuper(2,2)*nabla_v(2,1)
     .             +gsuper(2,3)*nabla_v(3,1) )

        t22 =
     .       0.25*( rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k)
     .             +rvy(i,jp,k)*vy(i,j,k) + rvy(i,j,k)*vy(i,jp,k))
     .       -vis*( gsuper(2,1)*nabla_v(1,2)
     .             +gsuper(2,2)*nabla_v(2,2)
     .             +gsuper(2,3)*nabla_v(3,2) )

        t23 =
     .       0.25*( rvy(i,jp,k)*vz(i,j,k) + rvy(i,j,k)*vz(i,jp,k)
     .             +rvz(i,jp,k)*vy(i,j,k) + rvz(i,j,k)*vy(i,jp,k))
     .       -vis*( gsuper(2,1)*nabla_v(1,3)
     .             +gsuper(2,2)*nabla_v(2,3)
     .             +gsuper(2,3)*nabla_v(3,3) )

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine eom_advc_y

c     eom_advc_z
c     #############################################################
      subroutine eom_advc_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
        real(8)    :: jac,ijac,vis

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

        vis = 2./(1./nuu(i,j,kp) + 1./nuu(i,j,k ))

        t31 =
     .       0.25*( rvz(i,j,kp)*vx(i,j,k) + rvz(i,j,k)*vx(i,j,kp)
     .             +rvx(i,j,kp)*vz(i,j,k) + rvx(i,j,k)*vz(i,j,kp) )
     .       -vis*( gsuper(3,1)*nabla_v(1,1)
     .             +gsuper(3,2)*nabla_v(2,1)
     .             +gsuper(3,3)*nabla_v(3,1) )

        t32 =
     .       0.25*( rvz(i,j,kp)*vy(i,j,k) + rvz(i,j,k)*vy(i,j,kp)
     .             +rvy(i,j,kp)*vz(i,j,k) + rvy(i,j,k)*vz(i,j,kp) )
     .       -vis*( gsuper(3,1)*nabla_v(1,2)
     .             +gsuper(3,2)*nabla_v(2,2)
     .             +gsuper(3,3)*nabla_v(3,2) )

        t33 =
     .       0.25*( rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp)
     .             +rvz(i,j,kp)*vz(i,j,k) + rvz(i,j,k)*vz(i,j,kp) )
     .       -vis*( gsuper(3,1)*nabla_v(1,3)
     .             +gsuper(3,2)*nabla_v(2,3)
     .             +gsuper(3,3)*nabla_v(3,3) )

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine eom_advc_z

c     eom_force_x
c     #############################################################
      subroutine eom_force_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for force term in EOM.
c     In the call sequence:
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
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,ptot

c     Begin program

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
     .               +gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac0 = gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

       !Recall p=a_p nT
        ptot = a_p*0.5*jac*(rho(ip,j,k)*tmp(i ,j,k)
     .                     +rho(i ,j,k)*tmp(ip,j,k))
cc          ptot = jac*(rho(ip,j,k)*tmp(ip,j,k)
cc     .               +rho(i ,j,k)*tmp(i ,j,k))
     .        +jac*(bx(ip,j,k)*bx_cov(i ,j,k)*ijacp
     .             +bx(i ,j,k)*bx_cov(ip,j,k)*ijac0
     .             +by(ip,j,k)*by_cov(i ,j,k)*ijac0
     .             +by(i ,j,k)*by_cov(ip,j,k)*ijacp
     .             +bz(ip,j,k)*bz_cov(i ,j,k)*ijacp
     .             +bz(i ,j,k)*bz_cov(ip,j,k)*ijac0)*0.25
cc     .        +(bx(ip,j,k)*bx_cov(i ,j,k)
cc     .         +bx(i ,j,k)*bx_cov(ip,j,k)
cc     .         +by(ip,j,k)*by_cov(i ,j,k)
cc     .         +by(i ,j,k)*by_cov(ip,j,k)
cc     .         +bz(ip,j,k)*bz_cov(i ,j,k)
cc     .         +bz(i ,j,k)*bz_cov(ip,j,k))*0.25

        t11 = -( bx(ip,j,k)*bx(i ,j,k)*ijacp*ijac0 )*jac**2
     .        +gsuper(1,1)*ptot

        t12 = -0.5*( bx(ip,j,k)*by(i ,j,k)*ijacp
     .              +bx(i ,j,k)*by(ip,j,k)*ijac0 )*jac
     .        +gsuper(1,2)*ptot

        t13 = -0.5*( bx(ip,j,k)*bz(i ,j,k)*ijacp*ijac0
     .              +bx(i ,j,k)*bz(ip,j,k)*ijacp*ijac0)*jac**2
     .        +gsuper(1,3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine eom_force_x

c     eom_force_y
c     #############################################################
      subroutine eom_force_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for force term in EOM.
c     In the call sequence:
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
        real(8)    :: jac,ijac,ptot

c     Begin program

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))

        !Recall p=2nT
cc          ptot = a_p*0.5*jac*(rho(i,jp,k)*tmp(i,jp,k)
cc     .                       +rho(i,j ,k)*tmp(i,j ,k))
        ptot = a_p*0.5*jac*(rho(i,jp,k)*tmp(i,j ,k)
     .                     +rho(i,j ,k)*tmp(i,jp,k))
     .         +(bx(i,jp,k)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,jp,k)
     .          +by(i,jp,k)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,jp,k)
     .          +bz(i,jp,k)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,jp,k))*0.25

        t21 = -0.5*( by(i,jp,k)*bx(i,j ,k)
     .              +by(i,j ,k)*bx(i,jp,k) )
     .        +gsuper(2,1)*ptot

        t22 = -0.5*( by(i,jp,k)*by(i,j ,k)
     .              +by(i,j ,k)*by(i,jp,k) )
     .        +gsuper(2,2)*ptot

        t23 = -0.5*( by(i,jp,k)*bz(i,j ,k)
     .              +by(i,j ,k)*bz(i,jp,k) )
     .        +gsuper(2,3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine eom_force_y

c     eom_force_z
c     #############################################################
      subroutine eom_force_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for force term in EOM.
c     In the call sequence:
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
        real(8)    :: jac,ijac,ptot

c     Begin program

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))

        !Recall p=2nT
cc          ptot = a_p*0.5*jac*(rho(i,j,kp)*tmp(i,j,kp)
cc     .                       +rho(i,j,k )*tmp(i,j,k ))
        ptot = a_p*0.5*jac*(rho(i,j,kp)*tmp(i,j,k )
     .                     +rho(i,j,k )*tmp(i,j,kp))
     .         +(bx(i,j,kp)*bx_cov(i,j,k)+bx(i,j,k)*bx_cov(i,j,kp)
     .          +by(i,j,kp)*by_cov(i,j,k)+by(i,j,k)*by_cov(i,j,kp)
     .          +bz(i,j,kp)*bz_cov(i,j,k)+bz(i,j,k)*bz_cov(i,j,kp))*0.25

        t31 = -0.5*( bz(i,j,kp)*bx(i,j,k )
     .              +bz(i,j,k )*bx(i,j,kp) )
     .        +gsuper(3,1)*ptot

        t32 = -0.5*( bz(i,j,kp)*by(i,j,k )
     .              +bz(i,j,k )*by(i,j,kp) )
     .        +gsuper(3,2)*ptot

        t33 = -0.5*( bz(i,j,kp)*bz(i,j,k )
     .              +bz(i,j,k )*bz(i,j,kp) )
     .        +gsuper(3,3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine eom_force_z

c     si_op
c     #####################################################################
      function si_op(i,j,k,nx,ny,nz,igx,igy,igz) result(cnv)

      implicit none

      integer(4) :: i,j,k,igx,igy,igz,nx,ny,nz
      real(8)    :: cnv(3)

      real(8)    :: psib(3),psit(3)

      !EM part
      psib = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                ,si_tnsr_x,si_tnsr_y,si_tnsr_z,vol=.false.)

      !Pressure part
      call psi_t(i,j,k,igx,igy,igz,psit)

      !SI operator
      cnv = dt**2*k_si*(0*psib + psit)

      end function si_op

c     si_tnsr_x
c     #############################################################
      subroutine si_tnsr_x(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer(4) :: ip,igrid
        real(8)    :: jac,jac0,jacp,gsuper(3,3)
     .               ,acnv(3),acnvp(3),b0cov(3),b0cnv(3),scalar_prod

c     Begin program

        igrid = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
     .               +gmetric%grid(igrid)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))

cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igrid)%jac(ip,j,k)
          jac0 = gmetric%grid(igrid)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        if (flag /= 0) then
          acnv=curl_bxv(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv_dt,b_n,1)
        else
          acnv=curl_bxv(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv_dt,b_n,0)
        endif

        b0cnv(1) = 0.5*(b_n(ip,j,k,1)/jacp+b_n(i,j,k,1)/jac0)*jac
        b0cnv(2) = 0.5*(b_n(ip,j,k,2)     +b_n(i,j,k,2))
        b0cnv(3) = 0.5*(b_n(ip,j,k,3)/jacp+b_n(i,j,k,3)/jac0)*jac

        if (flag /= 0) then
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=1)
        else
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=0)
        endif

        scalar_prod = dot_product(acnv,b0cov)

        t11 =( acnv(1)*b0cnv(1)
     .        +acnv(1)*b0cnv(1)
     .        -gsuper(1,1)*scalar_prod )

        t12 =( acnv(1)*b0cnv(2)
     .        +acnv(2)*b0cnv(1)
     .        -gsuper(1,2)*scalar_prod )

        t13 =( acnv(1)*b0cnv(3)
     .        +acnv(3)*b0cnv(1)
     .        -gsuper(1,3)*scalar_prod )

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        if (flag /= 0) then
          t11 = t11/jac
          if (.not.alt_eom) t12 = t12/jac
          t13 = t13/jac
        endif

c     End program

      end subroutine si_tnsr_x

c     si_tnsr_y
c     #############################################################
      subroutine si_tnsr_y(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EM SI operator
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer(4) :: jp,igrid
        real(8)    :: jac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .                ,scalar_prod

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          acnv=curl_bxv(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv_dt,b_n,2)
        else
          acnv=curl_bxv(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv_dt,b_n,0)
        endif

        b0cnv = 0.5*(b_n(i,jp,k,:)+b_n(i,j,k,:))

        if (flag /= 0) then
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=2)
        else
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=0)
        endif

        scalar_prod = dot_product(acnv,b0cov)

        t21 =( acnv(2)*b0cnv(1)
     .        +acnv(1)*b0cnv(2)
     .        -gsuper(2,1)*scalar_prod )

        t22 =( acnv(2)*b0cnv(2)
     .        +acnv(2)*b0cnv(2)
     .        -gsuper(2,2)*scalar_prod )

        t23 =( acnv(2)*b0cnv(3)
     .        +acnv(3)*b0cnv(2)
     .        -gsuper(2,3)*scalar_prod )

        if (flag /= 0) then
          t21 = t21/jac
          if (.not.alt_eom) t22 = t22/jac
          t23 = t23/jac
        endif


c     End program

      end subroutine si_tnsr_y

c     si_tnsr_z
c     #############################################################
      subroutine si_tnsr_z(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                   ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EM SI operator
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer(4) :: kp,igrid
        real(8)    :: jac,gsuper(3,3),acnv(3),b0cov(3),b0cnv(3)
     .               ,scalar_prod

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          acnv=curl_bxv(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv_dt,b_n,3)
        else
          acnv=curl_bxv(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid,dv_dt,b_n,0)
        endif

        b0cnv = 0.5*(b_n(i,j,kp,:)+b_n(i,j,k,:))

        if (flag /= 0) then
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=3)
        else
          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                        ,b0cov(1),b0cov(2),b0cov(3)
     .                        ,b0cnv(1),b0cnv(2),b0cnv(3)
     .                        ,.false.,half_elem=0)
        endif

        scalar_prod = dot_product(acnv,b0cov)

        t31 =( acnv(3)*b0cnv(1)
     .        +acnv(1)*b0cnv(3)
     .        -gsuper(3,1)*scalar_prod )

        t32 =( acnv(3)*b0cnv(2)
     .        +acnv(2)*b0cnv(3)
     .        -gsuper(3,2)*scalar_prod )

        t33 =( acnv(3)*b0cnv(3)
     .        +acnv(3)*b0cnv(3)
     .        -gsuper(3,3)*scalar_prod )

        if (flag /= 0) then
          t31 = t31/jac
          if (.not.alt_eom) t32 = t32/jac
          t33 = t33/jac
        endif

c     End program

      end subroutine si_tnsr_z

c     psi_t
c     #####################################################################
      subroutine psi_t(i,j,k,igx,igy,igz,psit)

      implicit none

      integer(4) :: i,j,k,igx,igy,igz
      real(8)    :: psit(3)

      integer(4) :: ip,im,jp,jm,kp,km,ig,jg,kg

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
     $             ,divip,divim,divjp,divjm,divkp,divkm
     .             ,cov(3),cnv(3),car(3)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      !Fluxes at faces for calculation of grad[div(dv_dt p_n)]
      flxip =( (dv_dt(ip,j ,k ,1)*p_n(ip,j ,k ,1)
     .         -dv_dt(i ,j ,k ,1)*p_n(i ,j ,k ,1))/dx(ig)
     .        +(dv_dt(ip,jp,k ,2)*p_n(ip,jp,k ,1)
     .         -dv_dt(ip,jm,k ,2)*p_n(ip,jm,k ,1)
     .         +dv_dt(i ,jp,k ,2)*p_n(i ,jp,k ,1)
     .         -dv_dt(i ,jm,k ,2)*p_n(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv_dt(ip,j ,kp,3)*p_n(ip,j ,kp,1)
     .         -dv_dt(ip,j ,km,3)*p_n(ip,j ,km,1)
     .         +dv_dt(i ,j ,kp,3)*p_n(i ,j ,kp,1)
     .         -dv_dt(i ,j ,km,3)*p_n(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacip)*2
      flxim =( (dv_dt(i ,j ,k ,1)*p_n(i ,j ,k ,1)
     .         -dv_dt(im,j ,k ,1)*p_n(im,j ,k ,1))/dx(ig-1)
     .        +(dv_dt(im,jp,k ,2)*p_n(im,jp,k ,1)
     .         -dv_dt(im,jm,k ,2)*p_n(im,jm,k ,1)
     .         +dv_dt(i ,jp,k ,2)*p_n(i ,jp,k ,1)
     .         -dv_dt(i ,jm,k ,2)*p_n(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv_dt(im,j ,kp,3)*p_n(im,j ,kp,1)
     .         -dv_dt(im,j ,km,3)*p_n(im,j ,km,1)
     .         +dv_dt(i ,j ,kp,3)*p_n(i ,j ,kp,1)
     .         -dv_dt(i ,j ,km,3)*p_n(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacim)*2

      flxjp =( (dv_dt(ip,jp,k ,1)*p_n(ip,jp,k ,1)
     .         -dv_dt(im,jp,k ,1)*p_n(im,jp,k ,1)
     .         +dv_dt(ip,j ,k ,1)*p_n(ip,j ,k ,1)
     .         -dv_dt(im,j ,k ,1)*p_n(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv_dt(i ,jp,k ,2)*p_n(i ,jp,k ,1)
     .         -dv_dt(i ,j ,k ,2)*p_n(i ,j ,k ,1))/dy(jg)
     .        +(dv_dt(i ,jp,kp,3)*p_n(i ,jp,kp,1)
     .         -dv_dt(i ,jp,km,3)*p_n(i ,jp,km,1)
     .         +dv_dt(i ,j ,kp,3)*p_n(i ,j ,kp,1)
     .         -dv_dt(i ,j ,km,3)*p_n(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacjp)*2
      flxjm =( (dv_dt(ip,jm,k ,1)*p_n(ip,jm,k ,1)
     .         -dv_dt(im,jm,k ,1)*p_n(im,jm,k ,1)
     .         +dv_dt(ip,j ,k ,1)*p_n(ip,j ,k ,1)
     .         -dv_dt(im,j ,k ,1)*p_n(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv_dt(i ,j ,k ,2)*p_n(i ,j ,k ,1)
     .         -dv_dt(i ,jm,k ,2)*p_n(i ,jm,k ,1))/dy(jg-1)
     .        +(dv_dt(i ,jm,kp,3)*p_n(i ,jm,kp,1)
     .         -dv_dt(i ,jm,km,3)*p_n(i ,jm,km,1)
     .         +dv_dt(i ,j ,kp,3)*p_n(i ,j ,kp,1)
     .         -dv_dt(i ,j ,km,3)*p_n(i ,j ,km,1))/dzh(kg)/4.)
     $        /(jac+jacjm)*2

      flxkp =( (dv_dt(ip,j ,kp,1)*p_n(ip,j ,kp,1)
     .         -dv_dt(im,j ,kp,1)*p_n(im,j ,kp,1)
     .         +dv_dt(ip,j ,k ,1)*p_n(ip,j ,k ,1)
     .         -dv_dt(im,j ,k ,1)*p_n(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv_dt(i ,jp,kp,2)*p_n(i ,jp,kp,1)
     .         -dv_dt(i ,jm,kp,2)*p_n(i ,jm,kp,1)
     .         +dv_dt(i ,jp,k ,2)*p_n(i ,jp,k ,1)
     .         -dv_dt(i ,jm,k ,2)*p_n(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv_dt(i ,j ,kp,3)*p_n(i ,j ,kp,1)
     .         -dv_dt(i ,j ,k ,3)*p_n(i ,j ,k ,1))/dz(kg) )
     $        /(jac+jackp)*2
      flxkm =( (dv_dt(ip,j ,km,1)*p_n(ip,j ,km,1)
     .         -dv_dt(im,j ,km,1)*p_n(im,j ,km,1)
     .         +dv_dt(ip,j ,k ,1)*p_n(ip,j ,k ,1)
     .         -dv_dt(im,j ,k ,1)*p_n(im,j ,k ,1))/dxh(ig)/4.
     .        +(dv_dt(i ,jp,km,2)*p_n(i ,jp,km,1)
     .         -dv_dt(i ,jm,km,2)*p_n(i ,jm,km,1)
     .         +dv_dt(i ,jp,k ,2)*p_n(i ,jp,k ,1)
     .         -dv_dt(i ,jm,k ,2)*p_n(i ,jm,k ,1))/dyh(jg)/4.
     .        +(dv_dt(i ,j ,k ,3)*p_n(i ,j ,k ,1)
     .         -dv_dt(i ,j ,km,3)*p_n(i ,j ,km,1))/dz(kg-1) )
     $        /(jac+jackm)*2

      !Fluxes at faces for calculation of grad[(gamma-1)*p_n*div(dv_dt)]

      !!Divergence at faces i+-1/2, etc.
      divip = (dv_dt(ip,j ,k,1)-dv_dt(i ,j ,k,1))/dx(ig)
     .       +(dv_dt(i ,jp,k,2)-dv_dt(i ,jm,k,2)
     .        +dv_dt(ip,jp,k,2)-dv_dt(ip,jm,k,2))/dyh(jg)/4.
     .       +(dv_dt(i ,j,kp,3)-dv_dt(i ,j,km,3)
     .        +dv_dt(ip,j,kp,3)-dv_dt(ip,j,km,3))/dzh(kg)/4.
      divim = (dv_dt(i ,j ,k,1)-dv_dt(im,j ,k,1))/dx(ig-1)
     .       +(dv_dt(i ,jp,k,2)-dv_dt(i ,jm,k,2)
     .        +dv_dt(im,jp,k,2)-dv_dt(im,jm,k,2))/dyh(jg)/4.
     .       +(dv_dt(i ,j,kp,3)-dv_dt(i ,j,km,3)
     .        +dv_dt(im,j,kp,3)-dv_dt(im,j,km,3))/dzh(kg)/4.

      divjp = (dv_dt(ip,j ,k,1)-dv_dt(im,j ,k,1)
     .        +dv_dt(ip,jp,k,1)-dv_dt(im,jp,k,1))/dxh(ig)/4.
     .       +(dv_dt(i ,jp,k,2)-dv_dt(i ,j ,k,2))/dy(jg)
     .       +(dv_dt(i,j ,kp,3)-dv_dt(i,j ,km,3)
     .        +dv_dt(i,jp,kp,3)-dv_dt(i,jp,km,3))/dzh(kg)/4.
      divjm = (dv_dt(ip,j ,k,1)-dv_dt(im,j ,k,1)
     .        +dv_dt(ip,jm,k,1)-dv_dt(im,jm,k,1))/dxh(ig)/4.
     .       +(dv_dt(i ,j ,k,2)-dv_dt(i ,jm,k,2))/dy(jg-1)
     .       +(dv_dt(i,j ,kp,3)-dv_dt(i,j ,km,3)
     .        +dv_dt(i,jm,kp,3)-dv_dt(i,jm,km,3))/dzh(kg)/4.

      divkp = (dv_dt(ip,j,k ,1)-dv_dt(im,j,k ,1)
     .        +dv_dt(ip,j,kp,1)-dv_dt(im,j,kp,1))/dxh(ig)/4.
     .       +(dv_dt(i,jp,k ,2)-dv_dt(i,jm,k ,2)
     .        +dv_dt(i,jp,kp,2)-dv_dt(i,jm,kp,2))/dyh(jg)/4.
     .       +(dv_dt(i,j ,kp,3)-dv_dt(i,j ,k ,3))/dz(kg)
      divkm = (dv_dt(ip,j,k ,1)-dv_dt(im,j,k ,1)
     .        +dv_dt(ip,j,km,1)-dv_dt(im,j,km,1))/dxh(ig)/4.
     .       +(dv_dt(i,jp,k ,2)-dv_dt(i,jm,k ,2)
     .        +dv_dt(i,jp,km,2)-dv_dt(i,jm,km,2))/dyh(jg)/4.
     .       +(dv_dt(i,j ,k ,3)-dv_dt(i,j ,km,3))/dz(kg-1)

      flxip = flxip
     .      + (gamma-1.)*(p_n(i,j,k,1)+p_n(ip,j,k,1))*divip/(jac+jacip)
      flxim = flxim
     .      + (gamma-1.)*(p_n(i,j,k,1)+p_n(im,j,k,1))*divim/(jac+jacim)

      flxjp = flxjp
     .      + (gamma-1.)*(p_n(i,j,k,1)+p_n(i,jp,k,1))*divjp/(jac+jacjp)
      flxjm = flxjm
     .      + (gamma-1.)*(p_n(i,j,k,1)+p_n(i,jm,k,1))*divjm/(jac+jacjm)

      flxkp = flxkp
     .      + (gamma-1.)*(p_n(i,j,k,1)+p_n(i,j,kp,1))*divkp/(jac+jackp)
      flxkm = flxkm
     .      + (gamma-1.)*(p_n(i,j,k,1)+p_n(i,j,km,1))*divkm/(jac+jackm)

      if (isSP(i,j,k,igx,igy,igz)) flxim = 0d0

      cov(1) = (flxip - flxim)/dxh(ig)
      cov(2) = (flxjp - flxjm)/dyh(jg)
      cov(3) = (flxkp - flxkm)/dzh(kg)

      call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                            ,cov(1),cov(2),cov(3)
     .                            ,cnv(1),cnv(2),cnv(3),.true.)
      psit = -cnv

      end subroutine psi_t

      end module nlfunction_setup

c module vectorOps
c #########################################################################
      module vectorOps

      contains

c     transformVector
c     ######################################################################
      subroutine transformVector(igx,igy,igz
     .                          ,imin,imax,jmin,jmax,kmin,kmax
     .                          ,arr1,arr2,arr3,covariant,to_cartsn)

c     ----------------------------------------------------------------------
c     Transforms vectors components in arrays arr1,arr2,arr3
c     from Cartesian to curvilinear (to_cartesian=.false.)
c     and viceversa, either with covariant (covariant=.true.) or 
c     contravariant curvilinear vectors.
c     ----------------------------------------------------------------------

      use grid

      implicit none

c     Input variables

        integer(4) :: imin,imax,jmin,jmax,kmin,kmax
        integer(4) :: igx,igy,igz
        logical    :: covariant,to_cartsn
        real(8)    :: arr1(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr2(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr3(imin:imax,jmin:jmax,kmin:kmax)

c     Local variables

        integer(4) :: i,j,k
        real(8)    :: vec(3)

c     Begin program

        if (to_cartsn) then

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call transformVectorToCartesian
     .               (i,j,k,igx,igy,igz
     .               ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .               ,covariant
     .               ,vec(1),vec(2),vec(3))

                arr1(i,j,k) = vec(1)
                arr2(i,j,k) = vec(2)
                arr3(i,j,k) = vec(3)
                
              enddo
            enddo
          enddo

        else

          do k=kmin,kmax
            do j=jmin,jmax
              do i=imin,imax

                call transformVectorToCurvilinear
     .               (i,j,k,igx,igy,igz
     .               ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .               ,covariant
     .               ,vec(1),vec(2),vec(3))

                arr1(i,j,k) = vec(1)
                arr2(i,j,k) = vec(2)
                arr3(i,j,k) = vec(3)
                
              enddo
            enddo
          enddo

        endif

      end subroutine transformVector

      end module vectorOps
