c module equilibrium
c ######################################################################
      module equilibrium

        integer    :: IRHO,IVX,IVY,IVZ,IBX,IBY,IBZ,ITMP,IJX,IJY,IJZ
     .               ,IAX,IAY,IAZ,IBCNV,IBCOV,IJCNV,IJCOV

#if !defined(vec_pot)
        parameter    (IRHO=1,IVX=2,IVY=3,IVZ=4,IBX=5,IBY=6,IBZ=7,ITMP=8
     .               ,IJX=9,IJY=10,IJZ=11,IAX=12,IAY=13,IAZ=14
     .               ,IJCNV=1,IJCOV=2,IBCNV=3,IBCOV=4)
#else
        parameter    (IRHO=1,IVX=2,IVY=3,IVZ=4,IAX=5,IAY=6,IAZ=7,ITMP=8
     .               ,IBX=9,IBY=10,IBZ=11,IJX=12,IJY=13,IJZ=14
     .               ,IJCNV=1,IJCOV=2,IBCNV=3,IBCOV=4)
#endif

        real(8)    :: dlambda,rshear,vparflow,vperflow

        character(5)  :: equil
        character(20) :: equ_file,prt_file
        real(8)       :: eq_params(6)

        logical       :: dcon=.false.

      end module equilibrium

c module transport_params
c ######################################################################
      module transport_params

        use grid

        use equilibrium

        real(8) :: nu,eta,dd,chi,gamma,E0(3)=0d0,B0(3)=0d0,di,a_p
     .            ,heta=0d0,aa_eta,c_hyper

        integer :: nn_eta

        logical :: adiabatic,hall_new=.false.

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
      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

c     Local variables

      integer    :: ig,jg,kg
      real(8)    :: x1,y1,z1
      logical    :: cartsn

c     Begin program

c     Resistivity profile eta*(1 + aa*x^nn)
c       Coefficient aa is set so that res = 20*eta at wall
c       and nn so that res=??*eta at sing. surf. xs ~ 0.33

      select case (equil)
      case ('ppnch','ppnsl','ppnst')

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
        res = eta*(1. + aa_eta*grid_params%xx(ig)**nn_eta)

      case default

        res = eta

      end select

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
      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

c     Local variables

c     Begin program

      vis = nu

c     End program

      end function vis

c     hres
c     ###############################################################
      function hres(nx,ny,nz,igx,igy,igz)
c     ---------------------------------------------------------------
c     Finds hyper-resistivity coefficient (global on the grid)
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: hres
      integer    :: nx,ny,nz,igx,igy,igz

c     Local variables

      real(8)    :: kk,norm,idx,idy,idz
      integer    :: ig,jg,kg,i,j,k

c     Begin program

      kk = 0d0

cc      do k=1,nz
cc        do j=1,ny
cc          do i=1,nx
cc
cc            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc            !Maximum kk
cc            idx  = 1./grid_params%dx(ig)
cc            if (nx == 1) idx = 1d-2
cc            idy  = 1./grid_params%dy(jg)
cc            if (ny == 1) idy = 1d-2
cc            idz  = 1./grid_params%dz(kg)
cc            if (nz == 1) idz = 1d-2
cc
cc            norm = vectorNorm(i,j,k,igx,igy,igz,idx,idy,idz,.true.)
cc            kk   = max(kk,norm)
cc
cc          enddo
cc        enddo
cc      enddo
cc
cc      idx = grid_params%nxgl(igx)/(xmax-xmin)
cc      idy = grid_params%nxgl(igy)/(ymax-ymin)
cc      idz = grid_params%nxgl(igz)/(zmax-zmin)
cc
cc      kk = (idx**2 + idy**2 + idz**2)
cc
cc      hres = .05*di/kk
cc      hres = di/kk
      hres = c_hyper

c     End

      end function hres

      end module transport_params

c module grid_aliases
c ######################################################################
      module grid_aliases

        use grid

        real(8),pointer,dimension(:) :: xx,yy,zz,dxh,dyh,dzh,dx,dy,dz

        integer    :: igx,igy,igz,nx,ny,nz

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

        use variables

        use equilibrium

#if defined(vec_pot)
        real(8),target,allocatable,dimension(:,:,:,:) :: acov,acnv
     .                                                  ,vlap_a
        real(8),pointer,dimension(:,:,:,:) :: bcnv,bcov
#else
        real(8),target,allocatable,dimension(:,:,:,:) :: bcnv,bcov
#endif

        real(8),target,allocatable,dimension(:,:,:) :: eeta,nuu,divrgV

        real(8),target,allocatable,dimension(:,:,:,:) ::vcnv,vcov
     .                                                 ,pcnv,vecnv
     .                                                 ,vdummy,E_ni
     .                                                 ,vlap_v
     .                                                 ,jcov_0

        real(8),pointer,dimension(:,:,:,:) :: jcnv,jcov

        logical :: alt_eom

        logical :: nc_eom_jxb=.false.,nc_eom_gp=.false.,nc_eom_v=.false.
     .            ,solve_rho=.true.

c SI operator
        real(8) :: k_si=0d0
        real(8),target,allocatable,dimension(:,:,:,:) :: b_n,p_n,dv_dt
     .                                                  ,vcov_n
c SI operator

      contains

c     allocAuxVariables
c     ###################################################################
      subroutine allocAuxVariables

c     -------------------------------------------------------------------
c     Allocates auxiliary variable storage.
c     -------------------------------------------------------------------

c     Begin program

c     Allocate storage

        if (.not.associated(gv%aux)) then
#if defined(vec_pot)
          call allocateAuxStruct(0,4,gv%aux)
#else
          call allocateAuxStruct(0,2,gv%aux)
#endif
        endif

        gv%aux%vec_list(IJCNV)%cnv = .true.
        jcnv => gv%aux%vec_list(IJCNV)%vec

        gv%aux%vec_list(IJCOV)%cnv = .false.
        jcov => gv%aux%vec_list(IJCOV)%vec

#if defined(vec_pot)
        gv%aux%vec_list(IBCNV)%cnv = .true.
        bcnv => gv%aux%vec_list(IBCNV)%vec

        gv%aux%vec_list(IBCOV)%cnv = .false.
        bcov => gv%aux%vec_list(IBCOV)%vec
#endif

c     End program

      end subroutine allocAuxVariables

c     defineAuxBCs
c     ###################################################################
cc      subroutine defineAuxBCs(neq,bbcs)
      subroutine defineAuxBCs

c     -------------------------------------------------------------------
c     Defines BCs for auxiliary quantities
c     -------------------------------------------------------------------

c     Call variables

cc      integer    :: neq,bbcs(6,neq)

c     Local variables

      integer    :: ieq,bcsq(6)

c     Begin program

cc#if defined(vec_pot)
cc      !Magnetic field
cc      bcsq = bcond
cc      where (bcsq == DEF) bcsq = DIR
cc      gv%aux%vec_list(IBCNV)%bconds(:,1) = bcsq
cc
cccc      if (equil == 'spnch' .or. equil == 'rfp2') then
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cccc        gv%aux%vec_list(IBCNV)%bconds(:,2) = bcsq
cccc
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cccc        gv%aux%vec_list(IBCNV)%bconds(:,3) = bcsq
cccc      else
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -NEU !On covariant components
cc        gv%aux%vec_list(IBCNV)%bconds(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -NEU !On covariant components
cc        gv%aux%vec_list(IBCNV)%bconds(:,3) = bcsq
cccc      endif
cc
cc      gv%aux%vec_list(IBCOV)%bconds = gv%aux%vec_list(IBCNV)%bconds
cc
cc      !Current
cccc      if (equil == 'spnch' .or. equil == 'rfp2') then
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = DIR
cccc        bcnd(:,1) = bcsq
cccc
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = EXT
cccccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cccc        bcnd(:,2) = bcsq
cccc
cccc        bcsq = bcond
cccc        where (bcsq == DEF) bcsq = EXT
cccccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cccc        bcnd(:,3) = bcsq
cccc      else
cc        gv%aux%vec_list(IJCNV)%bconds(:,1) = bbcs(:,IAX)
cc        gv%aux%vec_list(IJCNV)%bconds(:,2) = bbcs(:,IAY)
cc        gv%aux%vec_list(IJCNV)%bconds(:,3) = bbcs(:,IAZ)
cccc      endif
cc
cc      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
cc#else
cc      !Current
cc      gv%aux%vec_list(IJCNV)%bconds = bbcs(:,IBX:IBZ)
cc      where (gv%aux%vec_list(IJCNV)%bconds == -NEU)
cccc        gv%aux%vec_list(IJCNV)%bconds = -DIR  !Use covariant components for tangential BCs
cc        gv%aux%vec_list(IJCNV)%bconds = -EQU  !Use covariant components for tangential BCs
cc      end where
cc
cc      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
cc#endif

#if defined(vec_pot)
      !Magnetic field
      bcsq = bcond
      where (bcsq == DEF) bcsq = DIR
      gv%aux%vec_list(IBCNV)%bconds(:,1) = bcsq

cc      if (equil == 'spnch' .or. equil == 'rfp2') then
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cc        gv%aux%vec_list(IBCNV)%bconds(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cc        gv%aux%vec_list(IBCNV)%bconds(:,3) = bcsq
cc      else
        bcsq = bcond
        where (bcsq == DEF) bcsq = -NEU !On covariant components
        gv%aux%vec_list(IBCNV)%bconds(:,2) = bcsq

        bcsq = bcond
        where (bcsq == DEF) bcsq = -NEU !On covariant components
        gv%aux%vec_list(IBCNV)%bconds(:,3) = bcsq
cc      endif

      gv%aux%vec_list(IBCOV)%bconds = gv%aux%vec_list(IBCNV)%bconds

      !Current
cc      if (equil == 'spnch' .or. equil == 'rfp2') then
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = DIR
cc        bcnd(:,1) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = EXT
cccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc        bcnd(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = EXT
cccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc        bcnd(:,3) = bcsq
cc      else
        gv%aux%vec_list(IJCNV)%bconds(:,1) = u_0%array_var(IAX)%bconds
        gv%aux%vec_list(IJCNV)%bconds(:,2) = u_0%array_var(IAY)%bconds
        gv%aux%vec_list(IJCNV)%bconds(:,3) = u_0%array_var(IAZ)%bconds
cc      endif

      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
#else
      !Current
      gv%aux%vec_list(IJCNV)%bconds(:,1) = u_0%array_var(IBX)%bconds
      gv%aux%vec_list(IJCNV)%bconds(:,2) = u_0%array_var(IBY)%bconds
      gv%aux%vec_list(IJCNV)%bconds(:,3) = u_0%array_var(IBZ)%bconds

      where (gv%aux%vec_list(IJCNV)%bconds == -NEU)
        gv%aux%vec_list(IJCNV)%bconds = -EQU  !Use covariant components for tangential BCs
      end where

      gv%aux%vec_list(IJCOV)%bconds = gv%aux%vec_list(IJCNV)%bconds
#endif

c     End program

      end subroutine defineAuxBCs

      end module auxiliaryVariables

c module auxPlotVariables
c ######################################################################
      module auxPlotVariables

#if defined(vec_pot)
        real(8),target,allocatable,dimension(:,:,:) ::
     .          ax_car,ay_car,az_car
     .         ,ax_cnv,ay_cnv,az_cnv
     .         ,bx,by,bz
#endif

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_car,by_car,bz_car
     .         ,jx_car,jy_car,jz_car
     .         ,vx_car,vy_car,vz_car
     .         ,divrgB,divrgJ,Pflux
     .         ,qfactor,lambda,qfactr2,lambda2,p_tot
     .         ,vpar,vpsi,vperp

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_cov,by_cov,bz_cov
     .         ,jx,jy,jz,jx_cov,jy_cov,jz_cov
     .         ,vx,vy,vz,vx_cov,vy_cov,vz_cov
     .         ,vex,vey,vez

      end module auxPlotVariables

c module operators
c ######################################################################
      module operators

        use grid

        use grid_aliases

        use error

        use transport_params

        real(8),pointer,dimension(:,:,:,:) :: vec,vec1,vec2
        real(8),pointer,dimension(:,:,:)   :: coef,sc1,sc2

        logical    :: solenoidal=.true.

        INTERFACE veclaplacian
          module procedure veclap_diff,veclap_ndiff,veclap2
        end INTERFACE

        INTERFACE veclap_cov
          module procedure veclap_cov_diff,veclap_cov_ndiff
        end INTERFACE

        INTERFACE laplacian
          module procedure lap_diff,lap_ndiff
        end INTERFACE

        INTERFACE curlcurl
          module procedure curlcurl_diff,curlcurl_ndiff
        end INTERFACE

        PRIVATE :: grad_vrtx,curl_vrtx,div_vrtx

      contains

c     div_vrtx
c     ##################################################################
      function div_vrtx(i,j,k,igx,igy,igz,vx,vy,vz,vrtx) result(div)

c     ------------------------------------------------------------------
c     Calculates div at vertex at vertex (i+1/2,j+1/2,k+1/2) from cnv
c     vector defined at surrounding cell centers.
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igx,igy,igz
      real(8) :: vx(2,2,2),vy(2,2,2),vz(2,2,2),div
      logical :: vrtx

c     Local variables

      integer :: ig,jg,kg,ip,i0,jp,j0,kp,k0
      real(8) :: dhx,dhy,dhz,fp,fm

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      if (vrtx) then
        dhx = dx(ig)
        dhy = dy(jg)
        dhz = dz(kg)
      else
        dhx = dxh(ig)
        dhy = dyh(jg)
        dhz = dzh(kg)
      endif

      ip = 2
      i0 = 1
      jp = 2
      j0 = 1
      kp = 2
      k0 = 1

      !X component
      fp = 0.25*(vx(ip,j0,k0)+vx(ip,jp,k0)
     .          +vx(ip,j0,kp)+vx(ip,jp,kp))
      fm = 0.25*(vx(i0,j0,k0)+vx(i0,jp,k0)
     .          +vx(i0,j0,kp)+vx(i0,jp,kp))

      div = (fp-fm)/dhx

      !Y component
      fp = 0.25*(vy(i0,jp,k0)+vy(ip,jp,k0)
     .          +vy(i0,jp,kp)+vy(ip,jp,kp))
      fm = 0.25*(vy(i0,j0,k0)+vy(ip,j0,k0)
     .          +vy(i0,j0,kp)+vy(ip,j0,kp))

      div = div + (fp-fm)/dhy

      !Z component
      fp = 0.25*(vz(ip,j0,kp)+vz(i0,j0,kp)
     .          +vz(ip,jp,kp)+vz(i0,jp,kp))
      fm = 0.25*(vz(ip,j0,k0)+vz(i0,j0,k0)
     .          +vz(ip,jp,k0)+vz(i0,jp,k0))

      div = div + (fp-fm)/dhz

      !Jacobian factor
      if (vrtx) then
        jac = 0.125*(gmetric%grid(igx)%jac(i  ,j  ,k  )
     .              +gmetric%grid(igx)%jac(i+1,j  ,k  )
     .              +gmetric%grid(igx)%jac(i  ,j+1,k  )
     .              +gmetric%grid(igx)%jac(i+1,j+1,k  ) 
     .              +gmetric%grid(igx)%jac(i  ,j  ,k+1)
     .              +gmetric%grid(igx)%jac(i+1,j  ,k+1)
     .              +gmetric%grid(igx)%jac(i  ,j+1,k+1)
     .              +gmetric%grid(igx)%jac(i+1,j+1,k+1))
      else
        jac = gmetric%grid(igx)%jac(i,j,k)
      endif

      if (isSP(i+1,j,k,igx,igy,igz)) then
        div = 0d0
      else
        div = div/jac
      endif

c     End program

      end function div_vrtx

c     grad_vrtx
c     ##################################################################
      function grad_vrtx(i,j,k,igx,igy,igz,phi,vrtx) result(grd)

c     ------------------------------------------------------------------
c     Calculates gradient at vertex (i+1/2,j+1/2,k+1/2) from scalar
c     defined at surrounding cell centers
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igx,igy,igz
      real(8) :: phi(2,2,2),grd(3)
      logical :: vrtx

c     Local variables

      integer :: ig,jg,kg,ip,i0,jp,j0,kp,k0,nx,ny,nz
      real(8) :: idhx,idhy,idhz,fp,fm

c     Begin program

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      nx = grid_params%nxgl(igx)
      ny = grid_params%nygl(igy)
      nz = grid_params%nzgl(igz)

      idhx = 0d0
      idhy = 0d0
      idhz = 0d0

      if (vrtx) then
        if (nx > 1) idhx = 1./dx(ig)
        if (ny > 1) idhy = 1./dy(jg)
        if (nz > 1) idhz = 1./dz(kg)
      else
        if (nx > 1) idhx = 1./dxh(ig)
        if (ny > 1) idhy = 1./dyh(jg)
        if (nz > 1) idhz = 1./dzh(kg)
      endif

cc      if (vrtx) then
cc        idhx = 1./dx(ig)
cc        idhy = 1./dy(jg)
cc        idhz = 1./dz(kg)
cc      else
cc        idhx = 1./dxh(ig)
cc        idhy = 1./dyh(jg)
cc        idhz = 1./dzh(kg)
cc      endif
cc
cc      if (vrtx) then
cc        ip = 1 + min(grid_params%nxgl(igx)-(i+grid_params%ilo(igx)-1),1)
cc        i0 = 2 - min(i+grid_params%ilo(igx)-1,1)
cc        jp = 1 + min(grid_params%nygl(igy)-(j+grid_params%jlo(igy)-1),1)
cc        j0 = 2 - min(j+grid_params%jlo(igy)-1,1)
cc        kp = 1 + min(grid_params%nzgl(igz)-(k+grid_params%klo(igz)-1),1)
cc        k0 = 2 - min(k+grid_params%klo(igz)-1,1)
cc
cc        write (*,*) i,j,k,"  ;  ",i0,ip,j0,jp,k0,kp
cc
cc      else
        ip = 2
        i0 = 1
        jp = 2
        j0 = 1
        kp = 2
        k0 = 1
cc      endif

      !X component
      fp = 0.25*(phi(ip,j0,k0)+phi(ip,jp,k0)
     .          +phi(ip,j0,kp)+phi(ip,jp,kp))
      fm = 0.25*(phi(i0,j0,k0)+phi(i0,jp,k0)
     .          +phi(i0,j0,kp)+phi(i0,jp,kp))

      grd(1) = (fp-fm)*idhx

      !Y component
      fp = 0.25*(phi(i0,jp,k0)+phi(ip,jp,k0)
     .          +phi(i0,jp,kp)+phi(ip,jp,kp))
      fm = 0.25*(phi(i0,j0,k0)+phi(ip,j0,k0)
     .          +phi(i0,j0,kp)+phi(ip,j0,kp))

      grd(2) = (fp-fm)*idhy

      !Z component
      fp = 0.25*(phi(ip,j0,kp)+phi(i0,j0,kp)
     .          +phi(ip,jp,kp)+phi(i0,jp,kp))
      fm = 0.25*(phi(ip,j0,k0)+phi(i0,j0,k0)
     .          +phi(ip,jp,k0)+phi(i0,jp,k0))

      grd(3) = (fp-fm)*idhz

c     End program

      end function grad_vrtx

c     curl_vrtx
c     ##################################################################
      function curl_vrtx(i,j,k,igx,igy,igz,vv,vrtx) result(cnv)

c     ------------------------------------------------------------------
c     Calculates curl at vertex at vertex (i+1/2,j+1/2,k+1/2) from cov
c     vector defined at surrounding cell centers.
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igx,igy,igz
      real(8) :: vv(2,2,2,3),cnv(3)
      logical :: vrtx

c     Local variables

      real(8) :: gvx(3),gvy(3),gvz(3)

c     Begin program

      gvx = grad_vrtx(i,j,k,igx,igy,igz,vv(:,:,:,1),vrtx)
      gvy = grad_vrtx(i,j,k,igx,igy,igz,vv(:,:,:,2),vrtx)
      gvz = grad_vrtx(i,j,k,igx,igy,igz,vv(:,:,:,3),vrtx)

      cnv(1) = gvz(2)-gvy(3)
      cnv(2) = gvx(3)-gvz(1)
      cnv(3) = gvy(1)-gvx(2)

c     End program

      end function curl_vrtx

c     div
c     ###############################################################
      function div(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az,he,sp)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field at cell centers in
c     general non-orthogonal geometry.
c     ---------------------------------------------------------------

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1),div
      integer   ,optional :: he
      logical,optional    :: sp

c     Local variables

      integer    :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km
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
     .         +ax(i  ,j  ,k  )/jac0)*jach     )/dxx
     .        +(ay(i  ,j+1,k  )-ay(i  ,j-1,k  ))/dyy
     .        +(az(i  ,j  ,k+1)-az(i  ,j  ,k-1))/dzz
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

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
     .             ,v1 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v2 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v3 (0:nx+1,0:ny+1,0:nz+1),advec

      logical,optional :: vol,upwind

c     Local variables

      integer    :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km,su
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

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)    :: phi(0:nx+1,0:ny+1,0:nz+1)
     .             ,v1 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v2 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v3 (0:nx+1,0:ny+1,0:nz+1),c_advec

      logical,optional :: sp,vol,upwind

c     Local variables

      integer    :: ig,jg,kg,igrid,half_elem,ip,im,jp,jm,kp,km,su
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
          if (isSP(i,j,k,igx,igy,igz)) jach = 0d0
          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
        elseif (i+grid_params%ilo(igx)-1 == grid_params%nxgl(igx)) then
          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))

          jach = 0.5*(jac+jacim)
          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
        endif

cc        if (i+grid_params%ilo(igx)-1 == 1) then
cc          jach = 0.5*(jac+jacip)
cc          flxip = 0.5*(v1(ip,j,k)*phi(i ,j,k)/jacip
cc     .               + v1(i ,j,k)*phi(ip,j,k)/jac  )*jach
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cc        elseif (i+grid_params%ilo(igx)-1 == 2) then
cc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc
cc          jach = 0.5*(jac+jacim)
cc          flxim = 0.5*(v1(im,j,k)*phi(i ,j,k)/jacim
cc     .               + v1(i ,j,k)*phi(im,j,k)/jac  )*jach
cc        else
cc          flxip = 0.5*(v1(ip,j,k)*phi(i,j,k) + v1(i,j,k)*phi(ip,j,k))
cc          flxim = 0.5*(v1(im,j,k)*phi(i,j,k) + v1(i,j,k)*phi(im,j,k))
cc        endif

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
      
      if (vw) c_advec = c_advec*gmetric%grid(igx)%dvol(i,j,k)

c     End 

      end function c_advec

c     flx_advec
c     ###############################################################
      function flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,v1,v2,v3,arr
     .                  ,method,sp,vol) result(advec)
      implicit none
c     ---------------------------------------------------------------
c     Calculates divergence of vector field (v.arr) at cell centers in
c     general non-orthogonal geometry.
c
c     method: specifies the advective scheme:
c        1, 7 -> upwind,
c        2 -> ZIP,
c        3 -> QUICK,
c        4 -> SMART,
c        5 -> smooth SMART,
c        6 -> centered high-order
c        8 -> Gamma
c     ---------------------------------------------------------------

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz,method

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1)
     .             ,v1 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v2 (0:nx+1,0:ny+1,0:nz+1)
     .             ,v3 (0:nx+1,0:ny+1,0:nz+1),advec

      logical,optional :: sp,vol

c     Local variables

      integer    :: igrid,half_elem,su

      integer    :: im,ip,ipp,imm,igm,igp,igpp,igmm
     .             ,jm,jp,jpp,jmm,jgm,jgp,jgpp,jgmm
     .             ,km,kp,kpp,kmm,kgm,kgp,kgpp,kgmm
     .             ,ig,jg,kg,igl,jgl,kgl

      real(8)    :: dxx,dyy,dzz,x0,y0,z0,jacip,jacim,jac,jach
      real(8)    :: vxp,vxm,vyp,vym,vzp,vzm
     .             ,qip,qim,qjp,qjm,qkp,qkm
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
      logical    :: vw,spoint

c     Begin program

c     Set defaults

      if (PRESENT(vol)) then
        vw = vol
      else
        vw = .true.
      endif

      if (PRESENT(sp)) then
        spoint = sp
      else
        spoint = .false.
      endif

c     General setup

      igrid = igx

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      ipp = ip+1
      imm = im-1
cc      if (bcond(1) == PER .and. method > 2) then
cc        wrap_x = .true.
cc        if (i == nx) ipp = 2
cc        if (i == 1 ) imm = nx-1
cc      else
cc        wrap_x = .false.
        ipp = min(ipp,nx+1)
        imm = max(imm,0)
cc      endif

      igp  = ig+1
      igm  = ig-1
      igpp = ig + (ipp - i)
      igmm = ig + (imm - i)
cc      if (i == nx .and. wrap_x)
cc     .     call getMGmap(1   ,j,k,igx,igy,igz,igpp,jg,kg)
cc      if (i == 1  .and. wrap_x)
cc     .     call getMGmap(nx-1,j,k,igx,igy,igz,igmm,jg,kg)

      jpp = jp+1
      jmm = jm-1
cc      if (bcond(3) == PER .and. method > 2) then
cc        wrap_y = .true.
cc        if (j == ny) jpp = 2
cc        if (j == 1 ) jmm = ny-1
cc      else
cc        wrap_y = .false.
        jpp = min(jpp,ny+1)
        jmm = max(jmm,0)
cc      endif

      jgp  = jg+1
      jgm  = jg-1
      jgpp = jg + (jpp-j)
      jgmm = jg + (jmm-j)
cc      if (j == ny .and. wrap_y)
cc     .     call getMGmap(i,1   ,k,igx,igy,igz,ig,jgpp,kg)
cc      if (j == 1  .and. wrap_y)                 
cc     .     call getMGmap(i,ny-1,k,igx,igy,igz,ig,jgmm,kg)

      kpp = kp+1
      kmm = km-1
cc      if (bcond(3) == PER .and. method > 2) then
cc        wrap_z = .true.
cc        if (k == nz) kpp = 2
cc        if (k == 1 ) kmm = nz-1
cc      else
cc        wrap_z = .false.
        kpp = min(kpp,nz+1)
        kmm = max(kmm,0)
cc      endif

      kgp  = kg+1
      kgm  = kg-1
      kgpp = kg + (kpp-k)
      kgmm = kg + (kmm-k)
cc      if (k == nz .and. wrap_z)
cc     .     call getMGmap(i,j,1   ,igx,igy,igz,ig,jg,kgpp)
cc      if (k == 1  .and. wrap_z)
cc     .     call getMGmap(i,j,nz-1,igx,igy,igz,ig,jg,kgmm)

      igl = i+grid_params%ilo(igx)-1
      jgl = j+grid_params%jlo(igy)-1
      kgl = k+grid_params%klo(igz)-1

c     Grid info

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

c     Find ZIP fluxes

      !X flux
      if     (igl < grid_params%nxgl(igx)) then
        jach  = 0.5*(jac+jacip)
        flxip = 0.5*(v1(ip,j,k)*arr(i ,j,k)/jacip
     .             + v1(i ,j,k)*arr(ip,j,k)/jac  )*jach

        if (isSP(i,j,k,igx,igy,igz).and.spoint) then
          flxim = 0d0
        else
          jach  = 0.5*(jac+jacim)
          flxim = 0.5*(v1(im,j,k)*arr(i ,j,k)/jacim
     .               + v1(i ,j,k)*arr(im,j,k)/jac  )*jach
        endif
      elseif (igl == grid_params%nxgl(igx)) then
        flxip = 0.5*(v1(ip,j,k)*arr(i,j,k) + v1(i,j,k)*arr(ip,j,k))

        jach  = 0.5*(jac+jacim)
        flxim = 0.5*(v1(im,j,k)*arr(i ,j,k)/jacim
     .             + v1(i ,j,k)*arr(im,j,k)/jac  )*jach
      endif

      !Y flux
      flxjp = 0.5*(v2(i,jp,k)*arr(i,j,k) + v2(i,j,k)*arr(i,jp,k))
      flxjm = 0.5*(v2(i,jm,k)*arr(i,j,k) + v2(i,j,k)*arr(i,jm,k))

      !Z flux
      flxkp = 0.5*(v3(i,j,kp)*arr(i,j,k) + v3(i,j,k)*arr(i,j,kp))
      flxkm = 0.5*(v3(i,j,km)*arr(i,j,k) + v3(i,j,k)*arr(i,j,km))

c     Find velocities at faces from ZIP fluxes

      if (method /= 2) then

        vxp = 2.*flxip/(arr(ip,j,k)+arr(i,j,k))
        vxm = 2.*flxim/(arr(im,j,k)+arr(i,j,k))

        vyp = 2.*flxjp/(arr(i,jp,k)+arr(i,j,k))
        vym = 2.*flxjm/(arr(i,jm,k)+arr(i,j,k))

        vzp = 2.*flxkp/(arr(i,j,kp)+arr(i,j,k))
        vzm = 2.*flxkm/(arr(i,j,km)+arr(i,j,k))

c     Find advected quantity at faces

        if (i == nx) then
          qip=inter(xx(igpp),arr(ipp,j,k),xx(igp),arr(ip,j,k)
     .             ,xx(ig  ),arr(i  ,j,k),xx(igm),arr(im,j,k)
     .             ,vxp,min(method,2))
cc        if (i == nx .and. (.not.wrap_x)) then
cc          qip=inter(xx(igpp),arr(ipp,j,k),xx(igp),arr(ip,j,k)
cc       .           ,xx(ig  ),arr(i  ,j,k),xx(igm),arr(im,j,k)
cc       .           ,vxp,min(method,2))
cc        elseif (i == nx .and. wrap_x) then
cc          qip=inter(xx(igp)+dx(igpp),arr(ipp,j,k),xx(igp),arr(ip,j,k)
cc       .           ,xx(ig )         ,arr(i  ,j,k),xx(igm),arr(im,j,k)
cc       .           ,vxp,method)
        else
          qip=inter(xx(igpp),arr(ipp,j,k),xx(igp),arr(ip,j,k)
     .             ,xx(ig  ),arr(i  ,j,k),xx(igm),arr(im,j,k)
     .             ,vxp,method)
        endif

        if (i == 1) then
          qim=inter(xx(igp),arr(ip,j,k),xx(ig  ),arr(i  ,j,k)
     .             ,xx(igm),arr(im,j,k),xx(igmm),arr(imm,j,k)
     .             ,vxm,min(method,2))
cc        if (i == 1 .and. (.not.wrap_x)) then
cc          qim=inter(xx(igp),arr(ip,j,k),xx(ig  ),arr(i  ,j,k)
cc       .           ,xx(igm),arr(im,j,k),xx(igmm),arr(imm,j,k)
cc       .           ,vxm,min(method,2))
cc        elseif (i == 1 .and. wrap_x) then
cc          qim=inter(xx(igp),arr(ip,j,k),xx(ig )         ,arr(i  ,j,k)
cc       .           ,xx(igm),arr(im,j,k),xx(igm)-dx(igmm),arr(imm,j,k)
cc       .           ,vxm,method)
        else
          qim=inter(xx(igp),arr(ip,j,k),xx(ig  ),arr(i  ,j,k)
     .             ,xx(igm),arr(im,j,k),xx(igmm),arr(imm,j,k)
     .             ,vxm,method)
        endif

        if (j == ny) then
          qjp=inter(yy(jgpp),arr(i,jpp,k),yy(jgp),arr(i,jp,k)
     .             ,yy(jg  ),arr(i,j  ,k),yy(jgm),arr(i,jm,k)
     .             ,vyp,min(method,2))
cc        if (j == ny .and. (.not.wrap_y)) then
cc          qjp=inter(yy(jgpp),arr(i,jpp,k),yy(jgp),arr(i,jp,k)
cc       .           ,yy(jg  ),arr(i,j  ,k),yy(jgm),arr(i,jm,k)
cc       .           ,vyp,min(method,2))
cc        elseif (j == ny .and. wrap_y) then
cc          qjp=inter(yy(jgp)+dy(jgpp),arr(i,jpp,k),yy(jgp),arr(i,jp,k)
cc       .           ,yy(jg )         ,arr(i,j  ,k),yy(jgm),arr(i,jm,k)
cc       .           ,vyp,method)
        else
          qjp=inter(yy(jgpp),arr(i,jpp,k),yy(jgp),arr(i,jp,k)
     .             ,yy(jg  ),arr(i,j  ,k),yy(jgm),arr(i,jm,k)
     .             ,vyp,method)
        endif

        if (j == 1) then
          qjm=inter(yy(jgp),arr(i,jp,k),yy(jg  ),arr(i,j  ,k)
     .             ,yy(jgm),arr(i,jm,k),yy(jgmm),arr(i,jmm,k)
     .             ,vym,min(method,2))
cc        if (j == 1 .and. (.not.wrap_y)) then
cc          qjm=inter(yy(jgp),arr(i,jp,k),yy(jg  ),arr(i,j  ,k)
cc       .           ,yy(jgm),arr(i,jm,k),yy(jgmm),arr(i,jmm,k)
cc       .           ,vym,min(method,2))
cc        elseif (j == 1 .and. wrap_y) then
cc          qjm=inter(yy(jgp),arr(i,jp,k),yy(jg )         ,arr(i,j  ,k)
cc       .           ,yy(jgm),arr(i,jm,k),yy(jgm)-dy(jgmm),arr(i,jmm,k)
cc       .           ,vym,method)
        else
          qjm=inter(yy(jgp),arr(i,jp,k),yy(jg  ),arr(i,j  ,k)
     .             ,yy(jgm),arr(i,jm,k),yy(jgmm),arr(i,jmm,k)
     .             ,vym,method)
        endif

        if (k == nz) then
          qkp=inter(zz(kgpp),arr(i,j,kpp),zz(kgp),arr(i,j,kp)
     .             ,zz(kg  ),arr(i,j,k  ),zz(kgm),arr(i,j,km)
     .             ,vzp,min(method,2))
cc        if (k == nz .and. (.not.wrap_z)) then
cc          qkp=inter(zz(kgpp),arr(i,j,kpp),zz(kgp),arr(i,j,kp)
cc       .           ,zz(kg  ),arr(i,j,k  ),zz(kgm),arr(i,j,km)
cc       .           ,vzp,min(method,2))
cc        elseif (k == nz .and. wrap_z) then
cc          qkp=inter(zz(kgp)+dz(kgpp),arr(i,j,kpp),zz(kgp),arr(i,j,kp)
cc       .           ,zz(kg )         ,arr(i,j,k  ),zz(kgm),arr(i,j,km)
cc       .           ,vzp,method)
        else
          qkp=inter(zz(kgpp),arr(i,j,kpp),zz(kgp),arr(i,j,kp)
     .             ,zz(kg  ),arr(i,j,k  ),zz(kgm),arr(i,j,km)
     .             ,vzp,method)
        endif

        if (k == 1) then
          qkm=inter(zz(kgp),arr(i,j,kp),zz(kg  ),arr(i,j,k  )
     .             ,zz(kgm),arr(i,j,km),zz(kgmm),arr(i,j,kmm)
     .             ,vzm,min(method,2))
cc        if (k == 1 .and. (.not.wrap_z)) then
cc          qkm=inter(zz(kgp),arr(i,j,kp),zz(kg  ),arr(i,j,k  )
cc       .           ,zz(kgm),arr(i,j,km),zz(kgmm),arr(i,j,kmm)
cc       .           ,vzm,min(method,2))
cc        elseif (j == 1 .and. wrap_z) then
cc          qkm=inter(zz(kgp),arr(i,j,kp),zz(kg )         ,arr(i,j,k  )
cc       .           ,zz(kgm),arr(i,j,km),zz(kgm)-dz(kgmm),arr(i,j,kmm)
cc       .           ,vzm,method)       
        else
          qkm=inter(zz(kgp),arr(i,j,kp),zz(kg  ),arr(i,j,k  )
     .             ,zz(kgm),arr(i,j,km),zz(kgmm),arr(i,j,kmm)
     .             ,vzm,method)
        endif

c     Setup fluxes (ZIP at boundaries for conservation)

        if (igl<grid_params%nxgl(igx).or.bcond(2)==PER) flxip=vxp*qip
        if (igl>1                    .or.bcond(1)==PER) flxim=vxm*qim

        if (jgl<grid_params%nygl(igy).or.bcond(4)==PER) flxjp=vyp*qjp
        if (jgl>1                    .or.bcond(3)==PER) flxjm=vym*qjm

        if (kgl<grid_params%nzgl(igz).or.bcond(6)==PER) flxkp=vzp*qkp
        if (kgl>1                    .or.bcond(5)==PER) flxkm=vzm*qkm

      endif

c     Find advective term

      advec =( (flxip - flxim)/dxx
     .       + (flxjp - flxjm)/dyy
     .       + (flxkp - flxkm)/dzz )/jac
      
      if (vw) advec = advec*gmetric%grid(igx)%dvol(i,j,k)

c     End 

      end function flx_advec

c     inter
c     ###############################################################
      function inter(x1,q1,x2,q2,x3,q3,x4,q4,vel,advect)

c     ---------------------------------------------------------------
c     This function computes the advection interpolation at control
c     volume face. Options:
c       * advect = 1 or 7 => First-order upwind
c       * advect = 2 => Centered
c       * advect = 3 => High-order upwind (QUICK)
c       * advect = 4 => Monotone high-order upwind (SMART)
c       * advect = 5 => SMART with smooth transition
c       * advect = 6 => Centered, high-order
c       * advect = 8 => Gamma (modification of SMART)
c
c     Convention:
c
c       ---x-----x--o--x-----x---->
c          q4    q3 ^  q2    q1
c                   ^
c           Location of face
c     ---------------------------------------------------------------

      use oned_int

      use math

      implicit none                !For safe fortran

c     Call variables

      real(8)    :: x1,x2,x3,x4,q1,q2,q3,q4,vel,inter
      integer    :: advect

c     Local variables

      real(8)    :: dx4,dx3,dx2,qt1,qt2,qt3,qt4,slp1,slp2,a,b,c,xp1,xp2
     .             ,curv

c     Begin program

      select case(advect)
      case(1,7)

c     Upwind

        if (vel.gt.0d0) then
          inter = q3
        else
          inter = q2
        endif

      case(2)

c     Centered

        inter = 0.5*(q2+q3)

      case(3)

c     QUICK [Leonard, Comput. Meth. Appl. Mech. Eng. 19, p. 59 (1979)]

        if (vel.gt.0d0) then
          curv = ((q2-q3)/(x2-x3)-(q3-q4)/(x3-x4))/0.5/(x2-x4)
        else
          curv = ((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/0.5/(x1-x3)
        endif

        inter = 0.5*(q2+q3) - 0.125*curv*(x2-x3)**2

cc        if (vel.gt.0d0) then
cc          inter = 3.*q2/8. + 3.*q3/4.- q4/8.
cc        else
cc          inter = 3.*q3/8. + 3.*q2/4.- q1/8.
cc        endif

      case(4)

c     SMART [Gaskell and Lau, Int. J. Num. Meth. Fluids 8, p 617-641 (1988)]

        slp1 = 1.5
        slp2 = .5
        if (vel.gt.0d0) then
          curv = ((q2-q3)/(x2-x3)-(q3-q4)/(x3-x4))/0.5/(x2-x4)
          qt1 = 0.5*(q2+q3) - 0.125*curv*(x2-x3)**2
cc          qt1 = 3.*q2/8. + 3.*q3/4.- q4/8.
          qt2 = slp1*q3 + (1.-slp1)*q4
          qt3 = slp2*q3 + (1.-slp2)*q2
          qt4 = fmed(q3,qt2,qt3)
          inter = fmed(q3,qt4,qt1)
        else
          curv = ((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/0.5/(x1-x3)
          qt1 = 0.5*(q2+q3) - 0.125*curv*(x2-x3)**2
cc          qt1 = 3.*q3/8. + 3.*q2/4.- q1/8.
          qt2 = slp1*q2 + (1.-slp1)*q1
          qt3 = slp2*q2 + (1.-slp2)*q3
          qt4 = fmed(q2,qt2,qt3)
          inter = fmed(q2,qt4,qt1)
        endif

      case(5)

c     Smooth SMART

        write (*,*) 'Need to fix this for non-uniform grids'
        stop

        xp1 = .3
        xp2 = .5
        if (vel.gt.0d0) then
          qt1  = (q3-q4)/(q2-q4)
          if (qt1.lt.xp1.and.qt1.ge.0.) then
            a = .25/xp1**3*(xp1-3.)
            b = 1./xp1**2*(9./8.-xp1/2.)
            qt2 = qt1*(a*qt1**2 + b*qt1 + 1.)
            inter = q4 + qt2*(q2-q4)
          elseif (qt1.ge.xp2.and.qt1.lt.1.) then
            a = (xp2/4.-.5)/(xp2-1.)**3
            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
            inter = q4 + qt2*(q2-q4)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
            inter = 3.*q2/8. + 3.*q3/4.- q4/8.
          else
            inter = q3
          endif
        else
          qt1  = (q2-q1)/(q3-q1)
          if (qt1.lt.xp1.and.qt1.ge.0.) then
            a = .25/xp1**3*(xp1-3.)
            b = 1./xp1**2*(9./8.-xp1/2.)
            c = 1.
            qt2 = qt1*(a*qt1**2 + b*qt1 + c)
            inter = q1 + qt2*(q3-q1)
          elseif (qt1.ge.xp2.and.qt1.lt.1.) then
            a = (xp2/4.-.5)/(xp2-1.)**3
            b = 1./8./(xp2-1.)**3*(-4.*xp2**2+7.*xp2+1)
            c = 1. - (-.5*xp2**2-3./8+9./8*xp2)/(xp2-1.)**3
            qt2 = (qt1 - 1.)*(a*qt1**2 + b*qt1 + c) + 1.
            inter = q1 + qt2*(q3-q1)
          elseif (qt1.ge.xp1.and.qt1.lt.xp2) then
            inter = 3.*q3/8. + 3.*q2/4.- q1/8.
          else
            inter = q2
          endif
        endif

      case(6)

c     Centered high-order

        inter = (-q4 + 9.*q3 + 9.*q2 - q1)/16.

cc      case(8)
cc
ccc     Gamma [Jasak, Weller, and Gosman, Int. J. Numer. Meth. Fluids, 32, p431 (1999)]
cc
cc        slp1 = 1.5
cc        slp2 = .5
cc        if (vel.gt.0d0) then
cc          curv = ((q2-q3)/(x2-x3)-(q3-q4)/(x3-x4))/0.5/(x2-x4)
cc          qt1 = 0.5*(q2+q3) - 0.125*curv*(x2-x3)**2
cccc          qt1 = 3.*q2/8. + 3.*q3/4.- q4/8.
cc          qt2 = slp1*q3 + (1.-slp1)*q4
cc          qt3 = slp2*q3 + (1.-slp2)*q2
cc          qt4 = fmed(q3,qt2,qt3)
cc          inter = fmed(q3,qt4,qt1)
cc        else
cc          curv = ((q1-q2)/(x1-x2)-(q2-q3)/(x2-x3))/0.5/(x1-x3)
cc          qt1 = 0.5*(q2+q3) - 0.125*curv*(x2-x3)**2
cccc          qt1 = 3.*q3/8. + 3.*q2/4.- q1/8.
cc          qt2 = slp1*q2 + (1.-slp1)*q1
cc          qt3 = slp2*q2 + (1.-slp2)*q3
cc          qt4 = fmed(q2,qt2,qt3)
cc          inter = fmed(q2,qt4,qt1)
cc        endif

      case default

        call pstop('inter','Advection scheme not available')

      end select

c     End

      end function inter

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

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian

      real(8)    :: dff(0:nx+1,0:ny+1,0:nz+1)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

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

      d_xy_ipjp = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                 +dff(ip,j ,k)*gsup(ip,j ,k,1,2)
     .                 +dff(i ,jp,k)*gsup(i ,jp,k,1,2)
     .                 +dff(ip,jp,k)*gsup(ip,jp,k,1,2))
      d_xy_ipjm = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                 +dff(ip,j ,k)*gsup(ip,j ,k,1,2)
     .                 +dff(i ,jm,k)*gsup(i ,jm,k,1,2)
     .                 +dff(ip,jm,k)*gsup(ip,jm,k,1,2))
      d_xy_imjp = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                 +dff(im,j ,k)*gsup(im,j ,k,1,2)
     .                 +dff(i ,jp,k)*gsup(i ,jp,k,1,2)
     .                 +dff(im,jp,k)*gsup(im,jp,k,1,2))
      d_xy_imjm = 0.25*(dff(i ,j ,k)*gsup(i ,j ,k,1,2)
     .                 +dff(im,j ,k)*gsup(im,j ,k,1,2)
     .                 +dff(i ,jm,k)*gsup(i ,jm,k,1,2)
     .                 +dff(im,jm,k)*gsup(im,jm,k,1,2))
                                                      
      d_xz_ipkp = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                 +dff(ip,j,k )*gsup(ip,j,k ,1,3)
     .                 +dff(i ,j,kp)*gsup(i ,j,kp,1,3)
     .                 +dff(ip,j,kp)*gsup(ip,j,kp,1,3))
      d_xz_ipkm = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                 +dff(ip,j,k )*gsup(ip,j,k ,1,3)
     .                 +dff(i ,j,km)*gsup(i ,j,km,1,3)
     .                 +dff(ip,j,km)*gsup(ip,j,km,1,3))
      d_xz_imkp = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                 +dff(im,j,k )*gsup(im,j,k ,1,3)
     .                 +dff(i ,j,kp)*gsup(i ,j,kp,1,3)
     .                 +dff(im,j,kp)*gsup(im,j,kp,1,3))
      d_xz_imkm = 0.25*(dff(i ,j,k )*gsup(i ,j,k ,1,3)
     .                 +dff(im,j,k )*gsup(im,j,k ,1,3)
     .                 +dff(i ,j,km)*gsup(i ,j,km,1,3)
     .                 +dff(im,j,km)*gsup(im,j,km,1,3))
                                                      
      d_yz_jpkp = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                 +dff(i,jp,k )*gsup(i,jp,k ,2,3)
     .                 +dff(i,j ,kp)*gsup(i,j ,kp,2,3)
     .                 +dff(i,jp,kp)*gsup(i,jp,kp,2,3))
      d_yz_jpkm = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                 +dff(i,jp,k )*gsup(i,jp,k ,2,3)
     .                 +dff(i,j ,km)*gsup(i,j ,km,2,3)
     .                 +dff(i,jp,km)*gsup(i,jp,km,2,3))
      d_yz_jmkp = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                 +dff(i,jm,k )*gsup(i,jm,k ,2,3)
     .                 +dff(i,j ,kp)*gsup(i,j ,kp,2,3)
     .                 +dff(i,jm,kp)*gsup(i,jm,kp,2,3))
      d_yz_jmkm = 0.25*(dff(i,j ,k )*gsup(i,j ,k ,2,3)
     .                 +dff(i,jm,k )*gsup(i,jm,k ,2,3)
     .                 +dff(i,j ,km)*gsup(i,j ,km,2,3)
     .                 +dff(i,jm,km)*gsup(i,jm,km,2,3))

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

      if (.not.vol_wgt) laplacian=laplacian
     .                           /gmetric%grid(igx)%dvol(i,j,k)

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

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),laplacian

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid

      real(8)    :: d_xx_ip,d_xx_im,d_yy_jp,d_yy_jm,d_zz_kp,d_zz_km
     .             ,d_xy_ipjp,d_xy_ipjm,d_xy_imjp,d_xy_imjm
     .             ,d_xz_ipkp,d_xz_ipkm,d_xz_imkp,d_xz_imkm
     .             ,d_yz_jpkp,d_yz_jpkm,d_yz_jmkp,d_yz_jmkm

      real(8),dimension(:,:,:,:,:),pointer :: gsup

      logical    :: vol_wgt

c     Begin program

      if (PRESENT(vol)) then
        vol_wgt = vol
      else
        vol_wgt = .true.
      endif

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

      if (coords == 'car') then

        laplacian = 
     .     dyh(jg)*dzh(kg)*( (arr(ip,j,k)-arr(i,j,k))/dx(ig)
     .                      +(arr(im,j,k)-arr(i,j,k))/dx(ig-1) )
     .    +dxh(ig)*dzh(kg)*( (arr(i,jp,k)-arr(i,j,k))/dy(jg)
     .                      +(arr(i,jm,k)-arr(i,j,k))/dy(jg-1) )
     .    +dxh(ig)*dyh(jg)*( (arr(i,j,kp)-arr(i,j,k))/dz(kg)
     .                      +(arr(i,j,km)-arr(i,j,k))/dz(kg-1) )
     .    +0.5*( dzh(kg)*( (arr(ip,jp,k)-arr(i,j,k))
     .                    +(arr(im,jm,k)-arr(i,j,k))
     .                    -(arr(ip,jm,k)-arr(i,j,k))
     .                    -(arr(im,jp,k)-arr(i,j,k)) )
     .          +dyh(jg)*( (arr(ip,j,kp)-arr(i,j,k))
     .                    +(arr(im,j,km)-arr(i,j,k))
     .                    -(arr(ip,j,km)-arr(i,j,k))
     .                    -(arr(im,j,kp)-arr(i,j,k)) )
     .          +dxh(ig)*( (arr(i,jp,kp)-arr(i,j,k))
     .                    +(arr(i,jm,km)-arr(i,j,k))
     .                    -(arr(i,jp,km)-arr(i,j,k))
     .                    -(arr(i,jm,kp)-arr(i,j,k)) ) )

      else

        gsup => gmetric%grid(igrid)%gsup

cc        if (isSP(i,j,k,igx,igy,igz)) then
cc          gsup(im,jm:jp,km:kp,:,:) = gsup(i,jm:jp,km:kp,:,:)
cc        endif

        d_xx_ip = 0.5*(gsup(i ,j,k,1,1)
     .                +gsup(ip,j,k,1,1))
        d_xx_im = 0.5*(gsup(i ,j,k,1,1)
     .                +gsup(im,j,k,1,1))

        d_yy_jp = 0.5*(gsup(i,j ,k,2,2)
     .                +gsup(i,jp,k,2,2))
        d_yy_jm = 0.5*(gsup(i,j ,k,2,2)
     .                +gsup(i,jm,k,2,2))

        d_zz_kp = 0.5*(gsup(i,j,k ,3,3)
     .                +gsup(i,j,kp,3,3))
        d_zz_km = 0.5*(gsup(i,j,k ,3,3)
     .                +gsup(i,j,km,3,3))

        d_xy_ipjp = 0.25*(gsup(i ,j ,k,1,2)
     .                   +gsup(ip,j ,k,1,2)
     .                   +gsup(i ,jp,k,1,2)
     .                   +gsup(ip,jp,k,1,2))
        d_xy_ipjm = 0.25*(gsup(i ,j ,k,1,2)
     .                   +gsup(ip,j ,k,1,2)
     .                   +gsup(i ,jm,k,1,2)
     .                   +gsup(ip,jm,k,1,2))
        d_xy_imjp = 0.25*(gsup(i ,j ,k,1,2)
     .                   +gsup(im,j ,k,1,2)
     .                   +gsup(i ,jp,k,1,2)
     .                   +gsup(im,jp,k,1,2))
        d_xy_imjm = 0.25*(gsup(i ,j ,k,1,2)
     .                   +gsup(im,j ,k,1,2)
     .                   +gsup(i ,jm,k,1,2)
     .                   +gsup(im,jm,k,1,2))

        d_xz_ipkp = 0.25*(gsup(i ,j,k ,1,3)
     .                   +gsup(ip,j,k ,1,3)
     .                   +gsup(i ,j,kp,1,3)
     .                   +gsup(ip,j,kp,1,3))
        d_xz_ipkm = 0.25*(gsup(i ,j,k ,1,3)
     .                   +gsup(ip,j,k ,1,3)
     .                   +gsup(i ,j,km,1,3)
     .                   +gsup(ip,j,km,1,3))
        d_xz_imkp = 0.25*(gsup(i ,j,k ,1,3)
     .                   +gsup(im,j,k ,1,3)
     .                   +gsup(i ,j,kp,1,3)
     .                   +gsup(im,j,kp,1,3))
        d_xz_imkm = 0.25*(gsup(i ,j,k ,1,3)
     .                   +gsup(im,j,k ,1,3)
     .                   +gsup(i ,j,km,1,3)
     .                   +gsup(im,j,km,1,3))

        d_yz_jpkp = 0.25*(gsup(i,j ,k ,2,3)
     .                   +gsup(i,jp,k ,2,3)
     .                   +gsup(i,j ,kp,2,3)
     .                   +gsup(i,jp,kp,2,3))
        d_yz_jpkm = 0.25*(gsup(i,j ,k ,2,3)
     .                   +gsup(i,jp,k ,2,3)
     .                   +gsup(i,j ,km,2,3)
     .                   +gsup(i,jp,km,2,3))
        d_yz_jmkp = 0.25*(gsup(i,j ,k ,2,3)
     .                   +gsup(i,jm,k ,2,3)
     .                   +gsup(i,j ,kp,2,3)
     .                   +gsup(i,jm,kp,2,3))
        d_yz_jmkm = 0.25*(gsup(i,j ,k ,2,3)
     .                   +gsup(i,jm,k ,2,3)
     .                   +gsup(i,j ,km,2,3)
     .                   +gsup(i,jm,km,2,3))

cc        if (isSP(i,j,k,igx,igy,igz)) then
cc          d_xx_im   = 0d0
cc          d_xy_imjp = 0d0
cc          d_xy_imjm = 0d0
cc          d_xz_imkp = 0d0
cc          d_xz_imkm = 0d0
cc        endif

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

      endif


      if (.not.vol_wgt) laplacian=laplacian
     .                           /gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function lap_ndiff

c     div_tensor
c     ###############################################################
      function div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                   ,tsrx,tsry,tsrz,vol) result (divt)

c     ---------------------------------------------------------------
c     Calculatesthe contravariant components of div(tensor) at cell
c     centers in general non-orthogonal coordinates, given by:
c
c         div(T)^i = partial_j(T^ji/jac) + T^lm/jac Gamma(i)_lm)
c
c     The routines tsrx,tsry,tsrz must provide the contravariant
c     tersor components T^lm (if flag /= 0, they actually provide
c     T^lm/jac).
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8)    :: divt(3)

      integer    :: i,j,k,igx,igy,igz,nx,ny,nz

      logical    :: alt_eom

      external   tsrx,tsry,tsrz

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid

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
      if (vol_wgt) divt=divt*gmetric%grid(igrid)%dvol(i,j,k)

c     End program

      end function div_tensor

c     div_tensor_cov
c     ###############################################################
      function div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                       ,tsrx,tsry,tsrz,vol) result (divt)

c     ---------------------------------------------------------------
c     Calculates covariant components of div(tensor) at cell centers in
c     general non-orthogonal coordinates, which is given by:
c
c         div(T)_i = jac^(-1)*[ partial_j(T^j_i)-T^l_m Gamma(m)_il) ]
c
c     The routinex tsrx,tsry,tsrz must provide mixed representation 
c     tensor components T^l_m.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8)    :: divt(3)

      integer    :: i,j,k,igx,igy,igz,nx,ny,nz

      external   tsrx,tsry,tsrz

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid

      real(8)    :: jac,dvol,dS1,dS2,dS3,dxx,dyy,dzz

cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

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

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      dvol = dxx*dyy*dzz

      if (coords /= 'car') then
        hess = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
      endif

      call tsrx(i ,j,k,nx,ny,nz,igx,igy,igz,t11p,t12p,t13p, 1)
      call tsrx(im,j,k,nx,ny,nz,igx,igy,igz,t11m,t12m,t13m,-1)
      if (coords /= 'car')
     .    call tsrx(i,j,k,nx,ny,nz,igx,igy,igz,t11o,t12o,t13o,0)

      call tsry(i,j ,k,nx,ny,nz,igx,igy,igz,t21p,t22p,t23p, 1)
      call tsry(i,jm,k,nx,ny,nz,igx,igy,igz,t21m,t22m,t23m,-1)
      if (coords /= 'car')
     .    call tsry(i,j,k,nx,ny,nz,igx,igy,igz,t21o,t22o,t23o,0)

      call tsrz(i,j,k ,nx,ny,nz,igx,igy,igz,t31p,t32p,t33p, 1)
      call tsrz(i,j,km,nx,ny,nz,igx,igy,igz,t31m,t32m,t33m,-1)
      if (coords /= 'car')
     .    call tsrz(i,j,k,nx,ny,nz,igx,igy,igz,t31o,t32o,t33o,0)

      msource = 0d0

      !Component 1
cc      flxip = t11p
cc      flxim = t11m
ccc
cc      flxjp = t21p
cc      flxjm = t21m
cc
cc      flxkp = t31p
cc      flxkm = t31m

      if (coords /= 'car') then
        msource =  (t11o*hess(1,1,1)+t12o*hess(2,1,1)+t13o*hess(3,1,1)
     .             +t21o*hess(1,2,1)+t22o*hess(2,2,1)+t23o*hess(3,2,1)
     .             +t31o*hess(1,3,1)+t32o*hess(2,3,1)+t33o*hess(3,3,1))
      endif

cc      divt(1) =  (flxip - flxim)/dxx
cc     .          +(flxjp - flxjm)/dyy
cc     .          +(flxkp - flxkm)/dzz - msource
      divt(1) =  (t11p - t11m)/dxx
     .          +(t21p - t21m)/dyy
     .          +(t31p - t31m)/dzz - msource

      !Component 2
cc      flxip = t12p
cc      flxim = t12m
cc
cc      flxjp = t22p
cc      flxjm = t22m
cc
cc      flxkp = t32p
cc      flxkm = t32m

      if (coords /= 'car') then
        msource=  (t11o*hess(1,1,2)+t12o*hess(2,1,2)+t13o*hess(3,1,2)
     .            +t21o*hess(1,2,2)+t22o*hess(2,2,2)+t23o*hess(3,2,2)
     .            +t31o*hess(1,3,2)+t32o*hess(2,3,2)+t33o*hess(3,3,2))
      endif

cc      divt(2) =  (flxip - flxim)/dxx
cc     .          +(flxjp - flxjm)/dyy
cc     .          +(flxkp - flxkm)/dzz  - msource
      divt(2) =  (t12p - t12m)/dxx
     .          +(t22p - t22m)/dyy
     .          +(t32p - t32m)/dzz  - msource

      !Component 3
cc      flxip = t13p
cc      flxim = t13m
cc
cc      flxjp = t23p
cc      flxjm = t23m
cc
cc      flxkp = t33p
cc      flxkm = t33m

      if (coords /= 'car') then
        msource = (t11o*hess(1,1,3)+t12o*hess(2,1,3)+t13o*hess(3,1,3)
     .            +t21o*hess(1,2,3)+t22o*hess(2,2,3)+t23o*hess(3,2,3)
     .            +t31o*hess(1,3,3)+t32o*hess(2,3,3)+t33o*hess(3,3,3))
      endif

cc      divt(3) =  (flxip - flxim)/dxx
cc     .          +(flxjp - flxjm)/dyy
cc     .          +(flxkp - flxkm)/dzz - msource
      divt(3) =  (t13p - t13m)/dxx
     .          +(t23p - t23m)/dyy
     .          +(t33p - t33m)/dzz - msource

      divt = divt/jac

      !Volume factor
      if (vol_wgt) divt=divt*gmetric%grid(igrid)%dvol(i,j,k)

c     End program

      end function div_tensor_cov

c     veclap2
c     ###############################################################
      function veclap2(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,vcov,vol)
     .         result (vlap)

c     ---------------------------------------------------------------
c     Calculates dvol*lap(vector) at cell centers in general non-orthog.
c     coordinates using
c                  lap(vector)=grad(div(v)-curl(curl(v)). 
c     Returns lap(vector) in contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vcnv(0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,vcov(0:nx+1,0:ny+1,0:nz+1,3)

      real(8)    :: vlap(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: icomp

      real(8)    :: cov(3)

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

      if (coords == 'car') then
        do icomp=1,3
          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,vcnv(:,:,:,icomp)
     .                           ,vol=vol_wgt)
        enddo
      else
        cov = graddiv(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,vol=vol_wgt)
        vlap = matmul(gmetric%grid(igx)%gsup(i,j,k,:,:),cov) !Xforms to cnv

        vlap = vlap - curlcurl(i,j,k,nx,ny,nz,igx,igy,igz
     .                        ,vcov,vol=vol_wgt)
      endif

c     End program

      end function veclap2

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

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,diff   (0:nx+1,0:ny+1,0:nz+1)

      real(8)    :: vlap(3)

      logical    :: alteom

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: icomp

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
      else
        vec  => vfield !Pointer passed to nabtensor routines
        coef => diff   !Pointer passed to nabtensor routines
        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
     $                   ,vol=vol_wgt)
        nullify(coef)
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

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)

      real(8)    :: vlap(3)

      logical    :: alteom

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: icomp

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
      else
        vec  => vfield !Pointer passed to nabtensor routines
        nullify(coef)
        vlap = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                   ,nabtensor_x,nabtensor_y,nabtensor_z
     $                   ,vol=vol_wgt)
      endif

c     End program

      end function veclap_ndiff

c     nabtensor_x
c     #############################################################
      subroutine nabtensor_x(i,j,k,nx,ny,nz,igx,igy,igz,alteom
     .                      ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for nabla(vec)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alteom

c     Local variables

        integer    :: ig,jg,kg,ip,igrid
        real(8)    :: x,y,z,jac,jac0,jacp,ijac
        real(8)    :: nabla_v_cnv(3,3),gsuper(3,3),dd

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
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),1,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0,cnv=.true.)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(ip,j,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t11 = dd*nabla_v_cnv(1,1)
        t12 = dd*nabla_v_cnv(1,2)
        t13 = dd*nabla_v_cnv(1,3)

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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alteom

c     Local variables

        integer    :: ig,jg,kg,jp,igrid
        real(8)    :: x,y,z,jac,ijac
        real(8)    :: nabla_v_cnv(3,3),gsuper(3,3),dd

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
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),2,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0,cnv=.true.)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,jp,k) + 1./coef(i,j ,k))
        else
          dd = 1d0
        endif

        t21 = dd*nabla_v_cnv(2,1)
        t22 = dd*nabla_v_cnv(2,2)
        t23 = dd*nabla_v_cnv(2,3)

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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alteom

c     Local variables

        integer    :: ig,jg,kg,kp,igrid
        real(8)    :: x,y,z,jac,ijac
        real(8)    :: nabla_v_cnv(3,3),gsuper(3,3),dd

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
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),3,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0,cnv=.true.)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,j,kp) + 1./coef(i,j,k ))
        else
          dd = 1d0
        endif

        t31 = dd*nabla_v_cnv(3,1)
        t32 = dd*nabla_v_cnv(3,2)
        t33 = dd*nabla_v_cnv(3,3)

        if (flag /= 0) then
          t31 = t31*ijac
          if (.not.alteom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine nabtensor_z

c     veclap_cov_ndiff
c     ###############################################################
      function veclap_cov_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,vfield
     .                         ,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates covariant components of dvol*lap(vfield) at cell
c     centers in general non-orthogonal coordinates. Vector vfield
c     in covariant representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)

      real(8)    :: vlap(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: icomp

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

      if (coords == 'car') then
        do icomp=1,3
          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,vfield(:,:,:,icomp),vol=vol_wgt)
        enddo
      else
        vec => vfield !Pointer passed to nabtensor routines
        nullify(coef)
        vlap = div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .             ,nabtensor_cov_x,nabtensor_cov_y,nabtensor_cov_z
     $             ,vol=vol_wgt)
      endif

c     End program

      end function veclap_cov_ndiff

c     veclap_cov_diff
c     ###############################################################
      function veclap_cov_diff(i,j,k,nx,ny,nz,igx,igy,igz,vfield,diff
     .                         ,vol) result (vlap)

c     ---------------------------------------------------------------
c     Calculates covariant components of dvol*lap(vfield) at cell
c     centers in general non-orthogonal coordinates. Vector vfield
c     in covariant representation.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8),target :: vfield (0:nx+1,0:ny+1,0:nz+1,3)
     .                 ,diff   (0:nx+1,0:ny+1,0:nz+1)

      real(8)    :: vlap(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: icomp

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

      if (coords == 'car') then
        do icomp=1,3
          vlap(icomp)=laplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                         ,vfield(:,:,:,icomp),diff
     .                         ,vol=vol_wgt)
        enddo
      else
        vec  => vfield !Pointer passed to nabtensor routines
        coef => diff   !Pointer passed to nabtensor routines
        vlap = div_tensor_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .             ,nabtensor_cov_x,nabtensor_cov_y,nabtensor_cov_z
     $             ,vol=vol_wgt)
        nullify(coef)
      endif

c     End program

      end function veclap_cov_diff

c     nabtensor_cov_x
c     #############################################################
      subroutine nabtensor_cov_x(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t11-t13 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13

c     Local variables

        integer    :: ig,jg,kg,ip,igrid
        real(8)    :: x,y,z,jac,jac0,jacp
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

        if (flag /= 0) then
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),1)
        else
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(ip,j,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t11 = dd*( gsuper(1,1)*nabla_v(1,1)
     .            +gsuper(1,2)*nabla_v(2,1)
     .            +gsuper(1,3)*nabla_v(3,1) )

        t12 = dd*( gsuper(1,1)*nabla_v(1,2)
     .            +gsuper(1,2)*nabla_v(2,2)
     .            +gsuper(1,3)*nabla_v(3,2) )

        t13 = dd*( gsuper(1,1)*nabla_v(1,3)
     .            +gsuper(1,2)*nabla_v(2,3)
     .            +gsuper(1,3)*nabla_v(3,3) )

c     End program

      end subroutine nabtensor_cov_x

c     nabtensor_cov_y
c     #############################################################
      subroutine nabtensor_cov_y(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t21-t23 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23

c     Local variables

        integer    :: ig,jg,kg,jp,igrid
        real(8)    :: x,y,z,jac
        real(8)    :: nabla_v(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
     .               +gmetric%grid(igrid)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),2)
        else
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,jp,k) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t21 = dd*( gsuper(2,1)*nabla_v(1,1)
     .            +gsuper(2,2)*nabla_v(2,1)
     .            +gsuper(2,3)*nabla_v(3,1) )

        t22 = dd*( gsuper(2,1)*nabla_v(1,2)
     .            +gsuper(2,2)*nabla_v(2,2)
     .            +gsuper(2,3)*nabla_v(3,2) )

        t23 = dd*( gsuper(2,1)*nabla_v(1,3)
     .            +gsuper(2,2)*nabla_v(2,3)
     .            +gsuper(2,3)*nabla_v(3,3) )

c     End program

      end subroutine nabtensor_cov_y

c     nabtensor_cov_z
c     #############################################################
      subroutine nabtensor_cov_z(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates mixed representation tensor components t31-t33 for
c     nabla(vec), [nab(vec))^l_m.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33

c     Local variables

        integer    :: ig,jg,kg,kp,igrid
        real(8)    :: x,y,z,jac
        real(8)    :: nabla_v(3,3),gsuper(3,3),dd

c     Begin program

        igrid = igx

        kp=k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
     .               +gmetric%grid(igrid)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),3)
        else
          nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                          ,vec(:,:,:,1)
     .                          ,vec(:,:,:,2)
     .                          ,vec(:,:,:,3),0)
        endif

        if (associated(coef)) then
          dd = 2./(1./coef(i,j,kp) + 1./coef(i,j,k))
        else
          dd = 1d0
        endif

        t31 = dd*( gsuper(3,1)*nabla_v(1,1)
     .            +gsuper(3,2)*nabla_v(2,1)
     .            +gsuper(3,3)*nabla_v(3,1) )

        t32 = dd*( gsuper(3,1)*nabla_v(1,2)
     .            +gsuper(3,2)*nabla_v(2,2)
     .            +gsuper(3,3)*nabla_v(3,2) )

        t33 = dd*( gsuper(3,1)*nabla_v(1,3)
     .            +gsuper(3,2)*nabla_v(2,3)
     .            +gsuper(3,3)*nabla_v(3,3) )

c     End program

      end subroutine nabtensor_cov_z

c     grad
c     ###############################################################
      subroutine grad(i,j,k,nx,ny,nz,igx,igy,igz,arr,gx,gy,gz)

c     ---------------------------------------------------------------
c     Calculates grad(A)) in general non-orthogonal coordinates,
c     preserving the SPD property. The vector grad(A) is covariant.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: arr(0:nx+1,0:ny+1,0:nz+1),gx,gy,gz

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid
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

      integer    :: i,j,k,comp,nx,ny,nz,igx,igy,igz
      integer   , optional :: he

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid,half_elem
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
c     Calculates curl(A) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns the contravariant component "comp".
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      real(8)    :: curl2

      integer    :: i,j,k,comp,nx,ny,nz,igx,igy,igz

      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
     .             ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid
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
c     Finds contravariant components (a1,a2,a3) of curl(bb x vv) at the
c     grid node (i,j,k). One sided derivatives are employed when half_elem=1
c     (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nx,ny,nz,half_elem,igx,igy,igz
        real(8)    :: cnv(3)
        real(8)    :: vv(0:nx+1,0:ny+1,0:nz+1,3)
     $               ,bb(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq
        integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

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

        if (coords == 'car') then
          ijacip  = 1d0
          ijacim  = 1d0
          ijacjp  = 1d0
          ijacjm  = 1d0
          ijackp  = 1d0
          ijackm  = 1d0
        endif
cc        sing_point = isSP(i,j,k,igx,igy,igz)

c     Exceptions

        select case(half_elem)
        case (1)
          idhx = 1./dx(ig)
          im = i

          if (coords /= 'car') then
            jacip  = gmetric%grid(igx)%jac(ip,j,k)
            jacim  = gmetric%grid(igx)%jac(i ,j,k)
            jacjp  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
     .                   +gmetric%grid(igx)%jac(i ,jp,k))
            jacjm  = 0.5*(gmetric%grid(igx)%jac(ip,jm,k)
     .                   +gmetric%grid(igx)%jac(i ,jm,k))
            jackp  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
     .                   +gmetric%grid(igx)%jac(i ,j,kp))
            jackm  = 0.5*(gmetric%grid(igx)%jac(ip,j,km)
     .                   +gmetric%grid(igx)%jac(i ,j,km))

            if (isSP(i+1,j,k,igx,igy,igz)) then
              jacjp = SP_flsv
              jacjm = SP_flsv
              jackp = SP_flsv
              jackm = SP_flsv
cc              write (*,*) 'DIAG -- curl_bxv'
cc              jacjp = 1d0!SP_flsv
cc              jacjm = 1d0!SP_flsv
cc              jackp = 1d0!SP_flsv
cc              jackm = 1d0!SP_flsv
            endif
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

          if (coords /= 'car') then
            jacip  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
     .                   +gmetric%grid(igx)%jac(ip,j ,k))
            jacim  = 0.5*(gmetric%grid(igx)%jac(im,jp,k)
     .                   +gmetric%grid(igx)%jac(im,j ,k))
            jacjp  = gmetric%grid(igx)%jac(i,jp,k)
            jacjm  = gmetric%grid(igx)%jac(i,j ,k)
            jackp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
     .                   +gmetric%grid(igx)%jac(i,j ,kp))
            jackm  = 0.5*(gmetric%grid(igx)%jac(i,jp,km)
     .                   +gmetric%grid(igx)%jac(i,j ,km))
          endif

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

          if (coords /= 'car') then
            jacip  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
     .                   +gmetric%grid(igx)%jac(ip,j,k ))
            jacim  = 0.5*(gmetric%grid(igx)%jac(im,j,kp)
     .                   +gmetric%grid(igx)%jac(im,j,k ))
            jacjp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
     .                   +gmetric%grid(igx)%jac(i,jp,k ))
            jacjm  = 0.5*(gmetric%grid(igx)%jac(i,jm,kp)
     .                   +gmetric%grid(igx)%jac(i,jm,k ))
            jackp  = gmetric%grid(igx)%jac(i,j,kp)
            jackm  = gmetric%grid(igx)%jac(i,j,km)
          endif

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

          if (coords /= 'car') then
            jacip  = gmetric%grid(igx)%jac(ip,j,k)
            jacim  = gmetric%grid(igx)%jac(im,j,k)
            jacjp  = gmetric%grid(igx)%jac(i,jp,k)
            jacjm  = gmetric%grid(igx)%jac(i,jm,k)
            jackp  = gmetric%grid(igx)%jac(i,j,kp)
            jackm  = gmetric%grid(igx)%jac(i,j,km)
          endif

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

        if (coords /= 'car') then
          ijacip  = 1d0/jacip
          ijacim  = 1d0/jacim
          ijacjp  = 1d0/jacjp
          ijacjm  = 1d0/jacjm
          ijackp  = 1d0/jackp
          ijackm  = 1d0/jackm
        endif

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

ccc     curl_bxv_vrtx
ccc     ###################################################################
cc      function curl_bxv_vrtx(i,j,k,nx,ny,nz,igx,igy,igz,vv,bb)
cc     .         result(cnv)
cc
ccc     -------------------------------------------------------------------
ccc     Finds contravariant components (a1,a2,a3) of curl(bb x vv) at the
ccc     grid node (i,j,k). One sided derivatives are employed when hex=1
ccc     (i+1/1), hey=1 (j+1/2), and hez=1 (k+1/2).
ccc     -------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,nx,ny,nz,hex,hey,hez,igx,igy,igz
cc        real(8)    :: cnv(3)
cc        real(8)    :: vv(0:nx+1,0:ny+1,0:nz+1,3)
cc     $               ,bb(0:nx+1,0:ny+1,0:nz+1,3)
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq
cc        integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
cc
cc        real(8)    :: idhx,idhy,idhz,a(3)
cc        real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc
cc        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
cc     .               ,jacp,jacm,jacph,jacmh,jach,jac0
cc     .               ,ijacip,ijacim,ijacjp,ijacjm,ijackp,ijackm
cc     .               ,ijacp,ijacm,ijacph,ijacmh,ijach,ijac0
cc
cc        real(8)    :: vxip,vxim,vxjp,vxjm,vxkp,vxkm
cc     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
cc     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
cc
cc        real(8)    :: bxip,bxim,bxjp,bxjm,bxkp,bxkm
cc     .               ,byip,byim,byjp,byjm,bykp,bykm
cc     .               ,bzip,bzim,bzjp,bzjm,bzkp,bzkm
cc
ccc     Begin program
cc
ccc     Defaults
cc
cc        ip = i+1
cc        im = i-1
cc        jp = j+1
cc        jm = j-1
cc        kp = k+1
cc        km = k-1
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        idhx = 0.5/grid_params%dxh(ig)
cc        idhy = 0.5/grid_params%dyh(jg)
cc        idhz = 0.5/grid_params%dzh(kg)
cc
cc        jac  = gmetric%grid(igx)%jac(i,j,k)
cc
cccc        sing_point = isSP(i,j,k,igx,igy,igz)
cc
ccc     Exceptions
cc
cc        !X
cccc        igp = ig+(ip-i)
cccc        igm = ig-(i-im)
cccc
cccc        dh(1)= (grid_params%xx(igp)-grid_params%xx(igm))
cc
cc
cc        if (hex == 1) then
cc          idhx = 1./dx(ig)
cc          im = i
cc        endif
cc
cc        jacip = face_add(ip,j,k,gmetric%grid(igx)%jac,1)
cc        jacim = face_add(im,j,k,gmetric%grid(igx)%jac,1)
cc
cc        jacjp = face_add(i,jp
cc          jacip  = gmetric%grid(igx)%jac(ip,j,k)
cc          jacim  = gmetric%grid(igx)%jac(i ,j,k)
cc          jacjp  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
cc     .                 +gmetric%grid(igx)%jac(i ,jp,k))
cc          jacjm  = 0.5*(gmetric%grid(igx)%jac(ip,jm,k)
cc     .                 +gmetric%grid(igx)%jac(i ,jm,k))
cc          jackp  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
cc     .                 +gmetric%grid(igx)%jac(i ,j,kp))
cc          jackm  = 0.5*(gmetric%grid(igx)%jac(ip,j,km)
cc     .                 +gmetric%grid(igx)%jac(i ,j,km))
cc
cc          if (isSP(i+1,j,k,igx,igy,igz)) then
cc            jacjp = SP_flsv
cc            jacjm = SP_flsv
cc            jackp = SP_flsv
cc            jackm = SP_flsv
cccc            write (*,*) 'DIAG -- curl_bxv'
cccc            jacjp = 1d0!SP_flsv
cccc            jacjm = 1d0!SP_flsv
cccc            jackp = 1d0!SP_flsv
cccc            jackm = 1d0!SP_flsv
cc          endif
cc
cc          vxip = vv(ip,j,k,1)
cc          vxim = vv(i ,j,k,1)
cc          vyip = vv(ip,j,k,2)
cc          vyim = vv(i ,j,k,2)
cc          vzip = vv(ip,j,k,3)
cc          vzim = vv(i ,j,k,3)
cc
cc          vxjp = 0.5*(vv(ip,jp,k,1)+vv(i,jp,k,1))
cc          vxjm = 0.5*(vv(ip,jm,k,1)+vv(i,jm,k,1))
cc          vyjp = 0.5*(vv(ip,jp,k,2)+vv(i,jp,k,2))
cc          vyjm = 0.5*(vv(ip,jm,k,2)+vv(i,jm,k,2))
cc          vzjp = 0.5*(vv(ip,jp,k,3)+vv(i,jp,k,3))
cc          vzjm = 0.5*(vv(ip,jm,k,3)+vv(i,jm,k,3))
cc
cc          vxkp = 0.5*(vv(ip,j,kp,1)+vv(i,j,kp,1))
cc          vxkm = 0.5*(vv(ip,j,km,1)+vv(i,j,km,1))
cc          vykp = 0.5*(vv(ip,j,kp,2)+vv(i,j,kp,2))
cc          vykm = 0.5*(vv(ip,j,km,2)+vv(i,j,km,2))
cc          vzkp = 0.5*(vv(ip,j,kp,3)+vv(i,j,kp,3))
cc          vzkm = 0.5*(vv(ip,j,km,3)+vv(i,j,km,3))
cc
cc          bxip = bb(ip,j,k,1)
cc          bxim = bb(i ,j,k,1)
cc          byip = bb(ip,j,k,2)
cc          byim = bb(i ,j,k,2)
cc          bzip = bb(ip,j,k,3)
cc          bzim = bb(i ,j,k,3)
cc
cc          bxjp = 0.5*(bb(ip,jp,k,1)+bb(i,jp,k,1))
cc          bxjm = 0.5*(bb(ip,jm,k,1)+bb(i,jm,k,1))
cc          byjp = 0.5*(bb(ip,jp,k,2)+bb(i,jp,k,2))
cc          byjm = 0.5*(bb(ip,jm,k,2)+bb(i,jm,k,2))
cc          bzjp = 0.5*(bb(ip,jp,k,3)+bb(i,jp,k,3))
cc          bzjm = 0.5*(bb(ip,jm,k,3)+bb(i,jm,k,3))
cc
cc          bxkp = 0.5*(bb(ip,j,kp,1)+bb(i,j,kp,1))
cc          bxkm = 0.5*(bb(ip,j,km,1)+bb(i,j,km,1))
cc          bykp = 0.5*(bb(ip,j,kp,2)+bb(i,j,kp,2))
cc          bykm = 0.5*(bb(ip,j,km,2)+bb(i,j,km,2))
cc          bzkp = 0.5*(bb(ip,j,kp,3)+bb(i,j,kp,3))
cc          bzkm = 0.5*(bb(ip,j,km,3)+bb(i,j,km,3))
cc
cc        case (2)
cc
cc          idhy = 1./dy(jg)
cc          jm = j
cc
cc          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
cc     .                 +gmetric%grid(igx)%jac(ip,j ,k))
cc          jacim  = 0.5*(gmetric%grid(igx)%jac(im,jp,k)
cc     .                 +gmetric%grid(igx)%jac(im,j ,k))
cc          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
cc          jacjm  = gmetric%grid(igx)%jac(i,j ,k)
cc          jackp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
cc     .                 +gmetric%grid(igx)%jac(i,j ,kp))
cc          jackm  = 0.5*(gmetric%grid(igx)%jac(i,jp,km)
cc     .                 +gmetric%grid(igx)%jac(i,j ,km))
cc
cc          vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1))*0.5
cc          vxim = (vv(im,j,k,1)+vv(im,jp,k,1))*0.5
cc          vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2))*0.5
cc          vyim = (vv(im,j,k,2)+vv(im,jp,k,2))*0.5
cc          vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3))*0.5
cc          vzim = (vv(im,j,k,3)+vv(im,jp,k,3))*0.5
cc
cc          vxjp = vv(i,jp,k,1)
cc          vxjm = vv(i,j ,k,1)
cc          vyjp = vv(i,jp,k,2)
cc          vyjm = vv(i,j ,k,2)
cc          vzjp = vv(i,jp,k,3)
cc          vzjm = vv(i,j ,k,3)
cc
cc          vxkp = (vv(i,j,kp,1)+vv(i,jp,kp,1))*0.5
cc          vxkm = (vv(i,j,km,1)+vv(i,jp,km,1))*0.5
cc          vykp = (vv(i,j,kp,2)+vv(i,jp,kp,2))*0.5
cc          vykm = (vv(i,j,km,2)+vv(i,jp,km,2))*0.5
cc          vzkp = (vv(i,j,kp,3)+vv(i,jp,kp,3))*0.5
cc          vzkm = (vv(i,j,km,3)+vv(i,jp,km,3))*0.5
cc
cc          bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1))*0.5
cc          bxim = (bb(im,j,k,1)+bb(im,jp,k,1))*0.5
cc          byip = (bb(ip,j,k,2)+bb(ip,jp,k,2))*0.5
cc          byim = (bb(im,j,k,2)+bb(im,jp,k,2))*0.5
cc          bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3))*0.5
cc          bzim = (bb(im,j,k,3)+bb(im,jp,k,3))*0.5
cc
cc          bxjp = bb(i,jp,k,1)
cc          bxjm = bb(i,j ,k,1)
cc          byjp = bb(i,jp,k,2)
cc          byjm = bb(i,j ,k,2)
cc          bzjp = bb(i,jp,k,3)
cc          bzjm = bb(i,j ,k,3)
cc
cc          bxkp = (bb(i,j,kp,1)+bb(i,jp,kp,1))*0.5
cc          bxkm = (bb(i,j,km,1)+bb(i,jp,km,1))*0.5
cc          bykp = (bb(i,j,kp,2)+bb(i,jp,kp,2))*0.5
cc          bykm = (bb(i,j,km,2)+bb(i,jp,km,2))*0.5
cc          bzkp = (bb(i,j,kp,3)+bb(i,jp,kp,3))*0.5
cc          bzkm = (bb(i,j,km,3)+bb(i,jp,km,3))*0.5
cc
cc        case (3)
cc          idhz = 1./dz(kg)
cc          km = k
cc
cc          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
cc     .                 +gmetric%grid(igx)%jac(ip,j,k ))
cc          jacim  = 0.5*(gmetric%grid(igx)%jac(im,j,kp)
cc     .                 +gmetric%grid(igx)%jac(im,j,k ))
cc          jacjp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
cc     .                 +gmetric%grid(igx)%jac(i,jp,k ))
cc          jacjm  = 0.5*(gmetric%grid(igx)%jac(i,jm,kp)
cc     .                 +gmetric%grid(igx)%jac(i,jm,k ))
cc          jackp  = gmetric%grid(igx)%jac(i,j,kp)
cc          jackm  = gmetric%grid(igx)%jac(i,j,km)
cc
cc          vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1))*0.5
cc          vxim = (vv(im,j,k,1)+vv(im,j,kp,1))*0.5
cc          vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2))*0.5
cc          vyim = (vv(im,j,k,2)+vv(im,j,kp,2))*0.5
cc          vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3))*0.5
cc          vzim = (vv(im,j,k,3)+vv(im,j,kp,3))*0.5
cc
cc          vxjp = (vv(i,jp,k,1)+vv(i,jp,kp,1))*0.5
cc          vxjm = (vv(i,jm,k,1)+vv(i,jm,kp,1))*0.5
cc          vyjp = (vv(i,jp,k,2)+vv(i,jp,kp,2))*0.5
cc          vyjm = (vv(i,jm,k,2)+vv(i,jm,kp,2))*0.5
cc          vzjp = (vv(i,jp,k,3)+vv(i,jp,kp,3))*0.5
cc          vzjm = (vv(i,jm,k,3)+vv(i,jm,kp,3))*0.5
cc
cc          vxkp = vv(i,j,kp,1)
cc          vxkm = vv(i,j,k ,1)
cc          vykp = vv(i,j,kp,2)
cc          vykm = vv(i,j,k ,2)
cc          vzkp = vv(i,j,kp,3)
cc          vzkm = vv(i,j,k ,3)
cc
cc          bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1))*0.5
cc          bxim = (bb(im,j,k,1)+bb(im,j,kp,1))*0.5
cc          byip = (bb(ip,j,k,2)+bb(ip,j,kp,2))*0.5
cc          byim = (bb(im,j,k,2)+bb(im,j,kp,2))*0.5
cc          bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3))*0.5
cc          bzim = (bb(im,j,k,3)+bb(im,j,kp,3))*0.5
cc
cc          bxjp = (bb(i,jp,k,1)+bb(i,jp,kp,1))*0.5
cc          bxjm = (bb(i,jm,k,1)+bb(i,jm,kp,1))*0.5
cc          byjp = (bb(i,jp,k,2)+bb(i,jp,kp,2))*0.5
cc          byjm = (bb(i,jm,k,2)+bb(i,jm,kp,2))*0.5
cc          bzjp = (bb(i,jp,k,3)+bb(i,jp,kp,3))*0.5
cc          bzjm = (bb(i,jm,k,3)+bb(i,jm,kp,3))*0.5
cc
cc          bxkp = bb(i,j,kp,1)
cc          bxkm = bb(i,j,k ,1)
cc          bykp = bb(i,j,kp,2)
cc          bykm = bb(i,j,k ,2)
cc          bzkp = bb(i,j,kp,3)
cc          bzkm = bb(i,j,k ,3)
cc
cc        case default
cc          
cc          jacip  = gmetric%grid(igx)%jac(ip,j,k)
cc          jacim  = gmetric%grid(igx)%jac(im,j,k)
cc          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
cc          jacjm  = gmetric%grid(igx)%jac(i,jm,k)
cc          jackp  = gmetric%grid(igx)%jac(i,j,kp)
cc          jackm  = gmetric%grid(igx)%jac(i,j,km)
cc
cc          vxip = vv(ip,j,k,1)
cc          vxim = vv(im,j,k,1)
cc          vyip = vv(ip,j,k,2)
cc          vyim = vv(im,j,k,2)
cc          vzip = vv(ip,j,k,3)
cc          vzim = vv(im,j,k,3)
cc
cc          vxjp = vv(i,jp,k,1)
cc          vxjm = vv(i,jm,k,1)
cc          vyjp = vv(i,jp,k,2)
cc          vyjm = vv(i,jm,k,2)
cc          vzjp = vv(i,jp,k,3)
cc          vzjm = vv(i,jm,k,3)
cc
cc          vxkp = vv(i,j,kp,1)
cc          vxkm = vv(i,j,km,1)
cc          vykp = vv(i,j,kp,2)
cc          vykm = vv(i,j,km,2)
cc          vzkp = vv(i,j,kp,3)
cc          vzkm = vv(i,j,km,3)
cc
cc          bxip = bb(ip,j,k,1)
cc          bxim = bb(im,j,k,1)
cc          byip = bb(ip,j,k,2)
cc          byim = bb(im,j,k,2)
cc          bzip = bb(ip,j,k,3)
cc          bzim = bb(im,j,k,3)
cc
cc          bxjp = bb(i,jp,k,1)
cc          bxjm = bb(i,jm,k,1)
cc          byjp = bb(i,jp,k,2)
cc          byjm = bb(i,jm,k,2)
cc          bzjp = bb(i,jp,k,3)
cc          bzjm = bb(i,jm,k,3)
cc
cc          bxkp = bb(i,j,kp,1)
cc          bxkm = bb(i,j,km,1)
cc          bykp = bb(i,j,kp,2)
cc          bykm = bb(i,j,km,2)
cc          bzkp = bb(i,j,kp,3)
cc          bzkm = bb(i,j,km,3)
cc
cc        end select
cc
ccc     Components
cc
cc        ijacip  = 1d0/jacip
cc        ijacim  = 1d0/jacim
cc        ijacjp  = 1d0/jacjp
cc        ijacjm  = 1d0/jacjm
cc        ijackp  = 1d0/jackp
cc        ijackm  = 1d0/jackm
cc
cccc        write (*,*) 'DIAG -- curl_bxv'
cccc        ijacip  = 1d0
cccc        ijacim  = 1d0
cccc        ijacjp  = 1d0
cccc        ijacjm  = 1d0
cccc        ijackp  = 1d0
cccc        ijackm  = 1d0
cc
cc        !component 1
cc
cc        flxjp = ( vyjp*bxjp-vxjp*byjp )*ijacjp
cc        flxjm = ( vyjm*bxjm-vxjm*byjm )*ijacjm
cc
cc        flxkp = ( vzkp*bxkp-vxkp*bzkp )*ijackp
cc        flxkm = ( vzkm*bxkm-vxkm*bzkm )*ijackm
cc
cc        cnv(1) =  (flxjp-flxjm)*idhy
cc     .           +(flxkp-flxkm)*idhz
cc
cc        !component 2
cc
cc        flxip = ( vxip*byip-vyip*bxip )*ijacip
cc        flxim = ( vxim*byim-vyim*bxim )*ijacim
cc
cc        flxkp = ( vzkp*bykp-vykp*bzkp )*ijackp
cc        flxkm = ( vzkm*bykm-vykm*bzkm )*ijackm
cc
cc        cnv(2) =  (flxip-flxim)*idhx
cc     .           +(flxkp-flxkm)*idhz
cc
cc        !component 3
cc
cc        flxip = ( vxip*bzip-vzip*bxip )*ijacip
cc        flxim = ( vxim*bzim-vzim*bxim )*ijacim
cc
cc        flxjp = ( vyjp*bzjp-vzjp*byjp )*ijacjp
cc        flxjm = ( vyjm*bzjm-vzjm*byjm )*ijacjm
cc
cc        cnv(3) =  (flxip-flxim)*idhx
cc     .           +(flxjp-flxjm)*idhy
cc
cc      contains
cc
cc      function face_add(i,j,k,array,dir) result (sum)
ccc     -----------------------------------------------------------------
ccc     This function adds contributions of 'array' at a given face from 
ccc     different points in the 27-point stencil, depending on the values 
ccc     of hex, hey, hez. The face is determined by the direction 'dir'
ccc     (=1,2,3) and the corresponding index (i,j,k).
ccc     -----------------------------------------------------------------
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,dir
cc        real(8)    :: sum,array(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc        sum = array(i,j,k)
cc
cc        select case(dir)
cc        case(1)
cc
cc          if (hey == 1) then
cc            sum = sum + array(i,jp,k)
cc            coeff = 0.5
cc          endif
cc
cc          if (hez == 1) then
cc            sum = sum + array(i,j,kp)
cc            coeff = 0.5
cc          endif
cc
cc          if (hey == 1 .and. hez == 1) then
cc            sum = sum + array(i,jp,kp)
cc            coeff = 0.25
cc          endif
cc
cc          sum = coeff*sum
cc
cc        case(2)
cc
cc          if (hex == 1) then
cc            sum = sum + array(ip,j,k)
cc            coeff = 0.5
cc          endif
cc
cc          if (hez == 1) then
cc            sum = sum + array(i,j,kp)
cc            coeff = 0.5
cc          endif
cc
cc          if (hex == 1 .and. hez == 1) then
cc            sum = sum + array(ip,j,kp)
cc            coeff = 0.25
cc          endif
cc
cc          sum = coeff*sum
cc
cc        case(3)
cc
cc          if (hex == 1) then
cc            sum = sum + array(ip,j,k)
cc            coeff = 0.5
cc          endif
cc
cc          if (hey == 1) then
cc            sum = sum + array(i,jp,k)
cc            coeff = 0.5
cc          endif
cc
cc          if (hex == 1 .and. hey == 1) then
cc            sum = sum + array(ip,jp,k)
cc            coeff = 0.25
cc          endif
cc
cc          sum = coeff*sum
cc
cc        end select
cc
cc      end function face_add
cc
cc      end function curl_bxv2

ccc     curlcurl_ndiff
ccc     ###############################################################
cc      function curlcurl_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,aa,vol)
cc     .         result(cnv)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates curl(curl(A)) in general non-orthogonal
ccc     coordinates, preserving the SPD property. The vector A is
ccc     covariant, and returns contravariant components.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cnv(3)
cc
cc      logical,optional,intent(IN) :: vol
cc
ccc     Local variables
cc
cc      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
cc
cc      real(8)    :: dwa,dwb,dwc,dwd,dwe,dwf,dwg,dwh,dwi,dwj,dwk
cc
cc      real(8),pointer,dimension(:,:,:,:,:) :: gsub
cc
cc      logical    :: vol_wgt
cc
cc      real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
cc     .             ,ay(0:nx+1,0:ny+1,0:nz+1)
cc     .             ,az(0:nx+1,0:ny+1,0:nz+1)
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
cc      vol_wgt = .true.
cc      if (PRESENT(vol)) vol_wgt = vol
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
cc      ax = aa(:,:,:,1)
cc      ay = aa(:,:,:,2)
cc      az = aa(:,:,:,3)
cc
cc      gsub => gmetric%grid(igx)%gsub
cc
cc      !Component 1
cc
cc      d_yy_kp = 0.5*(gsub(i,j,k,2,2)+gsub(i,j,kp,2,2))
cc      d_yy_km = 0.5*(gsub(i,j,k,2,2)+gsub(i,j,km,2,2))
cc
cc      d_zz_jp = 0.5*(gsub(i,j,k,3,3)+gsub(i,jp,k,3,3))
cc      d_zz_jm = 0.5*(gsub(i,j,k,3,3)+gsub(i,jm,k,3,3))
cc
cc      d_yz_jpkp = 0.25*(gsub(i,j,k ,2,3)+gsub(i,jp,k ,2,3)
cc     .                 +gsub(i,j,kp,2,3)+gsub(i,jp,kp,2,3))
cc      d_yz_jpkm = 0.25*(gsub(i,j,k ,2,3)+gsub(i,jp,k ,2,3)
cc     .                 +gsub(i,j,km,2,3)+gsub(i,jp,km,2,3))
cc      d_yz_jmkp = 0.25*(gsub(i,j,k ,2,3)+gsub(i,jm,k ,2,3)
cc     .                 +gsub(i,j,kp,2,3)+gsub(i,jm,kp,2,3))
cc      d_yz_jmkm = 0.25*(gsub(i,j,k ,2,3)+gsub(i,jm,k ,2,3)
cc     .                 +gsub(i,j,km,2,3)+gsub(i,jm,km,2,3))
cc
cc      d_xy_kp = 0.5*(gsub(i,j,k,1,2)+gsub(i,j,kp,1,2))
cc      d_xy_km = 0.5*(gsub(i,j,k,1,2)+gsub(i,j,km,1,2))
cc
cc      d_yz_kp = 0.5*(gsub(i,j,k,2,3)+gsub(i,j,kp,2,3))
cc      d_yz_km = 0.5*(gsub(i,j,k,2,3)+gsub(i,j,km,2,3))
cc
cc      d_xz_jp = 0.5*(gsub(i,j,k,1,3)+gsub(i,jp,k,1,3))
cc      d_xz_jm = 0.5*(gsub(i,j,k,1,3)+gsub(i,jm,k,1,3))
cc
cc      d_yz_jp = 0.5*(gsub(i,j,k,2,3)+gsub(i,jp,k,2,3))
cc      d_yz_jm = 0.5*(gsub(i,j,k,2,3)+gsub(i,jm,k,2,3))
cc
cc      dwa = -dxh(ig)*dyh(jg)
cc     .            *(d_yy_kp*(ax(i,j,kp)-ax(i,j,k ))/dz(kg  )
cc     .             -d_yy_km*(ax(i,j,k )-ax(i,j,km))/dz(kg-1) )
cc      dwb = -dxh(ig)*dzh(kg)
cc     .            *(d_zz_jp*(ax(i,jp,k)-ax(i,j ,k))/dy(jg  )
cc     .             -d_zz_jm*(ax(i,j ,k)-ax(i,jm,k))/dy(jg-1) )
cc      dwc = -dxh(ig)/2.
cc     .            *(-d_yz_jpkp*(ax(i,jp,kp)-ax(i,j ,k ))
cc     .              +d_yz_jmkm*(ax(i,j ,k )-ax(i,jm,km))
cc     .              +d_yz_jpkm*(ax(i,jp,km)-ax(i,j ,k ))
cc     .              -d_yz_jmkp*(ax(i,j ,k )-ax(i,jm,kp)))
cc      dwd =  dyh(jg)/4.
cc     .            *(d_yy_kp*(az(ip,j,kp)-az(im,j,kp)
cc     .                      +az(ip,j,k )-az(im,j,k ))
cc     .             -d_yy_km*(az(ip,j,km)-az(im,j,km)
cc     .                      +az(ip,j,k )-az(im,j,k )))
cc      dwe = -dxh(ig)/4.
cc     .            *(d_xy_kp*(az(i,jp,kp)-az(i,jm,kp)
cc     .                      +az(i,jp,k )-az(i,jm,k ))
cc     .             -d_xy_km*(az(i,jp,km)-az(i,jm,km)
cc     .                      +az(i,jp,k )-az(i,jm,k )))
cc      dwf =  dxh(ig)*dyh(jg)
cc     .            *(d_xy_kp*(ay(i,j,kp)-ay(i,j,k ))/dz(kg)
cc     .             -d_xy_km*(ay(i,j,k )-ay(i,j,km))/dz(kg-1))
cc      dwg = -dyh(jg)/4.
cc     .            *(d_yz_kp*(ay(ip,j,kp)-ay(im,j,kp)
cc     .                      +ay(ip,j,k )-ay(im,j,k ))
cc     .             -d_yz_km*(ay(ip,j,km)-ay(im,j,km)
cc     .                      +ay(ip,j,k )-ay(im,j,k )))
cc      dwh =  dzh(kg)/4.
cc     .            *(d_zz_jp*(ay(ip,jp,k)-ay(im,jp,k)
cc     .                      +ay(ip,j ,k)-ay(im,j ,k))
cc     .             -d_zz_jm*(ay(ip,jm,k)-ay(im,jm,k)
cc     .                      +ay(ip,j ,k)-ay(im,j ,k)))
cc      dwi =  dxh(ig)*dzh(kg)
cc     .            *(d_xz_jp*(az(i,jp,k)-az(i,j ,k))/dy(jg  )
cc     .             -d_xz_jm*(az(i,j ,k)-az(i,jm,k))/dy(jg-1))
cc      dwj = -dxh(ig)/4.
cc     .            *(d_xz_jp*(ay(i,jp,kp)-ay(i,jp,km)
cc     .                      +ay(i,j ,kp)-ay(i,j ,km))
cc     .             -d_xz_jm*(ay(i,jm,kp)-ay(i,jm,km)
cc     .                      +ay(i,j ,kp)-ay(i,j ,km)))
cc      dwk = -dzh(kg)/4.
cc     .            *(d_yz_jp*(az(ip,jp,k)-az(im,jp,k)
cc     .                      +az(ip,j ,k)-az(im,j ,k))
cc     .             -d_yz_jm*(az(ip,jm,k)-az(im,jm,k)
cc     .                      +az(ip,j ,k)-az(im,j ,k)))
cc
cc      cnv(1) = dwa+dwb+dwc+dwd+dwe+dwf+dwg+dwh+dwi+dwj+dwk
cc
cc      !Component 2
cc
cc      d_xx_kp = 0.5*(gsub(i,j,k,1,1)+gsub(i,j,kp,1,1))
cc      d_xx_km = 0.5*(gsub(i,j,k,1,1)+gsub(i,j,km,1,1))
cc
cc      d_zz_ip = 0.5*(gsub(i,j,k,3,3)+gsub(ip,j,k,3,3))
cc      d_zz_im = 0.5*(gsub(i,j,k,3,3)+gsub(im,j,k,3,3))
cc
cc      d_xz_ipkp = 0.25*(gsub(i,j,k ,1,3)+gsub(ip,j,k ,1,3)
cc     .                 +gsub(i,j,kp,1,3)+gsub(ip,j,kp,1,3))
cc      d_xz_ipkm = 0.25*(gsub(i,j,k ,1,3)+gsub(ip,j,k ,1,3)
cc     .                 +gsub(i,j,km,1,3)+gsub(ip,j,km,1,3))
cc      d_xz_imkp = 0.25*(gsub(i,j,k ,1,3)+gsub(im,j,k ,1,3)
cc     .                 +gsub(i,j,kp,1,3)+gsub(im,j,kp,1,3))
cc      d_xz_imkm = 0.25*(gsub(i,j,k ,1,3)+gsub(im,j,k ,1,3)
cc     .                 +gsub(i,j,km,1,3)+gsub(im,j,km,1,3))
cc
cc      d_xy_kp = 0.5*(gsub(i,j,k,1,2)+gsub(i,j,kp,1,2))
cc      d_xy_km = 0.5*(gsub(i,j,k,1,2)+gsub(i,j,km,1,2))
cc
cc      d_yz_ip = 0.5*(gsub(i,j,k,2,3)+gsub(ip,j,k,2,3))
cc      d_yz_im = 0.5*(gsub(i,j,k,2,3)+gsub(im,j,k,2,3))
cc
cc      d_xz_kp = 0.5*(gsub(i,j,k,1,3)+gsub(i,j,kp,1,3))
cc      d_xz_km = 0.5*(gsub(i,j,k,1,3)+gsub(i,j,km,1,3))
cc
cc      d_xz_ip = 0.5*(gsub(i,j,k,1,3)+gsub(ip,j,k,1,3))
cc      d_xz_im = 0.5*(gsub(i,j,k,1,3)+gsub(im,j,k,1,3))
cc
cc      dwa = -dzh(kg)*dyh(jg)
cc     .            *(d_zz_ip*(ay(ip,j,k)-ay(i ,j,k))/dx(ig  )
cc     .             -d_zz_im*(ay(i ,j,k)-ay(im,j,k))/dx(ig-1) )
cc      dwb = -dxh(ig)*dyh(jg)
cc     .            *(d_xx_kp*(ay(i,j,kp)-ay(i,j,k ))/dz(kg  )
cc     .             -d_xx_km*(ay(i,j,k )-ay(i,j,km))/dz(kg-1) )
cc      dwc = -dyh(jg)/2.
cc     .            *(-d_xz_ipkp*(ay(ip,j,kp)-ay(i ,j,k ))
cc     .              +d_xz_imkm*(ay(i ,j,k )-ay(im,j,km))
cc     .              +d_xz_ipkm*(ay(ip,j,km)-ay(i ,j,k ))
cc     .              -d_xz_imkp*(ay(i ,j,k )-ay(im,j,kp)))
cc      dwd =  dxh(ig)/4.
cc     .            *(d_xx_kp*(az(i,jp,kp)-az(i,jm,kp)
cc     .                      +az(i,jp,k )-az(i,jm,k ))
cc     .             -d_xx_km*(az(i,jp,km)-az(i,jm,km)
cc     .                      +az(i,jp,k )-az(i,jm,k )))
cc      dwe = -dyh(jg)/4.
cc     .            *(d_xy_kp*(az(ip,j,kp)-az(im,j,kp)
cc     .                      +az(ip,j,k )-az(im,j,k ))
cc     .             -d_xy_km*(az(ip,j,km)-az(im,j,km)
cc     .                      +az(ip,j,k )-az(im,j,k )))
cc      dwf =  dxh(ig)*dyh(jg)
cc     .            *(d_xy_kp*(ax(i,j,kp)-ax(i,j,k ))/dz(kg)
cc     .             -d_xy_km*(ax(i,j,k )-ax(i,j,km))/dz(kg-1))
cc      dwg = -dxh(ig)/4.
cc     .            *(d_xz_kp*(ax(i,jp,kp)-ax(i,jm,kp)
cc     .                      +ax(i,jp,k )-ax(i,jm,k ))
cc     .             -d_xz_km*(ax(i,jp,km)-ax(i,jm,km)
cc     .                      +ax(i,jp,k )-ax(i,jm,k )))
cc      dwh =  dzh(kg)/4.
cc     .            *(d_zz_ip*(ax(ip,jp,k)-ax(ip,jm,k)
cc     .                      +ax(i ,jp,k)-ax(i ,jm,k))
cc     .             -d_zz_im*(ax(im,jp,k)-ax(im,jm,k)
cc     .                      +ax(i ,jp,k)-ax(i ,jm,k)))
cc      dwi =  dyh(jg)*dzh(kg)
cc     .            *(d_yz_ip*(az(ip,j,k)-az(i ,j,k))/dx(ig  )
cc     .             -d_yz_im*(az(i ,j,k)-az(im,j,k))/dx(ig-1))
cc      dwj = -dzh(kg)/4.
cc     .            *(d_xz_ip*(az(ip,jp,k)-az(ip,jm,k)
cc     .                      +az(i ,jp,k)-az(i ,jm,k))
cc     .             -d_xz_im*(az(im,jp,k)-az(im,jm,k)
cc     .                      +az(i ,jp,k)-az(i ,jm,k)))
cc      dwk = -dyh(jg)/4.
cc     .            *(d_yz_ip*(ax(ip,j,kp)-ax(ip,j,km)
cc     .                      +ax(i ,j,kp)-ax(i ,j,km))
cc     .             -d_yz_im*(ax(im,j,kp)-ax(im,j,km)
cc     .                      +ax(i ,j,kp)-ax(i ,j,km)))
cc
cc      cnv(2) = dwa+dwb+dwc+dwd+dwe+dwf+dwg+dwh+dwi+dwj+dwk
cc
cc      !Component 3
cc
cc      d_xx_jp = 0.5*(gsub(i,j,k,1,1)+gsub(i,jp,k,1,1))
cc      d_xx_jm = 0.5*(gsub(i,j,k,1,1)+gsub(i,jm,k,1,1))
cc
cc      d_yy_ip = 0.5*(gsub(i,j,k,2,2)+gsub(ip,j,k,2,2))
cc      d_yy_im = 0.5*(gsub(i,j,k,2,2)+gsub(im,j,k,2,2))
cc
cc      d_xy_ipjp = 0.25*(gsub(i,j ,k,1,2)+gsub(ip,j ,k,1,2)
cc     .                 +gsub(i,jp,k,1,2)+gsub(ip,jp,k,1,2))
cc      d_xy_ipjm = 0.25*(gsub(i,j ,k,1,2)+gsub(ip,j ,k,1,2)
cc     .                 +gsub(i,jm,k,1,2)+gsub(ip,jm,k,1,2))
cc      d_xy_imjp = 0.25*(gsub(i,j ,k,1,2)+gsub(im,j ,k,1,2)
cc     .                 +gsub(i,jp,k,1,2)+gsub(im,jp,k,1,2))
cc      d_xy_imjm = 0.25*(gsub(i,j ,k,1,2)+gsub(im,j ,k,1,2)
cc     .                 +gsub(i,jm,k,1,2)+gsub(im,jm,k,1,2))
cc
cc      d_xy_jp = 0.5*(gsub(i,j,k,1,2)+gsub(i,jp,k,1,2))
cc      d_xy_jm = 0.5*(gsub(i,j,k,1,2)+gsub(i,jm,k,1,2))
cc
cc      d_yz_ip = 0.5*(gsub(i,j,k,2,3)+gsub(ip,j,k,2,3))
cc      d_yz_im = 0.5*(gsub(i,j,k,2,3)+gsub(im,j,k,2,3))
cc
cc      d_xy_ip = 0.5*(gsub(i,j,k,1,2)+gsub(ip,j,k,1,2))
cc      d_xy_im = 0.5*(gsub(i,j,k,1,2)+gsub(im,j,k,1,2))
cc
cc      d_xz_jp = 0.5*(gsub(i,j,k,1,3)+gsub(i,jp,k,1,3))
cc      d_xz_jm = 0.5*(gsub(i,j,k,1,3)+gsub(i,jm,k,1,3))
cc
cc      dwa = -dzh(kg)*dyh(jg)
cc     .            *(d_yy_ip*(az(ip,j,k)-az(i ,j,k))/dx(ig  )
cc     .             -d_yy_im*(az(i ,j,k)-az(im,j,k))/dx(ig-1) )
cc      dwb = -dxh(ig)*dzh(kg)
cc     .            *(d_xx_jp*(az(i,jp,k)-az(i,j ,k))/dy(jg  )
cc     .             -d_xx_jm*(az(i,j ,k)-az(i,jm,k))/dy(jg-1) )
cc      dwc = -dzh(kg)/2.
cc     .            *(-d_xy_ipjp*(az(ip,jp,k)-az(i ,j ,k))
cc     .              +d_xy_imjm*(az(i ,j ,k)-az(im,jm,k))
cc     .              +d_xy_ipjm*(az(ip,jm,k)-az(i ,j ,k))
cc     .              -d_xy_imjp*(az(i ,j ,k)-az(im,jp,k)))
cc      dwd =  dxh(ig)/4.
cc     .            *(d_xx_jp*(ay(i,jp,kp)-ay(i,jp,km)
cc     .                      +ay(i,j ,kp)-ay(i,j ,km))
cc     .             -d_xx_jm*(ay(i,jm,kp)-ay(i,jm,km)
cc     .                      +ay(i,j ,kp)-ay(i,j ,km)))
cc      dwe=  -dxh(ig)/4.
cc     .            *(d_xy_jp*(ax(i,jp,kp)-ax(i,jp,km)
cc     .                      +ax(i,j ,kp)-ax(i,j ,km))
cc     .             -d_xy_jm*(ax(i,jm,kp)-ax(i,jm,km)
cc     .                      +ax(i,j ,kp)-ax(i,j ,km)))
cc      dwf =  dxh(ig)*dzh(kg)
cc     .            *(d_xz_jp*(ax(i,jp,k)-ax(i,j ,k))/dy(jg)
cc     .             -d_xz_jm*(ax(i,j ,k)-ax(i,jm,k))/dy(jg-1))
cc      dwg = -dzh(kg)/4.
cc     .            *(d_xz_jp*(ay(ip,jp,k)-ay(im,jp,k)
cc     .                      +ay(ip,j ,k)-ay(im,j ,k))
cc     .             -d_xz_jm*(ay(ip,jm,k)-ay(im,jm,k)
cc     .                      +ay(ip,j ,k)-ay(im,j ,k)))
cc      dwh =  dyh(jg)/4.
cc     .            *(d_yy_ip*(ax(ip,j,kp)-ax(ip,j,km)
cc     .                      +ax(i ,j,kp)-ax(i ,j,km))
cc     .             -d_yy_im*(ax(im,j,kp)-ax(im,j,km)
cc     .                      +ax(i ,j,kp)-ax(i ,j,km)))
cc      dwi =  dyh(jg)*dzh(kg)
cc     .            *(d_yz_ip*(ay(ip,j,k)-ay(i ,j,k))/dx(ig  )
cc     .             -d_yz_im*(ay(i ,j,k)-ay(im,j,k))/dx(ig-1))
cc      dwj = -dzh(kg)/4.
cc     .            *(d_yz_ip*(ax(ip,jp,k)-ax(ip,jm,k)
cc     .                      +ax(i ,jp,k)-ax(i ,jm,k))
cc     .             -d_yz_im*(ax(im,jp,k)-ax(im,jm,k)
cc     .                      +ax(i ,jp,k)-ax(i ,jm,k)))
cc      dwk = -dyh(jg)/4.
cc     .            *(d_xy_ip*(ay(ip,j,kp)-ay(ip,j,km)
cc     .                      +ay(i ,j,kp)-ay(i ,j,km))
cc     .             -d_xy_im*(ay(im,j,kp)-ay(im,j,km)
cc     .                      +ay(i ,j,kp)-ay(i ,j,km)))
cc
cc      cnv(3) = dwa+dwb+dwc+dwd+dwe+dwf+dwg+dwh+dwi+dwj+dwk
cc
cc      if (.not.vol_wgt) cnv=cnv/gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End program
cc
cc      end function curlcurl_ndiff

c     curlcurl_ndiff
c     ###############################################################
      function curlcurl_ndiff(i,j,k,nx,ny,nz,igx,igy,igz,aa,vol)
     .         result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cnv(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: il,jl,kl,ii,jj,kk

      real(8)    :: crl(2,2,2,3),cov(3)

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

c     Curl at vertices

      do kl=1,2
        do jl=1,2
          do il=1,2
            ii = i+il-2
            jj = j+jl-2
            kk = k+kl-2
            cnv = curl_vrtx(ii,jj,kk,igx,igy,igz
     .                     ,aa(ii:ii+1,jj:jj+1,kk:kk+1,:)
     .                     ,.true.)

            call transformFromCurvToCurv(ii,jj,kk,igx,igy,igz
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3)
     .                                  ,.false.
     .                                  ,hex=1,hey=1,hez=1)

            crl(il,jl,kl,:) = cov
          enddo
        enddo
      enddo

c     Curl at cell center

      cnv = curl_vrtx(i,j,k,igx,igy,igz,crl,.false.)

      if (vol_wgt) cnv=cnv*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function curlcurl_ndiff

ccc     curlcurl_diff
ccc     ###############################################################
cc      function curlcurl_diff(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az,diff
cc     .                      ,vol)
cc     .         result(cnv)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates curl(diff*curl(A)) in general non-orthogonal
ccc     coordinates, preserving the SPD property. The vector A is
ccc     covariant, and returns contravariant components.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8)    :: ax  (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,ay  (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,az  (0:nx+1,0:ny+1,0:nz+1)
cc     .             ,diff(0:nx+1,0:ny+1,0:nz+1)
cc     .             ,cnv(3)
cc
cc      logical,optional,intent(IN) :: vol
cc
ccc     Local variables
cc
cc      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg
cc
cccc      real(8)    :: x0,xp,xm,y0,yp,ym,z0,zp,zm
cc      real(8)    :: dwa,dwb,dwc,dwd,dwe,dwf,dwg,dwh,dwi,dwj,dwk
cc
cc      real(8),pointer,dimension(:,:,:,:,:) :: gsub
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
cc      logical    :: vol_wgt
cc
ccc     Begin program
cc
cc      vol_wgt = .true.
cc      if (PRESENT(vol)) vol_wgt = vol
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
cc      gsub => gmetric%grid(igx)%gsub
cc
cc      !Component 1
cc
cc      d_yy_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,2)+gsub(i,j,kp,2,2))
cccc     .         *g_sub_elem(2,2,x0,y0,zp,.false.)
cc      d_yy_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,2)+gsub(i,j,km,2,2))
cc
cc      d_zz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,3,3)+gsub(i,jp,k,3,3))
cc      d_zz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,3,3)+gsub(i,jm,k,3,3))
cc
cc      d_yz_jpkp = 4./(1./diff(i,jp,kp)+1./diff(i,jp,k )
cc     .               +1./diff(i,j ,kp)+1./diff(i,j ,k ))
cc     .         *0.25*(gsub(i,j,k ,2,3)+gsub(i,jp,k ,2,3)
cc     .               +gsub(i,j,kp,2,3)+gsub(i,jp,kp,2,3))
cc      d_yz_jpkm = 4./(1./diff(i,jp,km)+1./diff(i,jp,k )
cc     .               +1./diff(i,j ,km)+1./diff(i,j ,k ))
cc     .         *0.25*(gsub(i,j,k ,2,3)+gsub(i,jp,k ,2,3)
cc     .               +gsub(i,j,km,2,3)+gsub(i,jp,km,2,3))
cc      d_yz_jmkp = 4./(1./diff(i,jm,kp)+1./diff(i,jm,k )
cc     .               +1./diff(i,j ,kp)+1./diff(i,j ,k ))
cc     .         *0.25*(gsub(i,j,k ,2,3)+gsub(i,jm,k ,2,3)
cc     .               +gsub(i,j,kp,2,3)+gsub(i,jm,kp,2,3))
cc      d_yz_jmkm = 4./(1./diff(i,jm,km)+1./diff(i,jm,k )
cc     .               +1./diff(i,j ,km)+1./diff(i,j ,k ))
cc     .         *0.25*(gsub(i,j,k ,2,3)+gsub(i,jm,k ,2,3)
cc     .               +gsub(i,j,km,2,3)+gsub(i,jm,km,2,3))
cc
cc      d_xy_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,2)+gsub(i,j,kp,1,2))
cc      d_xy_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,2)+gsub(i,j,km,1,2))
cc
cc      d_yz_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,3)+gsub(i,j,kp,2,3))
cc      d_yz_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,3)+gsub(i,j,km,2,3))
cc
cc      d_xz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,3)+gsub(i,jp,k,1,3))
cc      d_xz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,3)+gsub(i,jm,k,1,3))
cc
cc      d_yz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,3)+gsub(i,jp,k,2,3))
cc      d_yz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,3)+gsub(i,jm,k,2,3))
cc
cc      dwa = -dxh(ig)*dyh(jg)
cc     .            *(d_yy_kp*(ax(i,j,kp)-ax(i,j,k ))/dz(kg  )
cc     .             -d_yy_km*(ax(i,j,k )-ax(i,j,km))/dz(kg-1) )
cc      dwb = -dxh(ig)*dzh(kg)
cc     .            *(d_zz_jp*(ax(i,jp,k)-ax(i,j ,k))/dy(jg  )
cc     .             -d_zz_jm*(ax(i,j ,k)-ax(i,jm,k))/dy(jg-1) )
cc      dwc = -dxh(ig)/2.
cc     .            *(-d_yz_jpkp*(ax(i,jp,kp)-ax(i,j ,k ))
cc     .              +d_yz_jmkm*(ax(i,j ,k )-ax(i,jm,km))
cc     .              +d_yz_jpkm*(ax(i,jp,km)-ax(i,j ,k ))
cc     .              -d_yz_jmkp*(ax(i,j ,k )-ax(i,jm,kp)))
cc      dwd =  dyh(jg)/4.
cc     .            *(d_yy_kp*(az(ip,j,kp)-az(im,j,kp)
cc     .                      +az(ip,j,k )-az(im,j,k ))
cc     .             -d_yy_km*(az(ip,j,km)-az(im,j,km)
cc     .                      +az(ip,j,k )-az(im,j,k )))
cc      dwe = -dxh(ig)/4.
cc     .            *(d_xy_kp*(az(i,jp,kp)-az(i,jm,kp)
cc     .                      +az(i,jp,k )-az(i,jm,k ))
cc     .             -d_xy_km*(az(i,jp,km)-az(i,jm,km)
cc     .                      +az(i,jp,k )-az(i,jm,k )))
cc      dwf =  dxh(ig)*dyh(jg)
cc     .            *(d_xy_kp*(ay(i,j,kp)-ay(i,j,k ))/dz(kg)
cc     .             -d_xy_km*(ay(i,j,k )-ay(i,j,km))/dz(kg-1))
cc      dwg = -dyh(jg)/4.
cc     .            *(d_yz_kp*(ay(ip,j,kp)-ay(im,j,kp)
cc     .                      +ay(ip,j,k )-ay(im,j,k ))
cc     .             -d_yz_km*(ay(ip,j,km)-ay(im,j,km)
cc     .                      +ay(ip,j,k )-ay(im,j,k )))
cc      dwh =  dzh(kg)/4.
cc     .            *(d_zz_jp*(ay(ip,jp,k)-ay(im,jp,k)
cc     .                      +ay(ip,j ,k)-ay(im,j ,k))
cc     .             -d_zz_jm*(ay(ip,jm,k)-ay(im,jm,k)
cc     .                      +ay(ip,j ,k)-ay(im,j ,k)))
cc      dwi =  dxh(ig)*dzh(kg)
cc     .            *(d_xz_jp*(az(i,jp,k)-az(i,j ,k))/dy(jg  )
cc     .             -d_xz_jm*(az(i,j ,k)-az(i,jm,k))/dy(jg-1))
cc      dwj = -dxh(ig)/4.
cc     .            *(d_xz_jp*(ay(i,jp,kp)-ay(i,jp,km)
cc     .                      +ay(i,j ,kp)-ay(i,j ,km))
cc     .             -d_xz_jm*(ay(i,jm,kp)-ay(i,jm,km)
cc     .                      +ay(i,j ,kp)-ay(i,j ,km)))
cc      dwk = -dzh(kg)/4.
cc     .            *(d_yz_jp*(az(ip,jp,k)-az(im,jp,k)
cc     .                      +az(ip,j ,k)-az(im,j ,k))
cc     .             -d_yz_jm*(az(ip,jm,k)-az(im,jm,k)
cc     .                      +az(ip,j ,k)-az(im,j ,k)))
cc
cc      cnv(1) = dwa+dwb+dwc+dwd+dwe+dwf+dwg+dwh+dwi+dwj+dwk
cc
cc      !Component 2
cc
cc      d_xx_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,1)+gsub(i,j,kp,1,1))
cc      d_xx_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,1)+gsub(i,j,km,1,1))
cc
cc      d_zz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,3,3)+gsub(ip,j,k,3,3))
cc      d_zz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,3,3)+gsub(im,j,k,3,3))
cc
cc      d_xz_ipkp = 4./(1./diff(ip,j,kp)+1./diff(ip,j,k )
cc     .               +1./diff(i ,j,kp)+1./diff(i ,j,k ))
cc     .         *0.25*(gsub(i,j,k ,1,3)+gsub(ip,j,k ,1,3)
cc     .               +gsub(i,j,kp,1,3)+gsub(ip,j,kp,1,3))
cc      d_xz_ipkm = 4./(1./diff(ip,j,km)+1./diff(ip,j,k )
cc     .               +1./diff(i ,j,km)+1./diff(i ,j,k ))
cc     .         *0.25*(gsub(i,j,k ,1,3)+gsub(ip,j,k ,1,3)
cc     .               +gsub(i,j,km,1,3)+gsub(ip,j,km,1,3))
cc      d_xz_imkp = 4./(1./diff(im,j,kp)+1./diff(im,j,k )
cc     .               +1./diff(i ,j,kp)+1./diff(i ,j,k ))
cc     .         *0.25*(gsub(i,j,k ,1,3)+gsub(im,j,k ,1,3)
cc     .               +gsub(i,j,kp,1,3)+gsub(im,j,kp,1,3))
cc      d_xz_imkm = 4./(1./diff(im,j,km)+1./diff(im,j,k )
cc     .               +1./diff(i ,j,km)+1./diff(i ,j,k ))
cc     .         *0.25*(gsub(i,j,k ,1,3)+gsub(im,j,k ,1,3)
cc     .               +gsub(i,j,km,1,3)+gsub(im,j,km,1,3))
cc
cc      d_xy_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,2)+gsub(i,j,kp,1,2))
cc      d_xy_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,2)+gsub(i,j,km,1,2))
cc
cc      d_yz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,3)+gsub(ip,j,k,2,3))
cc      d_yz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,3)+gsub(im,j,k,2,3))
cc
cc      d_xz_kp = 2./(1./diff(i,j,kp)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,3)+gsub(i,j,kp,1,3))
cc      d_xz_km = 2./(1./diff(i,j,km)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,3)+gsub(i,j,km,1,3))
cc
cc      d_xz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,3)+gsub(ip,j,k,1,3))
cc      d_xz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,3)+gsub(im,j,k,1,3))
cc
cc      dwa = -dzh(kg)*dyh(jg)
cc     .            *(d_zz_ip*(ay(ip,j,k)-ay(i ,j,k))/dx(ig  )
cc     .             -d_zz_im*(ay(i ,j,k)-ay(im,j,k))/dx(ig-1) )
cc      dwb = -dxh(ig)*dyh(jg)
cc     .            *(d_xx_kp*(ay(i,j,kp)-ay(i,j,k ))/dz(kg  )
cc     .             -d_xx_km*(ay(i,j,k )-ay(i,j,km))/dz(kg-1) )
cc      dwc = -dyh(jg)/2.
cc     .            *(-d_xz_ipkp*(ay(ip,j,kp)-ay(i ,j,k ))
cc     .              +d_xz_imkm*(ay(i ,j,k )-ay(im,j,km))
cc     .              +d_xz_ipkm*(ay(ip,j,km)-ay(i ,j,k ))
cc     .              -d_xz_imkp*(ay(i ,j,k )-ay(im,j,kp)))
cc      dwd =  dxh(ig)/4.
cc     .            *(d_xx_kp*(az(i,jp,kp)-az(i,jm,kp)
cc     .                      +az(i,jp,k )-az(i,jm,k ))
cc     .             -d_xx_km*(az(i,jp,km)-az(i,jm,km)
cc     .                      +az(i,jp,k )-az(i,jm,k )))
cc      dwe = -dyh(jg)/4.
cc     .            *(d_xy_kp*(az(ip,j,kp)-az(im,j,kp)
cc     .                      +az(ip,j,k )-az(im,j,k ))
cc     .             -d_xy_km*(az(ip,j,km)-az(im,j,km)
cc     .                      +az(ip,j,k )-az(im,j,k )))
cc      dwf =  dxh(ig)*dyh(jg)
cc     .            *(d_xy_kp*(ax(i,j,kp)-ax(i,j,k ))/dz(kg)
cc     .             -d_xy_km*(ax(i,j,k )-ax(i,j,km))/dz(kg-1))
cc      dwg = -dxh(ig)/4.
cc     .            *(d_xz_kp*(ax(i,jp,kp)-ax(i,jm,kp)
cc     .                      +ax(i,jp,k )-ax(i,jm,k ))
cc     .             -d_xz_km*(ax(i,jp,km)-ax(i,jm,km)
cc     .                      +ax(i,jp,k )-ax(i,jm,k )))
cc      dwh =  dzh(kg)/4.
cc     .            *(d_zz_ip*(ax(ip,jp,k)-ax(ip,jm,k)
cc     .                      +ax(i ,jp,k)-ax(i ,jm,k))
cc     .             -d_zz_im*(ax(im,jp,k)-ax(im,jm,k)
cc     .                      +ax(i ,jp,k)-ax(i ,jm,k)))
cc      dwi =  dyh(jg)*dzh(kg)
cc     .            *(d_yz_ip*(az(ip,j,k)-az(i ,j,k))/dx(ig  )
cc     .             -d_yz_im*(az(i ,j,k)-az(im,j,k))/dx(ig-1))
cc      dwj = -dzh(kg)/4.
cc     .            *(d_xz_ip*(az(ip,jp,k)-az(ip,jm,k)
cc     .                      +az(i ,jp,k)-az(i ,jm,k))
cc     .             -d_xz_im*(az(im,jp,k)-az(im,jm,k)
cc     .                      +az(i ,jp,k)-az(i ,jm,k)))
cc      dwk = -dyh(jg)/4.
cc     .            *(d_yz_ip*(ax(ip,j,kp)-ax(ip,j,km)
cc     .                      +ax(i ,j,kp)-ax(i ,j,km))
cc     .             -d_yz_im*(ax(im,j,kp)-ax(im,j,km)
cc     .                      +ax(i ,j,kp)-ax(i ,j,km)))
cc
cc      cnv(2) = dwa+dwb+dwc+dwd+dwe+dwf+dwg+dwh+dwi+dwj+dwk
cc
cc      !Component 3
cc
cc      d_xx_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,1)+gsub(i,jp,k,1,1))
cc      d_xx_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,1)+gsub(i,jm,k,1,1))
cc
cc      d_yy_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,2)+gsub(ip,j,k,2,2))
cc      d_yy_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,2)+gsub(im,j,k,2,2))
cc
cc      d_xy_ipjp = 4./(1./diff(ip,jp,k)+1./diff(ip,j,k )
cc     .               +1./diff(i ,jp,k)+1./diff(i ,j,k ))
cc     .         *0.25*(gsub(i,j ,k,1,2)+gsub(ip,j ,k,1,2)
cc     .               +gsub(i,jp,k,1,2)+gsub(ip,jp,k,1,2))
cc      d_xy_ipjm = 4./(1./diff(ip,jm,k)+1./diff(ip,j,k )
cc     .               +1./diff(i ,jm,k)+1./diff(i ,j,k ))
cc     .         *0.25*(gsub(i,j ,k,1,2)+gsub(ip,j ,k,1,2)
cc     .               +gsub(i,jm,k,1,2)+gsub(ip,jm,k,1,2))
cc      d_xy_imjp = 4./(1./diff(im,jp,k)+1./diff(im,j,k )
cc     .               +1./diff(i ,jp,k)+1./diff(i ,j,k ))
cc     .         *0.25*(gsub(i,j ,k,1,2)+gsub(im,j ,k,1,2)
cc     .               +gsub(i,jp,k,1,2)+gsub(im,jp,k,1,2))
cc      d_xy_imjm = 4./(1./diff(im,jm,k)+1./diff(im,j,k )
cc     .               +1./diff(i ,jm,k)+1./diff(i ,j,k ))
cc     .         *0.25*(gsub(i,j ,k,1,2)+gsub(im,j ,k,1,2)
cc     .               +gsub(i,jm,k,1,2)+gsub(im,jm,k,1,2))
cc
cc      d_xy_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,2)+gsub(i,jp,k,1,2))
cc      d_xy_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,2)+gsub(i,jm,k,1,2))
cc
cc      d_yz_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,3)+gsub(ip,j,k,2,3))
cc      d_yz_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,2,3)+gsub(im,j,k,2,3))
cc
cc      d_xy_ip = 2./(1./diff(ip,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,2)+gsub(ip,j,k,1,2))
cc      d_xy_im = 2./(1./diff(im,j,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,2)+gsub(im,j,k,1,2))
cc
cc      d_xz_jp = 2./(1./diff(i,jp,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,3)+gsub(i,jp,k,1,3))
cc      d_xz_jm = 2./(1./diff(i,jm,k)+1./diff(i,j,k))
cc     .         *0.5*(gsub(i,j,k,1,3)+gsub(i,jm,k,1,3))
cc
cc      dwa = -dzh(kg)*dyh(jg)
cc     .            *(d_yy_ip*(az(ip,j,k)-az(i ,j,k))/dx(ig  )
cc     .             -d_yy_im*(az(i ,j,k)-az(im,j,k))/dx(ig-1) )
cc      dwb = -dxh(ig)*dzh(kg)
cc     .            *(d_xx_jp*(az(i,jp,k)-az(i,j ,k))/dy(jg  )
cc     .             -d_xx_jm*(az(i,j ,k)-az(i,jm,k))/dy(jg-1) )
cc      dwc = -dzh(kg)/2.
cc     .            *(-d_xy_ipjp*(az(ip,jp,k)-az(i ,j ,k))
cc     .              +d_xy_imjm*(az(i ,j ,k)-az(im,jm,k))
cc     .              +d_xy_ipjm*(az(ip,jm,k)-az(i ,j ,k))
cc     .              -d_xy_imjp*(az(i ,j ,k)-az(im,jp,k)))
cc      dwd =  dxh(ig)/4.
cc     .            *(d_xx_jp*(ay(i,jp,kp)-ay(i,jp,km)
cc     .                      +ay(i,j ,kp)-ay(i,j ,km))
cc     .             -d_xx_jm*(ay(i,jm,kp)-ay(i,jm,km)
cc     .                      +ay(i,j ,kp)-ay(i,j ,km)))
cc      dwe=  -dxh(ig)/4.
cc     .            *(d_xy_jp*(ax(i,jp,kp)-ax(i,jp,km)
cc     .                      +ax(i,j ,kp)-ax(i,j ,km))
cc     .             -d_xy_jm*(ax(i,jm,kp)-ax(i,jm,km)
cc     .                      +ax(i,j ,kp)-ax(i,j ,km)))
cc      dwf =  dxh(ig)*dzh(kg)
cc     .            *(d_xz_jp*(ax(i,jp,k)-ax(i,j ,k))/dy(jg)
cc     .             -d_xz_jm*(ax(i,j ,k)-ax(i,jm,k))/dy(jg-1))
cc      dwg = -dzh(kg)/4.
cc     .            *(d_xz_jp*(ay(ip,jp,k)-ay(im,jp,k)
cc     .                      +ay(ip,j ,k)-ay(im,j ,k))
cc     .             -d_xz_jm*(ay(ip,jm,k)-ay(im,jm,k)
cc     .                      +ay(ip,j ,k)-ay(im,j ,k)))
cc      dwh =  dyh(jg)/4.
cc     .            *(d_yy_ip*(ax(ip,j,kp)-ax(ip,j,km)
cc     .                      +ax(i ,j,kp)-ax(i ,j,km))
cc     .             -d_yy_im*(ax(im,j,kp)-ax(im,j,km)
cc     .                      +ax(i ,j,kp)-ax(i ,j,km)))
cc      dwi =  dyh(jg)*dzh(kg)
cc     .            *(d_yz_ip*(ay(ip,j,k)-ay(i ,j,k))/dx(ig  )
cc     .             -d_yz_im*(ay(i ,j,k)-ay(im,j,k))/dx(ig-1))
cc      dwj = -dzh(kg)/4.
cc     .            *(d_yz_ip*(ax(ip,jp,k)-ax(ip,jm,k)
cc     .                      +ax(i ,jp,k)-ax(i ,jm,k))
cc     .             -d_yz_im*(ax(im,jp,k)-ax(im,jm,k)
cc     .                      +ax(i ,jp,k)-ax(i ,jm,k)))
cc      dwk = -dyh(jg)/4.
cc     .            *(d_xy_ip*(ay(ip,j,kp)-ay(ip,j,km)
cc     .                      +ay(i ,j,kp)-ay(i ,j,km))
cc     .             -d_xy_im*(ay(im,j,kp)-ay(im,j,km)
cc     .                      +ay(i ,j,kp)-ay(i ,j,km)))
cc
cc      cnv(3) = dwa+dwb+dwc+dwd+dwe+dwf+dwg+dwh+dwi+dwj+dwk
cc
cc      if (.not.vol_wgt) cnv=cnv/gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End program
cc
cc      end function curlcurl_diff

c     curlcurl_diff
c     ###############################################################
      function curlcurl_diff(i,j,k,nx,ny,nz,igx,igy,igz,aa,diff,vol)
     .         result(cnv)

c     ---------------------------------------------------------------
c     Calculates curl(D*curl(A)) in general non-orthogonal
c     coordinates, preserving the SPD property. The vector A is
c     covariant, and returns contravariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cnv(3)
     .             ,diff(0:nx+1,0:ny+1,0:nz+1)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: il,jl,kl,ii,jj,kk

      real(8)    :: crl(2,2,2,3),cov(3),dff

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

c     Curl at vertices

      do kl=1,2
        do jl=1,2
          do il=1,2
            ii = i+il-2
            jj = j+jl-2
            kk = k+kl-2
            cnv = curl_vrtx(ii,jj,kk,igx,igy,igz
     .                     ,aa(ii:ii+1,jj:jj+1,kk:kk+1,:)
     .                     ,.true.)

            call transformFromCurvToCurv(ii,jj,kk,igx,igy,igz
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3)
     .                                  ,.false.
     .                                  ,hex=1,hey=1,hez=1)

            dff = dff + diff(ii,jj,kk)

            crl(il,jl,kl,:) = cov
          enddo
        enddo
      enddo

      crl = crl*0.125*dff

c     Curl at cell center

      cnv = curl_vrtx(i,j,k,igx,igy,igz,crl,.false.)

      if (vol_wgt) cnv=cnv*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function curlcurl_diff

ccc     graddiv
ccc     ###############################################################
cc      function graddiv(i,j,k,nx,ny,nz,igx,igy,igz,aa,vol) result(cov)
cc
ccc     ---------------------------------------------------------------
ccc     Calculates grad(div(A)) in general non-orthogonal
ccc     coordinates, preserving the self-adjoint property. The vector A
ccc     is contravariant, and returns covariant components.
ccc     ---------------------------------------------------------------
cc
cc      implicit none           !For safe fortran
cc
ccc     Call variables
cc
cc      integer    :: i,j,k,nx,ny,nz,igx,igy,igz
cc
cc      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cov(3)
cc
cc      logical,optional,intent(IN) :: vol
cc
ccc     Local variables
cc
cc      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg,igrid
cc
cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc     .             ,divip,divim,divjp,divjm,divkp,divkm
cc     .             ,jacip,jacim,jacjp,jacjm,jackp,jackm
cc     .             ,coeff,jac
cc
cc      logical    :: vol_wgt
cc
ccc     Begin program
cc
cc      igrid = igx
cc
cc      vol_wgt = .true.
cc      if (PRESENT(vol)) vol_wgt = vol
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
cc      jacip  = gmetric%grid(igrid)%jac(ip,j,k)
cc      jacim  = gmetric%grid(igrid)%jac(im,j,k)
cc      jacjp  = gmetric%grid(igrid)%jac(i,jp,k)
cc      jacjm  = gmetric%grid(igrid)%jac(i,jm,k)
cc      jackp  = gmetric%grid(igrid)%jac(i,j,kp)
cc      jackm  = gmetric%grid(igrid)%jac(i,j,km)
cc      jac    = gmetric%grid(igrid)%jac(i,j,k )
cc
ccc     Fluxes at faces for calculation of grad(div(dv))
cc
cc      !!Divergence at faces i+-1/2, etc.
cc      divip = (aa(ip,j ,k,1)-aa(i ,j ,k,1))/dx(ig)
cc     .       +(aa(i ,jp,k,2)-aa(i ,jm,k,2)
cc     .        +aa(ip,jp,k,2)-aa(ip,jm,k,2))/dyh(jg)*0.25
cc     .       +(aa(i ,j,kp,3)-aa(i ,j,km,3)
cc     .        +aa(ip,j,kp,3)-aa(ip,j,km,3))/dzh(kg)*0.25
cc      divim = (aa(i ,j ,k,1)-aa(im,j ,k,1))/dx(ig-1)
cc     .       +(aa(i ,jp,k,2)-aa(i ,jm,k,2)
cc     .        +aa(im,jp,k,2)-aa(im,jm,k,2))/dyh(jg)*0.25
cc     .       +(aa(i ,j,kp,3)-aa(i ,j,km,3)
cc     .        +aa(im,j,kp,3)-aa(im,j,km,3))/dzh(kg)*0.25
cc
cc      divjp = (aa(ip,j ,k,1)-aa(im,j ,k,1)
cc     .        +aa(ip,jp,k,1)-aa(im,jp,k,1))/dxh(ig)*0.25
cc     .       +(aa(i ,jp,k,2)-aa(i ,j ,k,2))/dy(jg)
cc     .       +(aa(i,j ,kp,3)-aa(i,j ,km,3)
cc     .        +aa(i,jp,kp,3)-aa(i,jp,km,3))/dzh(kg)*0.25
cc      divjm = (aa(ip,j ,k,1)-aa(im,j ,k,1)
cc     .        +aa(ip,jm,k,1)-aa(im,jm,k,1))/dxh(ig)*0.25
cc     .       +(aa(i ,j ,k,2)-aa(i ,jm,k,2))/dy(jg-1)
cc     .       +(aa(i,j ,kp,3)-aa(i,j ,km,3)
cc     .        +aa(i,jm,kp,3)-aa(i,jm,km,3))/dzh(kg)*0.25
cc
cc      divkp = (aa(ip,j,k ,1)-aa(im,j,k ,1)
cc     .        +aa(ip,j,kp,1)-aa(im,j,kp,1))/dxh(ig)*0.25
cc     .       +(aa(i,jp,k ,2)-aa(i,jm,k ,2)
cc     .        +aa(i,jp,kp,2)-aa(i,jm,kp,2))/dyh(jg)*0.25
cc     .       +(aa(i,j ,kp,3)-aa(i,j ,k ,3))/dz(kg)
cc      divkm = (aa(ip,j,k ,1)-aa(im,j,k ,1)
cc     .        +aa(ip,j,km,1)-aa(im,j,km,1))/dxh(ig)*0.25
cc     .       +(aa(i,jp,k ,2)-aa(i,jm,k ,2)
cc     .        +aa(i,jp,km,2)-aa(i,jm,km,2))/dyh(jg)*0.25
cc     .       +(aa(i,j ,k ,3)-aa(i,j ,km,3))/dz(kg-1)
cc
cc      flxip = divip/(jac+jacip)
cc      if (.not.isSP(i,j,k,igx,igy,igz)) then
cc        flxim = divim/(jac+jacim)
cc      else
cc        flxim = 0d0
cc      endif
cc
cc      flxjp = divjp/(jac+jacjp)
cc      flxjm = divjm/(jac+jacjm)
cc
cc      flxkp = divkp/(jac+jackp)
cc      flxkm = divkm/(jac+jackm)
cc
ccc     Transform to contravariant
cc
cc      cov(1) = (flxip - flxim)/dxh(ig)
cc      cov(2) = (flxjp - flxjm)/dyh(jg)
cc      cov(3) = (flxkp - flxkm)/dzh(kg)
cc
cc      if (vol_wgt) cov=cov*gmetric%grid(igx)%dvol(i,j,k)
cc
ccc     End program
cc
cc      end function graddiv

c     graddiv
c     ###############################################################
      function graddiv(i,j,k,nx,ny,nz,igx,igy,igz,aa,vol) result(cov)

c     ---------------------------------------------------------------
c     Calculates grad(div(A)) in general non-orthogonal
c     coordinates, preserving the self-adjoint property. The vector A
c     is contravariant, and returns covariant components.
c     ---------------------------------------------------------------

      implicit none           !For safe fortran

c     Call variables

      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

      real(8)    :: aa(0:nx+1,0:ny+1,0:nz+1,3),cov(3)

      logical,optional,intent(IN) :: vol

c     Local variables

      integer    :: il,jl,kl,ii,jj,kk

      real(8)    :: div(2,2,2)

      logical    :: vol_wgt

c     Begin program

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

c     Divergence at vertices

      do kl=1,2
        do jl=1,2
          do il=1,2
            ii = i+il-2
            jj = j+jl-2
            kk = k+kl-2

            div(il,jl,kl) = div_vrtx(ii,jj,kk,igx,igy,igz
     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,1)
     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,2)
     .                         ,aa(ii:ii+1,jj:jj+1,kk:kk+1,3)
     .                         ,.true.)
          enddo
        enddo
      enddo

c     Gradient at cell centers

      cov = grad_vrtx(i,j,k,igx,igy,igz,div,.false.)

      if (vol_wgt) cov=cov*gmetric%grid(igx)%dvol(i,j,k)

c     End program

      end function graddiv

c     fnabla_v
c     #############################################################
      function fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
     .                 ,half_elem,cnv)
     .         result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor components T_l^m of nabla(vec v), with v
c     contravariant, and fills tensor(l,m) at the following positions:
c       + half_elem=0 --> i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

        logical,optional :: cnv

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: tsrc(3,3),jacip,jacjp,jacjm,jackp,jackm
     .               ,jacipkp,jacipkm,jacipjp,jacipjm,jac

        logical    :: cntv

c     Begin program

        igrid = igx

        !Defaults

        if (PRESENT(cnv)) then
          cntv = cnv
        else
          cntv = .false.
        endif

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

          if (cntv) then
            gsuper = 0.5*(gmetric%grid(igrid)%gsup(i ,j,k,:,:)
     .                   +gmetric%grid(igrid)%gsup(ip,j,k,:,:))

            if (isSP(ip,j,k,igx,igy,igz)) gsuper = 0d0
          endif

          vxip = ax(ip,j,k)
          vxim = ax(i ,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(i ,j,k)
          vzip = az(ip,j,k)
          vzim = az(i ,j,k)

          vxjp = 0.5*(ax(i,jp,k)+ax(ip,jp,k))
          vxjm = 0.5*(ax(i,jm,k)+ax(ip,jm,k))
          vyjp = 0.5*(ay(i,jp,k)+ay(ip,jp,k))
          vyjm = 0.5*(ay(i,jm,k)+ay(ip,jm,k))
          vzjp = 0.5*(az(i,jp,k)+az(ip,jp,k))
          vzjm = 0.5*(az(i,jm,k)+az(ip,jm,k))
                                      
          vxkp = 0.5*(ax(i,j,kp)+ax(ip,j,kp))
          vxkm = 0.5*(ax(i,j,km)+ax(ip,j,km))
          vykp = 0.5*(ay(i,j,kp)+ay(ip,j,kp))
          vykm = 0.5*(ay(i,j,km)+ay(ip,j,km))
          vzkp = 0.5*(az(i,j,kp)+az(ip,j,kp))
          vzkm = 0.5*(az(i,j,km)+az(ip,j,km))

          tsrc = 0.5*(nabla_v_src(ip,j,k)
     .               +nabla_v_src(i ,j,k))

        case (2)
          dhy = dy(jg)

          if (cntv)
     .         gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j ,k,:,:)
     .                      +gmetric%grid(igrid)%gsup(i,jp,k,:,:))

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

          if (cntv)
     .         gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,k ,:,:)
     .                      +gmetric%grid(igrid)%gsup(i,j,kp,:,:))

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

          if(cntv) gsuper=gmetric%grid(igrid)%gsup(i,j,k ,:,:)

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

        tensor = tensor + tsrc

        if (cntv) tensor = matmul(gsuper,tensor)

c     End program

      contains

c     nabla_v_src
c     #############################################################
      function nabla_v_src(i,j,k) result(tensor)

c     -------------------------------------------------------------
c     Finds geometric source of tensor nabla(v) at cell (i,j,k)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k
        real(8)    :: tensor(3,3)

c     Local variables

        integer    :: l,m
        real(8)    :: hessian(3,3,3),cnv(3)

c     Begin program

        hessian = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

        cnv = (/ ax(i,j,k),ay(i,j,k),az(i,j,k) /)

        do l=1,3
          do m=1,3
            tensor(l,m) = cnv(1)*hessian(m,1,l)
     .                   +cnv(2)*hessian(m,2,l)
     .                   +cnv(3)*hessian(m,3,l)
     .                   -cnv(m)*(hessian(1,1,l)
     .                           +hessian(2,2,l)
     .                           +hessian(3,3,l))
          enddo
        enddo

c     End program

      end function nabla_v_src

      end function fnabla_v

ccc     fnabla_v
ccc     #############################################################
cc      function fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
cc     .                 ,half_elem,cnv)
cc     .         result(tensor)
ccc     -------------------------------------------------------------
ccc     Calculates the tensor components T_l^m of nabla(vec v) and
ccc     fills tensor(l,m) at the following positions:
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
cc        integer    :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
cc        real(8)    :: tensor(3,3)
cc        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
cc     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
cc     .               ,az(0:nx+1,0:ny+1,0:nz+1)
cc
cc        logical,optional :: cnv
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
cc        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
cc        real(8)    :: vxx,vyy,vzz
cc     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
cc     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
cc     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
cc        real(8)    :: tsrc(3,3),jacip,jacjp,jacjm,jackp,jackm
cc     .               ,jacipkp,jacipkm,jacipjp,jacipjm,jac,jach
cc     .               ,gsuper(3,3)
cc
cc        real(8), pointer, dimension(:,:,:) :: jb
cc
cc        logical    :: cntv
cc
ccc     Begin program
cc
cc        igrid = igx
cc
cc        !Defaults
cc
cc        if (PRESENT(cnv)) then
cc          cntv = cnv
cc        else
cc          cntv = .false.
cc        endif
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
cc        jb => gmetric%grid(igrid)%jac
cc
cc        jach = jb(i,j,k)
cc
cc        if(cntv) gsuper=gmetric%grid(igrid)%gsup(i,j,k ,:,:)
cc
cc        !Exceptions
cc
cc        select case(half_elem)
cc        case (1)
cc          dhx = dx(ig)
cc
cc          if (isSP(ip,j,k,igx,igy,igz)) then
cc            jach = 0d0
cc          else
cc            jach = 0.5*(jb(i,j,k)+jb(ip,j,k))
cc          endif
cc
cc          if (cntv) then
cc            gsuper = 0.5*(gmetric%grid(igrid)%gsup(i ,j,k,:,:)
cc     .                   +gmetric%grid(igrid)%gsup(ip,j,k,:,:))
cc
cc            if (isSP(ip,j,k,igx,igy,igz)) then
cc              gsuper(1,1) = 0d0
cc              gsuper(2,2) = 0d0
cc              gsuper(3,3) = 0d0
cccc              write (*,*) 'DIAG -- fnabla_v',gsuper
cc            endif
cc          endif
cc  
cc          vxip = ax(ip,j,k)/jb(ip,j,k)
cc          vxim = ax(i ,j,k)/jb(i ,j,k)
cc          vyip = ay(ip,j,k)/jb(ip,j,k)
cc          vyim = ay(i ,j,k)/jb(i ,j,k)
cc          vzip = az(ip,j,k)/jb(ip,j,k)
cc          vzim = az(i ,j,k)/jb(i ,j,k)
cc
cc          vxjp = 0.5*(ax(i,jp,k)/jb(i,jp,k)+ax(ip,jp,k)/jb(ip,jp,k))
cc          vxjm = 0.5*(ax(i,jm,k)/jb(i,jm,k)+ax(ip,jm,k)/jb(ip,jm,k))
cc          vyjp = 0.5*(ay(i,jp,k)/jb(i,jp,k)+ay(ip,jp,k)/jb(ip,jp,k))
cc          vyjm = 0.5*(ay(i,jm,k)/jb(i,jm,k)+ay(ip,jm,k)/jb(ip,jm,k))
cc          vzjp = 0.5*(az(i,jp,k)/jb(i,jp,k)+az(ip,jp,k)/jb(ip,jp,k))
cc          vzjm = 0.5*(az(i,jm,k)/jb(i,jm,k)+az(ip,jm,k)/jb(ip,jm,k))
cc                                      
cc          vxkp = 0.5*(ax(i,j,kp)/jb(i,j,kp)+ax(ip,j,kp)/jb(ip,j,kp))
cc          vxkm = 0.5*(ax(i,j,km)/jb(i,j,km)+ax(ip,j,km)/jb(ip,j,km))
cc          vykp = 0.5*(ay(i,j,kp)/jb(i,j,kp)+ay(ip,j,kp)/jb(ip,j,kp))
cc          vykm = 0.5*(ay(i,j,km)/jb(i,j,km)+ay(ip,j,km)/jb(ip,j,km))
cc          vzkp = 0.5*(az(i,j,kp)/jb(i,j,kp)+az(ip,j,kp)/jb(ip,j,kp))
cc          vzkm = 0.5*(az(i,j,km)/jb(i,j,km)+az(ip,j,km)/jb(ip,j,km))
cc
cc          tsrc = 0.5*(nabla_v_src(ip,j,k)
cc     .               +nabla_v_src(i ,j,k))
cc
cc        case (2)
cc          dhy = dy(jg)
cc
cc          jach = 0.5*(jb(i,j,k)+jb(i,jp,k))
cc
cc          if (cntv) gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j ,k,:,:)
cc     .                           +gmetric%grid(igrid)%gsup(i,jp,k,:,:))
cc
cc          vxip = 0.5*(ax(ip,j,k)/jb(ip,j,k)+ax(ip,jp,k)/jb(ip,jp,k))
cc          vxim = 0.5*(ax(im,j,k)/jb(im,j,k)+ax(im,jp,k)/jb(im,jp,k))
cc          vyip = 0.5*(ay(ip,j,k)/jb(ip,j,k)+ay(ip,jp,k)/jb(ip,jp,k))
cc          vyim = 0.5*(ay(im,j,k)/jb(im,j,k)+ay(im,jp,k)/jb(im,jp,k))
cc          vzip = 0.5*(az(ip,j,k)/jb(ip,j,k)+az(ip,jp,k)/jb(ip,jp,k))
cc          vzim = 0.5*(az(im,j,k)/jb(im,j,k)+az(im,jp,k)/jb(im,jp,k))
cc
cc          vxjp = ax(i,jp,k)/jb(i,jp,k)
cc          vxjm = ax(i,j ,k)/jb(i,j ,k)
cc          vyjp = ay(i,jp,k)/jb(i,jp,k)
cc          vyjm = ay(i,j ,k)/jb(i,j ,k)
cc          vzjp = az(i,jp,k)/jb(i,jp,k)
cc          vzjm = az(i,j ,k)/jb(i,j ,k)
cc
cc          vxkp = 0.5*(ax(i,j,kp)/jb(i,j,kp)+ax(i,jp,kp)/jb(i,jp,kp))
cc          vxkm = 0.5*(ax(i,j,km)/jb(i,j,km)+ax(i,jp,km)/jb(i,jp,km))
cc          vykp = 0.5*(ay(i,j,kp)/jb(i,j,kp)+ay(i,jp,kp)/jb(i,jp,kp))
cc          vykm = 0.5*(ay(i,j,km)/jb(i,j,km)+ay(i,jp,km)/jb(i,jp,km))
cc          vzkp = 0.5*(az(i,j,kp)/jb(i,j,kp)+az(i,jp,kp)/jb(i,jp,kp))
cc          vzkm = 0.5*(az(i,j,km)/jb(i,j,km)+az(i,jp,km)/jb(i,jp,km))
cc
cc          tsrc = 0.5*(nabla_v_src(i,jp,k)
cc     .               +nabla_v_src(i,j ,k))
cc
cc        case (3)
cc          dhz = dz(kg)
cc
cc          jach = 0.5*(jb(i,j,k)+jb(i,j,kp))
cc
cc          if (cntv) gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,k ,:,:)
cc     .                           +gmetric%grid(igrid)%gsup(i,j,kp,:,:))
cc
cc          vxip = 0.5*(ax(ip,j,k)/jb(ip,j,k)+ax(ip,j,kp)/jb(ip,j,kp))
cc          vxim = 0.5*(ax(im,j,k)/jb(im,j,k)+ax(im,j,kp)/jb(im,j,kp))
cc          vyip = 0.5*(ay(ip,j,k)/jb(ip,j,k)+ay(ip,j,kp)/jb(ip,j,kp))
cc          vyim = 0.5*(ay(im,j,k)/jb(im,j,k)+ay(im,j,kp)/jb(im,j,kp))
cc          vzip = 0.5*(az(ip,j,k)/jb(ip,j,k)+az(ip,j,kp)/jb(ip,j,kp))
cc          vzim = 0.5*(az(im,j,k)/jb(im,j,k)+az(im,j,kp)/jb(im,j,kp))
cc
cc          vxjp = 0.5*(ax(i,jp,k)/jb(i,jp,k)+ax(i,jp,kp)/jb(i,jp,kp))
cc          vxjm = 0.5*(ax(i,jm,k)/jb(i,jm,k)+ax(i,jm,kp)/jb(i,jm,kp))
cc          vyjp = 0.5*(ay(i,jp,k)/jb(i,jp,k)+ay(i,jp,kp)/jb(i,jp,kp))
cc          vyjm = 0.5*(ay(i,jm,k)/jb(i,jm,k)+ay(i,jm,kp)/jb(i,jm,kp))
cc          vzjp = 0.5*(az(i,jp,k)/jb(i,jp,k)+az(i,jp,kp)/jb(i,jp,kp))
cc          vzjm = 0.5*(az(i,jm,k)/jb(i,jm,k)+az(i,jm,kp)/jb(i,jm,kp))
cc
cc          vxkp = ax(i,j,kp)/jb(i,j,kp)
cc          vxkm = ax(i,j,k )/jb(i,j,k )
cc          vykp = ay(i,j,kp)/jb(i,j,kp)
cc          vykm = ay(i,j,k )/jb(i,j,k )
cc          vzkp = az(i,j,kp)/jb(i,j,kp)
cc          vzkm = az(i,j,k )/jb(i,j,k )
cc
cc          tsrc = 0.5*(nabla_v_src(i,j,kp)
cc     .               +nabla_v_src(i,j,k ))
cc
cc        case default
cc
cc          vxip = ax(ip,j,k)/jb(ip,j,k)
cc          vxim = ax(im,j,k)/jb(im,j,k)
cc          vyip = ay(ip,j,k)/jb(ip,j,k)
cc          vyim = ay(im,j,k)/jb(im,j,k)
cc          vzip = az(ip,j,k)/jb(ip,j,k)
cc          vzim = az(im,j,k)/jb(im,j,k)
cc                                      
cc          vxjp = ax(i,jp,k)/jb(i,jp,k)
cc          vxjm = ax(i,jm,k)/jb(i,jm,k)
cc          vyjp = ay(i,jp,k)/jb(i,jp,k)
cc          vyjm = ay(i,jm,k)/jb(i,jm,k)
cc          vzjp = az(i,jp,k)/jb(i,jp,k)
cc          vzjm = az(i,jm,k)/jb(i,jm,k)
cc                                      
cc          vxkp = ax(i,j,kp)/jb(i,j,kp)
cc          vxkm = ax(i,j,km)/jb(i,j,km)
cc          vykp = ay(i,j,kp)/jb(i,j,kp)
cc          vykm = ay(i,j,km)/jb(i,j,km)
cc          vzkp = az(i,j,kp)/jb(i,j,kp)
cc          vzkm = az(i,j,km)/jb(i,j,km)
cc
cc          tsrc = nabla_v_src(i,j,k)
cc
cc        end select
cc
cc        idhx = 1d0/dhx
cc        idhy = 1d0/dhy
cc        idhz = 1d0/dhz
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
cc      ! Version I (requires uncommenting code in nabla_v_src)
cc
cccc        tensor = jach*tensor
cccc        if (cntv) tensor = matmul(gsuper,tensor)
cccc
cccc        tensor = tensor + tsrc
cc
cc      ! Version II
cc
cc        tensor = jach*tensor + tsrc
cc
cc        if (cntv) tensor = matmul(gsuper,tensor)
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
ccc     Call variables
cc
cc        integer    :: i,j,k
cc        real(8)    :: tensor(3,3)
cc
ccc     Local variables
cc
cc        integer    :: l,m,n
cc        real(8)    :: hessian(3,3,3),gsuper(3,3)
cc
ccc     Begin program
cc
cc        hessian = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)
cc
cc        !Code for Version I
cccc        gsuper  = gmetric%grid(igrid)%gsup (i,j,k,:,:)
cccc        if (cntv) then
cccc          do m=1,3
cccc            do n=1,3
cccc              do l=1,3
cccc                hessian(l,n,m)=dot_product(gsuper(n,:),hessian(l,:,m))
cccc              enddo
cccc            enddo
cccc          enddo
cccc        endif
cc
cc        do l=1,3
cc          do m=1,3
cc            tensor(l,m) = ax(i,j,k)*hessian(m,l,1)
cc     .                   +ay(i,j,k)*hessian(m,l,2)
cc     .                   +az(i,j,k)*hessian(m,l,3)
cc          enddo
cc        enddo
cc
ccc     End program
cc
cc      end function nabla_v_src
cc
ccccc     nabla_v_src
ccccc     #############################################################
cccc      function nabla_v_src(i,j,k) result(tensor)
cccc
ccccc     -------------------------------------------------------------
ccccc     Finds geometric source of tensor nabla(v) at cell (i,j,k)
ccccc     -------------------------------------------------------------
cccc
cccc        implicit none
cccc
ccccc     Call variables
cccc
cccc        integer    :: i,j,k
cccc        real(8)    :: tensor(3,3)
cccc
ccccc     Local variables
cccc
cccc        integer    :: l,m,vxh,vyh,vzh
cccc        real(8)    :: hessian(3,3,3)
cccc
ccccc     Begin program
cccc
cccc        select case(half_elem)
cccc        case (1)
cccc
cccc          hessian  = 0.5*(gmetric%grid(igrid)%Gamma(i ,j,k,:,:,:)
cccc     .                   +gmetric%grid(igrid)%Gamma(ip,j,k,:,:,:))
cccc
cccc          vxh = 0.5*(ax(ip,j,k)+ax(i,j,k))
cccc          vyh = 0.5*(ay(ip,j,k)+ay(i,j,k))
cccc          vzh = 0.5*(az(ip,j,k)+az(i,j,k))
cccc
cccc        case (2)
cccc
cccc          hessian  = 0.5*(gmetric%grid(igrid)%Gamma(i,j ,k,:,:,:)
cccc     .                   +gmetric%grid(igrid)%Gamma(i,jp,k,:,:,:))
cccc
cccc          vxh = 0.5*(ax(i,jp,k)+ax(i,j,k))
cccc          vyh = 0.5*(ay(i,jp,k)+ay(i,j,k))
cccc          vzh = 0.5*(az(i,jp,k)+az(i,j,k))
cccc
cccc        case (3)
cccc
cccc          hessian  = 0.5*(gmetric%grid(igrid)%Gamma(i,j,k ,:,:,:)
cccc     .                   +gmetric%grid(igrid)%Gamma(i,j,kp,:,:,:))
cccc
cccc          vxh = 0.5*(ax(i,j,kp)+ax(i,j,k))
cccc          vyh = 0.5*(ay(i,j,kp)+ay(i,j,k))
cccc          vzh = 0.5*(az(i,j,kp)+az(i,j,k))
cccc
cccc        case default
cccc
cccc          hessian  = gmetric%grid(igrid)%Gamma(i,j,k ,:,:,:)
cccc
cccc          vxh = ax(i,j,k)
cccc          vyh = ay(i,j,k)
cccc          vzh = az(i,j,k)
cccc
cccc        end select
cccc
cccc        do l=1,3
cccc          do m=1,3
cccc            tensor(l,m) = vxh*hessian(m,1,l)
cccc     .                   +vyh*hessian(m,2,l)
cccc     .                   +vzh*hessian(m,3,l)
cccc          enddo
cccc        enddo
cccc
ccccc     End program
cccc
cccc      end function nabla_v_src
cc
cc      end function fnabla_v

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

        integer    :: i,j,k,hex,hey,hez,nx,ny,nz,igx,igy,igz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid
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

c     fnabla_v_cov
c     #############################################################
      function fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz,ax,ay,az
     .                     ,half_elem) result(tensor)
c     -------------------------------------------------------------
c     Calculates the tensor components T_lm of nabla(vec v), with
c     v covariant, and fills tensor(l,m) at the following positions:
c       + half_elem=0 --> i,j,k
c       + half_elem=1 --> i+1/2,j,k
c       + half_elem=2 --> i,j+1/2,k
c       + half_elem=3 --> i,j,k+1/2
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,half_elem,nx,ny,nz,igx,igy,igz
        real(8)    :: tensor(3,3)
        real(8)    :: ax(0:nx+1,0:ny+1,0:nz+1)
     .               ,ay(0:nx+1,0:ny+1,0:nz+1)
     .               ,az(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igrid,jpp,jmm
        real(8)    :: dhx,dhy,dhz,idhx,idhy,idhz
        real(8)    :: vxx,vyy,vzz
     .               ,vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
        real(8)    :: tsrc(3,3),jacip,jacjp,jacjm,jackp,jackm
     .               ,jacipkp,jacipkm,jacipjp,jacipjm,jac
        real(8)    :: cf1,cf2,cf3,cf4

        logical    :: spoint

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

        spoint = bcSP()
cc        spoint = isSP(i,j,k,igx,igy,igz)
cc        spoint = .false.

        !Exceptions

        !!Coeffs for 4th order first derivative
        cf1 = 1./8.
        cf2 = 9./8.

        !!Coeffs for 4th order second derivative
        cf3 =  1./12.
        cf4 = 15./12.

        select case(half_elem)
        case (1)
          dhx = dx(ig)

          vxip = ax(ip,j,k)
          vxim = ax(i ,j,k)
          vyip = ay(ip,j,k)
          vyim = ay(i ,j,k)
          vzip = az(ip,j,k)
          vzim = az(i ,j,k)

          if (spoint) then
c FIX PARALLEL for PER BC
            !4th order
            jpp = jp+1
            jmm = jm-1
            if (j == ny) jpp = 2
            if (j == 1 ) jmm = ny-1

            vxjp = 0.5*(-cf1*(ax(ip,jpp,k)+ax(i,jpp,k))
     .                  +cf2*(ax(ip,jp ,k)+ax(i,jp ,k))
     .                  +cf2*(ax(ip,j  ,k)+ax(i,j  ,k))
     .                  -cf1*(ax(ip,jm ,k)+ax(i,jm ,k)))
            vxjm = 0.5*(-cf1*(ax(ip,jp ,k)+ax(i,jp ,k))
     .                  +cf2*(ax(ip,j  ,k)+ax(i,j  ,k))
     .                  +cf2*(ax(ip,jm ,k)+ax(i,jm ,k))
     .                  -cf1*(ax(ip,jmm,k)+ax(i,jmm,k)))

            vyjp = 0.5*(-cf1*(ay(ip,jpp,k)+ay(i,jpp,k))
     .                  +cf2*(ay(ip,jp ,k)+ay(i,jp ,k))
     .                  +cf2*(ay(ip,j  ,k)+ay(i,j  ,k))
     .                  -cf1*(ay(ip,jm ,k)+ay(i,jm ,k)))
            vyjm = 0.5*(-cf1*(ay(ip,jp ,k)+ay(i,jp ,k))
     .                  +cf2*(ay(ip,j  ,k)+ay(i,j  ,k))
     .                  +cf2*(ay(ip,jm ,k)+ay(i,jm ,k))
     .                  -cf1*(ay(ip,jmm,k)+ay(i,jmm,k)))

            vzjp = 0.5*(-cf1*(az(ip,jpp,k)+az(i,jpp,k))
     .                  +cf2*(az(ip,jp ,k)+az(i,jp ,k))
     .                  +cf2*(az(ip,j  ,k)+az(i,j  ,k))
     .                  -cf1*(az(ip,jm ,k)+az(i,jm ,k)))
            vzjm = 0.5*(-cf1*(az(ip,jp ,k)+az(i,jp ,k))
     .                  +cf2*(az(ip,j  ,k)+az(i,j  ,k))
     .                  +cf2*(az(ip,jm ,k)+az(i,jm ,k))
     .                  -cf1*(az(ip,jmm,k)+az(i,jmm,k)))
          else
          !2nd order
            vxjp = 0.5*(ax(ip,jp,k)+ax(i,jp,k))
            vxjm = 0.5*(ax(ip,jm,k)+ax(i,jm,k))
            vyjp = 0.5*(ay(ip,jp,k)+ay(i,jp,k))
            vyjm = 0.5*(ay(ip,jm,k)+ay(i,jm,k))
            vzjp = 0.5*(az(ip,jp,k)+az(i,jp,k))
            vzjm = 0.5*(az(ip,jm,k)+az(i,jm,k))
          endif

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

          if (spoint) then
            !4th order
c FIX PARALLEL for PER BC
            jpp = jp+1
            if (j == ny) jpp = 2
            if (j == 0 ) jm  = ny-1

            vxjp = cf4*ax(i,jp,k)-cf3*ax(i,jpp,k)
            vxjm = cf4*ax(i,j ,k)-cf3*ax(i,jm ,k)
            vyjp = cf4*ay(i,jp,k)-cf3*ay(i,jpp,k)
            vyjm = cf4*ay(i,j ,k)-cf3*ay(i,jm ,k)
            vzjp = cf4*az(i,jp,k)-cf3*az(i,jpp,k)
            vzjm = cf4*az(i,j ,k)-cf3*az(i,jm ,k)

            tsrc = 0.5*(-cf1*nabla_v_src(i,jpp,k)
     .                  +cf2*nabla_v_src(i,jp ,k)
     .                  +cf2*nabla_v_src(i,j  ,k)
     .                  -cf1*nabla_v_src(i,jm ,k))
          else
            !2nd order
            vxjp = ax(i,jp,k)
            vxjm = ax(i,j ,k)
            vyjp = ay(i,jp,k)
            vyjm = ay(i,j ,k)
            vzjp = az(i,jp,k)
            vzjm = az(i,j ,k)

            tsrc = 0.5*(nabla_v_src(i,jp,k)
     .                 +nabla_v_src(i,j ,k))
          endif

          vxkp = 0.5*(ax(i,j,kp)+ax(i,jp,kp))
          vxkm = 0.5*(ax(i,j,km)+ax(i,jp,km))
          vykp = 0.5*(ay(i,j,kp)+ay(i,jp,kp))
          vykm = 0.5*(ay(i,j,km)+ay(i,jp,km))
          vzkp = 0.5*(az(i,j,kp)+az(i,jp,kp))
          vzkm = 0.5*(az(i,j,km)+az(i,jp,km))

        case (3)
          dhz = dz(kg)

          vxip = 0.5*(ax(ip,j,k)+ax(ip,j,kp))
          vxim = 0.5*(ax(im,j,k)+ax(im,j,kp))
          vyip = 0.5*(ay(ip,j,k)+ay(ip,j,kp))
          vyim = 0.5*(ay(im,j,k)+ay(im,j,kp))
          vzip = 0.5*(az(ip,j,k)+az(ip,j,kp))
          vzim = 0.5*(az(im,j,k)+az(im,j,kp))

          if (spoint) then
c FIX PARALLEL for PER BC
          !4th order
            jpp = jp+1
            jmm = jm-1
            if (j == ny) jpp = 2
            if (j == 1 ) jmm = ny-1

            vxjp = 0.5*(-cf1*(ax(i,jpp,kp)+ax(i,jpp,k))
     .                  +cf2*(ax(i,jp ,kp)+ax(i,jp ,k))
     .                  +cf2*(ax(i,j  ,kp)+ax(i,j  ,k))
     .                  -cf1*(ax(i,jm ,kp)+ax(i,jm ,k)))
            vxjm = 0.5*(-cf1*(ax(i,jp ,kp)+ax(i,jp ,k))
     .                  +cf2*(ax(i,j  ,kp)+ax(i,j  ,k))
     .                  +cf2*(ax(i,jm ,kp)+ax(i,jm ,k))
     .                  -cf1*(ax(i,jmm,kp)+ax(i,jmm,k)))

            vyjp = 0.5*(-cf1*(ay(i,jpp,kp)+ay(i,jpp,k))
     .                  +cf2*(ay(i,jp ,kp)+ay(i,jp ,k))
     .                  +cf2*(ay(i,j  ,kp)+ay(i,j  ,k))
     .                  -cf1*(ay(i,jm ,kp)+ay(i,jm ,k)))
            vyjm = 0.5*(-cf1*(ay(i,jp ,kp)+ay(i,jp ,k))
     .                  +cf2*(ay(i,j  ,kp)+ay(i,j  ,k))
     .                  +cf2*(ay(i,jm ,kp)+ay(i,jm ,k))
     .                  -cf1*(ay(i,jmm,kp)+ay(i,jmm,k)))

            vzjp = 0.5*(-cf1*(az(i,jpp,kp)+az(i,jpp,k))
     .                  +cf2*(az(i,jp ,kp)+az(i,jp ,k))
     .                  +cf2*(az(i,j  ,kp)+az(i,j  ,k))
     .                  -cf1*(az(i,jm ,kp)+az(i,jm ,k)))
            vzjm = 0.5*(-cf1*(az(i,jp ,kp)+az(i,jp ,k))
     .                  +cf2*(az(i,j  ,kp)+az(i,j  ,k))
     .                  +cf2*(az(i,jm ,kp)+az(i,jm ,k))
     .                  -cf1*(az(i,jmm,kp)+az(i,jmm,k)))
          else
          !2nd order
            vxjp = 0.5*(ax(i,jp,k)+ax(i,jp,kp))
            vxjm = 0.5*(ax(i,jm,k)+ax(i,jm,kp))
            vyjp = 0.5*(ay(i,jp,k)+ay(i,jp,kp))
            vyjm = 0.5*(ay(i,jm,k)+ay(i,jm,kp))
            vzjp = 0.5*(az(i,jp,k)+az(i,jp,kp))
            vzjm = 0.5*(az(i,jm,k)+az(i,jm,kp))
          endif

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

          if (spoint) then
c FIX PARALLEL for PER BC
          !4th order
            jpp = jp+1
            jmm = jm-1
            if (j == ny) jpp = 2
            if (j == 1 ) jmm = ny-1

            vxjp = (-cf1*(ax(i,jpp,k))
     .              +cf2*(ax(i,jp ,k))
     .              +cf2*(ax(i,j  ,k))
     .              -cf1*(ax(i,jm ,k)))
            vxjm = (-cf1*(ax(i,jp ,k))
     .              +cf2*(ax(i,j  ,k))
     .              +cf2*(ax(i,jm ,k))
     .              -cf1*(ax(i,jmm,k)))

            vyjp = (-cf1*(ay(i,jpp,k))
     .              +cf2*(ay(i,jp ,k))
     .              +cf2*(ay(i,j  ,k))
     .              -cf1*(ay(i,jm ,k)))
            vyjm = (-cf1*(ay(i,jp ,k))
     .              +cf2*(ay(i,j  ,k))
     .              +cf2*(ay(i,jm ,k))
     .              -cf1*(ay(i,jmm,k)))

            vzjp = (-cf1*(az(i,jpp,k))
     .              +cf2*(az(i,jp ,k))
     .              +cf2*(az(i,j  ,k))
     .              -cf1*(az(i,jm ,k)))
            vzjm = (-cf1*(az(i,jp ,k))
     .              +cf2*(az(i,j  ,k))
     .              +cf2*(az(i,jm ,k))
     .              -cf1*(az(i,jmm,k)))
          else
          !2nd order
            vxjp = ax(i,jp,k)
            vxjm = ax(i,jm,k)
            vyjp = ay(i,jp,k)
            vyjm = ay(i,jm,k)
            vzjp = az(i,jp,k)
            vzjm = az(i,jm,k)
          endif

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

        integer    :: i,j,k
        real(8)    :: tensor(3,3)

        real(8)    :: hessian(3,3,3)

c     Begin program

        hessian = gmetric%grid(igrid)%Gamma(i,j,k,:,:,:)

      ! l = 1, m = 1
        tensor(1,1) =+ ax(i,j,k)* hessian(1,1,1)
     .               + ay(i,j,k)* hessian(2,1,1)
     .               + az(i,j,k)* hessian(3,1,1)

      ! l = 1, m = 2
        tensor(1,2) =+ ax(i,j,k)* hessian(1,1,2)
     .               + ay(i,j,k)* hessian(2,1,2)
     .               + az(i,j,k)* hessian(3,1,2)

      ! l = 1, m = 3
        tensor(1,3) =+ ax(i,j,k)* hessian(1,1,3)
     .               + ay(i,j,k)* hessian(2,1,3)
     .               + az(i,j,k)* hessian(3,1,3)

      ! l = 2, m = 1
        tensor(2,1) =+ ax(i,j,k)* hessian(1,2,1)
     .               + ay(i,j,k)* hessian(2,2,1)
     .               + az(i,j,k)* hessian(3,2,1)

      ! l = 2, m = 2
        tensor(2,2) =+ ax(i,j,k)* hessian(1,2,2)
     .               + ay(i,j,k)* hessian(2,2,2)
     .               + az(i,j,k)* hessian(3,2,2)

      ! l = 2, m = 3
        tensor(2,3) =+ ax(i,j,k)* hessian(1,2,3)
     .               + ay(i,j,k)* hessian(2,2,3)
     .               + az(i,j,k)* hessian(3,2,3)

      ! l = 3, m = 1
        tensor(3,1) =+ ax(i,j,k)* hessian(1,3,1)
     .               + ay(i,j,k)* hessian(2,3,1)
     .               + az(i,j,k)* hessian(3,3,1)

      ! l = 3, m = 2
        tensor(3,2) =+ ax(i,j,k)* hessian(1,3,2)
     .               + ay(i,j,k)* hessian(2,3,2)
     .               + az(i,j,k)* hessian(3,3,2)

      ! l = 3, m = 3
        tensor(3,3) =+ ax(i,j,k)* hessian(1,3,3)
     .               + ay(i,j,k)* hessian(2,3,3)
     .               + az(i,j,k)* hessian(3,3,3)

c     End program

      end function nabla_v_src

      end function fnabla_v_cov

c     btensor_x
c     #############################################################
      subroutine btensor_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp

c     Begin program

        ip = i+1
        if (flag == 0) ip = i
cc        if (flag == 0) then
cc          t11 = 0d0
cc          t12 = 0d0
cc          t13 = 0d0
cc          return
cc        endif

        jacp = gmetric%grid(igx)%jac(ip,j,k)
        jac0 = gmetric%grid(igx)%jac(i ,j,k)
        jac  = 0.5*(jacp+jac0)

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = 0d0

        ijac  = 1d0/jac
        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

cc        if (isSP(i+1,j,k,igx,igy,igz)) ijac = 1d0/SP_flsv

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

cc        if (isSP(i+1,j,k,igx,igy,igz)) then
cc          write (*,*) t12,t13
cc        endif

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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp

c     Begin program

        jp = j+1
        if (flag == 0) jp = j
cc        if (flag == 0) then
cc          t21 = 0d0
cc          t22 = 0d0
cc          t23 = 0d0
cc          return
cc        endif

        jacp = gmetric%grid(igx)%jac(i,jp,k)
        jac0 = gmetric%grid(igx)%jac(i,j ,k)
        jac  = 0.5*(jacp+jac0)
cc        jac  = 1d0

        ijac  = 1d0/jac
        ijacp = 1d0/jacp
        ijac0 = 1d0/jac0

        t22 = 0d0

        if (solenoidal) then
          t21 = 0.5*( vec1(i,j ,k,2)*vec2(i,j ,k,1)*ijac0
     .               +vec1(i,jp,k,2)*vec2(i,jp,k,1)*ijacp)*jac
     .         -0.5*( vec1(i,j ,k,1)*vec2(i,j ,k,2)*ijac0
     .               +vec1(i,jp,k,1)*vec2(i,jp,k,2)*ijacp)*jac

          t23 = 0.5*( vec1(i,j ,k,2)*vec2(i,j ,k,3)*ijac0
     .               +vec1(i,jp,k,2)*vec2(i,jp,k,3)*ijacp)*jac
     .         -0.5*( vec1(i,j ,k,3)*vec2(i,j ,k,2)*ijac0
     .               +vec1(i,jp,k,3)*vec2(i,jp,k,2)*ijacp)*jac
        else
          t21 = 0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,1)*ijac0
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,1)*ijacp)*jac
     .         -0.5*( vec1(i,jp,k,1)*vec2(i,j ,k,2)*ijacp
     .               +vec1(i,j ,k,1)*vec2(i,jp,k,2)*ijac0)*jac

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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp,ijac,ijacp,ijac0

c     Begin program

        kp = k+1
        if (flag == 0) kp = k
cc        if (flag == 0) then
cc          t31 = 0d0
cc          t32 = 0d0
cc          t33 = 0d0
cc          return
cc        endif

        jacp = gmetric%grid(igx)%jac(i,j,kp)
        jac0 = gmetric%grid(igx)%jac(i,j,k )
        jac  = 0.5*(jacp+jac0)

        ijac  = 1d0/jac
        ijacp = 1d0/jacp
        ijac0 = 1d0/jac0

        t33 = 0d0

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

        if (flag /= 0) then
          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine btensor_z

c     curl_inv
c     ###############################################################
      subroutine curl_inv(nx,ny,nz,igx,igy,igz,b,a,global)

c     ---------------------------------------------------------------
c     Finds COVARIANT vector potential from CONTRAVARIANT vector
c     field components. Employs gauge A1 = 0d0. In SP systems,
c     it integrates at faces in r, and then averages to nodes.
c     
c     Integrals done by Numerical Recipes p. 128, 4.1.12 -- actually
c     one bit better -- extrapolates quadratically to f(-1) and uses
c     the internal formula, rather than using midpoint for the first
c     point. Thus even the first point has third order accuracy.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      integer    :: nx,ny,nz,igx,igy,igz
      real(8)    :: a(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,b(0:nx+1,0:ny+1,0:nz+1,3)
      logical,optional :: global

c     Local variables

      integer    :: i,j,k,igrid,ii,bcnd(6,3),ig,jg,kg,ivar
      real(8)    :: cm,c0,cp,dxx,dyy,dvol,cov(3)

      logical    :: spoint,glbl

c     Begin program

      if (PRESENT(global)) then
        glbl = global
      else
        glbl = .false.
      endif

      igrid = igx

c     Consistency check

      if (igrid /= 1 .and. glbl) then
        call pstop('curl_inv','Inconsistent grid level')
      endif

c     New code (SP ready, parallel-ready, version I)

      a(:,:,:,1)=0d0  !Gauge

      a(0,0,:,3)=0d0

      i = 0

cc      spoint = isSP(1,1,1,igx,igy,igz)  !Check whether current domain contains SP
      spoint = .false.

      do j=1,ny+1
        if (.not.glbl) then
          call getMGmap(i,j,1,igx,igy,igz,ig,jg,kg)
          dyy = grid_params%dy(jg-1)
        else
          dyy = grid_params%yg(j)-grid_params%yg(j-1)
        endif

        if (spoint) then  !r=0 face
          a(i,j,:,3)=a(i,j-1,:,3) + dyy*2.5d-1*(b(i  ,j-1,:,1)
     .                                         +b(i  ,j  ,:,1)
     .                                         +b(i+1,j-1,:,1)
     .                                         +b(i+1,j  ,:,1))
        else
          a(i,j,:,3)=a(i,j-1,:,3) + dyy*  5d-1*(b(i,j-1,:,1)
     .                                         +b(i,j  ,:,1))
        endif
      enddo

      a(0,:,:,2)=0d0

      do i=1,nx+1
        if (.not.glbl) call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)

        if (spoint) then
          if (.not.glbl) then
            dxx = grid_params%dxh(ig-1)
          else
cc            dxx = 0.5*(grid_params%xg(i+1)-grid_params%xg(i-1))
            dxx = grid_params%xg(i)-grid_params%xg(i-1)
          endif

          a(i,:,:,3) = a(i-1,:,:,3) - dxx*b(i,:,:,2)

          if (isSP(i,1,1,igx,igy,igz)) then
            a(i,:,:,2) = a(i-1,:,:,2) + 5d-1*dxx*b(i,:,:,3)  !Factor of 1/2 due to geometry
          else
            a(i,:,:,2) = a(i-1,:,:,2) +      dxx*b(i,:,:,3)
          endif

        else

          if (.not.glbl) then
            dxx = grid_params%dx(ig-1)
          else
            dxx = grid_params%xg(i)-grid_params%xg(i-1)
          endif

          a(i,:,:,3) = a(i-1,:,:,3) - dxx*5d-1*(b(i-1,:,:,2)
     .                                         +b(i  ,:,:,2))
          a(i,:,:,2) = a(i-1,:,:,2) + dxx*5d-1*(b(i-1,:,:,3)
     .                                         +b(i  ,:,:,3))

        endif
      enddo

      !Average from radial faces to nodes in SP coordinate systems
      if (spoint) a(1:nx+1,:,:,2:3) = 0.5*(a(1:nx+1,:,:,2:3)
     .                                    +a(0:nx  ,:,:,2:3))

      return

c     New code (SP-ready, NOT parallel-ready, version II)

cc      a(:,:,:,1)=0d0  !Gauge
cc
cc      if (bcSP()) then
cc
cc        a(0,0,:,3)=0d0
cc
cc        i = 0
cc
cc        do k=0,nz+1
cc          do j=1,ny+1
cc            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc            dyy = grid_params%dy(jg-1)
cc
cc            a(i,j,k,3)=a(i,j-1,k,3) + dyy*0.25*(b(i  ,j-1,k,1)
cc     .                                         +b(i  ,j  ,k,1)
cc     .                                         +b(i+1,j-1,k,1)
cc     .                                         +b(i+1,j  ,k,1))
cc          enddo
cc        enddo
cc
cc        a(0,:,:,2)=0d0
cc
cc        do k=0,nz+1
cc          do j=0,ny+1
cc            do i=1,nx+1
cc
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              dxx = grid_params%dxh(ig-1)
cc
cc              a(i,j,k,3)=a(i-1,j,k,3) -dxx*b(i,j,k,2)
cc
cc              if (isSP(i,j,k,igx,igy,igz)) then
cc                a(i,j,k,2)=a(i-1,j,k,2) +0.5*dxx*b(i,j,k,3) !Factor of 1/2 due to geometry
cc              else
cc                a(i,j,k,2)=a(i-1,j,k,2) +    dxx*b(i,j,k,3)
cc              endif
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc        !Average from radial faces to centers
cc        a(1:nx+1,:,:,2:3) = 0.5*(a(1:nx+1,:,:,2:3)+a(0:nx,:,:,2:3))
cc
cc        !Boundary conditions: topological constraints (SP, PER)
cc
cc        bcnd(:,1) = bcond
cc        bcnd(:,2) = bcond
cc        bcnd(:,3) = bcond
cc
cc        call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IBX
cc     .            ,is_vec=.true.,is_cnv=.false.,iorder=2)
cc
cc      else
cc
cc        a(0,0,:,3)=0d0
cc
cc        i = 0
cc
cc        do k=0,nz+1
cc          do j=1,ny+1
cc            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc            dyy = grid_params%dy(jg-1)
cc
cc            a(i,j,k,3)=a(i,j-1,k,3) + dyy*0.5*(b(i,j-1,k,1)
cc     .                                        +b(i,j  ,k,1))
cc
cc          enddo
cc        enddo
cc
cc        a(0,:,:,2)=0d0
cc
cc        do k=0,nz+1
cc          do j=0,ny+1
cc            do i=1,nx+1
cc
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              dxx = grid_params%dx(ig-1)
cc
cc              a(i,j,k,3)=a(i-1,j,k,3) -dxx*0.5*(b(i-1,j,k,2)
cc     .                                         +b(i  ,j,k,2))
cc              a(i,j,k,2)=a(i-1,j,k,2) +dxx*0.5*(b(i-1,j,k,3)
cc     .                                         +b(i  ,j,k,3))
cc
cc            enddo
cc          enddo
cc        enddo
cc
cc        bcnd(:,1) = bcond
cc        bcnd(:,2) = bcond
cc        bcnd(:,3) = bcond
cc        
cc        call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IBX
cc     .              ,is_vec=.true.,is_cnv=.false.,iorder=2)

c     JOHN FINN's code (assumes uniform mesh)

cc        cm=-0.0833333333333
cc        c0= 0.6666666666667
cc        cp= 0.4166666666667
cc
cc        a(1,1,:,3) = 0d0
cc
cc        do k=1,nz
cc          do j=2,ny
cc            call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,y1,z1)
cc
cc            dyy = grid_params%dyh(jg)
cc
cccc            if(j.eq.2) then
cccc              a(1,j,k,3)=a(1,j-1,k,3)
cccc     $          +dyy*(cp*b(1,j-1,k,1)+c0*b(1,j  ,k,1)+cm*b(1,j+1,k,1))
cccc            else
cc              a(1,j,k,3)=a(1,j-1,k,3)
cc     .          +dyy*(cm*b(1,j-2,k,1)+c0*b(1,j-1,k,1)+cp*b(1,j  ,k,1))
cccc            endif
cc          enddo
cc        enddo
cc
cc        a(1,:,:,2) = 0d0
cc
cc        do k=1,nz
cc          do j=1,ny
cc            do i=2,nx
cc
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz
cc     .                                    ,ig,jg,kg,x1,y1,z1)
cc
cc              dxx = grid_params%dxh(ig)
cc
cccc              if(i.eq.2) then
cccc                a(i,j,k,3)=a(i-1,j,k,3)
cccc     .            -dxx*(cp*b(i-1,j,k,2)+c0*b(i,j,k,2)+cm*b(i+1,j,k,2))
cccc                a(i,j,k,2)=a(i-1,j,k,2)
cccc     .            +dxx*(cp*b(i-1,j,k,3)+c0*b(i,j,k,3)+cm*b(i+1,j,k,3))
cccc              else
cc                a(i,j,k,3)=a(i-1,j,k,3)
cc     .            -dxx*(cm*b(i-2,j,k,2)+c0*b(i-1,j,k,2)+cp*b(i,j,k,2))
cc                a(i,j,k,2)=a(i-1,j,k,2)
cc     .            +dxx*(cm*b(i-2,j,k,3)+c0*b(i-1,j,k,3)+cp*b(i,j,k,3))
cccc              endif
cc            enddo
cc          enddo
cc        enddo
cc
ccc     Boundary conditions: topological constraints (PER) and extrapolation
cc
cc        bcnd(:,1) = bcond
cc        bcnd(:,2) = bcond
cc        bcnd(:,3) = bcond
cc        where (bcnd == DEF) bcnd = EXT
cc        
cc        call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IBX
cc     .              ,is_vec=.true.,is_cnv=.false.,iorder=2)
cc
cc      endif


c     End program

      end subroutine curl_inv

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
cc        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
cc        real(8)    :: t11,t12,t13
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ip,igrid
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
cc        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
cc        real(8)    :: t21,t22,t23
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: jp,igrid
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
cc        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
cc        real(8)    :: t31,t32,t33
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: kp,igrid
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

c     eom_f
c     ###################################################################
      function eom_f(i,j,k,nx,ny,nz,igx,igy,igz,rr,tt,bb,jj,vec_cnv
     $              ,nc_eom_jxb,nc_eom_gp,alt_eom,vol)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of the EOM force term
c     (jxb)-grad(rho*tmp) at cell centers. Vectors are taken as
c     contravariant if vec_cnv=.true., and covariant otherwise
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nx,ny,nz,igx,igy,igz
        real(8)    :: cnv(3)
        real(8),target :: jj(0:nx+1,0:ny+1,0:nz+1,3)
     $                   ,bb(0:nx+1,0:ny+1,0:nz+1,3)
     $                   ,rr(0:nx+1,0:ny+1,0:nz+1)
     $                   ,tt(0:nx+1,0:ny+1,0:nz+1)
        logical    :: vec_cnv,nc_eom_jxb,nc_eom_gp,alt_eom
        logical,optional,intent(IN) :: vol

c     Local variables

        integer    :: igrid,ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: cov(3),jac,cnv1(3),dxx,dyy,dzz
        logical    :: vol_wgt

c     Begin program

        igrid = igx

        vol_wgt = .true.
        if (PRESENT(vol)) vol_wgt = vol

        !Lorentz force
        if (nc_eom_jxb) then
          cnv(1) = jj(i,j,k,2)*bb(i,j,k,3)
     .           - jj(i,j,k,3)*bb(i,j,k,2)

          cnv(2) = jj(i,j,k,3)*bb(i,j,k,1)
     .           - jj(i,j,k,1)*bb(i,j,k,3)

          cnv(3) = jj(i,j,k,1)*bb(i,j,k,2)
     .           - jj(i,j,k,2)*bb(i,j,k,1)

          if (vec_cnv) then
            jac = gmetric%grid(igrid)%jac(i,j,k)

            cov = cnv/jac

            call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                                  ,cov(1),cov(2),cov(3)
     .                                  ,cnv(1),cnv(2),cnv(3),.true.)
          endif

          if (vol_wgt) cnv = cnv*gmetric%grid(igrid)%dvol(i,j,k)

        else
          if (.not.vec_cnv)
     .      call pstop('eom_f','Magnetic field should be contravariant')

          vec1 => bb
          cnv =-div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,eom_jxb_x,eom_jxb_y,eom_jxb_z,vol=vol_wgt)
          nullify(vec1)
        endif

        !Pressure force
        cnv1 = cnv
        if (nc_eom_gp) then
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1

          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

          dxx = dxh(ig)
          dyy = dyh(jg)
          dzz = dzh(kg)

          cov(1) = 0.5*( rr(i,j,k)*(tt(ip,j,k)-tt(im,j,k))/dxx
     .                  +tt(i,j,k)*(rr(ip,j,k)-rr(im,j,k))/dxx)

          cov(2) = 0.5*( rr(i,j,k)*(tt(i,jp,k)-tt(i,jm,k))/dyy
     .                  +tt(i,j,k)*(rr(i,jp,k)-rr(i,jm,k))/dyy)
          
          cov(3) = 0.5*( rr(i,j,k)*(tt(i,j,kp)-tt(i,j,km))/dzz
     .                  +tt(i,j,k)*(rr(i,j,kp)-rr(i,j,km))/dzz)

cc          cov(1)=0.5*(rr(ip,j,k)*tt(ip,j,k)-rr(im,j,k)*tt(im,j,k))/dxx
cc          cov(2)=0.5*(rr(i,jp,k)*tt(i,jp,k)-rr(i,jm,k)*tt(i,jm,k))/dyy
cc          cov(3)=0.5*(rr(i,j,kp)*tt(i,j,kp)-rr(i,j,km)*tt(i,j,km))/dzz

          cov = cov*a_p    !Multiply by alpha_p=1 + Ti/Te

          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                                ,cov(1),cov(2),cov(3)
     .                                ,cnv(1),cnv(2),cnv(3),.true.)

          if (vol_wgt) cnv = cnv*gmetric%grid(igrid)%dvol(i,j,k)

        else
          sc1 => rr
          sc2 => tt
          cnv = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,eom_grad_p_x,eom_grad_p_y,eom_grad_p_z
     .                    ,vol=vol_wgt)
          nullify(sc1,sc2)
        endif

        cnv = cnv1 - cnv

c     End program

      end function eom_f

c     eom_jxb_x
c     #############################################################
      subroutine eom_jxb_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for jxb force term in EOM.
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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip
        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3)
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,pmag

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

        cnv0 = vec1(i,j,k,:)
        call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                              ,cov0(1),cov0(2),cov0(3)
     .                              ,cnv0(1),cnv0(2),cnv0(3),.false.)
        cnvp = vec1(ip,j,k,:)
        call transformFromCurvToCurv(ip,j,k,igx,igy,igz
     .                              ,covp(1),covp(2),covp(3)
     .                              ,cnvp(1),cnvp(2),cnvp(3),.false.)

        pmag = jac*(cnvp(1)*cov0(1)*ijacp
     .             +cnv0(1)*covp(1)*ijac0
     .             +cnvp(2)*cov0(2)*ijac0
     .             +cnv0(2)*covp(2)*ijacp
     .             +cnvp(3)*cov0(3)*ijacp
     .             +cnv0(3)*covp(3)*ijac0)*0.25

cc        pmag = jac*(bx(ip,j,k)*bx_cov(i ,j,k)*ijacp
cc     .             +bx(i ,j,k)*bx_cov(ip,j,k)*ijac0
cc     .             +by(ip,j,k)*by_cov(i ,j,k)*ijac0
cc     .             +by(i ,j,k)*by_cov(ip,j,k)*ijacp
cc     .             +bz(ip,j,k)*bz_cov(i ,j,k)*ijacp
cc     .             +bz(i ,j,k)*bz_cov(ip,j,k)*ijac0)*0.25
cccc     .        +(bx(ip,j,k)*bx_cov(i ,j,k)
cccc     .         +bx(i ,j,k)*bx_cov(ip,j,k)
cccc     .         +by(ip,j,k)*by_cov(i ,j,k)
cccc     .         +by(i ,j,k)*by_cov(ip,j,k)
cccc     .         +bz(ip,j,k)*bz_cov(i ,j,k)
cccc     .         +bz(i ,j,k)*bz_cov(ip,j,k))*0.25

        t11 = -( vec1(ip,j,k,1)*vec1(i ,j,k,1)*ijacp*ijac0 )*jac**2
     .        +gsuper(1,1)*pmag

        t12 = -0.5*( vec1(ip,j,k,1)*vec1(i ,j,k,2)*ijacp
     .              +vec1(i ,j,k,1)*vec1(ip,j,k,2)*ijac0 )*jac
     .        +gsuper(1,2)*pmag

        t13 = -0.5*( vec1(ip,j,k,1)*vec1(i ,j,k,3)*ijacp*ijac0
     .              +vec1(i ,j,k,1)*vec1(ip,j,k,3)*ijacp*ijac0)*jac**2
     .        +gsuper(1,3)*pmag

        if (flag /= 0) then
          ijac = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine eom_jxb_x

c     eom_jxb_y
c     #############################################################
      subroutine eom_jxb_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for jxb force term in EOM.
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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp
        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3)
        real(8)    :: jac,ijac,pmag

c     Begin program

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))

        cnv0 = vec1(i,j,k,:)
        call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                              ,cov0(1),cov0(2),cov0(3)
     .                              ,cnv0(1),cnv0(2),cnv0(3),.false.)
        cnvp = vec1(i,jp,k,:)
        call transformFromCurvToCurv(i,jp,k,igx,igy,igz
     .                              ,covp(1),covp(2),covp(3)
     .                              ,cnvp(1),cnvp(2),cnvp(3),.false.)

        pmag = (cnvp(1)*cov0(1)
     .         +cnv0(1)*covp(1)
     .         +cnvp(2)*cov0(2)
     .         +cnv0(2)*covp(2)
     .         +cnvp(3)*cov0(3)
     .         +cnv0(3)*covp(3))*0.25

        t21 = -0.5*( vec1(i,jp,k,2)*vec1(i,j ,k,1)
     .              +vec1(i,j ,k,2)*vec1(i,jp,k,1) )
     .        +gsuper(2,1)*pmag

        t22 = -0.5*( vec1(i,jp,k,2)*vec1(i,j ,k,2)
     .              +vec1(i,j ,k,2)*vec1(i,jp,k,2) )
     .        +gsuper(2,2)*pmag

        t23 = -0.5*( vec1(i,jp,k,2)*vec1(i,j ,k,3)
     .              +vec1(i,j ,k,2)*vec1(i,jp,k,3) )
     .        +gsuper(2,3)*pmag

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine eom_jxb_y

c     eom_jxb_z
c     #############################################################
      subroutine eom_jxb_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for jxb force term in EOM.
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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp
        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3)
        real(8)    :: jac,ijac,pmag

c     Begin program

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))

        cnv0 = vec1(i,j,k,:)
        call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                              ,cov0(1),cov0(2),cov0(3)
     .                              ,cnv0(1),cnv0(2),cnv0(3),.false.)
        cnvp = vec1(i,j,kp,:)
        call transformFromCurvToCurv(i,j,kp,igx,igy,igz
     .                              ,covp(1),covp(2),covp(3)
     .                              ,cnvp(1),cnvp(2),cnvp(3),.false.)

        pmag = (cnvp(1)*cov0(1)
     .         +cnv0(1)*covp(1)
     .         +cnvp(2)*cov0(2)
     .         +cnv0(2)*covp(2)
     .         +cnvp(3)*cov0(3)
     .         +cnv0(3)*covp(3))*0.25

        t31 = -0.5*( vec1(i,j,kp,3)*vec1(i,j,k ,1)
     .              +vec1(i,j,k ,3)*vec1(i,j,kp,1) )
     .        +gsuper(3,1)*pmag

        t32 = -0.5*( vec1(i,j,kp,3)*vec1(i,j,k ,2)
     .              +vec1(i,j,k ,3)*vec1(i,j,kp,2) )
     .        +gsuper(3,2)*pmag

        t33 = -0.5*( vec1(i,j,kp,3)*vec1(i,j,k ,3)
     .              +vec1(i,j,k ,3)*vec1(i,j,kp,3) )
     .        +gsuper(3,3)*pmag

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine eom_jxb_z

c     eom_grad_p_x
c     #############################################################
      subroutine eom_grad_p_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip
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

       !Recall p=a_p nT
        ptot = a_p*0.5*jac*(sc1(ip,j,k)*sc2(i ,j,k)
     .                     +sc1(i ,j,k)*sc2(ip,j,k))

        t11 = gsuper(1,1)*ptot

        t12 = gsuper(1,2)*ptot

        t13 = gsuper(1,3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine eom_grad_p_x

c     eom_grad_p_y
c     #############################################################
      subroutine eom_grad_p_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp
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
        ptot = a_p*0.5*jac*(sc1(i,jp,k)*sc2(i,j ,k)
     .                     +sc1(i,j ,k)*sc2(i,jp,k))

        t21 = gsuper(2,1)*ptot

        t22 = gsuper(2,2)*ptot

        t23 = gsuper(2,3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine eom_grad_p_y

c     eom_grad_p_z
c     #############################################################
      subroutine eom_grad_p_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp
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
        ptot = a_p*0.5*jac*(sc1(i,j,kp)*sc2(i,j,k )
     .                     +sc1(i,j,k )*sc2(i,j,kp))

        t31 = gsuper(3,1)*ptot

        t32 = gsuper(3,2)*ptot

        t33 = gsuper(3,3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine eom_grad_p_z

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
#if !defined(vec_pot)
        real(8),pointer,dimension(:,:,:):: rho,px,py,pz,bx,by,bz,tmp
#else
        real(8),pointer,dimension(:,:,:):: rho,px,py,pz,ax,ay,az,tmp
#endif

        integer    :: advect

        real(8) :: max_dv_dt

        logical :: sym_st=.false.

      contains

ccc     eom_advc_x
ccc     #############################################################
cc      subroutine eom_advc_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
ccc     This routine assumes vec1->v (velocity), vec2->p (momentum).
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t11,t12,t13
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,ip,im,ipp,igpp,igp,igm
cc        real(8)    :: x,y,z,vxp,qip,ph,nabla_v_cnv(3,3)
cc        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,vis
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
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cc
cc        ijac0 = 1d0/jac0
cc        ijacp = 1d0/jacp
cc
ccc       Advective part
cc
cc        t11 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,1)
cc     .              +vec2(i ,j,k,1)*vec1(ip,j,k,1)
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,1)
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,1))
cc     .            *ijacp*ijac0*jac**2
cc
cc        t12 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,2)*ijacp
cc     .              +vec2(i ,j,k,1)*vec1(ip,j,k,2)*ijac0
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp)
cc     .            *jac
cc
cc        t13 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,3)
cc     .              +vec2(i ,j,k,1)*vec1(ip,j,k,3)
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,3)
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,3))
cc     .            *ijacp*ijac0*jac**2
cc
cccc        if (flag /= 0 .and. advect /= 2) then
cccc          ipp = min(ip+1,nx+1)
cccc          im  = max(i -1,0   )
cccc
cccc          igp  = ig + (ip  - i)
cccc          igm  = ig + (im  - i)
cccc          igpp = ig + (ipp - i)
cccc        endif
cccc
cccc        !Advection of px
cccc        t11 = 0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,1)
cccc     .             +vec1(i ,j,k,1)*vec2(ip,j,k,1))
cccc     .       *ijacp*ijac0*jac**2
cccc
cccc        ph = 0.5*(vec2(i,j,k,1)+vec2(ip,j,k,1))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t11    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vxp = t11/ph
cccc          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cccc
cccc          if (i == nx .or. i == 0) then
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,1),xx(igp),vec2(ip,j,k,1)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,1),xx(igm),vec2(im,j,k,1)
cccc     .               ,vxp,min(advect,2))
cccc          else
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,1),xx(igp),vec2(ip,j,k,1)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,1),xx(igm),vec2(im,j,k,1)
cccc     .               ,vxp,advect)
cccc          endif
cccc
cccc          t11 = vxp*qip
cccc        endif
cccc
cccc        !Advection of py
cccc        t12 = 0.5*( vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0
cccc     .            + vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp)*jac
cccc        ph = 0.5*(vec2(i,j,k,2)+vec2(ip,j,k,2))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t12    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vxp = t12/ph
cccc          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cccc
cccc          if (i == nx .or. i == 0) then
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,2),xx(igp),vec2(ip,j,k,2)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,2),xx(igm),vec2(im,j,k,2)
cccc     .               ,vxp,min(advect,2))
cccc          else
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,2),xx(igp),vec2(ip,j,k,2)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,2),xx(igm),vec2(im,j,k,2)
cccc     .               ,vxp,advect)
cccc          endif
cccc
cccc          t12 = vxp*qip
cccc        endif
cccc
cccc        !Advection of pz
cccc        t13 = 0.5*( vec1(i ,j,k,1)*vec2(ip,j,k,3)
cccc     .            + vec1(ip,j,k,1)*vec2(i ,j,k,3))
cccc     .           *ijacp*ijac0*jac**2
cccc        ph = 0.5*(vec2(i,j,k,3)+vec2(ip,j,k,3))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t13    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vxp = t13/ph
cccc          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cccc
cccc          if (i == nx .or. i == 0) then
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,3),xx(igp),vec2(ip,j,k,3)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,3),xx(igm),vec2(im,j,k,3)
cccc     .               ,vxp,min(advect,2))
cccc          else
cccc            qip=inter(xx(igpp),vec2(ipp,j,k,3),xx(igp),vec2(ip,j,k,3)
cccc     .               ,xx(ig  ),vec2(i  ,j,k,3),xx(igm),vec2(im,j,k,3)
cccc     .               ,vxp,advect)
cccc          endif
cccc
cccc          t13 = vxp*qip
cccc        endif
cc
ccc       Diffusive part
cc
cc        !Harmonic average for calculation of viscosity coeff. at faces
cc        if (nu /= 0d0) then
cc          vis = 2.*nuu(ip,j,k)*nuu(i,j,k)/(nuu(ip,j,k) + nuu(i,j,k))
cc
cccc          if (flag /= 0) then
cccc            nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cccc     .                        ,vec1(:,:,:,1)
cccc     .                        ,vec1(:,:,:,2)
cccc     .                        ,vec1(:,:,:,3),1)
cccc          else
cccc            nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cccc     .                        ,vec1(:,:,:,1)
cccc     .                        ,vec1(:,:,:,2)
cccc     .                        ,vec1(:,:,:,3),0)
cccc          endif
cccc
cccc          if (sym_st) then
cccc            nabla_v_cnv = matmul(gsuper,nabla_v)
cccc
cccc            nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
cccc
cccc            t11 = t11 - vis*nabla_v_cnv(1,1)
cccc            t12 = t12 - vis*nabla_v_cnv(1,2)
cccc            t13 = t13 - vis*nabla_v_cnv(1,3)
cccc          else
cccc            t11 = t11 - vis*( gsuper(1,1)*nabla_v(1,1)  
cccc     .                       +gsuper(1,2)*nabla_v(2,1)  
cccc     .                       +gsuper(1,3)*nabla_v(3,1) )
cccc
cccc            t12 = t12 - vis*( gsuper(1,1)*nabla_v(1,2)  
cccc     .                       +gsuper(1,2)*nabla_v(2,2)  
cccc     .                       +gsuper(1,3)*nabla_v(3,2) )
cccc
cccc            t13 = t13 - vis*( gsuper(1,1)*nabla_v(1,3)  
cccc     .                       +gsuper(1,2)*nabla_v(2,3)  
cccc     .                       +gsuper(1,3)*nabla_v(3,3) )
cccc          endif
cc
cc          if (flag /= 0) then
cc            nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                      ,vec1(:,:,:,1)
cc     .                      ,vec1(:,:,:,2)
cc     .                      ,vec1(:,:,:,3),1,cnv=.true.)
cc          else
cc            nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                      ,vec1(:,:,:,1)
cc     .                      ,vec1(:,:,:,2)
cc     .                      ,vec1(:,:,:,3),0,cnv=.true.)
cc          endif
cc
cc          if (sym_st) then
cc            nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
cc          endif
cc
cc          t11 = t11 - vis*nabla_v_cnv(1,1)
cc          t12 = t12 - vis*nabla_v_cnv(1,2)
cc          t13 = t13 - vis*nabla_v_cnv(1,3)
cc
cc        endif
cc
ccc       Postprocessing
cc
cc        if (flag /= 0) then
cc          ijac  = 1d0/jac
cc
cc          t11 = t11*ijac
cc          if (.not.alt_eom) t12 = t12*ijac
cc          t13 = t13*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine eom_advc_x
cc
ccc     eom_advc_y
ccc     #############################################################
cc      subroutine eom_advc_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t21,t22,t23
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,jp,jm,jpp,jgpp,jgp,jgm
cc        real(8)    :: x,y,z,vyp,qjp,ph,nabla_v_cnv(3,3)
cc        real(8)    :: jac,ijac,vis
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
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
ccc       Advective part
cc
cc        t21 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,1)
cc     .              +vec2(i,j ,k,2)*vec1(i,jp,k,1)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,1)
cc     .              +vec1(i,jp,k,2)*vec2(i,j ,k,1))
cc
cc        t22 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,2)
cc     .              +vec2(i,j ,k,2)*vec1(i,jp,k,2)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,2)
cc     .              +vec1(i,jp,k,2)*vec2(i,j ,k,2))
cc
cc        t23 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,3)
cc     .              +vec2(i,j ,k,2)*vec1(i,jp,k,3)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,3)
cc     .              +vec1(i,jp,k,2)*vec2(i,j ,k,3))
cc
cccc        if (flag /= 0 .and. advect /= 2) then
cccc          jpp = min(jp+1,ny+1)
cccc          jm  = max(j-1 ,0   )
cccc
cccc          jgp  = jg + (jp  - j)
cccc          jgm  = jg + (jm  - j)
cccc          jgpp = jg + (jpp - j)
cccc        endif
cccc
cccc        t21 = 0.5*(vec1(i,j ,k,2)*vec1(i,jp,k,1)
cccc     .            +vec1(i,jp,k,2)*vec1(i,j ,k,1))
cccc
cccc        ph = 0.5*(vec2(i,j,k,1)+vec2(i,jp,k,1))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t21    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vyp = t21/ph
cccc          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cccc
cccc          if (j == ny .or. j == 0) then
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,1),yy(jgp),vec2(i,jp,k,1)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,1),yy(jgm),vec2(i,jm,k,1)
cccc     .               ,vyp,min(advect,2))
cccc          else
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,1),yy(jgp),vec2(i,jp,k,1)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,1),yy(jgm),vec2(i,jm,k,1)
cccc     .               ,vyp,advect)
cccc          endif
cccc
cccc          t21 = vyp*qjp
cccc        endif
cccc
cccc        t22 = 0.5*( vec1(i,j ,k,2)*vec2(i,jp,k,2)
cccc     .            + vec1(i,jp,k,2)*vec2(i,j ,k,2))
cccc
cccc        ph = 0.5*(vec2(i,j,k,2)+vec2(i,jp,k,2))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t22    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vyp = t22/ph
cccc          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cccc
cccc          if (j == ny .or. j == 0) then
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,2),yy(jgp),vec2(i,jp,k,2)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,2),yy(jgm),vec2(i,jm,k,2)
cccc     .               ,vyp,min(advect,2))
cccc          else
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,2),yy(jgp),vec2(i,jp,k,2)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,2),yy(jgm),vec2(i,jm,k,2)
cccc     .               ,vyp,advect)
cccc          endif
cccc
cccc          t22 = vyp*qjp
cccc        endif
cccc
cccc        t23 = 0.5*( vec1(i,j ,k,2)*vec2(i,jp,k,3)
cccc     .            + vec1(i,jp,k,2)*vec2(i,j ,k,3))
cccc
cccc        ph = 0.5*(vec2(i,j,k,3)+vec2(i,jp,k,3))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t23    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vyp = t23/ph
cccc          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cccc
cccc          if (j == ny .or. j == 0) then
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,3),yy(jgp),vec2(i,jp,k,3)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,3),yy(jgm),vec2(i,jm,k,3)
cccc     .               ,vyp,min(advect,2))
cccc          else
cccc            qjp=inter(yy(jgpp),vec2(i,jpp,k,3),yy(jgp),vec2(i,jp,k,3)
cccc     .               ,yy(jg  ),vec2(i,j  ,k,3),yy(jgm),vec2(i,jm,k,3)
cccc     .               ,vyp,advect)
cccc          endif
cccc
cccc          t23 = vyp*qjp
cccc        endif
cc
ccc       Diffusive part
cc
cc        if (nu /= 0d0) then
cc          vis = 2.*nuu(i,jp,k)*nuu(i,j,k)/(nuu(i,jp,k) + nuu(i,j,k))
cc
cccc          if (flag /= 0) then
cccc            nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cccc     .                      ,vec1(:,:,:,1)
cccc     .                      ,vec1(:,:,:,2)
cccc     .                      ,vec1(:,:,:,3),2)
cccc          else
cccc            nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cccc     .                      ,vec1(:,:,:,1)
cccc     .                      ,vec1(:,:,:,2)
cccc     .                      ,vec1(:,:,:,3),0)
cccc          endif
cccc
cccc          if (sym_st) then
cccc            nabla_v_cnv = matmul(gsuper,nabla_v)
cccc
cccc            nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
cccc
cccc            t21 = t21 - vis*nabla_v_cnv(2,1)
cccc            t22 = t22 - vis*nabla_v_cnv(2,2)
cccc            t23 = t23 - vis*nabla_v_cnv(2,3)
cccc          else
cccc            t21 = t21 - vis*( gsuper(2,1)*nabla_v(1,1)  
cccc     .                       +gsuper(2,2)*nabla_v(2,1)  
cccc     .                       +gsuper(2,3)*nabla_v(3,1) )
cccc
cccc            t22 = t22 - vis*( gsuper(2,1)*nabla_v(1,2)  
cccc     .                       +gsuper(2,2)*nabla_v(2,2)  
cccc     .                       +gsuper(2,3)*nabla_v(3,2) )
cccc
cccc            t23 = t23 - vis*( gsuper(2,1)*nabla_v(1,3)  
cccc     .                       +gsuper(2,2)*nabla_v(2,3)  
cccc     .                       +gsuper(2,3)*nabla_v(3,3) )
cccc          endif
cc
cc         if (flag /= 0) then
cc            nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                      ,vec1(:,:,:,1)
cc     .                      ,vec1(:,:,:,2)
cc     .                      ,vec1(:,:,:,3),2,cnv=.true.)
cc          else
cc            nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                      ,vec1(:,:,:,1)
cc     .                      ,vec1(:,:,:,2)
cc     .                      ,vec1(:,:,:,3),0,cnv=.true.)
cc          endif
cc
cc          if (sym_st) then
cc            nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
cc          endif
cc
cc          t21 = t21 - vis*nabla_v_cnv(2,1)
cc          t22 = t22 - vis*nabla_v_cnv(2,2)
cc          t23 = t23 - vis*nabla_v_cnv(2,3)
cc
cc        endif
cc
ccc       Postprocessing
cc
cc        if (flag /= 0) then
cc          ijac = 1d0/jac
cc
cc          t21 = t21*ijac
cc          if (.not.alt_eom) t22 = t22*ijac
cc          t23 = t23*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine eom_advc_y
cc
ccc     eom_advc_z
ccc     #############################################################
cc      subroutine eom_advc_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t31,t32,t33
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,kp,km,kpp,kgpp,kgp,kgm
cc        real(8)    :: x,y,z,vzp,qkp,ph,nabla_v_cnv(3,3)
cc        real(8)    :: jac,ijac,vis
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
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
ccc       Advective part
cc
cc        t31 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,1)
cc     .              +vec2(i,j,k ,3)*vec1(i,j,kp,1)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,1)
cc     .              +vec1(i,j,kp,3)*vec2(i,j,k ,1))
cc
cc        t32 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,2)
cc     .              +vec2(i,j,k ,3)*vec1(i,j,kp,2)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,2)
cc     .              +vec1(i,j,kp,3)*vec2(i,j,k ,2))
cc
cc        t33 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,3)
cc     .              +vec2(i,j,k ,3)*vec1(i,j,kp,3)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,3)
cc     .              +vec1(i,j,kp,3)*vec2(i,j,k ,3))
cc
cccc        if (flag /= 0 .and. advect /= 2) then
cccc          kpp = min(kp+1,nx+1)
cccc          km  = max(k -1,0   )
cccc
cccc          kgp  = kg + (kp  - k)
cccc          kgm  = kg + (km  - k)
cccc          kgpp = kg + (kpp - k)
cccc        endif
cccc
cccc        t31 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,kp,1)
cccc     .            + vec1(i,j,kp,3)*vec2(i,j,k ,1) )
cccc
cccc        ph = 0.5*(vec2(i,j,k,1)+vec2(i,j,kp,1))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t31    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vzp = t31/ph
cccc          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cccc
cccc          if (k == nz .or. k == 0) then
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,1),zz(kgp),vec2(i,j,kp,1)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,1),zz(kgm),vec2(i,j,km,1)
cccc     .               ,vzp,min(advect,2))
cccc          else
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,1),zz(kgp),vec2(i,j,kp,1)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,1),zz(kgm),vec2(i,j,km,1)
cccc     .               ,vzp,advect)
cccc          endif
cccc
cccc          t31 = vzp*qkp
cccc        endif
cccc
cccc        t32 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,kp,2)
cccc     .            + vec1(i,j,kp,3)*vec2(i,j,k ,2) )
cccc
cccc        ph = 0.5*(vec2(i,j,k,2)+vec2(i,j,kp,2))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t32    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vzp = t32/ph
cccc          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cccc
cccc          if (k == nz .or. k == 0) then
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,2),zz(kgp),vec2(i,j,kp,2)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,2),zz(kgm),vec2(i,j,km,2)
cccc     .               ,vzp,min(advect,2))
cccc          else
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,2),zz(kgp),vec2(i,j,kp,2)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,2),zz(kgm),vec2(i,j,km,2)
cccc     .               ,vzp,advect)
cccc          endif
cccc
cccc          t32 = vzp*qkp
cccc        endif
cccc
cccc        t33 = 0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,3)
cccc     .            + vec1(i,j,k ,3)*vec2(i,j,kp,3) )
cccc
cccc        ph = 0.5*(vec2(i,j,k,3)+vec2(i,j,kp,3))
cccc
cccc        if (      flag   /= 0
cccc     .      .and. advect /= 2
cccc     .      .and. t33    /= 0d0
cccc     .      ) then
cccccc     .      .and. abs(ph) > 1d-5) then
cccc
cccccc          vzp = t33/ph
cccc          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cccc
cccc          if (k == nz .or. k == 0) then
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,3),zz(kgp),vec2(i,j,kp,3)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,3),zz(kgm),vec2(i,j,km,3)
cccc     .               ,vzp,min(advect,2))
cccc          else
cccc            qkp=inter(zz(kgpp),vec2(i,j,kpp,3),zz(kgp),vec2(i,j,kp,3)
cccc     .               ,zz(kg  ),vec2(i,j,k  ,3),zz(kgm),vec2(i,j,km,3)
cccc     .               ,vzp,advect)
cccc          endif
cccc
cccc          t33 = vzp*qkp
cccc        endif
cc
ccc       Diffusive part
cc
cc        if (nu /= 0d0) then
cc          vis = 2.*nuu(i,j,kp)*nuu(i,j,k)/(nuu(i,j,kp) + nuu(i,j,k))
cc
cccc          if (flag /= 0) then
cccc            nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cccc     .                      ,vec1(:,:,:,1)
cccc     .                      ,vec1(:,:,:,2)
cccc     .                      ,vec1(:,:,:,3),3)
cccc          else
cccc            nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cccc     .                      ,vec1(:,:,:,1)
cccc     .                      ,vec1(:,:,:,2)
cccc     .                      ,vec1(:,:,:,3),0)
cccc          endif
cccc
cccc          if (sym_st) then
cccc            nabla_v_cnv = matmul(gsuper,nabla_v)
cccc
cccc            nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
cccc
cccc            t31 = t31 - vis*nabla_v_cnv(3,1)
cccc            t32 = t32 - vis*nabla_v_cnv(3,2)
cccc            t33 = t33 - vis*nabla_v_cnv(3,3)
cccc          else
cccc            t31 = t31 - vis*( gsuper(3,1)*nabla_v(1,1)
cccc     .                       +gsuper(3,2)*nabla_v(2,1)
cccc     .                       +gsuper(3,3)*nabla_v(3,1) )
cccc
cccc            t32 = t32 - vis*( gsuper(3,1)*nabla_v(1,2)
cccc     .                       +gsuper(3,2)*nabla_v(2,2)
cccc     .                       +gsuper(3,3)*nabla_v(3,2) )
cccc
cccc            t33 = t33 - vis*( gsuper(3,1)*nabla_v(1,3)
cccc     .                       +gsuper(3,2)*nabla_v(2,3)
cccc     .                       +gsuper(3,3)*nabla_v(3,3) )
cccc          endif
cc
cc          if (flag /= 0) then
cc            nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                      ,vec1(:,:,:,1)
cc     .                      ,vec1(:,:,:,2)
cc     .                      ,vec1(:,:,:,3),3,cnv=.true.)
cc          else
cc            nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                      ,vec1(:,:,:,1)
cc     .                      ,vec1(:,:,:,2)
cc     .                      ,vec1(:,:,:,3),0,cnv=.true.)
cc          endif
cc
cc          if (sym_st) then
cc            nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
cc          endif
cc
cc          t31 = t31 - vis*nabla_v_cnv(3,1)
cc          t32 = t32 - vis*nabla_v_cnv(3,2)
cc          t33 = t33 - vis*nabla_v_cnv(3,3)
cc
cc        endif
cc
ccc       Postprocessing
cc
cc        if (flag /= 0) then
cc          ijac = 1d0/jac
cc
cc          t31 = t31*ijac
cc          if (.not.alt_eom) t32 = t32*ijac
cc          t33 = t33*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine eom_advc_z

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
c     This routine assumes vec1->v (velocity), vec2->p (momentum).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip,im,ipp,igpp,igp,igm
        real(8)    :: x,y,z,vxp,qip,ph,nabla_v_cnv(3,3)
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

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

c       Advective part

        t11 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,1)
     .              +vec2(i ,j,k,1)*vec1(ip,j,k,1)
     .              +vec1(i ,j,k,1)*vec2(ip,j,k,1)
     .              +vec1(ip,j,k,1)*vec2(i ,j,k,1))
     .            *ijacp*ijac0*jac**2

        t12 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,2)*ijacp
     .              +vec2(i ,j,k,1)*vec1(ip,j,k,2)*ijac0
     .              +vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0
     .              +vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp)
     .            *jac

        t13 = 0.25*( vec2(ip,j,k,1)*vec1(i ,j,k,3)
     .              +vec2(i ,j,k,1)*vec1(ip,j,k,3)
     .              +vec1(i ,j,k,1)*vec2(ip,j,k,3)
     .              +vec1(ip,j,k,1)*vec2(i ,j,k,3))
     .            *ijacp*ijac0*jac**2

cc        if (flag /= 0 .and. advect /= 2) then
cc          ipp = min(ip+1,nx+1)
cc          im  = max(i -1,0   )
cc
cc          igp  = ig + (ip  - i)
cc          igm  = ig + (im  - i)
cc          igpp = ig + (ipp - i)
cc        endif
cc
cc        !Advection of px
cc        t11 = 0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,1)
cc     .             +vec1(i ,j,k,1)*vec2(ip,j,k,1))
cc     .       *ijacp*ijac0*jac**2
cc
cc        ph = 0.5*(vec2(i,j,k,1)+vec2(ip,j,k,1))
cc
cc        if (      flag   /= 0
cc     .      .and. advect /= 2
cc     .      .and. t11    /= 0d0
cc     .      ) then
cccc     .      .and. abs(ph) > 1d-5) then
cc
cccc          vxp = t11/ph
cc          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cc
cc          if (i == nx .or. i == 0) then
cc            qip=inter(xx(igpp),vec2(ipp,j,k,1),xx(igp),vec2(ip,j,k,1)
cc     .               ,xx(ig  ),vec2(i  ,j,k,1),xx(igm),vec2(im,j,k,1)
cc     .               ,vxp,min(advect,2))
cc          else
cc            qip=inter(xx(igpp),vec2(ipp,j,k,1),xx(igp),vec2(ip,j,k,1)
cc     .               ,xx(ig  ),vec2(i  ,j,k,1),xx(igm),vec2(im,j,k,1)
cc     .               ,vxp,advect)
cc          endif
cc
cc          t11 = vxp*qip
cc        endif
cc
cc        !Advection of py
cc        t12 = 0.5*( vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0
cc     .            + vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp)*jac
cc        ph = 0.5*(vec2(i,j,k,2)+vec2(ip,j,k,2))
cc
cc        if (      flag   /= 0
cc     .      .and. advect /= 2
cc     .      .and. t12    /= 0d0
cc     .      ) then
cccc     .      .and. abs(ph) > 1d-5) then
cc
cccc          vxp = t12/ph
cc          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cc
cc          if (i == nx .or. i == 0) then
cc            qip=inter(xx(igpp),vec2(ipp,j,k,2),xx(igp),vec2(ip,j,k,2)
cc     .               ,xx(ig  ),vec2(i  ,j,k,2),xx(igm),vec2(im,j,k,2)
cc     .               ,vxp,min(advect,2))
cc          else
cc            qip=inter(xx(igpp),vec2(ipp,j,k,2),xx(igp),vec2(ip,j,k,2)
cc     .               ,xx(ig  ),vec2(i  ,j,k,2),xx(igm),vec2(im,j,k,2)
cc     .               ,vxp,advect)
cc          endif
cc
cc          t12 = vxp*qip
cc        endif
cc
cc        !Advection of pz
cc        t13 = 0.5*( vec1(i ,j,k,1)*vec2(ip,j,k,3)
cc     .            + vec1(ip,j,k,1)*vec2(i ,j,k,3))
cc     .           *ijacp*ijac0*jac**2
cc        ph = 0.5*(vec2(i,j,k,3)+vec2(ip,j,k,3))
cc
cc        if (      flag   /= 0
cc     .      .and. advect /= 2
cc     .      .and. t13    /= 0d0
cc     .      ) then
cccc     .      .and. abs(ph) > 1d-5) then
cc
cccc          vxp = t13/ph
cc          vxp = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cc
cc          if (i == nx .or. i == 0) then
cc            qip=inter(xx(igpp),vec2(ipp,j,k,3),xx(igp),vec2(ip,j,k,3)
cc     .               ,xx(ig  ),vec2(i  ,j,k,3),xx(igm),vec2(im,j,k,3)
cc     .               ,vxp,min(advect,2))
cc          else
cc            qip=inter(xx(igpp),vec2(ipp,j,k,3),xx(igp),vec2(ip,j,k,3)
cc     .               ,xx(ig  ),vec2(i  ,j,k,3),xx(igm),vec2(im,j,k,3)
cc     .               ,vxp,advect)
cc          endif
cc
cc          t13 = vxp*qip
cc        endif

c       Postprocessing

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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp,jm,jpp,jgpp,jgp,jgm
        real(8)    :: x,y,z,vyp,qjp,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,ijac,vis

c     Begin program

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c       Advective part

        t21 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,1)
     .              +vec2(i,j ,k,2)*vec1(i,jp,k,1)
     .              +vec1(i,j ,k,2)*vec2(i,jp,k,1)
     .              +vec1(i,jp,k,2)*vec2(i,j ,k,1))

        t22 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,2)
     .              +vec2(i,j ,k,2)*vec1(i,jp,k,2)
     .              +vec1(i,j ,k,2)*vec2(i,jp,k,2)
     .              +vec1(i,jp,k,2)*vec2(i,j ,k,2))

        t23 = 0.25*( vec2(i,jp,k,2)*vec1(i,j ,k,3)
     .              +vec2(i,j ,k,2)*vec1(i,jp,k,3)
     .              +vec1(i,j ,k,2)*vec2(i,jp,k,3)
     .              +vec1(i,jp,k,2)*vec2(i,j ,k,3))

cc        if (flag /= 0 .and. advect /= 2) then
cc          jpp = min(jp+1,ny+1)
cc          jm  = max(j-1 ,0   )
cc
cc          jgp  = jg + (jp  - j)
cc          jgm  = jg + (jm  - j)
cc          jgpp = jg + (jpp - j)
cc        endif
cc
cc        t21 = 0.5*(vec1(i,j ,k,2)*vec1(i,jp,k,1)
cc     .            +vec1(i,jp,k,2)*vec1(i,j ,k,1))
cc
cc        ph = 0.5*(vec2(i,j,k,1)+vec2(i,jp,k,1))
cc
cc        if (      flag   /= 0
cc     .      .and. advect /= 2
cc     .      .and. t21    /= 0d0
cc     .      ) then
cccc     .      .and. abs(ph) > 1d-5) then
cc
cccc          vyp = t21/ph
cc          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cc
cc          if (j == ny .or. j == 0) then
cc            qjp=inter(yy(jgpp),vec2(i,jpp,k,1),yy(jgp),vec2(i,jp,k,1)
cc     .               ,yy(jg  ),vec2(i,j  ,k,1),yy(jgm),vec2(i,jm,k,1)
cc     .               ,vyp,min(advect,2))
cc          else
cc            qjp=inter(yy(jgpp),vec2(i,jpp,k,1),yy(jgp),vec2(i,jp,k,1)
cc     .               ,yy(jg  ),vec2(i,j  ,k,1),yy(jgm),vec2(i,jm,k,1)
cc     .               ,vyp,advect)
cc          endif
cc
cc          t21 = vyp*qjp
cc        endif
cc
cc        t22 = 0.5*( vec1(i,j ,k,2)*vec2(i,jp,k,2)
cc     .            + vec1(i,jp,k,2)*vec2(i,j ,k,2))
cc
cc        ph = 0.5*(vec2(i,j,k,2)+vec2(i,jp,k,2))
cc
cc        if (      flag   /= 0
cc     .      .and. advect /= 2
cc     .      .and. t22    /= 0d0
cc     .      ) then
cccc     .      .and. abs(ph) > 1d-5) then
cc
cccc          vyp = t22/ph
cc          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cc
cc          if (j == ny .or. j == 0) then
cc            qjp=inter(yy(jgpp),vec2(i,jpp,k,2),yy(jgp),vec2(i,jp,k,2)
cc     .               ,yy(jg  ),vec2(i,j  ,k,2),yy(jgm),vec2(i,jm,k,2)
cc     .               ,vyp,min(advect,2))
cc          else
cc            qjp=inter(yy(jgpp),vec2(i,jpp,k,2),yy(jgp),vec2(i,jp,k,2)
cc     .               ,yy(jg  ),vec2(i,j  ,k,2),yy(jgm),vec2(i,jm,k,2)
cc     .               ,vyp,advect)
cc          endif
cc
cc          t22 = vyp*qjp
cc        endif
cc
cc        t23 = 0.5*( vec1(i,j ,k,2)*vec2(i,jp,k,3)
cc     .            + vec1(i,jp,k,2)*vec2(i,j ,k,3))
cc
cc        ph = 0.5*(vec2(i,j,k,3)+vec2(i,jp,k,3))
cc
cc        if (      flag   /= 0
cc     .      .and. advect /= 2
cc     .      .and. t23    /= 0d0
cc     .      ) then
cccc     .      .and. abs(ph) > 1d-5) then
cc
cccc          vyp = t23/ph
cc          vyp = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cc
cc          if (j == ny .or. j == 0) then
cc            qjp=inter(yy(jgpp),vec2(i,jpp,k,3),yy(jgp),vec2(i,jp,k,3)
cc     .               ,yy(jg  ),vec2(i,j  ,k,3),yy(jgm),vec2(i,jm,k,3)
cc     .               ,vyp,min(advect,2))
cc          else
cc            qjp=inter(yy(jgpp),vec2(i,jpp,k,3),yy(jgp),vec2(i,jp,k,3)
cc     .               ,yy(jg  ),vec2(i,j  ,k,3),yy(jgm),vec2(i,jm,k,3)
cc     .               ,vyp,advect)
cc          endif
cc
cc          t23 = vyp*qjp
cc        endif

c       Postprocessing

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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp,km,kpp,kgpp,kgp,kgm
        real(8)    :: x,y,z,vzp,qkp,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,ijac,vis

c     Begin program

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c       Advective part

        t31 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,1)
     .              +vec2(i,j,k ,3)*vec1(i,j,kp,1)
     .              +vec1(i,j,k ,3)*vec2(i,j,kp,1)
     .              +vec1(i,j,kp,3)*vec2(i,j,k ,1))

        t32 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,2)
     .              +vec2(i,j,k ,3)*vec1(i,j,kp,2)
     .              +vec1(i,j,k ,3)*vec2(i,j,kp,2)
     .              +vec1(i,j,kp,3)*vec2(i,j,k ,2))

        t33 = 0.25*( vec2(i,j,kp,3)*vec1(i,j,k ,3)
     .              +vec2(i,j,k ,3)*vec1(i,j,kp,3)
     .              +vec1(i,j,k ,3)*vec2(i,j,kp,3)
     .              +vec1(i,j,kp,3)*vec2(i,j,k ,3))

cc        if (flag /= 0 .and. advect /= 2) then
cc          kpp = min(kp+1,nx+1)
cc          km  = max(k -1,0   )
cc
cc          kgp  = kg + (kp  - k)
cc          kgm  = kg + (km  - k)
cc          kgpp = kg + (kpp - k)
cc        endif
cc
cc        t31 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,kp,1)
cc     .            + vec1(i,j,kp,3)*vec2(i,j,k ,1) )
cc
cc        ph = 0.5*(vec2(i,j,k,1)+vec2(i,j,kp,1))
cc
cc        if (      flag   /= 0
cc     .      .and. advect /= 2
cc     .      .and. t31    /= 0d0
cc     .      ) then
cccc     .      .and. abs(ph) > 1d-5) then
cc
cccc          vzp = t31/ph
cc          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cc
cc          if (k == nz .or. k == 0) then
cc            qkp=inter(zz(kgpp),vec2(i,j,kpp,1),zz(kgp),vec2(i,j,kp,1)
cc     .               ,zz(kg  ),vec2(i,j,k  ,1),zz(kgm),vec2(i,j,km,1)
cc     .               ,vzp,min(advect,2))
cc          else
cc            qkp=inter(zz(kgpp),vec2(i,j,kpp,1),zz(kgp),vec2(i,j,kp,1)
cc     .               ,zz(kg  ),vec2(i,j,k  ,1),zz(kgm),vec2(i,j,km,1)
cc     .               ,vzp,advect)
cc          endif
cc
cc          t31 = vzp*qkp
cc        endif
cc
cc        t32 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,kp,2)
cc     .            + vec1(i,j,kp,3)*vec2(i,j,k ,2) )
cc
cc        ph = 0.5*(vec2(i,j,k,2)+vec2(i,j,kp,2))
cc
cc        if (      flag   /= 0
cc     .      .and. advect /= 2
cc     .      .and. t32    /= 0d0
cc     .      ) then
cccc     .      .and. abs(ph) > 1d-5) then
cc
cccc          vzp = t32/ph
cc          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cc
cc          if (k == nz .or. k == 0) then
cc            qkp=inter(zz(kgpp),vec2(i,j,kpp,2),zz(kgp),vec2(i,j,kp,2)
cc     .               ,zz(kg  ),vec2(i,j,k  ,2),zz(kgm),vec2(i,j,km,2)
cc     .               ,vzp,min(advect,2))
cc          else
cc            qkp=inter(zz(kgpp),vec2(i,j,kpp,2),zz(kgp),vec2(i,j,kp,2)
cc     .               ,zz(kg  ),vec2(i,j,k  ,2),zz(kgm),vec2(i,j,km,2)
cc     .               ,vzp,advect)
cc          endif
cc
cc          t32 = vzp*qkp
cc        endif
cc
cc        t33 = 0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,3)
cc     .            + vec1(i,j,k ,3)*vec2(i,j,kp,3) )
cc
cc        ph = 0.5*(vec2(i,j,k,3)+vec2(i,j,kp,3))
cc
cc        if (      flag   /= 0
cc     .      .and. advect /= 2
cc     .      .and. t33    /= 0d0
cc     .      ) then
cccc     .      .and. abs(ph) > 1d-5) then
cc
cccc          vzp = t33/ph
cc          vzp = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cc
cc          if (k == nz .or. k == 0) then
cc            qkp=inter(zz(kgpp),vec2(i,j,kpp,3),zz(kgp),vec2(i,j,kp,3)
cc     .               ,zz(kg  ),vec2(i,j,k  ,3),zz(kgm),vec2(i,j,km,3)
cc     .               ,vzp,min(advect,2))
cc          else
cc            qkp=inter(zz(kgpp),vec2(i,j,kpp,3),zz(kgp),vec2(i,j,kp,3)
cc     .               ,zz(kg  ),vec2(i,j,k  ,3),zz(kgm),vec2(i,j,km,3)
cc     .               ,vzp,advect)
cc          endif
cc
cc          t33 = vzp*qkp
cc        endif

c       Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine eom_advc_z

c     eom_diff_x
c     #############################################################
      subroutine eom_diff_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
c     This routine assumes vec1->v (velocity), vec2->p (momentum).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip,im,ipp,igpp,igp,igm
        real(8)    :: x,y,z,vxp,qip,ph,nabla_v_cnv(3,3)
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

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

c     Diffusive part

        !Harmonic average for calculation of viscosity coeff. at faces
        vis = 2.*nuu(ip,j,k)*nuu(i,j,k)/(nuu(ip,j,k) + nuu(i,j,k))

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),1,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),0,cnv=.true.)
        endif

        if (sym_st) then
          nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
        endif

        t11 = -vis*nabla_v_cnv(1,1)
        t12 = -vis*nabla_v_cnv(1,2)
        t13 = -vis*nabla_v_cnv(1,3)

c       Postprocessing

        if (flag /= 0) then
          ijac  = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine eom_diff_x

c     eom_diff_y
c     #############################################################
      subroutine eom_diff_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp,jm,jpp,jgpp,jgp,jgm
        real(8)    :: x,y,z,vyp,qjp,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,ijac,vis

c     Begin program

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c       Diffusive part

        vis = 2.*nuu(i,jp,k)*nuu(i,j,k)/(nuu(i,jp,k) + nuu(i,j,k))

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),2,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),0,cnv=.true.)
        endif

        if (sym_st) then
          nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
        endif

        t21 = -vis*nabla_v_cnv(2,1)
        t22 = -vis*nabla_v_cnv(2,2)
        t23 = -vis*nabla_v_cnv(2,3)

c       Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine eom_diff_y

c     eom_diff_z
c     #############################################################
      subroutine eom_diff_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp,km,kpp,kgpp,kgp,kgm
        real(8)    :: x,y,z,vzp,qkp,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,ijac,vis

c     Begin program

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c       Diffusive part

        vis = 2.*nuu(i,j,kp)*nuu(i,j,k)/(nuu(i,j,kp) + nuu(i,j,k))

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),3,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),0,cnv=.true.)
        endif

        if (sym_st) then
          nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
        endif

        t31 = -vis*nabla_v_cnv(3,1)
        t32 = -vis*nabla_v_cnv(3,2)
        t33 = -vis*nabla_v_cnv(3,3)

c       Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine eom_diff_z

c     si_op
c     #####################################################################
      function si_op(i,j,k,nx,ny,nz,igx,igy,igz) result(cnv)

      implicit none

      integer    :: i,j,k,igx,igy,igz,nx,ny,nz
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

        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ip,igrid
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

        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: jp,igrid
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

        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: kp,igrid
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

      integer    :: i,j,k,igx,igy,igz
      real(8)    :: psit(3)

      integer    :: ip,im,jp,jm,kp,km,ig,jg,kg

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
     .                          ,arr1,arr2,arr3,in_repr,out_repr)

c     ----------------------------------------------------------------------
c     Transforms vectors components in arrays arr1,arr2,arr3 from
c     Cartesian to curvilinear and viceversa, or from curvilinear to
c     curvilinear, depending on the value of in_repr (which indicates the
c     representation of the input vector) and out_repr (which indicates
c     representation of the input vector). Possible representations are:
c     'car','cnv',and 'cov'.
c     ----------------------------------------------------------------------

      use grid

      implicit none

c     Input variables

        integer    :: imin,imax,jmin,jmax,kmin,kmax
        integer    :: igx,igy,igz
        real(8)    :: arr1(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr2(imin:imax,jmin:jmax,kmin:kmax)
     .               ,arr3(imin:imax,jmin:jmax,kmin:kmax)

        character(3) :: in_repr,out_repr

c     Local variables

        integer    :: i,j,k
        real(8)    :: vec(3)

        logical    :: covariant,cnvariant,cartesian

c     Begin program

        cartesian = (in_repr == 'car')
        covariant = (in_repr == 'cov')
        cnvariant = (.not. cartesian .and. .not. covariant)

        select case(out_repr)
        case('car')

          if (cartesian) return

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

        case('cnv')

          if (cnvariant) return

          if (cartesian) then
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformVectorToCurvilinear
     .                 (i,j,k,igx,igy,igz
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,covariant
     .                 ,vec(1),vec(2),vec(3))

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

                  call transformFromCurvToCurv
     .                 (i,j,k,igx,igy,igz
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,vec(1),vec(2),vec(3)
     .                 ,covariant,0)

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          endif

        case('cov')

          if (covariant) return

          if (cartesian) then
            do k=kmin,kmax
              do j=jmin,jmax
                do i=imin,imax

                  call transformVectorToCurvilinear
     .                 (i,j,k,igx,igy,igz
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,covariant
     .                 ,vec(1),vec(2),vec(3))

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

                  call transformFromCurvToCurv
     .                 (i,j,k,igx,igy,igz
     .                 ,vec(1),vec(2),vec(3)
     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
     .                 ,covariant,0)

                  arr1(i,j,k) = vec(1)
                  arr2(i,j,k) = vec(2)
                  arr3(i,j,k) = vec(3)
                  
                enddo
              enddo
            enddo
          endif

        end select

cc        if (to_cartsn) then
cc
cc          if (cartesian) return
cc
cc          do k=kmin,kmax
cc            do j=jmin,jmax
cc              do i=imin,imax
cc
cc                call transformVectorToCartesian
cc     .               (i,j,k,igx,igy,igz
cc     .               ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
cc     .               ,covariant
cc     .               ,vec(1),vec(2),vec(3))
cc
cc                arr1(i,j,k) = vec(1)
cc                arr2(i,j,k) = vec(2)
cc                arr3(i,j,k) = vec(3)
cc                
cc              enddo
cc            enddo
cc          enddo
cc
cc        else
cc
cc          if (cartesian) then
cc            do k=kmin,kmax
cc              do j=jmin,jmax
cc                do i=imin,imax
cc
cc                  call transformVectorToCurvilinear
cc     .                 (i,j,k,igx,igy,igz
cc     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
cc     .                 ,covariant
cc     .                 ,vec(1),vec(2),vec(3))
cc
cc                  arr1(i,j,k) = vec(1)
cc                  arr2(i,j,k) = vec(2)
cc                  arr3(i,j,k) = vec(3)
cc                  
cc                enddo
cc              enddo
cc            enddo
cc          else
cc            do k=kmin,kmax
cc              do j=jmin,jmax
cc                do i=imin,imax
cc
cc                  call transformFromCurvToCurv
cc     .                 (i,j,k,igx,igy,igz
cc     .                 ,arr1(i,j,k),arr2(i,j,k),arr3(i,j,k)
cc     .                 ,covariant
cc     .                 ,vec(1),vec(2),vec(3))
cc
cc                  arr1(i,j,k) = vec(1)
cc                  arr2(i,j,k) = vec(2)
cc                  arr3(i,j,k) = vec(3)
cc                  
cc                enddo
cc              enddo
cc            enddo
cc          endif
cc
cc        endif

      end subroutine transformVector

      end module vectorOps
