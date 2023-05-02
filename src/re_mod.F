      module re

      use nlfunction_setup

      use par_int

      use local_BCS_variables, ONLY: order_bc

      real(8),pointer,dimension(:,:,:),private :: nre_nm,nre_n,nre_np

      logical,private :: B_weight,do_3d
      real(8),private :: G_tol,Dstar,c_va
      integer,private :: sp_order

      logical :: re_pkg=.false.

      type(mg_array),target,private :: gv_exb,bbmag
        
      contains

c     re_allocVariables
c     ###############################################################
      subroutine re_allocVariables

      implicit none

      integer :: igrid,nx,ny,nz
      
      igrid = 1

      nx = gv%gparams%nxv(igrid)
      ny = gv%gparams%nyv(igrid)
      nz = gv%gparams%nzv(igrid)

      allocate(nre_nm (0:nx+1,0:ny+1,0:nz+1)
     $        ,nre_n  (0:nx+1,0:ny+1,0:nz+1)
     $        ,nre_np (0:nx+1,0:ny+1,0:nz+1))

      nre_nm = 0d0 ; nre_n = 0d0 ; nre_np = 0d0
      
      end subroutine re_allocVariables

c     re_deallocVariables
c     ###############################################################
      subroutine re_deallocVariables

      implicit none

      deallocate(nre_nm
     $          ,nre_n 
     $          ,nre_np)
      
      end subroutine re_deallocVariables

c     re_init
c     ###############################################################
      subroutine re_init(nre)

        implicit none

        real(8) :: nre(0:,0:,0:)

        nre_nm = nre
        nre_n  = nre
        nre_np = nre
        
      end subroutine re_init

c     re_updt
c     ###############################################################
      subroutine re_updt(vnm,vn,vnp,iout)

        implicit none

        !Call variables
        type(var_array),pointer :: vnm,vn,vnp
        integer :: iout

        !Save old timestep solutions
        nre_nm = nre_n
        nre_n  = nre_np
        
        !Eulerian update for nre_np
        call re_Eul(vnm,vn,vnp,iout)

        !Lagrangian update for nre_np
        call re_lag(vnm,vn,vnp,iout)
        
      end subroutine re_updt

c     re_Eul
c     ###############################################################
      subroutine re_Eul(v_nm,v_n,v_np,iout)

        implicit none

c       Call variables

        type(var_array),pointer :: v_nm,v_n,v_np
        integer :: iout

c       Local variables

        integer :: nn,nx,ny,nz,igrid,bcs(6,3),bc_nre(6),guess,order

        real(8),allocatable,dimension(:,:) :: x,y
        real(8),allocatable,dimension(:,:,:) :: s0
        real(8),allocatable,dimension(:,:,:,:) :: E

c       Begin program
        
        igrid = 1

        nx = gv%gparams%nxv(igrid)
        ny = gv%gparams%nyv(igrid)
        nz = gv%gparams%nzv(igrid)

        nn = nx*ny*nz

        call allocateMGArray(gv%gparams,3,gv_exb)

        allocate(E(0:nx+1,0:ny+1,0:nz+1,3))
        allocate(s0(0:nx+1,0:ny+1,0:nz+1))

c       Prepare solver coefficients

        E = XformToCnv(gv%gparams,igrid,E_ni)
        gv_exb%grid(igrid)%array =
     $       crossproduct(gv%gparams,igrid,E,bcnv,.true.)
        gv_exb%grid(igrid)%array(:,:,:,1) =
     $       gv_exb%grid(igrid)%array(:,:,:,1)
     $      /vectornorm  (gv%gparams,igrid,bcnv,.false.)
        gv_exb%grid(igrid)%array(:,:,:,2) =
     $       gv_exb%grid(igrid)%array(:,:,:,2)
     $      /vectornorm  (gv%gparams,igrid,bcnv,.false.)
        gv_exb%grid(igrid)%array(:,:,:,3) =
     $       gv_exb%grid(igrid)%array(:,:,:,3)
     $      /vectornorm  (gv%gparams,igrid,bcnv,.false.)

        bcs(:,1) = v_np%array_var(IVX)%bconds
        bcs(:,2) = v_np%array_var(IVY)%bconds
        bcs(:,3) = v_np%array_var(IVZ)%bconds

        order = 0
        call restrictMGArray(IVX,3,gv_exb,bcs,igrid,order
     .                      ,order_bc=order_bc)

c       Store time-n status and solve for time-(n+1)

        nre_np = -(cn*nre_n + cnm*nre_nm)/cnp
        
c       Form RE rhs
c$$$        call formRErhs(nn,x,y,rhssi,igrid)

c       Solve n_re system

        allocate(x(nn,1),y(nn,1))
        
        x(:,1) = reshape(nre_np(1:nx,1:ny,1:nz),shape(x(:,1)))

        guess = 1
        bc_nre = v_np%array_var(IRHO)%bconds
        call cSolver(gv%gparams,1,nn,y,x,bc_nre
     .              ,igrid,iout,guess,nre_mtvc,vol_wgt
     .              ,tol           = rtol !pc_tol
     .              ,gm_driver     = .true. !pc_gm_driver
     .              ,mg_debug      = mg_debug
     .              ,mg_galerkin   = mg_glrkin
     .              ,mg_order_res  = mg_ores
     .              ,mg_order_prol = mg_oprol
     .              ,mg_gm_coarse_solve  = mg_gm_coarse
     .              ,mg_vcyc       = mg_vcyc
     .              ,mg_smooth     = sm_type
     .              ,sm_it         = sm_iter
     .              ,sm_omega      = sm_omega
     .              ,sm_ncolors    = sm_ncolors
     .              ,sm_line_relax = sm_line_relax
     .              ,sm_zebra_relax= sm_zebra_relax
     .              ,sm_zebra_it   = sm_zebra_it
     .              )

c       Return new-time RE density (with BCs)

        nre_np = reshape(x(:,1),shape(nre_n(1:nx,1:ny,1:nz)))

        s0 = 0d0
        call setBC(gv%gparams,IRHO,nx,ny,nz,nre_np,s0
     .            ,v_np%array_var(IRHO)%bconds
     .            ,igrid,igrid,igrid,iorder=order_bc)

c       Garbage collection
        
        deallocate(x,y,E,s0)
        call deallocateMGArray(gv_exb)

      end subroutine re_Eul

c     nre_mtvc
c     ####################################################################
      subroutine nre_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
c     --------------------------------------------------------------------
c     This subroutine calculates, for given x, y = A.x  matrix-free
c     for the continuity equation.
c     In call:
c      * gpos: vector index of position on the numerical grid
c            + If gpos = i + nx*(j-1) + ny*nx*(k-1), then only 
c              surrounding stencil is filled (9-pt stencil in 2D
c              , 27-pt stencil in 3D).
c            + If gpos = 0, all the grid is considered.
c            + If gpos < 0, all grid is mapped, but operations are 
c              restricted to stencil of abs(gpos) (useful for
c              matrix-light GS)
c      * neq: number of coupled equations
c      * ntot: total number of unknowns: neq*nx*ny*nz
c      * x(ntot): input vector
c      * y(ntot): output vector
c      * igrid: grid level
c      * bcnf: boundary conditions on x vector.
c     --------------------------------------------------------------------

      use matvec

      implicit none

c     Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)
      type(grid_mg_def),pointer :: g_def

c     Local variables

      integer    :: isig,ijk,ijkg,nnx,nny,nnz,ip,im,jp,jm,kp,km
      integer    :: imin,imax,jmin,jmax,kmin,kmax,i,j,k,ig,jg,kg
 
      real(8),allocatable,dimension(:,:,:,:) :: drho
      real(8),pointer    ,dimension(:,:,:,:) :: v0_cnv

      real(8),pointer    ,dimension(:,:,:) :: jac

      real(8)    :: upwind,slap,dvol

      logical    :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c     Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nnx = g_def%nxv(igrid)
      nny = g_def%nyv(igrid)
      nnz = g_def%nzv(igrid)

      jac => g_def%gmetric%grid(igrid)%jac

c     Find limits for loops

      call limits(mg_ctx_l,abs(gpos),igrid
     .           ,imin,imax,jmin,jmax,kmin,kmax)

c     Map vector x to array for processing

      allocate(drho(0:nnx+1,0:nny+1,0:nnz+1,neq))

      drho = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped
      !For finding the diagonal, gpos > 0
      call fillArray(g_def,igrid,neq,x,drho,bcnd,order_bc,gpos=gpos
     .              ,icomp=(/IRHO/))

c     Map velocity components

      call point_to_MG_array(nnx,nny,nnz,igrid,v0_cnv,gv_exb)
 
c     Calculate matrix-vector product

      do k = kmin,kmax
        do j = jmin,jmax
          do i = imin,imax

            call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            ijk = i + nnx*(j-1) + nnx*nny*(k-1)

            if (MG_get_res_vol_wgt()) then
              dvol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)
            else
              dvol = 1d0
            endif

            !Conservative advection
            upwind = flx_advec(g_def,i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                       ,v0_cnv,drho(:,:,:,1),UPWD)

c$$$            ! Diffusion
c$$$            if (dd > 0d0) then
c$$$              slap = dd*lap(g_def,i,j,k,nnx,nny,nnz,igrid,igrid,igrid
c$$$     .                     ,drho(:,:,:,1))
c$$$            else
              slap = 0d0
c$$$            endif

            y(ijk)=dvol*(cnp/dt*drho(i,j,k,1) + alpha*(upwind - slap))
     .                 *jac(i,j,k)
          enddo
        enddo
      enddo

c     End program

      deallocate(drho)
      nullify(v0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine nre_mtvc

c     re_lag
c     ###############################################################
      subroutine re_lag(vnm,vn,vnp,iout)

        implicit none

        !Call variables
        type(var_array),pointer :: vnm,vn,vnp
        integer :: iout

        !Local variables
        integer :: igrid,nx,ny,nz

        real(8) :: rkernel
        
        real(8),allocatable,dimension(:,:,:) :: tmp_nm,tmp_n,tmp_np,arr
     $                                         ,bmag_nm,bmag_n,bmag_np
        real(8),allocatable,dimension(:,:,:,:) :: bb,bcnv_nm,bcnv_n
     $                                           ,bcnv_np

        !Begin program
        
        if (iout > 0 .and. my_rank == 0) then
          write (*,*)
          write (*,*) 'Performing RE integration in postproc...'
        endif

        igrid = 1

        nx = gv%gparams%nxv(igrid)
        ny = gv%gparams%nyv(igrid)
        nz = gv%gparams%nzv(igrid)

        allocate(arr    (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_nm (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_n  (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_np (0:nx+1,0:ny+1,0:nz+1)
     $          ,bmag_n (0:nx+1,0:ny+1,0:nz+1)
     $          ,bmag_nm(0:nx+1,0:ny+1,0:nz+1)
     $          ,bmag_np(0:nx+1,0:ny+1,0:nz+1)
     $          ,bb     (0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bcnv_n (0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bcnv_nm(0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bcnv_np(0:nx+1,0:ny+1,0:nz+1,3))

        arr = 0d0
        
        !Setup lagrangian integration
        G_tol = 1d-8       !Green's function tol
        sp_order = 2       !Spline order of interpolation
        B_weight = .false.
        do_3d    = .true.
        Dstar = 0.15d0     !For stability
        c_va  = 1d2        !c/v_A parameter
        
        !Set orbit tols to Green's function tol
        call set_orbit_tols(G_tol,G_tol)

        !Find orbit flow
        bcnv_nm(:,:,:,1) = vnm%array_var(IBX)%array
        bcnv_nm(:,:,:,2) = vnm%array_var(IBY)%array
        bcnv_nm(:,:,:,3) = vnm%array_var(IBZ)%array
        bcnv_n (:,:,:,1) = vn %array_var(IBX)%array
        bcnv_n (:,:,:,2) = vn %array_var(IBY)%array
        bcnv_n (:,:,:,3) = vn %array_var(IBZ)%array
        bcnv_np(:,:,:,1) = vnp%array_var(IBX)%array
        bcnv_np(:,:,:,2) = vnp%array_var(IBY)%array
        bcnv_np(:,:,:,3) = vnp%array_var(IBZ)%array
        
        bmag_nm = vectorNorm(gv%gparams,igrid,bcnv_nm,.false.)
        bmag_n  = vectorNorm(gv%gparams,igrid,bcnv_n ,.false.)
        bmag_np = vectorNorm(gv%gparams,igrid,bcnv_np,.false.)

        bmag_nm = sqrt(bmag_nm)
        bmag_n  = sqrt(bmag_n )
        bmag_np = sqrt(bmag_np)
        
        bb(:,:,:,1) = c_va*bb(:,:,:,1)/bmag_np
        bb(:,:,:,2) = c_va*bb(:,:,:,2)/bmag_np
        bb(:,:,:,3) = c_va*bb(:,:,:,3)/bmag_np

        !Initialize Lagrangian step (!advection-diff mode with diff coeff Dstar)
        rkernel = -1d0/sqrt(Dstar) 
        call lag_init(gv%gparams,igrid,bb,rkernel,iout
     .               ,G_tol,sp_order,B_weight,do_3d)

        !Find temp from pressure (automatically takes care of jacobian factor)
        tmp_nm = nre_nm/bmag_nm
        tmp_n  = nre_n /bmag_n
        tmp_np = nre_np/bmag_np

        call lag_update_bdf2(gv%gparams,igrid,IRHO
     .                      ,vnp%array_var(IRHO)%bconds
     $                      ,tmp_nm,tmp_n,tmp_np,arr,iout
c$$$     .                      ,lsode_ext=lsode_3d_ext
     .                      ,no_global_splines=.true.)

        !Destroy Lagrangian step
        call lag_destroy

        !Return RE density
        nre_np = nre_np*bmag_np

        
        !Check and correct for particle conservation

        !Impose BCs
        call setBC(gv%gparams,IRHO,nx,ny,nz
     .            ,nre_np,arr
     .            ,vnp%array_var(IRHO)%bconds
     .            ,igrid,igrid,igrid,iorder=order_bc)

        !Garbage collection
        deallocate(tmp_nm,tmp_n,tmp_np,bb,bmag_nm,bmag_n,bmag_np
     $            ,bcnv_n,bcnv_nm,bcnv_np,arr)
        
      end subroutine re_lag

c$$$c lsode_3d_ext
c$$$c #########################################################################
c$$$      subroutine lsode_3d_ext(neq,t,y,ydot)
c$$$
c$$$c -------------------------------------------------------------------------
c$$$c     Driver routine for LSODE integration: ydot = ff(y,t)
c$$$c -------------------------------------------------------------------------
c$$$
c$$$      use par_int
c$$$      use kernels       
c$$$      use grid_splines
c$$$     
c$$$      implicit none
c$$$
c$$$c Call variables
c$$$
c$$$      integer :: neq
c$$$      real(8) :: t,y(neq),ydot(neq)
c$$$
c$$$c Local variables
c$$$
c$$$      integer :: ierror
c$$$      real(8) :: x1,y1,z1,flds(size(sp_coef,4)),B2,iB2
c$$$      logical :: exl1,exl2,exl3,exl4
c$$$      logical :: exl1a,exl2a
c$$$      
c$$$      real(8) :: gammaRi,muR
c$$$      real(8) :: Bstar, Brem, rratio, tmp, epsinvbyq,dt_e,dt2_e
c$$$      real(8) :: x2,y2,z2
c$$$      real(8) :: Ersource,ALDrsink,ALDpsink,smlr_ratio
c$$$
c$$$c Begin program
c$$$
c$$$      if (neq /= 9) then
c$$$        call pstop('lsode_3d_bdf2','Wrong number of equations')
c$$$      endif
c$$$
c$$$c Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)
c$$$c      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1
c$$$c     .                  ,B2=B2)
c$$$c
c$$$c       ydot(1:5)= 0d0       
c$$$      !Prepare variables for Lagrangian integral
c$$$      dt_e = dt*iepsilon_r
c$$$      dt2_e = (dt+dtold)*iepsilon_r
c$$$
c$$$      ! ratio of minor and major radii 
c$$$      rratio = flux_r*imajor_r
c$$$      smlr_ratio = flux_r*iminor_r
c$$$
c$$$      ! qsafe 
c$$$      qsafe = 2.1d0+2d0*(smlr_ratio)*(smlr_ratio)
c$$$
c$$$      !  1/epsilon, epsilon is small << 1
c$$$      epsinvbyq =1d0/qsafe      !v=1/epsilon used in kernels  
c$$$
c$$$      ! magnetic compnent remaining in dtheta/dt
c$$$      Brem  = 1d0/sqrt(1d0+(rratio*epsinvbyq)*(rratio*epsinvbyq))
c$$$      if( toroidB) Brem  = 1d0  ! Bmag = Btoroidal 
c$$$
c$$$      tmp   = 1d0/(1d0+rratio*cos(y(1)))
c$$$      Bstar = tmp/Brem
c$$$      
c$$$      ! Lorentz factor inverse
c$$$      gammaRi = 1d0/sqrt(1d0+y(2)*y(2)+y(3)*y(3))
c$$$
c$$$      ! Relativisitic magnetic moment
c$$$      muR = 0.5*y(3)*y(3)/Bstar             ! pp^2/2B* 
c$$$
c$$$      if (iEorbit) then 
c$$$      ! contains electric field and radiation damping  
c$$$      ! compute electric field
c$$$        Ersource=  Evecr*tmp*epsilon_r
c$$$      ! compute radiation sink
c$$$        ALDrsink=-sigma_r*gammaRi*abs(y(2))*y(3)*y(3)*epsilon_r 
c$$$        ALDpsink=-sigma_r*gammaRi*y(3)*(1+y(3)*y(3))*epsilon_r
c$$$
c$$$       ! d(Theta)/dt  
c$$$        ydot(1) = epsinvbyq*y(2)*gammaRi*Brem      
c$$$
c$$$       ! d(pr)/dt    
c$$$        ydot(2) =-epsinvbyq*rratio*sin(y(1))*muR*gammaRi*tmp*tmp  
c$$$     .           +Ersource + ALDrsink 
c$$$               
c$$$       ! d(pp)/dt
c$$$        ydot(3) = 0.5d0*y(3)*rratio*sin(y(1))*tmp*ydot(1)
c$$$     .           + ALDpsink
c$$$
c$$$      else   ! iEorbit false 
c$$$      ! electric field and radiation lumped
c$$$      ! with collisional sources
c$$$
c$$$        if (linBfield.eqv..false.) then
c$$$     
c$$$       ! d(Theta)/dt  
c$$$        ydot(1) = epsinvbyq*y(2)*gammaRi*Brem      
c$$$
c$$$       ! d(pr)/dt    
c$$$        ydot(2) =-epsinvbyq*rratio*sin(y(1))*muR*gammaRi*tmp*tmp
c$$$
c$$$       ! d(pp)/dt
c$$$        ydot(3) = 0.5d0*y(3)*rratio*sin(y(1))*tmp*ydot(1)
c$$$      
c$$$        else  
c$$$
c$$$        ydot(1) = epsinvbyq*y(2)*gammaRi     
c$$$
c$$$       ! d(pr)/dt    
c$$$        ydot(2) = -epsinvbyq*rratio
c$$$     .              *sin(y(1))*0.5d0*y(3)**2*gammaRi
c$$$
c$$$       ! d(pp)/dt
c$$$        ydot(3) = epsinvbyq*rratio
c$$$     .            *sin(y(1))*0.5d0*y(3)*y(2)*gammaRi 
c$$$      
c$$$        endif
c$$$      endif
c$$$  
c$$$      ydot(4) = 1d0     !tkernel_local(-y(5),dt_e,exl1)
c$$$
c$$$      ydot(5) = 1d0     !skernel_local(-y(5),dt_e,exl1)
c$$$
c$$$c Call additional integrators
c$$$
c$$$      x2=y(1)
c$$$      y2=y(2) 
c$$$      z2=y(3)
c$$$
c$$$      call chk_pos(x2,y2,z2,ierr=ierror,warning=.false.) 
c$$$      
c$$$      if (ierror /= ORB_OK) then
c$$$         flds(1:3)=0d0
c$$$c         print *,"outside domain",x2,y2,z2
c$$$      else
c$$$         if (use_local_splines) then
c$$$           !flds=0d0
c$$$           if (sp_order.eq.1) then
c$$$             call evalFldsLocalLinear(x2,y2,z2,sp_coef,flds)
c$$$c            print *,my_rank,'exact',x2*y2*z2          !exp(-0.5d0*(y2*y2+z2*z2)/1d-2)
c$$$c            print *,my_rank,'interpolated-2',flds(1) 
c$$$c            print *,x2,y2,z2
c$$$                
c$$$           else  
c$$$            call evalFldsLocalQuadratic(x2,y2,z2,sp_coef,flds)
c$$$            if (flds(1).lt.0d0 .or. flds(2).lt.0d0) then 
c$$$             call evalFldsLocalLinear(x2,y2,z2,sp_coef,flds)
c$$$            endif
c$$$
c$$$c            print *,my_rank,'exact',y2**2*x2**2*z2**2   
c$$$c            print *,my_rank,'interpolated-2',flds(1) 
c$$$c            print *,x2,y2,z2
c$$$                  
c$$$
c$$$           endif 
c$$$         else
c$$$           call evalFlds(x2,y2,z2,sp_coef,flds)
c$$$c           print *,my_rank,'exact',x2**2*y2*z2 +y2**2*x2**2*z2**2   
c$$$c           print *,my_rank,'exact',exp(-0.5d0*(y2*y2+z2*z2)/1d-2)
c$$$c           print *,my_rank,'interpolated-2',flds(1) 
c$$$  
c$$$         endif
c$$$      endif
c$$$
c$$$      !Enforce positive semi-definiteness of interpolation
c$$$c$$$      flds(3) = max(flds(3)-cn*flds(1)-cnm*flds(2),0d0) !Limit cnp*Tnp (from source)
c$$$      if (.not.use_local_splines) then 
c$$$        flds(1) = max(flds(1),0d0)                        !Limit Tn
c$$$        flds(2) = max(flds(2),0d0)                        !Limit Tnm
c$$$      endif
c$$$ 
c$$$      !Tn
c$$$      ydot(6) = flds(1)*tkernel_3d(y(5),dt_e,exl1)
c$$$
c$$$      !Tnm
c$$$      if (sm_flag.eq.0) then 
c$$$        ydot(7) =0d0
c$$$        exl2=.true.
c$$$      else
c$$$        ydot(7) = flds(2)*tkernel_3d(y(5),dt2_e,exl2)
c$$$      endif
c$$$
c$$$      !Src @ dt
c$$$c$$$      flds(3) = flds(3)+cn*flds(1)+cnm*flds(2)        !Recover source
c$$$      ydot(8) = flds(3)*skernel_3d(y(5),dt_e,exl3)
c$$$
c$$$      !Src @ dt+dtold
c$$$      if (sm_flag.eq.0) then 
c$$$         ydot(9)=0d0
c$$$         exl4=.true.
c$$$      else
c$$$         ydot(9) = flds(3)*skernel_3d(y(5),dt2_e,exl4)
c$$$      endif
c$$$
c$$$      !Transform back to arc-length integral
c$$$c      ydot(6:neq) = ydot(6:neq)*ydot(5)
c$$$
c$$$c Are we done?
c$$$
c$$$      call set_orbit_exit(exl1.and.exl2.and.exl3.and.exl4)
c$$$
c$$$c End program
c$$$
c$$$      end subroutine lsode_3d_ext
      
      end module re
