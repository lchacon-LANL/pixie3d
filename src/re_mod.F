      module re

      use xdraw_io

      use pixie3d_io
      
      use nlfunction_setup

c$$$      use grid
c$$$      
c$$$      use parameters
c$$$
c$$$      use variables
c$$$
c$$$      use mg_solver

      use par_int

      use timestepping

      integer,private :: order_bc=1

      real(8),pointer,dimension(:,:,:),private :: nre_nm=>null()
     $                                           ,nre_n =>null()
     $                                           ,nre_np=>null()
     $                                           ,nre_0 =>null()
      real(8),pointer,dimension(:,:,:,:),private :: j__re=>null()

      logical :: re_diag=.true.

      real(8),private :: c__va
      
      !Setup lagrangian integration
      logical,private :: do_3d=.true.
     $                  ,local_interp=.true.
      real(8),private :: G_tol=1e-10      !Green's function tol
     $                  ,Dstar=0.15d0     !For stability
      integer,private :: sp_order=2       !Spline order of interpolation

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

      allocate(nre_nm(0:nx+1,0:ny+1,0:nz+1)
     $        ,nre_n (0:nx+1,0:ny+1,0:nz+1)
     $        ,nre_np(0:nx+1,0:ny+1,0:nz+1)
     $        ,nre_0 (0:nx+1,0:ny+1,0:nz+1)
     $        ,j__re (0:nx+1,0:ny+1,0:nz+1,3))

      nre_nm = 0d0 ; nre_n = 0d0 ; nre_np = 0d0 ; nre_0 = 0d0

      j__re = 0d0
      
      end subroutine re_allocVariables

c     re_deallocVariables
c     ###############################################################
      subroutine re_deallocVariables

      implicit none

      deallocate(nre_nm
     $          ,nre_n 
     $          ,nre_np
     $          ,nre_0
     $          ,j__re)
      
      end subroutine re_deallocVariables

c     init_re
c     ###############################################################
      subroutine init_re(nre,jre)

        implicit none

        real(8) :: nre(0:,0:,0:),jre(0:,0:,0:,:)

        nre_nm = nre
        nre_n  = nre
        nre_np = nre
        nre_0  = nre

        if (re_diag) then
          j__re = jre
          call dump_re_diag(.true.)
        endif
        
      end subroutine init_re

c     get_nre
c     ###############################################################
      subroutine get_nre(nre)

        implicit none

        real(8) :: nre(0:,0:,0:)

        nre = nre_np
 
      end subroutine get_nre

c     re_updt
c     ###############################################################
      subroutine re_updt(c_va,vnm,vn,vnp,Ecov,iout)

        implicit none

        !Call variables
        type(var_array),pointer :: vnm,vn,vnp
        integer :: iout
        real(8) :: c_va

        real(8),pointer,dimension(:,:,:,:) :: Ecov

        if (iout > 0 .and. my_rank == 0) then
          write (*,*)
          write (*,*) 'Performing RE integration in postproc...'
        endif
 
        !Eulerian update for nre_np
        call re_Eul(c_va,vnm,vn,vnp,Ecov,iout,.false.)

        !Lagrangian update for nre_np
        call re_lag(c_va,vnm,vn,vnp,iout)

      end subroutine re_updt

c     re_Eul
c     ###############################################################
      subroutine re_Eul(c_va,v_nm,v_n,v_np,Ecov,iout,full_flow)

        implicit none

c       Call variables

        type(var_array),pointer :: v_nm,v_n,v_np
        integer :: iout
        real(8) :: c_va

        real(8),pointer,dimension(:,:,:,:) :: Ecov
        
        logical :: full_flow
        
c       Local variables

        real(8) :: dvol
        
        integer :: nn,nx,ny,nz,igrid,bcs(6,3),bc_nre(6),guess,order
     $            ,i,j,k,ii

        real(8),allocatable,dimension(:,:) :: x,y
        real(8),allocatable,dimension(:,:,:) :: ib2
        real(8),allocatable,dimension(:,:,:,:) :: b_cnv,b_cov

c       Begin program

        if (iout > 0 .and. my_rank == 0) then
          write (*,*)
          call info('  Eulerian RE integration...')
        endif

        igrid = 1

        nx = gv%gparams%nxv(igrid)
        ny = gv%gparams%nyv(igrid)
        nz = gv%gparams%nzv(igrid)

        nn = nx*ny*nz

        call allocateMGArray(gv%gparams,3,gv_exb)

        allocate(b_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     $          ,b_cov(0:nx+1,0:ny+1,0:nz+1,3)
     $          ,ib2  (0:nx+1,0:ny+1,0:nz+1)  )
        
        allocate(x(nn,1),y(nn,1))

       !Save old timestep solutions
        nre_nm = nre_n
        nre_n  = nre_np
      
c       Prepare solver coefficients

        b_cnv(:,:,:,1) = v_np%array_var(IBX)%array
        b_cnv(:,:,:,2) = v_np%array_var(IBY)%array
        b_cnv(:,:,:,3) = v_np%array_var(IBZ)%array

        b_cov = XformToCov(gv%gparams,igrid,b_cnv)

        ib2 = 1d0/vectornorm(gv%gparams,igrid,b_cnv,.false.)

        gv_exb%grid(igrid)%array =
     $       crossproduct(gv%gparams,igrid,Ecov,b_cov,.false.)
        gv_exb%grid(igrid)%array(:,:,:,1) =
     $       gv_exb%grid(igrid)%array(:,:,:,1)*ib2
        gv_exb%grid(igrid)%array(:,:,:,2) =
     $       gv_exb%grid(igrid)%array(:,:,:,2)*ib2
        gv_exb%grid(igrid)%array(:,:,:,3) =
     $       gv_exb%grid(igrid)%array(:,:,:,3)*ib2

        if (full_flow) then
          ib2 = sqrt(ib2)
          gv_exb%grid(igrid)%array(:,:,:,1) =
     $       gv_exb%grid(igrid)%array(:,:,:,1)
     $      +c_va*b_cnv(:,:,:,1)*ib2
          gv_exb%grid(igrid)%array(:,:,:,2) =
     $       gv_exb%grid(igrid)%array(:,:,:,2)
     $      +c_va*b_cnv(:,:,:,2)*ib2
          gv_exb%grid(igrid)%array(:,:,:,3) =
     $       gv_exb%grid(igrid)%array(:,:,:,3)
     $      +c_va*b_cnv(:,:,:,3)*ib2
        endif
        
        bcs(:,1) = v_np%array_var(IVX)%bconds
        bcs(:,2) = v_np%array_var(IVY)%bconds
        bcs(:,3) = v_np%array_var(IVZ)%bconds

        order = 0
        call restrictMGArray(IVX,3,gv_exb,bcs,igrid,order
     .                      ,order_bc=order_bc)

c       Solve n_re system for new-time n_re

        bc_nre = v_np%array_var(IRHO)%bconds

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              y(ii,1) = -(cn*nre_n(i,j,k) + cnm*nre_nm(i,j,k))/dt
     $         - cnfactor*flx_advec(gv%gparams,i,j,k,nx,ny,nz
     $                       ,igrid,igrid,igrid,gv_exb%grid(igrid)%array
     $                       ,nre_n,UPWD)

              dvol = gv%gparams%gmetric%grid(igrid)%jac(i,j,k)
              if (vol_wgt) then
                dvol = dvol*gv%gparams%gmetric%grid(igrid)%dlvol(i,j,k)
              endif
              
              y(ii,1) = y(ii,1)*dvol

              x(ii,1) = nre_n(i,j,k)
            enddo
          enddo
        enddo

        guess = 1
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

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              nre_np(i,j,k) = x(ii,1)
            enddo
          enddo
        enddo

        !Impose BCs on perturbation
        nre_np = nre_np - nre_0
        call setBC_ez(gv%gparams,igrid,nre_np,bc_nre,iorder=order_bc)
        nre_np = nre_np + nre_0

c       Diagnostics

        if (my_rank == 0.and.iout > 0) call info('  Done!')

        if (TS_chk_dump().and.re_diag) then
          j__re(:,:,:,1) = gv_exb%grid(igrid)%array(:,:,:,1)*nre_np
          j__re(:,:,:,2) = gv_exb%grid(igrid)%array(:,:,:,2)*nre_np
          j__re(:,:,:,3) = gv_exb%grid(igrid)%array(:,:,:,3)*nre_np
          call dump_re_diag(.false.)
        endif
                
c       Garbage collection
        
        deallocate(x,y,b_cnv,b_cov,ib2)
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
              dvol = g_def%gmetric%grid(igrid)%dlvol(i,j,k)*jac(i,j,k)
            else
              dvol = jac(i,j,k)
            endif

            !Conservative advection
            upwind = flx_advec(g_def,i,j,k,nnx,nny,nnz,igrid,igrid,igrid
     .                        ,v0_cnv,drho(:,:,:,1),UPWD)

c$$$            ! Diffusion
c$$$            if (dd > 0d0) then
c$$$              slap = dd*lap(g_def,i,j,k,nnx,nny,nnz,igrid,igrid,igrid
c$$$     .                     ,drho(:,:,:,1))
c$$$            else
              slap = 0d0
c$$$            endif

            y(ijk)=dvol*(cnp/dt*drho(i,j,k,1) + alpha*(upwind - slap))

          enddo
        enddo
      enddo

c     End program

      deallocate(drho)
      nullify(v0_cnv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine nre_mtvc

c     re_Lag
c     ###############################################################
      subroutine re_Lag(c_va,vnm,vn,vnp,iout)

        implicit none

        !Call variables
        type(var_array),pointer :: vnm,vn,vnp
        integer :: iout
        real(8) :: c_va

        !Local variables
        integer :: igrid,nx,ny,nz

        real(8) :: rkernel,renum_before,renum_after
        
        real(8),allocatable,dimension(:,:,:) :: tmp_nm,tmp_n,tmp_np,ieps
     $                                         ,bmag_nm,bmag_n,bmag_np
        real(8),allocatable,dimension(:,:,:,:) :: bcnv_nm,bcnv_n,bcnv_np

        !Begin program
        
        igrid = 1

        nx = gv%gparams%nxv(igrid)
        ny = gv%gparams%nyv(igrid)
        nz = gv%gparams%nzv(igrid)

        allocate(ieps   (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_nm (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_n  (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_np (0:nx+1,0:ny+1,0:nz+1)
     $          ,bmag_n (0:nx+1,0:ny+1,0:nz+1)
     $          ,bmag_nm(0:nx+1,0:ny+1,0:nz+1)
     $          ,bmag_np(0:nx+1,0:ny+1,0:nz+1)
     $          ,bcnv_n (0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bcnv_nm(0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bcnv_np(0:nx+1,0:ny+1,0:nz+1,3))
        
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
        
        !Initialize Lagrangian step (!advection-diff mode with diff coeff Dstar)
        rkernel = -1d0/sqrt(Dstar) 
        call lag_init(gv%gparams,igrid,bcnv_np,rkernel,iout
     .               ,G_tol,sp_order,.not.local_interp,.false.,do_3d)

        !Find temp from pressure (automatically takes care of jacobian factor)
        tmp_nm = nre_nm/bmag_nm
        tmp_n  = nre_n /bmag_n
        tmp_np = nre_np/bmag_np

        ieps = c_va
        
        call lag_update_bdf2(gv%gparams,igrid,IRHO
     .                      ,vnp%array_var(IRHO)%bconds
     $                      ,tmp_nm,tmp_n,tmp_np,ieps,iout)

        !Destroy Lagrangian step
        call lag_destroy

        !Return RE density
        nre_np = tmp_np*bmag_np
        
        !Check and correct for particle conservation?
c$$$        renum_before = integral(gv%gparams,igrid,nre_n )
c$$$        renum_after  = integral(gv%gparams,igrid,nre_np)
c$$$        nre_np = nre_np*renum_before/renum_after
        
        !Impose BCs on perturbation
        nre_np = nre_np - nre_0
        call setBC_ez(gv%gparams,igrid,nre_np,vnp%array_var(IRHO)%bconds
     $               ,iorder=order_bc)
        nre_np = nre_np + nre_0

        !Garbage collection
        deallocate(tmp_nm,tmp_n,tmp_np,bmag_nm,bmag_n,bmag_np
     $            ,bcnv_n,bcnv_nm,bcnv_np,ieps)
        
      end subroutine re_Lag

c     re_Lag_liu
c     ###############################################################
      subroutine re_Lag_liu(c_va,vnm,vn,vnp,Ecov,iout)

        implicit none

        !Call variables
        type(var_array),pointer :: vnm,vn,vnp
        integer :: iout
        real(8) :: c_va

        real(8),pointer,dimension(:,:,:,:) :: Ecov
        
        !Local variables
        integer :: igrid,nx,ny,nz,bcnd(6,3),bc_nre(6)

        real(8) :: rkernel,renum_before,renum_after
        
        real(8),allocatable,dimension(:,:,:) :: tmp_nm,tmp_n,tmp_np,ieps
     $                                        ,ibmag_nm,ibmag_n,ibmag_np
        real(8),allocatable,dimension(:,:,:,:) :: bb,bcnv_nm,bcnv_n
     $                                           ,bcnv_np,v_exb,j_cnv
     $                                           ,bcov_np

        !Begin program

        c__va = c_va

        if (iout > 0 .and. my_rank == 0) then
          write (*,*)
          call info('  Lagrangian RE integration...')
        endif

        igrid = 1

        nx = gv%gparams%nxv(igrid)
        ny = gv%gparams%nyv(igrid)
        nz = gv%gparams%nzv(igrid)

        bc_nre = vnp%array_var(IRHO)%bconds

       !Save old timestep solutions
        nre_nm = nre_n
        nre_n  = nre_np

        !Allocate variables
        allocate(ieps    (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_nm  (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_n   (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_np  (0:nx+1,0:ny+1,0:nz+1)
     $          ,ibmag_n (0:nx+1,0:ny+1,0:nz+1)
     $          ,ibmag_nm(0:nx+1,0:ny+1,0:nz+1)
     $          ,ibmag_np(0:nx+1,0:ny+1,0:nz+1)
     $          ,j_cnv   (0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bb      (0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bcnv_n  (0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bcnv_nm (0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bcnv_np (0:nx+1,0:ny+1,0:nz+1,3)
     $          ,bcov_np (0:nx+1,0:ny+1,0:nz+1,3))
        
        !Set orbit tols to Green's function tol
        call set_orbit_tols(G_tol,G_tol)

        !Find required vectors
        bcnv_nm(:,:,:,1) = vnm%array_var(IBX)%array
        bcnv_nm(:,:,:,2) = vnm%array_var(IBY)%array
        bcnv_nm(:,:,:,3) = vnm%array_var(IBZ)%array
        bcnv_n (:,:,:,1) = vn %array_var(IBX)%array
        bcnv_n (:,:,:,2) = vn %array_var(IBY)%array
        bcnv_n (:,:,:,3) = vn %array_var(IBZ)%array
        bcnv_np(:,:,:,1) = vnp%array_var(IBX)%array
        bcnv_np(:,:,:,2) = vnp%array_var(IBY)%array
        bcnv_np(:,:,:,3) = vnp%array_var(IBZ)%array
        
        ibmag_nm = vectorNorm(gv%gparams,igrid,bcnv_nm,.false.)
        ibmag_n  = vectorNorm(gv%gparams,igrid,bcnv_n ,.false.)
        ibmag_np = vectorNorm(gv%gparams,igrid,bcnv_np,.false.)

        ibmag_nm = 1d0/sqrt(ibmag_nm)
        ibmag_n  = 1d0/sqrt(ibmag_n )
        ibmag_np = 1d0/sqrt(ibmag_np)
        
        bb(:,:,:,1) = bcnv_np(:,:,:,1)*ibmag_np
        bb(:,:,:,2) = bcnv_np(:,:,:,2)*ibmag_np
        bb(:,:,:,3) = bcnv_np(:,:,:,3)*ibmag_np

        j_cnv = curl(gv%gparams,igrid,XformToCov(gv%gparams,igrid,bb))

        !Impose BCs on current
        bcnd = gv%aux%vec_list(IJCNV)%bconds
        where(    abs(bcnd) == DIR
     $       .or. abs(bcnd) == NEU
     $       .or. abs(bcnd) == EQU ) bcnd = EXT
        call setBC_ez(gv%gparams,igrid,j_cnv,.true.,bcnd,1)
c diag ****
c$$$        write (*,*) bcnd
c$$$        call createDrawInCfile(9,'jj.bin','j','t','x','y'
c$$$     $                        ,(/'jj_1','jj_2','jj_3'
c$$$     $                          ,'jj_1','jj_2','jj_3'
c$$$     $                          ,'jj_1','jj_2','jj_3'/)
c$$$     $                        ,'-c -X0 -L57','drawjj.in')
c$$$        open(unit=110,file='jj.bin',form='unformatted'
c$$$     .      ,status='unknown')
c$$$        call contour(j_cnv(0:nxd+1,0:nyd,0,1),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,0,110)
c$$$        call contour(j_cnv(0:nxd+1,0:nyd,0,2),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(j_cnv(0:nxd+1,0:nyd,0,3),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(j_cnv(0:nxd+1,0:nyd,1,1),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(j_cnv(0:nxd+1,0:nyd,1,2),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(j_cnv(0:nxd+1,0:nyd,1,3),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(j_cnv(0:nxd+1,0:nyd,2,1),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(j_cnv(0:nxd+1,0:nyd,2,2),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(j_cnv(0:nxd+1,0:nyd,2,3),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        close(110)
c$$$        stop
c diag ****

        !Find source
        nre_np =(Ecov(:,:,:,1)*j_cnv(:,:,:,1)
     $          +Ecov(:,:,:,2)*j_cnv(:,:,:,2)
     $          +Ecov(:,:,:,3)*j_cnv(:,:,:,3))
     $          /gv%gparams%gmetric%grid(igrid)%jac

        nre_np = nre_np*ibmag_np
        
        !Find ExB flow
        bcov_np = XformToCov(gv%gparams,igrid,bcnv_np)
        v_exb = crossproduct(gv%gparams,igrid,Ecov,bcov_np,.false.)
        v_exb(:,:,:,1) = v_exb(:,:,:,1)*ibmag_np*ibmag_np
        v_exb(:,:,:,2) = v_exb(:,:,:,2)*ibmag_np*ibmag_np
        v_exb(:,:,:,3) = v_exb(:,:,:,3)*ibmag_np*ibmag_np

        !Find orbit flow
        bb = bb + v_exb/c_va
c diag ****
c$$$        call createDrawInCfile(9,'bb.bin','bb','t','x','y'
c$$$     $     ,(/'bb_1','bb_2','bb_3'
c$$$     $       ,'bb_1','bb_2','bb_3'
c$$$     $       ,'bb_1','bb_2','bb_3'/)
c$$$     $     ,'-c -X0 -L57','drawbb.in')
c$$$        open(unit=110,file='bb.bin',form='unformatted'
c$$$     .      ,status='unknown')
c$$$        call contour(bb (0:nxd+1,0:nyd,0,1),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,0,110)
c$$$        call contour(bb (0:nxd+1,0:nyd,0,2),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(bb (0:nxd+1,0:nyd,0,3),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(bb (0:nxd+1,0:nyd,1,1),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(bb (0:nxd+1,0:nyd,1,2),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(bb (0:nxd+1,0:nyd,1,3),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(bb (0:nxd+1,0:nyd,2,1),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(bb (0:nxd+1,0:nyd,2,2),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(bb (0:nxd+1,0:nyd,2,3),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        close(110)
c$$$        stop
c diag ****
        
        !Initialize Lagrangian step (advection-diff mode with diff coeff Dstar)
        rkernel = -1d0/sqrt(Dstar) 
        call lag_init(gv%gparams,igrid,bb,rkernel,iout,G_tol,sp_order
     .               ,.not.local_interp,.false.,do_3d)

        !Find temp from pressure (automatically takes care of jacobian factor)
        tmp_nm = log(nre_nm*ibmag_nm)
        tmp_n  = log(nre_n *ibmag_n )
        tmp_np = (dt*nre_np - cn*tmp_n - cnm*tmp_nm)/cnp !Integrate source

c diag ****
c$$$        write (*,*) sum(tmp_nm**2),sum(tmp_n**2),sum(tmp_np**2)
c$$$        call createDrawInCfile(9,'tmp.bin','j_re','t','x','y'
c$$$     $     ,(/'tnm_1','tnn_2','tnp_3'
c$$$     $       ,'tnm_1','tnn_2','tnp_3'
c$$$     $       ,'tnm_1','tnn_2','tnp_3'/)
c$$$     $     ,'-c -X0 -L57','drawtmp.in')
c$$$        open(unit=110,file='tmp.bin',form='unformatted'
c$$$     .      ,status='unknown')
c$$$        call contour(tmp_nm(0:nxd+1,0:nyd,0),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,0,110)
c$$$        call contour(tmp_n (0:nxd+1,0:nyd,0),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(tmp_np(0:nxd+1,0:nyd,0),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(tmp_nm(0:nxd+1,0:nyd,1),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(tmp_n (0:nxd+1,0:nyd,1),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(tmp_np(0:nxd+1,0:nyd,1),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(tmp_nm(0:nxd+1,0:nyd,2),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(tmp_n (0:nxd+1,0:nyd,2),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        call contour(tmp_np(0:nxd+1,0:nyd,2),nxd+2,nyd+1
c$$$     .              ,0d0,xmax,0d0,ymax,1,110)
c$$$        close(110)
c$$$        stop
c diag ****
        ieps = 1d0 !c_va  !Not used, but needs to be <> 0
        
        call lag_update_bdf2(gv%gparams,igrid,IRHO,bc_nre
     .                      ,tmp_nm,tmp_n,tmp_np,ieps,iout
     .                      ,lsode_ext=lsode_3d_ext)

        !Destroy Lagrangian step
        call lag_destroy

        !Compute new-time RE density
        nre_np = exp(tmp_np)/ibmag_np
        
        !Check and correct for particle conservation?
c$$$        renum_before = integral(gv%gparams,igrid,nre_n )
c$$$        renum_after  = integral(gv%gparams,igrid,nre_np)
c$$$        nre_np = nre_np*renum_before/renum_after
        
        !Impose BCs on perturbation
        nre_np = nre_np - nre_0
        call setBC_ez(gv%gparams,igrid,nre_np,bc_nre,iorder=order_bc)
        nre_np = nre_np + nre_0

        !Diagnostics
        if (my_rank == 0.and.iout > 0) call info('  Done!')

        if (TS_chk_dump().and.re_diag) then
          j__re(:,:,:,1) = c_va*bb(:,:,:,1)*nre_np
          j__re(:,:,:,2) = c_va*bb(:,:,:,2)*nre_np
          j__re(:,:,:,3) = c_va*bb(:,:,:,3)*nre_np
          call dump_re_diag(.false.)
        endif
        
        !Garbage collection
        deallocate(tmp_nm,tmp_n,tmp_np,bb,ibmag_nm,ibmag_n,ibmag_np
     $            ,bcnv_n,bcnv_nm,bcnv_np,ieps,j_cnv,bcov_np)

      end subroutine re_Lag_liu
      
c     lsode_3d_ext
c     ##################################################################
      subroutine lsode_3d_ext(neq,t,y,ydot)

c     ------------------------------------------------------------------
c     Driver routine for LSODE integration: ydot = ff(y,t)
c     ------------------------------------------------------------------

      use par_int

      implicit none

c     Call variables

      integer :: neq
      real(8) :: t,y(neq),ydot(neq)

c     Local variables

      integer :: ierror
      real(8) :: x1,y1,z1,flds(size(sp_coef,4)),B2,iB2,i_eps
      logical :: exl1,exl2,exl3,exl4

c     Begin program

      if (neq /= 9) then
        call pstop('lsode_3d_bdf2','Wrong number of equations')
      endif

c     Call LSODE orbit evaluation function (extracting logical coords x1,y1,z1)

      call lsode_orbit_f(5,t,y(1:5),ydot(1:5),xx1=x1,yy1=y1,zz1=z1)

c     Call additional integrators

      i_eps = c__va
      
      call evalFlds(x1,y1,z1,sp_coef,flds)

      !Tn
      ydot(6) = flds(1)*tkernel_3d(y(5),i_eps*dt,exl1)

      !Tnm
      ydot(7) = flds(2)*tkernel_3d(y(5),i_eps*(dt+dtold),exl2)

      !Src @ dt
      ydot(8) = flds(3)*skernel_3d(y(5),i_eps*dt,exl3)

      !Src @ dt+dtold
      ydot(9) = flds(3)*skernel_3d(y(5),i_eps*(dt+dtold),exl4)

      !Transform back to arc-length integral (only accounting for Jacobian, not |B|)
      ydot(6:neq) = ydot(6:neq)*ydot(4)

c     Are we done?

      call set_orbit_exit(exl1.and.exl2.and.exl3.and.exl4)

c     End program

      end subroutine lsode_3d_ext

c     dump_re_diag
c     ######################################################################
      subroutine dump_re_diag(init)

      implicit none

      !Call variables
      logical :: init

      !Local variables
      integer :: ierr,efu,i,j,funit
      real(8) :: dh

      real(8),pointer,dimension(:,:,:) :: nre_g,jre1,jre2,jre3
      real(8),pointer,dimension(:,:,:,:) :: jre_g

      type (graph_var_def),pointer,dimension(:) :: nre_graph

      type (cont_info) :: m_conf

      !Begin program
      
      allocate(nre_g(0:nxdp,0:nydp,0:nzdp))
      allocate(jre_g(0:nxdp,0:nydp,0:nzdp,3))
      allocate(jre1 (0:nxdp,0:nydp,0:nzdp))
      allocate(jre2 (0:nxdp,0:nydp,0:nzdp))
      allocate(jre3 (0:nxdp,0:nydp,0:nzdp))

      if (init) then
        call find_global(nre_0,nre_g)
        call find_global(j__re,jre_g)
      else
        call find_global(nre_np-nre_0,nre_g)
        call find_global(j__re       ,jre_g)
      endif

      jre1 = jre_g(:,:,:,1)
      jre2 = jre_g(:,:,:,2)
      jre3 = jre_g(:,:,:,3)
      
      if (my_rank == 0) then

        allocate(nre_graph(4))
        nre_graph(1)%descr = "nre"
        nre_graph(1)%array => nre_g
        nre_graph(2)%descr = "jre^1"
        nre_graph(2)%array => jre1
        nre_graph(3)%descr = "jre^2"
        nre_graph(3)%array => jre2
        nre_graph(4)%descr = "jre^3"
        nre_graph(4)%array => jre3

        m_conf%plane  = 3
        m_conf%coords = (/1,1,1/)

        funit = find_unit(110)

        if (init) then
          open(unit=funit,file='re.bin'
     .        ,form='unformatted',status='unknown')

          call createDrawInMfile(4,'re.bin',.false.,'RE qtys','t','x'
     $          ,'v',(/1,2,3,4/),nre_graph%descr,'-L50','drawre.in')

          allocate(xi(0:nxdp,0:nydp,1,3))
          dh = (xmax-xmin)/nxd
          do i=0,nxdp
            xi(i,:,:,1) = xmin + (i-0.5)*dh
          enddo
          dh = (ymax-ymin)/nyd
          do i=0,nydp
            xi(:,i,:,2) = ymin + (i-0.5)*dh
          enddo
          do i=1,1
            xi(:,:,i,3) = (i-0.5)
          enddo
          call M_contour(funit,4,nre_graph,0,m_conf,xi=xi
     $                  ,imax_g=nxd,imin_g=1
     $                  ,jmax_g=nyd,jmin_g=1
     $                  ,kmax_g=1  ,kmin_g=1)
          deallocate(xi)
        else
          open(unit=funit,file='re.bin'
     .         ,form='unformatted',status='old',access='append')
        endif

        call M_contour(funit,4,nre_graph,1,m_conf
     $                ,imax_g=nxd,imin_g=1
     $                ,jmax_g=nyd,jmin_g=1
     $                ,kmax_g=1  ,kmin_g=1)
        close(funit)

        deallocate(nre_graph)
      endif

c diag ****
c$$$      do i=0,nxdp
c$$$        write (*,*) 1d0*i,nre_g(i,10,1)
c$$$      enddo
c diag ****
      
      deallocate(nre_g,jre1,jre2,jre3,jre_g)

#if defined(petsc)
      call MPI_Barrier(MPI_COMM_WORLD,mpierr)
#endif

      end subroutine dump_re_diag
      
      end module re
