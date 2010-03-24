c postProcessSol
c####################################################################
      subroutine postProcessSol(v_arr)

c--------------------------------------------------------------------
c     Postprocess solution (e.g. filtering). In this instance, we
c     do divergence cleaning of the magnetic field.
c--------------------------------------------------------------------

      use nlfunction_setup

      use imposeBCinterface

      use precond_variables, only: pc_B_solve,pc_type,pc_divclean

      use xdraw_io

      implicit none

c Call variables

      type(var_array),pointer :: v_arr

c Local variables

      integer :: i,j,k,igx,igy,igz,nx,ny,nz
      real(8) :: mag
      real(8),allocatable,dimension(:,:,:,:) :: db_cnv

c Interfaces

#if !defined(vec_pot) && !defined(samrai)
      INTERFACE
        subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine imposeBoundaryConditions
      END INTERFACE
#endif

c Begin program

      !Setup grid info
      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Postprocess adiabatic equation for density update

      if (adiabatic .and. (.not.solve_rho)) then
        if (my_rank == 0) write (*,*) 'Postprocessing density'
        v_arr%array_var(IRHO)%array(1:nx,1:ny,1:nz) =
     .       u_0  %array_var(IRHO)%array(1:nx,1:ny,1:nz)
     .     *(v_arr%array_var(ITMP)%array(1:nx,1:ny,1:nz)
     .      /u_0  %array_var(ITMP)%array(1:nx,1:ny,1:nz))**(1./gamma)
      endif

c Divergence-clean magnetic field when PC is on

#if !defined(vec_pot) && !defined(samrai)

      if (     (.not.pc_B_solve)
     .    .or.  pc_divclean
     .    .or. (pc_type == 'id')
     .    .or. (nxd==1.and.nyd==1)
     .    .or. (nxd==1.and.nzd==1)
     .    .or. (nyd==1.and.nzd==1)) return

      if (ilevel > 0 .and. my_rank == 0) then
        write (*,*)
        write (*,*) 'Performing divergence cleaning in postproc...'
      endif

c Find db=B^(n+1)-B^n

      allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3))

      db_cnv(:,:,:,1) = v_arr%array_var(IBX)%array
     .                 -u_n  %array_var(IBX)%array
      db_cnv(:,:,:,2) = v_arr%array_var(IBY)%array
     .                 -u_n  %array_var(IBY)%array
      db_cnv(:,:,:,3) = v_arr%array_var(IBZ)%array
     .                 -u_n  %array_var(IBZ)%array

c Clean db

cc      call grad_div_clean(db_cnv)

      call vlap_div_clean(db_cnv)

cc      call curl_div_clean(db_cnv)

c Find new B^(n+1)

      v_arr%array_var(IBX)%array = db_cnv(:,:,:,1)
     .                +u_n%array_var(IBX)%array
      v_arr%array_var(IBY)%array = db_cnv(:,:,:,2)
     .                +u_n%array_var(IBY)%array
      v_arr%array_var(IBZ)%array = db_cnv(:,:,:,3)
     .                +u_n%array_var(IBZ)%array

c End program

      deallocate(db_cnv)

      contains

c     grad_div_clean
c     #################################################################
      subroutine grad_div_clean(db)

c     -----------------------------------------------------------------
c     Finds dB=dB* + grad(phi).
c     -----------------------------------------------------------------

      use matvec, rhs_bc => rhs, bcs_pc => bcs

      implicit none

c     Call variables

      real(8) :: db(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      integer :: ii,bcs(6)
      real(8) :: rhs(ntotdp),x(ntotdp),phi(0:nx+1,0:ny+1,0:nz+1)
     .          ,mag,dvol

      logical :: vwgt

      external lap_mtvc

c     Begin program

      vwgt = vol_wgt

c     Find rhs=-div(db)

      rhs = 0d0

      do k=1,nz
        do j=1,ny
          do i=1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            if (vwgt) then
              dvol = gmetric%grid(igx)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            rhs(ii)= -dvol*div(i,j,k,nx,ny,nz,igx,igy,igz
     $                   ,db(:,:,:,1)
     $                   ,db(:,:,:,2)
     $                   ,db(:,:,:,3))
          enddo
        enddo
      enddo

c diag ***** DIAGONAL GENERATION FOR XDRAW PLOTTING
cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .    ,status='replace')
cc      call contour(rhs,nx,ny,0d0,1d0,0d0,1d0,0,110)
c diag *****

c     Solve for phi=lap^-1(rhs)

      bcs = bcond
      where (bcond == DEF) bcs = NEU  !B.n must be zero at boundary!!

      x = 0d0

      !Solve Poisson equation
      call cSolver(1,ntotdp,rhs,x,bcs
     .            ,1,ilevel,0,lap_mtvc,vwgt
     .            ,tol           = 1d-5
!     .            ,dg            = lap_diag
     .            ,gm_driver     = .not.vwgt
     .            ,cg_driver     = vwgt
     .            ,ks_it         = 20
     .            ,mg_debug      = .false.
     .            ,mg_galerkin   = .false.
     .            ,mg_order_res  = 0  !To avoid null spaces
     .            ,mg_order_prol = 0  !To avoid null spaces
     .            ,mg_coarse_grid_size = 2
     .            ,mg_gm_coarse_solve  = .false.
     .            ,mg_vcyc       = 1
     .            ,mg_smooth     = 'gs'
     .            ,sm_omega      = 1d0
     .            ,sm_it         = 3
     .            ,sm_ncolors    = sm_ncolors  !9-pt stencil
     .            ,sm_line_relax = mg_line_relax
     .            ,sm_zebra_relax= mg_zebra_relax
     .            )

      !Map solution to array and apply BCs
      call mapMGVectorToArray(0,1,x,nx,ny,nz,phi,igx,.false.)

      call setMGBC(0,1,nx,ny,nz,igx,phi,bcs,iorder=bc_order)

c diag ***** DIAGONAL GENERATION FOR XDRAW PLOTTING
cc      call contour(phi(1:nx,1:ny,1),nx,ny,0d0,1d0,0d0,1d0,1,110)
cc      close(110)
cc      pause
cc      stop
c diag *****

c     Find dB=-grad(phi) (in cnv representation and w/o ghost cells)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            db(i,j,k,:) = db(i,j,k,:)
     .           + matmul(gmetric%grid(igx)%gsup(i,j,k,:,:)
     .                   ,grad2(i,j,k,nx,ny,nz,igx,igy,igz,phi,0,0,0))
          enddo
        enddo
      enddo

      end subroutine grad_div_clean

c     curl_div_clean
c     #################################################################
      subroutine curl_div_clean(db)

c     -----------------------------------------------------------------
c     Finds B^(n+1)=B^n + curl(dA), with dA=inv_curl(dB)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: db(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      real(8) :: da(0:nx+1,0:ny+1,0:nz+1,3)

c     Find dA from dB (w/ ghost cells)

      call curl_inv(nx,ny,nz,igx,igy,igz,db,da)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            db(i,j,k,:)=curl(i,j,k,nx,ny,nz,igx,igy,igz,da)
          enddo
        enddo
      enddo

      end subroutine curl_div_clean

c     vlap_div_clean
c     #################################################################
      subroutine vlap_div_clean(db)

c     -----------------------------------------------------------------
c     Finds curl(dA)=B^(n+1)-B^n, with dA=inv_lap(dJ). 
c     -----------------------------------------------------------------

      use matvec, rhs_bc => rhs, bcs_pc => bcs

      implicit none

c     Call variables

      real(8) :: db(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      integer :: ii,icomp,bcsa(6,3),bcsb(6,3)

      real(8) :: da     (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,da_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,db_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,rhs(ntotdp,3),x(ntotdp,3),dvol

      logical :: vwgt

      external lap_mtvc,avlap_mtvc

c     Begin program

      vwgt = .true.

c     Find rhs of vec_lap solve: rhs=vol*dJ (with dJ=curl(dB); dB already has BCs)

      bcsb(:,1) = u_0%array_var(IBX)%bconds
      bcsb(:,2) = u_0%array_var(IBY)%bconds
      bcsb(:,3) = u_0%array_var(IBZ)%bconds

      db_cov = db_cnv
      call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,db_cov(:,:,:,1)
     .                    ,db_cov(:,:,:,2)
     .                    ,db_cov(:,:,:,3)
     .                    ,'cnv','cov')

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            if (vwgt) then
              dvol = gmetric%grid(igx)%dvol(i,j,k)
            else
              dvol = 1d0
            endif

            rhs(ii,:)=-dvol*curl(i,j,k,nx,ny,nz,igx,igy,igz,db_cov)
          enddo
        enddo
      enddo

c     Define BCs for dA

      bcsa = bcsb
      where (bcsa == -NEU) bcsa = -DIR  !Use covariant components for tangential BCs

c     Allocate diagonal matrix

      if (.not.allocated(vlap_diag)) then
        allocate(vlap_diag(3,3*ntotd2p))

        call find_mf_diag(3,1,3*ntotdp,avlap_mtvc,igx,bcsa
     .                   ,vlap_diag,sm_ncolors,vwgt,.false.)
      endif

c     Invert vector laplacian to find dA=-vlap^-1(dJ)

      x = 0d0

      !Solve Poisson equation
      call cSolver(3,ntotdp,rhs,x,bcsa,igx
     .            ,max(0,ilevel-1),1,avlap_mtvc,vwgt
     .            ,dg            = vlap_diag
     .            ,gm_sing_matrix= .false.
     .            ,ks_it         = 20
     .            ,mg_order_res  = 0  !To ensure conservation (avoid null space)
     .            ,mg_order_prol = 1
     .            ,mg_coarse_grid_size = mg_coarse_size
     .            ,mg_gm_coarse_solve  = mg_gm_coarse
     .            ,mg_smooth     = sm_type
     .            ,sm_omega      = sm_omega
     .            ,sm_ncolors    = sm_ncolors
     .            ,sm_line_relax = mg_line_relax
     .            !Original parameters
     .            ,tol           = rtol
     .            ,gm_driver     = .true.
     .            ,mg_vcyc       = mg_vcyc
     .            ,sm_it         = sm_iter
     .            ,sm_zebra_relax= mg_zebra_relax
     .            !MG convergence test
cc     .            ,cvrg_tst      = .true.
cc     .            ,gm_driver     = .false.
cc     .            ,mg_vcyc       = 20
cc     .            ,tol           = 1d-10
cc     .            ,sm_it         = 3
cc     .            ,sm_zebra_relax= .false.
     .            )

c     Find new dB (with BCs)

cc      call fillArray(igx,IAX,3,x,da,bcsa,bc_order
      call fillArray(igx,IAX,3,x,da,bcsa,0
     .              ,is_vec=.true.,is_cnv=.true.,arr_cov=da_cov)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            db(i,j,k,:)=curl(i,j,k,nx,ny,nz,igx,igy,igz,da_cov)
          enddo
        enddo
      enddo

      call setMGBC(0,3,nx,ny,nz,igx,db,bcsb
     .            ,icomp=IBX,is_vec=.true.,is_cnv=.true.
     .            ,iorder=bc_order)

c diag ***** DIAGONAL GENERATION FOR XDRAW PLOTTING
cc      write (*,*) 'DIAG -- postProcessSol'
cc
cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .    ,status='replace')
cc      call contour(da_cov(:,:,1,1),nx+2,ny+2,0d0,1d0,0d0,1d0,0,110)
cc      call contour(da_cov(:,:,1,2),nx+2,ny+2,0d0,1d0,0d0,1d0,1,110)
cc      call contour(da_cov(:,:,1,3),nx+2,ny+2,0d0,1d0,0d0,1d0,1,110)
cc      call contour(db    (:,:,1,1),nx+2,ny+2,0d0,1d0,0d0,1d0,1,110)
cc      call contour(db    (:,:,1,2),nx+2,ny+2,0d0,1d0,0d0,1d0,1,110)
cc      call contour(db    (:,:,1,3),nx+2,ny+2,0d0,1d0,0d0,1d0,1,110)
cc      close (110)
cc      stop
c diag *****

c     Store helical flux (for resistivity)

      psi = psi + da_cov(:,:,:,3)

      end subroutine vlap_div_clean

#endif

      end subroutine postProcessSol
