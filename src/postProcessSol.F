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

      use B_tools

      implicit none

c Call variables

      type(var_array),pointer :: v_arr

c Local variables

      integer :: i,j,k,igx,igy,igz,nx,ny,nz,ieq
      real(8) :: mag,acc
      real(8),allocatable,dimension(:,:,:,:) :: db_cnv,da_cov

c Interfaces

#if !defined(vec_pot) && !defined(samrai)
      INTERFACE
        subroutine applyBC (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine applyBC
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

cc      if (solve_prs .and. (.not.solve_rho)) then
cc        if (my_rank == 0) write (*,*) 'Postprocessing density'
cccc        v_arr%array_var(IRHO)%array(1:nx,1:ny,1:nz) =
cccc     .       u_0  %array_var(IRHO)%array(1:nx,1:ny,1:nz)
cccc     .     *(v_arr%array_var(ITMP)%array(1:nx,1:ny,1:nz)
cccc     .      /u_0  %array_var(ITMP)%array(1:nx,1:ny,1:nz))**(1./gamma)
cc        v_arr%array_var(IRHO)%array =
cc     .       u_0  %array_var(IRHO)%array
cc     .     *(v_arr%array_var(ITMP)%array
cc     .      /u_0  %array_var(ITMP)%array)**(1./gamma)
cc      endif

c Add random noise in EOM (to break symmetry)

      if (noise) then
        if (nc_eom_v) then
          do k=1,nz
            do j=1,ny
              do i=1,nx
                do ieq=IVX,IVZ
                  call random_number(acc) !Between 0 and 1
                  acc = noise_lev*sqrt(dt)*(1-2.*acc)
                  v_arr%array_var(ieq)%array=v_arr%array_var(ieq)%array
     .                                      +acc
                enddo
              enddo
            enddo
          enddo
        else
          do k=1,nz
            do j=1,ny
              do i=1,nx
                do ieq=IVX,IVZ
                  call random_number(acc) !Between 0 and 1
                  acc = noise_lev*sqrt(dt)*(1-2.*acc)
                  v_arr%array_var(ieq)%array=v_arr%array_var(ieq )%array
     .                                  +acc*v_arr%array_var(IRHO)%array
                enddo
              enddo
            enddo
          enddo
        endif
      endif

c Divergence-clean magnetic field when PC is on

#if !defined(vec_pot) && !defined(samrai)

      if (     (.not.pc_B_solve)
     .    .or.  pc_divclean
     .    .or. (.not.solenoidal)
     .    .or. (pc_type == 'id')
     .    .or. (nxd==1.and.nyd==1)
     .    .or. (nxd==1.and.nzd==1)
     .    .or. (nyd==1.and.nzd==1)) return

      if (ilevel > 0 .and. my_rank == 0) then
        write (*,*)
        write (*,*) 'Performing divergence cleaning in postproc...'
      endif

c Find db=B^(n+1)-B^n

      allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,da_cov(0:nx+1,0:ny+1,0:nz+1,3))

      db_cnv(:,:,:,1) = v_arr%array_var(IBX)%array
     .                 -u_n  %array_var(IBX)%array
      db_cnv(:,:,:,2) = v_arr%array_var(IBY)%array
     .                 -u_n  %array_var(IBY)%array
      db_cnv(:,:,:,3) = v_arr%array_var(IBZ)%array
     .                 -u_n  %array_var(IBZ)%array

c Clean db

cc      call grad_div_clean(db_cnv)

      call vlap_div_clean(ilevel,nx,ny,nz
     .                   ,igx,igy,igz,db_cnv,dacov=da_cov
     .                   )!,diag_plots=.true.)

cc      call curl_div_clean(db_cnv)

c Find new B^(n+1)

      v_arr%array_var(IBX)%array = db_cnv(:,:,:,1)
     .                +u_n%array_var(IBX)%array
      v_arr%array_var(IBY)%array = db_cnv(:,:,:,2)
     .                +u_n%array_var(IBY)%array
      v_arr%array_var(IBZ)%array = db_cnv(:,:,:,3)
     .                +u_n%array_var(IBZ)%array

c Update BCs of all quantities

      call applyBC (v_arr,igx,igy,igz)

c Update helical flux

      psi = psi + da_cov(:,:,:,3)

c End program

      deallocate(db_cnv,da_cov)

      contains

c     grad_div_clean
c     #################################################################
      subroutine grad_div_clean(db)

c     -----------------------------------------------------------------
c     Finds dB=dB* + grad(phi).
c     -----------------------------------------------------------------

      use matvec, bcs_pc => bcs

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

            rhs(ii)= -dvol*div(i,j,k,nx,ny,nz,igx,igy,igz,db)
          enddo
        enddo
      enddo

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

c     Find dB=-grad(phi) (in cnv representation and w/o ghost cells)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            db(i,j,k,:) = db(i,j,k,:)
     .           + matmul(gmetric%grid(igx)%gsup(i,j,k,:,:)
     .                   ,grad(i,j,k,nx,ny,nz,igx,igy,igz,phi))
          enddo
        enddo
      enddo

      end subroutine grad_div_clean

#endif

      end subroutine postProcessSol

c preProcessSol
c####################################################################
      subroutine preProcessSol(v1,v2)

c--------------------------------------------------------------------
c     Preprocess solution (place holder)
c--------------------------------------------------------------------

      use variables

      implicit none

c Call variables

      type(var_array),pointer :: v1,v2

c End program

      end subroutine preProcessSol
