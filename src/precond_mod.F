c  module matvec
c ###################################################################
      module matvec

        use grid

        use setMGBC_interface

        use mg_internal

        use constants

      end module matvec

c module precond_variables
c ######################################################################
      module precond_variables

        use timeStepping

        use precond_setup

        use transport_params

        use equilibrium

        use constants

        use auxiliaryVariables

        use operators

        use variables

        use grid

        use mgarraySetup

#if defined(petsc)
        integer(4) :: jit=1
#else
        use newtongm
#endif

        integer(4) :: i,j,k,ig,jg,kg,ieq
        real(8)    :: x1,y1,z1
        logical    :: cartsn,covariant,to_cartsn,to_cnv

        integer(4) :: ntotd2p,si_it

cc        real(8), allocatable, dimension(:,:,:,:):: p0

        real(8), allocatable, dimension(:,:)    :: rho_diag,tmp_diag
     .                                            ,b_diag,v_diag
     .                                            ,mstar_diag

        real(8), allocatable, dimension(:,:)    :: mgj0cnv,mgadvdiffV0

        real(8), allocatable, dimension(:,:,:)  :: mgnablaV0

        real(8), allocatable, dimension(:)      :: mgdivV0

        integer(4), allocatable, dimension(:,:) :: bcs

        type (var_array),target :: varray

        type (mg_array ),target :: gp0,gb0,gv0,grho0,gb0_cov

        real(8),pointer,dimension(:,:,:):: rho,px,py,pz,bx,by,bz,tmp

        logical :: form_diag=.true.,gm_smooth=.false.,mg_glrkin=.false.
     .            ,si_car=.true.,line_relax=.false.,rat_limit=.true.

cc        INTERFACE ASSIGNMENT (=)
cc          module procedure scatter,gather
cc        END INTERFACE

      contains

c     allocPrecVariables
c     ###################################################################
      subroutine allocPrecVariables

c     -------------------------------------------------------------------
c     Allocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer(4) :: alloc_stat

c     Begin program

        ntotd2p = 2*ntotd/neqd

        allocate (mgj0cnv    (ntotd2p,3))
        allocate (mgadvdiffV0(ntotd2p,3))
        allocate (mgnablaV0  (ntotd2p,3,3))
        allocate (mgdivV0    (ntotd2p))

        allocate (bcs(6,neqd+3))

        allocate (rho_diag  (1,  ntotd2p)
     .           ,tmp_diag  (1,  ntotd2p)
     .           ,  b_diag  (3,3*ntotd2p)
     .           ,  v_diag  (3,3*ntotd2p)
     .           ,mstar_diag(3,3*ntotd2p),STAT=alloc_stat)

        call allocateMGArray(1,gp0)
        call allocateMGArray(1,grho0)
        call allocateMGArray(3,gv0)
        call allocateMGArray(3,gb0)
        call allocateMGArray(3,gb0_cov)

c     End program

      end subroutine allocPrecVariables

c     deallocPrecVariables
c     ###################################################################
      subroutine deallocPrecVariables

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

c     Begin program

        deallocate (mgj0cnv)
        deallocate (mgdivV0)
        deallocate (mgadvdiffV0)
        deallocate (mgnablaV0)
        deallocate (bcs)

        call deallocateMGArray(gp0)
        call deallocateMGArray(grho0)
        call deallocateMGArray(gb0)
        call deallocateMGArray(gb0_cov)
        call deallocateMGArray(gv0)

        call deallocateDerivedType(varray)

c     End program

      end subroutine deallocPrecVariables

c     gather
c     ###################################################################
      subroutine gather(xout,xin)

c     -------------------------------------------------------------------
c     Gathers magnitudes from different arrays to vector.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(IN),dimension(:,:) :: xin
        real(8),intent(OUT),dimension(:)  :: xout

c     Local variables

        integer(4) :: i,j,k,ii,iii,ieq
        real(8)    :: dvol

c     Begin program

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)
                xout(iii) = xin(ii,ieq)
              enddo
            enddo
          enddo
        enddo

c     End program

      end subroutine gather

c     scatter
c     ###################################################################
      subroutine scatter(xout,xin)

c     -------------------------------------------------------------------
c     Scatters magnitudes in vector to different arrays indexed by the
c     variable identifier ieq.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8),intent(IN)  :: xin(:)
        real(8),intent(OUT) :: xout(:,:)

c     Local variables

        integer(4) :: i,j,k,ii,iii,ieq

c     Begin program

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              do ieq=1,neqd
                iii = ieq + neqd*(ii-1)
                xout(ii,ieq) = xin(iii)
              enddo

            enddo
          enddo
        enddo

c     End program

      end subroutine scatter

c     XtoCar
c     ###################################################################
      subroutine XtoCar(igrid,x,xout,cov)

c     -------------------------------------------------------------------
c     Gathers magnitudes from different arrays to vector.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: igrid
        real(8),intent(INOUT),dimension(:,:) :: x
        real(8),intent(OUT),optional,dimension(:,:) :: xout
        logical, optional :: cov

c     Local variables

        integer(4) :: i,j,k,ii,nx,ny,nz
        real(8)    :: cnv(3)
        logical    :: covar

c     Begin program

        if (PRESENT(cov)) then
          covar = cov
        else
          covar = .false.
        endif

        if (size(x,2) /= 3) then
          call pstop('XtoCar','Vector has wrong dimensions')
        endif

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (PRESENT(xout)) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                cnv(:) = x(ii,:)

                call transformVectorToCartesian
     .               (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3),covar
     .               ,xout(ii,1),xout(ii,2),xout(ii,3))

              enddo
            enddo
          enddo
        else
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                cnv(:) = x(ii,:)

                call transformVectorToCartesian
     .               (i,j,k,igrid,igrid,igrid
     .               ,cnv(1),cnv(2),cnv(3),covar
     .               ,x(ii,1),x(ii,2),x(ii,3))

              enddo
            enddo
          enddo
        endif

c     End program

      end subroutine XtoCar

c     XtoCurv
c     ###################################################################
      subroutine XtoCurv(igrid,x,xout,cov)

c     -------------------------------------------------------------------
c     Gathers magnitudes from different arrays to vector.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: igrid
        real(8),intent(INOUT),dimension(:,:) :: x
        real(8),intent(OUT),optional,dimension(:,:) :: xout
        logical, optional :: cov

c     Local variables

        integer(4) :: i,j,k,ii,nx,ny,nz
        real(8)    :: car(3)
        logical    :: covar

c     Begin program

        if (PRESENT(cov)) then
          covar = cov
        else
          covar = .false.
        endif

        if (size(x,2) /= 3) then
          call pstop('XtoCar','Vector has wrong dimensions')
        endif

        nx = grid_params%nxv(igrid)
        ny = grid_params%nyv(igrid)
        nz = grid_params%nzv(igrid)

        if (PRESENT(xout)) then
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                call transformVectorToCurvilinear
     .               (i,j,k,igrid,igrid,igrid
     .               ,x(ii,1),x(ii,2),x(ii,3)
     .               ,covar
     .               ,xout(ii,1),xout(ii,2),xout(ii,3))
              enddo
            enddo
          enddo
        else
          do k = 1,nz
            do j = 1,ny
              do i = 1,nx
                ii  = i + nx*(j-1) + nx*ny*(k-1)

                car = x(ii,:)

                call transformVectorToCurvilinear
     .               (i,j,k,igrid,igrid,igrid
     .               ,car(1),car(2),car(3)
     .               ,covar
     .               ,x(ii,1),x(ii,2),x(ii,3))
              enddo
            enddo
          enddo
        endif

c     End program

      end subroutine XtoCurv

c     findCoeffs
c     ###################################################################
      subroutine findCoeffs

c     -------------------------------------------------------------------
c     Finds coefficients for linearized systems in preconditioner
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer(4) :: order,nxx,nyy,nzz,igrid,ii,ivar,igr,i,j
        real(8)    :: dvol,veclap(3)

        real(8), allocatable, dimension(:,:,:,:) :: vector,vel
        real(8), allocatable, dimension(:,:,:,:,:) :: tensor

c     Begin program

        nxx = grid_params%nxv(igx)
        nyy = grid_params%nyv(igy)
        nzz = grid_params%nzv(igz)

        nx  = nxx
        ny  = nyy
        nz  = nzz

        igrid = igx

        order = 2

c     Store density in all grids (w/o BCs)

        grho0%grid(igrid)%array(:,:,:,1) = rho

        call restrictMGArray(IRHO,1,grho0,bcs(:,IRHO),igrid,order)

c     Store magnetic field cnv components in all grids (w/ BCs)

cc        call restrict_B(.true.)

        gb0%grid(igrid)%array(:,:,:,1) = bx
        gb0%grid(igrid)%array(:,:,:,2) = by
        gb0%grid(igrid)%array(:,:,:,3) = bz

        call restrictMGArray(IBX,3,gb0,bcs(:,IBX:IBZ),igrid,order)

c     Store magnetic field cov components in all grids (w/ BCs)

        gb0_cov%grid(igrid)%array(:,:,:,1) = bx_cov
        gb0_cov%grid(igrid)%array(:,:,:,2) = by_cov
        gb0_cov%grid(igrid)%array(:,:,:,3) = bz_cov

        call restrictMGArray(IBX,3,gb0_cov,bcs(:,IBX:IBZ),igrid,order
     .                      ,iscnv=.false.)

c     Store ion velocity components in all grids (w/ BCs)

        gv0%grid(igrid)%array(:,:,:,1) = vx
        gv0%grid(igrid)%array(:,:,:,2) = vy
        gv0%grid(igrid)%array(:,:,:,3) = vz

        call restrictMGArray(IVX,3,gv0,bcs(:,IVX:IVZ),igrid,order)

c     Store current components in all grids (w/o BCs)

cc        call restrictArrayToMGVector(1,nxx,nyy,nzz,jx,mgj0cnv(:,1),igrid
cc     .                  ,order,.false.)
cc        call restrictArrayToMGVector(1,nxx,nyy,nzz,jy,mgj0cnv(:,2),igrid
cc     .                  ,order,.false.)
cc        call restrictArrayToMGVector(1,nxx,nyy,nzz,jz,mgj0cnv(:,3),igrid
cc     .                  ,order,.false.)

c     Find auxiliary quantities and store them in all grids

        !Velocity divergence (w/o BCs)
        allocate(divrgV(0:nxx+1,0:nyy+1,0:nzz+1))

        do k=1,nzz
          do j=1,nyy
            do i=1,nxx
              divrgV(i,j,k) = div(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                           ,vx,vy,vz)
cc              divrgV(i,j,k) = 0d0
            enddo
          enddo
        enddo

        call restrictArrayToMGVector(1,nxx,nyy,nzz,divrgV,mgdivV0,igrid
     .                            ,order,.false.)

        deallocate(divrgV)

        !pressure (w/ BCs)
        gp0%grid(igrid)%array(:,:,:,1) = a_p*rho*tmp

        call restrictMGArray(IRHO,1,gp0,bcs(:,IRHO),igrid,order)

        !v0/dt+theta(v0.grad(v0)-mu*veclap(v0)) and nabla_v0
        allocate(vector(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .          ,vel   (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .          ,tensor(0:nxx+1,0:nyy+1,0:nzz+1,3,3))

        vel (:,:,:,1) = vx
        vel (:,:,:,2) = vy
        vel (:,:,:,3) = vz

        do k = 1,nzz
          do j = 1,nyy
            do i = 1,nxx
              ii  = i + nxx*(j-1) + nxx*nyy*(k-1)

              jac    = gmetric%grid(igrid)%jac(i,j,k)
              nabla_v= fnabla_v(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                         ,vx,vy,vz,0)
              veclap = nuu(i,j,k)
     $                *veclaplacian(i,j,k,nxx,nyy,nzz,igrid,igrid,igrid
     .                             ,vel,alt_eom,vol=.false.)

              do ivar=1,3
                vector(i,j,k,ivar) = 
     .                            gv0%grid(igrid)%array(i,j,k,ivar)/dt
     .                          + alpha*vx(i,j,k)*nabla_v(1,ivar)/jac
     .                          + alpha*vy(i,j,k)*nabla_v(2,ivar)/jac
     .                          + alpha*vz(i,j,k)*nabla_v(3,ivar)/jac
     .                          - alpha*veclap(ivar)
              enddo

              tensor(i,j,k,:,:) = nabla_v
            enddo
          enddo
        enddo

        do ivar=1,3
          call restrictArrayToMGVector(1,nxx,nyy,nzz
     .                              ,vector(:,:,:,ivar)
     .                              ,mgadvdiffV0(:,ivar)
     .                              ,igrid,order,.false.)
        enddo

        do i=1,3
          do j=1,3
            call restrictArrayToMGVector(1,nxx,nyy,nzz
     .                              ,tensor(:,:,:,i,j)
     .                              ,mgnablaV0(:,i,j)
     .                              ,igrid,order,.false.)
          enddo
        enddo

        deallocate(vector,vel,tensor)

c     End program

      contains

c     restrict_B
c     ###############################################################
cc      subroutine restrictCurvVec(ax,ay,az,gvec,order)
      subroutine restrict_B(cartesian)

c     ---------------------------------------------------------------
c     This routines restricts B Cartesian components, and then finds
c     corresponding cov, cnv components with BCs.
c     ---------------------------------------------------------------

      use setMGBC_interface

      use vectorOps

      implicit none

c     Call variables

cc      integer(4) :: order
      logical    :: cartesian

c     Local variables

      integer(4) :: igr,bcmod(6,3),nxc,nyc,nzc,neq

c     Begin program

      neq = 3

      igr = 1

      nxc = grid_params%nxv(igr)
      nyc = grid_params%nyv(igr)
      nzc = grid_params%nzv(igr)

c     Find cartesian coordinates

      gb0%grid(1)%array(:,:,:,1) = bx
      gb0%grid(1)%array(:,:,:,2) = by
      gb0%grid(1)%array(:,:,:,3) = bz

      if (cartesian)
     .   call transformVector(igr,igr,igr,0,nxc+1,0,nyc+1,0,nzc+1
     .                       ,gb0%grid(1)%array(:,:,:,1)
     .                       ,gb0%grid(1)%array(:,:,:,2)
     .                       ,gb0%grid(1)%array(:,:,:,3)
     .                       ,.false.,.true.)

c     Restrict components (w/o BCs)

      call restrictMGArray(0,neq,gb0,bcs(:,IBX:IBZ),igr,order)
cc        call restrictMGArray(IBX,3,gb0,bcs(:,IBX:IBZ),igrid,order
cc     .                      ,iscnv=.true.,isvec=.false.
cc     .                      ,result_is_vec=.true.)

c     Impose BCs, and return cov, cnv components

      bcmod = bcs(:,IBX:IBZ)
      where (bcmod == EQU)
        bcmod = EXT
      end where

      do igr=1,grid_params%ngrid
        nxc = grid_params%nxv(igr)
        nyc = grid_params%nyv(igr)
        nzc = grid_params%nzv(igr)

        call setMGBC(0,neq,nxc,nyc,nzc,igr
     .              ,gb0%grid(igr)%array,bcmod       !Input/output cnv array
     .              ,arr_cov=gb0_cov%grid(igr)%array !Input/output cov array
     .              ,icomp=IBX
     .              ,is_cnv=.true.          !This informs that input is in cnv array
     .              ,is_vec=.false.         !This informs that input is Cartesian
     .              ,result_is_vec=.true.   !This returns curvilinear components
     .              ,iorder=order)
      enddo

c     End program

      end subroutine restrict_B

      end subroutine findCoeffs

c     cSolver
c #   ##################################################################
      subroutine cSolver(neq,ntotp,b,x,bcnd,igrid,out,guess
     .                  ,matvec,dg,ncolors,line_relax,gm_coarse_solve
     .                  ,gm_driver,cvrg_tst,order_res,order_prol,smooth
     .                  ,omega,mvcyc,mu_cycle,galerkin)
c     ---------------------------------------------------------------
c     This subroutine solves a coupled system of neq equations. 
c     In call sequence:
c       * neq: number of coupled equations
c       * ntotp: number of mesh points
c       * b: rhs
c       * x: solution
c       * bcnd: boundary condition defs.
c       * igrid: MG grid level (igrid=1 is finest level)
c       * out: level of output information
c       * guess: whether a non-trivial initial guess is provided
c               (iguess=1) or not (iguess=0)
c       * matvec (external): matrix-vector product definition.
c
c     Optional variables:
c       * dg: matrix neq*neq diagonal block (for stationary its).
c       * ncolors: number of colors in grid (for GS).
c       * line_relax: whether we want line relaxation
c       * gm_smooth: whether we want GMRES as a smoother
c       * gm_driver: whether we want GMRES as a driver
c       * cvrg_tst: whether we want to perform a convergence test
c       * order_res: restriction order
c       * order_prol: restriction order
c       * smooth: determines smoothing (jb,gs,gm; default jb)
c       * omega: damping for JB smoother
c       * mvcyc: maximum number of MG V-cycles (default 1)
c       * mu_cycle: number of mu_cycles in MG (default 1)
c       * galerkin: whether to do Galerkin coarsening of operator
c                  (true) or rediscretization (default: false)
c     ---------------------------------------------------------------

        use mlsolverSetup

        implicit none

c     Call variables

        integer(4) :: neq,ntotp,igrid,bcnd(6,neq),out,guess
        real(8)    :: x(ntotp,neq),b(ntotp,neq)

        real(8)   ,optional,intent(IN) :: dg(neq,2*neq*ntotp),omega
        integer(4),optional,intent(IN) :: ncolors,order_res,order_prol
     .                                   ,mvcyc,mu_cycle
        logical   ,optional,intent(IN) :: line_relax,gm_coarse_solve
     .                                   ,cvrg_tst,gm_driver,galerkin

        character(2),optional :: smooth

        external   matvec

c     Local variables

        integer(4) :: ntot,iorder_res,iorder_prol,vcyc,mu_cyc
        real(8)    :: xi(ntotp*neq),bi(ntotp*neq),omg
        real(8), target :: ddg(neq,2*neq*ntotp)
        logical    :: gm_crse_solve,cvrg_test,lrelax,gm_drvr,glkin

        character(2) :: smth

c     Begin program

c     Process optional arguments

        if (PRESENT(dg)) ddg = dg

        if (PRESENT(gm_coarse_solve)) then
          gm_crse_solve = gm_coarse_solve
        else
          gm_crse_solve = .false.
        endif

        if (PRESENT(gm_driver)) then
          gm_drvr = gm_driver
        else
          gm_drvr = .false.
        endif

        if (PRESENT(cvrg_tst)) then
          cvrg_test = cvrg_tst
        else
          cvrg_test = .false.
        endif

        if (PRESENT(line_relax)) then
          lrelax = line_relax
        else
          lrelax = .false.
        endif

        if(PRESENT(order_res)) then
          iorder_res = order_res
        else
          iorder_res = 2
        endif

        if(PRESENT(order_prol)) then
          iorder_prol = order_prol
        else
          iorder_prol = 2
        endif

        if(PRESENT(mvcyc)) then
          vcyc = mvcyc
        else
          vcyc = maxvcyc
        endif

        if(PRESENT(mu_cycle)) then
          mu_cyc = mu_cycle
        else
cc          if (bcSP()) then
cc            mu_cyc = 2
cc          else
            mu_cyc = 1
cc          endif
        endif

        if(PRESENT(smooth)) then
          smth = smooth
        else
          smth = 'jb'
        endif

        if(PRESENT(omega)) then
          omg = omega
        else
          omg = 0.7
        endif

        if(PRESENT(galerkin)) then
          glkin = galerkin
        else
          glkin = .false.
        endif

c     Gather variables for coupled solve

        if (.not.cvrg_test) then
          do i=1,ntotp
            do ieq=1,neq
              xi(neq*(i-1)+ieq) = x(i,ieq)
              bi(neq*(i-1)+ieq) = b(i,ieq)
            enddo
          enddo
        else
          bi = 0d0
          call random_number(xi)
          guess = 1
        endif

c     Solve coupled MG

c       Initialize solver

        call solverInit

c       Check whether GMRES is driving solve

        if (gm_drvr.and.vcyc > 1) then
          call solverOptionsInit

          solverOptions%iter    = vcyc
          solverOptions%tol     = mgtol

          solverOptions%stp_test = 1 
          solverOptions%krylov_subspace = solverOptions%iter 
          call assembleSolverHierarchy('gm')
        endif

c diag ****
cc        write (*,*) 'DIAG -- cSolver'
cc
cc        call solverOptionsInit
cc
cc        solverOptions%iter    = maxvcyc
cc        solverOptions%tol     = mgtol
cc
cc        if (PRESENT(dg)) solverOptions%diag => ddg
cc
cc        select case (smth)
cc        case('gs')
cc          solverOptions%omega   = 1d0
cc          solverOptions%ncolors = ncolors
cc          call assembleSolverHierarchy('gs')
cc        case('gm')
cc          solverOptions%stp_test = 1 
cc          solverOptions%krylov_subspace = solverOptions%iter 
cc          call assembleSolverHierarchy('gm')
cc
cc          solverOptions%iter    = 1
cc          solverOptions%omega   = omg
cc          call assembleSolverHierarchy('jb')
cccc          call assembleSolverHierarchy('id')
cc        case default !JB
cc          solverOptions%omega   = omg
cc          call assembleSolverHierarchy('jb')
cc        end select
c diag ****

c       MG solver options

        call solverOptionsInit

        solverOptions%tol      = mgtol
        if (gm_drvr .and. vcyc > 1) then
          solverOptions%vcyc     = 1
        else
          solverOptions%vcyc     = vcyc
        endif
cc        if (np > 1) then
cc          solverOptions%igridmin = 1
cc        else
          solverOptions%igridmin = 2
cc        endif
        solverOptions%orderres = iorder_res
        solverOptions%orderprol= iorder_prol
        solverOptions%mg_mu    = mu_cyc
        solverOptions%vol_res  = vol_wgt
        solverOptions%galerkin = glkin

        if (PRESENT(dg)) solverOptions%diag => ddg

        if (PRESENT(ncolors)) solverOptions%ncolors  = ncolors

        !Coarse solve
        if (gm_crse_solve) then
        !With GMRES as coarse solver (defined below)
          if (         (gm_drvr  .and. smth/='gm')
     .       .or.((.not.gm_drvr) .and. smth=='gm') ) then
            solverOptions%mg_coarse_solver_depth = 4
          elseif (gm_drvr .and. smth=='gm') then
            solverOptions%mg_coarse_solver_depth = 5
          endif
        else
        !Use smoother as coarse solver
          solverOptions%mg_coarse_solver_depth = 0
        endif

        !Vertex relaxation
cc        solverOptions%vertex_based_relax = .true.

        !Plane/line relaxation
        if (lrelax) then
          solverOptions%igridmin        = 3
          solverOptions%orderres        = 0
          solverOptions%orderprol       = 2
          solverOptions%mg_mu           = 1    
          solverOptions%mg_line_relax   = .true.
          solverOptions%mg_line_nsweep  = 1
          solverOptions%mg_line_vcyc    = 100
          solverOptions%mg_line_tol     = 1d-1
          solverOptions%mg_line_omega   = 1d0
          solverOptions%mg_line_x       = .true.
          solverOptions%mg_line_y       = .true.
          solverOptions%mg_line_z       = .true.
          solverOptions%mg_line_solve   = "mg"
          solverOptions%mg_line_coarse_solver_depth = 0
        endif

        call assembleSolverHierarchy('mg')

c       MG smoother

        call solverOptionsInit

        solverOptions%iter    = nsweep
        solverOptions%tol     = mgtol

        if (lrelax) then
          solverOptions%iter  = 2
        endif

        select case (smth)
        case('gs')
          solverOptions%omega   = 1d0
          solverOptions%ncolors = ncolors
          call assembleSolverHierarchy('gs')
        case('gm')
          solverOptions%stp_test        = 1 
          solverOptions%krylov_subspace = solverOptions%iter 
          call assembleSolverHierarchy('gm')

          if (gm_smooth) then
            call assembleSolverHierarchy('id')
          else
            call solverOptionsInit
            solverOptions%iter    = 1
            solverOptions%omega   = omg
            call assembleSolverHierarchy('jb')
          endif
        case default !JB
          solverOptions%omega   = omg
          call assembleSolverHierarchy('jb')
        end select

c       MG coarsest grid solve

        call solverOptionsInit

        solverOptions%tol             = 1d-5
        solverOptions%krylov_subspace = 100
        solverOptions%iter            = 100
        solverOptions%stp_test        = 1 
        solverOptions%omega           = 1d0

        call assembleSolverHierarchy('gm')
        call assembleSolverHierarchy('id') !GMRES preconditioner

c       Coarsest grid solve for inner line/plane MG

        call solverOptionsInit

        solverOptions%omega   = 1d0
        solverOptions%iter    = 100
        solverOptions%tol     = 1d-4

        call assembleSolverHierarchy('gs')

c       Invoke solver

        ntot=neq*ntotp
        call getSolver(neq,ntot,bi,xi,matvec,igrid,bcnd,guess,out,1)

c       Get output data

cc      call getSolverOptions(1)

c       Kill solver

        call solverKill

c     Scatter solution for output

        do i = 1,ntotp
          do ieq=1,neq
            x(i,ieq) = xi(neq*(i-1)+ieq)
          enddo
        enddo

c     End program

      end subroutine cSolver

c     findVectorPotential
c     ###############################################################
      subroutine findVectorPotential(nx,ny,nz,igx,igy,igz,b,a)

c     ---------------------------------------------------------------
c     Finds COVARIANT vector potential from CONTRAVARIANT magnetic
c     field components. Employs gauge A1 = 0d0. In SP systems,
c     it integrates at faces in r, and then averages to nodes.
c     
c     Integrals done by Numerical Recipes p. 128, 4.1.12 -- actually
c     one bit better -- extrapolates quadratically to f(-1) and uses
c     the internal formula, rather than using midpoint for the first
c     point. Thus even the first point has third order accuracy.
c     ---------------------------------------------------------------

      use setMGBC_interface

      implicit none

c     Call variables

      integer(4) :: nx,ny,nz,igx,igy,igz
      real(8)    :: a(0:nx+1,0:ny+1,0:nz+1,3)
     .             ,b(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      integer(4) :: i,j,k,igrid,ii,bcnd(6,3),ig,jg,kg,ivar
      real(8)    :: cm,c0,cp,dxx,dyy,dvol,cov(3)

c     Begin program

      igrid = igx

c     New code (SP ready, version I)

      a(:,:,:,1)=0d0  !Gauge

      a(0,0,:,3)=0d0

      i = 0

      do k=0,nz+1
        do j=1,ny+1
          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

          dyy = grid_params%dy(jg-1)

          if (bcSP()) then  !r=0 face
            a(i,j,k,3)=a(i,j-1,k,3) + dyy*0.25*(b(i  ,j-1,k,1)
     .                                         +b(i  ,j  ,k,1)
     .                                         +b(i+1,j-1,k,1)
     .                                         +b(i+1,j  ,k,1))
          else
            a(i,j,k,3)=a(i,j-1,k,3) + dyy*0.5*(b(i,j-1,k,1)
     .                                        +b(i,j  ,k,1))
          endif

        enddo
      enddo

      a(0,:,:,2)=0d0

      do k=0,nz+1
        do j=0,ny+1
          do i=1,nx+1

            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            if (bcSP()) then

              dxx = grid_params%dxh(ig-1)

              a(i,j,k,3)=a(i-1,j,k,3) -dxx*b(i,j,k,2)

              if (isSP(i,j,k,igx,igy,igz)) then
                a(i,j,k,2)=a(i-1,j,k,2) +0.5*dxx*b(i,j,k,3)  !Factor of 1/2 due to geometry
              else
                a(i,j,k,2)=a(i-1,j,k,2) +    dxx*b(i,j,k,3)
              endif

            else

              dxx = grid_params%dx(ig-1)

              a(i,j,k,3)=a(i-1,j,k,3) -dxx*0.5*(b(i-1,j,k,2)
     .                                         +b(i  ,j,k,2))
              a(i,j,k,2)=a(i-1,j,k,2) +dxx*0.5*(b(i-1,j,k,3)
     .                                         +b(i  ,j,k,3))

            endif
          enddo
        enddo
      enddo

      !Average from radial faces to nodes in SP coordinate systems
      if (bcSP())
     .   a(1:nx+1,:,:,2:3) = 0.5*(a(1:nx+1,:,:,2:3)+a(0:nx,:,:,2:3))

      !Impose BCs
      bcnd(:,1) = bcond
      bcnd(:,2) = bcond
      bcnd(:,3) = bcond
cc      where (bcnd == DEF) bcnd = EXT

      call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IBX
     .            ,is_vec=.true.,is_cnv=.false.,iorder=2)

      return

c     New code (SP-ready, version II)

      a(:,:,:,1)=0d0  !Gauge

      if (bcSP()) then

        a(0,0,:,3)=0d0

        i = 0

        do k=0,nz+1
          do j=1,ny+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            dyy = grid_params%dy(jg-1)

            a(i,j,k,3)=a(i,j-1,k,3) + dyy*0.25*(b(i  ,j-1,k,1)
     .                                         +b(i  ,j  ,k,1)
     .                                         +b(i+1,j-1,k,1)
     .                                         +b(i+1,j  ,k,1))
          enddo
        enddo

        a(0,:,:,2)=0d0

        do k=0,nz+1
          do j=0,ny+1
            do i=1,nx+1

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              dxx = grid_params%dxh(ig-1)

              a(i,j,k,3)=a(i-1,j,k,3) -dxx*b(i,j,k,2)

              if (isSP(i,j,k,igx,igy,igz)) then
                a(i,j,k,2)=a(i-1,j,k,2) +0.5*dxx*b(i,j,k,3) !Factor of 1/2 due to geometry
              else
                a(i,j,k,2)=a(i-1,j,k,2) +    dxx*b(i,j,k,3)
              endif

            enddo
          enddo
        enddo

        !Average from radial faces to centers
        a(1:nx+1,:,:,2:3) = 0.5*(a(1:nx+1,:,:,2:3)+a(0:nx,:,:,2:3))

        !Boundary conditions: topological constraints (SP, PER)

        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond

        call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IBX
     .            ,is_vec=.true.,is_cnv=.false.,iorder=2)

      else

        a(0,0,:,3)=0d0

        i = 0

        do k=0,nz+1
          do j=1,ny+1
            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

            dyy = grid_params%dy(jg-1)

            a(i,j,k,3)=a(i,j-1,k,3) + dyy*0.5*(b(i,j-1,k,1)
     .                                        +b(i,j  ,k,1))

          enddo
        enddo

        a(0,:,:,2)=0d0

        do k=0,nz+1
          do j=0,ny+1
            do i=1,nx+1

              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              dxx = grid_params%dx(ig-1)

              a(i,j,k,3)=a(i-1,j,k,3) -dxx*0.5*(b(i-1,j,k,2)
     .                                         +b(i  ,j,k,2))
              a(i,j,k,2)=a(i-1,j,k,2) +dxx*0.5*(b(i-1,j,k,3)
     .                                         +b(i  ,j,k,3))

            enddo
          enddo
        enddo

        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond
        
        call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IBX
     .              ,is_vec=.true.,is_cnv=.false.,iorder=2)


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

      endif


c     End program

      end subroutine findVectorPotential

c     vecA
c     ###################################################################
      function vecA(i,j,k,nx,ny,nz,igx,igy,igz,vv,vv_cov,bb,half_elem)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of vector
c     
c       A = curl(bb x vv + di Dv vv)
c
c     at the grid node (i,j,k). One sided derivatives are employed when
c     half_elem=1 (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer(4) :: i,j,k,nx,ny,nz,half_elem,igx,igy,igz
        real(8)    :: cnv(3)
        real(8)    :: vv    (0:nx+1,0:ny+1,0:nz+1,3)
     .               ,vv_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .               ,bb    (0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

        integer(4) :: ig,jg,kg
        real(8)    :: vxx,vyy,vzz,minv,idx,idy,idz

c     Begin program

c     Resistive MHD part

        cnv = curl_bxv(i,j,k,nx,ny,nz,igx,igy,igz,vv,bb,half_elem)

c     Hall MHD part

        if (di > 0d0) then

          !Hall term coefficient
cc          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc          idx  = 1./grid_params%dxh(ig)
cc          if (nxd == 1) idx = 0d0
cc          idy  = 1./grid_params%dyh(jg)
cc          if (nyd == 1) idy = 0d0
cc          idz  = 1./grid_params%dzh(kg)
cc          if (nzd == 1) idz = 0d0
cc
cc          vxx = abs(vv(i,j,k,1))
cc          vyy = abs(vv(i,j,k,2))
cc          vzz = abs(vv(i,j,k,3))
cc
cc          minv = cnp/dt
cc     .           +alpha*(scalarProduct(i,j,k,igx,igy,igz
cc     .                              ,idx,idy,idz,vxx,vyy,vzz)
cc     .           +nu*vectorNorm(i,j,k,igx,igy,igz,idx,idy,idz,.true.))

cc          minv = minv*di/alpha
          minv = di*cnp/dt/alpha

          !Hall term correction
          cnv = cnv
     .        + minv*curl(i,j,k,nx,ny,nz,igx,igy,igz
     .                   ,vv_cov(:,:,:,1)
     .                   ,vv_cov(:,:,:,2)
     .                   ,vv_cov(:,:,:,3),he=half_elem)
     .         
        endif

c     End program

      end function vecA

c     diag_v
c     #####################################################################
      function diag_v(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,da) result(psiv)

c     ---------------------------------------------------------------------
c     Find advection/diffusion component of SI operator
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer(4) :: i,j,k,igx,igy,igz,nxx,nyy,nzz
      real(8)    :: da(0:nxx+1,0:nyy+1,0:nzz+1,3)
     .             ,v0(0:nxx+1,0:nyy+1,0:nzz+1,3),psiv(3)

c     Local Variables

      integer(4) :: hex,hey,hez
      real(8)    :: upwind,mul,nabla_vv0(3,3),veclap(3)

c     Begin program

      hex = 1
      hey = 1
      hez = 1
      if (v0(i,j,k,1) > 0d0) hex = -1
      if (v0(i,j,k,2) > 0d0) hey = -1
      if (v0(i,j,k,3) > 0d0) hez = -1
cc      hex = 0
cc      hey = 0
cc      hez = 0

      nabla_v = fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                       ,da(:,:,:,1),da(:,:,:,2),da(:,:,:,3)
     .                       ,hex,hey,hez)

      nabla_vv0 =fnabla_v_upwd(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                        ,v0(:,:,:,1)
     .                        ,v0(:,:,:,2)
     .                        ,v0(:,:,:,3)
     .                        ,0,0,0)

      mul = vis(i,j,k,nxx,nyy,nzz,igx,igy,igz)

      veclap = veclaplacian(i,j,k,nxx,nyy,nzz,igx,igy,igz,da
     .                     ,alt_eom,vol=.false.)

      do ieq=1,3
        upwind =( v0(i,j,k,1)*nabla_v(1,ieq)
     .           +v0(i,j,k,2)*nabla_v(2,ieq)
     .           +v0(i,j,k,3)*nabla_v(3,ieq))/jac

        upwind = upwind
     .           +( da(i,j,k,1)*nabla_vv0(1,ieq)
     .             +da(i,j,k,2)*nabla_vv0(2,ieq)
     .             +da(i,j,k,3)*nabla_vv0(3,ieq))/jac

cc        if (.not.nc_eom_v) then
cc          hex = floor(sign(1d0,-mgadvdiffV0(ijkg,ieq)))
cccc          hex = 0
cc          upwind = upwind
cc     .          -dt*mgadvdiffV0(ijkg,ieq)*div_upwd(da,hex)/rho0(i,j,k,1)
cc        endif

        psiv(ieq) = cnp*da(i,j,k,ieq)/dt
     .            + alpha*upwind
     .            - alpha*mul*veclap(ieq)
      enddo

      end function diag_v

c     diagonalScaling
c     #####################################################################
      subroutine diagonalScaling(neq,ntotp,idiag,y,x,igrid)
c     ---------------------------------------------------------------
c     Performs x=idiag*y, where idiag contains the inverse of the
c     diagonal.
c     ---------------------------------------------------------------

cc      use mg_internal

      implicit none

c     Call variables

      integer(4) :: neq,ntotp,igrid
      real(8)    :: idiag(neq,*),y(ntotp,neq),x(ntotp,neq)

c     Local variables

      integer(4) :: ii,iii,iig,isig
      logical    :: fpointers

c     Begin program

      call allocPointers(neq,fpointers)

      isig = istart(igrid)

      do ii = 1,ntotp

        iii = neq*(ii-1)
        iig = iii + isig - 1

        x(ii,:) = matmul(transpose(idiag(:,iig+1:iig+neq)),y(ii,:))

      enddo

      call deallocPointers(fpointers)

c     End program

      end subroutine diagonalScaling

      end module precond_variables
