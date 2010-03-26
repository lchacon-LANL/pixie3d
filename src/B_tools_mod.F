c module B_tools
c ####################################################################
      module B_tools

        use setMGBC_interface

        use precond_variables

        use local_BCS_variables, ONLY: bc_order

        use xdraw_io

      contains

c     curl_inv
c     ###############################################################
      subroutine curl_inv(nx,ny,nz,igx,igy,igz,b,a,global,w_bcs)

c     ---------------------------------------------------------------
c     Finds COVARIANT vector potential "a" from CONTRAVARIANT vector
c     field components "b". Employs gauge A_1 = 0d0. In SP systems,
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
      logical,optional :: global,w_bcs

c     Local variables

      integer    :: i,j,k,igrid,ii,bcnd(6,3),ig,jg,kg,ivar
      real(8)    :: cm,c0,cp,dxx,dyy,dvol,cov(3)

      logical    :: spoint,glbl,do_bc

c     Begin program

      if (PRESENT(global)) then
        glbl = global
      else
        glbl = .false.
      endif

      if (PRESENT(w_bcs)) then
        do_bc = w_bcs
      else
        do_bc = (.not.glbl).or.(np == 1)
      endif

      igrid = igx

      spoint = bcSP()

c     Consistency check

      if (igrid /= 1 .and. glbl) then
        call pstop('curl_inv','Inconsistent grid level')
      endif

c     Gauge

      a(:,:,:,1)=0d0

c     Quadrant-based code (SP ready, parallel-ready)

      if (spoint) then  !Separate quadrants only in theta

cc        call line_int(1,nx,1,ny,b,a)

        !Lower-right quadrant
        call line_int(1,nx,ny/2-1,1
     .               ,b(0:nx+1,0:ny/2,0:nz+1,1:3)
     .               ,a(0:nx+1,0:ny/2,0:nz+1,1:3))

        !Upper-right quadrant
        call line_int(1,nx,ny/2+1,ny
     .               ,b(0:nx+1,ny/2:ny+1,0:nz+1,1:3)
     .               ,a(0:nx+1,ny/2:ny+1,0:nz+1,1:3))
      else
        !Lower-left quadrant
        call line_int(nx/2-1,1,ny/2-1,1
     .               ,b(0:nx/2,0:ny/2,0:nz+1,1:3)
     .               ,a(0:nx/2,0:ny/2,0:nz+1,1:3))

        !Lower-right quadrant
        call line_int(nx/2+1,nx,ny/2-1,1
     .               ,b(nx/2:nx+1,0:ny/2,0:nz+1,1:3)
     .               ,a(nx/2:nx+1,0:ny/2,0:nz+1,1:3))

        !Upper-left quadrant
        call line_int(nx/2-1,1,ny/2+1,ny
     .               ,b(0:nx/2,ny/2:ny+1,0:nz+1,1:3)
     .               ,a(0:nx/2,ny/2:ny+1,0:nz+1,1:3))

        !Upper-right quadrant
        call line_int(nx/2+1,nx,ny/2+1,ny
     .               ,b(nx/2:nx+1,ny/2:ny+1,0:nz+1,1:3)
     .               ,a(nx/2:nx+1,ny/2:ny+1,0:nz+1,1:3))
      endif

      !Boundary conditions: topological constraints (SP, PER)
      if (do_bc) then
        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond

        call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IAX
     .              ,is_vec=.true.,is_cnv=.false.,iorder=2)
      endif

      return

c     New code (SP ready, parallel-ready, version I)

      call line_int(1,nx,1,ny,b,a)

      !Boundary conditions: topological constraints (SP, PER)
      if (do_bc) then
        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond

        call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IAX
     .              ,is_vec=.true.,is_cnv=.false.,iorder=2)
      endif

cc      a(0,0,:,3)=0d0
cc
cc      i = 0
cc
cccc      spoint = isSP(1,1,1,igx,igy,igz)  !Check whether current domain contains SP
cc      spoint = .false.
cc
cc      do j=1,ny+1
cc        if (.not.glbl) then
cc          call getMGmap(i,j,1,igx,igy,igz,ig,jg,kg)
cc          dyy = grid_params%dy(jg-1)
cc        else
cc          dyy = grid_params%yg(j)-grid_params%yg(j-1)
cc        endif
cc
cc        if (spoint) then  !r=0 face
cc          a(i,j,:,3)=a(i,j-1,:,3) + dyy*2.5d-1*(b(i  ,j-1,:,1)
cc     .                                         +b(i  ,j  ,:,1)
cc     .                                         +b(i+1,j-1,:,1)
cc     .                                         +b(i+1,j  ,:,1))
cc        else
cc          a(i,j,:,3)=a(i,j-1,:,3) + dyy*  5d-1*(b(i,j-1,:,1)
cc     .                                         +b(i,j  ,:,1))
cc        endif
cc      enddo
cc
cc      a(0,:,:,2)=0d0
cc
cc      do i=1,nx+1
cc        if (.not.glbl) call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cc
cc        if (spoint) then
cc          if (.not.glbl) then
cc            dxx = grid_params%dxh(ig-1)
cc          else
cccc            dxx = 0.5*(grid_params%xg(i+1)-grid_params%xg(i-1))
cc            dxx = grid_params%xg(i)-grid_params%xg(i-1)
cc          endif
cc
cc          a(i,:,:,3) = a(i-1,:,:,3) - dxx*b(i,:,:,2)
cc
cc          if (isSP(i,1,1,igx,igy,igz)) then
cc            a(i,:,:,2) = a(i-1,:,:,2) + 5d-1*dxx*b(i,:,:,3)  !Factor of 1/2 due to geometry
cc          else
cc            a(i,:,:,2) = a(i-1,:,:,2) +      dxx*b(i,:,:,3)
cc          endif
cc
cc        else
cc
cc          if (.not.glbl) then
cc            dxx = grid_params%dx(ig-1)
cc          else
cc            dxx = grid_params%xg(i)-grid_params%xg(i-1)
cc          endif
cc
cc          a(i,:,:,3) = a(i-1,:,:,3) - dxx*5d-1*(b(i-1,:,:,2)
cc     .                                         +b(i  ,:,:,2))
cc          a(i,:,:,2) = a(i-1,:,:,2) + dxx*5d-1*(b(i-1,:,:,3)
cc     .                                         +b(i  ,:,:,3))
cc
cc        endif
cc      enddo
cc
cc      !Average from radial faces to nodes in SP coordinate systems
cc      if (spoint) a(1:nx+1,:,:,2:3) = 0.5*(a(1:nx+1,:,:,2:3)
cc     .                                    +a(0:nx  ,:,:,2:3))
cc

cc#if !defined(petsc)
cc      if (spoint) then
cc        !Impose BCs
cc        bcnd(:,1) = bcond
cc        bcnd(:,2) = bcond
cc        bcnd(:,3) = bcond
cc        where (bcnd /= SP) bcnd = DEF
cccc        where (bcnd == DEF) bcnd = EXT
cc
cc        call setMGBC(0,3,nx,ny,nz,igrid,da_cov,bcnd,icomp=IBX
cc     .              ,is_vec=.true.,is_cnv=.false.,iorder=2)
cccc        call setBC(IBX,3,nx,ny,nz,da_cov,da_cnv,vzeros,bcnd,igrid,igrid,igrid
cccc     .            ,is_cnv=.false.,iorder=2)
cc
cc      endif
cc
cc      !Boundary conditions: topological constraints (SP, PER)
cc      bcnd(:,1) = bcond
cc      bcnd(:,2) = bcond
cc      bcnd(:,3) = bcond
cc
cc      call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IAX
cc     .          ,is_vec=.true.,is_cnv=.false.,iorder=2)
cc
cc#endif

      return

ccc     New code (SP-ready, NOT parallel-ready, version II)
cc
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
cc
ccc     JOHN FINN's code (assumes uniform mesh)
cc
cc        cm=-0.0833333333333
cc        c0= 0.6666666666667
cc        cp= 0.4166666666667
cc
cc        a(1,1,:,3) = 0d0
cc
cc        do k=1,nz
cc          do j=2,ny
cc            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
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
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
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

      contains

c     line_int
c     #############################################################
      subroutine line_int(ilo,ihi,jlo,jhi,bb,aa)

c     -------------------------------------------------------------
c     Performs line integral to find vector potential from magnetic
c     field from starting point (ilo,jlo) in any given quadrant
c     defined by [ilo,ihi]x[jlo,jhi]. We do NOT assume ilo<ihi or
c     jlo<jhi.
c     -------------------------------------------------------------

c     Call variables

        integer :: ilo,ihi,jlo,jhi
        real(8) :: bb(min(ilo,ihi)-1:max(ilo,ihi)+1
     .               ,min(jlo,jhi)-1:max(jlo,jhi)+1,0:nz+1,1:3)
        real(8) :: aa(min(ilo,ihi)-1:max(ilo,ihi)+1
     .               ,min(jlo,jhi)-1:max(jlo,jhi)+1,0:nz+1,1:3)

c     Local variables

        integer :: i,im,ip,imin,imax,istep
     .            ,j,jm,jp,jmin,jmax,jstep

c     Begin program

        istep = -1
        if (ilo < ihi) istep = 1

        jstep = -1
        if (jlo < jhi) jstep = 1

        a(ilo-istep,jlo-jstep,:,3)=0d0

        i = ilo-istep

        jmin = jlo-min(jstep,0)
        jmax = jhi+max(jstep,0)

        do j=jmin,jmax,jstep
          jm = j-max(jstep,0)
          jp = j+min(jstep,0)

          if (.not.glbl) then
            call getMGmap(i,jm,1,igx,igy,igz,ig,jg,kg)
            dyy = sign(1,jstep)*grid_params%dy(jg)
          else
            dyy = grid_params%yg(jp)
     .           -grid_params%yg(jm)
          endif

cc          if (spoint) then  !r=0 face
cc            aa(i,j,:,3)=aa(i,j-1,:,3) + dyy*2.5d-1*(bb(i  ,j-1,:,1)
cc     .                                             +bb(i  ,j  ,:,1)
cc     .                                             +bb(i+1,j-1,:,1)
cc     .                                             +bb(i+1,j  ,:,1))
cc          else
            aa(i,jp,:,3)=aa(i,jm,:,3) + dyy*5d-1*(bb(i,jm,:,1)
     .                                           +bb(i,jp,:,1))
cc          endif
        enddo

        aa(ilo-istep,:,:,2)=0d0

        imin = ilo-min(istep,0)
        imax = ihi+max(istep,0)

        do i=imin,imax,istep

          im = i-max(istep,0)
          ip = i+min(istep,0)

cc          if (.not.glbl) call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cc
cc          if (spoint) then
cc            if (.not.glbl) then
cc              dxx = grid_params%dxh(ig-1)
cc            else
cccc              dxx = 0.5*(grid_params%xg(i+1)-grid_params%xg(i-1))
cc              dxx = grid_params%xg(i)-grid_params%xg(i-1)
cc            endif
cc
cc            a(i,:,:,3) = a(i-1,:,:,3) - dxx*b(i,:,:,2)
cc
cc            if (isSP(i,1,1,igx,igy,igz)) then
cc              a(i,:,:,2) = a(i-1,:,:,2) + 5d-1*dxx*b(i,:,:,3)  !Factor of 1/2 due to geometry
cc            else
cc              a(i,:,:,2) = a(i-1,:,:,2) +      dxx*b(i,:,:,3)
cc            endif
cc
cc          else

            if (.not.glbl) then
              call getMGmap(im,1,1,igx,igy,igz,ig,jg,kg)
              dxx = sign(1,istep)*grid_params%dx(ig)
            else
              dxx = grid_params%xg(ip)
     .             -grid_params%xg(im)
            endif

            aa(ip,:,:,3) = aa(im,:,:,3) - dxx*5d-1*(bb(im,:,:,2)
     .                                             +bb(ip,:,:,2))
            aa(ip,:,:,2) = aa(im,:,:,2) + dxx*5d-1*(bb(im,:,:,3)
     .                                             +bb(ip,:,:,3))

cc          endif
        enddo

      !Average from radial faces to nodes in SP coordinate systems
cc        if (spoint) aa(1:nx+1,:,:,2:3) = 0.5*(aa(1:nx+1,:,:,2:3)
cc     .                                       +aa(0:nx  ,:,:,2:3))

      end subroutine line_int

      end subroutine curl_inv

c     vlap_div_clean
c     #################################################################
      subroutine vlap_div_clean(iout,nx,ny,nz,igx,igy,igz,db,dbcov,dacov
     .                         ,w_bcs,diag_plots)

c     -----------------------------------------------------------------
c     Finds curl(dA)=B^(n+1)-B^n, with dA=inv_lap(dJ). 
c     -----------------------------------------------------------------

cc      use matvec, rhs_bc => rhs, bcs_pc => bcs

      implicit none

c     Call variables

      integer :: nx,ny,nz,igx,igy,igz,iout
      real(8) :: db(0:nx+1,0:ny+1,0:nz+1,3)

      real(8),optional,INTENT(IN) :: dbcov(0:nx+1,0:ny+1,0:nz+1,3)
      real(8),optional,INTENT(OUT):: dacov(0:nx+1,0:ny+1,0:nz+1,3)

      logical,optional :: w_bcs,diag_plots

c     Local variables

      integer :: i,j,k,ii,icomp,bcsa(6,3),bcsb(6,3)

      real(8) :: da     (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,da_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,db_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,rhs(ntotdp,3),x(ntotdp,3),dvol

      logical :: vwgt,wbcs,dplots

      external lap_mtvc,avlap_mtvc

c     Begin program

      vwgt = .true.

      if (PRESENT(w_bcs)) then
        wbcs = w_bcs
      else
        wbcs = .false.
      endif

      if (PRESENT(diag_plots)) then
        dplots = diag_plots
      else
        dplots = .false.
      endif

c     Find rhs of vec_lap solve: rhs=vol*dJ (with dJ=curl(dB); dB already has BCs)

      bcsb(:,1) = u_0%array_var(IBX)%bconds
      bcsb(:,2) = u_0%array_var(IBY)%bconds
      bcsb(:,3) = u_0%array_var(IBZ)%bconds

      if (PRESENT(dbcov)) then
        db_cov = dbcov
      else
        db_cov = db
        call transformVector(igx,igy,igz,0,nx+1,0,ny+1,0,nz+1
     .                    ,db_cov(:,:,:,1)
     .                    ,db_cov(:,:,:,2)
     .                    ,db_cov(:,:,:,3)
     .                    ,'cnv','cov')
      endif

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
      where (bcsa == -NEU) bcsa = EQU  !Use covariant components for tangential BCs

c     Allocate diagonal matrix

      if (.not.allocated(vlap_diag)) then
        allocate(vlap_diag(3,3*ntotd2p))

        call find_mf_diag(3,1,3*ntotdp,avlap_mtvc,igx,bcsa
     .                   ,vlap_diag,sm_ncolors,vwgt,.false.)
      endif

c     Invert vector laplacian to find dA=-vlap^-1(dJ)

      x = 0d0

      !Solve Poisson equation
      if (coords == 'car') then
        call cSolver(3,ntotdp,rhs,x,bcsa,igx
     .            ,iout,1,avlap_mtvc,vwgt
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

        call fillArray(igx,IAX,3,x,da,bcsa,0
     .                ,is_vec=.true.,is_cnv=.true.,arr_cov=da_cov)
      else
        call cSolver(3,ntotdp,rhs,x,bcsa,igx
     .            ,iout,1,avlap_mtvc,vwgt
     .            ,dg            = vlap_diag
     .            ,gm_sing_matrix= .false.
     .            ,ks_it         = 20
     .            ,mg_order_res  = 0  !To ensure conservation (avoid null space)
     .            ,mg_order_prol = 2
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
cc     .            ,tol           = 1d-4
cc     .            ,sm_it         = 3
cc     .            ,sm_zebra_relax= .false.
     .            )

        where (bcsa == EQU) bcsa =-EQU !Use covariant components for tangential BCs

        call fillArray(igx,IAX,3,x,da,bcsa,bc_order
     .                ,is_vec=.true.,is_cnv=.true.,arr_cov=da_cov)
      endif

c     Find new dB (with BCs)

      db = curl(nx,ny,nz,igx,igy,igz,da_cov)

      if (wbcs) then
        call setMGBC(0,3,nx,ny,nz,igx,db,bcsb
     .            ,icomp=IBX,is_vec=.true.,is_cnv=.true.
     .            ,iorder=bc_order)
      endif

c     Diagnostics

      if (dplots) then
        write (*,*) 'DIAG -- vlap_div_clean'

        open(unit=110,file='divclean.bin',form='unformatted'
     .      ,status='replace')
        call contour(da_cov(:,:,1,1),nx+2,ny+2
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,0,110)
        call contour(da_cov(:,:,1,2),nx+2,ny+2
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,1,110)
        call contour(da_cov(:,:,1,3),nx+2,ny+2
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,1,110)
        call contour(db    (:,:,1,1),nx+2,ny+2
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,1,110)
        call contour(db    (:,:,1,2),nx+2,ny+2
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,1,110)
        call contour(db    (:,:,1,3),nx+2,ny+2
     .                ,grid_params%gxmin,grid_params%gxmax
     .                ,grid_params%gymin,grid_params%gymax
     .                ,1,110)
        close (110)
      endif

c     Store helical flux (for resistivity)

      if (PRESENT(dacov)) dacov = da_cov

      end subroutine vlap_div_clean

      end module B_tools
