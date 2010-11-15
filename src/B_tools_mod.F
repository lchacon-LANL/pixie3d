c module B_tools
c ####################################################################
      module B_tools

        use setMGBC_interface

        use precond_variables

        use local_BCS_variables, ONLY: bc_order

        use app_iosetup

        use xdraw_io

        integer :: bc_order_a

        logical :: is_a_cnv=.true.

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

        if (ny == 1) then
          call line_int(1,nx,1,ny,b,a)
        else
          !Lower-right quadrant
          call line_int(1,nx,ny/2-1,1
     .               ,b(0:nx+1,0:ny/2,0:nz+1,1:3)
     .               ,a(0:nx+1,0:ny/2,0:nz+1,1:3))

          !Upper-right quadrant
          call line_int(1,nx,ny/2+1,ny
     .               ,b(0:nx+1,ny/2:ny+1,0:nz+1,1:3)
     .               ,a(0:nx+1,ny/2:ny+1,0:nz+1,1:3))
        endif
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

        call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=(/IAX/)
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

        call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=(/IAX/)
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
cc      call setMGBC(0,3,nx,ny,nz,igrid,a,bcnd,icomp=IBX
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

c     curl_div_clean
c     #################################################################
      subroutine curl_div_clean(nx,ny,nz,igx,igy,igz,db)

c     -----------------------------------------------------------------
c     Finds B^(n+1)=B^n + curl(dA), with dA=inv_curl(dB)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,igx,igy,igz
      real(8) :: db(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      real(8) :: da_cov(0:nx+1,0:ny+1,0:nz+1,3)

c     Find dA from dB (w/ ghost cells)

      call curl_inv(nx,ny,nz,igx,igy,igz,db,da_cov)

      db = curl(igx,da_cov)

      end subroutine curl_div_clean

c     vlap_div_clean
c     #################################################################
      subroutine vlap_div_clean(iout,nx,ny,nz,igx,igy,igz,db,dbcov,dacov
     .                         ,w_bcs,diag_plots,pc_call)

c     -----------------------------------------------------------------
c     Finds curl(dA)=B^(n+1)-B^n, with dA=inv_lap(dJ). 
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,igx,igy,igz,iout
      real(8) :: db(0:nx+1,0:ny+1,0:nz+1,3)

      real(8),optional,INTENT(IN) :: dbcov(0:nx+1,0:ny+1,0:nz+1,3)
      real(8),optional,INTENT(OUT):: dacov(0:nx+1,0:ny+1,0:nz+1,3)

      logical,optional :: w_bcs,diag_plots,pc_call

c     Local variables

      integer :: i,j,k,ii,icomp,bcsa(6,3),bcsb(6,3)
     .          ,order_res,order_prol

      real(8) :: da_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,da_cnv (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,db_cov (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,rhs(ntotdp,3),x(ntotdp,3),dvol,etol

      logical :: vwgt,wbcs,dplots,pc,sing_mat

      character(20) :: dcf

c     Begin program

      vwgt = .true.

      if (PRESENT(diag_plots)) then
        dplots = diag_plots
      else
        dplots = .false.
      endif

      if (PRESENT(w_bcs)) then
        wbcs = w_bcs
      else
        wbcs = dplots
      endif

      if (PRESENT(pc_call)) then
        pc = pc_call
      else
        pc = .false.
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
            if (.not.is_a_cnv)
     .        rhs(ii,:)=matmul(gmetric%grid(igx)%gsub(i,j,k,:,:)
     .                        ,rhs(ii,:))
          enddo
        enddo
      enddo

c     Define BCs for dA

      bcsa = bcsb
      where (bcond == FSYM)
        bcsa(:,1) = FSYM !Recover symmetry BCs
        bcsa(:,2) = FSYM !Recover symmetry BCs
        bcsa(:,3) = FSYM !Recover symmetry BCs
      end where
      where (bcsa == -NEU) bcsa = -EQU  !Use covariant components for tangential BCs

      call vector_symm_bc(bcsa,-1)

c     Setup solver parameters

      if (coords == 'car') then
        order_res  = 0  !To ensure conservation (avoid null space)
        order_prol = 1

        bc_order_a = 0

        sing_mat   = .true.
      else
        order_res  = 0
        order_prol = 2

        bc_order_a = 1

        sing_mat   = .false.
      endif

c     Allocate diagonal matrix

      if (.not.allocated(vlap_diag)) then
        allocate(vlap_diag(3,3*ntotd2p))

        call find_mf_diag(3,1,3*ntotdp,avlap_mtvc,igx,bcsa
     .                   ,vlap_diag,sm_ncolors,vwgt,.false.)
      endif

c     Invert vector laplacian to find dA=-vlap^-1(dJ)

      if (pc) then
        etol = 0.5*pc_tol
      else
        etol = 0.5*rtol
      endif

      x = 0d0

      call cSolver(3,ntotdp,rhs,x,bcsa,igx
     .            ,iout,1,avlap_mtvc,vwgt
     .            ,tol           = etol
     .            ,iters         = it_divcl
     .            ,dg            = vlap_diag
     .            ,gm_driver     = .not.pc
     .            ,gm_sing_matrix= sing_mat
     .            ,ks_it         = 30
     .            ,mg_vcyc       = mg_vcyc
     .            ,mg_order_res  = order_res
     .            ,mg_order_prol = order_prol
     .            ,mg_coarse_grid_size = mg_coarse_size
     .            ,mg_gm_coarse_solve  = mg_gm_coarse
cc     .            ,mg_gm_coarse_solve  = .true.
     .            ,sm_ncolors    = sm_ncolors
     .            ,sm_line_relax = mg_line_relax
     .            !Zebra relax
cc     .            ,mg_smooth     = 'gs'
cc     .            ,sm_omega      = 0.9d0
cc     .            ,sm_it         = 2
cc     .            ,sm_zebra_relax= .true.
     .            !Standard relax
     .            ,mg_smooth     = sm_type
     .            ,sm_omega      = sm_omega
     .            ,sm_it         = sm_iter
     .            ,sm_zebra_relax= mg_zebra_relax
     .            !SI relax
cc     .            ,mg_smooth     = sm_type_si
cc     .            ,sm_omega      = sm_omega_si
cc     .            ,sm_it         = sm_iter_si
cc     .            ,sm_zebra_relax= mg_zebra_relax_si
     .            )

c diag ****
cc      !Calculate residual
cc      call avlap_mtvc(0,3,3*ntotdp,x,rhs,igx,bcsa)
cc      x = rhs
c diag ****

      if (is_a_cnv) then
        call fillArray(igx,IAX,3,x,da_cnv,bcsa,bc_order_a
     .                ,is_vec=.true.,is_cnv=is_a_cnv
     .                ,arr_cov=da_cov)
      else
        call fillArray(igx,IAX,3,x,da_cov,bcsa,bc_order_a
     .                ,is_vec=.true.,is_cnv=is_a_cnv)
      endif

c     Find new dB (with BCs)

      db = curl(igx,da_cov)

      if (wbcs) then
        call setMGBC(0,3,nx,ny,nz,igx,db,bcsb
     .              ,icomp=(/IBX/),is_vec=.true.
     .              ,is_cnv=.true.,iorder=bc_order)
      endif

c     Diagnostics

      if (dplots) then
        dcf='divcl_p'//trim(int2char(my_rank))//'.bin'

        call createDrawInCfile(6,trim(dcf),'Solution'
     .          ,'t','x','y',(/ 'Ax','Ay','Az','Bx','By','Bz'/)
     .          ,'-c -X0 -L57'
     .          ,'drawdc_p'//trim(int2char(my_rank))//'.in')

        open(unit=110,file=trim(dcf),form='unformatted'
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

        call pstop('vlap_div_clean','Dumped graphics')
      endif

c     Store vector potential

      if (PRESENT(dacov)) dacov = da_cov

      end subroutine vlap_div_clean

c     avlap_mtvc
c     ###############################################################
      subroutine avlap_mtvc(gpos,neq,ntot,x,y,igrid,bcnd)
c     ---------------------------------------------------------------
c     This subroutine calculates lap(vector) in contravariant form.
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
c     ---------------------------------------------------------------

      use matvec

      implicit none

c     Call variables

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer    :: imin,imax,jmin,jmax,kmin,kmax,igx,igy,igz,ieq
      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

      real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm,jac
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
     .             ,veclap(3),cnv(3),cov(3),car(3)
      real(8)    :: upwind,etal

      real(8),allocatable,dimension(:,:,:,:) :: vv,v0

      logical    :: fpointers,is_cnv

c     Begin program

      is_cnv = is_a_cnv

      call allocPointers(neq,fpointers)

      isig = MGgrid%istartp(igrid)

      nxx = grid_params%nxv(igrid)
      nyy = grid_params%nyv(igrid)
      nzz = grid_params%nzv(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

c     Find limits for loops

      call limits(abs(gpos),nxx,nyy,nzz
     .           ,igrid,imin,imax,jmin,jmax,kmin,kmax)

c     Map vector x to array for processing

      allocate(vv(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      vv = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call mapMGVectorToArray(max(0,gpos),neq,x,nxx,nyy,nzz,vv,igrid
     .                       ,.false.)

      call setMGBC(max(0,gpos),neq,nxx,nyy,nzz,igrid,vv,bcnd
     .            ,icomp=(/IAX/),is_cnv=is_cnv,is_vec=.true.
     .            ,result_is_vec=.true.,iorder=bc_order_a)

c     Calculate matrix-vector product

      if (is_a_cnv) then
        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

              veclap = veclaplacian(i,j,k,nxx,nyy,nzz
     .                             ,igrid,igrid,igrid,vv
     .                             ,vol=volf)
cc              veclap =-curlcurl(i,j,k,nxx,nyy,nzz
cc       .                           ,igrid,igrid,igrid,vv
cc       .                           ,vol=volf)

              do ieq=1,3
                y(neq*(ijk-1)+ieq) = veclap(ieq)
              enddo

            enddo
          enddo
        enddo
      else
        do k = kmin,kmax
          do j = jmin,jmax
            do i = imin,imax
              ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

              veclap = veclap_cov(i,j,k,nxx,nyy,nzz
     .                           ,igrid,igrid,igrid,vv
     .                           ,vol=volf)

              do ieq=1,3
                y(neq*(ijk-1)+ieq) = veclap(ieq)
              enddo

            enddo
          enddo
        enddo
      endif

c     End program

      deallocate(vv)

      call deallocPointers(fpointers)

      end subroutine avlap_mtvc

      end module B_tools
