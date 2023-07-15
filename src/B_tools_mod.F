c module B_tools
c ####################################################################
      module B_tools

        use setMGBC_interface

        use PC_var

        use local_BCS_variables, ONLY:order_bc
     .                               ,default_A_BCs
     .                               ,default_B_BCs

        use app_iosetup

        use xdraw_io

        integer,private :: bc_order_a

      contains

c     curl_inv
c     ###############################################################
      subroutine curl_inv(g_def,igrid,b,a,global,w_bcs,order,sp)

c     ---------------------------------------------------------------
c     Finds COVARIANT vector potential "a" from CONTRAVARIANT vector
c     field components "b". Employs gauge A_1 = 0d0. In SP systems,
c     it integrates at faces in r, and then averages to nodes.
c     
c     Integrals done by high-order collocated integration formulas.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer    :: igrid
      real(8)    :: a(0:,0:,0:,:)
     .             ,b(0:,0:,0:,:)
      logical,optional :: global,w_bcs,sp
      integer,optional :: order

c     Local variables

      integer    :: nx,ny,nz,igx,igy,igz,lorder
      integer    :: i,j,k,ii,bcnd(6,3),ig,jg,kg,ivar
      real(8)    :: dxx,dyy,dvol,cov(3)

      logical    :: spoint,glbl,do_bc

c     Begin program

      igx = igrid ; igy = igrid ; igz = igrid

      nx = size(a,1)-2
      ny = size(a,2)-2
      nz = size(a,3)-2

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

      if (PRESENT(sp)) then
        spoint = sp
      else
        spoint = bcSP()
      endif

      if (PRESENT(order)) then
        lorder = order
      else
        lorder = 3
      endif

c     Consistency check

      if (igrid /= 1 .and. glbl) then
        call pstop('curl_inv','Inconsistent grid level')
      endif

c     Gauge

      a(:,:,:,1)=0d0

c     Quadrant-based code (SP ready, parallel-ready)

      if (spoint) then  !Separate quadrants only in theta

        if (ny == 1) then
          call line_int(1,nx,1,ny,b,a,lorder)
        else
          !Lower-right quadrant
          call line_int(1,nx,ny/2-1,1
     .               ,b(0:nx+1,0:ny/2,0:nz+1,1:3)
     .               ,a(0:nx+1,0:ny/2,0:nz+1,1:3),lorder)

          !Upper-right quadrant
          call line_int(1,nx,ny/2+1,ny
     .               ,b(0:nx+1,ny/2:ny+1,0:nz+1,1:3)
     .               ,a(0:nx+1,ny/2:ny+1,0:nz+1,1:3),lorder)
        endif
      else
        !Lower-left quadrant
        call line_int(nx/2-1,1,ny/2-1,1
     .               ,b(0:nx/2,0:ny/2,0:nz+1,1:3)
     .               ,a(0:nx/2,0:ny/2,0:nz+1,1:3),lorder)

        !Lower-right quadrant
        call line_int(nx/2+1,nx,ny/2-1,1
     .               ,b(nx/2:nx+1,0:ny/2,0:nz+1,1:3)
     .               ,a(nx/2:nx+1,0:ny/2,0:nz+1,1:3),lorder)

        !Upper-left quadrant
        call line_int(nx/2-1,1,ny/2+1,ny
     .               ,b(0:nx/2,ny/2:ny+1,0:nz+1,1:3)
     .               ,a(0:nx/2,ny/2:ny+1,0:nz+1,1:3),lorder)

        !Upper-right quadrant
        call line_int(nx/2+1,nx,ny/2+1,ny
     .               ,b(nx/2:nx+1,ny/2:ny+1,0:nz+1,1:3)
     .               ,a(nx/2:nx+1,ny/2:ny+1,0:nz+1,1:3),lorder)
      endif

      !Boundary conditions: topological constraints (SP, PER)
      if (do_bc) then
        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond

        call setMGBC(g_def,0,3,nx,ny,nz,igrid,a,bcnd,icomp=(/IAX/)
     .              ,is_vec=.true.,is_cnv=.false.,iorder=min(lorder,2))
      endif

c     End program

      contains

c     line_int
c     #############################################################
      subroutine line_int(ilo,ihi,jlo,jhi,bb,aa,lorder)

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
        integer :: lorder
        
c     Local variables

        real(8) :: r1,r2,r3,jac
        integer :: i,im,ip,ipp,imm,im3,ip3,imin,imax,istep
     .            ,j,jm,jp,jpp,jmm,jm3,jp3,jmin,jmax,jstep

        !3rd order
        real(8),parameter :: cm=-1d0/12d0       !Upwinded
     .                      ,c0=2d0/3d0 
     .                      ,cp=5d0/12d0
     .                      ,b1=6250d0/23040d0  !@ SP
     .                      ,b2=185d0 /23040d0
     .                      ,b3=13d0  /23040d0

        !4th order
        real(8),parameter :: ccm=-1d0/24d0  !Centered
     .                      ,cc0=13d0/24d0
     .                      ,smm= 9d0/24d0  !One-sided
     .                      ,sm =19d0/24d0
     .                      ,sp =-5d0/24d0
     .                      ,spp= 1d0/24d0
     .                      ,a1=512667246d0/1.8579456d9 !@ SP
     .                      ,a2=20337324d0 /1.8579456d9
     .                      ,a3=3062796d0  /1.8579456d9
     .                      ,a4=388149d0   /1.8579456d9
     .                      ,a5=26021d0    /1.8579456d9

c     Begin program

        istep = -1
        if (ilo < ihi) istep = 1

        jstep = -1
        if (jlo < jhi) jstep = 1

        aa(ilo-istep,jlo-jstep,:,3)=0d0

        i = ilo-istep

        jmin = jlo-min(jstep,0)
        jmax = jhi+max(jstep,0)

        do j=jmin,jmax,jstep
          jp  = j  + min (0,jstep)
          jm  = j  - max (0,jstep)
          jpp = jp + sign(1,jstep)
          jp3 = jp + sign(2,jstep)
          jmm = jm - sign(1,jstep)
          jm3 = jm - sign(2,jstep)

          if (.not.glbl) then
            call getMGmap(g_def,i,jm,1,igx,igy,igz,ig,jg,kg)
            dyy = sign(1,jstep)*g_def%dy(jg)
          else
            dyy = g_def%yg(jp)
     .           -g_def%yg(jm)
          endif

          select case(lorder)
          case(2)
            aa(i,jp,:,3)=aa(i,jm,:,3) + dyy*5d-1*(bb(i,jm,:,1)
     .                                           +bb(i,jp,:,1))
          case(3)
            if(j == jmin) then
              aa(:,jp,:,3)=aa(:,jm,:,3)
     .                      +dyy*(cp*bb(:,jm ,:,1)
     .                           +c0*bb(:,jp ,:,1)
     .                           +cm*bb(:,jpp,:,1))
            else                              
              aa(:,jp,:,3)=aa(:,jm,:,3)       
     .                      +dyy*(cm*bb(:,jmm,:,1)
     .                           +c0*bb(:,jm ,:,1)
     .                           +cp*bb(:,jp ,:,1))
            endif
          case(4)
            if(j == jmin) then
              aa(:,jp,:,3)=aa(:,jm,:,3)
     .                      +dyy*(smm*bb(:,jm ,:,1)
     .                           +sm *bb(:,jp ,:,1)
     .                           +sp *bb(:,jpp,:,1)
     .                           +spp*bb(:,jp3,:,1))
            elseif(j == jmax) then
              aa(:,jp,:,3)=aa(:,jm,:,3)
     .                      +dyy*(spp*bb(:,jm3,:,1)
     .                           +sp *bb(:,jmm,:,1)
     .                           +sm *bb(:,jm ,:,1)
     .                           +smm*bb(:,jp ,:,1))
            else                              
              aa(:,jp,:,3)=aa(:,jm,:,3)       
     .                      +dyy*(ccm*(bb(:,jmm,:,1)
     .                                +bb(:,jpp,:,1))
     .                           +cc0*(bb(:,jm ,:,1)
     .                                +bb(:,jp ,:,1)))
            endif
          end select
          
        enddo

        aa(ilo-istep,:,:,2)=0d0

        imin = ilo-min(istep,0)
        imax = ihi+max(istep,0)

        do i=imin,imax,istep

          ip  = i  + min (0,istep)
          im  = i  - max (0,istep)
          ipp = ip + sign(1,istep)
          imm = im - sign(1,istep)
          ip3 = ip + sign(2,istep)
          im3 = im - sign(2,istep)

          if (spoint.and.i==1) then
            if (.not.glbl) then
              call getMGmap(g_def,im,1,1,igx,igy,igz,ig,jg,kg)
              dxx = sign(1,istep)*g_def%dxh(ig)
            else
              dxx = g_def%xg(ip)
     .             -g_def%xg(im)
            endif

            select case(lorder)
            case(1)
              aa(ip,:,:,3) =  - 0.25*dxx*bb(ip,:,:,2) !A_3=0 @ r=0
              aa(ip,:,:,2) =  + 0.25*dxx*bb(ip,:,:,3) !A_2=0 @ r=0
            case(2)
              !A_3=0 @ r=0
              aa(ip,:,:,3) =-dxx*(51*bb(ip  ,:,:,2)
     .                              -bb(ip+1,:,:,2))/192d0
              !A_2=0 @ r=0
              aa(ip,:,:,2) = dxx*(51*bb(ip  ,:,:,3)
     .                              -bb(ip+1,:,:,3))/192d0
            case(3)
              !A_3=0 @ r=0
              aa(ip,:,:,3) =-dxx*(b1*bb(ip  ,:,:,2)
     .                           -b2*bb(ip+1,:,:,2)
     .                           +b3*bb(ip+2,:,:,2))
              !A_2=0 @ r=0
              aa(ip,:,:,2) = dxx*(b1*bb(ip  ,:,:,3)
     .                           -b2*bb(ip+1,:,:,3)
     .                           +b3*bb(ip+2,:,:,3))
            case(4)
              aa(ip,:,:,3) =-dxx*(a1*bb(ip  ,:,:,2)
     .                           -a2*bb(ip+1,:,:,2)
     .                           +a3*bb(ip+2,:,:,2)
     .                           -a4*bb(ip+3,:,:,2)
     .                           +a5*bb(ip+4,:,:,2))
              aa(ip,:,:,2) = dxx*(a1*bb(ip  ,:,:,3)
     .                           -a2*bb(ip+1,:,:,3)
     .                           +a3*bb(ip+2,:,:,3)
     .                           -a4*bb(ip+3,:,:,3)
     .                           +a5*bb(ip+4,:,:,3))
            end select

          else

            if (.not.glbl) then
              call getMGmap(g_def,im,1,1,igx,igy,igz,ig,jg,kg)
              dxx = sign(1,istep)*g_def%dx(ig)
            else
              dxx = g_def%xg(ip)
     .             -g_def%xg(im)
            endif

            select case(lorder)
            case(2)
              aa(ip,:,:,3) = aa(im,:,:,3)
     .                     - dxx*5d-1*(bb(im,:,:,2)
     .                                +bb(ip,:,:,2))
              aa(ip,:,:,2) = aa(im,:,:,2)
     .                     + dxx*5d-1*(bb(im,:,:,3)
     .                                +bb(ip,:,:,3))
            case(3)
              if(i == imin) then
                aa(ip,:,:,3)=aa(im,:,:,3)
     .                      -dxx*(cp*bb(im ,:,:,2)
     .                           +c0*bb(ip ,:,:,2)
     .                           +cm*bb(ipp,:,:,2))
                aa(ip,:,:,2)=aa(im,:,:,2)
     .                      +dxx*(cp*bb(im ,:,:,3)
     .                           +c0*bb(ip ,:,:,3)
     .                           +cm*bb(ipp,:,:,3))
              else
                aa(ip,:,:,3)=aa(im,:,:,3)
     .                      -dxx*(cm*bb(imm,:,:,2)
     .                           +c0*bb(im ,:,:,2)
     .                           +cp*bb(ip ,:,:,2))
                aa(ip,:,:,2)=aa(im,:,:,2)
     .                      +dxx*(cm*bb(imm,:,:,3)
     .                           +c0*bb(im ,:,:,3)
     .                           +cp*bb(ip ,:,:,3))
              endif
            case(4)
              if(i == imin) then
                aa(ip,:,:,3)=aa(im,:,:,3)
     .                      -dxx*(smm*bb(im ,:,:,2)
     .                           +sm *bb(ip ,:,:,2)
     .                           +sp *bb(ipp,:,:,2)
     .                           +spp*bb(ip3,:,:,2))
                aa(ip,:,:,2)=aa(im,:,:,2)
     .                      +dxx*(smm*bb(im ,:,:,3)
     .                           +sm *bb(ip ,:,:,3)
     .                           +sp *bb(ipp,:,:,3)
     .                           +spp*bb(ip3,:,:,3))
              elseif(i == imax) then
                aa(ip,:,:,3)=aa(im,:,:,3)
     .                      -dxx*(spp*bb(im3,:,:,2)
     .                           +sp *bb(imm,:,:,2)
     .                           +sm *bb(im ,:,:,2)
     .                           +smm*bb(ip ,:,:,2))
                aa(ip,:,:,2)=aa(im,:,:,2)
     .                      +dxx*(spp*bb(im3,:,:,3)
     .                           +sp *bb(imm,:,:,3)
     .                           +sm *bb(im ,:,:,3)
     .                           +smm*bb(ip ,:,:,3))
              else
                aa(ip,:,:,3)=aa(im,:,:,3)
     .                      -dxx*(ccm*(bb(imm,:,:,2)
     .                                +bb(ipp,:,:,2))
     .                           +cc0*(bb(im ,:,:,2)
     .                                +bb(ip ,:,:,2)))
                aa(ip,:,:,2)=aa(im,:,:,2)
     .                      +dxx*(ccm*(bb(imm,:,:,3)
     .                                +bb(ipp,:,:,3))
     .                           +cc0*(bb(im ,:,:,3)
     .                                +bb(ip ,:,:,3)))
              endif
            end select
            
          endif
        enddo

      end subroutine line_int

      end subroutine curl_inv

c     curl_div_clean
c     #################################################################
      subroutine curl_div_clean(g_def,igx,db,diag_plots,global,order)

c     -----------------------------------------------------------------
c     Finds B^(n+1)=B^n + curl(dA), with dA=inv_curl(dB)
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igx
      real(8) :: db(0:,0:,0:,:)

      logical :: diag_plots
      logical,optional :: global
      integer,optional :: order

c     Local variables

      real(8) :: da_cov(0:size(db,1)-1
     $                 ,0:size(db,2)-1
     $                 ,0:size(db,3)-1
     $                 ,size(db,4))

c     Find dA from dB (w/ ghost cells)

      call curl_inv(g_def,igx,db,da_cov,global=global,order=order)

c     Diagnostics (call before db = curl(A))

      if (diag_plots) call B_dump(g_def,igx,da_cov,db)

c     Find dB

      db = curl(g_def,igx,da_cov)  !This step is NOT global!!

      end subroutine curl_div_clean

c     vlap_curl_inv
c     ###############################################################
      subroutine vlap_curl_inv(g_def,igrid,db_cnv,da_cov,bcsb,bcsa
     .                        ,iout,dbcov,da0cov,w_pre_bcs,gm_it,rtol)

c     ---------------------------------------------------------------
c     Finds COVARIANT vector potential "a" from CONTRAVARIANT vector
c     field components "b". Employs gauge div(A) = 0d0. Uses vector
c     Laplacian.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid,iout,bcsb(6,3),bcsa(6,3)
      real(8) :: db_cnv(0:,0:,0:,:)
      real(8) :: da_cov(0:,0:,0:,:)

      integer,optional :: gm_it
      real(8),optional :: rtol
      real(8),optional,INTENT(IN) :: dbcov(0:size(db_cnv,1)-1
     .                                    ,0:size(db_cnv,2)-1
     .                                    ,0:size(db_cnv,3)-1
     .                                    ,  size(db_cnv,4))
      real(8),optional,INTENT(IN) :: da0cov(0:size(db_cnv,1)-1
     .                                     ,0:size(db_cnv,2)-1
     .                                     ,0:size(db_cnv,3)-1
     .                                     ,  size(db_cnv,4))

      logical,optional :: w_pre_bcs

c     Local variables

      integer :: nx,ny,nz,igx,igy,igz,i,j,k,ii,icomp,order_res
     .          ,order_prol,gmits

      real(8) :: da_cnv (0:size(db_cnv,1)-1
     .                  ,0:size(db_cnv,2)-1
     .                  ,0:size(db_cnv,3)-1,size(db_cnv,4))
     .          ,dlapa0 (0:size(db_cnv,1)-1
     .                  ,0:size(db_cnv,2)-1
     .                  ,0:size(db_cnv,3)-1,size(db_cnv,4))
     .          ,db_cov (0:size(db_cnv,1)-1
     .                  ,0:size(db_cnv,2)-1
     .                  ,0:size(db_cnv,3)-1,size(db_cnv,4))
     .          ,dj_cnv (0:size(db_cnv,1)-1
     .                  ,0:size(db_cnv,2)-1
     .                  ,0:size(db_cnv,3)-1,size(db_cnv,4))
     .          ,rhs(ntotdp,3),x(ntotdp,3),dvol,etol
     $          ,lxmin,lxmax,lymin,lymax,lzmin,lzmax

      logical :: vwgt,wbcs,w_prebcs,sing_mat

c     Begin program
            
      vwgt = .true.

      igx = igrid
      igy = igrid
      igz = igrid

      nx = size(db_cnv,1)-2
      ny = size(db_cnv,2)-2
      nz = size(db_cnv,3)-2

      if (PRESENT(w_pre_bcs)) then
        w_prebcs = w_pre_bcs
      else
        w_prebcs = .false.
      endif

      if (PRESENT(rtol)) then
        etol = rtol
      else
        etol = 1d-7
      endif

c     Find Lap(A0) for Dirichlet BCs

      dlapa0 = 0d0

      if (PRESENT(da0cov)) then
        bc_order_a = 1          !To ensure consistency with Dirichlet BC provided
        !Extrapolate from face to ghost cell
        call setMGBC(g_def,0,3,nx,ny,nz,igrid,dlapa0,bcsa
     .              ,icomp=(/IAX/),arr0=da0cov,is_vec=.true.
     .              ,is_cnv=.false.,iorder=bc_order_a)
c     =================
c$$$        dlapa0 = veclap_cov(g_def,igrid,dlapa0) !Returns cov
c$$$        dlapa0 = XformToCnv(g_def,igrid,dlapa0)
c     =================
        dlapa0 = XformToCnv(g_def,igrid,dlapa0)
        dlapa0 = veclap(g_def,igrid,dlapa0)  !Returns cnv
      else
        bc_order_a = 2
      endif

c     Find rhs of vec_lap solve: rhs=-vol*dJ (with dJ=curl(dB); dB already has BCs)

      if (w_prebcs) then
        call setMGBC(g_def,0,3,nx,ny,nz,igx,db_cnv,bcsb
     .              ,arr_cov=db_cov,icomp=(/IBX/),is_vec=.true.
     .              ,is_cnv=.true.,iorder=order_bc)
      else
        if (PRESENT(dbcov)) then
          db_cov = dbcov
        else
          db_cov = XformVector(g_def,igx,db_cnv,.false.)
        endif
      endif

      dj_cnv = curl(g_def,igrid,db_cov) + dlapa0

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            ii  = i + nx*(j-1) + nx*ny*(k-1)

            if (vwgt) then
              dvol = g_def%gmetric%grid(igx)%dlvol(i,j,k)
            else
              dvol = 1d0
            endif

            rhs(ii,:) = -dvol*dj_cnv(i,j,k,:)
          enddo
        enddo
      enddo

c     Setup solver parameters

      sing_mat = .not.have_map()

      if (sing_mat) then
        order_res  = 0         !To ensure conservation (avoid null space)
        order_prol = 1
      else
        order_res  = 0
        order_prol = 1
      endif
      
c     Invert vector laplacian to find dA=vlap^-1(dJ)

      x = 0d0

      call cSolver(g_def,3,ntotdp,rhs,x,bcsa,igx
     .            ,iout,1,avlap_mtvc,vwgt
     .            ,tol           = etol
     .            ,iters         = gmits
     .            ,gm_driver     = .true.
     .            ,sing_matrix   = sing_mat
     .            ,ks_it         = 100
     .            ,mg_vcyc       = 1
     .            ,mg_order_res  = order_res
     .            ,mg_order_prol = order_prol
cc     .            ,mg_coarse_grid_size = mg_coarse_size  !Deprecated
     .            ,mg_gm_coarse_solve  = mg_gm_coarse
cc     .            ,mg_gm_coarse_solve  = .true.
     .            ,sm_ncolors    = sm_ncolors
     .            ,sm_line_relax = sm_line_relax
     .            !Zebra relax
cc     .            ,mg_smooth     = 'gs'
cc     .            ,sm_omega      = 0.9d0
cc     .            ,sm_it         = 2
cc     .            ,sm_zebra_relax= .true.
     .            !Standard relax
     .            ,mg_smooth     = sm_type
     .            ,sm_omega      = sm_omega
     .            ,sm_it         = sm_iter
     .            ,sm_zebra_relax= sm_zebra_relax
     .            !SI relax
cc     .            ,mg_smooth     = sm_type_si
cc     .            ,sm_omega      = sm_omega_si
cc     .            ,sm_it         = sm_iter_si
cc     .            ,sm_zebra_relax= mg_zebra_relax_si
     .            )
      
      call fillArray(g_def,igx,3,x,da_cnv,bcsa,bc_order_a
     .              ,is_vec=.true.,is_cnv=.true.,icomp=(/IAX/)
     .              ,arr_cov=da_cov,arr0=da0cov)

      if (PRESENT(gm_it)) gm_it = gmits

c     End program

      end subroutine vlap_curl_inv
      
c     vlap_div_clean
c     #################################################################
      subroutine vlap_div_clean(g_def,iout,nx,ny,nz,igx,igy,igz,db,dacov
     .              ,da0cov,dbcov,w_bcs,w_pre_bcs,diag_plots,gm_it,rtol)

c     -----------------------------------------------------------------
c     Finds curl(dA)=B^(n+1)-B^n, with dA=inv_lap(dJ). 
c     On input:
c       * iout: level of output info
c       * nx,ny,nz: mesh size
c       * igx,igy,igz: grid level size
c       * db: input cnv magnetic field 
c       * dacov,da0cov,dbcov (optional): cov A0 and B (input) and A (output)
c       * w_bcs (optional): whether output has BCs
c       * w_pre_bcs (optional): whether input requires BCs
c       * diag_plots (optional): whether we want plots of answer
c       * gm_it,rtol (optional): number of GMRES iterations and conv. tol.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: nx,ny,nz,igx,igy,igz,iout
      real(8) :: db(0:nx+1,0:ny+1,0:nz+1,3)

      integer,optional :: gm_it
      real(8),optional :: rtol
      real(8),optional,INTENT(IN) :: dbcov (0:nx+1,0:ny+1,0:nz+1,3)
     $                              ,da0cov(0:nx+1,0:ny+1,0:nz+1,3)
      real(8),optional,INTENT(OUT):: dacov (0:nx+1,0:ny+1,0:nz+1,3)

      logical,optional :: w_bcs,w_pre_bcs,diag_plots

c     Local variables

      integer :: bcsa(6,3),bcsb(6,3)

      real(8) :: da_cov (0:nx+1,0:ny+1,0:nz+1,3)
      logical :: dplots,wbcs
      
c     Begin program

      if (my_rank == 0.and.iout > 0) then
        write (*,*)
        call info('  Div-cleaning B-field...')
      endif

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

c     Define BCs

      bcsb(:,1) = u_0%array_var(IBX)%bconds
      bcsb(:,2) = u_0%array_var(IBY)%bconds
      bcsb(:,3) = u_0%array_var(IBZ)%bconds
c$$$      bcsb = gv%aux%vec_list(IBCNV_0)%bconds

      call default_A_BCs(bcsa)

c     Find dA

      call vlap_curl_inv(g_def,igx,db,da_cov,bcsb,bcsa,iout
     .                  ,w_pre_bcs=w_pre_bcs,da0cov=da0cov
     .                  ,gm_it=gm_it,rtol=rtol)

c     Diagnostics

      if (dplots) call B_dump(g_def,igx,da_cov,db)

c     Find new dB (with BCs)

      db = curl(g_def,igx,da_cov)

      if (wbcs) then
        call setMGBC(g_def,0,3,nx,ny,nz,igx,db,bcsb
     .              ,icomp=(/IBX/),is_vec=.true.
     .              ,is_cnv=.true.,iorder=order_bc)
      endif

c     Store vector potential

      if (PRESENT(dacov)) dacov = da_cov

      if (my_rank == 0.and.iout > 0) then
        call info('  Done!')
        write (*,*)
      endif

      end subroutine vlap_div_clean

c     B_div_clean
c     #################################################################
      subroutine B_div_clean(g_def,iout,igrid,db,meth,dacov
     .                      ,w_bcs,w_pre_bcs,diag_plots,gm_it,rtol)

c     -----------------------------------------------------------------
c     Finds curl(dA)=B^(n+1)-B^n, with dA=inv_lap(dJ). 
c     On input:
c       * iout: level of output info
c       * igrid: grid level size
c       * db: input cnv magnetic field 
c       * dacov (optional): cov A (output)
c       * meth: divergence cleaning method ("vlap","line")
c       * w_bcs (optional): whether output has BCs
c       * w_pre_bcs (optional): whether input requires BCs
c       * diag_plots (optional): whether we want plots of answer
c       * gm_it,rtol (optional): number of GMRES iterations and conv. tol.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: iout,igrid
      real(8) :: db(0:,0:,0:,:)

      character(4) :: meth

      integer,optional :: gm_it
      real(8),optional :: rtol
      real(8),optional,INTENT(OUT):: dacov(0:,0:,0:,:)

      logical,optional :: w_bcs,w_pre_bcs,diag_plots

c     Local variables

      integer :: bcsb(6,3),nx,ny,nz,its

      real(8) :: da_cov (0:size(db,1)-1
     .                  ,0:size(db,2)-1
     .                  ,0:size(db,3)-1,3)
      logical :: dplots,wbcs

c     Begin program

      if (my_rank == 0.and.iout > 0) then
        write (*,*)
        write (*,*) 'Div-cleaning B-field...'
        write (*,*)
      endif

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

c     Find vector potential

      its = 0

      call findVectorPotential(g_def,iout,igrid,db,da_cov,meth
     .                        ,its=its,rtol=rtol)

      if (my_rank == 0.and.iout > 0.and.(.not.PRESENT(gm_it))) then
        write (*,'(a,i3,a)') meth//' divergence cleaning converged in'
     .                      ,its,' its'
      endif

      if (PRESENT(gm_it)) gm_it = its

c     Diagnostics

      if (dplots) call B_dump(g_def,igrid,da_cov,db)

c     Find new dB (with BCs)

      db = curl(g_def,igrid,da_cov)

      if (wbcs) then
        call default_B_BCs(bcsb)

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        call setMGBC(g_def,0,3,nx,ny,nz,igrid,db,bcsb
     .              ,icomp=(/IBX/),is_vec=.true.
     .              ,is_cnv=.true.,iorder=order_bc)
      endif

c     Store vector potential

      if (PRESENT(dacov)) dacov = da_cov

      end subroutine B_div_clean

c     avlap_mtvc
c     ###############################################################
      subroutine avlap_mtvc(g_def,gpos,neq,ntot,x,y,igrid,bcnd)
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

      type(grid_mg_def),pointer :: g_def

      integer    :: neq,ntot,igrid,gpos,bcnd(6,neq)
      real(8)    :: x(ntot),y(ntot)

c     Local variables

      integer    :: i,j,k,ig,jg,kg,isig,ip,im,jp,jm,kp,km,nxx,nyy,nzz
      integer    :: imin,imax,jmin,jmax,kmin,kmax,igx,igy,igz,ieq
      integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

      real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm,jac
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,vol
     .             ,vlap(3),cnv(3),cov(3),car(3)
      real(8)    :: upwind,etal

      real(8),pointer,dimension(:,:,:,:) :: vv

      logical    :: fpointers

      type(mg_ctx),pointer :: mg_ctx_l => null()

c     Begin program

      call allocMGPointers(neq,g_def,mg_ctx_l,fpointers)

      isig = g_def%istartp(igrid)

      nxx = g_def%nxv(igrid)
      nyy = g_def%nyv(igrid)
      nzz = g_def%nzv(igrid)

      igx = igrid
      igy = igrid
      igz = igrid

c     Find limits for loops

      call limits(mg_ctx_l,abs(gpos)
     .           ,igrid,imin,imax,jmin,jmax,kmin,kmax)

c     Map vector x to array for processing

      allocate(vv(0:nxx+1,0:nyy+1,0:nzz+1,neq))

      vv = 0d0

      !For GS, gpos < 0 so that the whole vector x is mapped and BCs are filled
      !For finding the diagonal, gpos > 0
      call fillArray(g_def,igrid,neq,x,vv,bcnd,bc_order_a,gpos=gpos
     .              ,icomp=(/IAX/),is_cnv=.true.,is_vec=.true.
     .              ,result_is_curv=.true.)

c     Calculate matrix-vector product

       do k = kmin,kmax
         do j = jmin,jmax
           do i = imin,imax
             ijk    = i + nxx*(j-1) + nxx*nyy*(k-1)

             vlap = veclap    (g_def,i,j,k,nxx,nyy,nzz
cc             vlap = veclap_cnv(g_def,i,j,k,nxx,nyy,nzz
     .                        ,igrid,igrid,igrid,vv
     .                        ,vol=MG_get_res_vol_wgt())

             do ieq=1,3
               y(neq*(ijk-1)+ieq) = vlap(ieq)
             enddo
           enddo
         enddo
       enddo

c     End program

      deallocate(vv)

      call deallocMGPointers(fpointers,mg_ctx_l)

      end subroutine avlap_mtvc

ccc     grad_div_clean
ccc     #################################################################
cc      subroutine grad_div_clean(db)
cc
ccc     -----------------------------------------------------------------
ccc     Finds dB=dB* + grad(phi).
ccc     -----------------------------------------------------------------
cc
cc      use matvec, bcs_pc => bcs
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      real(8) :: db(0:nx+1,0:ny+1,0:nz+1,3)
cc
ccc     Local variables
cc
cc      integer :: ii,bcs(6)
cc      real(8) :: rhs(ntotdp),x(ntotdp),phi(0:nx+1,0:ny+1,0:nz+1)
cc     .          ,mag,dvol
cc
cc      logical :: vwgt
cc
cc      external lap_mtvc
cc
ccc     Begin program
cc
cc      vwgt = vol_wgt
cc
ccc     Find rhs=-div(db)
cc
cc      rhs = 0d0
cc
cc      do k=1,nz
cc        do j=1,ny
cc          do i=1,nx
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc            if (vwgt) then
cc              dvol = gmetric%grid(igx)%dvol(i,j,k)
cc            else
cc              dvol = 1d0
cc            endif
cc
cc            rhs(ii)= -dvol*div(i,j,k,nx,ny,nz,igx,igy,igz,db)
cc          enddo
cc        enddo
cc      enddo
cc
ccc     Solve for phi=lap^-1(rhs)
cc
cc      bcs = bcond
cc      where (bcond == DEF) bcs = NEU  !B.n must be zero at boundary!!
cc
cc      x = 0d0
cc
cc      !Solve Poisson equation
cc      call cSolver(1,ntotdp,rhs,x,bcs
cc     .            ,1,ilevel,0,lap_mtvc,vwgt
cc     .            ,tol           = 1d-5
cc!     .            ,dg            = lap_diag
cc     .            ,gm_driver     = .not.vwgt
cc     .            ,cg_driver     = vwgt
cc     .            ,ks_it         = 20
cc     .            ,mg_debug      = .false.
cc     .            ,mg_galerkin   = .false.
cc     .            ,mg_order_res  = 0  !To avoid null spaces
cc     .            ,mg_order_prol = 0  !To avoid null spaces
cc     .            ,mg_coarse_grid_size = 2
cc     .            ,mg_gm_coarse_solve  = .false.
cc     .            ,mg_vcyc       = 1
cc     .            ,mg_smooth     = 'gs'
cc     .            ,sm_omega      = 1d0
cc     .            ,sm_it         = 3
cc     .            ,sm_ncolors    = sm_ncolors  !9-pt stencil
cc     .            ,sm_line_relax = sm_line_relax
cc     .            ,sm_zebra_relax= sm_zebra_relax
cc     .            )
cc
cc      !Map solution to array and apply BCs
cc      call mapMGVectorToArray(0,1,x,nx,ny,nz,phi,igx,.false.)
cc
cc      call setMGBC(0,1,nx,ny,nz,igx,phi,bcs,iorder=bc_order)
cc
ccc     Find dB=-grad(phi) (in cnv representation and w/o ghost cells)
cc
cc      do k = 1,nz
cc        do j = 1,ny
cc          do i = 1,nx
cc            db(i,j,k,:) = db(i,j,k,:)
cc     .           + matmul(gmetric%grid(igx)%gsup(i,j,k,:,:)
cc     .                   ,grad(i,j,k,nx,ny,nz,igx,igy,igz,phi))
cc          enddo
cc        enddo
cc      enddo
cc
cc      end subroutine grad_div_clean

c     B_dump
c     ###############################################################
      subroutine B_dump(g_def,igx,da_cov,db)

c     ---------------------------------------------------------------
c     Dumps magnetic and vector potential fields
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igx
      real(8) :: da_cov(0:,0:,0:,:)
     .          ,db    (0:,0:,0:,:)

c     Local variables

      integer :: nx,ny,nz
      real(8) :: ddb(0:size(db,1)-1
     $              ,0:size(db,2)-1
     $              ,0:size(db,3)-1
     $              ,size(db,4))

      character(20) :: dcf

c     Begin program

      nx = g_def%nxv(igx)
      ny = g_def%nyv(igx)
      nz = g_def%nzv(igx)

      dcf='divcl_p'//trim(int2char(my_rank))//'.bin'

      call createDrawInCfile(9,trim(dcf),'Solution','t','x','y'
     .        ,(/'Bx_in ','By_in ','Bz_in '
     $          ,'Ax_out','Ay_out','Az_out'
     $          ,'Bx_dff','By_dff','Bz_dff'/)
     .        ,'-c -X0 -L57','drawdc_p'//trim(int2char(my_rank))//'.in')

      open(unit=110,file=trim(dcf),form='unformatted'
     .    ,status='replace')

      call contour(db   (:,:,1,1),nx+2,ny+2
     .            ,g_def%gxmin,g_def%gxmax
     .            ,g_def%gymin,g_def%gymax
     .            ,0,110)
      call contour(db   (:,:,1,2),nx+2,ny+2
     .            ,g_def%gxmin,g_def%gxmax
     .            ,g_def%gymin,g_def%gymax
     .            ,1,110)
      call contour(db   (:,:,1,3),nx+2,ny+2
     .            ,g_def%gxmin,g_def%gxmax
     .            ,g_def%gymin,g_def%gymax
     .            ,1,110)

      call contour(da_cov(:,:,1,1),nx+2,ny+2
     .            ,g_def%gxmin,g_def%gxmax
     .            ,g_def%gymin,g_def%gymax
     .            ,1,110)
      call contour(da_cov(:,:,1,2),nx+2,ny+2
     .            ,g_def%gxmin,g_def%gxmax
     .            ,g_def%gymin,g_def%gymax
     .            ,1,110)
      call contour(da_cov(:,:,1,3),nx+2,ny+2
     .            ,g_def%gxmin,g_def%gxmax
     .            ,g_def%gymin,g_def%gymax
     .            ,1,110)

      ddb = curl(g_def,igx,da_cov)!-db

      call contour(ddb    (:,:,1,1),nx+2,ny+2
     .            ,g_def%gxmin,g_def%gxmax
     .            ,g_def%gymin,g_def%gymax
     .            ,1,110)
      call contour(ddb    (:,:,1,2),nx+2,ny+2
     .            ,g_def%gxmin,g_def%gxmax
     .            ,g_def%gymin,g_def%gymax
     .            ,1,110)
      call contour(ddb    (:,:,1,3),nx+2,ny+2
     .            ,g_def%gxmin,g_def%gxmax
     .            ,g_def%gymin,g_def%gymax
     .            ,1,110)
      close (110)

      write (*,*) "Dumped graphics in B_dump"
cc      call pstop('B_dump','Dumped graphics')

c     End program

      end subroutine B_dump

c     findVectorPotential
c     ################################################################
      subroutine findVectorPotential(g_def,iout,igrid,b_cnv,a_cov,meth
     .                              ,its,rtol,sp,lorder)

c     ----------------------------------------------------------------
c     Finds COVARIANT vector potential from CONTRAVARIANT magnetic
c     field components. There are two methods:
c       1) meth='line': performs line integrals
c       2) meth='vlap': performs vector Laplacian inversion
c     Other input variables include:
c     * sp: whether we are dealing with SP geometry
c     * lorder: order of line integrals
c
c     On input:
c       * b_cnv -> cnv magnetic field
c     On output:
c       * a_cov -> cov vector potential
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: iout,igrid
      real(8) :: a_cov(0:,0:,0:,:)
     .          ,b_cnv(0:,0:,0:,:)
      character(4) :: meth

      integer,optional,INTENT(OUT) :: its
      real(8),optional,INTENT(IN)  :: rtol
      logical,optional :: sp
      integer,optional :: lorder

c     Local variables

      integer :: nxg,nyg,nzg,gm_its
      integer :: i,j,k,ii,nx,ny,nz,ig,jg,kg
     .          ,bcnd(6,3),bcs_b(6,3),bcs_a(6,3),bcsq(6)
      real(8) :: a0_cov(ilom:ihip,jlom:jhip,klom:khip,3)
     .          ,a_cnv (ilom:ihip,jlom:jhip,klom:khip,3)
     .          ,db_cnv(ilom:ihip,jlom:jhip,klom:khip,3)
     .          ,vzeros(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: bx0,bxl,by0,byl,bz0,bzl,x1,y1,z1,r_tol
      logical :: spoint

      real(8),pointer,dimension(:,:,:,:) :: ag,bg
      
c     Begin program

      nx = size(a_cov,1)-2
      ny = size(a_cov,2)-2
      nz = size(a_cov,3)-2

c     Gather global magnetic field

      nxg = g_def%nxgl(igrid)  !Only communicate along X for line integral
      nyg = ny
      nzg = nz

      allocate(ag(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .        ,bg(0:nxg+1,0:nyg+1,0:nzg+1,3))

      call find_global(b_cnv,bg
#if defined(petsc)
     .                ,mpi_comm=g_def%MPI_COMM_X
#endif
     .                )

      if (my_rank == 0.and.(.not.test)) then
        write (*,*)
        call info('  ***********************')
        call info('  Vector potential solver')
        call info('  ***********************')
        write (*,*)
      endif

c     Find A BCs (fluxes) via line integrals

      call curl_inv(g_def,igrid,bg,ag,global=.true.,order=lorder,sp=sp)

      !Restrict solution to current domain (when appropriate)
      call XferGlobalArray(ag,a_cov)

      deallocate(ag,bg)

      !Set topological BCs
      bcnd(:,1) = bcond
      bcnd(:,2) = bcond
      bcnd(:,3) = bcond

      vzeros = 0d0
      a_cnv = XformToCnv(gv%gparams,igrid,a_cov)
      call setBC(gv%gparams,IAX,3,nx,ny,nz,a_cnv,a_cov,vzeros,bcnd
     .          ,igrid,igrid,igrid,is_cnv=.false.,iorder=order_bc)

c diag ****
c$$$      if (test) then
c$$$        call XferGlobalArray(b,db_cnv)
c$$$        call checkCurl  (a_cov,db_cnv)
c$$$      endif
c diag ****

c$$$      !Diagnostics
c$$$      call XferGlobalArray(b,db_cnv)
c$$$      call B_dump(gv%gparams,igrid,a_cov,db_cnv)

c     Correct A covariant

      if (meth == 'vlap') then

        a0_cov = a_cov

        !BC defs (not yet defined in application)
        call default_A_BCs(bcs_a)
        call default_B_BCs(bcs_b)

        !Find B and substract BCs
        db_cnv = b_cnv - curl_mesh(gv%gparams,igrid,a0_cov)

        !Find covariant vector potential update with div(A)=0
        call vlap_curl_inv(gv%gparams,igrid,db_cnv,a_cov
     .                    ,bcs_b,bcs_a,iout,gm_it=gm_its
     .                    ,rtol=rtol,w_pre_bcs=.true.)

        if (.not.PRESENT(its)) then
          if(my_rank==0.and.iout==0) then
            write (*,'(a,i3)') '  A-solve GMRES its=',gm_its
          endif
        else
          its = gm_its
        endif

        !Add BCs
        a_cov = a_cov + a0_cov

c$$$        !Diagnostics
c$$$        call XferGlobalArray(b,db_cnv)
c$$$        call B_dump(gv%gparams,igrid,a_cov,db_cnv)

      endif

c     End program

      contains

c     XferGlobalArray
c     ####################################################################
      subroutine XferGlobalArray(ag,al)

        implicit none
      
        real(8),dimension(0:,0:,0:,:) :: ag,al

        integer :: igm,igp,jgm,jgp,kgm,kgp

        if (nxg == nxd) then
          igm = ilomg ; igp = ihipg
        else
          igm = ilom  ; igp = ihip
        endif

        if (nyg == nyd) then
          jgm = jlomg ; jgp = jhipg
        else
          jgm = jlom  ; jgp = jhip
        endif

        if (nzg == nzd) then
          kgm = klomg ; kgp = khipg
        else
          kgm = klom  ; kgp = khip
        endif

        al(ilom:ihip,jlom:jhip,klom:khip,:)
     .     = ag(igm:igp,jgm:jgp,kgm:kgp,:)

      end subroutine XferGlobalArray

c     checkCurl
c     ####################################################################
      subroutine checkCurl(a_cov,b_cnv)

      implicit none

      integer :: order

      real(8) :: a_cov(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: b_cnv(ilom:ihip,jlom:jhip,klom:khip,3)

      real(8) :: b_cnv2(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: b_cov(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: b_cov2(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: diva (ilom:ihip,jlom:jhip,klom:khip)
      real(8) :: a_cnv(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: j_cnv(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: j_cnv2(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: j_cnv3(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: j_cnv4(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: j_cnv5(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: j_cov (ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: j_cov2(ilom:ihip,jlom:jhip,klom:khip,3)
      real(8) :: j_cov3(ilom:ihip,jlom:jhip,klom:khip,3)

      integer :: j=2,k=1
      
c     Begin program

      a_cnv = XformToCnv(gv%gparams,igrid,a_cov)

      do i=0,nx
        write(*,*) "Acnvdiag",gv%gparams%xx(i+1)
     .                    ,a_cnv (i,j,k,2)  !Benchmark
     .                    ,a_cnv (i,j,k,3)  !Benchmark
        write(*,*) "Acovdiag",gv%gparams%xx(i+1)
     .                    ,a_cov (i,j,k,2)  !Benchmark
     .                    ,a_cov (i,j,k,3)  !curl(A)
      enddo

      !Div(A) diagnostic
      diva = div(gv%gparams,igrid,a_cnv)
      write (*,*) "divA=",sqrt(sum(diva**2))

      !Curl benchmark
      b_cov = XformToCov(gv%gparams,igrid,b_cnv)
      j_cnv = curl(gv%gparams,igrid,b_cov)

      !Curl-veclap diagnostic
      order = 2

c$$$      b_cnv2 = curl(gv%gparams,igrid,a_cov,order=4)
c$$$      b_cov2 = XformToCov(gv%gparams,igrid,b_cnv2)
      b_cnv2 = curl_cell2vrtx_mesh(gv%gparams,igrid,a_cov)
      b_cov2 = XformToCov_vrtx(gv%gparams,igrid,b_cnv2)

      do i=0,nx
        write(*,*) "Bcnvdiag",gv%gparams%xx(i+1)
     .                    ,b_cnv (i,j,k,2)  !Benchmark
     .                    ,b_cnv2(i,j,k,2)  !curl(A)
     .                    ,b_cnv (i,j,k,3)  !Benchmark
     .                    ,b_cnv2(i,j,k,3)  !curl(A)
        write(*,*) "Bcovdiag",gv%gparams%xx(i+1)
     .                    ,b_cov (i,j,k,2)  !Benchmark
     .                    ,b_cov2(i,j,k,2)  !curl(A)
     .                    ,b_cov (i,j,k,3)  !Benchmark
     .                    ,b_cov2(i,j,k,3)  !curl(A)
      enddo

      !Currents
cc      j_cnv2= curl_vrtx2cell_mesh(gv%gparams,igrid,b_cov2)
cc      j_cnv2= curl  (gv%gparams,igrid,b_cov2,order=order)
      j_cnv2=-veclap    (gv%gparams,igrid,a_cnv)
cc      j_cnv3=-veclap_cnv(gv%gparams,igrid,a_cnv)
      j_cnv4=-veclap_cov(gv%gparams,igrid,a_cov)
      j_cnv4= XformToCnv(gv%gparams,igrid,j_cnv4)
      j_cnv5= curlcurl_vrtx(gv%gparams,igrid,a_cov)

      do i=1,nx
        write(*,*) "Jcnvdiag",gv%gparams%xx(i+1)
     .                    ,j_cnv (i,j,k,2)  !Benchmark
     .                    ,j_cnv2(i,j,k,2)  
cc     .                    ,j_cnv3(i,j,k,2)  
     .                    ,j_cnv4(i,j,k,2)  
     .                    ,j_cnv5(i,j,k,2)  
     .                    ,j_cnv (i,j,k,3)  !Benchmark
     .                    ,j_cnv2(i,j,k,3)  
cc     .                    ,j_cnv3(i,j,k,3)  
     .                    ,j_cnv4(i,j,k,3)  
     .                    ,j_cnv5(i,j,k,3)  
      enddo

      j_cov = XformToCov(gv%gparams,igrid,j_cnv)
      j_cov2=-veclap_cov(gv%gparams,igrid,a_cov)
      j_cov3= XformToCov(gv%gparams,igrid,j_cnv3)
      do i=1,nx
        write(*,*) "Jcovdiag",gv%gparams%xx(i+1)
     .                    ,j_cov (i,j,k,2)  !Benchmark
     .                    ,j_cov2(i,j,k,2)  !veclap_cov
     .                    ,j_cov3(i,j,k,2)  !veclap_cnv
     .                    ,j_cov (i,j,k,3)  !Benchmark
     .                    ,j_cov2(i,j,k,3)  !veclap_cov
     .                    ,j_cov3(i,j,k,3)  !veclap_cnv
      enddo

      stop

      end subroutine checkCurl

      end subroutine findVectorPotential

      end module B_tools

c fixDivB
c ######################################################################
      subroutine fixDivB(patch_var)
      use variables
      use auxiliaryVariables 
      implicit none

      type(patch), TARGET :: patch_var
      integer :: nbcgrps
      real(8),pointer :: B_cnv(:,:,:,:)
      real(8),pointer :: B_cov(:,:,:,:)

      INTERFACE
        subroutine solen_bc(g_def,array,dim,loc,face,safe)
          use grid_def_st
          type(grid_mg_def),pointer :: g_def
          integer :: dim,loc
          logical :: face,safe
          real(8) :: array(0:,0:,0:,:)
        end subroutine solen_bc
      END INTERFACE

c Get some basic properties

      gv => patch_var

      call setVectorDimensions

      B_cnv => patch_var%aux%vec_list(IBCNV)%vec
      B_cov => patch_var%aux%vec_list(IBCOV)%vec

      nxbc = gv%gparams%nxv(1)
      nybc = gv%gparams%nyv(1)
      nzbc = gv%gparams%nzv(1)
      iimnm = gv%gparams%ilo(1)
      iimxp = gv%gparams%ihi(1)
      jjmnm = gv%gparams%jlo(1)
      jjmxp = gv%gparams%jhi(1)
      kkmnm = gv%gparams%klo(1)
      kkmxp = gv%gparams%khi(1)
      igxbc = 1
      igybc = 1
      igzbc = 1

c Correct B cnv

      call solen_bc(gv%gparams,B_cnv,1,0,.false.,.false.)
      call solen_bc(gv%gparams,B_cnv,2,0,.false.,.false.)
      call solen_bc(gv%gparams,B_cnv,3,0,.false.,.false.)
      call solen_bc(gv%gparams,B_cnv,1,1,.false.,.false.)
      call solen_bc(gv%gparams,B_cnv,2,1,.false.,.false.)
      call solen_bc(gv%gparams,B_cnv,3,1,.false.,.false.)

c Correct B cov

      call XformVector_BC(gv%gparams,igxbc,B_cnv,B_cov,.false.)

      end subroutine fixDivB
