      MODULE efit

      use io

      use xdraw_io

      use grid

      use local_BCS_variables, ONLY: default_B_BCs,bcond,setMGBC,order_bc

      IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     WRITTEN 10/19/18 BY L. CHACON AS PART OF THE PIXIE3D PROJECT (c)
!     
!     PURPOSE: COMPUTES AND STORES COORDINATES AND FIELDS FROM EFIT DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     VARIABLE DECLARATIONS
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      !EFIT variables
      integer :: nbbbs,limitr,idum,nw,nh,kvtor,nmass,neqdsk
   
      real(8), allocatable,dimension(:,:) :: psirz
      real(8), allocatable,dimension(:)   :: fpol,pres,ffprim,pressw,pwprim,dmion,rhovn &
                                            ,pprime,qpsi,rbbbs,zbbbs,rlim,zlim
      real(8) :: rdim,zdim,rcentr,rleft,zmid
      real(8) :: rmaxis,zmaxis,simag,sibry,bcentr
      real(8) :: current,xdum,rvtor

      !SLATEC spline variables
      integer :: kx,kz,nxs,nzs,dim,flg,sorder
      real(8),dimension(:)  ,allocatable :: tx,tz,tps,work,xs,zs,ps &
                                           ,fpol_coef,pres_coef &
                                           ,ffprim_coef,pprime_coef,qpsi_coef
      real(8),dimension(:,:),allocatable :: psi_coef

      !Module variables
      real(8) :: r_max,r_min,z_max,LL

!!$      REAL(8), DIMENSION(:,:,:), ALLOCATABLE ::                     &
!!$     &             sqrtg,                                               &!sqrt(g): Jacobian on half grid
!!$     &             gss, gsu, gsv, guu, guv, gvv,                        &!symmetric elements of lower metric tensor (full mesh)
!!$     &             hss, hsu, hsv, huu, huv, hvv                          !symmetric elements of upper metric tensor (half mesh)
!!$
!!$      REAL(8), DIMENSION(:,:), ALLOCATABLE ::                       &
!!$     &             rmnc_spline, zmns_spline, bsupumnc_spline,           &
!!$     &             bsupvmnc_spline, bsubvmnc_spline
!!$
!!$      REAL(8), DIMENSION(:,:,:), ALLOCATABLE ::                     &
!!$     &             rmnc_i, zmns_i, bsupumncf_i, bsupvmncf_i, bsubvmncf_i

!!$      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: rr,zz

!!$      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: bsupuijcf,          &
!!$     &  bsupvijcf, bsubvijcf, bsupsijsf, presijf
!!$
!!$!     LOCAL (PRIVATE) HELPER ROUTINES (NOT ACCESSIBLE FROM OUTSIDE THIS MODULE)
!!$!
!!$!LC 2/2/07      PRIVATE spline_fourier_modes, loadrzl_vmec,       &
!!$      PRIVATE spline_fourier_modes,convert_to_full_mesh, vmec_load_rz

      CONTAINS

!     find_RZ
!     ################################################################################
      SUBROUTINE find_RZ(g_def,igrid,i,j,k,RR,ZZ)

        IMPLICIT NONE

        !Call variables
        
        type(grid_mg_def),pointer :: g_def
        integer :: i,j,k,igrid
        real(8) :: RR,ZZ

        !Local variables
        
        integer :: ig,jg,kg
        real(8) :: x1,y1,z1
        
        !Begin program
        
        call getCartesianCoordinates(g_def,i,j,k      &
     &            ,igrid,igrid,igrid,ig,jg,kg,x1,y1,z1)

        RR = sqrt(x1**2+y1**2)
        ZZ = z1

      END SUBROUTINE find_RZ
        
!     read_efit_file
!     ################################################################################
      SUBROUTINE read_efit_file(efit_file,istat)

        IMPLICIT NONE

        !Call variables
        character*(*) :: efit_file
        integer :: istat
        
        !Local variables
        character(10) :: case(6)

        integer :: i,j
        !Begin program
        neqdsk = find_unit(123)
        
        open (unit=neqdsk,file=trim(efit_file),status="old")
   
        read (neqdsk,2000,IOSTAT=istat) (case(i),i=1,6),idum,nw,nh
        if (istat /= 0) return
        
        allocate(psirz(nw,nh),fpol(nw),pres(nw),ffprim(nw),pprime(nw),qpsi(nw) &
                ,pressw(nw),pwprim(nw),dmion(nw),rhovn(nw))
   
        read (neqdsk,2020,IOSTAT=istat) rdim,zdim,rcentr,rleft,zmid
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) rmaxis,zmaxis,simag,sibry,bcentr
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) current,simag,xdum,rmaxis,xdum
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) zmaxis,xdum,sibry,xdum,xdum
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) (fpol(i),i=1,nw)
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) (pres(i),i=1,nw)
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) (ffprim(i),i=1,nw)
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) (pprime(i),i=1,nw)
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) ((psirz(i,j),i=1,nw),j=1,nh)
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) (qpsi(i),i=1,nw)
        if (istat /= 0) return

        read (neqdsk,2022,IOSTAT=istat) nbbbs,limitr
        if (istat /= 0) return

        allocate(rbbbs(nbbbs),zbbbs(nbbbs),rlim(limitr),zlim(limitr))

        read (neqdsk,2020,IOSTAT=istat) (rbbbs(i),zbbbs(i),i=1,nbbbs)
        if (istat /= 0) return
        read (neqdsk,2020,IOSTAT=istat) (rlim(i),zlim(i),i=1,limitr)
        if (istat /= 0) return

        read (neqdsk,2024,IOSTAT=istat) kvtor,rvtor,nmass
        if (istat /= 0) return
        if (kvtor.gt.0) then
          read (neqdsk,2020,IOSTAT=istat) (pressw(i),i=1,nw)
          if (istat /= 0) return
          read (neqdsk,2020,IOSTAT=istat) (pwprim(i),i=1,nw)
          if (istat /= 0) return
        endif
        if (nmass.gt.0) then
          read (neqdsk,2020,IOSTAT=istat) (dmion(i),i=1,nw)
          if (istat /= 0) return
        endif
        read (neqdsk,2020,IOSTAT=istat) (rhovn(i),i=1,nw)
        if (istat /= 0) return

!!        call read_chk()
        close (neqdsk)
   
2000    format (6a8,3i4)
2020    format (5e16.9)
2022    format (2i5)
2024    format (i5,e16.9,i5)

      contains

        subroutine read_chk

          !Plot Psi(R,Z)
          call createDrawInCfile(1,"efit_psi.bin",'Solution','t','x','y' &
                              ,(/'Psi'/),'-c -X0 -L57','drawpsi.in')

          open(unit=110,file="efit_psi.bin",form='unformatted',status='unknown')

          call contour(psirz,nw,nh,rleft,rleft+rdim,zmid-0.5*zdim,zmid+0.5*zdim &
                      ,0,110)
          close(110)

          !Plot flux boundary (gnuplot)
          open(unit=110,file="efit_bdrys.txt",status='unknown')
          do i=1,nbbbs
             write(110,2020) rbbbs(i),zbbbs(i)
          enddo
          close(110)

          !Plot limiter boundary (gnuplot)
          open(unit=110,file="efit_limits.txt",status='unknown')
          do i=1,limitr
             write(110,2020) rlim(i),zlim(i)
          enddo
          close(110)

          write(*,2000) case,idum,nw,nh

          write(*,2020) rdim,zdim,rcentr,rleft,zmid
          write(*,2020) rmaxis,zmaxis,simag,sibry,bcentr
          write(*,2020) current,simag,xdum,rmaxis,xdum
          write(*,2020) zmaxis,xdum,sibry,xdum,xdum
          write(*,2020) (fpol(i),i=1,nw)
          write(*,2020) (pres(i),i=1,nw)
          write(*,2020) (ffprim(i),i=1,nw)
          write(*,2020) (pprime(i),i=1,nw)
          write(*,2020) ((psirz(i,j),i=1,nw),j=1,nh)
          write(*,2020) (qpsi(i),i=1,nw)

          write(*,2022) nbbbs,limitr

          write(*,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
          write(*,2020) (rlim(i),zlim(i),i=1,limitr)
          write(*,2024) kvtor,rvtor,nmass

          if (kvtor.gt.0) then
            write(*,2020) (pressw(i),i=1,nw)
            write(*,2020) (pwprim(i),i=1,nw)
          endif
          if (nmass.gt.0) then
            write(*,2020) (dmion(i),i=1,nw)
          endif
          write(*,2020) (rhovn(i),i=1,nw)

2000      format (6a8,3i4)
2020      format (5e16.9)
2022      format (2i5)
2024      format (i5,e16.9,i5)

          stop "EFIT DIAG"

        end subroutine read_chk
        
      end SUBROUTINE read_efit_file
        
!     efit_init
!     ################################################################################
      SUBROUTINE efit_init(efit_file)

        IMPLICIT NONE

!     -----------------------------------------------------------------------------
!     LOADS VALUES FOR MESHES AND FIELDS TO BE USED
!     -----------------------------------------------------------------------------

!     Call variables

        CHARACTER*(*)       :: efit_file

!     Local variables

        INTEGER :: istat,ig,kg,alloc_stat
        REAL(8) :: dr,dz

        real(8),allocatable,dimension(:) :: q

!     Begin program

!     Initialize variables

!     READ-IN DATA FROM EFIT

        CALL read_efit_file(efit_file, istat)
        IF (istat .ne. 0) STOP 'Read-efit error in efit_init'
        
!     FIND DOMAIN LIMITS AND SETUP GEOMETRY

        z_max = maxval(zbbbs)
        r_max = maxval(rbbbs)
        r_min = minval(rbbbs)

        LL = 0.5*(r_max-r_min) !Dimensionalization length

        write (*,*) "a (m)=",LL
        write (*,*) "R_max (m)=",r_max
        write (*,*) "R_min (m)=",r_min
        write (*,*) "Z_max (m)=",z_max

        r_max = r_max/LL
        r_min = r_min/LL
        z_max = z_max/LL

        write (*,*) "R_max/a=",r_max
        write (*,*) "R_min/a=",r_min
        write (*,*) "Z_max/a=",z_max

        rdim = rdim/LL
        rleft = rleft/LL

        zdim = zdim/LL
        zmid = zmid/LL

!     SPLINE PSI ON R-Z MESH

        sorder = 3
        nxs = nw
        nzs = nh

        flg = 0 !Let spline routine find interpolation knots
        kx = min(sorder+1,nxs-1)
        kz = min(sorder+1,nzs-1)

        dim = nxs*nzs + max(2*kx*(nxs+1),2*kz*(nzs+1))

        allocate(work(dim) ,stat=alloc_stat)
        allocate(tx(nxs+kx),stat=alloc_stat)
        allocate(tz(nzs+kz),stat=alloc_stat)

        !Knots
        allocate(xs(nxs),zs(nzs),stat=alloc_stat)

        dr = rdim/(nxs-1)
        do ig = 1,nxs
          xs(ig) = rleft+dr*(ig-1)
        enddo

        dz = zdim/(nzs-1)
        do kg = 1,nzs
          zs(kg) = zmid-0.5*zdim + dz*(kg-1)
        enddo

        !Coeffs
        allocate(psi_coef(nxs,nzs),stat=alloc_stat)

        call db2ink(xs,nxs,zs,nzs,psirz,nxs,kx,kz,tx,tz,psi_coef,work,flg)

!     SPLINE 1D PSI-DEPENDENT ARRAYS

        allocate(ps(nxs),stat=alloc_stat)

        allocate(fpol_coef(nxs) &
                ,pres_coef(nxs) &
                ,ffprim_coef(nxs)&
                ,pprime_coef(nxs)&
                ,qpsi_coef(nxs),stat=alloc_stat)

        allocate(q((2*kx-1)*nxs))

        !Define poloidal flux domain
        dr = (sibry-simag)/(nxs-1)
        do ig = 1,nxs
          ps(ig) = simag+dr*(ig-1)
        enddo

        allocate(fpol_coef  (nxs) &
                ,pres_coef  (nxs) &
                ,ffprim_coef(nxs) &
                ,pprime_coef(nxs) &
                ,qpsi_coef  (nxs) &
                ,stat=alloc_stat)

        allocate(tps(nxs+kx),stat=alloc_stat)

        call dbknot(ps,nxs,kx,tps)
        call dbintk(ps,fpol  ,tps,nxs,kx,fpol_coef  ,q,work)
        call dbintk(ps,pres  ,tps,nxs,kx,pres_coef  ,q,work)
        call dbintk(ps,ffprim,tps,nxs,kx,ffprim_coef,q,work)
        call dbintk(ps,pprime,tps,nxs,kx,pprime_coef,q,work)
        call dbintk(ps,qpsi  ,tps,nxs,kx,qpsi_coef  ,q,work)

        deallocate(q)
      
!     Destroy storage
        
        deallocate(xs,zs,ps)
         
      END SUBROUTINE efit_init

!     efit_cleanup
!     ################################################################################
      SUBROUTINE efit_cleanup

        IMPLICIT NONE

        INTEGER :: istat

        deallocate(psirz,fpol,pres,ffprim,pprime,qpsi,rbbbs,zbbbs,rlim,zlim &
                  ,pressw,pwprim,dmion,rhovn,STAT=istat)
   
        DEALLOCATE(work,tx,tz,tps,psi_coef,fpol_coef,ffprim_coef,pprime_coef &
                  ,qpsi_coef,stat=istat)

        IF (istat .ne. 0) STOP 'Deallocation error in efit_cleanup'

      END SUBROUTINE efit_cleanup

      END MODULE efit

!     efit_map
!     #################################################################
      subroutine efit_map(equ_file,g_def)

!     -----------------------------------------------------------------
!     Give Cartesian coordinates of each logical mesh point at grid
!     level (igrid).
!     -----------------------------------------------------------------

      use efit

      use equilibrium

      implicit none

!     Input variables

      type(grid_mg_def),pointer :: g_def
      character(*) :: equ_file

!     Local variables

      integer :: nxg,nyg,nzg
      real(8) :: scale

!     Begin program

      if (my_rank == 0) write (*,'(a)') ' Creating EFIT map'

      !Read equilibrium file

      call efit_init(equ_file)

      !Setup geometry
      
      coords = 'tor'

      scale = 1.1
      gparams(1) = 0.5*(r_max+r_min)   !Major radius
      gparams(2) = scale               !Horizontal elliptical radius
      gparams(3) = z_max*scale         !Vertical elliptical radius

      write (*,*) "R/a =",0.5*(r_max+r_min)
      write (*,*) "Z/a =",z_max
      write (*,*) "Domain limits=",xmax,xmin,ymax,ymin,zmax,zmin
      
      nxg = g_def%nglx
      nyg = g_def%ngly
      nzg = g_def%nglz

      call createGrid(nxg,nyg,nzg,xmin,xmax,ymin,ymax,zmin,zmax,g_def)

!     End program

      end subroutine efit_map

!     efit_equ
!     #################################################################
      subroutine efit_equ(iout,igrid,nx,ny,nz,bb,prs,rho,gam,equ_file)

!     -----------------------------------------------------------------
!     Give equilibrium fields at each logical mesh point in grid
!     level (igrid).
!     -----------------------------------------------------------------

        use efit
        
        use equilibrium

        implicit none

!      Call variables

        integer :: igrid,iout,nx,ny,nz
        real(8) :: gam
        real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3) :: bb
        real(8),dimension(0:nx+1,0:ny+1,0:nz+1)   :: prs,rho
        character(*) :: equ_file

!!        logical :: dcon,divcl

!     Local variables

        integer :: i,j,k,igl,jgl,kgl,ig,jg,kg,istat,its,inbv,bcsb(6,3)
        real(8) :: max_prs,RR,ZZ,psi_max,psi_min,BB0,iB0,ip0,ipsi0,a1,a2,x1,y1,z1

        real(8),allocatable, dimension(:,:,:) :: bsub3,psi

        integer :: udcon=1111

        real(8)  :: db2val,dbvalu
        external :: db2val,dbvalu

!     Begin program

        if (my_rank == 0) then
           write (*,*)
           write (*,*) 'Reading EFIT equilibrium...'
        endif

!!$!     Get GLOBAL limits (VMEC operates on global domain)
!!$
!!$        nxg = gv%gparams%nxgl(igrid)
!!$        nyg = gv%gparams%nygl(igrid)
!!$        nzg = gv%gparams%nzgl(igrid)
!!$
!!$!     DCON dump
!!$
!!$        if (dcon .and. my_rank == 0) then
!!$
!!$          write (*,*)
!!$          write (*,*) 'Dumping DCON file with VMEC solution...'
!!$
!!$          allocate(pflx(ns_i),ff_i(ns_i),q_i(ns_i))
!!$
!!$          !VMEC POLOIDAL flux surface positions (integral of 2*pi*jac*B^2)
!!$          ds = 1d0/(ns_i-1)
!!$          ppi = acos(-1d0)
!!$          pflx(1) = 0d0
!!$          do i = 2,ns_i
!!$            pflx(i) = pflx(i-1) + ppi*ds*sqrtg(i,1,1)*(bsupuijcf(i,1,1)+bsupuijcf(i-1,1,1))
!!$          enddo
!!$
!!$          !F factor (R*Bt=B_3)
!!$          ff_i = bsubvijcf(:,1,1)
!!$
!!$          !q-profile
!!$          q_i = 1d0/iotaf_i
!!$
!!$          !DCON dump
!!$          open(unit=udcon,file='pixie-dcon.bin',form='unformatted',status='unknown')
!!$
!!$          !Poloidal plane size
!!$          write (udcon) nxg,nyg
!!$
!!$          write (udcon) pflx      !Poloidal flux
!!$          write (udcon) ff_i      !R*B_t (flux function)
!!$          write (udcon) presf_i   !Pressure
!!$          write (udcon) q_i       !Q-profile
!!$
!!$          !R(psi,theta), Z(psi,theta)
!!$
!!$          k = 1  !Fix poloidal plane
!!$          do j=1,nyg/2+1   !Cycle in poloidal angle (only half-plane)
!!$            write (udcon) rr(:,j,k)
!!$            write (udcon) zz(:,j,k)
!!$          enddo
!!$
!!$          close (udcon)
!!$
!!$          deallocate(pflx,ff_i,q_i)
!!$        endif

        a1 = gv%gparams%params(2)  !Minor radius
        a2 = gv%gparams%params(3)  !Major radius
        
!     Interpolate poloidal flux and associated qtys

        allocate(psi  (0:nx+1,0:ny+1,0:nz+1) &
                ,bsub3(0:nx+1,0:ny+1,0:nz+1))

        inbv = 1

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1

              call find_RZ(gv%gparams,igrid,i,j,k,RR,ZZ)
!!$              call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)
!!$
!!$              !Find toroidal coordinates
!!$              r1  = gv%gparams%xx(ig)
!!$              th1 = gv%gparams%yy(jg)
!!$              ph1 = gv%gparams%zz(kg)
!!$
!!$              !Find R,Z coordinates
!!$              RR = gv%gparams%params(1) + a1*r1*cos(th1)
!!$              ZZ = a2*r1*sin(th1)

              !Interpolate poloidal flux
              psi(i,j,k) = db2val(RR,ZZ,0,0,tx,tz,nxs,nzs  &
                                 ,kx,kz,psi_coef,work)

              !Interpolate pressure
              if (psi(i,j,k) > sibry) then
                prs(i,j,k)   = dbvalu(tps,pres_coef,nxs,kx,0,sibry,inbv,work)
                bsub3(i,j,k) = dbvalu(tps,fpol_coef,nxs,kx,0,sibry,inbv,work)
              else
                prs(i,j,k)   = dbvalu(tps,pres_coef,nxs,kx,0,psi(i,j,k),inbv,work)
                bsub3(i,j,k) = dbvalu(tps,fpol_coef,nxs,kx,0,psi(i,j,k),inbv,work)
              endif

              if (prs(i,j,k) < 0d0) prs(i,j,k) = 0d0

            enddo
          enddo
        enddo

!     Check poloidal flux limits
        
        open(unit=110,file="efit_psi.txt",status='unknown')

        do j=0,ny+1
           do i=0,nx+1

              call find_RZ(gv%gparams,igrid,i,j,1,RR,ZZ)
!!$              call getMGmap(gv%gparams,i,j,1,igrid,igrid,igrid,ig,jg,kg)
!!$
!!$              !Find toroidal coordinates
!!$              r1  = gv%gparams%xx(ig)
!!$              th1 = gv%gparams%yy(jg)
!!$              ph1 = gv%gparams%zz(kg)
!!$
!!$              !Find R,Z coordinates
!!$              RR = gv%gparams%params(1) + a1*r1*cos(th1)
!!$              ZZ = a2*r1*sin(th1)

              write (110,*) RR,ZZ,psi(i,j,1)
           enddo
           write (110,*)
        enddo

        close(110)

!!$        !Plot Psi(R,Z)
!!$        call createDrawInCfile(3,"efit_p3d.bin",'Solution','t','x','y' &
!!$                            ,(/'Psi','Prs','B_3'/),'-c -X0 -L57','drawp3d.in')
!!$
!!$        open(unit=110,file="efit_p3d.bin",form='unformatted',status='unknown')
!!$
!!$        call contour(psi  (:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,0,110)
!!$        call contour(prs  (:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
!!$        call contour(bsub3(:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
!!$        close(110)

!!$        write (*,*) "Magnetic axis coords=",rmaxis/LL,zmaxis/LL
!!$        psi_min = db2val(rmaxis/LL,zmaxis/LL,0,0,tx,tz,nxs,nzs  &
!!$     &                  ,kx,kz,psi_coef,work)
!!$
!!$        psi_max = db2val(r_max,0d0,0,0,tx,tz,nxs,nzs  &
!!$     &                  ,kx,kz,psi_coef,work)
!!$
!!$        write (*,*) "dpsi_min",psi_min-simag
!!$        write (*,*) "dpsi_max =",psi_max-sibry
!!$
!!$        STOP "efit_equ"

!!$        psi_min = simag
!!$        psi_max = sibry

!     Normalization magnetic field (toroidal at magnetic axis) and pressure

        BB0 = dbvalu(tps,fpol_coef,nxs,kx,0,simag,inbv,work)/rmaxis

        iB0 = 1d0/BB0
        
        ip0 = 2*iB0*iB0*(4*pi)*1d-7

        ipsi0 = iB0/LL

!     Normalized pressure

        prs = prs*ip0

!     Find magnetic field components

        !Poloidal B (cnv)
        bb = 0d0
        bb(:,:,:,3) = psi*ipsi0  !Cov A

        bb = curl(gv%gparams,igrid,bb)

        !Toroidal B
        bb(:,:,:,3) = bsub3*iB0  !Covariant

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              bb(i,j,k,3) = (bb(i,j,k,3)                                        &
     &           -(gv%gparams%gmetric%grid(igrid)%gsub(i,j,k,3,1)*bb(i,j,k,1)   &
     &            +gv%gparams%gmetric%grid(igrid)%gsub(i,j,k,3,2)*bb(i,j,k,2))) &
     &            /gv%gparams%gmetric%grid(igrid)%gsub(i,j,k,3,3)       !Xform to cnv

            enddo
          enddo
        enddo

        !BCs
        bcsb(:,1) = bcond
        bcsb(:,2) = bcond
        bcsb(:,3) = bcond
        call default_B_BCs(bcsb)
        where (bcsb == -NEU) bcsb = -EXT  !Extrapolate tangential components
        
        call setMGBC(gv%gparams,0,3,nx,ny,nz,igrid,bb,bcsb &
     &              ,icomp=(/IBX/),is_vec=.true.           &
     &              ,is_cnv=.true.,iorder=order_bc)

!     Find normalized density

        max_prs = maxval(prs)
        max_prs = pmax(max_prs)
        
        if (max_prs == 0d0) then
          rho = 1d0
        else
          rho = (abs(prs/max_prs)**(1d0/gam))  !Forces positive density
        endif

!!$        !Diagnostics plots
!!$        call createDrawInCfile(7,"efit_p3d.bin",'Solution','t','x','y' &
!!$                            ,(/'Psi','Prs','B_3','B^1','B^2','B^3','rho'/) &
!!$                            ,'-c -X0 -L57','drawp3d.in')
!!$
!!$        open(unit=110,file="efit_p3d.bin",form='unformatted',status='unknown')
!!$
!!$        call contour(psi  (:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,0,110)
!!$        call contour(prs  (:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
!!$        call contour(bsub3(:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
!!$        call contour(bb   (:,:,1,1),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
!!$        call contour(bb   (:,:,1,2),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
!!$        call contour(bb   (:,:,1,3),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
!!$        call contour(rho  (:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
!!$
!!$        close(110)

!     Free work space

        DEALLOCATE(psi,bsub3,stat=istat)
        
        call efit_cleanup

!     End program

      end subroutine efit_equ

!!$!     vmec_init_fields
!!$!     #################################################################
!!$      SUBROUTINE vmec_init_fields
!!$
!!$!     -----------------------------------------------------------------
!!$!     Gives contravariant magnetic fields components and pressure
!!$!     at grid level igrid
!!$!     -----------------------------------------------------------------
!!$
!!$       USE stel_kinds
!!$       USE stel_constants            
!!$       USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,                &
!!$    &          mpol=>mpol_i, ntor=>ntor_i, nuv=>nuv_i, mnmax=>mnmax_i,        &
!!$    &          ohs=>ohs_i, nfp=>nfp_i
!!$       USE vmec_mod, ONLY: bsupumncf_i,bsupvmncf_i,bsubvmncf_i,presf_i        &
!!$                          ,bsupuijcf  ,bsupvijcf  ,bsubvijcf  ,bsupsijsf      &
!!$                          ,presijf
!!$       USE fourier, ONLY:  toijsp
!!$
!!$       IMPLICIT NONE
!!$
!!$!      Local variables
!!$
!!$       INTEGER:: istat, jk
!!$
!!$!      Begin program
!!$
!!$!      Allocate variables
!!$
!!$       ALLOCATE (bsupsijsf(ns,ntheta,nzeta)  &
!!$                ,bsupuijcf(ns,ntheta,nzeta)  &
!!$                ,bsupvijcf(ns,ntheta,nzeta)  &
!!$                ,bsubvijcf(ns,ntheta,nzeta), stat=istat)
!!$       IF (istat .ne. 0) STOP 'Allocation error in vmec_init_fields'
!!$
!!$       ALLOCATE(presijf(ns,ntheta,nzeta), stat=istat)
!!$       IF (istat .ne. 0) STOP 'Allocation error in vmec_init_fields'
!!$
!!$!      Fill GLOBAL variables (at VMEC integer mesh -- PIXIE3D's half mesh)
!!$
!!$!LC12/14/06       bsupsijsh = zero                       ! It is zero in VMEC (flux surfaces)
!!$       bsupsijsf = zero                                ! It is zero in VMEC (flux surfaces)
!!$       CALL toijsp (bsupumncf_i, bsupuijcf, 0, 0, 0, 0)
!!$       CALL toijsp (bsupvmncf_i, bsupvijcf, 0, 0, 0, 0)
!!$       CALL toijsp (bsubvmncf_i, bsubvijcf, 0, 0, 0, 0)
!!$
!!$       DO jk = 1, ns
!!$         presijf(jk,:,:) =  presf_i(jk)                ! Init pressure
!!$       ENDDO
!!$
!!$      END SUBROUTINE vmec_init_fields 
