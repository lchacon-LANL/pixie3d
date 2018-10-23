      MODULE efit_mod

      use io

      use xdraw_io

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

      integer :: nxe,nye,nze

      !EFIT variables
      integer :: nbbbs,limitr,idum,nw,nh,kvtor,nmass,neqdsk
   
      real(8), allocatable,dimension(:,:) :: psirz
      real(8), allocatable,dimension(:)   :: fpol,pres,ffprim,pressw,pwprim,dmion,rhovn &
                                            ,pprime,qpsi,rbbbs,zbbbs,rlim,zlim
      real(8) :: rdim,zdim,rcentr,rleft,zmid
      real(8) :: rmaxis,zmaxis,simag,sibry,bcentr
      real(8) :: current,xdum,rvtor

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

      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: rr,zz

!!$      REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: bsupuijcf,          &
!!$     &  bsupvijcf, bsubvijcf, bsupsijsf, presijf
!!$
!!$!     LOCAL (PRIVATE) HELPER ROUTINES (NOT ACCESSIBLE FROM OUTSIDE THIS MODULE)
!!$!
!!$!LC 2/2/07      PRIVATE spline_fourier_modes, loadrzl_vmec,       &
!!$      PRIVATE spline_fourier_modes,convert_to_full_mesh, vmec_load_rz

      CONTAINS

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

        call read_chk()

        close (neqdsk)
   
2000    format (6a8,3i4)
2020    format (5e16.9)
2022    format (2i5)
2024    format (i5,e16.9,i5)

      contains

        subroutine read_chk

          !Plot Psi(R,Z)
          call createDrawInCfile(1,"efit_psi.bin",'Solution','t','x','y' &
                              ,(/'Psi'/),'-c -X0 -L57','drawefit.in')

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


!!$          write(*,2000) case,idum,nw,nh
!!$
!!$          write(*,2020) rdim,zdim,rcentr,rleft,zmid
!!$          write(*,2020) rmaxis,zmaxis,simag,sibry,bcentr
!!$          write(*,2020) current,simag,xdum,rmaxis,xdum
!!$          write(*,2020) zmaxis,xdum,sibry,xdum,xdum
!!$          write(*,2020) (fpol(i),i=1,nw)
!!$          write(*,2020) (pres(i),i=1,nw)
!!$          write(*,2020) (ffprim(i),i=1,nw)
!!$          write(*,2020) (pprime(i),i=1,nw)
!!$          write(*,2020) ((psirz(i,j),i=1,nw),j=1,nh)
!!$          write(*,2020) (qpsi(i),i=1,nw)
!!$
!!$          write(*,2022) nbbbs,limitr
!!$
!!$          write(*,2020) (rbbbs(i),zbbbs(i),i=1,nbbbs)
!!$          write(*,2020) (rlim(i),zlim(i),i=1,limitr)
!!$          write(*,2024) kvtor,rvtor,nmass
!!$
!!$          if (kvtor.gt.0) then
!!$            write(*,2020) (pressw(i),i=1,nw)
!!$            write(*,2020) (pwprim(i),i=1,nw)
!!$          endif
!!$          if (nmass.gt.0) then
!!$            write(*,2020) (dmion(i),i=1,nw)
!!$          endif
!!$          write(*,2020) (rhovn(i),i=1,nw)

2000      format (6a8,3i4)
2020      format (5e16.9)
2022      format (2i5)
2024      format (i5,e16.9,i5)

          stop "EFIT DIAG"

        end subroutine read_chk
        
      end SUBROUTINE read_efit_file
        
!     efit_init
!     ################################################################################
      SUBROUTINE efit_init(nx,ny,nz,efit_file)

        IMPLICIT NONE

!     -----------------------------------------------------------------------------
!     LOADS VALUES FOR MESHES AND FIELDS TO BE USED
!     -----------------------------------------------------------------------------

!     Call variables

        CHARACTER*(*)       :: efit_file
        INTEGER, INTENT(in) :: nx,ny,nz

!     Local variables

        INTEGER :: istat, ntype, imesh, js, ig, kg
        REAL(8) :: t1, t2, dr, dz
        INTEGER :: m, n, mn
        LOGICAL :: load_metrics

        !SLATEC spline variables
        integer :: kx,kz,nxs,nzs,dim,flg,sorder,alloc_stat
        real(8),dimension(:)  ,allocatable :: tx,tz,work,xs,zs
        real(8),dimension(:,:),allocatable :: psi_coef

!     Begin program

!     Initialize variables

        nxe = nx
        nye = ny
        nze = nz

!     READ-IN DATA FROM EFIT

        CALL read_efit_file(efit_file, istat)
        IF (istat .ne. 0) STOP 'Read-efit error in efit_init'

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

        call db2ink(xs,nxs,zs,nzs,psirz,nxs,nzs,kx,kz,tx,tz,psi_coef,work,flg)

!     INVERT PSI MAP: R(psi,theta), Z(psi,theta)
      
!     SPLINE R AND Z DATA
        

!!$!       Spline 1-D arrays: careful -> convert phipf VMEC and multiply
!!$!       by ds-vmec/ds-island, since phipf_i = d(PHI)/ds-island
!!$
!!$        CALL Spline_OneD_Array (iotaf_vmec, iotaf_i, istat)
!!$        CALL Spline_OneD_Array (phipf_vmec, phipf_i, istat)
!!$        presf_vmec = mu0 * presf_vmec
!!$        CALL Spline_OneD_Array (presf_vmec, presf_i, istat)
!!$        DO js = 1, ns_i
!!$          phipf_i(js) = 2 * hs_i*(js-1) * phipf_i(js) / (2*pi)
!!$        END DO
!!$          
!!$!     CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON "POLAR" MESH
!!$!     AND COMPUTE METRIC ELEMENTS AND JACOBIAN
!!$
!!$        CALL efit_load_rz(istat)
!!$        IF (istat .ne. 0) STOP 'Load R,Z error in EFIT_INIT'

        !Destroy storage
        deallocate(psi_coef,xs,zs)
        call efit_cleanup
        
      END SUBROUTINE efit_init

!!$!       efit_load_rz
!!$!       ################################################################################
!!$        SUBROUTINE efit_load_rz(istat)
!!$
!!$!       ------------------------------------------------------------------------------
!!$!       Computes R-Z map in real space from VMEC equilibrium
!!$!
!!$!       Modified by L. Chacon 11/30/06 from routines by S. Hirschman and R. Sanchez
!!$!       ------------------------------------------------------------------------------
!!$
!!$          USE Fourier, ONLY: toijsp
!!$!          USE dump_output
!!$
!!$!       Call variables
!!$
!!$          INTEGER, INTENT(out)     :: istat
!!$          LOGICAL :: load_metrics
!!$
!!$!       Local variables
!!$
!!$          REAL(8), DIMENSION(ns_vmec)  :: fac1, fac2
!!$          INTEGER                          :: mpol_save, ntor_save,iparity,ipol
!!$
!!$!       Begin program
!!$
!!$!       LOAD FOURIER FACTORS USING VMEC mpol,ntor AND ISLAND nu,nv
!!$!       IDEA: USE VMEC mpol, ntor VALUES AND CALL fixarray
!!$      
!!$          CALL fixarray
!!$
!!$!       FIRST, REPACK SPLINED ARRAYS FROM VMEC-ORDERING TO ISLAND-ORDERING
!!$!       AND SWAP THE N->-N MODES TO BE CONSISTENT WITH mu+nv ISLAND ARGUMENT
!!$
!!$          CALL repack (istat)
!!$          IF (istat .ne. 0) STOP 'REPACK error in VMEC_LOAD_RZ'
!!$
!!$!       COMPUTE AND STORE R, Z
!!$
!!$          iparity = 0
!!$          CALL toijsp (rmnc_i,rr(:,:,1:nv_i),0,0,iparity,0)
!!$      
!!$          iparity = 1
!!$          CALL toijsp (zmns_i,zz(:,:,1:nv_i),0,0,iparity,0)
!!$
!!$!       ENFORCE PERIODIC BCs ALONG TOROIDAL DIRECTION
!!$
!!$          rr(:,:,0     ) = rr(:,:,nv_i)
!!$          rr(:,:,nv_i+1) = rr(:,:,1   )
!!$
!!$          zz(:,:,0     ) = zz(:,:,nv_i)
!!$          zz(:,:,nv_i+1) = zz(:,:,1   )
!!$
!!$!       DEBUG: CHECK IF CORRECT
!!$!
!!$!        istat = (ns_i-1)/2 + 1    !rho = 1/2 => s = 1/4
!!$!        CALL dump_special(r1_i,z1_i,ru_i,zu_i,rv_i,zv_i,istat)
!!$
!!$!       LOAD METRIC ELEMENTS
!!$
!!$          if (load_metrics) then
!!$            call vmec_load_metrics(istat)
!!$            IF (istat .ne. 0) STOP 'Error in VMEC_LOAD_METRICS'
!!$          endif
!!$
!!$!       CLEAN-UP EXTRA ARRAYS
!!$
!!$          DEALLOCATE (rmnc_i, zmns_i, stat=istat)
!!$          IF (lwout_opened) CALL read_wout_deallocate
!!$
!!$        END SUBROUTINE efit_load_rz

!       efit_cleanup
!       ################################################################################
        SUBROUTINE efit_cleanup

          IMPLICIT NONE

          INTEGER :: istat

          deallocate(psirz,fpol,pres,ffprim,pprime,qpsi,rbbbs,zbbbs,rlim,zlim &
                    ,pressw,pwprim,dmion,rhovn,STAT=istat)
   
        END SUBROUTINE efit_cleanup

!!$!       init_metric_elements
!!$!       ################################################################################
!!$        SUBROUTINE init_metric_elements(ns_in,mpol_in,ntor_in,wout_file)
!!$
!!$          IMPLICIT NONE
!!$!-----------------------------------------------
!!$!   D u m m y   A r g u m e n t s
!!$!-----------------------------------------------
!!$          CHARACTER*(*)            :: wout_file
!!$          INTEGER, INTENT(in)      :: ns_in, mpol_in, ntor_in
!!$          INTEGER                  :: istat, ntype, imesh, js
!!$          REAL(8), DIMENSION(:,:), ALLOCATABLE ::                       &
!!$     &                            bsupumnc_v,  bsupvmnc_v
!!$          REAL(8) :: t1, t2
!!$          INTEGER     :: m, n, mn
!!$!-----------------------------------------------
!!$!
!!$!     LOADS VALUES FOR MESHES TO BE USED IN VMECPP ISLAND SOLVER
!!$!     (GENERALLY DIFFERENT FROM VMEC MESHES. USE SUBSCRIPT _i FOR ISLAND VARIABLES)
!!$!     READS wout_file FROM VMEC TO GET FOURIER COMPONENTS ON VMEC MESH
!!$!     SPLINES THE VMEC COMPONENTS TO RADIAL ISLAND MESH (PHI->SQRT(PHI))
!!$!     CALLS Fourier MODULE fixarray TO LOAD TRIG ARRAYS FOR COMPUTING ISLAND METRICS
!!$!     COMPUTES gij (sub/sup) AND JACOBIAN ON ISLAND MESHES IN REAL SPACE
!!$!
!!$          ns_i = ns_in
!!$          nsh  = ns_i-1
!!$          ntor_i = ntor_in
!!$          mpol_i = mpol_in
!!$
!!$! Set number of points == number of modes for now! (May want mid-points for flux conservation)
!!$!SPH : # POINTS???
!!$          nu_i = mpol_i + 3
!!$          nv_i = 2*ntor_i + 2
!!$          IF (ntor_i .eq. 0) nv_i = 1
!!$          nuv_i = nu_i*nv_i
!!$          mnmax_i = (mpol_i + 1)*(2*ntor_i + 1)             ! Added RS. Contains total number of modes.
!!$
!!$!
!!$!     READ-IN DATA FROM VMEC-PRODUCED WOUT FILE (LIBSTELL ROUTINE)
!!$!
!!$          CALL read_wout_file(wout_file, istat)
!!$          IF (istat .ne. 0) STOP 'Read-wout error in INIT_METRIC_ELEMENTS'
!!$
!!$          nfp_i = nfp_vmec
!!$          wb_i  = (4*pi*pi)*wb_vmec
!!$
!!$          IF (wb_i .eq. 0._dp) STOP 'wb_vmec = 0!'
!!$
!!$!
!!$!     Allocate space for splined arrays
!!$!
!!$          ALLOCATE(rmnc_spline(mnmax,ns_i), zmns_spline(mnmax,ns_i),        &
!!$     &             bsupumnc_v(mnmax,ns_vmec),  bsupvmnc_v(mnmax,ns_vmec),   &
!!$     &             bsupumnc_spline(mnmax,ns_i), bsupvmnc_spline(mnmax,ns_i),&
!!$     &             phipf_i(ns_i), iotaf_i(ns_i),presf_i(ns_i),              &
!!$     &             stat=istat)
!!$          IF (istat .ne. 0) STOP 'Allocation error 1 in INIT_METRIC_ELEMENTS'
!!$
!!$!
!!$!     SPLINE R, Z. L FOURIER COMPONENTS in s FROM ORIGINAL VMEC MESH (s ~ phi, ns_vmec points) 
!!$!     TO A "POLAR" MESH [s ~ sqrt(phi), ns_i POINTS] WITH BETTER AXIS RESOLUTION
!!$!
!!$          DO ntype = 1, 3
!!$             IF (ntype .eq. 1) THEN
!!$                istat = 0
!!$                CALL Spline_Fourier_Modes(rmnc_vmec, rmnc_spline, istat)   
!!$             ELSE IF (ntype .eq. 2) THEN
!!$                istat = 0
!!$                CALL Spline_Fourier_Modes(zmns_vmec, zmns_spline, istat)   
!!$             ELSE 
!!$                CALL Convert_From_Nyq(bsupumnc_v,bsupumnc_vmec)       !These have mnmax_nyq members
!!$                CALL Convert_To_Full_Mesh(bsupumnc_v)
!!$                istat = 1
!!$                CALL Spline_Fourier_Modes(bsupumnc_v, bsupumnc_spline, istat)   
!!$
!!$                CALL Convert_From_Nyq(bsupvmnc_v,bsupvmnc_vmec)       !These have mnmax_nyq members
!!$                CALL Convert_To_Full_Mesh(bsupvmnc_v)
!!$                istat = 1
!!$                CALL Spline_Fourier_Modes(bsupvmnc_v, bsupvmnc_spline, istat)   
!!$                DEALLOCATE(bsupumnc_v, bsupvmnc_v)
!!$             END IF
!!$
!!$             IF (istat .ne. 0) STOP 'Spline error in INIT_METRIC_ELEMENTS'
!!$          END DO
!!$   
!!$!
!!$!     Spline 1-D arrays: careful -> convert phipf VMEC and multiply
!!$!     by ds-vmec/ds-island, since phipf_i = d(PHI)/ds-island
!!$!
!!$          CALL Spline_OneD_Array (iotaf_vmec, iotaf_i, istat)
!!$          CALL Spline_OneD_Array (phipf_vmec, phipf_i, istat)
!!$          presf_vmec = mu0 * presf_vmec
!!$          CALL Spline_OneD_Array (presf_vmec, presf_i, istat)
!!$          DO js = 1, ns_i
!!$             phipf_i(js) = 2 * hs_i*(js-1) * phipf_i(js) / (2*pi)
!!$          END DO
!!$          
!!$!
!!$!     CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON "POLAR" MESH
!!$!     AND COMPUTE METRIC ELEMENTS AND JACOBIAN
!!$!
!!$          CALL LoadRZL_VMEC(istat)
!!$          IF (istat .ne. 0) STOP 'LoadRZL error in INIT_METRIC_ELEMENTS'
!!$
!!$        END SUBROUTINE init_metric_elements

!!$      SUBROUTINE Spline_Fourier_Modes(ymn_vmec, ymn_spline, istat)
!!$
!!$      use oned_int   !Added by L. Chacon 6/5/07
!!$
!!$!-----------------------------------------------
!!$!   D u m m y   A r g u m e n t s
!!$!-----------------------------------------------
!!$      INTEGER, INTENT(inout)     :: istat
!!$      REAL(8), DIMENSION(mnmax, ns_vmec), TARGET,                   &
!!$     &                     INTENT(in)  :: ymn_vmec
!!$      REAL(8), DIMENSION(mnmax, ns_i), TARGET,                      &
!!$     &                     INTENT(out) :: ymn_spline
!!$!-----------------------------------------------
!!$!   L o c a l   V a r i a b l e s
!!$!-----------------------------------------------
!!$      REAL(8), PARAMETER          :: one = 1
!!$      INTEGER                         :: js, modes, ntype, mp
!!$      REAL(8), DIMENSION(ns_vmec) :: y_vmec!, y_spline
!!$      REAL(8), DIMENSION(ns_vmec) :: snodes_vmec, y2_vmec, fac1
!!$      REAL(8), DIMENSION(ns_i)    :: snodes, fac2
!!$      REAL(8)                     :: hs_vmec, yp1, ypn, expm
!!$!-----------------------------------------------
!!$!
!!$!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM VMEC (s~phi) TO s~sqrt(phi) MESH
!!$!
!!$
!!$!
!!$!     1. Factors for taking out (putting back) [sqrt(s)]**m factor
!!$
!!$      IF (ns_vmec .le. 1) THEN
!!$         istat = 1
!!$         RETURN
!!$      END IF
!!$
!!$      fac1 = 0
!!$      hs_vmec = one/(ns_vmec-1)
!!$      DO js = 2, ns_vmec
!!$         fac1(js) = one/SQRT(hs_vmec*(js-1))
!!$      END DO
!!$
!!$      IF (ns_i .le. 1) THEN
!!$         istat = 2
!!$         RETURN
!!$      END IF
!!$
!!$      ohs_i = ns_i-1
!!$      hs_i = one/ohs_i
!!$      DO js = 1, ns_i
!!$         fac2(js) = hs_i*(js-1)
!!$      END DO
!!$
!!$!
!!$!     2. Set up s-nodes in original (vmec) radial spatial coordinates for the final (splined) mesh 
!!$!        using the mapping relation: s_vmec(sfinal) = sfinal**2 since sfinal = sqrt(phi)
!!$!
!!$
!!$      DO js = 1, ns_i
!!$         snodes(js) = fac2(js)**2
!!$      END DO
!!$
!!$!
!!$!     3. Set up "knots" on initial (vmec, svmec~phi) mesh 
!!$!
!!$      DO js = 1, ns_vmec
!!$         snodes_vmec(js) = hs_vmec*((js-1))
!!$      END DO
!!$
!!$      ymn_spline = 0d0   !L. Chacon, 2/20/07, to avoid definition error below
!!$
!!$      DO modes = 1, mnmax
!!$         y_vmec   = ymn_vmec(modes,:)
!!$!         y_spline = ymn_spline(modes,:)
!!$         mp = xm_vmec(modes)
!!$
!!$         IF (istat.eq.0 .and. mp.gt.0) THEN
!!$            IF (MOD(mp,2) .eq. 1) THEN 
!!$               expm = 1d0
!!$            ELSE
!!$               expm = 2d0
!!$            END IF
!!$            y_vmec = y_vmec*(fac1**expm)
!!$!LC 6/5/07            IF (mp .le. 2) y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
!!$            IF (mp .le. 2) then
!!$               call IntDriver1d(4,snodes_vmec(2:5),y_vmec(2:5)   &
!!$     &                         ,1,snodes_vmec(1:1),y_vmec(1:1),3)
!!$            ENDIF
!!$         END IF
!!$
!!$!
!!$!        4 and 5: spline vmec coefficients into snodes mesh (L. Chacon, 6/5/07)
!!$!
!!$         call IntDriver1d(ns_vmec,snodes_vmec,y_vmec   &
!!$     &                   ,ns_i   ,snodes     ,ymn_spline(modes,:),3)
!!$
!!$!
!!$!        RECOVER RADIAL DEPENDENCE
!!$!
!!$         IF (istat.eq.0 .and. mp.gt.0) THEN
!!$            ymn_spline(modes,:) = ymn_spline(modes,:)*(fac2**expm)
!!$         END IF
!!$
!!$      END DO
!!$
!!$!LC      stop
!!$      
!!$      istat = 0
!!$
!!$      END SUBROUTINE Spline_Fourier_Modes
!!$
!!$      SUBROUTINE Spline_OneD_Array (ns_vmec,y_vmec,ns_i,y_spline,istat)
!!$      IMPLICIT NONE
!!$!-----------------------------------------------
!!$!   D u m m y   A r g u m e n t s
!!$!-----------------------------------------------
!!$      INTEGER, INTENT(out) :: istat
!!$      INTEGER, INTENT(in)  :: ns_vmec,ns_i  
!!$      REAL(8), DIMENSION(ns_vmec), INTENT(in)  :: y_vmec
!!$      REAL(8), DIMENSION(ns_i)   , INTENT(out) :: y_spline
!!$!-----------------------------------------------
!!$!   L o c a l   V a r i a b l e s
!!$!-----------------------------------------------
!!$      REAL(8), PARAMETER          :: one = 1
!!$      INTEGER                         :: js, modes, ntype, mp
!!$      REAL(8), DIMENSION(ns_vmec) :: snodes_vmec, y2_vmec, fac1
!!$      REAL(8), DIMENSION(ns_i)    :: snodes, fac2
!!$      REAL(8)                     :: hs_vmec, yp1, ypn
!!$!-----------------------------------------------
!!$!
!!$!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM s~PSI TO s~sqrt(phi) MESH
!!$!
!!$
!!$!
!!$!     1. Factors for taking out (putting back) sqrt(s) factor for m-odd modes
!!$
!!$      IF (ns_vmec .le. 1) THEN
!!$         istat = 1
!!$         RETURN
!!$      END IF
!!$
!!$      fac1 = 0
!!$      hs_vmec = one/(ns_vmec-1)
!!$      DO js = 2, ns_vmec
!!$         fac1(js) = one/SQRT(hs_vmec*(js-1))
!!$      END DO
!!$
!!$      IF (ns_i .le. 1) THEN
!!$         istat = 2
!!$         RETURN
!!$      END IF
!!$
!!$      ohs_i = ns_i-1
!!$      hs_i = one/ohs_i
!!$      DO js = 1, ns_i
!!$         fac2(js) = hs_i*(js-1)
!!$      END DO
!!$
!!$!
!!$!     2. Set up s-nodes on final (splined) mesh [sfinal ~ sqrt(s_vmec)]
!!$!
!!$
!!$      DO js = 1, ns_i
!!$         snodes(js) = fac2(js)*fac2(js)
!!$      END DO
!!$
!!$!
!!$!     3. Set up "knots" on initial (vmec, svmec~phi) mesh 
!!$!
!!$      DO js = 1, ns_vmec
!!$         snodes_vmec(js) = hs_vmec*((js-1))
!!$      END DO
!!$
!!$
!!$!     4. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!!$!
!!$      yp1 = -1.e30_dp;  ypn = -1.e30_dp
!!$      CALL spline (snodes_vmec, y_vmec, ns_vmec, yp1, ypn, y2_vmec)
!!$
!!$!
!!$!     5. Interpolate onto snodes mesh
!!$!
!!$      DO js = 1, ns_i
!!$         CALL splint (snodes_vmec, y_vmec, y2_vmec, ns_vmec,         &
!!$     &                   snodes(js), y_spline(js))
!!$      END DO
!!$
!!$      END SUBROUTINE Spline_OneD_Array

!!$!       half_to_int
!!$!       ###########################################################################
!!$        subroutine half_to_int(igl,jgl,kgl,vmec_arr,local_val,order)
!!$
!!$!       ---------------------------------------------------------------------------
!!$!       Averages quantities from VMEC's full radial mesh (which is
!!$!       PIXIE's half radial mesh) to PIXIE's collocated mesh. At ghost
!!$!       cells, it extrapolates using a first-order formula.
!!$!       ---------------------------------------------------------------------------
!!$
!!$          use oned_int
!!$
!!$!       Call variables
!!$
!!$          integer :: igl,jgl,kgl
!!$          integer,optional :: order
!!$
!!$          real(8) :: vmec_arr(ns_i,nu_i,nv_i),local_val
!!$
!!$!       Local variables
!!$
!!$          integer :: ordr
!!$          real(8) :: pos(5),val(5)
!!$
!!$!       Begin program
!!$
!!$          if (PRESENT(order)) then
!!$             ordr = order
!!$          else
!!$             ordr = 2
!!$          endif
!!$
!!$          pos(1:4) = (/ 0d0,1d0,2d0,3d0 /)
!!$
!!$          if (igl == 0) then
!!$             pos(5) = -0.5d0
!!$             val(1:4) = vmec_arr(igl+1:igl+4,jgl,kgl)
!!$             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
!!$          elseif (igl == 1) then
!!$             pos(5) = 0.5d0
!!$             val(1:4) = vmec_arr(igl:igl+3,jgl,kgl)
!!$             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
!!$          elseif (igl == ns_i) then
!!$             pos(5) = 3.5d0
!!$             val(1:4) = vmec_arr(igl-3:igl,jgl,kgl)
!!$             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
!!$          elseif (igl == ns_i - 1) then
!!$             pos(5) = 2.5d0
!!$             val(1:4) = vmec_arr(igl-2:igl+1,jgl,kgl)
!!$             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
!!$          else
!!$             pos(5) = 1.5d0
!!$             val(1:4) = vmec_arr(igl-1:igl+2,jgl,kgl)
!!$             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
!!$          endif
!!$
!!$          local_val = val(5)
!!$
!!$        end subroutine half_to_int

      END MODULE efit_mod

!     efit_map
!     #################################################################
      subroutine efit_map(equ_file,g_def)

!     -----------------------------------------------------------------
!     Give Cartesian coordinates of each logical mesh point at grid
!     level (igrid).
!     -----------------------------------------------------------------

      use efit_mod
      use grid
      use equilibrium

      implicit none

!     Input variables

      type(grid_mg_def),pointer :: g_def
      character(*) :: equ_file

!     Local variables

      integer :: igrid,nx,ny,nz,alloc_stat,icomp,ierr
      integer :: nxg,nyg,nzg,i,j,k,igl,jgl,kgl,ig,jg,kg
      real(8) :: r1,r2,r3,th1,v1,ph1,sgn,ds,dth,dphi,rr1(1),zz1(1),ph0

      real(8),allocatable,dimension(:,:,:,:) :: xcar

      !SLATEC spline variables
      integer :: kx,ky,kz,nxs,nys,nzs,dim,flg,sorder
      real(8),dimension(:)    ,allocatable :: tx,ty,tz,work,xs,ys,zs
      real(8),dimension(:,:,:),allocatable :: rr_coef,zz_coef

      real(8)  :: db3val
      external    db3val

!     Begin program

!     Cycle grid levels

      do igrid=1,g_def%ngrid

!     Get LOCAL limits and allocate local map array

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        allocate(xcar(0:nx+1,0:ny+1,0:nz+1,3))

!     Get GLOBAL limits (EFIT operates on global domain)

        nxg = g_def%nxgl(igrid)
        nyg = g_def%nygl(igrid)
        nzg = g_def%nzgl(igrid)

        if (my_rank == 0) then
           write (*,'(a,i3,a,i3,a,i3,a,i3)') &
          ' Reading EFIT map on grid',igrid,', nx x ny x nz=',nxg,'x',nyg,'x',nzg
        endif

!     Read equilibrium file and setup arrays
!     [VMEC++ assumes solution is up-down symmetric wrt Z=0, 
!     and hence only gives theta=(0,pi); thereby the limit nyg/2+1 in theta]

        call efit_init(nxg,nyg,nzg,equ_file)

!     Spline EFIT global map

        !Prepare splines
        sorder = 2
        nxs = nxg+1
        nys = nyg/2+1
        nzs = nzg+2

        flg = 0 !Let spline routine find interpolation knots
        kx = min(sorder+1,nxs-1)
        ky = min(sorder+1,nys-1)
        kz = min(sorder+1,nzs-1)

        dim = nxs*nys*nzs + max(2*kx*(nxs+1),2*ky*(nys+1),2*kz*(nzs+1))

        allocate(work(dim) ,stat=alloc_stat)
        allocate(tx(nxs+kx),stat=alloc_stat)
        allocate(ty(nys+ky),stat=alloc_stat)
        allocate(tz(nzs+kz),stat=alloc_stat)

        !Global VMEC coordinates 
        allocate(xs(nxs),ys(nys),zs(nzs),stat=alloc_stat)

        ds = 1d0/nxg
        do ig = 1,nxs
          xs(ig) = ds*(ig-1)
        enddo

        dth = 2*pi/nyg
        do jg = 1,nys
          ys(jg) = dth*(jg-1)
        enddo

        dphi = 2*pi/nzg
        do kg = 1,nzs
          zs(kg) = dphi*(kg-2)
        enddo

        !Spline RR, ZZ
        allocate(rr_coef(nxs,nys,nzs),zz_coef(nxs,nys,nzs),stat=alloc_stat)

        call db3ink(xs,nxs,ys,nys,zs,nzs,rr,nxs,nys,kx,ky,kz,tx,ty,tz,rr_coef,work,flg)
        call db3ink(xs,nxs,ys,nys,zs,nzs,zz,nxs,nys,kx,ky,kz,tx,ty,tz,zz_coef,work,flg)

!     Transfer map (from GLOBAL in VMEC to LOCAL)

        ph0 = 0

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1

              call getMGmap(g_def,i,j,k,igrid,igrid,igrid,ig,jg,kg)

              !Find local coordinates
              r1  = g_def%xx(ig)
              th1 = g_def%yy(jg)
              ph1 = g_def%zz(kg)

              !Impose SP BCs
              if (r1 < 0d0) then
                 r1  = -r1
                 th1 = th1 + pi
              endif

              !Impose periodic BCs in theta, phi
              if (th1 < 0d0 ) th1 = th1 + 2*pi
              if (th1 > 2*pi) th1 = th1 - 2*pi

              if (ph1 < 0d0 ) ph1 = ph1 + 2*pi
              if (ph1 > 2*pi) ph1 = ph1 - 2*pi

              !Impose VMEC symmetries in theta, phi: R(th,phi) = R(2pi-th,2pi-phi)
              !                                      Z(th,phi) =-Z(2pi-th,2pi-phi)
              sgn = 1d0
              if (th1 > pi) then
                 th1 = 2*pi - th1
                 ph1 = 2*pi - ph1
                 sgn = -1d0
              endif

              !Map phi to VMEC toroidal coordinate, v
              v1  = mod(ph1,2*pi)

              !Extrapolate to ghost cell at psi=1 boundary (second-order)
              if (isBdry(g_def,i-1,igrid,2)) then
                 r1 = 0.5*(g_def%xx(ig-1)+g_def%xx(ig))
                 r2 =      g_def%xx(ig-1)
                 r3 =      g_def%xx(ig-2)
                 rr1(1) = 8./3.*db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                                ,kx,ky,kz,rr_coef,work)                &
     &                      -2.*db3val(r2,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                                ,kx,ky,kz,rr_coef,work)                &
     &                   +1./3.*db3val(r3,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                                ,kx,ky,kz,rr_coef,work)

                 zz1(1) = 8./3.*db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                                ,kx,ky,kz,zz_coef,work)                &
     &                      -2.*db3val(r2,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                                ,kx,ky,kz,zz_coef,work)                &
     &                   +1./3.*db3val(r3,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                                ,kx,ky,kz,zz_coef,work)
              else
                 rr1(1) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs    &
     &                          ,kx,ky,kz,rr_coef,work)

                 zz1(1) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs    &
     &                          ,kx,ky,kz,zz_coef,work)
              endif

              !Transform to Cartesian geometry (minus sign in phi to preserve a right-handed ref. sys.)
              ph1 =-g_def%zz(kg) + ph0
              xcar(i,j,k,1)=rr1(1)*cos(ph1)
              xcar(i,j,k,2)=rr1(1)*sin(ph1)
              xcar(i,j,k,3)=sgn*zz1(1)
            enddo
          enddo
        enddo

!     Fill grid metrics hierarchy

        call defineGridMetric(g_def,xcar=xcar,igr=igrid,ierr=ierr)

        if (check_grid.and.igrid==1) call checkGrid(g_def)

!     Free work space (to allow processing of different grid levels)

        deallocate(xcar,xs,ys,zs,work,tx,ty,tz,rr_coef,zz_coef)

        call efit_cleanup

        if (ierr /= 0) exit

      enddo

!     Check for errors, else set up gmetric MG hierarchy

      if (ierr /= 0) then
        if (my_rank == 0) write (*,*) ' >>> Discarding leftover VMEC meshes due to jac < 0'
        g_def%ngrid = igrid - 1
        call deallocateGridStructure(g_def%g_crse_def)
      else
        if (my_rank == 0) write (*,*) ' >>> Coarsening VMEC mesh hierarchy'
        call createMGMetricHierarchy(g_def)
      endif

!     End program

      end subroutine efit_map

!!$!     vmec_equ
!!$!     #################################################################
!!$      subroutine vmec_equ(iout,igrid,nx,ny,nz,bb,prs,rho,gam,equ_file &
!!$     &                   ,dcon,divcl)
!!$
!!$!     -----------------------------------------------------------------
!!$!     Give equilibrium fields at each logical mesh point in grid
!!$!     level (igrid).
!!$!     -----------------------------------------------------------------
!!$
!!$        use vmec_mod, pi_stel => pi
!!$        use grid
!!$        use equilibrium
!!$
!!$        implicit none
!!$
!!$!      Call variables
!!$
!!$        integer :: igrid,iout,nx,ny,nz
!!$        real(8) :: gam
!!$        real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3) :: bb
!!$        real(8),dimension(0:nx+1,0:ny+1,0:nz+1)   :: prs,rho
!!$        character(*) :: equ_file
!!$
!!$        logical :: dcon,divcl
!!$
!!$!     Local variables
!!$
!!$        integer :: nxg,nyg,nzg,i,j,k,igl,jgl,kgl,ig,jg,kg,istat,its
!!$        real(8) :: r1,z1,th1,ph1,v1,dphi,dth,sgn,jac,ds,dum1,dum2,ppi,max_rho
!!$
!!$        real(8),allocatable, dimension(:)     :: pflx,ff_i,q_i
!!$        real(8),allocatable, dimension(:,:,:) :: bsub3,prsg
!!$        real(8),allocatable, dimension(:,:,:,:) :: bsup
!!$
!!$        integer :: udcon=1111
!!$
!!$        logical :: enf_tor_flx_fn  !Whether we enforce flux functions in toroidal geom.
!!$
!!$        !SLATEC spline variables
!!$        integer :: kx,ky,kz,nxs,nys,nzs,dim,flg,sorder,alloc_stat
!!$        real(8),dimension(:)    ,allocatable :: tx,ty,tz,work,xs,ys,zs
!!$        real(8),dimension(:,:,:),allocatable :: b1_coef,b2_coef,b3_coef,prs_coef
!!$
!!$        real(8)  :: db3val
!!$        external :: db3val
!!$
!!$!     Begin program
!!$
!!$        if (my_rank == 0) then
!!$           write (*,*)
!!$           write (*,*) 'Reading VMEC solution...'
!!$        endif
!!$
!!$!     Get GLOBAL limits (VMEC operates on global domain)
!!$
!!$
!!$        nxg = gv%gparams%nxgl(igrid)
!!$        nyg = gv%gparams%nygl(igrid)
!!$        nzg = gv%gparams%nzgl(igrid)
!!$
!!$!     Read equilibrium file and setup arrays
!!$!     [Assumes solution is up-down symmetric wrt Z=0, 
!!$!     and hence only gives theta=(0,pi); thereby the limit nyg/2+1 in theta]
!!$
!!$        call vmec_init(nxg+1,nyg/2+1,nzg,equ_file,.true.)
!!$
!!$!     Setup equilibrium fields (at VMEC integer mesh -- PIXIE3D's half mesh)
!!$
!!$        call vmec_init_fields
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
!!$
!!$!     Find half-mesh PIXIE3D magnetic field components in GLOBAL mesh (flux functions)
!!$
!!$        allocate(bsup (0:nxg+1,0:nyg+1,0:nzg+1,3))
!!$        allocate(bsub3(0:nxg+1,0:nyg+1,0:nzg+1))
!!$        allocate(prsg (0:nxg+1,0:nyg+1,0:nzg+1))
!!$
!!$        do k=0,nzg+1
!!$          do j=0,nyg+1
!!$            do i=0,nxg+1
!!$               igl = i
!!$               jgl = j
!!$               kgl = k
!!$
!!$              !Periodic boundary in theta=0 (other boundary enforced by symmetry)
!!$              if (jgl == 0) jgl = nyg
!!$
!!$              !Periodic boundary in phi=0
!!$              if (kgl == 0) kgl = nzg
!!$
!!$              !Up-down symmetry in theta, phi: R(th,phi) = R(2pi-th,2pi-phi)
!!$              !                                Z(th,phi) =-Z(2pi-th,2pi-phi)
!!$              sgn = 1d0
!!$              if (jgl > nyg/2+1) then
!!$                 jgl = nyg + 2 - jgl
!!$
!!$                 kgl = nzg + 2 - kgl
!!$                 sgn = -1d0
!!$              endif
!!$
!!$              !Periodic boundary in phi=2*pi
!!$              if (kgl == nzg+1) kgl = 1
!!$
!!$              !Interpolate VMEC fields to PIXIE3D's integer radial mesh
!!$              call half_to_int(igl,jgl,kgl,bsupsijsf,bsup(i,j,k,1))
!!$              call half_to_int(igl,jgl,kgl,bsupuijcf,bsup(i,j,k,2))
!!$              call half_to_int(igl,jgl,kgl,bsupvijcf,bsup(i,j,k,3))
!!$              call half_to_int(igl,jgl,kgl,bsubvijcf,bsub3(i,j,k))
!!$              call half_to_int(igl,jgl,kgl,presijf  ,prsg (i,j,k))
!!$
!!$              !Find flux functions w/ PIXIE3D convention
!!$              jac = sqrtg(igl+1,jgl,kgl) !VMEC's half-mesh jacobian, including s=sqrt(psi) transformation
!!$              bsup(i,j,k,1) = jac*bsup(i,j,k,1)/(2.*gv%gparams%xg(igl)) !B1=B1'/(2s)
!!$              bsup(i,j,k,2) = jac*bsup(i,j,k,2)!Flux coordinate, B2=B2'
!!$              bsup(i,j,k,3) = jac*bsup(i,j,k,3)!Flux coordinate, B3=B3'
!!$            enddo
!!$          enddo
!!$        enddo
!!$
!!$
!!$        !Calculate maximum density (to normalize density later)
!!$        max_rho = maxval(prsg(1:nxg,1:nyg,1:nzg)**(1d0/gam))
!!$
!!$        !Clean flux functions (Tokamak case)
!!$        enf_tor_flx_fn = (nfp_i == 1)
!!$
!!$        divcl = .not.enf_tor_flx_fn !Whether to divergence clean
!!$        
!!$        if (enf_tor_flx_fn) then
!!$          do k=0,nzg+1
!!$            do i=0,nxg+1
!!$               igl = i
!!$               jgl = j
!!$               kgl = k
!!$
!!$              !Periodic boundary in theta=0 (other boundary enforced by symmetry)
!!$              if (jgl == 0) jgl = nyg
!!$
!!$              !Periodic boundary in phi=0
!!$              if (kgl == 0) kgl = nzg
!!$
!!$              !Up-down symmetry in theta, phi: R(th,phi) = R(2pi-th,2pi-phi)
!!$              !                                Z(th,phi) =-Z(2pi-th,2pi-phi)
!!$              sgn = 1d0
!!$              if (jgl > nyg/2+1) then
!!$                 jgl = nyg + 2 - jgl
!!$
!!$                 kgl = nzg + 2 - kgl
!!$                 sgn = -1d0
!!$              endif
!!$
!!$              !Periodic boundary in phi=2*pi
!!$              if (kgl == nzg+1) kgl = 1
!!$
!!$              dum1=sum(bsup (i,1:nyg,k,2))/nyg
!!$              dum2=sum(bsub3(i,1:nyg,k))/nyg
!!$
!!$              bsup (i,:,k,2) = dum1
!!$              bsub3(i,:,k)   = dum2
!!$             
!!$            enddo
!!$          enddo
!!$        endif
!!$
!!$!     Spline PIXIE3D global variables
!!$
!!$        !Prepare splines
!!$        sorder = 2
!!$        nxs = nxg+2
!!$        nys = nyg+2
!!$        nzs = nzg+2
!!$
!!$        flg = 0 !Let spline routine find interpolation knots
!!$        kx = min(sorder+1,nxs-1)
!!$        ky = min(sorder+1,nys-1)
!!$        kz = min(sorder+1,nzs-1)
!!$
!!$        dim = nxs*nys*nzs + max(2*kx*(nxs+1),2*ky*(nys+1),2*kz*(nzs+1))
!!$
!!$        allocate(work(dim) ,stat=alloc_stat)
!!$        allocate(tx(nxs+kx),stat=alloc_stat)
!!$        allocate(ty(nys+ky),stat=alloc_stat)
!!$        allocate(tz(nzs+kz),stat=alloc_stat)
!!$
!!$        !Global coordinates 
!!$        allocate(xs(nxs),ys(nys),zs(nzs),stat=alloc_stat)
!!$
!!$        !Radial half-mesh (PIXIE3D's full mesh)
!!$        ds = 1d0/nxg
!!$        do ig = 1,nxs
!!$!!           xs(ig) = ds*(ig-1)
!!$           xs(ig) = ds*(ig-1.5)
!!$        enddo
!!$
!!$        !Full angle meshes, starting at angle 0
!!$        dth = 2*pi/nyg
!!$        do jg = 1,nys
!!$           ys(jg) = dth*(jg-1)
!!$        enddo
!!$
!!$        dphi = 2*pi/nzg/nfp_i  !VMEC stores 1 of NFP periods in phi, ie., v=2*pi/nfp
!!$        do kg = 1,nzs
!!$           zs(kg) = dphi*(kg-2)
!!$        enddo
!!$
!!$
!!$        !Spline B-field components
!!$        allocate(b1_coef (nxs,nys,nzs)  &
!!$     &          ,b2_coef (nxs,nys,nzs)  &
!!$     &          ,b3_coef (nxs,nys,nzs)  &
!!$     &          ,prs_coef(nxs,nys,nzs),stat=alloc_stat)
!!$
!!$        call db3ink(xs,nxs,ys,nys,zs,nzs,bsup(:,:,:,1),nxs,nys,kx,ky,kz,tx,ty,tz,b1_coef,work,flg)
!!$        call db3ink(xs,nxs,ys,nys,zs,nzs,bsup(:,:,:,2),nxs,nys,kx,ky,kz,tx,ty,tz,b2_coef,work,flg)
!!$        if (enf_tor_flx_fn) then
!!$          call db3ink(xs,nxs,ys,nys,zs,nzs,bsub3,nxs,nys,kx,ky,kz,tx,ty,tz,b3_coef,work,flg)
!!$        else
!!$          call db3ink(xs,nxs,ys,nys,zs,nzs,bsup(:,:,:,3),nxs,nys,kx,ky,kz,tx,ty,tz,b3_coef,work,flg)
!!$        endif
!!$        call db3ink(xs,nxs,ys,nys,zs,nzs,prsg ,nxs,nys,kx,ky,kz,tx,ty,tz,prs_coef ,work,flg)
!!$
!!$!     Transfer variables (from GLOBAL in VMEC to LOCAL in PIXIE3D)
!!$
!!$        do k=0,nz+1
!!$          do j=0,ny+1
!!$            do i=0,nx+1
!!$
!!$              call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)
!!$
!!$              !Find global limits
!!$              call fromLocalToGlobalLimits(gv%gparams,igrid,i,j,k,igl,jgl,kgl)
!!$
!!$              !Find local coordinates
!!$              r1  = gv%gparams%xx(ig)
!!$              th1 = gv%gparams%yy(jg)
!!$              ph1 = gv%gparams%zz(kg)
!!$
!!$              !Impose SP BCs
!!$              if (r1 < 0d0) then
!!$                 r1  = -r1
!!$                 th1 = th1 + pi
!!$              endif
!!$
!!$              !Impose periodic BCs in theta, phi
!!$              if (th1 < 0d0 ) th1 = th1 + 2*pi
!!$              if (th1 > 2*pi) th1 = th1 - 2*pi
!!$
!!$              if (ph1 < 0d0 ) ph1 = ph1 + 2*pi
!!$              if (ph1 > 2*pi) ph1 = ph1 - 2*pi
!!$
!!$              !Map phi to VMEC toroidal coordinate, v
!!$              v1  = mod(ph1,2*pi/nfp_i)
!!$
!!$              bb(i,j,k,1) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
!!$     &                          ,kx,ky,kz,b1_coef,work)
!!$
!!$              bb(i,j,k,2) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
!!$     &                          ,kx,ky,kz,b2_coef,work)             !Flux coordinate
!!$
!!$              if (enf_tor_flx_fn) then
!!$                bb(i,j,k,3) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
!!$     &                            ,kx,ky,kz,b3_coef,work)             !Cov
!!$                bb(i,j,k,3) = (bb(i,j,k,3)                                        &
!!$     &                    -(gv%gparams%gmetric%grid(igrid)%gsub(i,j,k,3,1)*bb(i,j,k,1)   &
!!$     &                     +gv%gparams%gmetric%grid(igrid)%gsub(i,j,k,3,2)*bb(i,j,k,2))) &
!!$     &                    /gv%gparams%gmetric%grid(igrid)%gsub(i,j,k,3,3)       !Xform to cnv
!!$              else
!!$                bb(i,j,k,3) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
!!$     &                            ,kx,ky,kz,b3_coef,work)             !Cnv
!!$              endif
!!$
!!$
!!$              prs(i,j,k) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs &
!!$     &                           ,kx,ky,kz,prs_coef,work)
!!$
!!$              rho(i,j,k) = abs(prs(i,j,k))**(1d0/gam)  !Forces positive density
!!$            enddo
!!$          enddo
!!$        enddo
!!$
!!$        if (max_rho == 0d0) then
!!$          rho = 1d0
!!$        else
!!$          rho = rho/max_rho
!!$        endif
!!$
!!$!     Free work space
!!$
!!$        DEALLOCATE(bsup,bsub3,prsg,stat=istat)
!!$
!!$        DEALLOCATE(bsupuijcf,bsupvijcf,bsubvijcf,bsupsijsf,presijf,stat=istat)
!!$
!!$        DEALLOCATE(xs,ys,zs,work,tx,ty,tz,b1_coef,b2_coef,b3_coef,prs_coef,stat=istat)
!!$
!!$        IF (istat .ne. 0) STOP 'Deallocation error in vmec_equ'
!!$
!!$        call vmec_cleanup
!!$
!!$        call vmec_cleanup_metrics
!!$
!!$!     End program
!!$
!!$      end subroutine vmec_equ

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
