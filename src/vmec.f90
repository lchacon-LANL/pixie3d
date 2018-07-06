      MODULE vmec_mod

      USE stel_kinds
      USE stel_constants
      USE island_params
      USE read_wout_mod, ns_vmec=>ns, mpol_vmec=>mpol, ntor_vmec=>ntor, &
     &   rmnc_vmec=>rmnc, zmns_vmec=>zmns,                              &
     &   xm_vmec=>xm, xn_vmec=>xn, iotaf_vmec=>iotaf, phipf_vmec=>phipf,&
     &   presf_vmec=>presf, nfp_vmec=>nfp, bsupumnc_vmec=>bsupumnc,     &
     &   bsupvmnc_vmec=>bsupvmnc, wb_vmec=>wb,bsubvmnc_vmec=>bsubvmnc

      IMPLICIT NONE
!     
!     WRITTEN 06-27-06 BY S. P. HIRSHMAN AS PART OF THE VMEC++ PROJECT (c)
!     MODIFIED 11-28-06 BY L. CHACON AS PART OF THE PIXIE3D PROJECT (c)
!     
!     PURPOSE: COMPUTES AND STORES THE REAL-SPACE METRIC ELEMENTS, JACOBIAN BASED 
!     ON A SQRT(FLUX) SPLINED VMEC - COORDINATE SYSTEM
!

!     VARIABLE DECLARATIONS (3D for now; will convert to 2D or 1D as needed)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                     &
     &             sqrtg,                                               &!sqrt(g): Jacobian on half grid
     &             gss, gsu, gsv, guu, guv, gvv,                        &!symmetric elements of lower metric tensor (full mesh)
     &             hss, hsu, hsv, huu, huv, hvv                          !symmetric elements of upper metric tensor (half mesh)

      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
     &             rmnc_spline, zmns_spline, bsupumnc_spline,           &
     &             bsupvmnc_spline, bsubvmnc_spline

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                     &
     &             rmnc_i, zmns_i, bsupumncf_i, bsupvmncf_i, bsubvmncf_i

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: rr,zz

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: bsupuijcf,          &
     &  bsupvijcf, bsubvijcf, bsupsijsf, presijf

!     LOCAL (PRIVATE) HELPER ROUTINES (NOT ACCESSIBLE FROM OUTSIDE THIS MODULE)
!
!LC 2/2/07      PRIVATE spline_fourier_modes, loadrzl_vmec,       &
      PRIVATE spline_fourier_modes,convert_to_full_mesh, vmec_load_rz

      CONTAINS

!       vmec_init
!       ################################################################################
        SUBROUTINE vmec_init(ns_in,nu_in,nv_in,wout_file,load_metrics)

          IMPLICIT NONE

!       -----------------------------------------------------------------------------
!       LOADS VALUES FOR MESHES TO BE USED IN VMECPP ISLAND SOLVER
!       (GENERALLY DIFFERENT FROM VMEC MESHES. USE SUBSCRIPT _i FOR ISLAND VARIABLES)
!       READS wout_file FROM VMEC TO GET FOURIER COMPONENTS ON VMEC MESH
!       SPLINES THE VMEC COMPONENTS TO RADIAL ISLAND MESH (PHI->SQRT(PHI))
!       CALLS Fourier MODULE fixarray TO LOAD TRIG ARRAYS FOR COMPUTING ISLAND METRICS
!       COMPUTES gij (sub/sup) AND JACOBIAN ON ISLAND MESHES IN REAL SPACE
!
!       Modified by L. Chacon 11/30/06, from routines by S. Hirschman and R. Sanchez
!       -----------------------------------------------------------------------------

!       Call variables

          CHARACTER*(*)            :: wout_file
          INTEGER, INTENT(in)      :: ns_in, nu_in, nv_in

!       Local variables

          INTEGER                  :: istat, ntype, imesh, js
          REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
     &                            bsupumnc_v,  bsupvmnc_v,  bsubvmnc_v
          REAL(rprec) :: t1, t2
          INTEGER     :: m, n, mn
          LOGICAL     :: load_metrics

!       Begin program

!       Initialize variables

          ns_i = ns_in
          nu_i = nu_in
          nv_i = nv_in

          nsh   = ns_i-1
          nuv_i = nu_i*nv_i

!       READ-IN DATA FROM VMEC-PRODUCED WOUT FILE (LIBSTELL ROUTINE)

          CALL read_wout_file(wout_file, istat)
          IF (istat .ne. 0) STOP 'Read-wout error in vmec_init'

          ntor_i = ntor_vmec   !Take from VMEC equilibrium file
          mpol_i = mpol_vmec   !Take from VMEC equilibrium file

          mnmax_i = (mpol_i + 1)*(2*ntor_i + 1)   ! Added RS. Contains total number of modes.

          nfp_i = nfp_vmec
          wb_i  = (4*pi*pi)*wb_vmec

          IF (wb_i .eq. 0._dp) STOP 'wb_vmec = 0!'

!       Allocate space for splined arrays

          ALLOCATE(rr(ns_i,nu_i,0:nv_i+1),zz(ns_i,nu_i,0:nv_i+1))

          ALLOCATE(rmnc_spline(mnmax,ns_i), zmns_spline(mnmax,ns_i),        &
     &             bsupumnc_v(mnmax,ns_vmec),bsupumnc_spline(mnmax,ns_i),   &
     &             bsupvmnc_v(mnmax,ns_vmec),bsupvmnc_spline(mnmax,ns_i),   &
     &             bsubvmnc_v(mnmax,ns_vmec),bsubvmnc_spline(mnmax,ns_i),   &
     &             phipf_i(ns_i), iotaf_i(ns_i),presf_i(ns_i),              &
     &             stat=istat)
          IF (istat .ne. 0) STOP 'Allocation error 1 in VMEC_INIT'

!       SPLINE R, Z. L FOURIER COMPONENTS in s FROM ORIGINAL VMEC MESH (s ~ phi, ns_vmec points) 
!       TO A "POLAR" MESH [s ~ sqrt(phi), ns_i POINTS] WITH BETTER AXIS RESOLUTION

          DO ntype = 1, 3
             IF (ntype .eq. 1) THEN
                istat = 0
                CALL Spline_Fourier_Modes(rmnc_vmec, rmnc_spline, istat)
             ELSE IF (ntype .eq. 2) THEN
                istat = 0
                CALL Spline_Fourier_Modes(zmns_vmec, zmns_spline, istat)   
              ELSE 
                CALL Convert_From_Nyq(bsupumnc_v,bsupumnc_vmec)       !These have mnmax_nyq members
                CALL Convert_To_Full_Mesh(bsupumnc_v)
!!                istat = 1
                istat = 0    !Changed by L. Chacon, 6/5/07
                CALL Spline_Fourier_Modes(bsupumnc_v, bsupumnc_spline, istat)   

!!L. Chacon 9/10/07:  Added cov component of B_tor (flux function)
                CALL Convert_From_Nyq(bsupvmnc_v,bsupvmnc_vmec)       !These have mnmax_nyq members
                CALL Convert_From_Nyq(bsubvmnc_v,bsubvmnc_vmec)       !These have mnmax_nyq members
                CALL Convert_To_Full_Mesh(bsupvmnc_v)
                CALL Convert_To_Full_Mesh(bsubvmnc_v)
!!L. Chacon 6/5/07                istat = 1
                istat = 0    !Changed by L. Chacon, 6/5/07
                CALL Spline_Fourier_Modes(bsupvmnc_v, bsupvmnc_spline, istat)   
                CALL Spline_Fourier_Modes(bsubvmnc_v, bsubvmnc_spline, istat)   
                DEALLOCATE(bsupumnc_v, bsupvmnc_v, bsubvmnc_v)
             END IF

             IF (istat .ne. 0) STOP 'Spline error in VMEC_INIT'
          END DO

!       Spline 1-D arrays: careful -> convert phipf VMEC and multiply
!       by ds-vmec/ds-island, since phipf_i = d(PHI)/ds-island

          CALL Spline_OneD_Array (iotaf_vmec, iotaf_i, istat)
          CALL Spline_OneD_Array (phipf_vmec, phipf_i, istat)
          presf_vmec = mu0 * presf_vmec
          CALL Spline_OneD_Array (presf_vmec, presf_i, istat)
          DO js = 1, ns_i
             phipf_i(js) = 2 * hs_i*(js-1) * phipf_i(js) / (2*pi)
          END DO
          
!     CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON "POLAR" MESH
!     AND COMPUTE METRIC ELEMENTS AND JACOBIAN

          CALL vmec_load_rz(load_metrics,istat)
          IF (istat .ne. 0) STOP 'Load R,Z error in VMEC_INIT'

        END SUBROUTINE vmec_init

!       vmec_load_rz
!       ################################################################################
        SUBROUTINE vmec_load_rz(load_metrics,istat)

!       ------------------------------------------------------------------------------
!       Computes R-Z map in real space from VMEC equilibrium
!
!       Modified by L. Chacon 11/30/06 from routines by S. Hirschman and R. Sanchez
!       ------------------------------------------------------------------------------

          USE Fourier, ONLY: toijsp
!          USE dump_output

!       Call variables

          INTEGER, INTENT(out)     :: istat
          LOGICAL :: load_metrics

!       Local variables

          REAL(rprec), DIMENSION(ns_vmec)  :: fac1, fac2
          INTEGER                          :: mpol_save, ntor_save,iparity,ipol

!       Begin program

!       LOAD FOURIER FACTORS USING VMEC mpol,ntor AND ISLAND nu,nv
!       IDEA: USE VMEC mpol, ntor VALUES AND CALL fixarray
      
          CALL fixarray

!       FIRST, REPACK SPLINED ARRAYS FROM VMEC-ORDERING TO ISLAND-ORDERING
!       AND SWAP THE N->-N MODES TO BE CONSISTENT WITH mu+nv ISLAND ARGUMENT

          CALL repack (istat)
          IF (istat .ne. 0) STOP 'REPACK error in VMEC_LOAD_RZ'

!       COMPUTE AND STORE R, Z

          iparity = 0
          CALL toijsp (rmnc_i,rr(:,:,1:nv_i),0,0,iparity,0)
      
          iparity = 1
          CALL toijsp (zmns_i,zz(:,:,1:nv_i),0,0,iparity,0)

!       ENFORCE PERIODIC BCs ALONG TOROIDAL DIRECTION

          rr(:,:,0     ) = rr(:,:,nv_i)
          rr(:,:,nv_i+1) = rr(:,:,1   )

          zz(:,:,0     ) = zz(:,:,nv_i)
          zz(:,:,nv_i+1) = zz(:,:,1   )

!       DEBUG: CHECK IF CORRECT
!
!        istat = (ns_i-1)/2 + 1    !rho = 1/2 => s = 1/4
!        CALL dump_special(r1_i,z1_i,ru_i,zu_i,rv_i,zv_i,istat)

!       LOAD METRIC ELEMENTS

          if (load_metrics) then
            call vmec_load_metrics(istat)
            IF (istat .ne. 0) STOP 'Error in VMEC_LOAD_METRICS'
          endif

!       CLEAN-UP EXTRA ARRAYS

          DEALLOCATE (rmnc_i, zmns_i, stat=istat)
          IF (lwout_opened) CALL read_wout_deallocate

        END SUBROUTINE vmec_load_rz

!       vmec_load_metrics
!       ################################################################################
        SUBROUTINE vmec_load_metrics(istat)

!       ------------------------------------------------------------------------------
!       Loads metric elements in VMEC's half mesh (PIXIE3D integer mesh) in real space
!       from VMEC equilibrium
!
!       Modified by L. Chacon 01/31/07, from routines by S. Hirschman and R. Sanchez
!       ------------------------------------------------------------------------------

          USE Fourier, ONLY: toijsp
!          USE dump_output

!       Call variables

          INTEGER, INTENT(out)     :: istat
          REAL(8), DIMENSION(:,:,:), ALLOCATABLE ::ru_i,zu_i,rv_i,zv_i

!       Local variables

!!$          REAL(rprec), DIMENSION(ns_vmec)  :: fac1, fac2
          INTEGER                          :: iparity

!       Begin program

!       COMPUTE AND STORE ANGULAR DERIVATIVES OF R, Z

          ALLOCATE (ru_i(ns_i, nu_i, nv_i), zu_i(ns_i, nu_i, nv_i),         &
     &              rv_i(ns_i, nu_i, nv_i), zv_i(ns_i, nu_i, nv_i),         &
     &              stat = istat)

          IF (istat .ne. 0) STOP 'Allocation failed in VMEC_LOAD_METRICS'

          iparity = 0
          CALL toijsp (rmnc_i, ru_i, 1, 0, iparity, 0)
          CALL toijsp (rmnc_i, rv_i, 0, 1, iparity, 0)
      
          iparity = 1
          CALL toijsp (zmns_i, zu_i, 1, 0, iparity, 0)
          CALL toijsp (zmns_i, zv_i, 0, 1, iparity, 0)

!       DEBUG: CHECK IF CORRECT
!
!        istat = (ns_i-1)/2 + 1    !rho = 1/2 => s = 1/4
!        CALL dump_special(r1_i,z1_i,ru_i,zu_i,rv_i,zv_i,istat)

!       COMPUTE HALF-MESH LOWER/UPPER METRIC ELEMENTS AND JACOBIAN

          CALL vmec_half_mesh_metrics (rr(:,:,1:nv_i),ru_i,rv_i,zz(:,:,1:nv_i),zu_i,zv_i)

!       CLEAN-UP EXTRA ARRAYS

          DEALLOCATE (ru_i,zu_i,rv_i,zv_i,stat=istat)

        END SUBROUTINE vmec_load_metrics

!       vmec_cleanup
!       ################################################################################
        SUBROUTINE vmec_cleanup

          USE Fourier, ONLY: orthonorm

          IMPLICIT NONE

          INTEGER         :: istat

          DEALLOCATE(rr,zz,stat=istat)
          IF (istat .ne. 0) STOP 'Deallocation error 1 in vmec_cleanup'

          DEALLOCATE(phipf_i,iotaf_i,presf_i,stat=istat)
          IF (istat .ne. 0) STOP 'Deallocation error 2 in vmec_cleanup'

          DEALLOCATE(bsupumncf_i,bsupvmncf_i,bsubvmncf_i,stat=istat)
          IF (istat .ne. 0) STOP 'Deallocation error 3 in vmec_cleanup'

          call DEALLOC_FIXARRAY

          DEALLOCATE(orthonorm,stat=istat)
          IF (istat .ne. 0) STOP 'Deallocation error 4 in vmec_cleanup'

        END SUBROUTINE vmec_cleanup

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
!!$          REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
!!$     &                            bsupumnc_v,  bsupvmnc_v
!!$          REAL(rprec) :: t1, t2
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

      SUBROUTINE Spline_Fourier_Modes(ymn_vmec, ymn_spline, istat)

      use oned_int   !Added by L. Chacon 6/5/07

!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(inout)     :: istat
      REAL(rprec), DIMENSION(mnmax, ns_vmec), TARGET,                   &
     &                     INTENT(in)  :: ymn_vmec
      REAL(rprec), DIMENSION(mnmax, ns_i), TARGET,                      &
     &                     INTENT(out) :: ymn_spline
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER          :: one = 1
      INTEGER                         :: js, modes, ntype, mp
      REAL(rprec), DIMENSION(ns_vmec) :: y_vmec!, y_spline
      REAL(rprec), DIMENSION(ns_vmec) :: snodes_vmec, y2_vmec, fac1
      REAL(rprec), DIMENSION(ns_i)    :: snodes, fac2
      REAL(rprec)                     :: hs_vmec, yp1, ypn, expm
!-----------------------------------------------
!
!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM VMEC (s~phi) TO s~sqrt(phi) MESH
!

!
!     1. Factors for taking out (putting back) [sqrt(s)]**m factor

      IF (ns_vmec .le. 1) THEN
         istat = 1
         RETURN
      END IF

      fac1 = 0
      hs_vmec = one/(ns_vmec-1)
      DO js = 2, ns_vmec
         fac1(js) = one/SQRT(hs_vmec*(js-1))
      END DO

      IF (ns_i .le. 1) THEN
         istat = 2
         RETURN
      END IF

      ohs_i = ns_i-1
      hs_i = one/ohs_i
      DO js = 1, ns_i
         fac2(js) = hs_i*(js-1)
      END DO

!
!     2. Set up s-nodes in original (vmec) radial spatial coordinates for the final (splined) mesh 
!        using the mapping relation: s_vmec(sfinal) = sfinal**2 since sfinal = sqrt(phi)
!

      DO js = 1, ns_i
         snodes(js) = fac2(js)**2
      END DO

!
!     3. Set up "knots" on initial (vmec, svmec~phi) mesh 
!
      DO js = 1, ns_vmec
         snodes_vmec(js) = hs_vmec*((js-1))
      END DO

      ymn_spline = 0d0   !L. Chacon, 2/20/07, to avoid definition error below

      DO modes = 1, mnmax
         y_vmec   = ymn_vmec(modes,:)
!         y_spline = ymn_spline(modes,:)
         mp = xm_vmec(modes)

         IF (istat.eq.0 .and. mp.gt.0) THEN
            IF (MOD(mp,2) .eq. 1) THEN 
               expm = 1d0
            ELSE
               expm = 2d0
            END IF
            y_vmec = y_vmec*(fac1**expm)
!LC 6/5/07            IF (mp .le. 2) y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
            IF (mp .le. 2) then
               call IntDriver1d(4,snodes_vmec(2:5),y_vmec(2:5)   &
     &                         ,1,snodes_vmec(1:1),y_vmec(1:1),3)
            ENDIF
         END IF

!
!        4. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!
!!$         yp1 = -1.e30_dp;  ypn = -1.e30_dp
!!$         CALL spline (snodes_vmec, y_vmec, ns_vmec, yp1, ypn, y2_vmec)
!!$
!!$!
!!$!        5. Interpolate onto snodes mesh
!!$!
!!$         DO js = 1, ns_i
!!$            CALL splint (snodes_vmec, y_vmec, y2_vmec, ns_vmec,         &
!!$     &                   snodes(js), y_spline(js))
!!$         END DO

!
!        4 and 5: spline vmec coefficients into snodes mesh (L. Chacon, 6/5/07)
!
         call IntDriver1d(ns_vmec,snodes_vmec,y_vmec   &
     &                   ,ns_i   ,snodes     ,ymn_spline(modes,:),3)

!
!        PRINT OUT FOR CHECKING (Modified by L. Chacon 6/5/07)
!
!!$         IF (xn_vmec(modes) .eq. 0) THEN
!!$             WRITE (*, *) mp
!!$             WRITE (*, *) 'index   Spline_knots     Vmec positions'
!!$             DO js = 1, ns_vmec
!!$                WRITE (*, '(i3,1p2e14.5)') js,snodes_vmec(js), y_vmec(js)
!!$             END DO
!!$             WRITE (*, *)
!!$             WRITE (*, *) 'index   Radial nodes     Splined positions'
!!$             DO js = 1, ns_i
!!$                WRITE (*, '(i3,1p2e14.5)') js,snodes(js), y_spline(js)
!!$             END DO
!!$             WRITE (*, *)
!!$         END IF
!
!        RECOVER RADIAL DEPENDENCE
!
         IF (istat.eq.0 .and. mp.gt.0) THEN
            ymn_spline(modes,:) = ymn_spline(modes,:)*(fac2**expm)
         END IF

      END DO

!LC      stop
      
      istat = 0

      END SUBROUTINE Spline_Fourier_Modes

!!$      SUBROUTINE Spline_Fourier_Modes(ymn_vmec, ymn_spline, istat)
!!$!-----------------------------------------------
!!$!   D u m m y   A r g u m e n t s
!!$!-----------------------------------------------
!!$      INTEGER, INTENT(inout)     :: istat
!!$      REAL(rprec), DIMENSION(mnmax, ns_vmec), TARGET,                   &
!!$     &                     INTENT(in)  :: ymn_vmec
!!$      REAL(rprec), DIMENSION(mnmax, ns_i), TARGET,                      &
!!$     &                     INTENT(out) :: ymn_spline
!!$!-----------------------------------------------
!!$!   L o c a l   V a r i a b l e s
!!$!-----------------------------------------------
!!$      REAL(rprec), PARAMETER          :: one = 1
!!$      INTEGER                         :: js, modes, ntype, mp
!!$      REAL(rprec), DIMENSION(:), POINTER :: y_vmec, y_spline
!!$      REAL(rprec), DIMENSION(ns_vmec) :: snodes_vmec, y2_vmec, fac1
!!$      REAL(rprec), DIMENSION(ns_i)    :: snodes, fac2
!!$      REAL(rprec)                     :: hs_vmec, yp1, ypn
!!$!-----------------------------------------------
!!$!
!!$!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM VMEC (s~phi) TO s~sqrt(phi) MESH
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
!!$      ymn_spline = 0d0   !L. Chacon, 11/29/06, to avoid definition error below
!!$
!!$      DO modes = 1, mnmax
!!$         y_vmec => ymn_vmec(modes,:)
!!$         y_spline => ymn_spline(modes,:)
!!$         mp = xm_vmec(modes)
!!$
!!$         IF (istat.eq.0 .and. MOD(mp,2).eq.1) THEN
!!$            y_vmec = y_vmec*fac1
!!$            y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
!!$         END IF
!!$
!!$!
!!$!        4. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!!$!
!!$         yp1 = -1.e30_dp;  ypn = -1.e30_dp
!!$         CALL spline (snodes_vmec, y_vmec, ns_vmec, yp1, ypn, y2_vmec)
!!$
!!$!
!!$!        5. Interpolate onto snodes mesh
!!$!
!!$         DO js = 1, ns_i
!!$            CALL splint (snodes_vmec, y_vmec, y2_vmec, ns_vmec,         &
!!$     &                   snodes(js), y_spline(js))
!!$         END DO
!!$
!!$         IF (istat.eq.0 .and. MOD(mp,2).eq.1) THEN
!!$            y_spline = y_spline*fac2
!!$         END IF
!!$
!!$!
!!$!        PRINT OUT FOR CHECKING
!!$!
!!$         IF (xn_vmec(modes) .eq. 0) THEN
!!$!             WRITE (33, *) mp
!!$             DO js = 1, ns_vmec
!!$!                WRITE (33, '(1p2e14.4)') snodes_vmec(js), y_vmec(js)
!!$             END DO
!!$!                WRITE (33, *)
!!$             DO js = 1, ns_i
!!$!                WRITE (33, '(1p2e14.4)') snodes(js), y_spline(js)
!!$             END DO
!!$!                WRITE (33, *)
!!$         END IF
!!$
!!$      END DO
!!$
!!$      istat = 0
!!$
!!$      END SUBROUTINE Spline_Fourier_Modes
      

      SUBROUTINE Spline_OneD_Array (y_vmec, y_spline, istat)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(out)     :: istat
      REAL(rprec), DIMENSION(ns_vmec), INTENT(in)  :: y_vmec
      REAL(rprec), DIMENSION(ns_i)   , INTENT(out) :: y_spline
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER          :: one = 1
      INTEGER                         :: js, modes, ntype, mp
      REAL(rprec), DIMENSION(ns_vmec) :: snodes_vmec, y2_vmec, fac1
      REAL(rprec), DIMENSION(ns_i)    :: snodes, fac2
      REAL(rprec)                     :: hs_vmec, yp1, ypn
!-----------------------------------------------
!
!     CALL LIBRARY SPLINE ROUTINES TO SPLINE FROM VMEC (s~phi) TO s~sqrt(phi) MESH
!

!
!     1. Factors for taking out (putting back) sqrt(s) factor for m-odd modes

      IF (ns_vmec .le. 1) THEN
         istat = 1
         RETURN
      END IF

      fac1 = 0
      hs_vmec = one/(ns_vmec-1)
      DO js = 2, ns_vmec
         fac1(js) = one/SQRT(hs_vmec*(js-1))
      END DO

      IF (ns_i .le. 1) THEN
         istat = 2
         RETURN
      END IF

      ohs_i = ns_i-1
      hs_i = one/ohs_i
      DO js = 1, ns_i
         fac2(js) = hs_i*(js-1)
      END DO

!
!     2. Set up s-nodes on final (splined) mesh [sfinal ~ sqrt(s_vmec)]
!

      DO js = 1, ns_i
         snodes(js) = fac2(js)*fac2(js)
      END DO

!
!     3. Set up "knots" on initial (vmec, svmec~phi) mesh 
!
      DO js = 1, ns_vmec
         snodes_vmec(js) = hs_vmec*((js-1))
      END DO


!     4. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!
      yp1 = -1.e30_dp;  ypn = -1.e30_dp
      CALL spline (snodes_vmec, y_vmec, ns_vmec, yp1, ypn, y2_vmec)

!
!     5. Interpolate onto snodes mesh
!
      DO js = 1, ns_i
         CALL splint (snodes_vmec, y_vmec, y2_vmec, ns_vmec,         &
     &                   snodes(js), y_spline(js))
      END DO

      END SUBROUTINE Spline_OneD_Array

      SUBROUTINE Convert_From_Nyq(ymn_o, ymn_vmec)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(in)          :: ymn_vmec(mnmax_nyq,ns_vmec)
      REAL(rprec), INTENT(out)         :: ymn_o(mnmax,ns_vmec)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER  :: m, n, m_o, n_o, mn, mn_o
!-----------------------------------------------
      ymn_o = zero

      DO mn = 1, mnmax_nyq
         m = xm_nyq(mn);  n = xn_nyq(mn)
         DO mn_o = 1, mnmax
            IF (xm_vmec(mn_o).eq.m .and. xn_vmec(mn_o).eq.n) THEN
               ymn_o(mn_o, :) = ymn_vmec(mn, :)
               CYCLE
            END IF
         END DO
      END DO

      END SUBROUTINE Convert_From_Nyq
      
      SUBROUTINE Convert_To_Full_Mesh(ymn_vmec)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(inout)         :: ymn_vmec(mnmax,ns_vmec)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER             :: p5 = 0.5_dp
      INTEGER                            :: modes, js, mp
      REAL(rprec), DIMENSION(ns_vmec)    :: y_vmec
      REAL(rprec), DIMENSION(ns_vmec)    :: fac1, fac2
      REAL(rprec)                        :: temp
!-----------------------------------------------

!
!     CONVERTS INPUT ARRAY y_vmec FROM HALF-RADIAL MESH TO FULL MESH, EXTRAPOLATING
!     TO THE ORIGIN AND EDGE POINTS. STORE INTERPOLATED VALUES IN y_vmec
!

!
!     COMPUTE, STORE FACTORS FOR INTERPOLATING m ODD MODES
!
      fac1 = 0;  fac2 = 0
      DO js = 2, ns_vmec-1
         temp = REAL(js-1,dp)
         fac1(js) = p5*SQRT(temp/(js-1-p5))
         fac2(js) = p5*SQRT(temp/(js-p5))
      END DO

      DO modes = 1, mnmax
         y_vmec = ymn_vmec(modes,:)
         mp = xm_vmec(modes)

!
!     EXTRAPOLATE ORIGIN POINT FIRST
!
         IF (mp .eq. 0) THEN
            y_vmec(1) = p5*(3*y_vmec(2) - y_vmec(3))
         ELSE
            y_vmec(1) = 0
         END IF   

         IF (MOD(mp,2) .eq. 1) THEN
            DO js = 2, ns_vmec-1
               y_vmec(js) = fac1(js)*y_vmec(js) + fac2(js)*y_vmec(js+1)
            END DO
         ELSE
            DO js = 2, ns_vmec-1
               y_vmec(js) = p5*(y_vmec(js) + y_vmec(js+1))
            END DO
         END IF

!
!     EXTRAPOLATE EDGE POINT LAST
!
         y_vmec(ns_vmec) = 2*y_vmec(ns_vmec) - y_vmec(ns_vmec-1)
         ymn_vmec(modes,:) = y_vmec

      END DO

      END SUBROUTINE


!!$      SUBROUTINE LoadRZL_VMEC(istat)
!!$      USE Fourier, ONLY: toijsp
!!$!      USE dump_output
!!$!-----------------------------------------------
!!$!   D u m m y   A r g u m e n t s
!!$!-----------------------------------------------
!!$      INTEGER, INTENT(out)     :: istat
!!$!-----------------------------------------------
!!$!   L o c a l   V a r i a b l e s
!!$!-----------------------------------------------
!!$      REAL(rprec), DIMENSION(ns_vmec)  :: fac1, fac2
!!$      INTEGER                          :: mpol_save, ntor_save, iparity
!!$      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE    ::                  &
!!$     &   r1_i, z1_i, ru_i, zu_i, rv_i, zv_i
!!$      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: temp
!!$!-----------------------------------------------
!!$!
!!$!     1. LOAD FOURIER FACTORS USING VMEC mpol,ntor AND ISLAND nu,nv
!!$!        IDEA: USE VMEC mpol, ntor VALUES AND CALL fixarray
!!$!     2. COMPUTE METRIC ELEMENTS, JACOBIAN, LAMBDA ON ISLAND MESH
!!$      
!!$      CALL fixarray
!!$
!!$!
!!$!     COMPUTE R, Z (Metric elements, jacobian), LAMBDA (For initializing B)
!!$!
!!$!     FIRST, REPACK SPLINED ARRAYS FROM VMEC-ORDERING TO ISLAND-ORDERING
!!$!     AND SWAP THE N->-N MODES TO BE CONSISTENT WITH mu+nv ISLAND ARGUMENT
!!$!
!!$      CALL repack (istat)
!!$
!!$!
!!$!     COMPUTE AND STORE R, Z AND THEIR ANGULAR DERIVATIVES, AND LAMBDA (NEED FOR VECTOR POT)
!!$!
!!$      ALLOCATE (r1_i(ns_i, nu_i, nv_i), z1_i(ns_i, nu_i, nv_i),         &
!!$     &          ru_i(ns_i, nu_i, nv_i), zu_i(ns_i, nu_i, nv_i),         &
!!$     &          rv_i(ns_i, nu_i, nv_i), zv_i(ns_i, nu_i, nv_i),         &
!!$     &          stat = istat)
!!$
!!$      IF (istat .ne. 0) STOP 'Allocation failed in LoadRZL_VMEC'
!!$
!!$      iparity = 0
!!$      CALL toijsp (rmnc_i, r1_i, 0, 0, iparity, 0)
!!$      CALL toijsp (rmnc_i, ru_i, 1, 0, iparity, 0)
!!$      CALL toijsp (rmnc_i, rv_i, 0, 1, iparity, 0)
!!$      
!!$      iparity = 1
!!$      CALL toijsp (zmns_i, z1_i, 0, 0, iparity, 0)
!!$      CALL toijsp (zmns_i, zu_i, 1, 0, iparity, 0)
!!$      CALL toijsp (zmns_i, zv_i, 0, 1, iparity, 0)
!!$
!!$!
!!$!     DEBUG: CHECK IF CORRECT
!!$!
!!$      istat = (ns_i-1)/2 + 1    !rho = 1/2 => s = 1/4
!!$!      CALL dump_special(r1_i,z1_i,ru_i,zu_i,rv_i,zv_i,istat)
!!$
!!$!
!!$!     COMPUTE HALF-MESH LOWER/UPPER METRIC ELEMENTS AND JACOBIAN
!!$!
!!$      CALL half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)
!!$
!!$!
!!$!     CLEAN-UP EXTRA ARRAYS
!!$!
!!$      DEALLOCATE (rmnc_i, zmns_i, stat=istat)
!!$      IF (lwout_opened) CALL read_wout_deallocate
!!$
!!$!
!!!!$!     CONVERT UPPER METRIC ELEMENTS TO FULL MESH
!!!!$      ALLOCATE (temp(ns_i, nuv_i))
!!!!$      temp = hss
!!!!$      CALL to_full_mesh_sp (temp, hss)
!!!!$      temp = hsu
!!!!$      CALL to_full_mesh_sp (temp, hsu)
!!!!$      temp = hsv
!!!!$      CALL to_full_mesh_sp (temp, hsv)
!!!!$      temp = huu
!!!!$      CALL to_full_mesh_sp (temp, huu)
!!!!$      huu(3,:) = (9._dp)*huu(4,:)/4     !  huu ~ 1/rho**2
!!!!$      huu(2,:) = 4*huu(3,:)
!!!!$      temp = huv
!!!!$      CALL to_full_mesh_sp (temp, huv)
!!!!$      temp = hvv
!!!!$      CALL to_full_mesh_sp (temp, hvv)
!!!!$
!!!!$      DEALLOCATE(temp)
!!$
!!$      END SUBROUTINE LoadRZL_VMEC


      SUBROUTINE repack (istat)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(out)     :: istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                  :: modes, js, m, n, ntype, n1
      REAL, PARAMETER          :: nfactor2 = 1._dp/nfactor
!-----------------------------------------------
!
!     The splined arrays (rmnc_spline, etc) are all in VMEC ordering
!
!         rmnc_spline(1:mnmax,ns_i)
!
!     Now pack them so they can be used by ISLAND Fourier routines
!
!         rmnc_i(ns_i, 0:mpol, -ntor:ntor)
!
!     NOTE: mpol_i == mpol_vmec, ntor_i == ntor_vmec here
!

      ALLOCATE(rmnc_i(ns_i,0:mpol_i,-ntor_i:ntor_i),                    & 
     &         zmns_i(ns_i,0:mpol_i,-ntor_i:ntor_i),                    &
     &         bsupumncf_i(ns_i,0:mpol_i,-ntor_i:ntor_i),               &
     &         bsupvmncf_i(ns_i,0:mpol_i,-ntor_i:ntor_i),               &
     &         bsubvmncf_i(ns_i,0:mpol_i,-ntor_i:ntor_i),               &
     &         stat=istat)

      IF (istat .ne. 0) STOP 'Allocation error in REPACK'
      rmnc_i = 0;  zmns_i = 0;  bsupumncf_i = 0; bsupvmncf_i = 0 ; bsubvmncf_i =0

!
!     LOAD n>=0 ONLY FOR M=0
!

      IF (MAXVAL(xm_vmec).gt.mpol_i .or.                                &
     &    MAXVAL(ABS(xn_vmec))/nfp_vmec.gt.ntor_i)                      & 
     &    PRINT *,' You should increase number of modes to at least VMEC values!'
      

      DO js = 1,ns_i
         DO modes = 1,mnmax
            m = xm_vmec(modes)
            n = xn_vmec(modes)/nfp_vmec
            IF (m.gt.mpol_i .or. ABS(n).gt.ntor_i) CYCLE
            IF (m .eq. 0) THEN
               n1 = ABS(n)
               rmnc_i(js,m,n1) = rmnc_i(js,m,n1)                        &
      &                        + rmnc_spline(modes,js)
               zmns_i(js,m,n1) = zmns_i(js,m,n1)                        &
      &                        -SIGN(1,n)*zmns_spline(modes,js)
               bsupumncf_i(js,m,n1) = bsupumncf_i(js,m,n1)              &
      &                        + bsupumnc_spline(modes,js)
               bsupvmncf_i(js,m,n1) = bsupvmncf_i(js,m,n1)              &
      &                        + bsupvmnc_spline(modes,js)
               bsubvmncf_i(js,m,n1) = bsubvmncf_i(js,m,n1)              &
      &                        + bsubvmnc_spline(modes,js)
            ELSE
               rmnc_i(js,m,-n) = rmnc_spline(modes,js)
               zmns_i(js,m,-n) = zmns_spline(modes,js)
               bsupumncf_i(js,m,-n) = bsupumnc_spline(modes,js)
               bsupvmncf_i(js,m,-n) = bsupvmnc_spline(modes,js)
               bsubvmncf_i(js,m,-n) = bsubvmnc_spline(modes,js)
            ENDIF
         END DO
      END DO
!
! RS (10/03/06) = Add FACTOR for n=0/m=0 normalization    
! VMEC harmonics with m!=0 (or n!=0) VMEC are NFACTOR times bigger than the ones in VMEC++.                                 
!
      rmnc_i = nfactor2*rmnc_i
      rmnc_i(:,0,0) = rmnc_i(:,0,0)*nfactor
      zmns_i = nfactor2*zmns_i
      zmns_i(:,0,0) = zmns_i(:,0,0)*nfactor
      bsupumncf_i = nfactor2*bsupumncf_i
      bsupumncf_i(:,0,0) = bsupumncf_i(:,0,0)*nfactor
      bsupvmncf_i = nfactor2*bsupvmncf_i
      bsubvmncf_i = nfactor2*bsubvmncf_i
      bsupvmncf_i(:,0,0) = bsupvmncf_i(:,0,0)*nfactor
      bsubvmncf_i(:,0,0) = bsubvmncf_i(:,0,0)*nfactor

      DEALLOCATE (rmnc_spline, zmns_spline, bsupumnc_spline,                &
     &            bsupvmnc_spline, bsubvmnc_spline, stat=istat)

      END SUBROUTINE repack

!!$      SUBROUTINE half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)
!!$!-----------------------------------------------
!!$!   D u m m y   A r g u m e n t s
!!$!-----------------------------------------------
!!$      REAL(rprec), DIMENSION(1:ns_i,nuv_i), INTENT(in) ::               &
!!$     &   r1_i, ru_i, rv_i, z1_i, zu_i, zv_i
!!$!-----------------------------------------------
!!$!   L o c a l   V a r i a b l e s
!!$!-----------------------------------------------
!!$      REAL(rprec), PARAMETER :: p5 = 1.d0/2.d0, zero = 0
!!$      INTEGER  :: js, lk, js1, istat
!!$      REAL(rprec)            :: mintest, maxtest, r1, s1, t1, eps
!!$      REAL(rprec), DIMENSION(nuv_i) ::                                  &
!!$     &      r12, ru12, rv12, zu12, zv12, rs12, zs12, det
!!$!-----------------------------------------------
!!$!
!!$!     ALLOCATE METRIC ELEMENT ARRAYS
!!$!
!!$      ALLOCATE(sqrtg(ns_i,nuv_i),                                       &
!!$     &         gss(ns_i,nuv_i), gsu(ns_i,nuv_i),                      &
!!$     &         gsv(ns_i,nuv_i), guu(ns_i,nuv_i),                      &
!!$     &         guv(ns_i,nuv_i), gvv(ns_i,nuv_i),                      &
!!$     &         hss(ns_i,nuv_i), hsu(ns_i,nuv_i),                      &
!!$     &         hsv(ns_i,nuv_i), huu(ns_i,nuv_i),                      &
!!$     &         huv(ns_i,nuv_i), hvv(ns_i,nuv_i),                      &
!!$     &         stat=istat)
!!$
!!$!     COMPUTE ALL ON THE HALF MESH
!!$      IF (istat .ne. 0) STOP 'Allocation error in VMEC_HALF_MESH_METRICS'
!!$      gss(1,:) = 0;  gsu(1,:) = 0;  gsv(1,:) = 0
!!$      guu(1,:) = 0;  guv(1,:) = 0;  gvv(1,:) = 0
!!$      sqrtg(1,:) = 0
!!$
!!$      !Do inner domain
!!$      DO js = 2, ns_i
!!$         js1 = js-1
!!$         r12 = p5*(r1_i(js,:) + r1_i(js1,:))
!!$         rs12 = (r1_i(js,:) - r1_i(js1,:))*ohs_i
!!$         ru12= p5*(ru_i(js,:) + ru_i(js1,:))
!!$         rv12= p5*(rv_i(js,:) + rv_i(js1,:))
!!$         zs12 = (z1_i(js,:) - z1_i(js1,:))*ohs_i
!!$         zu12= p5*(zu_i(js,:) + zu_i(js1,:))
!!$         zv12= p5*(zv_i(js,:) + zv_i(js1,:))
!!$         guu(js,:) = ru12*ru12 + zu12*zu12
!!$         guv(js,:) = ru12*rv12 + zu12*zv12
!!$         gvv(js,:) = r12*r12 + rv12*rv12 + zv12*zv12
!!$         gsu(js,:) = rs12*ru12 + zs12*zu12
!!$         gsv(js,:) = rs12*rv12 + zs12*zv12
!!$         gss(js,:) = rs12*rs12 + zs12*zs12
!!$         sqrtg(js,:) = r12*(ru12*zs12 - rs12*zu12)
!!$      END DO
!!$
!!$      mintest = MINVAL(sqrtg(2:ns_i+1,:))
!!$      maxtest = MAXVAL(sqrtg(2:ns_i+1,:))
!!$
!!$      IF (mintest*maxtest .le. zero) STOP "Jacobian changed sign!"
!!$
!!$!     Compute upper metric elements (inverse of lower matrix: det = sqrtg**2?)
!!$!
!!$      hss(1,:) = 0;  hsu(1,:) = 0;  hsv(1,:) = 0
!!$      huu(1,:) = 0;  huv(1,:) = 0;  hvv(1,:) = 0
!!$
!!$      DO js = 2, ns_i
!!$         det(:) = one/(  gss(js,:)*guu(js,:)*gvv(js,:)   &
!!$     &                +2*gsu(js,:)*guv(js,:)*gsv(js,:)   &
!!$     &                -  gsv(js,:)*guu(js,:)*gsv(js,:)   &
!!$     &                -  gss(js,:)*guv(js,:)*guv(js,:)   &
!!$     &                -  gvv(js,:)*gsu(js,:)*gsu(js,:))
!!$         hss(js,:) = det(:)*                                            &
!!$     &              (guu(js,:)*gvv(js,:) - guv(js,:)*guv(js,:))
!!$         hsu(js,:) = det(:)*                                            &
!!$     &              (guv(js,:)*gsv(js,:) - gsu(js,:)*gvv(js,:))
!!$         hsv(js,:) = det(:)*                                            &
!!$     &              (gsu(js,:)*guv(js,:) - gsv(js,:)*guu(js,:))
!!$         huu(js,:) = det(:)*                                            &
!!$     &              (gss(js,:)*gvv(js,:) - gsv(js,:)*gsv(js,:))
!!$         huv(js,:) = det(:)*                                            &
!!$     &              (gsv(js,:)*gsu(js,:) - gss(js,:)*guv(js,:))
!!$         hvv(js,:) = det(:)*                                            &
!!$     &              (gss(js,:)*guu(js,:) - gsu(js,:)*gsu(js,:))
!!$
!!$      END DO
!!$
!!$!
!!$!     CONFIRM ACCURACY OF INVERSE
!!$!
!!$!     FIRST ROW
!!$     
!!$      maxtest = 0
!!$      eps = 0.01_dp*SQRT(EPSILON(eps))
!!$
!!$      DO js = 2, ns_i
!!$        r12 = gss(js,:)*hss(js,:) + gsu(js,:)*hsu(js,:)                 &
!!$     &       +gsv(js,:)*hsv(js,:)
!!$        ru12 = gss(js,:)*hsu(js,:) + gsu(js,:)*huu(js,:)                &
!!$     &       +gsv(js,:)*huv(js,:) 
!!$        rs12 = gss(js,:)*hsv(js,:) + gsu(js,:)*huv(js,:)                &
!!$     &       +gsv(js,:)*hvv(js,:) 
!!$        r1 = MAXVAL(ABS(r12-1))
!!$        s1 = MAXVAL(ABS(ru12))
!!$        t1 = MAXVAL(ABS(rs12))
!!$        maxtest = MAX(maxtest,r1,s1,t1)
!!$      END DO
!!$    
!!$      IF (maxtest .gt. eps) then
!!$         write(*,*) "Error1 in metric elements",maxtest,r1,s1,t1
!!$         stop
!!$      ENDIF
!!$
!!$!     SECOND ROW
!!$     
!!$      maxtest = 0
!!$
!!$      DO js = 2, ns_i
!!$        r12 = gsu(js,:)*hss(js,:) + guu(js,:)*hsu(js,:)                 &
!!$     &       +guv(js,:)*hsv(js,:)
!!$        ru12 = gsu(js,:)*hsu(js,:) + guu(js,:)*huu(js,:)                &
!!$     &       +guv(js,:)*huv(js,:) 
!!$        rs12 = gsu(js,:)*hsv(js,:) + guu(js,:)*huv(js,:)                &
!!$     &       +guv(js,:)*hvv(js,:) 
!!$        r1 = MAXVAL(ABS(r12))
!!$        s1 = MAXVAL(ABS(ru12-1))
!!$        t1 = MAXVAL(ABS(rs12))
!!$        maxtest = MAX(maxtest,r1,s1,t1)
!!$      END DO
!!$
!!$      IF (maxtest .gt. eps) then
!!$         write(*,*) "Error2 in metric elements",maxtest,r1,s1,t1
!!$         stop
!!$      ENDIF
!!$
!!$!     THIRD ROW
!!$     
!!$      maxtest = 0
!!$
!!$      DO js = 2, ns_i
!!$        r12 = gsv(js,:)*hss(js,:) + guv(js,:)*hsu(js,:)                 &
!!$     &       +gvv(js,:)*hsv(js,:)
!!$        ru12 = gsv(js,:)*hsu(js,:) + guv(js,:)*huu(js,:)                &
!!$     &       +gvv(js,:)*huv(js,:) 
!!$        rs12 = gsv(js,:)*hsv(js,:) + guv(js,:)*huv(js,:)                &
!!$     &       +gvv(js,:)*hvv(js,:) 
!!$        r1 = MAXVAL(ABS(r12))
!!$        s1 = MAXVAL(ABS(ru12))
!!$        t1 = MAXVAL(ABS(rs12-1))
!!$        maxtest = MAX(maxtest,r1,s1,t1)
!!$      END DO
!!$
!!$      IF (maxtest .gt. eps) then
!!$         write(*,*) "Error3 in metric elements",maxtest,r1,s1,t1
!!$         stop
!!$      ENDIF
!!$
!!$      END SUBROUTINE half_mesh_metrics

      SUBROUTINE vmec_half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)

      use oned_int, det1 => det   !Added by L. Chacon 6/19/17
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(1:ns_i,nu_i,nv_i), INTENT(in) ::               &
     &   r1_i, ru_i, rv_i, z1_i, zu_i, zv_i
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 1.d0/2.d0, zero = 0
      INTEGER :: js, us, vs, lk, js1, istat
      REAL(rprec) :: mintest, maxtest, r1, s1, t1, eps
      REAL(rprec), DIMENSION(nu_i,nv_i) ::                                  &
     &      r12, ru12, rv12, zu12, zv12, rs12, zs12, det                    &
     &     ,guu0,guv0,gvv0,gsu0,gsv0,gss0,sqrtg0                        
!!$     &     ,r12_0,ru12_0,rv12_0,zu12_0,zv12_0,rs12_0,zs12_0
      REAL(rprec),dimension(1:ns_i) :: xs

!-----------------------------------------------
!
!     ALLOCATE METRIC ELEMENT ARRAYS
!
      ALLOCATE(sqrtg(ns_i+1,nu_i,nv_i),                                       &
     &         gss(ns_i+1,nu_i,nv_i), gsu(ns_i+1,nu_i,nv_i),                  &
     &         gsv(ns_i+1,nu_i,nv_i), guu(ns_i+1,nu_i,nv_i),                  &
     &         guv(ns_i+1,nu_i,nv_i), gvv(ns_i+1,nu_i,nv_i),                  &
     &         hss(ns_i+1,nu_i,nv_i), hsu(ns_i+1,nu_i,nv_i),                  &
     &         hsv(ns_i+1,nu_i,nv_i), huu(ns_i+1,nu_i,nv_i),                  &
     &         huv(ns_i+1,nu_i,nv_i), hvv(ns_i+1,nu_i,nv_i),                  &
     &         stat=istat)

!     COMPUTE ALL ON THE HALF MESH
      IF (istat .ne. 0) STOP 'Allocation error in VMEC_HALF_MESH_METRICS'
!!$ LC 02/01/2007
!!$      gss(1,:) = 0;  gsu(1,:) = 0;  gsv(1,:) = 0
!!$      guu(1,:) = 0;  guv(1,:) = 0;  gvv(1,:) = 0
!!$      sqrtg(1,:) = 0

      !Do inner domain
      DO js = 2, ns_i
         js1 = js-1
         r12 = p5*(r1_i(js,:,:) + r1_i(js1,:,:))
         rs12 = (r1_i(js,:,:) - r1_i(js1,:,:))*ohs_i
         ru12= p5*(ru_i(js,:,:) + ru_i(js1,:,:))
         rv12= p5*(rv_i(js,:,:) + rv_i(js1,:,:))
         zs12 = (z1_i(js,:,:) - z1_i(js1,:,:))*ohs_i
         zu12= p5*(zu_i(js,:,:) + zu_i(js1,:,:))
         zv12= p5*(zv_i(js,:,:) + zv_i(js1,:,:))
         guu(js,:,:) = ru12*ru12 + zu12*zu12
         guv(js,:,:) = ru12*rv12 + zu12*zv12
         gvv(js,:,:) = r12*r12 + rv12*rv12 + zv12*zv12
         gsu(js,:,:) = rs12*ru12 + zs12*zu12
         gsv(js,:,:) = rs12*rv12 + zs12*zv12
         gss(js,:,:) = rs12*rs12 + zs12*zs12
!!$         sqrtg(js,:,:) = r12*(ru12*zs12 - rs12*zu12)
         sqrtg(js,:,:) = sqrt((  gss(js,:,:)*guu(js,:,:)*gvv(js,:,:)   &
     &                        +2*gsu(js,:,:)*guv(js,:,:)*gsv(js,:,:)   &
     &                        -  gsv(js,:,:)*guu(js,:,:)*gsv(js,:,:)   &
     &                        -  gss(js,:,:)*guv(js,:,:)*guv(js,:,:)   &
     &                        -  gvv(js,:,:)*gsu(js,:,:)*gsu(js,:,:)))
      END DO

      !Extrapolate (first-order) at singular point (L. Chacon, 02/01/2007)
      js = 1   !SINGULAR POINT

!!$      !Find derivatives at boundary face (dx/du == 0d0 at SP)
!!$      r12_0  = r1_i(js,:,:)
!!$      rs12_0 = (r1_i(js+1,:,:) - r1_i(js,:,:))*ohs_i
!!$      ru12_0 = ru_i(js,:,:)
!!$      rv12_0 = rv_i(js,:,:)
!!$      zs12_0 = (z1_i(js+1,:,:) - z1_i(js,:,:))*ohs_i
!!$      zu12_0 = zu_i(js,:,:)
!!$      zv12_0 = zv_i(js,:,:)
!!$
!!$      !Extrapolate derivatives to ghost cell 
!!$      r12 = p5*(r1_i(js+1,:,:) + r1_i(js,:,:))
!!$      rs12 = (r1_i(js+1,:,:) - r1_i(js,:,:))*ohs_i
!!$      ru12= p5*(ru_i(js+1,:,:) + ru_i(js,:,:))
!!$      rv12= p5*(rv_i(js+1,:,:) + rv_i(js,:,:))
!!$      zs12 = (z1_i(js+1,:,:) - z1_i(js,:,:))*ohs_i
!!$      zu12= p5*(zu_i(js+1,:,:) + zu_i(js,:,:))
!!$      zv12= p5*(zv_i(js+1,:,:) + zv_i(js,:,:))
!!$
!!$      r12  = 2*r12_0  - r12 
!!$      rs12 = 2*rs12_0 - rs12
!!$      ru12 = 2*ru12_0 - ru12
!!$      rv12 = 2*rv12_0 - rv12
!!$      zs12 = 2*zs12_0 - zs12
!!$      zu12 = 2*zu12_0 - zu12
!!$      zv12 = 2*zv12_0 - zv12 
!!$
!!$      !Find metrics at ghost cell
!!$      guu(js,:,:)  = ru12*ru12 + zu12*zu12
!!$      guv(js,:,:)  = ru12*rv12 + zu12*zv12
!!$      gvv(js,:,:)  = r12*r12 + rv12*rv12 + zv12*zv12
!!$      gsu(js,:,:)  = rs12*ru12 + zs12*zu12
!!$      gsv(js,:,:)  = rs12*rv12 + zs12*zv12
!!$      gss(js,:,:)  = rs12*rs12 + zs12*zs12

      !Find derivatives at boundary face (dx/du == 0d0 at SP)
      r12  = r1_i(js,:,:)
      rs12 = (r1_i(js+1,:,:) - r1_i(js,:,:))*ohs_i
      ru12 = ru_i(js,:,:)
      rv12 = rv_i(js,:,:)
      zs12 = (z1_i(js+1,:,:) - z1_i(js,:,:))*ohs_i
      zu12 = zu_i(js,:,:)
      zv12 = zv_i(js,:,:)

      !Find metrics at boundary face  (dx/du == 0d0 at SP)
      guu0  = ru12*ru12 + zu12*zu12
      guv0  = ru12*rv12 + zu12*zv12
      gvv0  = r12*r12 + rv12*rv12 + zv12*zv12
      gsu0  = rs12*ru12 + zs12*zu12
      gsv0  = rs12*rv12 + zs12*zv12
      gss0  = rs12*rs12 + zs12*zs12

!!      sqrtg0= 0d0
!!      sqrtg(js,:,:) = 2*sqrtg0 - sqrtg(js+1,:,:)

      !Extrapolate metrics to ghost cell
      guu(js,:,:) = 2*guu0 - guu(js+1,:,:)
      guv(js,:,:) = 2*guv0 - guv(js+1,:,:)
      gvv(js,:,:) = 2*gvv0 - gvv(js+1,:,:)
      gsu(js,:,:) = 2*gsu0 - gsu(js+1,:,:)
      gsv(js,:,:) = 2*gsv0 - gsv(js+1,:,:)
      gss(js,:,:) = 2*gss0 - gss(js+1,:,:)

      sqrtg(js,:,:) = -sqrt( abs(gss(js,:,:)*guu(js,:,:)*gvv(js,:,:)   &
     &                        +2*gsu(js,:,:)*guv(js,:,:)*gsv(js,:,:)   &
     &                        -  gsv(js,:,:)*guu(js,:,:)*gsv(js,:,:)   &
     &                        -  gss(js,:,:)*guv(js,:,:)*guv(js,:,:)   &
     &                        -  gvv(js,:,:)*gsu(js,:,:)*gsu(js,:,:)))

      !Extrapolate (third-order) at outer radial boundary (L. Chacon, 06/19/2017)
      do js=2,ns_i+1
        xs(js-1) = js
      enddo
      
      do vs=1,nv_i
        do us=1,nu_i
          call IntDriver1d(ns_i-1,xs(1:ns_i-1),sqrtg(2:ns_i,us,vs)  &
     &                    ,1     ,xs(ns_i)    ,sqrtg(ns_i+1:ns_i+1,us,vs),3)
        enddo
      enddo
      
      js = ns_i
      
      !Find derivatives at boundary face
      r12  = r1_i(js,:,:)
      rs12 =(r1_i(js,:,:) - r1_i(js-1,:,:))*ohs_i
      ru12 = ru_i(js,:,:)
      rv12 = rv_i(js,:,:)
      zs12 =(z1_i(js,:,:) - z1_i(js-1,:,:))*ohs_i
      zu12 = zu_i(js,:,:)
      zv12 = zv_i(js,:,:)

      !Find metrics at boundary face
      guu0 = ru12*ru12 + zu12*zu12
      guv0 = ru12*rv12 + zu12*zv12
      gvv0 = r12*r12 + rv12*rv12 + zv12*zv12
      gsu0 = rs12*ru12 + zs12*zu12
      gsv0 = rs12*rv12 + zs12*zv12
      gss0 = rs12*rs12 + zs12*zs12

!!!$      sqrtg0 = sqrt((  gss0*guu0*gvv0   &
!!!$     &              +2*gsu0*guv0*gsv0   &
!!!$     &              -  gsv0*guu0*gsv0   &
!!!$     &              -  gss0*guv0*guv0   &
!!!$     &              -  gvv0*gsu0*gsu0))
!!!$
!!!$      sqrtg(js+1,:,:) = 2*sqrtg0 - sqrtg(js,:,:)
      
      !Extrapolate metrics to ghost cell
      guu(js+1,:,:) = 2*guu0 - guu(js,:,:)
      guv(js+1,:,:) = 2*guv0 - guv(js,:,:)
      gvv(js+1,:,:) = 2*gvv0 - gvv(js,:,:)
      gsu(js+1,:,:) = 2*gsu0 - gsu(js,:,:)
      gsv(js+1,:,:) = 2*gsv0 - gsv(js,:,:)
      gss(js+1,:,:) = 2*gss0 - gss(js,:,:)

!!$      sqrtg(js+1,:,:) = sqrt((  gss(js+1,:,:)*guu(js+1,:,:)*gvv(js+1,:,:)   &
!!$     &                       +2*gsu(js+1,:,:)*guv(js+1,:,:)*gsv(js+1,:,:)   &
!!$     &                       -  gsv(js+1,:,:)*guu(js+1,:,:)*gsv(js+1,:,:)   &
!!$     &                       -  gss(js+1,:,:)*guv(js+1,:,:)*guv(js+1,:,:)   &
!!$     &                       -  gvv(js+1,:,:)*gsu(js+1,:,:)*gsu(js+1,:,:)))

      !Test
      mintest = MINVAL(sqrtg(2:ns_i+1,:,:))
      maxtest = MAXVAL(sqrtg(2:ns_i+1,:,:))

!!$      write (*,*) 'DIAG -- vmec_half_mesh_metrics;'
!!$      write (*,*) sqrtg(1:ns_i+1,1)
!!$      stop
!!$      write (*,*) 'DIAG -- half_mesh_metrics:',mintest,maxtest
      IF (mintest*maxtest .le. zero) STOP "Jacobian changed sign!"

!     Compute upper metric elements (inverse of lower matrix: det = sqrtg**2?)
!
!!$ LC 02/01/2007
!!$      hss(1,:,:) = 0;  hsu(1,:,:) = 0;  hsv(1,:,:) = 0
!!$      huu(1,:,:) = 0;  huv(1,:,:) = 0;  hvv(1,:,:) = 0
!!$
!!$      DO js = 2, ns_i
      DO js = 1, ns_i+1
         det(:,:) = one/(  gss(js,:,:)*guu(js,:,:)*gvv(js,:,:)   &
     &                 + 2*gsu(js,:,:)*guv(js,:,:)*gsv(js,:,:)   &
     &                 -   gsv(js,:,:)*guu(js,:,:)*gsv(js,:,:)   &
     &                 -   gss(js,:,:)*guv(js,:,:)*guv(js,:,:)   &
     &                 -   gvv(js,:,:)*gsu(js,:,:)*gsu(js,:,:))
         hss(js,:,:) = det(:,:)*                                            &
     &              (guu(js,:,:)*gvv(js,:,:) - guv(js,:,:)*guv(js,:,:))
         hsu(js,:,:) = det(:,:)*                                            &
     &              (guv(js,:,:)*gsv(js,:,:) - gsu(js,:,:)*gvv(js,:,:))
         hsv(js,:,:) = det(:,:)*                                            &
     &              (gsu(js,:,:)*guv(js,:,:) - gsv(js,:,:)*guu(js,:,:))
         huu(js,:,:) = det(:,:)*                                            &
     &              (gss(js,:,:)*gvv(js,:,:) - gsv(js,:,:)*gsv(js,:,:))
         huv(js,:,:) = det(:,:)*                                            &
     &              (gsv(js,:,:)*gsu(js,:,:) - gss(js,:,:)*guv(js,:,:))
         hvv(js,:,:) = det(:,:)*                                            &
     &              (gss(js,:,:)*guu(js,:,:) - gsu(js,:,:)*gsu(js,:,:))

      END DO

!
!     CONFIRM ACCURACY OF INVERSE
!
!     FIRST ROW
     
      maxtest = 0
      eps = 0.01_dp*SQRT(EPSILON(eps))

      DO js = 1, ns_i+1
        r12 = gss(js,:,:)*hss(js,:,:) + gsu(js,:,:)*hsu(js,:,:)                 &
     &       +gsv(js,:,:)*hsv(js,:,:)
        ru12 = gss(js,:,:)*hsu(js,:,:) + gsu(js,:,:)*huu(js,:,:)                &
     &       +gsv(js,:,:)*huv(js,:,:) 
        rs12 = gss(js,:,:)*hsv(js,:,:) + gsu(js,:,:)*huv(js,:,:)                &
     &       +gsv(js,:,:)*hvv(js,:,:) 
        r1 = MAXVAL(ABS(r12-1))
        s1 = MAXVAL(ABS(ru12))
        t1 = MAXVAL(ABS(rs12))
        maxtest = MAX(maxtest,r1,s1,t1)
      END DO
    
      IF (maxtest .gt. eps) then
         write(*,*) "Error1 in metric elements",maxtest,r1,s1,t1
         stop
      ENDIF

!     SECOND ROW
     
      maxtest = 0

      DO js = 1, ns_i+1
        r12 = gsu(js,:,:)*hss(js,:,:) + guu(js,:,:)*hsu(js,:,:)                 &
     &       +guv(js,:,:)*hsv(js,:,:)
        ru12 = gsu(js,:,:)*hsu(js,:,:) + guu(js,:,:)*huu(js,:,:)                &
     &       +guv(js,:,:)*huv(js,:,:) 
        rs12 = gsu(js,:,:)*hsv(js,:,:) + guu(js,:,:)*huv(js,:,:)                &
     &       +guv(js,:,:)*hvv(js,:,:) 
        r1 = MAXVAL(ABS(r12))
        s1 = MAXVAL(ABS(ru12-1))
        t1 = MAXVAL(ABS(rs12))
        maxtest = MAX(maxtest,r1,s1,t1)
      END DO

      IF (maxtest .gt. eps) then
         write(*,*) "Error2 in metric elements",maxtest,r1,s1,t1
         stop
      ENDIF

!     THIRD ROW
     
      maxtest = 0

      DO js = 1, ns_i+1
        r12 = gsv(js,:,:)*hss(js,:,:) + guv(js,:,:)*hsu(js,:,:)                 &
     &       +gvv(js,:,:)*hsv(js,:,:)
        ru12 = gsv(js,:,:)*hsu(js,:,:) + guv(js,:,:)*huu(js,:,:)                &
     &       +gvv(js,:,:)*huv(js,:,:) 
        rs12 = gsv(js,:,:)*hsv(js,:,:) + guv(js,:,:)*huv(js,:,:)                &
     &       +gvv(js,:,:)*hvv(js,:,:) 
        r1 = MAXVAL(ABS(r12))
        s1 = MAXVAL(ABS(ru12))
        t1 = MAXVAL(ABS(rs12-1))
        maxtest = MAX(maxtest,r1,s1,t1)
      END DO

      IF (maxtest .gt. eps) then
         write(*,*) "Error3 in metric elements",maxtest,r1,s1,t1
         stop
      ENDIF

      END SUBROUTINE vmec_half_mesh_metrics

      SUBROUTINE vmec_cleanup_metrics
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER         :: istat
!-----------------------------------------------

      DEALLOCATE(gss, gsu, gsv, guu, guv, gvv,                          &
     &           hss, hsu, hsv, huu, huv, hvv,                          &
     &           sqrtg, stat=istat)
	
      END SUBROUTINE vmec_cleanup_metrics

      
      SUBROUTINE FIXARRAY
!
!     This subroutine computes the cosine-sine factors that will be needed when moving between
!     Fourier and real space. All normalizations are contained in the poloidal quantities 
!     used in the Fourier to real space transformation: SINMUI, COSMUI
!
!     Fourier representations are assumed to have STELLARATOR symmetry:
!         1. COSINE (iparity=0):      C(u, v) = sum_u sum_v (C_mn COS(mu + n*nfp*v))
!         2. SINE   (iparity=1):      S(u, v) = sum_u sum_v (S_mn SIN(mu + n*nfp*v))
!
!     The number of collocation points have been set initially equal to the number of modes:
!       theta_j = j*pi/M, (j=0,...,M);   zeta_k = k*2*pi/(2N+1), k = 0,...., 2N
!       where M = mpol - 1 and N = ntor.
!
        USE stel_kinds
        USE Fourier, ONLY: orthonorm
        IMPLICIT NONE
        INTEGER(iprec) :: jk, m, n, lk, istat, ntheta, nzeta, ntheta1
        INTEGER(iprec) :: mpol, ntor
        REAL(rprec)    :: arg, cosi, sini
        
        ntheta  = nu_i; nzeta   = nv_i
        ntheta1 = ntheta - 1
        mpol = mpol_i;  ntor = ntor_i

        ALLOCATE(cosmu(0:mpol,1:ntheta), cosmum(0:mpol,1:ntheta),      &
     &    cosmui(0:mpol,1:ntheta), sinmu(0:mpol,1:ntheta),             &
     &    sinmum(0:mpol,1:ntheta), sinmui(0:mpol,1:ntheta),            &
     &    cosnv(0:ntor,1:nzeta), cosnvn(0:ntor,1:nzeta),               &
     &    sinnv(0:ntor,1:nzeta), sinnvn(0:ntor,1:nzeta), stat=istat)

        IF (istat .ne. 0) STOP 'Allocation error in fixarray'
        
        dnorm_i = 1.d0/((ntheta-1)*nzeta)                         ! Norm for Fourier-Space -> UV-Space
        DO jk = 1, ntheta
          DO m = 0, mpol
            arg = m*pi*(jk - 1)/ntheta1                           ! Poloidal collocation points*m
            cosi = COS(arg)
            sini = SIN(arg)
!            IF (m .ne. 0) THEN
!              cosi = cosi*nfactor; sini = sini*nfactor            ! RS (10/03/06): Added n=0/m=0 normalization
!            ENDIF
!            IF (m .eq. 0) THEN
!              cosi = cosi/nfactor; sini = sini/nfactor            ! RS (10/03/06): Added n=0/m=0 normalization
!            ENDIF
            cosmu(m, jk)  = cosi                                  ! Used in UV-SPace to Fourier-space
            cosmum(m, jk)  = m*cosi                               ! Used in UV-SPace to Fourier-space
            cosmui(m, jk) = dnorm_i*cosi                          ! Used in Fourier-Space to UV-space
            sinmu(m, jk)  = sini                                  ! Used in UV-SPace to Fourier-space
            sinmum(m, jk)  = m*sini                               ! Used in UV-SPace to Fourier-space
            sinmui(m, jk) = dnorm_i*sini                          ! Used in Fourier-Space to UV-space
            IF (jk == 1 .or. jk == ntheta) then
              cosmui(m, jk) = 0.5d0*cosmui(m, jk)                 ! 1st and last poloidal points are divided by 2!
              sinmui(m, jk) = 0.5d0*sinmui(m, jk)
            ENDIF
          ENDDO
        ENDDO

!
!FOR NOW, WE ARE IN 1 FIELD PERIOD
        DO lk = 1, nzeta
          DO n = 0, ntor
            arg = n*twopi*(lk - 1)/nzeta                          ! Toroidal collocation points*ABS(n*nfp_i); 
            cosi = COS(arg)                                       ! sign of "n" taken care later through ISIGN
            sini = SIN(arg)
!            IF (n .ne. 0) THEN
!              cosi = cosi*nfactor; sini = sini*nfactor           ! RS (10/03/06): Added n=0/m=0 normalization
!            ENDIF
!            IF (n .eq. 0) THEN
!              cosi = cosi/nfactor; sini = sini/nfactor           ! RS (10/03/06): Added n=0/m=0 normalization
!            ENDIF
            cosnv(n, lk)  = cosi                                  ! All normalization factors goes in poloidal quantities
            cosnvn(n, lk)  = n*nfp_i*cosi                         ! NFP is number of field periods
            sinnv(n, lk)  = sini
            sinnvn(n, lk)  = n*nfp_i*sini
          ENDDO
        ENDDO

!COMPUTE orthonorm FACTOR FOR cos(mu+nv), sin(mu+nv) BASIS
        ALLOCATE (orthonorm(0:mpol,-ntor:ntor), stat=istat)
        IF (istat .ne. 0) STOP 'Allocate orthonorm failed in FIXARRAY!'
        orthonorm = nfactor
        orthonorm(0,0) = 1

        END SUBROUTINE FIXARRAY


        SUBROUTINE DEALLOC_FIXARRAY
        USE stel_kinds
        IMPLICIT NONE
        INTEGER(kind=iprec):: istat

        DEALLOCATE(cosmu, cosmum, cosmui, sinmu, sinmum, sinmui,            &
     &    cosnv, cosnvn, sinnv, sinnvn, stat=istat) 
        IF (istat .ne. 0) STOP 'Allocation error in dealloc_fixarray'
         
        END SUBROUTINE DEALLOC_FIXARRAY

!       half_to_int
!       ###########################################################################
        subroutine half_to_int(igl,jgl,kgl,vmec_arr,local_val,order)

!       ---------------------------------------------------------------------------
!       Averages quantities from VMEC's full radial mesh (which is
!       PIXIE's half radial mesh) to PIXIE's collocated mesh. At ghost
!       cells, it extrapolates using a first-order formula.
!       ---------------------------------------------------------------------------

          use oned_int

!       Call variables

          integer :: igl,jgl,kgl
          integer,optional :: order

          real(rprec) :: vmec_arr(ns_i,nu_i,nv_i),local_val

!       Local variables

          integer :: ordr
          real(8) :: pos(5),val(5)

!       Begin program

          if (PRESENT(order)) then
             ordr = order
          else
             ordr = 2
          endif

          pos(1:4) = (/ 0d0,1d0,2d0,3d0 /)

          if (igl == 0) then
             pos(5) = -0.5d0
             val(1:4) = vmec_arr(igl+1:igl+4,jgl,kgl)
             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
          elseif (igl == 1) then
             pos(5) = 0.5d0
             val(1:4) = vmec_arr(igl:igl+3,jgl,kgl)
             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
          elseif (igl == ns_i) then
             pos(5) = 3.5d0
             val(1:4) = vmec_arr(igl-3:igl,jgl,kgl)
             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
          elseif (igl == ns_i - 1) then
             pos(5) = 2.5d0
             val(1:4) = vmec_arr(igl-2:igl+1,jgl,kgl)
             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
          else
             pos(5) = 1.5d0
             val(1:4) = vmec_arr(igl-1:igl+2,jgl,kgl)
             call IntDriver1d(4,pos(1:4),val(1:4),1,pos(5),val(5),ordr)
          endif

          local_val = val(5)

        end subroutine half_to_int

      END MODULE vmec_mod

!     vmec_map
!     #################################################################
      subroutine vmec_map(equ_file,g_def)

!     -----------------------------------------------------------------
!     Give Cartesian coordinates of each logical mesh point at grid
!     level (igrid).
!     -----------------------------------------------------------------

      use vmec_mod, pi_stel => pi
      use grid
      use equilibrium

      implicit none

!     Input variables

      type(grid_mg_def),pointer :: g_def
!!$      logical :: metrics
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

!!$      load_metrics = metrics  !Set flag for metric elements routine

!     Cycle grid levels

!!$      igrid = 1
      do igrid=1,g_def%ngrid

!     Get LOCAL limits and allocate local map array

        nx = g_def%nxv(igrid)
        ny = g_def%nyv(igrid)
        nz = g_def%nzv(igrid)

        allocate(xcar(0:nx+1,0:ny+1,0:nz+1,3))

!     Get GLOBAL limits (VMEC operates on global domain)

        nxg = g_def%nxgl(igrid)
        nyg = g_def%nygl(igrid)
        nzg = g_def%nzgl(igrid)

        if (my_rank == 0) then
           write (*,'(a,i3,a,i3,a,i3,a,i3)') &
          ' Reading VMEC map on grid',igrid,', nx x ny x nz=',nxg,'x',nyg,'x',nzg
        endif

!     Read equilibrium file and setup arrays
!     [VMEC++ assumes solution is up-down symmetric wrt Z=0, 
!     and hence only gives theta=(0,pi); thereby the limit nyg/2+1 in theta]

        call vmec_init(nxg+1,nyg/2+1,nzg,equ_file,.false.)

!     Spline VMEC global map

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

        if( mod(nzg,nfp_i) /= 0) then
          if (my_rank == 0) write (*,*) 'Number of VMEC periods',nfp_i
          call pstop('vmec_map','Toroidal mesh cannot accommodate # periods')
        endif

        dphi = 2*pi/nzg/nfp_i  !VMEC stores 1 of NFP periods in phi, ie., v=2*pi/nfp
        do kg = 1,nzs
          zs(kg) = dphi*(kg-2)
        enddo

        !Spline RR, ZZ
        allocate(rr_coef(nxs,nys,nzs),zz_coef(nxs,nys,nzs),stat=alloc_stat)

        call db3ink(xs,nxs,ys,nys,zs,nzs,rr,nxs,nys,kx,ky,kz,tx,ty,tz,rr_coef,work,flg)
        call db3ink(xs,nxs,ys,nys,zs,nzs,zz,nxs,nys,kx,ky,kz,tx,ty,tz,zz_coef,work,flg)

!     Transfer map (from GLOBAL in VMEC to LOCAL)

!!$        ph0 = g_def%zg(1)  !Reference phi=0 plane at phi(k=1)
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
              v1  = mod(ph1,2*pi/nfp_i)

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

        call vmec_cleanup

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

      end subroutine vmec_map

!!$!     vmec_map
!!$!     #################################################################
!!$      subroutine vmec_map(metrics)
!!$
!!$!     -----------------------------------------------------------------
!!$!     Give Cartesian coordinates of each logical mesh point at grid
!!$!     level (igrid).
!!$!     -----------------------------------------------------------------
!!$
!!$      use vmec_mod, pi_stel => pi
!!$      use grid
!!$      use equilibrium
!!$
!!$      implicit none
!!$
!!$!     Input variables
!!$
!!$      logical :: metrics
!!$
!!$!     Local variables
!!$
!!$      integer    :: igrid,nx,ny,nz,alloc_stat,icomp,nxs
!!$      integer    :: nxg,nyg,nzg,i,j,k,igl,jgl,kgl,ig,jg,kg
!!$      real(8)    :: r1,r2,r3,th1,v1,ph1,sgn,ds,dth,dphi,rr1(1),zz1(1)
!!$
!!$      real(8),dimension(:),allocatable :: xs
!!$
!!$      real(8),allocatable,dimension(:,:,:,:) :: xcar
!!$
!!$!     Begin program
!!$
!!$      load_metrics = metrics  !Set flag for metric elements routine
!!$
!!$      nullify(gmetric)
!!$
!!$!     Cycle grid levels
!!$
!!$      do igrid=1,gv%gparams%ngrid
!!$
!!$        if (my_rank == 0) then
!!$           write (*,'(a,i3)') ' Reading VMEC map on grid',igrid
!!$        endif
!!$
!!$!     Get LOCAL limits and allocate local map array
!!$
!!$        nx = gv%gparams%nxv(igrid)
!!$        ny = gv%gparams%nyv(igrid)
!!$        nz = gv%gparams%nzv(igrid)
!!$
!!$        allocate(xcar(0:nx+1,0:ny+1,0:nz+1,3))
!!$
!!$!     Get GLOBAL limits (VMEC operates on global domain)
!!$
!!$        nxg = gv%gparams%nxgl(igrid)
!!$        nyg = gv%gparams%nygl(igrid)
!!$        nzg = gv%gparams%nzgl(igrid)
!!$
!!$        nxs = nxg+1
!!$
!!$        allocate(xs(nxs),stat=alloc_stat)
!!$
!!$        ds = 1d0/nxg
!!$        do ig = 1,nxs
!!$           xs(ig) = ds*(ig-1)
!!$        enddo
!!$
!!$!     Read equilibrium file and setup arrays
!!$!     [VMEC++ assumes solution is up-down symmetric wrt Z=0, 
!!$!     and hence only gives theta=(0,pi); thereby the limit nyg/2+1 in theta]
!!$
!!$        call vmec_init(nxg+1,nyg/2+1,nzg,equ_file)
!!$
!!$!     Transfer map (from GLOBAL in VMEC to LOCAL)
!!$
!!$        do k=0,nz+1
!!$          do j=0,ny+1
!!$            do i=0,nx+1
!!$
!!$              call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)
!!$
!!$              !Find local coordinates
!!$              r1  = gv%gparams%xx(ig)
!!$              th1 = gv%gparams%yy(jg)
!!$              ph1 = gv%gparams%zz(kg)
!!$
!!$              !Find global limits
!!$              call fromLocalToGlobalLimits(i,j,k,igl,jgl,kgl,igrid,igrid,igrid)
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
!!$              !Perform interpolation (r,th,phi) -> (R,Z)
!!$              call IntDriver1d(nxs,xs,rr(:,jgl,kgl),1,(/r1/),rr1,3)
!!$              call IntDriver1d(nxs,xs,zz(:,jgl,kgl),1,(/r1/),zz1,3)
!!$
!!$              !Transform to Cartesian geometry (minus sign in phi to preserve a right-handed ref. sys.)
!!$              xcar(i,j,k,1)=rr1(1)*cos(-ph1)
!!$              xcar(i,j,k,2)=rr1(1)*sin(-ph1)
!!$              xcar(i,j,k,3)=sgn*zz1(1)
!!$            enddo
!!$          enddo
!!$        enddo
!!$
!!$!     Fill grid metrics hierarchy
!!$
!!$        call defineGridMetric(gv%gparams,xcar=xcar,igr=igrid)
!!$
!!$!     Free work space (to allow processing of different grid levels)
!!$
!!$        deallocate(xcar,xs)
!!$
!!$        call vmec_cleanup
!!$
!!$      enddo
!!$
!!$!     Set up gmetric pointer
!!$
!!$      gmetric => gv%gparams%gmetric
!!$
!!$!     End program
!!$
!!$      end subroutine vmec_map

      !THIS ROUTINE IS INCORRECT, AS IT DOES NOT TAKE INTO ACCOUNT APPROPRIATELY
      !THE CHANGE OF RADIAL COORDINATES BETWEEN VMEC AND PIXIE3D (PSI -> SQRT(PSI))
!!$!     vmec_metrics
!!$!     #################################################################
!!$      subroutine vmec_metrics(igrid,nx,ny,nz,jac,gsub,gsup)
!!$
!!$!     -----------------------------------------------------------------
!!$!     Give metrics at each logical mesh point in grid level (igrid).
!!$!     -----------------------------------------------------------------
!!$
!!$        use vmec_mod
!!$        use grid
!!$
!!$        implicit none
!!$
!!$!     Input variables
!!$
!!$        integer    :: igrid,nx,ny,nz
!!$        real(8)    :: jac (0:nx+1,0:ny+1,0:nz+1)            &
!!$     &               ,gsub(0:nx+1,0:ny+1,0:nz+1,3,3)        &
!!$     &               ,gsup(0:nx+1,0:ny+1,0:nz+1,3,3)
!!$
!!$!     Local variables
!!$
!!$        integer    :: nxg,nyg,nzg,i,j,k,igl,jgl,kgl,ig,jg,kg,jkg,m,l
!!$        real(8)    :: sgn,ds
!!$
!!$!     Begin program
!!$
!!$!     Get GLOBAL limits (VMEC operates on global domain)
!!$
!!$        nxg = gv%gparams%nxgl(igrid)
!!$        nyg = gv%gparams%nygl(igrid)
!!$        nzg = gv%gparams%nzgl(igrid)
!!$
!!$!     Transfer metrics (from GLOBAL in VMEC to LOCAL)
!!$
!!$        do k=0,nz+1
!!$          do j=0,ny+1
!!$            do i=0,nx+1
!!$
!!$              call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)
!!$
!!$              !Find global limits
!!$              call fromLocalToGlobalLimits(i,j,k,igl,jgl,kgl,igrid,igrid,igrid)
!!$
!!$              !Periodic boundary in theta (other boundary enforced by symmetry)
!!$              if (jgl == 0) jgl = nyg
!!$
!!$              !Up-down symmetry in theta (temporary, until VMEC interface is fixed)
!!$              sgn = 1d0
!!$              if (jgl > nyg/2+1) then
!!$                 jgl = nyg + 2 - jgl
!!$                 sgn = -1d0
!!$              endif
!!$
!!$              !Periodic boundary in phi
!!$              if (kgl == 0) kgl = nzg
!!$              if (kgl == nzg+1) kgl = 1
!!$
!!$              !Assign metric elements
!!$              jac(i,j,k)      = sqrtg(igl+1,jgl,kgl)
!!$
!!$              gsub(i,j,k,1,1) = gss(igl+1,jgl,kgl)
!!$              gsub(i,j,k,1,2) = gsu(igl+1,jgl,kgl)
!!$              gsub(i,j,k,1,3) = gsv(igl+1,jgl,kgl)
!!$              gsub(i,j,k,2,2) = guu(igl+1,jgl,kgl)
!!$              gsub(i,j,k,2,3) = guv(igl+1,jgl,kgl)
!!$              gsub(i,j,k,3,3) = gvv(igl+1,jgl,kgl)
!!$
!!$              gsup(i,j,k,1,1) = hss(igl+1,jgl,kgl)
!!$              gsup(i,j,k,1,2) = hsu(igl+1,jgl,kgl)
!!$              gsup(i,j,k,1,3) = hsv(igl+1,jgl,kgl)
!!$              gsup(i,j,k,2,2) = huu(igl+1,jgl,kgl)
!!$              gsup(i,j,k,2,3) = huv(igl+1,jgl,kgl)
!!$              gsup(i,j,k,3,3) = hvv(igl+1,jgl,kgl)
!!$
!!$              !Enforce symmetry
!!$              do m=1,3
!!$                 do l=1,m
!!$                  gsub(i,j,k,m,l) = gsub(i,j,k,l,m)
!!$                  gsup(i,j,k,m,l) = gsup(i,j,k,l,m)
!!$                enddo
!!$              enddo
!!$
!!$              !Transform to PIXIE3D's convention
!!$              gsup(i,j,k,:,:) = gsup(i,j,k,:,:)*jac(i,j,k)
!!$              gsub(i,j,k,:,:) = gsub(i,j,k,:,:)/jac(i,j,k)
!!$            enddo
!!$          enddo
!!$        enddo
!!$
!!$!     Free work space (to allow multiple calls to vmec_map for different grid levels)
!!$
!!$        call vmec_cleanup_metrics
!!$
!!$!     End program
!!$
!!$      end subroutine vmec_metrics

!     vmec_equ
!     #################################################################
      subroutine vmec_equ(iout,igrid,nx,ny,nz,bb,prs,rho,gam,equ_file &
     &                   ,dcon,divcl)

!     -----------------------------------------------------------------
!     Give equilibrium fields at each logical mesh point in grid
!     level (igrid).
!     -----------------------------------------------------------------

        use vmec_mod, pi_stel => pi
        use grid
        use equilibrium

        implicit none

!      Call variables

        integer :: igrid,iout,nx,ny,nz
        real(8) :: gam
        real(8),dimension(0:nx+1,0:ny+1,0:nz+1,3) :: bb
        real(8),dimension(0:nx+1,0:ny+1,0:nz+1)   :: prs,rho
        character(*) :: equ_file

        logical :: dcon,divcl

!     Local variables

        integer :: nxg,nyg,nzg,i,j,k,igl,jgl,kgl,ig,jg,kg,istat,its
        real(8) :: r1,z1,th1,ph1,v1,dphi,dth,sgn,jac,ds,dum1,dum2,ppi,max_rho

        real(8),allocatable, dimension(:)     :: pflx,ff_i,q_i
        real(8),allocatable, dimension(:,:,:) :: bsub3,prsg
        real(8),allocatable, dimension(:,:,:,:) :: bsup

        integer :: udcon=1111

        logical :: enf_tor_flx_fn  !Whether we enforce flux functions in toroidal geom.

        !SLATEC spline variables
        integer :: kx,ky,kz,nxs,nys,nzs,dim,flg,sorder,alloc_stat
        real(8),dimension(:)    ,allocatable :: tx,ty,tz,work,xs,ys,zs
        real(8),dimension(:,:,:),allocatable :: b1_coef,b2_coef,b3_coef,prs_coef

        real(8)  :: db3val
        external :: db3val

!     Begin program

        if (my_rank == 0) then
           write (*,*)
           write (*,*) 'Reading VMEC solution...'
        endif

!     Get GLOBAL limits (VMEC operates on global domain)

!!$        nx  = gv%gparams%nxv(igrid)
!!$        ny  = gv%gparams%nyv(igrid)
!!$        nz  = gv%gparams%nzv(igrid)

        nxg = gv%gparams%nxgl(igrid)
        nyg = gv%gparams%nygl(igrid)
        nzg = gv%gparams%nzgl(igrid)

!     Read equilibrium file and setup arrays
!     [Assumes solution is up-down symmetric wrt Z=0, 
!     and hence only gives theta=(0,pi); thereby the limit nyg/2+1 in theta]

!!$        load_metrics = .true.  !Whether to use VMEC metrics to define B components

        call vmec_init(nxg+1,nyg/2+1,nzg,equ_file,.true.)

!     Setup equilibrium fields (at VMEC integer mesh -- PIXIE3D's half mesh)

        call vmec_init_fields

!     DCON dump

        if (dcon .and. my_rank == 0) then

          write (*,*)
          write (*,*) 'Dumping DCON file with VMEC solution...'

          allocate(pflx(ns_i),ff_i(ns_i),q_i(ns_i))

          !VMEC POLOIDAL flux surface positions (integral of 2*pi*jac*B^2)
          ds = 1d0/(ns_i-1)
          ppi = acos(-1d0)
          pflx(1) = 0d0
          do i = 2,ns_i
            pflx(i) = pflx(i-1) + ppi*ds*sqrtg(i,1,1)*(bsupuijcf(i,1,1)+bsupuijcf(i-1,1,1))
          enddo

          !F factor (R*Bt=B_3)
          ff_i = bsubvijcf(:,1,1)

          !q-profile
          q_i = 1d0/iotaf_i

          !DCON dump
          open(unit=udcon,file='pixie-dcon.bin',form='unformatted',status='unknown')

          !Poloidal plane size
          write (udcon) nxg,nyg

          write (udcon) pflx      !Poloidal flux
          write (udcon) ff_i      !R*B_t (flux function)
          write (udcon) presf_i   !Pressure
          write (udcon) q_i       !Q-profile

          !R(psi,theta), Z(psi,theta)

          k = 1  !Fix poloidal plane
          do j=1,nyg/2+1   !Cycle in poloidal angle (only half-plane)
            write (udcon) rr(:,j,k)
            write (udcon) zz(:,j,k)
          enddo

          close (udcon)

          deallocate(pflx,ff_i,q_i)
        endif

!     Find half-mesh PIXIE3D magnetic field components in GLOBAL mesh (flux functions)

        allocate(bsup (0:nxg+1,0:nyg+1,0:nzg+1,3))
        allocate(bsub3(0:nxg+1,0:nyg+1,0:nzg+1))
        allocate(prsg (0:nxg+1,0:nyg+1,0:nzg+1))

        do k=0,nzg+1
          do j=0,nyg+1
            do i=0,nxg+1
               igl = i
               jgl = j
               kgl = k

              !Periodic boundary in theta=0 (other boundary enforced by symmetry)
              if (jgl == 0) jgl = nyg

              !Periodic boundary in phi=0
              if (kgl == 0) kgl = nzg

              !Up-down symmetry in theta, phi: R(th,phi) = R(2pi-th,2pi-phi)
              !                                Z(th,phi) =-Z(2pi-th,2pi-phi)
              sgn = 1d0
              if (jgl > nyg/2+1) then
                 jgl = nyg + 2 - jgl

                 kgl = nzg + 2 - kgl
                 sgn = -1d0
              endif

              !Periodic boundary in phi=2*pi
              if (kgl == nzg+1) kgl = 1

              !Interpolate VMEC fields to PIXIE3D's integer radial mesh
              call half_to_int(igl,jgl,kgl,bsupsijsf,bsup(i,j,k,1))
              call half_to_int(igl,jgl,kgl,bsupuijcf,bsup(i,j,k,2))
              call half_to_int(igl,jgl,kgl,bsupvijcf,bsup(i,j,k,3))
              call half_to_int(igl,jgl,kgl,bsubvijcf,bsub3(i,j,k))
              call half_to_int(igl,jgl,kgl,presijf  ,prsg (i,j,k))

              !Find flux functions w/ PIXIE3D convention
              jac = sqrtg(igl+1,jgl,kgl) !VMEC's half-mesh jacobian, including s=sqrt(psi) transformation
              bsup(i,j,k,1) = jac*bsup(i,j,k,1)/(2.*gv%gparams%xg(igl)) !B1=B1'/(2s)
              bsup(i,j,k,2) = jac*bsup(i,j,k,2)!Flux coordinate, B2=B2'
              bsup(i,j,k,3) = jac*bsup(i,j,k,3)!Flux coordinate, B3=B3'
            enddo
          enddo
        enddo

!!$          write (*,*) 'jac',sqrtg(:,1,1)
!!$          write (*,*) 'b1',bsup(:,1,1,1)
!!$          write (*,*) 'b2',bsup(:,1,1,2)
!!$          write (*,*) 'b3',bsup(:,1,1,3)
!!$          stop

        !Calculate maximum density (to normalize density later)
        max_rho = maxval(prsg(1:nxg,1:nyg,1:nzg)**(1d0/gam))

        !Clean flux functions (Tokamak case)
        enf_tor_flx_fn = (nfp_i == 1)

        divcl = .not.enf_tor_flx_fn !Whether to divergence clean
        
        if (enf_tor_flx_fn) then
          do k=0,nzg+1
            do i=0,nxg+1
               igl = i
               jgl = j
               kgl = k

              !Periodic boundary in theta=0 (other boundary enforced by symmetry)
              if (jgl == 0) jgl = nyg

              !Periodic boundary in phi=0
              if (kgl == 0) kgl = nzg

              !Up-down symmetry in theta, phi: R(th,phi) = R(2pi-th,2pi-phi)
              !                                Z(th,phi) =-Z(2pi-th,2pi-phi)
              sgn = 1d0
              if (jgl > nyg/2+1) then
                 jgl = nyg + 2 - jgl

                 kgl = nzg + 2 - kgl
                 sgn = -1d0
              endif

              !Periodic boundary in phi=2*pi
              if (kgl == nzg+1) kgl = 1

              dum1=sum(bsup (i,1:nyg,k,2))/nyg
              dum2=sum(bsub3(i,1:nyg,k))/nyg

              bsup (i,:,k,2) = dum1
              bsub3(i,:,k)   = dum2
             
            enddo
          enddo
        endif

!     Spline PIXIE3D global variables

        !Prepare splines
        sorder = 2
        nxs = nxg+2
        nys = nyg+2
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

        !Global coordinates 
        allocate(xs(nxs),ys(nys),zs(nzs),stat=alloc_stat)

        !Radial half-mesh (PIXIE3D's full mesh)
        ds = 1d0/nxg
        do ig = 1,nxs
!!           xs(ig) = ds*(ig-1)
           xs(ig) = ds*(ig-1.5)
        enddo

        !Full angle meshes, starting at angle 0
        dth = 2*pi/nyg
        do jg = 1,nys
           ys(jg) = dth*(jg-1)
        enddo

        dphi = 2*pi/nzg/nfp_i  !VMEC stores 1 of NFP periods in phi, ie., v=2*pi/nfp
        do kg = 1,nzs
           zs(kg) = dphi*(kg-2)
        enddo

!!$        xs = gv%gparams%xg(0:nxg+1)
!!$        ys = gv%gparams%yg(0:nyg+1)
!!$        zs = gv%gparams%zg(0:nzg+1)

        !Spline B-field components
        allocate(b1_coef (nxs,nys,nzs)  &
     &          ,b2_coef (nxs,nys,nzs)  &
     &          ,b3_coef (nxs,nys,nzs)  &
     &          ,prs_coef(nxs,nys,nzs),stat=alloc_stat)

        call db3ink(xs,nxs,ys,nys,zs,nzs,bsup(:,:,:,1),nxs,nys,kx,ky,kz,tx,ty,tz,b1_coef,work,flg)
        call db3ink(xs,nxs,ys,nys,zs,nzs,bsup(:,:,:,2),nxs,nys,kx,ky,kz,tx,ty,tz,b2_coef,work,flg)
        if (enf_tor_flx_fn) then
          call db3ink(xs,nxs,ys,nys,zs,nzs,bsub3,nxs,nys,kx,ky,kz,tx,ty,tz,b3_coef,work,flg)
        else
          call db3ink(xs,nxs,ys,nys,zs,nzs,bsup(:,:,:,3),nxs,nys,kx,ky,kz,tx,ty,tz,b3_coef,work,flg)
        endif
        call db3ink(xs,nxs,ys,nys,zs,nzs,prsg ,nxs,nys,kx,ky,kz,tx,ty,tz,prs_coef ,work,flg)

!     Transfer variables (from GLOBAL in VMEC to LOCAL in PIXIE3D)

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1

              call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)

              !Find global limits
              call fromLocalToGlobalLimits(gv%gparams,igrid,i,j,k,igl,jgl,kgl)

              !Find local coordinates
              r1  = gv%gparams%xx(ig)
              th1 = gv%gparams%yy(jg)
              ph1 = gv%gparams%zz(kg)

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

              !Map phi to VMEC toroidal coordinate, v
              v1  = mod(ph1,2*pi/nfp_i)

              bb(i,j,k,1) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                          ,kx,ky,kz,b1_coef,work)

              bb(i,j,k,2) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                          ,kx,ky,kz,b2_coef,work)             !Flux coordinate

              if (enf_tor_flx_fn) then
                bb(i,j,k,3) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                            ,kx,ky,kz,b3_coef,work)             !Cov
                bb(i,j,k,3) = (bb(i,j,k,3)                                        &
     &                    -(gv%gparams%gmetric%grid(igrid)%gsub(i,j,k,3,1)*bb(i,j,k,1)   &
     &                     +gv%gparams%gmetric%grid(igrid)%gsub(i,j,k,3,2)*bb(i,j,k,2))) &
     &                    /gv%gparams%gmetric%grid(igrid)%gsub(i,j,k,3,3)       !Xform to cnv
              else
                bb(i,j,k,3) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs  &
     &                            ,kx,ky,kz,b3_coef,work)             !Cnv
              endif

!This code does not work for hard stellarator cases
!!$              b3(i,j,k) = gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,3,2)          &
!!$     &                   /gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,2,2)*b2(i,j,k)&
!!$     &                  +(gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,3,3)          &
!!$     &                   -gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,3,2)**2       &
!!$     &                   /gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,2,2))         &
!!$     &                  *db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs   &
!!$     &                         ,kx,ky,kz,bsub3_coef,work)

              prs(i,j,k) = db3val(r1,th1,v1,0,0,0,tx,ty,tz,nxs,nys,nzs &
     &                           ,kx,ky,kz,prs_coef,work)

              rho(i,j,k) = sign(1d0,prs(i,j,k))*abs(prs(i,j,k))**(1d0/gam)
            enddo
          enddo
        enddo

        if (max_rho == 0d0) then
          rho = 1d0
        else
          rho = rho/max_rho
        endif

!     Free work space

        DEALLOCATE(bsup,bsub3,prsg,stat=istat)

        DEALLOCATE(bsupuijcf,bsupvijcf,bsubvijcf,bsupsijsf,presijf,stat=istat)

        DEALLOCATE(xs,ys,zs,work,tx,ty,tz,b1_coef,b2_coef,b3_coef,prs_coef,stat=istat)

        IF (istat .ne. 0) STOP 'Deallocation error in vmec_equ'

        call vmec_cleanup

        call vmec_cleanup_metrics

!     End program

      end subroutine vmec_equ

!!$!     vmec_equ
!!$!     #################################################################
!!$      subroutine vmec_equ(igrid,nx,ny,nz,b1,b2,b3,prs,rho,gam)
!!$
!!$!     -----------------------------------------------------------------
!!$!     Give equilibrium fields at each logical mesh point in grid
!!$!     level (igrid).
!!$!     -----------------------------------------------------------------
!!$
!!$        use vmec_mod, pi_stel => pi
!!$        use app_iosetup
!!$        use grid
!!$        use equilibrium
!!$
!!$        implicit none
!!$
!!$!      Call variables
!!$
!!$        integer :: igrid,nx,ny,nz
!!$        real(8) :: gam
!!$        real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: b1,b2,b3,prs,rho
!!$
!!$!     Local variables
!!$
!!$        integer :: nxg,nyg,nzg,i,j,k,igl,jgl,kgl,ig,jg,kg,istat
!!$        real(8) :: r1,z1,th1,ph1,v1,dphi,dth,sgn,jac,ds,dum1,dum2,ppi,max_rho
!!$
!!$        real(8),allocatable, dimension(:)     :: pflx,ff_i,q_i
!!$
!!$        integer :: udcon=1111
!!$
!!$        logical :: enf_tor_flx_fn  !Whether we enforce flux functions in toroidal geom.
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
!!$        nxg = gv%gparams%nxgl(igrid)
!!$        nyg = gv%gparams%nygl(igrid)
!!$        nzg = gv%gparams%nzgl(igrid)
!!$
!!$!     Read equilibrium file and setup arrays
!!$!     [VMEC++ assumes solution is up-down symmetric wrt Z=0, 
!!$!     and hence only gives theta=(0,pi); thereby the limit nyg/2+1 in theta]
!!$
!!$        load_metrics = .true.  !Whether to use VMEC metrics to define B components
!!$
!!$        call vmec_init(nxg+1,nyg/2+1,nzg,equ_file)
!!$
!!$!     Setup equilibrium fields
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
!!$!     Transfer variables (from GLOBAL in VMEC to LOCAL in PIXIE3D)
!!$
!!$        do k=0,nz+1
!!$          do j=0,ny+1
!!$            do i=0,nx+1
!!$
!!$              call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)
!!$
!!$              !Find global limits
!!$              call fromLocalToGlobalLimits(i,j,k,igl,jgl,kgl,igrid,igrid,igrid)
!!$
!!$              !Singular point boundary
!!$              if (igl == 0) then
!!$                jgl = mod(jgl+nyg/2,nyg)
!!$                igl = 1
!!$              endif
!!$
!!$              !Periodic boundary in theta (other boundary enforced by symmetry)
!!$              if (jgl == 0) jgl = nyg
!!$
!!$              !Up-down symmetry in theta (temporary, until SIESTA interface is fixed)
!!$              sgn = 1d0
!!$              if (jgl > nyg/2+1) then
!!$                 jgl = nyg + 2 - jgl
!!$                 sgn = -1d0
!!$              endif
!!$
!!$              !Periodic boundary in phi
!!$              if (kgl == 0) kgl = nzg
!!$              if (kgl == nzg+1) kgl = 1
!!$
!!$              !Average to our integer radial mesh (half-mesh in VMEC)
!!$              call half_to_int(igl,jgl,kgl,bsupsijsf,b1 (i,j,k))
!!$              call half_to_int(igl,jgl,kgl,bsupuijcf,b2 (i,j,k))
!!$!!              call half_to_int(igl,jgl,kgl,bsupvijcf,b3 (i,j,k))
!!$              call half_to_int(igl,jgl,kgl,bsubvijcf,b3 (i,j,k))   !Use covariant component (flux function)
!!$              call half_to_int(igl,jgl,kgl,presijf  ,prs(i,j,k))
!!$
!!$              jac = sqrtg(igl+1,jgl,kgl)
!!$
!!$              !Transform b1, b2 to PIXIE's contravariant representation (b3 left covariant)
!!$              b1(i,j,k) = jac*b1(i,j,k)*0.5/gv%gparams%xx(ig)  !This correction comes because here
!!$                                                                !the variable is x1=sqrt(s), s-> VMEC.
!!$              b2(i,j,k) = jac*b2(i,j,k)   !Flux coordinate
!!$
!!$!!              b3(i,j,k) = gv%gparams%gmetric%grid(igrid)%jac(i,j,k)*b3(i,j,k)
!!$!!            Cov component               !Flux coordinate
!!$!!              b3(i,j,k) = gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,3,2)            &
!!$!!     &                   /gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,2,2)*b2(i,j,k)  &
!!$!!     &                  +(gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,3,3)            &
!!$!!     &                   -gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,3,2)**2         &
!!$!!     &                   /gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,2,2))*b3(i,j,k)
!!$
!!$              rho(i,j,k) = max(prs(i,j,k),0d0)**(1d0/gam)
!!$
!!$            enddo
!!$          enddo
!!$        enddo
!!$
!!$        max_rho = maxval(rho)
!!$        rho = rho/max_rho
!!$
!!$!       Clean flux functions (Tokamak case)
!!$        
!!$        enf_tor_flx_fn = (nfp_i == 1)
!!$        if (enf_tor_flx_fn.and.(np == 1)) then
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
!!$              dum1=sum(b2(i,1:nyg,k))/nyg
!!$              dum2=sum(b3(i,1:nyg,k))/nyg
!!$
!!$              b2(i,:,k) = dum1
!!$              b3(i,:,k) = dum2
!!$             
!!$            enddo
!!$          enddo
!!$        endif
!!$
!!$!  Get cnv b3 component
!!$
!!$        do k=0,nz+1
!!$          do j=0,ny+1
!!$            do i=0,nx+1
!!$              b3(i,j,k) = gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,3,2)            &
!!$     &                   /gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,2,2)*b2(i,j,k)  &
!!$     &                  +(gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,3,3)            &
!!$     &                   -gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,3,2)**2         &
!!$     &                   /gv%gparams%gmetric%grid(igrid)%gsup(i,j,k,2,2))*b3(i,j,k)
!!$            enddo
!!$          enddo
!!$        enddo
!!$
!!$!     Free work space
!!$
!!$        call vmec_cleanup
!!$
!!$        call vmec_cleanup_metrics
!!$
!!$!     End program
!!$
!!$      end subroutine vmec_equ

!     vmec_init_fields
!     #################################################################
      SUBROUTINE vmec_init_fields

!     -----------------------------------------------------------------
!     Gives contravariant magnetic fields components and pressure
!     at grid level igrid
!     -----------------------------------------------------------------

       USE stel_kinds
       USE stel_constants            
       USE island_params, ns=>ns_i, ntheta=>nu_i, nzeta=>nv_i,                &
    &          mpol=>mpol_i, ntor=>ntor_i, nuv=>nuv_i, mnmax=>mnmax_i,        &
    &          ohs=>ohs_i, nfp=>nfp_i
       USE vmec_mod, ONLY: bsupumncf_i,bsupvmncf_i,bsubvmncf_i,presf_i        &
                          ,bsupuijcf  ,bsupvijcf  ,bsubvijcf  ,bsupsijsf      &
                          ,presijf
       USE fourier, ONLY:  toijsp

       IMPLICIT NONE

!      Local variables

       INTEGER:: istat, jk

!      Begin program

!      Allocate variables

       ALLOCATE (bsupsijsf(ns,ntheta,nzeta)  &
                ,bsupuijcf(ns,ntheta,nzeta)  &
                ,bsupvijcf(ns,ntheta,nzeta)  &
                ,bsubvijcf(ns,ntheta,nzeta), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation error in vmec_init_fields'

       ALLOCATE(presijf(ns,ntheta,nzeta), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation error in vmec_init_fields'

!      Fill GLOBAL variables (at VMEC integer mesh -- PIXIE3D's half mesh)

!LC12/14/06       bsupsijsh = zero                       ! It is zero in VMEC (flux surfaces)
       bsupsijsf = zero                                ! It is zero in VMEC (flux surfaces)
       CALL toijsp (bsupumncf_i, bsupuijcf, 0, 0, 0, 0)
       CALL toijsp (bsupvmncf_i, bsupvijcf, 0, 0, 0, 0)
       CALL toijsp (bsubvmncf_i, bsubvijcf, 0, 0, 0, 0)

       DO jk = 1, ns
         presijf(jk,:,:) =  presf_i(jk)                ! Init pressure
       ENDDO

      END SUBROUTINE vmec_init_fields 
