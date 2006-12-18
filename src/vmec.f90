      MODULE vmec_mod

      USE stel_kinds
      USE stel_constants
      USE island_params
      USE read_wout_mod, ns_vmec=>ns, mpol_vmec=>mpol, ntor_vmec=>ntor, &
     &   rmnc_vmec=>rmnc, zmns_vmec=>zmns,                              &
     &   xm_vmec=>xm, xn_vmec=>xn, iotaf_vmec=>iotaf, phipf_vmec=>phipf,&
     &   presf_vmec=>presf, nfp_vmec=>nfp, bsupumnc_vmec=>bsupumnc,     &
     &   bsupvmnc_vmec=>bsupvmnc, wb_vmec=>wb

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
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
     &             sqrtg,                                               &!sqrt(g): Jacobian on half grid
     &             gss, gsu, gsv, guu, guv, gvv,                        &!symmetric elements of lower metric tensor (full mesh)
     &             hss, hsu, hsv, huu, huv, hvv                          !symmetric elements of upper metric tensor (half mesh)

      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
     &             rmnc_spline, zmns_spline, bsupumnc_spline,           &
     &             bsupvmnc_spline

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                     &
     &             rmnc_i, zmns_i, bsupumncf_i, bsupvmncf_i

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: rr,zz

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE:: bsupuijcf,           &
     &  bsupvijcf, bsupsijsf, presijf

!     LOCAL (PRIVATE) HELPER ROUTINES (NOT ACCESSIBLE FROM OUTSIDE THIS MODULE)
!
      PRIVATE spline_fourier_modes, loadrzl_vmec,       &
     &        convert_to_full_mesh, vmec_load_rz

      CONTAINS

!       vmec_init_equ
!       ################################################################################
        SUBROUTINE vmec_init_equ(ns_in,nu_in,nv_in,wout_file)

          IMPLICIT NONE

!       -----------------------------------------------------------------------------
!          LOADS VALUES FOR MESHES TO BE USED IN VMECPP ISLAND SOLVER
!          (GENERALLY DIFFERENT FROM VMEC MESHES. USE SUBSCRIPT _i FOR ISLAND VARIABLES)
!          READS wout_file FROM VMEC TO GET FOURIER COMPONENTS ON VMEC MESH
!          SPLINES THE VMEC COMPONENTS TO RADIAL ISLAND MESH (PHI->SQRT(PHI))
!          CALLS Fourier MODULE fixarray TO LOAD TRIG ARRAYS FOR COMPUTING ISLAND METRICS
!          COMPUTES gij (sub/sup) AND JACOBIAN ON ISLAND MESHES IN REAL SPACE
!
!          Created by L. Chacon 11/30/06
!       -----------------------------------------------------------------------------

!       Call variables

          CHARACTER*(*)            :: wout_file
          INTEGER, INTENT(in)      :: ns_in, nu_in, nv_in

!       Local variables

          INTEGER                  :: istat, ntype, imesh, js
          REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
     &                            bsupumnc_v,  bsupvmnc_v
          REAL(rprec) :: t1, t2
          INTEGER     :: m, n, mn

!       Begin program

!       Initialize variables

          ns_i = ns_in
          nu_i = nu_in
          nv_i = nv_in

          nsh   = ns_i-1
          nuv_i = nu_i*nv_i

!       READ-IN DATA FROM VMEC-PRODUCED WOUT FILE (LIBSTELL ROUTINE)

          CALL read_wout_file(wout_file, istat)
          IF (istat .ne. 0) STOP 'Read-wout error in INIT_METRIC_ELEMENTS'

          ntor_i = ntor_vmec   !Take from VMEC equilibrium file
          mpol_i = mpol_vmec   !Take from VMEC equilibrium file

          mnmax_i = (mpol_i + 1)*(2*ntor_i + 1)             ! Added RS. Contains total number of modes.

          nfp_i = nfp_vmec
          wb_i  = (4*pi*pi)*wb_vmec

          IF (wb_i .eq. 0._dp) STOP 'wb_vmec = 0!'

!       Allocate space for splined arrays

          ALLOCATE(rr(ns_i,nu_i,nv_i),zz(ns_i,nu_i,nv_i))

          ALLOCATE(rmnc_spline(mnmax,ns_i), zmns_spline(mnmax,ns_i),        &
     &             bsupumnc_v(mnmax,ns_vmec),  bsupvmnc_v(mnmax,ns_vmec),   &
     &             bsupumnc_spline(mnmax,ns_i), bsupvmnc_spline(mnmax,ns_i),&
     &             phipf_i(ns_i), iotaf_i(ns_i),presf_i(ns_i),              &
     &             stat=istat)
          IF (istat .ne. 0) STOP 'Allocation error 1 in INIT_METRIC_ELEMENTS'

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
                istat = 1
                CALL Spline_Fourier_Modes(bsupumnc_v, bsupumnc_spline, istat)   

                CALL Convert_From_Nyq(bsupvmnc_v,bsupvmnc_vmec)       !These have mnmax_nyq members
                CALL Convert_To_Full_Mesh(bsupvmnc_v)
                istat = 1
                CALL Spline_Fourier_Modes(bsupvmnc_v, bsupvmnc_spline, istat)   
                DEALLOCATE(bsupumnc_v, bsupvmnc_v)
             END IF

             IF (istat .ne. 0) STOP 'Spline error in INIT_METRIC_ELEMENTS'
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

          CALL vmec_load_rz(istat)
          IF (istat .ne. 0) STOP 'LoadRZL error in VMEC_INIT_EQU'

        END SUBROUTINE vmec_init_equ

!       vmec_load_rz
!       ################################################################################
        SUBROUTINE vmec_load_rz(istat)

!       ------------------------------------------------------------------------------
!       Computes R-Z map in real space from VMEC equilibrium
!
!       Created by L. Chacon 11/30/06
!       ------------------------------------------------------------------------------

          USE Fourier, ONLY: toijsp
!          USE dump_output

!       Call variables

          INTEGER, INTENT(out)     :: istat

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
          CALL toijsp (rmnc_i,rr,0,0,iparity,0)
      
          iparity = 1
          CALL toijsp (zmns_i,zz,0,0,iparity,0)

!       DEBUG: CHECK IF CORRECT
!
!        istat = (ns_i-1)/2 + 1    !rho = 1/2 => s = 1/4
!        CALL dump_special(r1_i,z1_i,ru_i,zu_i,rv_i,zv_i,istat)

!       CLEAN-UP EXTRA ARRAYS

          DEALLOCATE (rmnc_i, zmns_i, stat=istat)
          IF (lwout_opened) CALL read_wout_deallocate

        END SUBROUTINE vmec_load_rz

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

          DEALLOCATE(bsupumncf_i,bsupvmncf_i,stat=istat)
          IF (istat .ne. 0) STOP 'Deallocation error 3 in vmec_cleanup'

          call DEALLOC_FIXARRAY

          DEALLOCATE(orthonorm,stat=istat)
          IF (istat .ne. 0) STOP 'Deallocation error 4 in vmec_cleanup'

        END SUBROUTINE vmec_cleanup

!       init_metric_elements
!       ################################################################################
        SUBROUTINE init_metric_elements(ns_in,mpol_in,ntor_in,wout_file)

          IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
          CHARACTER*(*)            :: wout_file
          INTEGER, INTENT(in)      :: ns_in, mpol_in, ntor_in
          INTEGER                  :: istat, ntype, imesh, js
          REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
     &                            bsupumnc_v,  bsupvmnc_v
          REAL(rprec) :: t1, t2
          INTEGER     :: m, n, mn
!-----------------------------------------------
!
!     LOADS VALUES FOR MESHES TO BE USED IN VMECPP ISLAND SOLVER
!     (GENERALLY DIFFERENT FROM VMEC MESHES. USE SUBSCRIPT _i FOR ISLAND VARIABLES)
!     READS wout_file FROM VMEC TO GET FOURIER COMPONENTS ON VMEC MESH
!     SPLINES THE VMEC COMPONENTS TO RADIAL ISLAND MESH (PHI->SQRT(PHI))
!     CALLS Fourier MODULE fixarray TO LOAD TRIG ARRAYS FOR COMPUTING ISLAND METRICS
!     COMPUTES gij (sub/sup) AND JACOBIAN ON ISLAND MESHES IN REAL SPACE
!
          ns_i = ns_in
          nsh  = ns_i-1
          ntor_i = ntor_in
          mpol_i = mpol_in

! Set number of points == number of modes for now! (May want mid-points for flux conservation)
!SPH : # POINTS???
          nu_i = mpol_i + 3
          nv_i = 2*ntor_i + 2
          IF (ntor_i .eq. 0) nv_i = 1
          nuv_i = nu_i*nv_i
          mnmax_i = (mpol_i + 1)*(2*ntor_i + 1)             ! Added RS. Contains total number of modes.

!
!     READ-IN DATA FROM VMEC-PRODUCED WOUT FILE (LIBSTELL ROUTINE)
!
          CALL read_wout_file(wout_file, istat)
          IF (istat .ne. 0) STOP 'Read-wout error in INIT_METRIC_ELEMENTS'

          nfp_i = nfp_vmec
          wb_i  = (4*pi*pi)*wb_vmec

          IF (wb_i .eq. 0._dp) STOP 'wb_vmec = 0!'

!
!     Allocate space for splined arrays
!
          ALLOCATE(rmnc_spline(mnmax,ns_i), zmns_spline(mnmax,ns_i),        &
     &             bsupumnc_v(mnmax,ns_vmec),  bsupvmnc_v(mnmax,ns_vmec),   &
     &             bsupumnc_spline(mnmax,ns_i), bsupvmnc_spline(mnmax,ns_i),&
     &             phipf_i(ns_i), iotaf_i(ns_i),presf_i(ns_i),              &
     &             stat=istat)
          IF (istat .ne. 0) STOP 'Allocation error 1 in INIT_METRIC_ELEMENTS'

!
!     SPLINE R, Z. L FOURIER COMPONENTS in s FROM ORIGINAL VMEC MESH (s ~ phi, ns_vmec points) 
!     TO A "POLAR" MESH [s ~ sqrt(phi), ns_i POINTS] WITH BETTER AXIS RESOLUTION
!
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
                istat = 1
                CALL Spline_Fourier_Modes(bsupumnc_v, bsupumnc_spline, istat)   

                CALL Convert_From_Nyq(bsupvmnc_v,bsupvmnc_vmec)       !These have mnmax_nyq members
                CALL Convert_To_Full_Mesh(bsupvmnc_v)
                istat = 1
                CALL Spline_Fourier_Modes(bsupvmnc_v, bsupvmnc_spline, istat)   
                DEALLOCATE(bsupumnc_v, bsupvmnc_v)
             END IF

             IF (istat .ne. 0) STOP 'Spline error in INIT_METRIC_ELEMENTS'
          END DO
   
!
!     Spline 1-D arrays: careful -> convert phipf VMEC and multiply
!     by ds-vmec/ds-island, since phipf_i = d(PHI)/ds-island
!
          CALL Spline_OneD_Array (iotaf_vmec, iotaf_i, istat)
          CALL Spline_OneD_Array (phipf_vmec, phipf_i, istat)
          presf_vmec = mu0 * presf_vmec
          CALL Spline_OneD_Array (presf_vmec, presf_i, istat)
          DO js = 1, ns_i
             phipf_i(js) = 2 * hs_i*(js-1) * phipf_i(js) / (2*pi)
          END DO
          
!
!     CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON "POLAR" MESH
!     AND COMPUTE METRIC ELEMENTS AND JACOBIAN
!
          CALL LoadRZL_VMEC(istat)
          IF (istat .ne. 0) STOP 'LoadRZL error in INIT_METRIC_ELEMENTS'

        END SUBROUTINE init_metric_elements

      
      SUBROUTINE Spline_Fourier_Modes(ymn_vmec, ymn_spline, istat)
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
      REAL(rprec), DIMENSION(:), POINTER :: y_vmec, y_spline
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

      ymn_spline = 0d0   !L. Chacon, 11/29/06, to avoid definition error below

      DO modes = 1, mnmax
         y_vmec => ymn_vmec(modes,:)
         y_spline => ymn_spline(modes,:)
         mp = xm_vmec(modes)

         IF (istat.eq.0 .and. MOD(mp,2).eq.1) THEN
            y_vmec = y_vmec*fac1
            y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
         END IF

!
!        4. Initialize spline for each mode amplitude (factor out sqrt(s) factor for odd-m)
!
         yp1 = -1.e30_dp;  ypn = -1.e30_dp
         CALL spline (snodes_vmec, y_vmec, ns_vmec, yp1, ypn, y2_vmec)

!
!        5. Interpolate onto snodes mesh
!
         DO js = 1, ns_i
            CALL splint (snodes_vmec, y_vmec, y2_vmec, ns_vmec,         &
     &                   snodes(js), y_spline(js))
         END DO

         IF (istat.eq.0 .and. MOD(mp,2).eq.1) THEN
            y_spline = y_spline*fac2
         END IF

!
!        PRINT OUT FOR CHECKING
!
         IF (xn_vmec(modes) .eq. 0) THEN
!             WRITE (33, *) mp
             DO js = 1, ns_vmec
!                WRITE (33, '(1p2e14.4)') snodes_vmec(js), y_vmec(js)
             END DO
!                WRITE (33, *)
             DO js = 1, ns_i
!                WRITE (33, '(1p2e14.4)') snodes(js), y_spline(js)
             END DO
!                WRITE (33, *)
         END IF

      END DO

      istat = 0

      END SUBROUTINE Spline_Fourier_Modes
      

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


      SUBROUTINE LoadRZL_VMEC(istat)
      USE Fourier, ONLY: toijsp
!      USE dump_output
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(out)     :: istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns_vmec)  :: fac1, fac2
      INTEGER                          :: mpol_save, ntor_save, iparity
      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE    ::                  &
     &   r1_i, z1_i, ru_i, zu_i, rv_i, zv_i
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE :: temp
!-----------------------------------------------
!
!     1. LOAD FOURIER FACTORS USING VMEC mpol,ntor AND ISLAND nu,nv
!        IDEA: USE VMEC mpol, ntor VALUES AND CALL fixarray
!     2. COMPUTE METRIC ELEMENTS, JACOBIAN, LAMBDA ON ISLAND MESH
      
      CALL fixarray

!
!     COMPUTE R, Z (Metric elements, jacobian), LAMBDA (For initializing B)
!
!     FIRST, REPACK SPLINED ARRAYS FROM VMEC-ORDERING TO ISLAND-ORDERING
!     AND SWAP THE N->-N MODES TO BE CONSISTENT WITH mu+nv ISLAND ARGUMENT
!
      CALL repack (istat)

!
!     COMPUTE AND STORE R, Z AND THEIR ANGULAR DERIVATIVES, AND LAMBDA (NEED FOR VECTOR POT)
!
      ALLOCATE (r1_i(ns_i, nu_i, nv_i), z1_i(ns_i, nu_i, nv_i),         &
     &          ru_i(ns_i, nu_i, nv_i), zu_i(ns_i, nu_i, nv_i),         &
     &          rv_i(ns_i, nu_i, nv_i), zv_i(ns_i, nu_i, nv_i),         &
     &          stat = istat)

      IF (istat .ne. 0) STOP 'Allocation failed in LoadRZL_VMEC'

      iparity = 0
      CALL toijsp (rmnc_i, r1_i, 0, 0, iparity, 0)
      CALL toijsp (rmnc_i, ru_i, 1, 0, iparity, 0)
      CALL toijsp (rmnc_i, rv_i, 0, 1, iparity, 0)
      
      iparity = 1
      CALL toijsp (zmns_i, z1_i, 0, 0, iparity, 0)
      CALL toijsp (zmns_i, zu_i, 1, 0, iparity, 0)
      CALL toijsp (zmns_i, zv_i, 0, 1, iparity, 0)

!
!     DEBUG: CHECK IF CORRECT
!
      istat = (ns_i-1)/2 + 1    !rho = 1/2 => s = 1/4
!      CALL dump_special(r1_i,z1_i,ru_i,zu_i,rv_i,zv_i,istat)

!
!     COMPUTE HALF-MESH LOWER/UPPER METRIC ELEMENTS AND JACOBIAN
!
      CALL half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)

!
!     CLEAN-UP EXTRA ARRAYS
!
      DEALLOCATE (rmnc_i, zmns_i, stat=istat)
      DEALLOCATE (r1_i, z1_i, ru_i, zu_i, rv_i, zv_i, stat=istat)
      IF (lwout_opened) CALL read_wout_deallocate

!
!     CONVERT UPPER METRIC ELEMENTS TO FULL MESH
!!$      ALLOCATE (temp(ns_i, nuv_i))
!!$      temp = hss
!!$      CALL to_full_mesh_sp (temp, hss)
!!$      temp = hsu
!!$      CALL to_full_mesh_sp (temp, hsu)
!!$      temp = hsv
!!$      CALL to_full_mesh_sp (temp, hsv)
!!$      temp = huu
!!$      CALL to_full_mesh_sp (temp, huu)
!!$      huu(3,:) = (9._dp)*huu(4,:)/4     !  huu ~ 1/rho**2
!!$      huu(2,:) = 4*huu(3,:)
!!$      temp = huv
!!$      CALL to_full_mesh_sp (temp, huv)
!!$      temp = hvv
!!$      CALL to_full_mesh_sp (temp, hvv)

      DEALLOCATE(temp)

      END SUBROUTINE LoadRZL_VMEC


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
     &         stat=istat)

      IF (istat .ne. 0) STOP 'Allocation error in REPACK'
      rmnc_i = 0;  zmns_i = 0;  bsupumncf_i = 0; bsupvmncf_i = 0

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
            ELSE
               rmnc_i(js,m,-n) = rmnc_spline(modes,js)
               zmns_i(js,m,-n) = zmns_spline(modes,js)
               bsupumncf_i(js,m,-n) = bsupumnc_spline(modes,js)
               bsupvmncf_i(js,m,-n) = bsupvmnc_spline(modes,js)
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
      bsupvmncf_i(:,0,0) = bsupvmncf_i(:,0,0)*nfactor

      DEALLOCATE (rmnc_spline, zmns_spline, bsupumnc_spline,                &
     &            bsupvmnc_spline, stat=istat)

      END SUBROUTINE repack


      SUBROUTINE half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(1:ns_i,nuv_i), INTENT(in) ::               &
     &   r1_i, ru_i, rv_i, z1_i, zu_i, zv_i
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      REAL(rprec), PARAMETER :: p5 = 1.d0/2.d0, zero = 0
      INTEGER  :: js, lk, js1, istat
      REAL(rprec)            :: mintest, maxtest, r1, s1, t1, eps
      REAL(rprec), DIMENSION(nuv_i) ::                                  &
     &      r12, ru12, rv12, zu12, zv12, rs12, zs12, det
!-----------------------------------------------
!
!     ALLOCATE METRIC ELEMENT ARRAYS
!
      ALLOCATE(sqrtg(ns_i,nuv_i),                                       &
     &         gss(ns_i,nuv_i), gsu(ns_i,nuv_i),                        &
     &         gsv(ns_i,nuv_i), guu(ns_i,nuv_i),                        &
     &         guv(ns_i,nuv_i), gvv(ns_i,nuv_i),                        &
     &         hss(ns_i,nuv_i), hsu(ns_i,nuv_i),                        & 
     &         hsv(ns_i,nuv_i), huu(ns_i,nuv_i),                        &
     &         huv(ns_i,nuv_i), hvv(ns_i,nuv_i),                        &
     &         stat=istat)

!     COMPUTE ALL ON THE HALF MESH
      IF (istat .ne. 0) STOP 'Allocation error in HALF_MESH_METRICS'
      gss(1,:) = 0;  gsu(1,:) = 0;  gsv(1,:) = 0
      guu(1,:) = 0;  guv(1,:) = 0;  gvv(1,:) = 0
      sqrtg(1,:) = 0

      DO js = 2, ns_i
         js1 = js-1
         r12 = p5*(r1_i(js,:) + r1_i(js1,:))
         rs12 = (r1_i(js,:) - r1_i(js1,:))*ohs_i
         ru12= p5*(ru_i(js,:) + ru_i(js1,:))
         rv12= p5*(rv_i(js,:) + rv_i(js1,:))
         zs12 = (z1_i(js,:) - z1_i(js1,:))*ohs_i
         zu12= p5*(zu_i(js,:) + zu_i(js1,:))
         zv12= p5*(zv_i(js,:) + zv_i(js1,:))
         guu(js,:) = ru12*ru12 + zu12*zu12
         guv(js,:) = ru12*rv12 + zu12*zv12
         gvv(js,:) = r12*r12 + rv12*rv12 + zv12*zv12
         gsu(js,:) = rs12*ru12 + zs12*zu12
         gsv(js,:) = rs12*rv12 + zs12*zv12
         gss(js,:) = rs12*rs12 + zs12*zs12
         sqrtg(js,:) = r12*(ru12*zs12 - rs12*zu12)
      END DO

      mintest = MINVAL(sqrtg(2:,:))
      maxtest = MAXVAL(sqrtg(2:,:))

      IF (mintest*maxtest .le. zero) STOP "Jacobian changed sign!"

      hss(1,:) = 0;  hsu(1,:) = 0;  hsv(1,:) = 0
      huu(1,:) = 0;  huv(1,:) = 0;  hvv(1,:) = 0
!
!     Compute upper metric elements (inverse of lower matrix: det = sqrtg**2?)
!
      DO js = 2, ns_i
         det(:) = one/(sqrtg(js,:)*sqrtg(js,:))
         hss(js,:) = det(:)*                                            &
     &              (guu(js,:)*gvv(js,:) - guv(js,:)*guv(js,:))
         hsu(js,:) = det(:)*                                            &
     &              (guv(js,:)*gsv(js,:) - gsu(js,:)*gvv(js,:))
         hsv(js,:) = det(:)*                                            &
     &              (gsu(js,:)*guv(js,:) - gsv(js,:)*guu(js,:))
         huu(js,:) = det(:)*                                            &
     &              (gss(js,:)*gvv(js,:) - gsv(js,:)*gsv(js,:))
         huv(js,:) = det(:)*                                            &
     &              (gsv(js,:)*gsu(js,:) - gss(js,:)*guv(js,:))
         hvv(js,:) = det(:)*                                            &
     &              (gss(js,:)*guu(js,:) - gsu(js,:)*gsu(js,:))

      END DO

!
!     CONFIRM ACCURACY OF INVERSE
!
!     FIRST ROW
     
      maxtest = 0
      eps = 0.01_dp*SQRT(EPSILON(eps))

      DO js = 2, ns_i
        r12 = gss(js,:)*hss(js,:) + gsu(js,:)*hsu(js,:)                 &
     &       +gsv(js,:)*hsv(js,:)
        ru12 = gss(js,:)*hsu(js,:) + gsu(js,:)*huu(js,:)                &
     &       +gsv(js,:)*huv(js,:) 
        rs12 = gss(js,:)*hsv(js,:) + gsu(js,:)*huv(js,:)                &
     &       +gsv(js,:)*hvv(js,:) 
        r1 = MAXVAL(ABS(r12)-1)
        s1 = MAXVAL(ABS(ru12))
        t1 = MAXVAL(ABS(rs12))
        maxtest = MAX(maxtest,r1,s1,t1)
      END DO

    
      IF (maxtest .gt. eps) STOP "Error1 in metric elements"

!     SECOND ROW
     
      maxtest = 0

      DO js = 2, ns_i
        r12 = gsu(js,:)*hss(js,:) + guu(js,:)*hsu(js,:)                 &
     &       +guv(js,:)*hsv(js,:)
        ru12 = gsu(js,:)*hsu(js,:) + guu(js,:)*huu(js,:)                &
     &       +guv(js,:)*huv(js,:) 
        rs12 = gsu(js,:)*hsv(js,:) + guu(js,:)*huv(js,:)                &
     &       +guv(js,:)*hvv(js,:) 
        r1 = MAXVAL(ABS(r12))
        s1 = MAXVAL(ABS(ru12-1))
        t1 = MAXVAL(ABS(rs12))
        maxtest = MAX(maxtest,r1,s1,t1)
      END DO

      IF (maxtest .gt. eps) STOP "Error2 in metric elements"

!     THIRD ROW
     
      maxtest = 0

      DO js = 2, ns_i
        r12 = gsv(js,:)*hss(js,:) + guv(js,:)*hsu(js,:)                 &
     &       +gvv(js,:)*hsv(js,:)
        ru12 = gsv(js,:)*hsu(js,:) + guv(js,:)*huu(js,:)                &
     &       +gvv(js,:)*huv(js,:) 
        rs12 = gsv(js,:)*hsv(js,:) + guv(js,:)*huv(js,:)                &
     &       +gvv(js,:)*hvv(js,:) 
        r1 = MAXVAL(ABS(r12))
        s1 = MAXVAL(ABS(ru12))
        t1 = MAXVAL(ABS(rs12-1))
        maxtest = MAX(maxtest,r1,s1,t1)
      END DO

      IF (maxtest .gt. eps) STOP "Error3 in metric elements"

      END SUBROUTINE half_mesh_metrics


      SUBROUTINE vmec_cleanup_metric_elements
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER         :: istat
!-----------------------------------------------
	
!     Note: sqrtg is deallocated in init_bcovar and stored in jacob variable

      DEALLOCATE(gss, gsu, gsv, guu, guv, gvv,                          &
     &           hss, hsu, hsv, huu, huv, hvv,                          &
     &           phipf_i, iotaf_i, presf_i,                             &
     &           stat=istat)
	
      END SUBROUTINE vmec_cleanup_metric_elements

      
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

!       full_to_half
!       ###########################################################################
        subroutine full_to_half(igl,jgl,kgl,vmec_arr,local_val)

!       ---------------------------------------------------------------------------
!       Averages quantities from VMEC's full mesh (which is PIXIE's half mesh)
!       to PIXIE's collocated mesh. At ghost cells, it extrapolates using a 
!       first-order formula.
!       ---------------------------------------------------------------------------

!       Call variables

          integer(4)  :: igl,jgl,kgl
          real(rprec) :: vmec_arr(ns_i,nu_i,nv_i),local_val

!       Local variables

          if (igl == 0) then
             local_val = 1.5*vmec_arr(1,jgl,kgl)-0.5*vmec_arr(2,jgl,kgl)
          elseif (igl == ns_i) then
             local_val = 1.5*vmec_arr(igl,jgl,kgl)-0.5*vmec_arr(igl-1,jgl,kgl)
          else
             local_val = 0.5*(vmec_arr(igl,jgl,kgl)+vmec_arr(igl+1,jgl,kgl))
          endif

        end subroutine full_to_half

      END MODULE vmec_mod

!     vmec_map
!     #################################################################
      subroutine vmec_map(igrid,nx,ny,nz,xcar)

!     -----------------------------------------------------------------
!     Give Cartesian coordinates of each logical mesh point at grid
!     level (igrid).
!     -----------------------------------------------------------------

        use vmec_mod
        use equilibrium
        use grid

        implicit none

!     Input variables

        integer(4) :: igrid,nx,ny,nz
        real(8)    :: xcar(0:nx+1,0:ny+1,0:nz+1,3)

!     Local variables

        integer(4) :: nxg,nyg,nzg,i,j,k,igl,jgl,kgl,ig,jg,kg
        real(8)    :: r1,z1,ph1,sgn

!     Begin program

!     Get GLOBAL limits (VMEC operates on global domain)

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

!     Read equilibrium file and setup arrays
!     [VMEC++ assumes solution is up-down symmetric wrt Z=0, 
!     and hence only gives theta=(0,pi); thereby the limit nyg/2+1 in theta]

        call vmec_init_equ(nxg+1,nyg/2+1,nzg,equ_file)

!!$        write (*,*) 'DIAG--vmec_map',nxg,nyg,nzg
!!$        do k=1,nzg
!!$           write (*,*) 'slice=',k
!!$           do j=1,nyg/2+1
!!$              write (*,*) 'rr',j,k,rr(:,j,k)
!!$              write (*,*) 'zz',j,k,zz(:,j,k)
!!$           enddo
!!$           write (*,*)
!!$        enddo
!!$        stop

!     Transfer map (from GLOBAL in VMEC to LOCAL)

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              !Find global limits
              call fromLocalToGlobalLimits(i,j,k,igl,jgl,kgl,igrid,igrid,igrid)

              !Ensures we don't step over physical periodic boundaries
!!$              if (    (jgl < 1 .or. jgl > nyg)                       &
!!$                  .or.(kgl < 1 .or. kgl > nzg) ) cycle

              !Periodic boundary in theta (other boundary enforced by symmetry)
              if (jgl == 0) jgl = nyg

              !Up-down symmetry in theta (temporary, until VMEC++ interface is fixed)
              sgn = 1d0
              if (jgl > nyg/2+1) then
                 jgl = nyg + 2 - jgl
                 sgn = -1d0
              endif

              !Periodic boundary in phi
              if (kgl == 0) kgl = nzg
              if (kgl == nzg+1) kgl = 1

              !Average to our radial mesh (half-mesh in VMEC)
              !(except for radial ghost cells, where we extrapolate the boundary value)
              call full_to_half(igl,jgl,kgl,rr,r1)
              call full_to_half(igl,jgl,kgl,zz,z1)

              ph1 = -grid_params%zz(kg)  !Minus sign to preserve a right-handed ref. sys.

              !Transform to Cartesian geometry
              xcar(i,j,k,1)=r1*cos(ph1)
              xcar(i,j,k,2)=r1*sin(ph1)
              xcar(i,j,k,3)=sgn*z1

            enddo
          enddo
        enddo

!!$        write (*,*) 'DIAG--vmec_map',nxg,nyg,nzg
!!$        do k=1,nz
!!$           write (*,*) 'slice=',k
!!$           do j=1,ny
!!$              write (*,*) 'X',j,k,xcar(:,j,k,1)
!!$              write (*,*) 'Y',j,k,xcar(:,j,k,2)
!!$              write (*,*) 'Z',j,k,xcar(:,j,k,3)
!!$              write (*,*)
!!$           enddo
!!$           write (*,*)
!!$        enddo
!!$        stop

!     Free work space (to allow multiple calls to vmec_map for different grid levels)

        call vmec_cleanup

!     End program

      end subroutine vmec_map

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
       USE vmec_mod, ONLY: bsupumncf_i,bsupvmncf_i,presf_i                    &
                          ,bsupuijcf  ,bsupvijcf  ,bsupsijsf,presijf
       USE fourier, ONLY:  toijsp

       IMPLICIT NONE

!      Local variables

       INTEGER:: istat, jk

!      Begin program

!      Allocate variables

       ALLOCATE (bsupsijsf(ns,ntheta,nzeta)  &
                ,bsupuijcf(ns,ntheta,nzeta)  &
                ,bsupvijcf(ns,ntheta,nzeta), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation error in init_bcovar'

       ALLOCATE(presijf(ns,ntheta,nzeta), stat=istat)
       IF (istat .ne. 0) STOP 'Allocation error in init_fields'

!      Fill GLOBAL variables (at VMEC integer mesh -- PIXIE3D's half mesh)

!LC12/14/06       bsupsijsh = zero                       ! It is zero in VMEC (flux surfaces)
       bsupsijsf = zero                                ! It is zero in VMEC (flux surfaces)
       CALL toijsp (bsupumncf_i, bsupuijcf, 0, 0, 0, 0)
       CALL toijsp (bsupvmncf_i, bsupvijcf, 0, 0, 0, 0)

       DO jk = 1, ns
         presijf(jk,:,:) =  presf_i(jk)                ! Init pressure
       ENDDO

      END SUBROUTINE vmec_init_fields 

!     vmec_equ
!     #################################################################
      subroutine vmec_equ(igrid,nx,ny,nz,b1,b2,b3,prs)

!     -----------------------------------------------------------------
!     Give Cartesian coordinates of each logical mesh point at grid
!     level (igrid).
!     -----------------------------------------------------------------

        use vmec_mod
        use equilibrium
        use grid

        implicit none

!      Call variables

        integer(4) :: igrid,nx,ny,nz
        real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: b1,b2,b3,prs

!     Local variables

        integer(4) :: nxg,nyg,nzg,i,j,k,igl,jgl,kgl,ig,jg,kg,istat
        real(8)    :: r1,z1,ph1,sgn

!     Begin program

!     Get GLOBAL limits (VMEC operates on global domain)

        nxg = grid_params%nxgl(igrid)
        nyg = grid_params%nygl(igrid)
        nzg = grid_params%nzgl(igrid)

!     Read equilibrium file and setup arrays
!     [VMEC++ assumes solution is up-down symmetric wrt Z=0, 
!     and hence only gives theta=(0,pi); thereby the limit nyg/2+1 in theta]

        call vmec_init_equ(nxg+1,nyg/2+1,nzg,equ_file)

!     Setup equilibrium fields

        call vmec_init_fields

!     Transfer variables (from GLOBAL in VMEC to LOCAL)

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              !Find global limits
              call fromLocalToGlobalLimits(i,j,k,igl,jgl,kgl,igrid,igrid,igrid)

              !Ensures we don't step over physical periodic boundaries
!!$              if (    (jgl < 1 .or. jgl > nyg)                       &
!!$                  .or.(kgl < 1 .or. kgl > nzg) ) cycle

              !Periodic boundary in theta (other boundary enforced by symmetry)
              if (jgl == 0) jgl = nyg

              !Up-down symmetry in theta (temporary, until VMEC++ interface is fixed)
              sgn = 1d0
              if (jgl > nyg/2+1) then
                 jgl = nyg + 2 - jgl
                 sgn = -1d0
              endif

              !Periodic boundary in phi
              if (kgl == 0) kgl = nzg
              if (kgl == nzg+1) kgl = 1

              !Average to our radial mesh (half-mesh in VMEC)
              call full_to_half(igl,jgl,kgl,bsupsijsf,b1 (i,j,k))
              call full_to_half(igl,jgl,kgl,bsupuijcf,b2 (i,j,k))
              call full_to_half(igl,jgl,kgl,bsupvijcf,b3 (i,j,k))
              call full_to_half(igl,jgl,kgl,presijf  ,prs(i,j,k))

              !Transform to PIXIE's contravariant representation
              b1(i,j,k) = gmetric%grid(igrid)%jac(i,j,k)*b1(i,j,k)  &
                          *0.5/grid_params%xx(ig)   !This correction comes because here
                                                    ! the variable is x1=sqrt(s), s-> VMEC.
              b2(i,j,k) = gmetric%grid(igrid)%jac(i,j,k)*b2(i,j,k)
              b3(i,j,k) = gmetric%grid(igrid)%jac(i,j,k)*b3(i,j,k)
            enddo
          enddo
        enddo

!     Free work space (to allow multiple calls to vmec_map for different grid levels)

        DEALLOCATE(bsupuijcf,bsupvijcf, bsupsijsf, presijf,stat=istat)
        IF (istat .ne. 0) STOP 'Deallocation error in vmec_equ'

        call vmec_cleanup

!     End program

      end subroutine vmec_equ
