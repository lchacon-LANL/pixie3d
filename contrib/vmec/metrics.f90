      MODULE metrics
      USE stel_kinds
      USE stel_constants
      USE island_params
      USE read_wout_mod, ns_vmec=>ns, mpol_vmec=>mpol, ntor_vmec=>ntor, &
     &   rmnc_vmec=>rmnc, zmns_vmec=>zmns, lmns_vmec=>lmns,             &
     &   xm_vmec=>xm, xn_vmec=>xn, iotaf_vmec=>iotaf, phipf_vmec=>phipf,&
     &   presf_vmec=>presf, nfp_vmec=>nfp, wb_vmec=>wb, wp_vmec=>wp
      IMPLICIT NONE
!     
!     WRITTEN 06-27-06 BY S. P. HIRSHMAN AS PART OF THE VMEC++ PROJECT (c)
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
     &             gss, gsu, gsv, guu, guv, gvv,                        &!symmetric elements of lower metric tensor (half mesh)
	 &             hss, hsu, hsv, huu, huv, hvv                          !symmetric elements of upper metric tensor (full mesh)

      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
     &             rmnc_spline, zmns_spline, lmns_spline

      REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE ::                     &
     &             rmnc_i, zmns_i, lmns_i, jbsupumncf_i, jbsupvmncf_i

      REAL(rprec) :: rmax, rmin, zmax, zmin

      REAL(rprec), DIMENSION(:,:), ALLOCATABLE:: gssf, guuf, gvvf,      &! Lower Metric elements (full)
     &   gsuf, gsvf, guvf   
!
!     LOCAL (PRIVATE) HELPER ROUTINES (NOT ACCESSIBLE FROM OUTSIDE THIS MODULE)
!
      PRIVATE spline_fourier_modes, loadrzl_vmec,       &
     &        convert_to_full_mesh

      CONTAINS

	  SUBROUTINE init_metric_elements(ns_in, mpol_in, ntor_in, wout_file)
	  IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER*(*)            :: wout_file
	  INTEGER, INTENT(in)      :: ns_in, mpol_in, ntor_in
	  INTEGER                  :: istat, ntype, imesh, js
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::                       &
     &                            jbsupumnc_v,  jbsupvmnc_v
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

!Set number of points == number of modes for now! (May want mid-points for flux conservation)
!SPH : # POINTS???

      nu_i = mpol_i+3  
      nv_i = 2*ntor_i+2
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

      WRITE (33, *) 'INITIAL VMEC PARAMETERS'
      WRITE (33, 100) wp_vmec/wb_vmec, phi(ns_vmec)
 100  FORMAT(' <BETA>: ',1pe12.4,' TFLUX: ',1pe12.4,/,                  &
     &       21('-'))
      
!
!     RS = add dump of RZ info to file for Line_tracing
!      
      IF (l_tracing) THEN
        OPEN(unit=77, file='tracing_VMEC.dat', status='unknown')
        WRITE(77, 165) 0, ns_vmec, mnmax, ntor_vmec, mpol_vmec
        DO mn = 1, mnmax
          WRITE(77, 166) INT(xm_vmec(mn)), INT(xn_vmec(mn)),            &
     &      rmnc_vmec(mn,1), zmns_vmec(mn,1), zero, zero
        ENDDO
        DO js = 2, ns_vmec
          DO mn = 1, mnmax
            WRITE(77, 167) rmnc_vmec(mn, js), zmns_vmec(mn, js),        &
     &        zero, zero
          ENDDO
        ENDDO
        CLOSE(unit=77)
      ENDIF
165   FORMAT(5(2x,i5))
166   FORMAT(2(2x,i5), 4(2x,1pe16.6))
167   FORMAT(4(2x,1pe16.6))
!
!     Allocate space for splined arrays
!
      ALLOCATE(rmnc_spline(mnmax,ns_i), zmns_spline(mnmax,ns_i),        &
     &         lmns_spline(mnmax,ns_i),                                 &
     &         phipf_i(ns_i), iotaf_i(ns_i), presf_i(ns_i),             &
     &         stat=istat)
      IF (istat .ne. 0) STOP 'Allocation error 1 in INIT_METRIC_ELEMENTS'
!
!     SPLINE R, Z. L FOURIER COMPONENTS in s FROM ORIGINAL VMEC MESH (s ~ phi, ns_vmec points) 
!     TO A "POLAR" MESH [s ~ sqrt(phi), ns_i POINTS] WITH BETTER AXIS RESOLUTION
!
      DO ntype = 1, 3
         IF (ntype .eq. 1) THEN
            istat = 0
            CALL Spline_Fourier_Modes(rmnc_vmec, rmnc_spline, istat)
            IF (istat .ne. 0) EXIT
         ELSE IF (ntype .eq. 2) THEN
            istat = 0
            CALL Spline_Fourier_Modes(zmns_vmec, zmns_spline, istat)   
            IF (istat .ne. 0) EXIT
         ELSE 
            CALL Convert_To_Full_Mesh(lmns_vmec)
            istat = 0
            CALL Spline_Fourier_Modes(lmns_vmec, lmns_spline, istat)   
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
!
!     Scale phipf_i and convert to sqrt(flux) mesh by multiplying by 2*s
!
      phipf_i = phipf_i / (2*pi)
      DO js = 1, ns_i
         phipf_i(js) = 2 * hs_i*(js-1) * phipf_i(js)
      END DO

!
!     CONSTRUCT R, Z, L REAL-SPACE ARRAYS ON "POLAR" MESH
!     AND COMPUTE METRIC ELEMENTS AND JACOBIAN
!
      CALL LoadRZL_VMEC(istat)
      IF (istat .ne. 0) STOP 'LoadRZL error in INIT_METRIC_ELEMENTS'

      CALL get_full_lower_metrics         ! RS: used to update current for resistive runs

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


      DO modes = 1, mnmax
         y_vmec => ymn_vmec(modes,:)
         y_spline => ymn_spline(modes,:)
         mp = xm_vmec(modes)

         IF (istat.eq.0 .and. mp.gt.0) THEN
            IF (MOD(mp,2) .eq. 1) THEN 
               expm = 1
            ELSE
               expm = 2
            END IF
            y_vmec = y_vmec*(fac1**expm)
            IF (mp .le. 2) y_vmec(1) = 2*y_vmec(2) - y_vmec(3)
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

         IF (istat.eq.0 .and. mp.gt.0) THEN
            y_spline = y_spline*(fac2**expm)
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
!     2. Set up s-nodes on final (splined) mesh [s_vmec(sfinal) ~ sfinal**2]
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
      fac1(1) = 0;  fac2(1) = 0
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
      USE dump_output
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

      rmax = MAXVAL(r1_i);      rmin = MINVAL(r1_i)
      zmax = MAXVAL(z1_i);      zmin = MINVAL(z1_i)
!
!     DEBUG: CHECK IF CORRECT
!
!      istat = (ns_i-1)/2 + 1    !rho = 1/2 => s = 1/4
!      CALL dump_special(r1_i,z1_i,ru_i,zu_i,rv_i,zv_i,istat)

!
!     COMPUTE HALF-MESH LOWER/UPPER METRIC ELEMENTS AND JACOBIAN
!
      CALL half_mesh_metrics (r1_i, ru_i, rv_i, z1_i, zu_i, zv_i)

!
!     CLEAN-UP EXTRA ARRAYS
!
!     IF(.NOT.l_tracing) DEALLOCATE (rmnc_i, zmns_i, stat=istat)    ! RS: Do not deallocate in L_TRACING
      DEALLOCATE (r1_i, z1_i, ru_i, zu_i, rv_i, zv_i, stat=istat)
!      IF (lwout_opened) CALL read_wout_deallocate

!
!     CONVERT UPPER METRIC ELEMENTS TO FULL MESH
      ALLOCATE (temp(ns_i, nuv_i))
      temp = hss
      CALL to_full_mesh_sp (temp, hss)
      temp = hsu
      CALL to_full_mesh_sp (temp, hsu)
      temp = hsv
      CALL to_full_mesh_sp (temp, hsv)
      temp = huu
      CALL to_full_mesh_sp (temp, huu)
      huu(3,:) = (9._dp)*huu(4,:)/4     !  huu ~ 1/rho**2
      huu(2,:) = 4*huu(3,:)
      temp = huv
      CALL to_full_mesh_sp (temp, huv)
      temp = hvv
      CALL to_full_mesh_sp (temp, hvv)

      DEALLOCATE(temp)

      END SUBROUTINE LoadRZL_VMEC


      SUBROUTINE repack (istat)
      USE fourier, ONLY: orthonorm
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(out)     :: istat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER                  :: modes, js, m, n, ntype, n1
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
     &         lmns_i(ns_i,0:mpol_i,-ntor_i:ntor_i),                    &
     &         stat=istat)

      IF (istat .ne. 0) STOP 'Allocation error in REPACK'
      rmnc_i = 0;  zmns_i = 0;  lmns_i = 0

!
!     LOAD n>=0 ONLY FOR M=0
!

      IF (MAXVAL(xm_vmec).gt.mpol_i .or.                                &
     &    MAXVAL(ABS(xn_vmec))/nfp_vmec.gt.ntor_i)                      & 
     &    PRINT *,'It is recommended to increase the number of modes',  &
     &            ' to at least VMEC values!'
      

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
               lmns_i(js,m,n1) = lmns_i(js,m,n1)                        &
      &                        -SIGN(1,n)*lmns_spline(modes,js)
            ELSE
               rmnc_i(js,m,-n) = rmnc_spline(modes,js)
               zmns_i(js,m,-n) = zmns_spline(modes,js)
               lmns_i(js,m,-n) = lmns_spline(modes,js)
            ENDIF
         END DO
      END DO
!
! RS (10/03/06)  Divide out "orthonorm" factor use in Island FFT routines
! 
      DO js = 1, ns_i
         rmnc_i(js,:,:) = rmnc_i(js,:,:)/orthonorm
         zmns_i(js,:,:) = zmns_i(js,:,:)/orthonorm
         lmns_i(js,:,:) = lmns_i(js,:,:)/orthonorm
      END DO

      DEALLOCATE (rmnc_spline, zmns_spline, lmns_spline, stat=istat)

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


      SUBROUTINE cleanup_metric_elements
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
	
      END SUBROUTINE cleanup_metric_elements

      
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
        USE Fourier, ONLY: orthonorm, tomnsp, toijsp
        IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: jk, m, n, lk, istat, ntheta, nzeta, ntheta1
        INTEGER :: mpol, ntor
        REAL(rprec)    :: arg, cosi, sini
        REAL(rprec), ALLOCATABLE, DIMENSION(:,:,:) :: testr, testi
!-----------------------------------------------
        
        ntheta = nu_i; nzeta = nv_i
        ntheta1 = ntheta - 1
        mpol = mpol_i;  ntor = ntor_i

        ALLOCATE(cosmu(0:mpol,1:ntheta), cosmum(0:mpol,1:ntheta),      &
     &    cosmui(0:mpol,1:ntheta), sinmu(0:mpol,1:ntheta),             &
     &    sinmum(0:mpol,1:ntheta), sinmui(0:mpol,1:ntheta),            &
     &    cosnv(0:ntor,1:nzeta), cosnvn(0:ntor,1:nzeta),               &
     &    sinnv(0:ntor,1:nzeta), sinnvn(0:ntor,1:nzeta), stat=istat)

        IF (istat .ne. 0) STOP 'Allocation error in fixarray'
        
!        dnorm_i = 2.d0/((ntheta-1)*nzeta)                         ! Norm for Fourier-Space -> UV-Space

        dnorm_i = one/((ntheta-1)*nzeta)                          ! ST:Norm for Fourier-Space -> UV-Space
        DO jk = 1, ntheta
          DO m = 0, mpol
            arg = m*pi*(jk-1)/ntheta1                             ! Poloidal collocation points*m
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
            IF (jk == 1 .or. jk == ntheta) THEN
               cosmui(m, jk) = 0.5_dp*cosmui(m, jk)               ! 1st and last poloidal points are divided by 2!
               sinmui(m, jk) = 0.5_dp*sinmui(m, jk)
            ENDIF
          ENDDO
        ENDDO

!
!FOR NOW, WE ARE IN 1 FIELD PERIOD
!
        DO lk = 1, nzeta
          DO n = 0, ntor
            arg = n*twopi*(lk-1)/nzeta                            ! Toroidal collocation points*ABS(n*nfp_i); 
            cosi = COS(arg)                                       ! sign of "n" taken care later through ISIGN
            sini = SIN(arg)
            
!            IF (n .ne. 0) THEN
!              cosi = cosi*nfactor; sini = sini*nfactor           ! RS (10/03/06): Added n=0/m=0 normalization
!            ENDIF
!            IF (n .eq. 0) THEN
!              cosi = cosi/nfactor; sini = sini/nfactor           ! RS (10/03/06): Added n=0/m=0 normalization
!            ENDIF

            cosnv( n, lk)  = cosi                                  ! All normalization factors goes in poloidal quantities
            cosnvn(n, lk)  = n*nfp_i*cosi                         ! NFP is number of field periods
            sinnv( n, lk)  = sini
            sinnvn(n, lk)  = n*nfp_i*sini
          ENDDO
        ENDDO

!COMPUTE orthonorm FACTOR FOR cos(mu+nv), sin(mu+nv) BASIS (NOT cos(mu)cos(nv) basis!)
!BE SURE TO FIX THIS WHEN DOING ASYMMETRIC CASE!

        ALLOCATE (orthonorm(0:mpol,-ntor:ntor), stat=istat)
        IF (istat .ne. 0) STOP 'Allocate ORTHONORM failed in FIXARRAY!'
        orthonorm = nfactor
        orthonorm(0,0) = 1
        IF (nu_i .eq. (mpol_i+1)) orthonorm(mpol,:) = 1


        END SUBROUTINE FIXARRAY


        SUBROUTINE DEALLOC_FIXARRAY
        USE stel_kinds
        IMPLICIT NONE
        INTEGER(kind=iprec):: istat

        DEALLOCATE(cosmu, cosmum, cosmui, sinmu, sinmum, sinmui,            &
     &    cosnv, cosnvn, sinnv, sinnvn, stat=istat) 
        IF (istat .ne. 0) STOP 'Allocation error in dealloc_fixarray'
         
        END SUBROUTINE DEALLOC_FIXARRAY
        
        
        SUBROUTINE get_full_lower_metrics
        USE stel_kinds
        USE island_params, ns=>ns_i, nuv=>nuv_i
        IMPLICIT NONE  
        REAL(rprec), DIMENSION(1:ns,1:nuv):: det  
        INTEGER:: istat   
!
!     This subroutine gets the lower metric elements on the full mesh
!     avoiding interpolations. It is only called if l_RESISTIVE = .TRUE.
!
        det = hss*huu*hvv + 2*hsu*huv*hsv  - hsv*huu*hsv -              &
     &    hsu*hsu*hvv - hss*hsv*huv
        
        ALLOCATE(gssf(1:ns,1:nuv), guuf(1:ns,1:nuv),                    &
     &    gvvf(1:ns,1:nuv), gsuf(1:ns,1:nuv),                           &
     &    gsvf(1:ns,1:nuv), guvf(1:ns,1:nuv), stat = istat)
        IF (istat .ne. 0) STOP 'Problem alloc. in LOWER_METRIC_FULL'
        
        gssf = (huu*hvv - huv*huv)/det
        guuf = (hss*hvv - hsv*hsv)/det
        gvvf = (hss*huu - hsu*hsu)/det
        gsuf = (huv*hsv - hvv*hsu)/det
        gsvf = (hsu*huv - huu*hsv)/det
        guvf = (hsv*hsu - hss*huv)/det
                    
        END SUBROUTINE get_full_lower_metrics
       
        
        SUBROUTINE dealloc_full_lower_metrics
        USE stel_kinds
        IMPLICIT NONE
        INTEGER(kind=iprec):: istat
       
        DEALLOCATE(gssf, guuf, gvvf, gsuf, gsvf, guvf, stat = istat)
        IF (istat .ne. 0) STOP 'Problem dealloc. in LOWER_METRIC_FULL'
        
        END SUBROUTINE dealloc_full_lower_metrics

      END MODULE metrics
