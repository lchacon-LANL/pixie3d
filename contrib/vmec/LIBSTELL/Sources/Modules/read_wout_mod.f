      MODULE read_wout_mod
!
!     USE READ_WOUT_MOD to include variables dynamically allocated
!     in this module
!     Call DEALLOCATE_READ_WOUT to free this memory when it is no longer needed
! 
!     Reads in output from VMEC equilibrium code(s), contained in wout file
!
!     Contained subroutines:
!
!     read_wout_file      wrapper alias called to read/open wout file      
!     read_wout_text      called by read_wout_file to read text file wout
!     read_wout_nc        called by read_wout_file to read netcdf file wout
!
!     Post-processing routines
!
!     mse_pitch           user-callable function to compute mse pitch angle
!                         for the computed equilibrium
!
      USE stel_kinds
      USE mgrid_mod

      IMPLICIT NONE
!DEC$ IF DEFINED (NETCDF)
C-----------------------------------------------
C   L O C A L   P A R A M E T E R S
C-----------------------------------------------
! Variable names (vn_...) : put eventually into library, used by read_wout too...
      CHARACTER(LEN=*), PARAMETER ::
     1  vn_version = 'version_',
     2  vn_extension = 'input_extension', vn_mgrid = 'mgrid_file',
     3  vn_magen = 'wb', vn_therm = 'wp', vn_gam = 'gamma', 
     4  vn_maxr = 'rmax_surf', vn_minr = 'rmin_surf', 
     5  vn_maxz = 'zmax_surf', vn_fp = 'nfp',
     6  vn_radnod = 'ns', vn_polmod = 'mpol', vn_tormod = 'ntor', 
     7  vn_maxmod = 'mnmax', vn_maxit = 'niter', vn_actit = 'itfsq',
     8  vn_asym = 'lasym', vn_recon = 'lrecon', vn_free = 'lfreeb',
     9  vn_error = 'ier_flag', vn_aspect = 'aspect', 
     1  vn_maxmod_nyq = 'mnmax_nyq',
     1  vn_beta = 'betatotal', vn_pbeta = 'betapol',
     2  vn_tbeta = 'betator', vn_abeta = 'betaxis',
     3  vn_b0 = 'b0', vn_rbt0 = 'rbtor0', vn_rbt1 = 'rbtor', 
     4  vn_sgs = 'signgs', vn_lar = 'IonLarmor', vn_modB = 'volavgB',
     5  vn_ctor = 'ctor', vn_amin = 'Aminor_p', vn_Rmaj = 'Rmajor_p', 
     6  vn_vol = 'volume_p', vn_am = 'am', vn_ai = 'ai', vn_ac = 'ac',
     6  vn_mse = 'imse', vn_thom = 'itse',
     7  vn_pmod = 'xm', vn_tmod = 'xn', vn_pmod_nyq = 'xm_nyq', 
     7  vn_tmod_nyq = 'xn_nyq',
     8  vn_racc = 'raxis_cc', vn_zacs = 'zaxis_cs', 
     1  vn_racs = 'raxis_cs', vn_zacc = 'zaxis_cc', vn_iotaf = 'iotaf',
     2  vn_presf = 'presf', vn_phi = 'phi', vn_phipf = 'phipf', 
     3  vn_jcuru = 'jcuru', vn_jcurv = 'jcurv', vn_iotah = 'iotas',
     4  vn_mass = 'mass', vn_presh = 'pres', vn_betah = 'beta_vol', 
     5  vn_buco = 'buco', vn_bvco = 'bvco', vn_vp = 'vp',
     6  vn_specw = 'specw', vn_phip = 'phips', vn_jdotb = 'jdotb', 
     6  vn_overr = 'over_r',
     7  vn_bgrv = 'bdotgradv', vn_merc = 'DMerc', vn_mshear = 'DShear',
     8  vn_mwell = 'DWell', vn_mcurr = 'DCurr', vn_mgeo = 'DGeod', 
     9  vn_equif = 'equif', vn_fsq = 'fsqt', vn_wdot = 'wdot', 
     1  vn_ftolv = 'ftolv', vn_fsql= 'fsql', vn_fsqr = 'fsqr', 
     1  vn_fsqz = 'fsqz', 
     1  vn_extcur = 'extcur', vn_curlab = 'curlabel', vn_rmnc = 'rmnc',
     2  vn_zmns = 'zmns', vn_lmns = 'lmns', vn_gmnc = 'gmnc', 
     3  vn_bmnc = 'bmnc', vn_bsubumnc = 'bsubumnc', 
     4  vn_bsubvmnc = 'bsubvmnc', vn_bsubsmns = 'bsubsmns', 
     5  vn_bsupumnc = 'bsupumnc', vn_bsupvmnc = 'bsupvmnc', 
     6  vn_rmns = 'rmns', vn_zmnc = 'zmnc',
     7  vn_lmnc = 'lmnc', vn_gmns = 'gmns', vn_bmns = 'bmns', 
     8  vn_bsubumns = 'bsubumns', vn_bsubvmns = 'bsubvmns', 
     9  vn_bsubsmnc = 'bsubsmnc', vn_bsupumns = 'bsupumns', 
     1  vn_bsupvmns = 'bsupvmns', 
     2  vn_rbc = 'rbc', vn_zbs = 'zbs', vn_rbs = 'rbs', vn_zbc = 'zbc'

! Long names (ln_...)
      CHARACTER(LEN=*), PARAMETER ::
     1  ln_version = 'VMEC Version',
     2  ln_extension = 'Input file extension',
     3  ln_mgrid = 'MGRID file',
     4  ln_magen = 'Magnetic Energy', ln_therm = 'Thermal Energy',
     5  ln_gam = 'Gamma', ln_maxr = 'Maximum R', ln_minr = 'Minimum R',
     6  ln_maxz = 'Maximum Z', ln_fp = 'Field Periods',
     7  ln_radnod = 'Radial nodes', ln_polmod = 'Poloidal modes',
     8  ln_tormod = 'Toroidal modes', ln_maxmod = 'Fourier modes',
     8  ln_maxmod_nyq = 'Fourier modes (Nyquist)',
     9  ln_maxit = 'Max iterations', ln_actit = 'Actual iterations',
     1  ln_asym = 'Asymmetry', ln_recon = 'Reconstruction',
     1  ln_free = 'Free boundary',
     2  ln_error = 'Error flag', ln_aspect = 'Aspect ratio',
     3  ln_beta = 'Total beta', ln_pbeta = 'Poloidal beta',
     4  ln_tbeta = 'Toroidal beta', ln_abeta = 'Beta axis',
     5  ln_b0 = 'RB-t over R axis', ln_rbt0 = 'RB-t axis',
     6  ln_rbt1 = 'RB-t edge', ln_sgs = 'Sign jacobian',
     7  ln_lar = 'Ion Larmor radius', ln_modB = 'avg mod B',
     8  ln_ctor = 'Toroidal current', ln_amin = 'minor radius',
     9  ln_Rmaj = 'major radius', ln_vol = 'Plasma volume',
     1  ln_mse = 'Number of MSE points', 
     1  ln_thom = 'Number of Thompson scattering points',
     1  ln_am = 'polynomial coeficients from mass(s)',
     1  ln_ac = 'polynomial coeficients from <J>(s)',
     1  ln_ai = 'polynomial coeficients from iota(s)',
     4  ln_pmod = 'Poloidal mode numbers', 
     5  ln_tmod = 'Toroidal mode numbers', 
     4  ln_pmod_nyq = 'Poloidal mode numbers (Nyquist)', 
     5  ln_tmod_nyq = 'Toroidal mode numbers (Nyquist)', 
     5  ln_racc = 'raxis (cosnv)', ln_racs = 'raxis (sinnv)',
     6  ln_zacs = 'zaxis (sinnv)', ln_zacc = 'zaxis (cosnv)',
     7  ln_iotaf = 'iota on full mesh',
     8  ln_presf = 'pressure on full mesh', 
     8  ln_phi = 'Toroidal flux on full mesh',
     9  ln_phipf = 'd(phi)/ds: Toroidal flux deriv on full mesh', 
     9  ln_jcuru = 'j dot gradu full',
     1  ln_jcurv = 'j dot gradv full', ln_iotah = 'iota half',
     2  ln_mass = 'mass half', ln_presh = 'pressure half',
     3  ln_betah = 'beta half', ln_buco = 'bsubu half',
     4  ln_bvco = 'bsubv half', ln_vp = 'volume deriv half',
     5  ln_specw = 'Spectral width half',
     6  ln_phip = 'tor flux deriv over 2pi half',
     7  ln_jdotb = 'J dot B', ln_bgrv = 'B dot grad v',
     8  ln_merc = 'Mercier criterion', ln_mshear = 'Shear Mercier',
     9  ln_mwell = 'Well Mercier', ln_mcurr = 'Current Mercier',
     1  ln_mgeo = 'Geodesic Mercier', ln_equif='Average force balance',
     1  ln_fsq = 'Residual decay',
     2  ln_wdot = 'Wdot decay', ln_extcur = 'External coil currents', 
     2  ln_fsqr = 'Residual decay - radial',
     2  ln_fsqz = 'Residual decay - vertical',
     2  ln_fsql = 'Residual decay - hoop',
     2  ln_ftolv = 'Residual decay - requested',
     3  ln_curlab = 'External current names', 

     3  ln_rmnc = 'cosmn component of cylindrical R, full mesh',
     4  ln_zmns = 'sinmn component of cylindrical Z, full mesh',
     4  ln_lmns = 'sinmn component of lambda, half mesh',  
     5  ln_gmnc = 'cosmn component of jacobian, half mesh',  
     6  ln_bmnc = 'cosmn component of mod-B, half mesh',  
     6  ln_bsubumnc = 'cosmn covariant u-component of B, half mesh', 
     6  ln_bsubvmnc = 'cosmn covariant v-component of B, half mesh', 
     7  ln_bsubsmns = 'sinmn covariant s-component of B, full mesh', 

     7  ln_bsupumnc = 'BSUPUmnc half',
     8  ln_bsupvmnc = 'BSUPVmnc half', 

     3  ln_rmns = 'sinmn component of cylindrical R, full mesh',
     4  ln_zmnc = 'cosmn component of cylindrical Z, full mesh',
     4  ln_lmnc = 'cosmn component of lambda, half mesh',  
     5  ln_gmns = 'sinmn component of jacobian, half mesh',  
     6  ln_bmns = 'sinmn component of mod-B, half mesh',  
     6  ln_bsubumns = 'sinmn covariant u-component of B, half mesh', 
     6  ln_bsubvmns = 'sinmn covariant v-component of B, half mesh', 
     7  ln_bsubsmnc = 'cosmn covariant s-component of B, full mesh', 

     4  ln_bsupumns = 'BSUPUmns half', ln_bsupvmns = 'BSUPVmns half',
     6  ln_rbc = 'Initial boundary R cos(mu-nv) coefficients', 
     7  ln_zbs = 'Initial boundary Z sin(mu-nv) coefficients',  
     8  ln_rbs = 'Initial boundary R sin(mu-nv) coefficients', 
     9  ln_zbc = 'Initial boundary Z cos(mu-nv) coefficients'
!DEC$ ENDIF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nfp, ns, mpol, ntor, mnmax, mnmax_nyq, itfsq, niter, 
     1    iasym, ireconstruct, ierr_vmec, imse, itse, nstore_seq, 
     2    isnodes, ipnodes, imatch_phiedge, isigng, mnyq, nnyq, ntmax
      REAL(rprec) :: wb, wp, gamma, pfac, rmax_surf, rmin_surf,
     1    zmax_surf, aspect, betatot, betapol, betator, betaxis, b0,
     2    tswgt, msewgt, flmwgt, bcwgt, phidiam, version_,
     3    delphid, IonLarmor, VolAvgB,
     3    fsql, fsqr, fsqz, ftolv,
     4    Aminor, Rmajor, Volume, RBtor, RBtor0, Itor
      REAL(rprec), ALLOCATABLE :: rzl_local(:,:,:,:)
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    rmnc, zmns, lmns, rmns, zmnc, lmnc, bmnc, gmnc, bsubumnc,
     2    bsubvmnc, bsubsmns, bsupumnc, bsupvmnc, currvmnc, 
     3    currumnc, bbc, raxis, zaxis
      REAL(rprec), DIMENSION(:,:), ALLOCATABLE ::
     1    bmns, gmns, bsubumns, bsubvmns, bsubsmnc, 
     2    bsupumns, bsupvmns, currumns, currvmns 
      REAL(rprec), DIMENSION(:), ALLOCATABLE ::
     1   iotas, iotaf, presf, phipf, mass, pres, beta_vol, xm, xn, 
     2   xm_nyq, xn_nyq, phip, buco, bvco, phi, vp, overr, jcuru, jcurv, 
     3   specw, jdotb, bdotgradv, fsqt, wdot, am, ac, ai,
     3   Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, extcur,
     4   sknots, ystark, y2stark, pknots, ythom, y2thom,
     5   anglemse, rmid, qmid, shear, presmid, alfa, curmid, rstark,
     6   qmeas, datastark, rthom, datathom, dsiobt
      LOGICAL :: lasym, lthreed, lwout_opened=.false.
      CHARACTER :: mgrid_file*200, input_extension*100

      INTEGER, PARAMETER :: norm_term_flag=0, 
     1   bad_jacobian_flag=1, more_iter_flag=2, jac75_flag=4

!     OVERLOAD SUBROUTINE READ_WOUT_FILE TO ACCEPT BOTH UNIT NO. (OPENED EXTERNALLY)
!     OR FILENAME (HANDLE OPEN/CLOSE HERE)
      INTERFACE read_wout_file
          MODULE PROCEDURE readw_and_open, readw_only
      END INTERFACE

!DEC$ IF DEFINED (NETCDF)
      PRIVATE :: read_wout_text, read_wout_nc
!DEC$ ELSE
      PRIVATE :: read_wout_text
!DEC$ ENDIF
      PRIVATE :: norm_term_flag, bad_jacobian_flag, 
     1           more_iter_flag, jac75_flag

      CONTAINS

      SUBROUTINE readw_and_open(file_or_extension, ierr, iopen)
      USE safe_open_mod
!DEC$ IF DEFINED(IFORT)
      USE IFPORT
!DEC$ ENDIF      
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: ierr
      INTEGER, OPTIONAL :: iopen
      CHARACTER(LEN=*), INTENT(in) :: file_or_extension
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, PARAMETER :: iunit_init = 10
      INTEGER :: iunit
      LOGICAL :: isnc
      CHARACTER(len=LEN_TRIM(file_or_extension)+10) :: filename
C-----------------------------------------------
!
!     THIS SUBROUTINE READS THE WOUT FILE CREATED BY THE VMEC CODE
!     AND STORES THE DATA IN THE READ_WOUT MODULE
!
!     FIRST, CHECK IF THIS IS A FULLY-QUALIFIED PATH NAME
!     MAKE SURE wout IS NOT EMBEDDED IN THE NAME (PERVERSE USER...)
!
      filename = 'wout'
      CALL parse_extension(filename, file_or_extension, isnc)
      CALL flush(6)
!SPH  IF (.not.isnc) STOP 'ISNC ERR IN READ_WOUT_MOD'
      IF (isnc) THEN
!DEC$ IF DEFINED (NETCDF)
         CALL read_wout_nc(filename, ierr)
!DEC$ ELSE
         PRINT *, "NETCDF wout file can not be opened on this platform"
         ierr = -100
!DEC$ ENDIF
      ELSE
         iunit = iunit_init
         CALL safe_open (iunit, ierr, filename, 'old', 'formatted')
         IF (ierr .eq. 0) CALL read_wout_text(iunit, ierr)
         CLOSE(unit=iunit)
      END IF

      IF (PRESENT(iopen)) iopen = ierr
      lwout_opened = (ierr .eq. 0)

      END SUBROUTINE readw_and_open


      SUBROUTINE readw_only(iunit, ierr, iopen)
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(out):: ierr
      INTEGER, OPTIONAL :: iopen
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat
      CHARACTER(LEN=256) :: vmec_version
      LOGICAL :: exfile
C-----------------------------------------------
!
!     User opened the file externally and has a unit number, iunit
!
      ierr = 0

      INQUIRE(unit=iunit, exist=exfile, name=vmec_version,iostat=istat)
      IF (istat.ne.0 .or. .not.exfile) THEN
        PRINT *,' In READ_WOUT_FILE, Unit = ',iunit,
     1          ' File = ',TRIM(vmec_version),' DOES NOT EXIST'
        IF (PRESENT(iopen)) iopen = -1
        ierr = -1
        RETURN
      ELSE
        IF (PRESENT(iopen)) iopen = 0
      END IF

      CALL read_wout_text(iunit, ierr)
      lwout_opened = (ierr .eq. 0)

      END SUBROUTINE readw_only


      SUBROUTINE read_wout_text(iunit, ierr)
      USE stel_constants, ONLY: mu0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: iunit, ierr
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: eps_w = 1.e-4_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: istat(15), i, j, k, js, m, n, n1, mn, nparts_in
      CHARACTER(LEN=256) :: vmec_version
      LOGICAL            :: lcurr
C-----------------------------------------------
!
!     THIS SUBROUTINE READS THE TEXT FILE WOUT CREATED BY THE VMEC CODE
!     AND STORES THE INFORMATION IN THE read_WOUT MODULE
!
!     CALL read_WOUT_FILE - GENERIC INTERFACE - CAN BE CALLED WITH EITHER UNIT NO. OR FILENAME
!
!     RMNC, ZMNS: FULL-GRID
!     LMNS      : HALF-GRID
!
      istat = 0
      ierr = 0
      nextcur = 0

      READ (iunit, '(a)', iostat=istat(2), err=1000) vmec_version

      i = INDEX(vmec_version,'=')
      IF (i .ge. 0) THEN
         READ(vmec_version(i+1:len_trim(vmec_version)),*) version_
      ELSE
         version_ = -1.0
      END IF

      ierr_vmec = norm_term_flag

      IF (version_ .le. (5.10 + eps_w)) THEN
         READ (iunit, *, iostat=istat(2), err=1000) wb, wp, gamma,
     1      pfac, nfp, ns,
     1      mpol, ntor, mnmax, itfsq, niter, iasym, ireconstruct
      ELSE
         IF (version_ .lt. 6.54) THEN
            READ (iunit, *, iostat=istat(2), err=1000) wb, wp, gamma,
     1         pfac, rmax_surf, rmin_surf
         ELSE
            READ (iunit, *, iostat=istat(2), err=1000) wb, wp, gamma,
     1         pfac, rmax_surf, rmin_surf, zmax_surf
         END IF
         IF (version_ .le. (8.0+eps_w)) THEN
            READ (iunit, *, iostat=istat(2), err=1000) nfp, ns, mpol, 
     1      ntor, mnmax, itfsq, niter, iasym, ireconstruct, ierr_vmec
            mnmax_nyq = mnmax
         ELSE
            READ (iunit, *, iostat=istat(2), err=1000) nfp, ns, mpol, 
     1      ntor, mnmax, mnmax_nyq, itfsq, niter, iasym, ireconstruct, 
     2      ierr_vmec
         END IF
      END IF

      lfreeb = (iasym .gt. 0)

      IF (version_ .gt. (6.20+eps_w)) THEN
         READ (iunit, *, iostat=istat(1), err=1000)imse, itse, nbsets,
     1      nobd, nextcur, nstore_seq
      ELSE
         READ (iunit, *, iostat=istat(1), err=1000)imse, itse, nbsets,
     1        nobd, nextcur
         nstore_seq = 100
      END IF

      IF (ierr_vmec.ne.norm_term_flag .and. ierr_vmec.ne.more_iter_flag)
     1   GOTO 1000

      IF (nextcur .gt. nigroup) istat(15) = -1

      IF (ALLOCATED(xm)) CALL read_wout_deallocate

      ALLOCATE (xm(mnmax), xn(mnmax), xm_nyq(mnmax_nyq), 
     1  xn_nyq(mnmax_nyq),rmnc(mnmax,ns), zmns(mnmax,ns),
     2  lmns(mnmax,ns), bmnc(mnmax_nyq,ns), gmnc(mnmax_nyq,ns), 
     3  bsubumnc(mnmax_nyq,ns), bsubvmnc(mnmax_nyq,ns), 
     4  bsubsmns(mnmax_nyq,ns), bsupumnc(mnmax_nyq,ns), 
     5  bsupvmnc(mnmax_nyq,ns), currvmnc(mnmax_nyq,ns), 
     5  iotas(ns), mass(ns), pres(ns), beta_vol(ns), phip(ns), 
     6  buco(ns), bvco(ns), phi(ns), iotaf(ns), presf(ns), phipf(ns),
     5  vp(ns), overr(ns), jcuru(ns), jcurv(ns), specw(ns), Dmerc(ns),
     6  Dshear(ns), Dwell(ns), Dcurr(ns), Dgeod(ns), equif(ns),
     7  extcur(nextcur), curlabel(nextcur), raxis(0:ntor,2),
     8  zaxis(0:ntor,2), jdotb(ns), bdotgradv(ns), 
     8  am(0:10), ac(0:10),  ai(0:10),
     9  fsqt(nstore_seq), wdot(nstore_seq), stat = istat(6))

      IF (lasym) 
     1   ALLOCATE (rmns(mnmax,ns), zmnc(mnmax,ns), lmnc(mnmax,ns), 
     2             bmns(mnmax_nyq,ns), gmns(mnmax_nyq,ns), 
     3             bsubumns(mnmax_nyq,ns), 
     3             bsubvmns(mnmax_nyq,ns), bsubsmnc(mnmax_nyq,ns), 
     4             bsupumns(mnmax_nyq,ns), bsupvmns(mnmax_nyq,ns), 
     5             stat=istat(6))

      fsqt = 0; wdot = 0; raxis = 0; zaxis = 0

      IF (nbsets .gt. 0) READ (iunit, *, iostat=istat(4), err=1000)
     1   (nbfld(i),i=1,nbsets)
      READ (iunit, '(a)', iostat=istat(5), err=1000) mgrid_file

      DO js = 1, ns
         DO mn = 1, mnmax
            IF(js .eq. 1) THEN
               READ (iunit, *, iostat=istat(7), err=1000) m, n
               xm(mn) = REAL(m,rprec)
               xn(mn) = REAL(n,rprec)
            END IF
            IF (version_ .le. (6.20+eps_w)) THEN
              READ (iunit, 730, iostat=istat(8), err=1000)
     1        rmnc(mn,js), zmns(mn,js), lmns(mn,js),
     2        bmnc(mn,js), gmnc(mn,js), bsubumnc(mn,js), 
     3        bsubvmnc(mn,js), bsubsmns(mn,js), 
     4        bsupumnc(mn,js), bsupvmnc(mn,js),
     5        currvmnc(mn,js)
            ELSE IF (version_ .le. (8.0+eps_w)) THEN
              READ (iunit, *, iostat=istat(8), err=1000)
     1        rmnc(mn,js), zmns(mn,js), lmns(mn,js),
     2        bmnc(mn,js), gmnc(mn,js), bsubumnc(mn,js), 
     3        bsubvmnc(mn,js), bsubsmns(mn,js), 
     4        bsupumnc(mn,js), bsupvmnc(mn,js),
     5        currvmnc(mn,js)
            ELSE
              READ (iunit, *, iostat=istat(8), err=1000)
     1        rmnc(mn,js), zmns(mn,js), lmns(mn,js)
            END IF
            IF (lasym) THEN
               IF (version_ .le. (8.0+eps_w)) THEN
                  READ (iunit, *, iostat=istat(8), err=1000)
     1            rmns(mn,js), zmnc(mn,js), lmnc(mn,js), 
     2            bmns(mn,js), gmns(mn,js), bsubumns(mn,js), 
     3            bsubvmns(mn,js), bsubsmnc(mn,js), 
     4            bsupumns(mn,js), bsupvmns(mn,js)
               ELSE
                  READ (iunit, *, iostat=istat(8), err=1000)
     1            rmns(mn,js), zmnc(mn,js), lmnc(mn,js)
               END IF
            END IF
            IF (js.eq.1 .and. m.eq.0) THEN
               n1 = ABS(n/nfp)
               IF (n1 .le. ntor) THEN
                  raxis(n1,1) = rmnc(mn,1)
                  zaxis(n1,1) = zmns(mn,1)
                  IF (lasym) THEN
                     raxis(n1,2) = rmns(mn,1)
                     zaxis(n1,2) = zmnc(mn,1)
                  END IF
               END IF
            END IF
         END DO

         IF (version_ .le. (8.0+eps_w)) CYCLE
         DO mn = 1, mnmax_nyq
            IF(js .eq. 1) THEN
               READ (iunit, *, iostat=istat(7), err=1000) m, n
               xm_nyq(mn) = REAL(m,rprec)
               xn_nyq(mn) = REAL(n,rprec)
            END IF
            READ (iunit, *, iostat=istat(8), err=1000)
     1        bmnc(mn,js), gmnc(mn,js), bsubumnc(mn,js), 
     2        bsubvmnc(mn,js), bsubsmns(mn,js), 
     3        bsupumnc(mn,js), bsupvmnc(mn,js)
            IF (lasym) THEN
               READ (iunit, *, iostat=istat(8), err=1000)
     2         bmns(mn,js), gmns(mn,js), bsubumns(mn,js), 
     3         bsubvmns(mn,js), bsubsmnc(mn,js), 
     4         bsupumns(mn,js), bsupvmns(mn,js)
            END IF
         END DO

      END DO

      mnyq = INT(MAXVAL(xm_nyq));  nnyq = INT(MAXVAL(ABS(xn_nyq)))/nfp

!
!     Read FULL AND HALF-MESH QUANTITIES 
!
!     NOTE: In version_ <= 6.00, mass, press were written out in INTERNAL (VMEC) units
!     and are therefore multiplied here by 1/mu0 to transform to pascals. Same is true
!     for ALL the currents (jcuru, jcurv, jdotb). Also, in version_ = 6.10 and
!     above, PHI is the true (physical) toroidal flux (has the sign of jacobian correctly
!     built into it)
!
      iotas(1) = 0; mass(1) = 0; pres(1) = 0; phip(1) = 0; 
      buco(1) = 0; bvco(1) = 0; vp(1) = 0; overr(1) = 0;  specw(1) = 1
      beta_vol(1) = 0

      IF (version_ .le. (6.05+eps_w)) THEN
         READ (iunit, 730, iostat=istat(9), err=1000)
     1     (iotas(js), mass(js), pres(js),
     2      phip(js), buco(js), bvco(js), phi(js), vp(js), overr(js),
     3      jcuru(js), jcurv(js), specw(js),js=2,ns)
         READ (iunit, 730, iostat=istat(10), err=1000)
     1      aspect, betatot, betapol, betator, betaxis, b0
      ELSE IF (version_ .le. (6.20+eps_w)) THEN
         READ (iunit, 730, iostat=istat(9), err=1000)
     1     (iotas(js), mass(js), pres(js), beta_vol(js),
     2      phip(js), buco(js), bvco(js), phi(js), vp(js), overr(js),
     3      jcuru(js), jcurv(js), specw(js),js=2,ns)
         READ (iunit, 730, iostat=istat(10), err=1000)
     1      aspect, betatot, betapol, betator, betaxis, b0
      ELSE IF (version_ .le. (6.95+eps_w)) THEN
         READ (iunit, *, iostat=istat(9), err=1000)
     1     (iotas(js), mass(js), pres(js), beta_vol(js),
     2      phip(js), buco(js), bvco(js), phi(js), vp(js), overr(js),
     3      jcuru(js), jcurv(js), specw(js),js=2,ns)
         READ (iunit, *, iostat=istat(10), err=1000)
     1      aspect, betatot, betapol, betator, betaxis, b0
      ELSE
         READ (iunit, *, iostat=istat(9), err=1000)
     1   (iotaf(js), presf(js), phipf(js), phi(js), 
     2   jcuru(js), jcurv(js), js=1,ns)
         READ (iunit, *, iostat=istat(9), err=1000)
     1   (iotas(js), mass(js), pres(js),
     2   beta_vol(js), phip(js), buco(js), bvco(js), vp(js),
     3   overr(js), specw(js),js=2,ns)
         READ (iunit, *, iostat=istat(10), err=1000)
     1      aspect, betatot, betapol, betator, betaxis, b0
      END IF


      IF (version_ .gt. (6.10+eps_w)) THEN
         READ (iunit, *, iostat=istat(10), err=1000) isigng
         READ (iunit, *, iostat=istat(10), err=1000) input_extension
         READ (iunit, *, iostat=istat(10), err=1000) IonLarmor,
     1     VolAvgB, RBtor0, RBtor, Itor, Aminor, Rmajor, Volume
      END IF


!-----------------------------------------------
!     MERCIER CRITERION
!-----------------------------------------------
      IF (version_.gt.(5.10+eps_w) .and. version_.lt.(6.20-eps_w)) THEN
         READ (iunit, 730, iostat=istat(11), err=1000)
     1      (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),
     2       Dgeod(js), equif(js), js=2,ns-1)
      ELSE IF (version_ .ge. (6.20-eps_w)) THEN
         READ (iunit, *, iostat=istat(11), err=1000)
     1      (Dmerc(js), Dshear(js), Dwell(js), Dcurr(js),
     2       Dgeod(js), equif(js), js=2,ns-1)
      END IF

      IF (nextcur .gt. 0) THEN
         IF (version_ .le. (6.20+eps_w)) THEN
            READ (iunit, 730, iostat=istat(12), err=1000)
     1      (extcur(i),i=1,nextcur)
         ELSE
            READ (iunit, *, iostat=istat(12), err=1000)
     1      (extcur(i),i=1,nextcur)
         END IF
         READ (iunit, *, iostat=istat(13)) lcurr
         IF (lcurr) READ (iunit, *, iostat=istat(13), err=1000)
     1      (curlabel(i),i=1,nextcur)
      END IF

      IF (version_ .le. (6.20+eps_w)) THEN
         READ (iunit, 730, iostat=istat(14))
     1     (fsqt(i), wdot(i), i=1,nstore_seq)
      ELSE
         READ (iunit, *, iostat=istat(14))
     1     (fsqt(i), wdot(i), i=1,nstore_seq)
      END IF

      IF ((version_.ge.6.20-eps_w) .and. (version_ .lt. (6.50-eps_w))
     1   .and. (istat(14).eq.0)) THEN
         READ (iunit, 730, iostat=istat(14), err=1000)
     1     (jdotb(js), bdotgradv(js),js=1,ns)
      ELSE IF (version_ .ge. (6.50-eps_w)) THEN
         READ (iunit, *, iostat=istat(14), err=1000)
     1     (jdotb(js), bdotgradv(js),js=1,ns)
      ELSE
         istat(14) = 0
      END IF

!
!     CONVERT FROM INTERNAL UNITS TO PHYSICAL UNITS IF NEEDED
!
      IF (version_ .le. (6.05+eps_w)) THEN
         mass = mass/mu0
         pres = pres/mu0
         jcuru = jcuru/mu0
         jcurv = jcurv/mu0
         jdotb = jdotb/mu0
         phi   = -phi
      END IF

!-----------------------------------------------
!     DATA AND MSE FITS
!-----------------------------------------------
      IF (ireconstruct .gt. 0) THEN

        n1 = MAXVAL(nbfld(:nbsets))
        ALLOCATE (sknots(isnodes), ystark(isnodes), y2stark(isnodes),
     1     pknots(ipnodes), ythom(ipnodes), y2thom(ipnodes),
     2     anglemse(2*ns), rmid(2*ns), qmid(2*ns), shear(2*ns),
     3     presmid(2*ns), alfa(2*ns), curmid(2*ns), rstark(imse),
     4     datastark(imse), rthom(itse), datathom(itse),
     5     dsiext(nobd), plflux(nobd), dsiobt(nobd), bcoil(n1,nbsets),
     6     plbfld(n1,nbsets), bbc(n1,nbsets))
         IF (imse.ge.2 .or. itse.gt.0) THEN
            READ (iunit, *) tswgt, msewgt
            READ (iunit, *) isnodes, (sknots(i),ystark(i),y2stark(i),
     1         i=1,isnodes)
            READ (iunit, *) ipnodes, (pknots(i), ythom(i),
     1         y2thom(i),i=1,ipnodes)
            READ(iunit, *)(anglemse(i),rmid(i),qmid(i),shear(i),
     1      presmid(i),alfa(i),curmid(i),i=1,2*ns-1)
            READ(iunit, *)(rstark(i),datastark(i),qmeas(i),i=1,imse)
            READ(iunit, *)(rthom(i),datathom(i),i=1,itse)
         END IF

         IF (nobd .gt. 0) THEN
            READ (iunit, *) (dsiext(i),plflux(i),dsiobt(i),i=1,nobd)
            READ (iunit, *) flmwgt
         END IF

         nbfldn = SUM(nbfld(:nbsets))
         IF (nbfldn .gt. 0) THEN
            DO n = 1, nbsets
               READ (iunit, *) (bcoil(i,n),plbfld(i,n),bbc(i,n),
     1            i=1,nbfld(n))
            END DO
            READ (iunit, *) bcwgt
         END IF

         READ (iunit, *) phidiam, delphid
!
!     READ Limiter & Prout plotting specs
!
         READ (iunit, *) nsets, nparts_in, nlim

         ALLOCATE (nsetsn(nsets))
         READ (iunit, *) (nsetsn(i),i=1,nsets)

         n1 = MAXVAL(nsetsn(:nsets))
         ALLOCATE (pfcspec(nparts_in,n1,nsets), limitr(nlim))

         READ (iunit, *) (((pfcspec(i,j,k),i=1,nparts_in),
     1      j=1,nsetsn(k)),k=1,nsets)

         READ (iunit, *) (limitr(i), i=1,nlim)

         m  = MAXVAL(limitr(:nlim))
         ALLOCATE (rlim(m,nlim), zlim(m,nlim))

         READ (iunit, *) ((rlim(i,j),zlim(i,j),i=1,limitr(j)),
     1      j=1,nlim)
         READ (iunit, *) nrgrid, nzgrid
         READ (iunit, *) tokid
         READ (iunit, *) rx1, rx2, zy1, zy2, condif
         READ (iunit, *) imatch_phiedge

      END IF

 1000 CONTINUE

      READ (iunit, iostat=ierr) mgrid_mode
      IF (ierr .ne. 0) THEN
         ierr = 0; mgrid_mode = 'N'
      END IF

      IF (istat(2) .ne. 0) ierr_vmec = 1

      DO m = 1,15
        IF (istat(m) .gt. 0) THEN
           PRINT *,' Error No. ',m,' in READ_WOUT, iostat = ',istat(m)
           ierr = m
           EXIT
        END IF
      END DO


  720 FORMAT(8i10)
  730 FORMAT(5e20.13)
  740 FORMAT(a)
  790 FORMAT(i5,/,(1p,3e12.4))

      END SUBROUTINE read_wout_text


!DEC$ IF DEFINED (NETCDF)
      SUBROUTINE read_wout_nc(filename, ierr)
      USE ezcdf
      USE vparams, ONLY: mu0
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(out) :: ierr
      CHARACTER(LEN=*), INTENT(in) :: filename
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: nwout, ierror, js
      INTEGER, DIMENSION(3)   :: dimlens
      REAL(rprec) :: ohs
      REAL(rprec), DIMENSION(:), ALLOCATABLE :: raxis_cc, raxis_cs,
     1                                          zaxis_cs, zaxis_cc
C-----------------------------------------------
! Open cdf File
      CALL cdf_open(nwout,filename,'r', ierr)
      IF (ierr .ne. 0) THEN
         PRINT *,' Error opening wout .nc file'
         RETURN
      END IF

! Be sure all arrays are deallocated
      CALL read_wout_deallocate

! Read in scalar variables
      CALL cdf_read(nwout, vn_error, ierr_vmec)
      
      IF (ierr_vmec.ne.norm_term_flag .and. ierr_vmec.ne.more_iter_flag)
     1   GOTO 1000

      CALL cdf_read(nwout, vn_version, version_)
      CALL cdf_read(nwout, vn_extension, input_extension)
      CALL cdf_read(nwout, vn_mgrid, mgrid_file)
      CALL cdf_read(nwout, vn_magen, wb)
      CALL cdf_read(nwout, vn_therm, wp)
      CALL cdf_read(nwout, vn_gam, gamma)
      CALL cdf_read(nwout, vn_maxr, rmax_surf)
      CALL cdf_read(nwout, vn_minr, rmin_surf)
      CALL cdf_read(nwout, vn_maxz, zmax_surf)
      CALL cdf_read(nwout, vn_fp, nfp)
      CALL cdf_read(nwout, vn_radnod, ns)
      CALL cdf_read(nwout, vn_polmod, mpol)
      CALL cdf_read(nwout, vn_tormod, ntor)
      CALL cdf_read(nwout, vn_maxmod, mnmax)
      mnmax_nyq = -1
      CALL cdf_read(nwout, vn_maxmod_nyq, mnmax_nyq)
      CALL cdf_read(nwout, vn_maxit, niter)
      CALL cdf_read(nwout, vn_actit, itfsq)
      CALL cdf_read(nwout, vn_asym, lasym)
      IF (lasym) iasym = 1
      CALL cdf_read(nwout, vn_recon, lrecon)
      IF (lrecon) ireconstruct = 1
      CALL cdf_read(nwout, vn_free, lfreeb)
      CALL cdf_read(nwout, vn_aspect, aspect)
      CALL cdf_read(nwout, vn_beta, betatot)
      CALL cdf_read(nwout, vn_pbeta, betapol)
      CALL cdf_read(nwout, vn_tbeta, betator)
      CALL cdf_read(nwout, vn_abeta, betaxis)
      CALL cdf_read(nwout, vn_b0, b0)
      CALL cdf_read(nwout, vn_rbt0, rbtor0)
      CALL cdf_read(nwout, vn_rbt1, rbtor)
      CALL cdf_read(nwout, vn_sgs, isigng)
      CALL cdf_read(nwout, vn_lar, IonLarmor)
      CALL cdf_read(nwout, vn_modB, volAvgB)
      CALL cdf_read(nwout, vn_ctor, Itor)
      CALL cdf_read(nwout, vn_amin, Aminor)
      CALL cdf_read(nwout, vn_rmaj, Rmajor)
      CALL cdf_read(nwout, vn_vol, volume)
      CALL cdf_read(nwout, vn_ftolv, ftolv)
      CALL cdf_read(nwout, vn_fsqr, fsqr)
      CALL cdf_read(nwout, vn_fsqz, fsqz)
      CALL cdf_read(nwout, vn_fsql, fsql)
      imse = -1
      IF (lrecon) THEN
         CALL cdf_read(nwout, vn_mse, imse)
         CALL cdf_read(nwout, vn_thom, itse)
      END IF
      CALL cdf_read(nwout, vn_nextcur, nextcur)

      mgrid_mode = 'N'
      CALL cdf_inquire(nwout, vn_mgmode, dimlens, ier=ierror)
      IF (ierror.eq.0) CALL cdf_read(nwout, vn_mgmode, mgrid_mode)
      IF (lfreeb) THEN
         CALL cdf_read(nwout, vn_flp, nobser)
         CALL cdf_read(nwout, vn_nobd, nobd)
         CALL cdf_read(nwout, vn_nbset, nbsets)
      END IF

! Inquire existence, dimensions of arrays for allocation
! 1D Arrays
      IF (lfreeb .and. nbsets.gt.0) THEN
         call cdf_read(nwout, vn_nbfld, nbfld)
      END IF

      CALL cdf_inquire(nwout, vn_pmod, dimlens)
      ALLOCATE (xm(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_tmod, dimlens)
      ALLOCATE (xn(dimlens(1)), stat = ierror)
      IF (mnmax_nyq .gt. 0) THEN
         CALL cdf_inquire(nwout, vn_pmod_nyq, dimlens)
         ALLOCATE (xm_nyq(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_tmod_nyq, dimlens)
         ALLOCATE (xn_nyq(dimlens(1)), stat = ierror)
      END IF

      CALL cdf_inquire(nwout, vn_racc, dimlens)
      ALLOCATE (raxis_cc(0:dimlens(1)-1), stat = ierror)
      CALL cdf_inquire(nwout, vn_zacs, dimlens)
      ALLOCATE (zaxis_cs(0:dimlens(1)-1), stat = ierror)
      IF (lasym) THEN
         CALL cdf_inquire(nwout, vn_racs, dimlens) 
         ALLOCATE (raxis_cs(0:dimlens(1)-1), stat = ierror)
         CALL cdf_inquire(nwout, vn_zacc, dimlens)
         ALLOCATE (zaxis_cc(0:dimlens(1)-1), stat = ierror)
      END IF

!  Profile coefficients, dimensioned from 0
      CALL cdf_inquire(nwout, vn_am, dimlens)
      ALLOCATE (am(0:dimlens(1)-1), stat = ierror)
      CALL cdf_inquire(nwout, vn_ac, dimlens)
      ALLOCATE (ac(0:dimlens(1)-1), stat = ierror)
      CALL cdf_inquire(nwout, vn_ai, dimlens)
      ALLOCATE (ai(0:dimlens(1)-1), stat = ierror)

      CALL cdf_inquire(nwout, vn_iotaf, dimlens)
      ALLOCATE (iotaf(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_presf, dimlens) 
      ALLOCATE (presf(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_phi, dimlens)
      ALLOCATE (phi(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_phipf, dimlens)
      ALLOCATE (phipf(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_jcuru, dimlens)
      ALLOCATE (jcuru(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_jcurv, dimlens)
      ALLOCATE (jcurv(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_iotah, dimlens)
      ALLOCATE (iotas(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mass, dimlens) 
      ALLOCATE (mass(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_presh, dimlens)
      ALLOCATE (pres(dimlens(1)), stat = ierror)
            CALL cdf_inquire(nwout, vn_am, dimlens)
      ALLOCATE (am(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_ac, dimlens)
      ALLOCATE (ac(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_ai, dimlens)
      ALLOCATE (ai(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_betah, dimlens)
      ALLOCATE (beta_vol(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_buco, dimlens) 
      ALLOCATE (buco(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bvco, dimlens)
      ALLOCATE (bvco(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_vp, dimlens)
      ALLOCATE (vp(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_specw, dimlens) 
      ALLOCATE (specw(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_phip, dimlens)
      ALLOCATE (phip(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_overr, dimlens)
      ALLOCATE (overr(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_jdotb, dimlens)
      ALLOCATE (jdotb(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bgrv, dimlens)
      ALLOCATE (bdotgradv(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_merc, dimlens)
      ALLOCATE (Dmerc(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mshear, dimlens)
      ALLOCATE (Dshear(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mwell, dimlens)
      ALLOCATE (Dwell(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mcurr, dimlens)
      ALLOCATE (Dcurr(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_mgeo, dimlens)
      ALLOCATE (Dgeod(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_equif, dimlens)
      ALLOCATE (equif(dimlens(1)), stat = ierror)

      CALL cdf_inquire(nwout, vn_fsq, dimlens)
      ALLOCATE (fsqt(dimlens(1)), stat = ierror)
      CALL cdf_inquire(nwout, vn_wdot, dimlens)
      ALLOCATE (wdot(dimlens(1)), stat = ierror)

      IF (nextcur .gt. 0) THEN
         CALL cdf_inquire(nwout, vn_extcur, dimlens)
         ALLOCATE (extcur(dimlens(1)), stat = ierror)
         CALL cdf_inquire(nwout, vn_curlab, dimlens)
         ALLOCATE (curlabel(dimlens(1)), stat = ierror)
      ENDIF

! 2D Arrays
      CALL cdf_inquire(nwout, vn_rmnc, dimlens)
      ALLOCATE (rmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_zmns, dimlens)
      ALLOCATE (zmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_lmns, dimlens)
      ALLOCATE (lmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_gmnc, dimlens)
      ALLOCATE (gmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bmnc, dimlens)
      ALLOCATE (bmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubumnc, dimlens)
      ALLOCATE (bsubumnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubvmnc, dimlens)
      ALLOCATE (bsubvmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubsmns, dimlens)
      ALLOCATE (bsubsmns(dimlens(1),dimlens(2)), stat = ierror)

!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
      CALL cdf_inquire(nwout, vn_bsupumnc, dimlens)
      ALLOCATE (bsupumnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsupvmnc, dimlens)
      ALLOCATE (bsupvmnc(dimlens(1),dimlens(2)), stat = ierror)

      IF (.NOT. lasym) GO TO 800

      CALL cdf_inquire(nwout, vn_rmns, dimlens)
      ALLOCATE (rmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_zmnc, dimlens)
      ALLOCATE (zmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_lmnc, dimlens)
      ALLOCATE (lmnc(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_gmns, dimlens)
      ALLOCATE (gmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bmns, dimlens)
      ALLOCATE (bmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubumns, dimlens)
      ALLOCATE (bsubumns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubvmns, dimlens)
      ALLOCATE (bsubvmns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsubsmnc, dimlens)
      ALLOCATE (bsubsmnc(dimlens(1),dimlens(2)), stat = ierror)

!     ELIMINATE THESE EVENTUALLY: DO NOT NEED THEM
      CALL cdf_inquire(nwout, vn_bsupumns, dimlens)
      ALLOCATE (bsupumns(dimlens(1),dimlens(2)), stat = ierror)
      CALL cdf_inquire(nwout, vn_bsupvmns, dimlens)
      ALLOCATE (bsupvmns(dimlens(1),dimlens(2)), stat = ierror)

 800  CONTINUE

! Read Arrays
      CALL cdf_read(nwout, vn_pmod, xm)
      CALL cdf_read(nwout, vn_tmod, xn)
      IF (mnmax_nyq .eq. -1) THEN
         mnmax_nyq = mnmax
         ALLOCATE (xm_nyq(mnmax_nyq), xn_nyq(mnmax_nyq))
         xm_nyq = xm;  xn_nyq = xn
      ELSE
         CALL cdf_read(nwout, vn_pmod_nyq, xm_nyq)
         CALL cdf_read(nwout, vn_tmod_nyq, xn_nyq)
      END IF

      mnyq = INT(MAXVAL(xm_nyq));  nnyq = INT(MAXVAL(ABS(xn_nyq)))/nfp

      CALL cdf_read(nwout, vn_racc, raxis_cc)
      CALL cdf_read(nwout, vn_zacs, zaxis_cs) 

      IF (SIZE(raxis_cc) .ne. ntor+1) 
     1   STOP 'WRONG SIZE(raxis_cc) in READ_WOUT_NC'
      ALLOCATE (raxis(0:ntor,2), zaxis(0:ntor,2))
      raxis(:,1) = raxis_cc(0:ntor);   zaxis(:,1) = zaxis_cs(0:ntor)
      raxis(:,2) = 0;                  zaxis(:,2) = 0
      DEALLOCATE (raxis_cc, zaxis_cs, stat=ierror)

      CALL cdf_read(nwout, vn_rmnc, rmnc)
      CALL cdf_read(nwout, vn_zmns, zmns)
      CALL cdf_read(nwout, vn_lmns, lmns)
      CALL cdf_read(nwout, vn_gmnc, gmnc)              !Half mesh
      CALL cdf_read(nwout, vn_bmnc, bmnc)              !Half mesh
      CALL cdf_read(nwout, vn_bsubumnc, bsubumnc)      !Half mesh
      CALL cdf_read(nwout, vn_bsubvmnc, bsubvmnc)      !Half mesh
      CALL cdf_read(nwout, vn_bsubsmns, bsubsmns)      !Full mesh
!     ELIMINATE  THESE EVENTUALLY: DON'T NEED THEM (can express in terms of lambdas)
      CALL cdf_read(nwout, vn_bsupumnc, bsupumnc)
      CALL cdf_read(nwout, vn_bsupvmnc, bsupvmnc)
      IF (lasym) THEN
         CALL cdf_read(nwout, vn_racs, raxis_cs)
         CALL cdf_read(nwout, vn_zacc, zaxis_cc) 
         raxis(:,2) = raxis_cs;   zaxis(:,2) = zaxis_cc
         DEALLOCATE (raxis_cs, zaxis_cc, stat=ierror)
         CALL cdf_read(nwout, vn_rmns, rmns)
         CALL cdf_read(nwout, vn_zmnc, zmnc)
         CALL cdf_read(nwout, vn_lmnc, lmnc)
         CALL cdf_read(nwout, vn_gmns, gmns)
         CALL cdf_read(nwout, vn_bmns, bmns) 
         CALL cdf_read(nwout, vn_bsubumns, bsubumns)
         CALL cdf_read(nwout, vn_bsubvmns, bsubvmns)
         CALL cdf_read(nwout, vn_bsubsmnc, bsubsmnc)
!     ELIMINATE THESE EVENTUALLY: DON'T NEED THEM
         CALL cdf_read(nwout, vn_bsupumns, bsupumns)
         CALL cdf_read(nwout, vn_bsupvmns, bsupvmns)
      END IF

      CALL cdf_read(nwout, vn_am, am)
      CALL cdf_read(nwout, vn_ac, ac)
      CALL cdf_read(nwout, vn_ai, ai)

      CALL cdf_read(nwout, vn_iotaf, iotaf) 
      CALL cdf_read(nwout, vn_presf, presf) 
      CALL cdf_read(nwout, vn_phi, phi) 
      CALL cdf_read(nwout, vn_phipf, phipf)
      CALL cdf_read(nwout, vn_jcuru, jcuru)
      CALL cdf_read(nwout, vn_jcurv, jcurv)
 
!     HALF-MESH quantities
!     NOTE: jdotb is in units_of_A (1/mu0 incorporated in jxbforce...)
!     prior to version 6.00, this was output in internal VMEC units...
      CALL cdf_read(nwout, vn_iotah, iotas)
      CALL cdf_read(nwout, vn_mass, mass) 
      CALL cdf_read(nwout, vn_presh, pres)
      CALL cdf_read(nwout, vn_betah, beta_vol)
      CALL cdf_read(nwout, vn_buco, buco)
      CALL cdf_read(nwout, vn_bvco, bvco) 
      CALL cdf_read(nwout, vn_vp, vp)
      CALL cdf_read(nwout, vn_specw, specw)
      CALL cdf_read(nwout, vn_phip, phip)
      CALL cdf_read(nwout, vn_jdotb, jdotb)
      CALL cdf_read(nwout, vn_bgrv, bdotgradv)

!     MERCIER_CRITERION
      CALL cdf_read(nwout, vn_merc, Dmerc)
      CALL cdf_read(nwout, vn_mshear, Dshear)
      CALL cdf_read(nwout, vn_mwell, Dwell)
      CALL cdf_read(nwout, vn_mcurr, Dcurr)
      CALL cdf_read(nwout, vn_mgeo, Dgeod)
      CALL cdf_read(nwout, vn_equif, equif)

      CALL cdf_read(nwout, vn_fsq, fsqt)
      CALL cdf_read(nwout, vn_wdot, wdot)
     
      IF (nextcur .gt. 0) THEN
         CALL cdf_read(nwout, vn_extcur, extcur)
         CALL cdf_read(nwout, vn_curlab, curlabel)
      ENDIF

 1000 CONTINUE

      CALL cdf_close(nwout, ierr)

!
!     COMPUTE CONTRAVARIANT CURRENT COMPONENTS (* JACOBIAN) currumn, currvmn IN AMPS 
!     ON THE FULL RADIAL MESH, WHERE JACOBIAN = SQRT(G)
!
!     CURRU = SQRT(G) * J dot grad(u)
!     CURRV = SQRT(G) * J dot grad(v)
!
      ohs = (ns-1)
     
      IF (.not.ALLOCATED(bsubumnc)) RETURN

      ALLOCATE (currumnc(mnmax_nyq,ns), currvmnc(mnmax_nyq,ns), 
     1          stat=ierror)

      DO js = 2, ns-1
         currumnc(:,js) = -xn_nyq(:)*bsubsmns(:,js)
     1          - ohs*(bsubvmnc(:,js+1) - bsubvmnc(:,js))
         currvmnc(:,js) = -xm_nyq(:)*bsubsmns(:,js)
     1          + ohs*(bsubumnc(:,js+1) - bsubumnc(:,js))
      END DO         

      WHERE (MOD(NINT(xm_nyq),2) .eq. 0)
         currvmnc(:,1) = 2*currvmnc(:,2) - currvmnc(:,3)
      ELSE WHERE 
          currvmnc(:,1) = 0
      END WHERE

      currumnc(:,1) = 2*currumnc(:,2) - currumnc(:,3)
      currumnc(:,ns) = 2*currumnc(:,ns-1) - currumnc(:,ns-2)
      currvmnc(:,ns) = 2*currvmnc(:,ns-1) - currvmnc(:,ns-2)
      currumnc = currumnc/mu0;   currvmnc = currvmnc/mu0

      IF (.not.lasym) RETURN

      ALLOCATE (currumns(mnmax_nyq,ns), currvmns(mnmax_nyq,ns), 
     1          stat=ierror)

      DO js = 2, ns-1
         currumns(:,js) = -xn_nyq(:)*bsubsmnc(:,js)
     1         -  ohs*(bsubvmns(:,js+1) - bsubvmns(:,js))
         currvmns(:,js) = -xm_nyq(:)*bsubsmnc(:,js)
     1          + ohs*(bsubumns(:,js+1) - bsubumns(:,js))
      END DO         

      currumns(:,1) = 2*currumns(:,2) - currumns(:,3)
      currumns(:,ns) = 2*currumns(:,ns-1) - currumns(:,ns-2)
      currvmns(:,1) = 2*currvmns(:,2) - currvmns(:,3)
      currvmns(:,ns) = 2*currvmns(:,ns-1) - currvmns(:,ns-2)
      currumns = currumns/mu0;   currvmns = currvmns/mu0
      IF (ierr. ne. 0) PRINT*,"in read_wout_nc ierr=",ierr
      IF (ierror. ne. 0) PRINT*,"in read_wout_nc ierror=",ierror

      END SUBROUTINE read_wout_nc
!DEC$ ENDIF

      SUBROUTINE read_wout_deallocate
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: istat(10)
!-----------------------------------------------
      istat = 0
      lwout_opened = .false.

      IF (ALLOCATED(extcur)) DEALLOCATE (extcur, curlabel,
     1         stat = istat(1))
      IF (ALLOCATED(overr)) DEALLOCATE (overr, stat = istat(2))

      IF (ALLOCATED(xm)) DEALLOCATE (xm, xn, xm_nyq, xn_nyq, 
     1  rmnc, zmns, lmns, bmnc, gmnc, bsubumnc, iotaf, presf, phipf,
     2  bsubvmnc, bsubsmns, bsupumnc, bsupvmnc, currvmnc, iotas, mass,
     3  pres, beta_vol, phip, buco, bvco, phi, vp, jcuru, am, ac, ai,
     4  jcurv, specw, Dmerc, Dshear, Dwell, Dcurr, Dgeod, equif, jdotb,
     5  bdotgradv, raxis, zaxis, fsqt, wdot, stat = istat(3))

      IF (ireconstruct.gt.0 .and. ALLOCATED(sknots)) DEALLOCATE (
     1    ystark, y2stark, pknots, anglemse, rmid, qmid, shear,
     2    presmid, alfa, curmid, rstark, datastark, rthom, datathom,
     3    ythom, y2thom, plflux, dsiobt, bcoil, plbfld, bbc, sknots,
     4    pfcspec, limitr, rlim, zlim, nsetsn, stat = istat(4))

      IF (ALLOCATED(rmns)) DEALLOCATE (rmns, zmnc, lmnc, 
     1    bmns, gmns, bsubumns, bsubvmns, bsubsmnc, 
     2    bsupumns, bsupvmns, currvmns, stat=istat(5))

      IF (ALLOCATED(currumnc)) DEALLOCATE (currumnc)
      IF (ALLOCATED(currumns)) DEALLOCATE (currumns)
      IF (ALLOCATED(rzl_local)) DEALLOCATE (rzl_local)

      IF (ANY(istat .ne. 0)) 
     1      STOP 'Deallocation error in read_wout_deallocate'

      END SUBROUTINE read_wout_deallocate


      SUBROUTINE tosuvspace (s_in, u_in, v_in, gsqrt,
     1     bsupu, bsupv)
      USE stel_constants, ONLY: zero, one
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(rprec), INTENT(in) :: s_in, u_in, v_in
      REAL(rprec), INTENT(out), OPTIONAL :: gsqrt, bsupu, bsupv
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: c1p5 = 1.5_dp
      INTEGER :: m, n, n1, mn, ipresent, jslo, jshi
      REAL(rprec) :: hs1, wlo, whi, wlo_odd, whi_odd
      REAL(rprec), DIMENSION(mnmax_nyq) :: gmnc1, gmns1, bsupumnc1,
     1   bsupumns1, bsupvmnc1, bsupvmns1, wmins, wplus
      REAL(rprec) :: cosu, sinu, cosv, sinv, tcosmn, tsinmn, sgn
      REAL(rprec) :: cosmu(0:mnyq), sinmu(0:mnyq),
     1               cosnv(0:nnyq), sinnv(0:nnyq)
      LOGICAL :: lgsqrt, lbsupu, lbsupv
C-----------------------------------------------
!
!     COMPUTE VARIOUS HALF-RADIAL GRID QUANTITIES AT THE INPUT POINT
!     (S, U, V) , WHERE 
!        S = normalized toroidal flux (0 - 1),
!        U = poloidal angle 
!        V = N*phi = toroidal angle * no. field periods
!
C-----------------------------------------------
      IF (s_in.lt.zero .or. s_in.gt.one) THEN
         WRITE(6, *)
     1   ' In tosuvspace, s(flux) must be between 0 and 1'
         RETURN
      END IF

      IF (.not.lwout_opened) THEN
         WRITE(6, *)
     1   ' tosuvspace can only be called AFTER opening wout file!'
         RETURN
      END IF
!
!     FIND INTERPOLATED s VALUE AND COMPUTE INTERPOLATION WEIGHTS wlo, whi
!     RECALL THAT THESE QUANTITIES ARE ON THE HALF-RADIAL GRID...
!     s-half(j) = (j-1.5)*hs, for j = 2,...ns
!
      hs1 = one/(ns-1)
      jslo = INT(c1p5 + s_in/hs1)
      jshi = jslo+1
      wlo = (hs1*(jshi-c1p5) - s_in)/hs1
      whi = 1 - wlo
      IF (jslo .eq. ns) THEN
!        USE Xhalf(ns+1) = 2*Xhalf(ns) - Xhalf(ns-1) FOR "GHOST" POINT VALUE 1/2hs OUTSIDE EDGE
!        THEN, X = wlo*Xhalf(ns) + whi*Xhalf(ns+1) == Xhalf(ns) + whi*(Xhalf(ns) - Xhalf(ns-1)) 
         jshi = jslo-1
         wlo = 1+whi; whi = -whi
      ELSE IF (jslo .eq. 1) THEN
         jslo = 2
      END IF

!
!     FOR ODD-m MODES X ~ SQRT(s), SO INTERPOLATE Xmn/SQRT(s)
! 
      whi_odd = whi*SQRT(s_in/(hs1*(jshi-c1p5)))
      IF (jslo .ne. 1) THEN
         wlo_odd = wlo*SQRT(s_in/(hs1*(jslo-c1p5)))
      ELSE
         wlo_odd = 0
         whi_odd = SQRT(s_in/(hs1*(jshi-c1p5)))
      END IF

      WHERE (MOD(NINT(xm_nyq(:)),2) .eq. 0)
         wmins = wlo
         wplus = whi
      ELSEWHERE
         wmins = wlo_odd
         wplus = whi_odd
      END WHERE

      ipresent = 0
      lgsqrt = PRESENT(gsqrt)
      IF (lgsqrt) THEN
         gsqrt = 0 ;  ipresent = ipresent+1
         gmnc1 = wmins*gmnc(:,jslo) + wplus*gmnc(:,jshi)
         IF (lasym)
     1   gmns1 = wmins*gmns(:,jslo) + wplus*gmns(:,jshi)
      END IF
      lbsupu = PRESENT(bsupu)
      IF (lbsupu) THEN
         bsupu = 0 ;  ipresent = ipresent+1
         bsupumnc1 = wmins*bsupumnc(:,jslo) + wplus*bsupumnc(:,jshi)
         IF (lasym)
     1   bsupumns1 = wmins*bsupumns(:,jslo) + wplus*bsupumns(:,jshi)
      END IF
      lbsupv = PRESENT(bsupv)
      IF (lbsupv) THEN
         bsupv = 0 ;  ipresent = ipresent+1
         bsupvmnc1 = wmins*bsupvmnc(:,jslo) + wplus*bsupvmnc(:,jshi)
         IF (lasym)
     1   bsupvmns1 = wmins*bsupvmns(:,jslo) + wplus*bsupvmns(:,jshi)
      END IF

      IF (ipresent .eq. 0) RETURN

!
!     SETUP TRIG ARRAYS
!
      cosu = COS(u_in);   sinu = SIN(u_in)
      cosv = COS(v_in);   sinv = SIN(v_in)

      cosmu(0) = 1;    sinmu(0) = 0
      cosnv(0) = 1;    sinnv(0) = 0
      DO m = 1, mnyq
         cosmu(m) = cosmu(m-1)*cosu - sinmu(m-1)*sinu
         sinmu(m) = sinmu(m-1)*cosu + cosmu(m-1)*sinu
      END DO

      DO n = 1, nnyq
         cosnv(n) = cosnv(n-1)*cosv - sinnv(n-1)*sinv
         sinnv(n) = sinnv(n-1)*cosv + cosnv(n-1)*sinv
      END DO

!
!     COMPUTE GSQRT, ... IN REAL SPACE
!     tcosmn = cos(mu - nv);  tsinmn = sin(mu - nv)
!
      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tcosmn = cosmu(m)*cosnv(n1) + sgn*sinmu(m)*sinnv(n1)   
         IF (lgsqrt) gsqrt = gsqrt + gmnc1(mn)*tcosmn
         IF (lbsupu) bsupu = bsupu + bsupumnc1(mn)*tcosmn
         IF (lbsupv) bsupv = bsupv + bsupvmnc1(mn)*tcosmn
      END DO

      IF (.not.lasym) GOTO 1000

      DO mn = 1, mnmax_nyq
         m = NINT(xm_nyq(mn));  n = NINT(xn_nyq(mn))/nfp
         n1 = ABS(n);   sgn = SIGN(1,n)
         tsinmn = sinmu(m)*cosnv(n1) - sgn*cosmu(m)*sinnv(n1)
         IF (lgsqrt) gsqrt = gsqrt + gmns1(mn)*tsinmn
         IF (lbsupu) bsupu = bsupu + bsupumns1(mn)*tsinmn
         IF (lbsupv) bsupv = bsupv + bsupvmns1(mn)*tsinmn
      END DO

 1000 CONTINUE

      END SUBROUTINE tosuvspace

      SUBROUTINE LoadRZL
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER     :: rcc, rss, zsc, zcs, rsc, rcs, zcc, zss
      INTEGER     :: mpol1, mn, m, n, n1
      REAL(rprec) :: sgn
C-----------------------------------------------
!
!     Arrays must be stacked (and ns,ntor,mpol ordering imposed)
!     as coefficients of cos(mu)*cos(nv), etc
!     Only need R, Z components(not lambda, for now anyhow)
!      
      IF (ALLOCATED(rzl_local)) RETURN

      lthreed = ANY(NINT(xn) .ne. 0)
      mpol1 = mpol-1
      rcc = 1;  zsc = 1
      IF (.not.lasym) THEN
         IF (lthreed) THEN
            ntmax = 2
            rss = 2;  zcs = 2
         ELSE
            ntmax = 1
         END IF
      ELSE
         IF (lthreed) THEN
            ntmax = 4
            rss = 2;  rsc = 3;  rcs = 4
            zcs = 2;  zcc = 3;  zss = 4
         ELSE
            ntmax = 2
            rsc = 2;  zcc = 2
         END IF
      END IF

      zsc = 1+ntmax; zcs = zcs+ntmax; zcc = zcc+ntmax
      ALLOCATE(rzl_local(ns,0:ntor,0:mpol1,2*ntmax))
      rzl_local = 0

      DO mn = 1, mnmax
         m = NINT(xm(mn));  n = NINT(xn(mn))/nfp; n1 = ABS(n)
         sgn = SIGN(1, n)
         rzl_local(:,n1,m,rcc) = rzl_local(:,n1,m,rcc) + rmnc(mn,:)
         rzl_local(:,n1,m,zsc) = rzl_local(:,n1,m,zsc) + zmns(mn,:)
         IF (lthreed) THEN
            rzl_local(:,n1,m,rss) = rzl_local(:,n1,m,rss) 
     1                            + sgn*rmnc(mn,:)
            rzl_local(:,n1,m,zcs) = rzl_local(:,n1,m,zcs)
     1                            - sgn*zmns(mn,:)
         END IF
         IF (lasym) THEN
            rzl_local(:,n1,m,rsc) = rzl_local(:,n1,m,rsc) 
     1                            + rmns(mn,:)
            rzl_local(:,n1,m,zcc) = rzl_local(:,n1,m,zcc) 
     1                            + zmnc(mn,:)
            IF (lthreed) THEN
                rzl_local(:,n1,m,rcs) = rzl_local(:,n1,m,rcs)
     1                                - sgn*rmns(mn,:)
                rzl_local(:,n1,m,zss) = rzl_local(:,n1,m,zss) 
     1                                + sgn*zmnc(mn,:)
            END IF
         END IF
      END DO

      END SUBROUTINE LoadRZL

      END MODULE read_wout_mod

