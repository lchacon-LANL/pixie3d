      MODULE efit

      use io

      use xdraw_io

      use grid

      use local_BCS_variables, ONLY: default_B_BCs,default_A_BCs,bcond,setMGBC,order_bc

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
      real(8) :: r_max,r_min,z_max,z_min,LL,iLL,psisgn=1d0,e_bp=0d0

      logical :: short_efit_file=.true.,efit_dbg=.false.

      !Physical constants
      real(8),parameter ::  mu_0 = (4*pi)*1d-7 &!Magnetic permeability (SI units)
                           ,eps_0= 8.8542d-12  &!Permittivity of free space
                           ,n_0  = 1d20        &!m^-3
                           ,kb   = 1.6022d-19  &!J/eV
                           ,mi   = 1.6726d-27  &!kg
                           ,me   = 9.1094d-31  &!kg
                           ,qe   = 1.6022d-19   !C
      
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
                 ,igrid,igrid,igrid,ig,jg,kg,x1,y1,z1)

        RR = sqrt(x1**2+y1**2)
        ZZ = z1 !+ zmaxis*iLL

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
        real(8) :: chk_simag,chk_sibry

        integer :: i,j
        !Begin program
        neqdsk = find_unit(123)
        
        open (unit=neqdsk,file=trim(efit_file),status="old")
   
        read (neqdsk,2000,IOSTAT=istat) (case(i),i=1,6),idum,nw,nh
        if (istat /= 0) then
          write (*,*) "EFIT read error 1"
          return
        endif
        
        allocate(psirz(nw,nh),fpol(nw),pres(nw),ffprim(nw),pprime(nw),qpsi(nw) &
                ,pressw(nw),pwprim(nw),dmion(nw),rhovn(nw))
   
        read (neqdsk,2020,IOSTAT=istat) rdim,zdim,rcentr,rleft,zmid
        if (istat /= 0) then
          write (*,*) "EFIT read error 2"
          return
        endif

        read (neqdsk,2020,IOSTAT=istat) rmaxis,zmaxis,simag,sibry,bcentr
        if (istat /= 0) then
          write (*,*) "EFIT read error 3"
          return
        endif
        
        read (neqdsk,2020,IOSTAT=istat) current,chk_simag,xdum,rmaxis,xdum
        if (istat /= 0) then
          write (*,*) "EFIT read error 4"
          return
        endif

        read (neqdsk,2020,IOSTAT=istat) zmaxis,xdum,chk_sibry,xdum,xdum
        if (istat /= 0) then
          write (*,*) "EFIT read error 5"
          return
        endif

        if (chk_sibry /= sibry.or.chk_simag/=simag) then
          write (*,*) "Inconsistent psi limits; FIX in EQU file"
          istat = 1
          return
        endif

        read (neqdsk,2020,IOSTAT=istat) (fpol(i),i=1,nw)
        if (istat /= 0) then
          write (*,*) "EFIT read error 6"
          return
        endif

        read (neqdsk,2020,IOSTAT=istat) (pres(i),i=1,nw)
        if (istat /= 0) then
          write (*,*) "EFIT read error 7"
          return
        endif

        read (neqdsk,2020,IOSTAT=istat) (ffprim(i),i=1,nw)
        if (istat /= 0) then
          write (*,*) "EFIT read error 8"
          return
        endif

        read (neqdsk,2020,IOSTAT=istat) (pprime(i),i=1,nw)
        if (istat /= 0) then
          write (*,*) "EFIT read error 9"
          return
        endif

        read (neqdsk,2020,IOSTAT=istat) ((psirz(i,j),i=1,nw),j=1,nh)
        if (istat /= 0) then
          write (*,*) "EFIT read error 10"
          return
        endif

        read (neqdsk,2020,IOSTAT=istat) (qpsi(i),i=1,nw)
        if (istat /= 0) then
          write (*,*) "EFIT read error 11"
          return
        endif

        read (neqdsk,2022,IOSTAT=istat) nbbbs,limitr
        if (istat /= 0) then
          write (*,*) "EFIT read error 12"
          return
        endif
        
        allocate(rbbbs(nbbbs),zbbbs(nbbbs),rlim(limitr),zlim(limitr))

        read (neqdsk,2020,IOSTAT=istat) (rbbbs(i),zbbbs(i),i=1,nbbbs)
        if (istat /= 0) then
          write (*,*) "EFIT read error 13"
          return
        endif

        read (neqdsk,2020,IOSTAT=istat) (rlim(i),zlim(i),i=1,limitr)
        if (istat /= 0) then
          write (*,*) "EFIT read error 14"
          return
        endif

        call clean_sptrx(nbbbs,rbbbs,zbbbs)
        
        if (.not.short_efit_file) then
          read (neqdsk,2024,IOSTAT=istat) kvtor,rvtor,nmass
          if (istat /= 0) then
            write (*,*) "EFIT read error 15"
            istat = 0 !Incomplete EFIT file; not a terminating error
            short_efit_file = .true.
            return
          endif

          if (kvtor.gt.0) then
            read (neqdsk,2020,IOSTAT=istat) (pressw(i),i=1,nw)
            if (istat /= 0) then
              write (*,*) "EFIT read error 16"
              return
            endif

            read (neqdsk,2020,IOSTAT=istat) (pwprim(i),i=1,nw)
            if (istat /= 0) then
              write (*,*) "EFIT read error 17"
              return
            endif
          endif
        
          if (nmass.gt.0) then
            read (neqdsk,2020,IOSTAT=istat) (dmion(i),i=1,nw)
            if (istat /= 0) then
              write (*,*) "EFIT read error 18"
              return
            endif
          endif
        
          read (neqdsk,2020,IOSTAT=istat) (rhovn(i),i=1,nw)
          if (istat /= 0) then
            write (*,*) "EFIT read error 19"
            return
          endif
        endif

        close (neqdsk)
        
        if (my_rank == 0) call read_chk(0)
   
2000    format (6a8,3i4)
2020    format (5e16.9)
2022    format (2i5)
2024    format (i5,e16.9,i5)

      contains

!       read_chk
!       #########################################################################
        subroutine read_chk(flag)

          integer :: flag
          
          !Plot Psi(R,Z)
          call createDrawInCfile(1,"efit_psi.bin",'Solution','t','x','y' &
                              ,(/'Psi'/),'-c -X0 -L57','drawefit_psi.in')

          open(unit=110,file="efit_psi.bin",form='unformatted',status='unknown')

          call contour(psirz,nw,nh,rleft,rleft+rdim,zmid-0.5*zdim,zmid+0.5*zdim &
                      ,0,110)
          close(110)

          !Plot flux boundary (gnuplot)
          open(unit=110,file="efit_psi_bdry.txt",status='unknown')
          do i=1,nbbbs
             write(110,*) rbbbs(i),zbbbs(i)
          enddo
          close(110)

          !Plot limiter boundary (gnuplot)
          open(unit=110,file="efit_limiter_bdry.txt",status='unknown')
          do i=1,limitr
             write(110,*) rlim(i),zlim(i)
          enddo
          close(110)

          !Plot q profile (gnuplot)
          open(unit=110,file="efit_q.txt",status='unknown')
          do i=1,nw
             write(110,*) 1d0*i/nw,abs(qpsi(i))
          enddo
          close(110)

          !Plot pressure profile (gnuplot)
          open(unit=110,file="efit_prs.txt",status='unknown')
          do i=1,nw
             write(110,*) 1d0*i/nw,abs(pres(i))
          enddo
          close(110)

          if (flag == 0) return
          
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

          if (.not.short_efit_file) then
            write(*,2024) kvtor,rvtor,nmass

            if (kvtor.gt.0) then
              write(*,2020) (pressw(i),i=1,nw)
              write(*,2020) (pwprim(i),i=1,nw)
            endif
            if (nmass.gt.0) then
              write(*,2020) (dmion(i),i=1,nw)
            endif
            write(*,2020) (rhovn(i),i=1,nw)
          endif
          
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
        
        if (my_rank == 0) then
          write (*,*) "Magnetic axis (R,Z)",rmaxis,zmaxis
          write (*,*) "Psi at magnetic axis",simag
          write (*,*) "Psi at magnetic bdry",sibry
        endif
        
!     FIND DOMAIN LIMITS AND SETUP GEOMETRY

        !Limiter box
        z_max = maxval(zlim)
        z_min = minval(zlim)
        r_max = maxval(rlim)
        r_min = minval(rlim)

        LL = 0.5*(r_max-r_min) !Dimensionalization length

        if (LL == 0d0) then
           LL = 1d0
           r_max = rleft+rdim
           r_min = rleft
           z_max = zmid + zdim/2d0
           z_min = zmid - zdim/2d0
        endif
        
        iLL = 1d0/LL
        
        if (my_rank == 0) then
          write (*,*) "Minor rad.  a (m)=",LL
          write (*,*) "Limiter R_max (m)=",r_max
          write (*,*) "Limiter R_min (m)=",r_min
          write (*,*) "Limiter Z_max (m)=",z_max
          write (*,*) "Limiter Z_min (m)=",z_min
        endif
        
        r_max = r_max*iLL
        r_min = r_min*iLL
        z_max = z_max*iLL
        z_min = z_min*iLL

        rdim  = rdim *iLL
        rleft = rleft*iLL

        zdim = zdim*iLL
        zmid = zmid*iLL

        if (my_rank == 0) then
          write (*,*) "Psi box R_max/a=",rleft+rdim
          write (*,*) "Psi box R_min/a=",rleft
          write (*,*) "Psi box Z_max/a=",zmid + 0.5*zdim
          write (*,*) "Psi box Z_min/a=",zmid - 0.5*zdim
        endif
        
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

        if (efit_dbg.and.my_rank==0) write (*,*) "Splining psirz..."

        psisgn = (sibry-simag)/abs(sibry-simag)
        psirz = (psirz-simag)*psisgn

        sibry = abs(sibry-simag)
        simag = 0d0
        
        if (efit_dbg.and.my_rank==0) write (*,*) "After Psi limits",simag,sibry

        call db2ink(xs,nxs,zs,nzs,psirz,nxs,kx,kz,tx,tz,psi_coef,work,flg)

        if (efit_dbg.and.my_rank==0) write (*,*) "Finished!"

!     SPLINE 1D PSI-DEPENDENT ARRAYS

        if (efit_dbg.and.my_rank==0) write (*,*) "Splining psi functions..."

        allocate(ps(nxs),stat=alloc_stat)

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

        allocate(q((2*kx-1)*nxs),stat=alloc_stat)

        call dbknot(ps,nxs,kx,tps)
        call dbintk(ps,fpol  ,tps,nxs,kx,fpol_coef  ,q,work)
        call dbintk(ps,pres  ,tps,nxs,kx,pres_coef  ,q,work)
        call dbintk(ps,ffprim,tps,nxs,kx,ffprim_coef,q,work)
        call dbintk(ps,pprime,tps,nxs,kx,pprime_coef,q,work)
        call dbintk(ps,qpsi  ,tps,nxs,kx,qpsi_coef  ,q,work)

        if (efit_dbg) write (*,*) "Finished!"

        deallocate(q)
        
      END SUBROUTINE efit_init

!     efit_cleanup
!     ################################################################################
      SUBROUTINE efit_cleanup

        IMPLICIT NONE

        INTEGER :: istat

        deallocate(psirz,fpol,pres,ffprim,pprime,qpsi,rbbbs,zbbbs,rlim,zlim &
                  ,pressw,pwprim,dmion,rhovn,STAT=istat)
   
        DEALLOCATE(work,tx,tz,tps,psi_coef,fpol_coef,ffprim_coef,pprime_coef &
                  ,qpsi_coef,pres_coef,stat=istat)
        
        deallocate(xs,zs,ps)

        IF (istat .ne. 0) STOP 'Deallocation error in efit_cleanup'

      END SUBROUTINE efit_cleanup

!     inside_sptrx
!     #################################################################
      function inside_sptrx(R,Z) result(inside)

!     -----------------------------------------------------------------
!     Determines whether a point is inside separatrix.
!     -----------------------------------------------------------------

        implicit none

!     Call variables

        real(8) :: R,Z
        logical :: inside
        
!     Local variables

        integer :: npt,ipt
        real(8) :: rz(2),dth,theta,a(2),b(2)

!     Begin program
        
        rz(1) = R/iLL
        rz(2) = Z/iLL

        npt = nbbbs

        !Order N algorithm based on angular integration
        !(approximates shape by a polygon)
        !Works for arbitrary curve shape
        theta = 0d0
        b = rz-(/rbbbs(1),zbbbs(1)/)
        do ipt=2,npt
           a = b
           b = rz-(/rbbbs(ipt),zbbbs(ipt)/)
           dth = atan2((a(1)*b(2)-a(2)*b(1)),dot_product(a,b))
           theta = theta + dth
        enddo
        a = b
        b = rz-(/rbbbs(1),zbbbs(1)/)
        dth = atan2((a(1)*b(2)-a(2)*b(1)),dot_product(a,b))
        theta = theta + dth
        
!!        write (*,*) npt,theta/(2*pi)
        inside = (abs(theta) > 0.98*2*pi)
!!        if (.not.inside) write (*,*) "diag",rz
        
      end function inside_sptrx

!     inside_lmtr
!     #################################################################
      function inside_lmtr(R,Z) result(inside)

!     -----------------------------------------------------------------
!     Determines whether a point is inside limiter.
!     -----------------------------------------------------------------

        implicit none

!     Call variables

        real(8) :: R,Z
        logical :: inside
        
!     Local variables

        integer :: npt,ipt
        real(8) :: rz(2),dth,theta,a(2),b(2)

!     Begin program
        
        rz(1) = R!/iLL
        rz(2) = Z!/iLL

        npt = size(rlim)

        !Order N algorithm based on angular integration
        !(approximates shape by a polygon)
        !Works for arbitrary curve shape
        theta = 0d0
        b = rz-(/rlim(1),zlim(1)/)
        do ipt=2,npt
           a = b
           b = rz-(/rlim(ipt),zlim(ipt)/)
           dth = atan2((a(1)*b(2)-a(2)*b(1)),dot_product(a,b))
           theta = theta + dth
        enddo
        a = b
        b = rz-(/rlim(1),zlim(1)/)
        dth = atan2((a(1)*b(2)-a(2)*b(1)),dot_product(a,b))
        theta = theta + dth

!!$        write (*,*) npt,theta/(2*pi)
        inside = (abs(theta) > 0.98*2*pi)

      end function inside_lmtr

!     clean_sptrx
!     #################################################################
      subroutine clean_sptrx(npt,rb,zb)

!     -----------------------------------------------------------------
!     Cleans up separatrix to ensure points are inside limiter.
!     -----------------------------------------------------------------
        
        implicit none

        integer :: npt
        real(8) :: rb(:),zb(:)

!       Local variables

        integer :: np,i
        real(8) :: mag,iLL
        real(8), allocatable,dimension(:) :: rl,zl

!       Begin program

        allocate(rl(npt),zl(npt))
        
        !First pass (select inside limiter)
        np = 0
        do i=1,npt
           if (inside_lmtr(rb(i),zb(i))) then
              np = np + 1
              zl(np) = zb(i)
              rl(np) = rb(i)
           endif
        enddo

        npt = np
        rb(1:np) = rl(1:np)
        zb(1:np) = zl(1:np)

        !Second pass (detect jumps in curve)
        iLL = 1d0/(maxval(rlim)-minval(rlim))
        mag = 0d0
        np = 0
        do i=1,npt-1
           np = np + 1
           mag = sqrt((rb(i+1)-rb(i))**2+(zb(i+1)-zb(i))**2)
           if (mag*iLL > 0.3.and.np < npt/2) then 
              np = 0 !large jump across points; line too short: cut out
           else
              zl(np) = zb(i)
              rl(np) = rb(i)
           endif
        enddo

        !Store new separatrix
        npt = np
        rb(1:np) = rl(1:np)
        zb(1:np) = zl(1:np)

        if (np < npt) then
           rb(np+1:)=0d0
           zb(np+1:)=0d0
        endif
        
        deallocate(rl,zl)
        
      end subroutine clean_sptrx
      
      END MODULE efit

!     efit_map
!     #################################################################
      subroutine efit_map(equ_file,g_def,coord)

!     -----------------------------------------------------------------
!     Give Cartesian coordinates of each logical mesh point at grid
!     level (igrid).
!     -----------------------------------------------------------------

      use efit

      use equilibrium

      implicit none

!     Input variables

      type(grid_mg_def),pointer :: g_def
      character(*) :: equ_file,coord

!     Local variables

      integer :: nxg,nyg,nzg
      real(8) :: scale

!     Begin program

      if (my_rank == 0) then
        write (*,*)
        write (*,'(a)') ' Creating EFIT map...'
        write (*,*)
        write (*,'(a)') ' EFIT Domain dimensions:'
        write (*,*)
      endif
      
      !Read equilibrium file

      call efit_init(equ_file)

      !Setup geometry

      select case(trim(coord))
      case('tor')

        coords = coord

        gparams(1) = 0.5*((1-gparams(1))*r_max+(1+gparams(1))*r_min)   !Major radius (biased)
!!        gparams(2) = scale               !Horizontal elliptical radius
        if (gparams(3) > 0d0) then
           gparams(3) = z_max*gparams(3)    !Vertical elliptical radius
        else
           gparams(3) = gparams(2)
        endif

      case('sh3')  !See grid_anal_map_mod.F for gparams legend

        coords = coord

        !Origin R-coord (set to magnetic axis if not provided)
        if (gparams(1) == 0d0) then
           gparams(1) = rmaxis*iLL         !Magnetic axis R-coord
        else
           if (gparams(1) > 0d0) then      !Do NOT shift to magnetic axis
              gparams(1) = gparams(1)*iLL  !Specified axis R-coord; backward compatibility
           else
              gparams(8) = -(abs(gparams(1))-rmaxis)*iLL  !R-coordinate bdry shift
              gparams(1) = rmaxis*iLL      !Magnetic axis R-coord
           endif
        endif

        !Origin Z-coord (set to magnetic axis if not provided)
        if (gparams(2) == 0d0) then
           gparams(2) = zmaxis*iLL         !Magnetic axis Z-coord
        else
           if (gparams(2) > 0d0) then      !Do NOT shift to magnetic axis
              gparams(2) = gparams(2)*iLL  !Specified axis Z-coord; backward compatibility
           else
              gparams(9) = -(abs(gparams(2))-zmaxis)*iLL  !Z-coordinate bdry shift
              gparams(2) = zmaxis*iLL      !Magnetic axis Z-coord
           endif
        endif

        !Minor radius (set if not provided)
        if (gparams(3) == 0d0) then
          gparams(3) = 1d0
        else
          gparams(3) = gparams(3)*iLL
        endif

        !Elongation of bdry (estimate if not provided)
        if (gparams(4) == 0d0) then
          gparams(4) = (maxval(zbbbs)-minval(zbbbs))/(maxval(rbbbs)-minval(rbbbs))
        endif
        
        ! gparams(5),delta, provided in input deck
        ! gparams(6), zeta, provided in input deck

        !Elongation at SP (estimate if not provided)
        if (gparams(7) == 0d0) then
          if (gparams(8)==0d0.or.gparams(9)==0d0) then
            gparams(7) = gparams(4)  !Origin NOT at magnetic axis: use bdry elongation
          else
            gparams(7) = zdim/rdim   !Estimate elongation of flux surfaces near mag-axis (crudely)
          endif
        endif

      case('tsq')

        coords = coord

        !Find box dimensions based on limiter
        if (gparams(1) == 0d0) gparams(1) = 0.5*(r_max+r_min)  !Domain center R-coord
        if (gparams(2) == 0d0) gparams(2) = 0.5*(z_max+z_min)  !Domain center Z-coord
        if (gparams(3) == 0d0) gparams(3) = r_max-r_min        !Horizontal dimension
        if (gparams(4) == 0d0) gparams(4) = z_max-z_min        !Vertical   dimension

      case default

        call pstop("efit_map","Coordinates unknown")

      end select
      
      if (my_rank == 0) then
        write (*,*) "R/a =",0.5*(r_max+r_min)
        write (*,*) "Z/a =",z_max
      endif
        
      nxg = g_def%nglx
      nyg = g_def%ngly
      nzg = g_def%nglz

      call destroyGrid(g_def)
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

!!$        logical :: dcon

!     Local variables

        integer :: i,j,k,ig,jg,kg,istat,inbv,bcs(6,3),nxg,nyg,nzg
        real(8) :: max_prs,RR,ZZ,BB0,iB0,ip0,ipsi0,x1,y1,z1

        real(8),allocatable, dimension(:,:,:) :: bsub3,psi
        real(8),allocatable, dimension(:,:,:,:) :: aa

        real(8)  :: db2val,dbvalu
        external :: db2val,dbvalu

        real(8) :: Zi,mi_mp,log_lamb,v_A,t_A,T_max,tau_e,tau_i,omega_ce,omega_ci &
                  ,n0_cgs,B0_cgs
        
!     Begin program

        if (my_rank == 0) then
           write (*,*)
           write (*,*) 'Reading EFIT equilibrium...'
           write (*,*)
        endif

!     Get GLOBAL limits

        nxg = gv%gparams%nxgl(igrid)
        nyg = gv%gparams%nygl(igrid)
        nzg = gv%gparams%nzgl(igrid)

!     Interpolate poloidal flux and associated qtys

        allocate(psi  (0:nx+1,0:ny+1,0:nz+1) &
                ,bsub3(0:nx+1,0:ny+1,0:nz+1) &
                ,aa   (0:nx+1,0:ny+1,0:nz+1,3))

        inbv = 1

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1

              call find_RZ(gv%gparams,igrid,i,j,k,RR,ZZ)

              !Interpolate poloidal flux
              psi(i,j,k) = db2val(RR,ZZ,0,0,tx,tz,nxs,nzs  &
                                 ,kx,kz,psi_coef,work)

              !Interpolate pressure
              if (  (nbbbs >0.and.inside_sptrx(RR,ZZ)) &
                .or.(nbbbs==0.and.(psi(i,j,k)<=sibry))) then
                if (psi(i,j,k) > sibry) then
                  call pstop("efit_equ","Error in separatrix boundary")
                endif
                prs  (i,j,k) = dbvalu(tps,pres_coef,nxs,kx,0,psi(i,j,k),inbv,work)
                bsub3(i,j,k) = dbvalu(tps,fpol_coef,nxs,kx,0,psi(i,j,k),inbv,work)
              else
                prs  (i,j,k) = dbvalu(tps,pres_coef,nxs,kx,0,sibry,inbv,work)
                bsub3(i,j,k) = dbvalu(tps,fpol_coef,nxs,kx,0,sibry,inbv,work)
              endif

              prs(i,j,k) = abs(prs(i,j,k))

            enddo
          enddo
        enddo

        !BCs
        aa(:,:,:,1:2) = 0d0
        aa(:,:,:,3  ) = psi
        
        call default_A_BCs(bcs)
        where (bcs == -EQU) bcs = -DEF  !Do nothing to tangential components

        call setMGBC(gv%gparams,0,3,nx,ny,nz,igrid,aa,bcs &
     &              ,icomp=(/IAX/),is_vec=.true.           &
     &              ,is_cnv=.false.,iorder=2)

!     Normalization constants for magnetic field (toroidal at magnetic axis) and pressure
        
        BB0 = abs(dbvalu(tps,fpol_coef,nxs,kx,0,simag,inbv,work))/rmaxis

        iB0 = 1d0/BB0
        ip0 = mu_0*iB0*iB0  !(B0^2/mu_0)^(-1)

        ipsi0 = iB0*iLL*(2*pi)**(-e_bp)
!!        ipsi0 = iB0*iLL*(2*pi)**(-1d0)
        
!     Normalized pressure

        prs = prs*ip0 + 1d-4 !Add pressure floor

!     Find magnetic field components

        !Poloidal B (cnv)
        aa = aa*ipsi0*iLL  !Cov A (R*A_T)
        bb = curl(gv%gparams,igrid,aa)

        !Toroidal B
        bb(:,:,:,3) = bsub3*ipsi0   !Cov B (R*B_T)

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
        call default_B_BCs(bcs)
        where (bcs == -NEU) bcs = -DEF  !Do nothing to tangential components

        call setMGBC(gv%gparams,0,3,nx,ny,nz,igrid,bb,bcs &
     &              ,icomp=(/IBX/),is_vec=.true.          &
     &              ,is_cnv=.true.,iorder=2)

!     Find normalized density

        max_prs = maxval(prs)
        max_prs = pmax(max_prs)
        
        if (max_prs == 0d0.or.(gam==0d0)) then
          rho = 1d0
        else
          rho = (abs(prs/max_prs)**(1d0/gam))  !Forces positive density
        endif

!     Find normalization constants

        n0_cgs = n_0*1d-6  !1/cm^3
        B0_cgs = BB0*1d4   !Gauss
        Zi = 1d0
        mi_mp = 1d0
        log_lamb = 1d1
        
        v_A = BB0/sqrt(mi*n_0*mu_0)  !m/s
        t_A = 1d0/(v_A*iLL) !s
        T_max = 0.5*max_prs/(kb*ip0*n_0) !eV
!!$        tau_e = (4*pi*eps_0)**(2)*3*sqrt(me)*(kb*T_max)**1.5/(4*sqrt(2*pi)*n_0*log_lamb*qe**4) !s (SI)
!!$        write (*,*) tau_e/(3.44d5*(T_max**1.5)/(n0_cgs*log_lamb))
        tau_e = 3.44d5*(T_max**1.5)/(n0_cgs*log_lamb) !s (cgs)
        tau_i = 2.08d7*(T_max**1.5)/(Zi**4*n0_cgs*log_lamb)*sqrt(mi_mp) !s (cgs)
        omega_ce = qe*BB0/me    !1/s (SI)
        omega_ci = qe*Zi*BB0/mi !1/s (SI)
        
        if (my_rank == 0) then
          write (*,'(a)') " Normalization constants:"
          write (*,*)
          write (*,'(a,1pe14.7,a)') " L_0 = ",1d0/iLL," m"
          write (*,'(a,1pe14.7,a)') " n_0 = ",n_0," 1/m^3"
          write (*,'(a,1pe14.7,a)') " B_0 = ",BB0," T"
          write (*,'(a,1pe14.7,a)') " p_0 = ",1d0/ip0," N/m^2"
          write (*,'(a,1pe14.7,a)') " v_A = ",v_A," m/s"
          write (*,'(a,1pe14.7,a)') " t_A = ",t_A," s"
          write (*,'(a,1pe14.7,a)') " Te max= ",T_max*1d-3," keV "
          write (*,*)
          write (*,'(a)') " Normalized parameters:"
          write (*,*)
          write (*,'(a,1pe14.7)') " beta       = ",max_prs
          write (*,'(a,1pe14.7)') " chi_par e  = ",3.2*kb*T_max*tau_e/(me*v_A)*iLL  !electrons
!!          write (*,'(a,1pe14.7)') " chi_perp e = ",4.7*kb*T_max/(me*omega_ce**2*tau_e*v_A)*iLL  !e
          write (*,'(a,1pe14.7)') " chi_perp i = ",2.0*kb*T_max/(mi*omega_ci**2*tau_i*v_A)*iLL  !ions
          write (*,'(a,1pe14.7)') " eta        = ",me/(n_0*qe*qe*tau_e*v_A*mu_0)*iLL
          write (*,*)
          write (*,'(a)') " Other relevant parameters:"
          write (*,*)
          write (*,'(a,1pe14.7,a)') " tau_ee   =",tau_e," s^(-1)"
          write (*,'(a,1pe14.7,a)') " tau_ii   =",tau_i," s^(-1)"
          write (*,'(a,1pe14.7,a)') " omega_ci =",omega_ci," rad/s"
          write (*,'(a,1pe14.7,a)') " omega_ce =",omega_ce," rad/s"
        endif

!     Check EFIT qtys

        call efit_chk(.false.)

!     Check GS equilibrium

        call GS_chk(.true.)

!     Free work space

        DEALLOCATE(psi,bsub3,aa,stat=istat)
        
        call efit_cleanup

!     End program

      contains

!     efit_chk
!     #################################################################
      subroutine efit_chk(dump)

        implicit none
        
        logical :: dump

        real(8) :: psi_mag,psi_bry
        
        if (dump) then
           !GNUPLOT Psi(R,Z), P(R,Z), B_3(R,Z)
           open(unit=110,file="efit_qtys.txt",status='unknown')

           k = 1
           do j=0,ny+1
              do i=0,nx+1
                 call find_RZ(gv%gparams,igrid,i,j,k,RR,ZZ)
                 write (110,*) RR,ZZ,psi(i,j,k),prs(i,j,k),bsub3(i,j,k)
              enddo
              write (110,*)
           enddo
           
           close(110)

           !XDRAW plot Psi, P, I=R*Bt
           call createDrawInCfile(6,"efit_qtys.bin",'Solution','t','x','y' &
                ,(/'Psi','Prs','B_3','B^1','B^2','B^3'/),'-c -X0 -L57'     &
                ,'drawefit_qtys.in')

           open(unit=110,file="efit_qtys.bin",form='unformatted',status='unknown')

           call contour(psi  (:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,0,110)
           call contour(prs  (:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
           call contour(bsub3(:,:,1),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
           call contour(bb   (:,:,1,1),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
           call contour(bb   (:,:,1,2),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)
           call contour(bb   (:,:,1,3),nx+2,ny+2,xmin,xmax,ymin,ymax,1,110)

           close(110)
        endif

        if (my_rank == 0) then
           write (*,*) "Magnetic axis normalized coords=",rmaxis*iLL,zmaxis*iLL
           psi_mag = db2val(rmaxis*iLL,zmaxis*iLL,0,0,tx,tz,nxs,nzs  &
     &                     ,kx,kz,psi_coef,work)

           if (nbbbs > 0) then
              psi_bry = db2val(rbbbs(1)*iLL,zbbbs(1)*iLL,0,0,tx,tz,nxs,nzs  &
     &                     ,kx,kz,psi_coef,work)
           else
              psi_bry = sibry
           endif
           
!!           write (*,*) simag,sibry,psi_mag,psi_bry
           write (*,*) "delta psi_mag =",(psi_mag-simag)/max(simag,sibry)
           write (*,*) "delta psi_bry =",(psi_bry-sibry)/max(simag,sibry)
        endif

        if (dump) STOP "EFIT CHECK"
        
      end subroutine efit_chk
      
!     GS_chk
!     #################################################################
      subroutine GS_chk(dump)

        implicit none

        logical :: dump
        
        real(8),allocatable, dimension(:,:,:,:) :: vv,jj
        real(8),allocatable, dimension(:,:,:)   :: mag1,mag2,mag3,mag1g,mag2g,mag3g

        real(8) :: gs_error,norm1,norm2
        
        allocate(vv  (0:nx+1,0:ny+1,0:nz+1,3) &
                ,jj  (0:nx+1,0:ny+1,0:nz+1,3) &
                ,mag1(0:nx+1,0:ny+1,0:nz+1)   &
                ,mag2(0:nx+1,0:ny+1,0:nz+1)   &
                ,mag3(0:nx+1,0:ny+1,0:nz+1))

        vv = XformToCov(gv%gparams,igrid,bb)
        jj = curl(gv%gparams,igrid,vv)

        vv = crossProduct(gv%gparams,igrid,jj,bb,.true.)
        mag1 = vectorNorm(gv%gparams,igrid,vv,.true.)
        norm1 = integral(gv%gparams,igrid,mag1,average=.true.)

        jj = grad(gv%gparams,igrid,prs)
        mag2 = vectorNorm(gv%gparams,igrid,jj,.true.)
        
        vv = vv - jj
        mag3 = vectorNorm(gv%gparams,igrid,vv,.true.)
        gs_error = integral(gv%gparams,igrid,mag3,average=.true.)

        if (my_rank == 0) then
          write (*,*)
          write (*,*) "Checking EFIT GS error..."
          write (*,*)
          write (*,*) "Relative GS error =",sqrt(gs_error/norm1)
        endif
        
        !Diagnostics plots
        if (dump) then
           allocate(mag1g(0:nxg+1,0:nyg+1,0:nzg+1))
           allocate(mag2g(0:nxg+1,0:nyg+1,0:nzg+1))
           allocate(mag3g(0:nxg+1,0:nyg+1,0:nzg+1))

           mag1 = sqrt(mag1)
           mag2 = sqrt(mag2)
           mag3 = sqrt(mag3)
        
           call find_global_nobc(mag1 (1:nx ,1:ny ,1:nz ) &
     &                          ,mag1g(1:nxg,1:nyg,1:nzg))
           call find_global_nobc(mag2 (1:nx ,1:ny ,1:nz ) &
     &                          ,mag2g(1:nxg,1:nyg,1:nzg))
           call find_global_nobc(mag3 (1:nx ,1:ny ,1:nz ) &
     &                          ,mag3g(1:nxg,1:nyg,1:nzg))

           if (my_rank == 0) then
              call createDrawInCfile(3,"efit_gschk.bin",'Solution','t','x','y' &
                            ,(/'JxB ','grdP','GS E'/),'-c -X0 -L57','drawefit_gschk.in')

              open(unit=110,file="efit_gschk.bin",form='unformatted',status='unknown')

              call contour(mag1g(1:nxg,1:nyg,1),nxg,nyg,xmin,xmax,ymin,ymax,0,110)
              call contour(mag2g(1:nxg,1:nyg,1),nxg,nyg,xmin,xmax,ymin,ymax,1,110)
              call contour(mag3g(1:nxg,1:nyg,1),nxg,nyg,xmin,xmax,ymin,ymax,1,110)

              close(110)
           endif
        endif
        
        DEALLOCATE(vv,jj,mag1,mag2,mag3,mag1g,mag2g,mag3g,stat=istat)
 
      end subroutine GS_chk

!!$!     dcon_dump
!!$!     #################################################################
!!$      subroutine dcon_dump
!!$
!!$        implicit none
!!$
!!$        real(8),allocatable, dimension(:) :: pflx,ff_i,q_i
!!$
!!$        REAL(rprec), DIMENSION(:,:,:), ALLOCATABLE :: rr,zz
!!$      
!!$        integer :: udcon=1111
!!$
!!$!     DCON dump (NOT WORKING)
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
!!$      end subroutine dcon_dump

      end subroutine efit_equ
