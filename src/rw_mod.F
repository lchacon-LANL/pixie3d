c module rw_bc
c ####################################################################
      module rw_bc

        use local_BCS_variables

        use var_io, ONLY: ilevel

        use app_iosetup, ONLY: test,critical
        
        use rw_bim_mod, only: rw_bim_symm_init

        !Resistive wall vars
        logical :: rw_bcs(6)      = .false.
     $            ,rw_fft_tupdate = .false. ! Advance Fourier RW modes in time
     $            ,rw_implicit_bc = .false. ! Use implicit RW BCs                 
     $            ,rw_exp_ssp_bc  = .true.  ! Use explicit SSP RW BCs             
     $            ,rw_exp_pc_bc   = .false. ! Use explicit PC BCs
     $            ,rw_exp_fix     = .true.  ! Use instability fix for small tau_w

        real(8) :: tau_w=1d30,r_v=1d0  !Ideal wall by default

        real(8),pointer,dimension(:,:,:,:) :: E_bc

        logical,private :: debug_io=.false.

        real(8),private :: dr_w

        integer,private :: m_max = 0

      contains

c     rw_allocBCVariables
c     ###################################################################
      subroutine rw_allocBCVariables

c     -------------------------------------------------------------------
c     Allocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: alloc_stat

        integer :: nx,ny,nz,nzg
      
c     Begin program

        nx = gv%gparams%nxv(1) 
        ny = gv%gparams%nyv(1)
        nz = gv%gparams%nzv(1)

        nzg = gv%gparams%nzgl(1)

        allocate(E_bc(0:nx+1,0:ny+1,0:nz+1,3))

        E_bc = 0d0

        select case(coords)
        case('tor','sha','sh2')
          if (nzd == 1) then
            call rw_bim_symm_init(gv%gparams,(ilevel>0)
     $                           ,(test.and.coords=="tor"))
            if (test.and.coords=="tor") stop
          else
            call critical('rw_alloc'
     $      ,'VACUUM RW response not implemented for non-symm tor geom')
          endif
        end select
      
c     End program

      end subroutine rw_allocBCVariables

c     rw_deallocBCVariables
c     ###################################################################
      subroutine rw_deallocBCVariables

c     -------------------------------------------------------------------
c     Deallocates preconditioner variables.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

c     Local variables

        integer :: alloc_stat

c     Begin program

        deallocate(E_bc,stat=alloc_stat)

c     End program

      end subroutine rw_deallocBCVariables

c     rw_Ebc
c     #################################################################
      subroutine rw_Ebc(g_def,igr,bcnv,E_0,extrp_b,newtime,dbg,bcnvp)

c     -----------------------------------------------------------------
c     Finds tangential cov electric field components according to
c     resistive wall BCs.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr
      
      real(8),dimension(0:,0:,0:,:) :: bcnv
      real(8),dimension(0:,0:,0:,:) :: E_0
      real(8),dimension(0:,0:,0:,:),optional :: bcnvp

      logical,optional :: extrp_b,dbg,newtime
      
c     Local variables

      integer :: i_bc,j,k,nx,ny,nz,nn,mm,my_rank_y,nxg,my_rank_yz

      real(8) :: bcovext  (0:size(bcnv,2)-1
     .                    ,0:size(bcnv,3)-1,3)
     .          ,bcovint  (0:size(bcnv,2)-1
     .                    ,0:size(bcnv,3)-1,3)
     .          ,bcnv1int (0:size(bcnv,2)-1
     .                    ,0:size(bcnv,3)-1)
     .          ,bcnv1intm(0:size(bcnv,2)-1
     .                    ,0:size(bcnv,3)-1)


      real(8) :: dum,dh,bcov_avg(2:3),send_buf(2),rec_buf(2),nfactor
      real(8) :: mag

      logical :: extrap,debug,new_time
      
c     Begin program

      nxg = g_def%nxgl(igr)
      if (g_def%ihi(igr) /= nxg) return
      
      nx = g_def%nxv(igr) 
      ny = g_def%nyv(igr)
      nz = g_def%nzv(igr)

      i_bc = nx + 1
      
      if (PRESENT(extrp_b)) then
        extrap = extrp_b
      else
        extrap = .false.
      endif

      if (PRESENT(dbg)) then
        debug = dbg
      else
        debug = .true.
      endif

      if (PRESENT(bcnvp)) then
        new_time = .true.
      elseif (PRESENT(newtime)) then
        new_time = newtime
      else
        new_time = rw_fft_tupdate
      endif

#if defined(petsc)
      call MPI_Comm_rank(g_def%MPI_COMM_YZ,my_rank_yz,mpierr)
      debug_io = debug.and.my_rank_yz==0.and.(.not.test)
#else
      my_rank_yz = 0
      debug_io = debug.and.(.not.test)
#endif

c     Compute INTERNAL B-field components @ RW

      bcnv1int  = bcnv1_plasma_bdry(bcnv,extrap)       !Previous time step

      if (PRESENT(bcnvp)) then                         !New time step
        bcovint   = bcov_plasma_bdry(bcnvp,.false.)
        bcnv1intm = bcnv1_plasma_int(bcnvp)
      else                                             !Old time step
        bcovint   = bcov_plasma_bdry(bcnv,extrap.or.rw_exp_fix)
        bcnv1intm = bcnv1_plasma_int(bcnv)
      endif

c     Compute EXTERNAL tangential B-field components @ RW (VACUUM response)

      select case(coords)
      case('hel')
       call vacuum_hel(g_def,igr,bcnv1int,bcnv1intm,bcovext
     $                ,new_time,debug)
c$$$      case('cyl')
      case('tor','sha','sh2')
        if (nzd == 1) then
          call vacuum_tor_symm(g_def,igr,bcnv1int,bcovext,dbg)
        else
          call critical('rw-bc'
     $      ,'VACUUM B-field response not implemented for non-symm tor')
        endif
      case default
        call critical('rw_bc'
     .      ,'VACUUM B-field response not implemented for '//coords)
      end select

c     Find electric field components @ r=1 (store in (nx+1))

      !Cnv from curl
      E_0(i_bc,:,:,1) = 0d0
      E_0(i_bc,:,:,2) = -(bcovext(:,:,3)-bcovint(:,:,3))/tau_w
      E_0(i_bc,:,:,3) =  (bcovext(:,:,2)-bcovint(:,:,2))/tau_w
      
      !Xform E_0 to Cov @ boundary
      do k=0,nz+1
        do j=0,ny+1
          E_0(i_bc,j,k,:)
     .      =XformToCov_ijk_face(g_def,i_bc-1,j,k,igr,E_0(i_bc,j,k,:),1)
        enddo
      enddo

c     Debug IO

      if (debug.and.(.not.test)) then
        mag = maxval(abs(E_0(i_bc,1:ny,1:nz,:)))
#if defined(petsc)
        mag = pmax(mag,mpi_comm=g_def%MPI_COMM_YZ)
#endif
        if (my_rank_yz==0) write (*,*) "max RW Ebc=",mag
      elseif (my_rank_yz==0) then
        write (*,*) "Computing RW BCs..."
      endif
      
c     End program

      contains

c     bcnv1_plasma_bdry
c     ####################################################################      
      function bcnv1_plasma_bdry(bcnv,extrap) result(bcnv1int)

c     --------------------------------------------------------------------
c     Compute Bcnv^1 @ wall on plasma side
c     --------------------------------------------------------------------

        implicit none

c     Call variables
        
        real(8),dimension(0:,0:,0:,:) :: bcnv
        real(8) :: bcnv1int(0:size(bcnv,2)-1
     .                     ,0:size(bcnv,3)-1)

        logical :: extrap
                
c     Local variables

        integer :: ig,jg,kg,i,j,k,ord
        real(8) :: xh(1)

c     Begin program

        if (extrap) then
          !Extrapolate to face from interior
          i = i_bc-1 ; ord = 2
          call getMGmap(g_def,i,1,1,igr,igr,igr,ig,jg,kg)
          xh = 0.5*(g_def%xx(ig)+g_def%xx(ig+1))
        else
          !Interpolate to face
          i = i_bc ; ord = 1
          call getMGmap(g_def,i,1,1,igr,igr,igr,ig,jg,kg)
          xh = 0.5*(g_def%xx(ig)+g_def%xx(ig-1))
        endif
        
        do k=0,nz+1
          do j=0,ny+1
            call IntDriver1d(4,g_def%xx(ig-3:ig),bcnv(i-3:i,j,k,1)
     .                        ,1,xh,bcnv1int(j,k),ord,0)
          enddo
        enddo

c$$$          bcnv1int = 0.5*bcnv(i_bc  ,:,:,1)
c$$$     .              +0.5*bcnv(i_bc-1,:,:,1)

      end function bcnv1_plasma_bdry

c     bcov_plasma_bdry
c     ####################################################################      
      function bcov_plasma_bdry(bcnv,extrap) result(bcovint)

c     --------------------------------------------------------------------
c     Compute Bcov @ wall on plasma side (for E-field computation)
c     --------------------------------------------------------------------

        implicit none

c     Call variables
        
        real(8),dimension(0:,0:,0:,:) :: bcnv
        real(8) :: bcovint (0:size(bcnv,2)-1
     .                     ,0:size(bcnv,3)-1,3)

        logical :: extrap
        
c     Local variables
        
        integer :: ig,jg,kg,i,j,k,ord,nyg,nzg
        real(8) :: bcov_avg(2:3),send_buf(2),rec_buf(2),bcovl(4,3),xh(1)

c     Begin program

        if (extrap) then
          i = i_bc-1 ; ord = 2
          call getMGmap(g_def,i,1,1,igr,igr,igr,ig,jg,kg)
          xh = 0.5*(g_def%xx(ig)+g_def%xx(ig+1))
        else
          !WARNING: tangential comps. seed RW instability
          !         for small tau_w for arbitrary dt
          i = i_bc ; ord = 1
          call getMGmap(g_def,i,1,1,igr,igr,igr,ig,jg,kg)
          xh = 0.5*(g_def%xx(ig)+g_def%xx(ig-1))
        endif
        
        do k=0,nz+1
          do j=0,ny+1
            bcovl(1,:) = XformToCov(g_def,i-3,j,k,igr,bcnv(i-3,j,k,:))
            bcovl(2,:) = XformToCov(g_def,i-2,j,k,igr,bcnv(i-2,j,k,:))
            bcovl(3,:) = XformToCov(g_def,i-1,j,k,igr,bcnv(i-1,j,k,:))
            bcovl(4,:) = XformToCov(g_def,i  ,j,k,igr,bcnv(i  ,j,k,:))

            call IntDriver1d(4,g_def%xx(ig-3:ig),bcovl(:,1)
     .                        ,1,xh,bcovint(j,k,1),ord,0)
            call IntDriver1d(4,g_def%xx(ig-3:ig),bcovl(:,2)
     .                        ,1,xh,bcovint(j,k,2),ord,0)
            call IntDriver1d(4,g_def%xx(ig-3:ig),bcovl(:,3)
     .                        ,1,xh,bcovint(j,k,3),ord,0)
          enddo
        enddo

c$$$          do k=0,nz+1
c$$$            do j=0,ny+1
c$$$              bcovint(j,k,:)
c$$$     .          =0.5*XformToCov(g_def,i_bc-1,j,k,igr,bcnv(i_bc-1,j,k,:))
c$$$     $          +0.5*XformToCov(g_def,i_bc  ,j,k,igr,bcnv(i_bc,  j,k,:))
c$$$            enddo
c$$$          enddo
        
c     Perform angular average: EXCLUDE m=n=0 mode (CRITICAL)

        bcov_avg(2) = sum(bcovint(1:ny,1:nz,2))
        bcov_avg(3) = sum(bcovint(1:ny,1:nz,3))

#if defined(petsc)
        send_buf = bcov_avg
        call MPI_Allreduce(send_buf,rec_buf,2,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,g_def%MPI_COMM_YZ,mpierr)
        bcov_avg = rec_buf
#endif

        nyg = g_def%nygl(igr)
        nzg = g_def%nzgl(igr)

        bcov_avg = bcov_avg/(nyg*nzg)
      
        !Subtract average
        bcovint(:,:,2) = bcovint(:,:,2) - bcov_avg(2)
        bcovint(:,:,3) = bcovint(:,:,3) - bcov_avg(3)

      end function bcov_plasma_bdry
      
c     bcnv1_plasma_int
c     ####################################################################      
      function bcnv1_plasma_int(bcnv) result(bcnv1intm)

c     --------------------------------------------------------------------
c     Store Bcnv inside plasma (dr from wall) for radial-B response
c     --------------------------------------------------------------------

        implicit none
      
        real(8),dimension(0:,0:,0:,:) :: bcnv
        real(8) :: bcnv1intm(0:size(bcnv,2)-1
     .                      ,0:size(bcnv,3)-1)

        dr_w = (g_def%xg(nxg)-g_def%xg(nxg-1))

c$$$        bcnvintm(:,:,:) = bcnv(i_bc-2,:,:,:)
c$$$        dr_w = 1.5*dr_w
        bcnv1intm(:,:) = bcnv(i_bc-1,:,:,1)
        dr_w = 0.5*dr_w

      end function bcnv1_plasma_int

      end subroutine rw_Ebc

c     vacuum_hel
c     ####################################################################
      subroutine vacuum_hel(g_def,igr,bcnvn_int,bcnvn_int_m,bcov_ext
     $                     ,new_time,dbg)
      
c     -----------------------------------------------------------------
c     Performs RW BC. Inputs normal internal B-component, outputs
c     tangential external B-components
c     -----------------------------------------------------------------

      use special_funcs
      
      implicit none

c     Call variables
        
      type(grid_mg_def),pointer :: g_def

      integer :: igr
      
      real(8) :: bcnvn_int(0:,0:),bcnvn_int_m(0:,0:),bcov_ext(0:,0:,:)

      logical :: dbg,new_time
      
c     Local variables
      
      integer :: m,n,jglobal,i,j,m1,mm2,mm1,nx,ny,nz,nyg
      real(8),dimension(:),pointer :: theta,mm,phir_m,phii_m
     .                               ,Brn_m,Bin_m,Br_m,Bi_m
     .                               ,gBr_m,gBi_m,Br,Bth,Bz
     .                               ,phir_m_np,phii_m_np

      real(8) :: dImv,dIm,Im,d2Im,dKmv,dKm,Km,d2Km,pv,pw,RR,k_h,rm
     .          ,efac,gamma,dum,m_h,dr
     .          ,bcov_ext_np(0:size(bcov_ext,1)-1
     .                      ,0:size(bcov_ext,2)-1
     .                      ,  size(bcov_ext,3)  )

      real(8),allocatable,dimension(:,:,:) :: Bn,Bnl,Bnml,Bnm

c     Begin program

      !Helical pitch modes
      m_h = g_def%params(1)
      k_h = g_def%params(2)

      nx = g_def%nxv(igr) 
      ny = g_def%nyv(igr)
      nz = g_def%nzv(igr)

      nyg = g_def%nygl(igr)

      allocate(Bnml(1,ny,1),Bnm(1,nyg,1))
      allocate( Bnl(1,ny,1), Bn(1,nyg,1))

      allocate(theta(nyg),mm(nyg),phir_m(nyg),phii_m(nyg)
     .        ,Br(nyg),Bth(nyg),Bz(nyg))

      allocate(Brn_m(nyg),Bin_m(nyg),Br_m(nyg),Bi_m(nyg))

      allocate(gBr_m(nyg),gBi_m(nyg))

      allocate(phir_m_np(nyg),phii_m_np(nyg))
      
c     Global angular mesh

      theta = g_def%yg(1:nyg)

c     FFT Bnormal at plasma wall (divided by J|n|=1/m* @ r=1)

      Bnl(1,:,1) = m_h*bcnvn_int(1:ny,1)  

      !Gather Bnormal in parallel along boundary
#if defined(petsc)
      call find_global_nobc(Bnl,Bn,mpi_comm=g_def%MPI_COMM_Y)
#else
      call find_global_nobc(Bnl,Bn)
#endif

      !Perform FFT
      call fft(nyg,Bn(1,:,1),theta,nyg,Brn_m,Bin_m,mm,0,0)

c$$$      if (m_max > 0) call filter_br()
      
c     Process plasma B-normal

      if (new_time) then
        !Transform normal to cov by dividing by J=r/m*, r=1
        Bnml(1,:,1) = m_h*bcnvn_int_m(1:ny,1)

        !Gather inner plasma B-field in parallel along boundary
#if defined(petsc)
        call find_global_nobc(Bnml,Bnm,mpi_comm=g_def%MPI_COMM_Y)
#else
        call find_global_nobc(Bnml,Bnm)
#endif
        !Perform FFT
        call fft(nyg,Bnm(1,:,1),theta,nyg,gBr_m,gBi_m,mm,0,0)
      endif
      
c     Find VACUUM potential solution Phi outside of resistive wall

      if (k_h == 0d0) then
        do m=1,nyg
          m1   = abs(nint(mm(m)))
          if (m1 == 0) then
            phir_m   (m) = 0d0 ; phii_m   (m) = 0d0
            phir_m_np(m) = 0d0 ; phii_m_np(m) = 0d0
            cycle
          endif

          rm = r_v**(2*mm(m))

          !Evolve normal B component
          if (new_time) then
            gamma = 1 + mm(m)*(rm+1d0)/(rm-1d0)
            efac  = 1d0/(1d0 + (gamma + 1d0/dr_w)*dt/tau_w)
            Br_m(m) = efac*(Brn_m(m)+dt/(tau_w*dr_w)*gBr_m(m))
            Bi_m(m) = efac*(Bin_m(m)+dt/(tau_w*dr_w)*gBi_m(m))
          endif
          
          !Find phi_nm
          dum = (rm+1d0)/((1d0-rm)*mm(m))
          if (new_time) then
            phir_m_np(m) = Br_m(m)*dum
            phii_m_np(m) = Bi_m(m)*dum
          endif
          phir_m   (m) = Brn_m(m)*dum
          phii_m   (m) = Bin_m(m)*dum
        enddo
      else
        do m=1,nyg
          m1   = abs(nint(mm(m)))
          if (m1 == 0) then
            phir_m   (m) = 0d0 ; phii_m   (m) = 0d0
            phir_m_np(m) = 0d0 ; phii_m_np(m) = 0d0
            cycle
          endif
          
          mm1  = abs(m1-1)
          mm2  = abs(m1-2)

          pw = m1*abs(k_h)/m_h  !Forces helical pitch, k_n = m x (k_n*/m*)
          pv = pw*r_v

          !Bessel function evals
c$$$          Km   = BESSK(m1 ,pw)                   
c$$$          dKm  =-BESSK(mm1,pw) - m1/pw*BESSK(m1,pw)
c$$$          d2Km =                   BESSK(mm2,pw)
c$$$     .         +  (2*m1-1)/ pw    *BESSK(mm1,pw)
c$$$     .         + m1*(m1+1)/(pw*pw)*BESSK(m1 ,pw)
c$$$          dKmv =-BESSK(mm1,pv) - m1/pv*BESSK(m1,pv)
c$$$          Im   = BESSI(m1 ,pw)
c$$$          dIm  = BESSI(mm1,pw) - m1/pw*BESSI(m1,pw)
c$$$          d2Im =                   BESSI(mm2,pw)
c$$$     .         -  (2*m1-1)/ pw    *BESSI(mm1,pw)
c$$$     .         + m1*(m1+1)/(pw*pw)*BESSI(m1 ,pw)
c$$$          dImv = BESSI(mm1,pv) - m1/pv*BESSI(m1,pv)
          Km   = BESSELK(1d0*m1 ,pw)                   
          dKm  =-BESSELK(1d0*mm1,pw) - m1/pw*BESSELK(1d0*m1,pw)
          d2Km =                   BESSELK(1d0*mm2,pw)
     .         +  (2*m1-1)/ pw    *BESSELK(1d0*mm1,pw)
     .         + m1*(m1+1)/(pw*pw)*BESSELK(1d0*m1 ,pw)
          dKmv =-BESSELK(1d0*mm1,pv) - m1/pv*BESSELK(1d0*m1,pv)
          Im   = BESSELI(1d0*m1 ,pw)
          dIm  = BESSELI(1d0*mm1,pw) - m1/pw*BESSELI(1d0*m1,pw)
          d2Im =                   BESSELI(1d0*mm2,pw)
     .         -  (2*m1-1)/ pw    *BESSELI(1d0*mm1,pw)
     .         + m1*(m1+1)/(pw*pw)*BESSELI(1d0*m1 ,pw)
          dImv = BESSELI(1d0*mm1,pv) - m1/pv*BESSELI(1d0*m1,pv)
          dum = (dImv*Km-dKmv*Im)/((dImv*dKm-dKmv*dIm)*pw)

          !Find phi_nm
          phir_m(m) = Brn_m(m)*dum
          phii_m(m) = Bin_m(m)*dum

          if (new_time) then
            !Evolve normal B component
            gamma =-pw*(dImv*d2Km-dKmv*d2Im)/(dImv*dKm -dKmv*dIm)
            efac  = 1d0/(1d0 + (gamma + 1d0/dr_w)*dt/tau_w)
            Br_m(m) = efac*(Brn_m(m)+dt/(tau_w*dr_w)*gBr_m(m))
            Bi_m(m) = efac*(Bin_m(m)+dt/(tau_w*dr_w)*gBi_m(m))

            phir_m_np(m) = Br_m(m)*dum
            phii_m_np(m) = Bi_m(m)*dum
          endif
        enddo

c$$$        k_h = -2d0 ; m_h = 1d0 ; r_v = 1.5d0
c$$$        
c$$$        do m1=1,4
c$$$          mm1  = abs(m1-1)
c$$$          mm2  = abs(m1-2)
c$$$
c$$$          pw = m1*abs(k_h)/m_h  !Forces helical pitch, k_n = m x (k_n*/m*)
c$$$          pv = pw*r_v
c$$$
c$$$          !Bessel function evals
c$$$          Km   = BESSK(m1 ,pw)                   
c$$$          dKm  =-BESSK(mm1,pw) - 1d0*m1/pw*BESSK(m1,pw)
c$$$          d2Km =                   BESSK(mm2,pw)
c$$$     .         +  (2d0*m1-1d0)/ pw    *BESSK(mm1,pw)
c$$$     .         + 1d0*m1*(1d0*m1+1d0)/(pw*pw)*BESSK(m1 ,pw)
c$$$          dKmv =-BESSK(mm1,pv) - 1d0*m1/pv*BESSK(m1,pv)
c$$$          Im   = BESSI(m1 ,pw)
c$$$          dIm  = BESSI(mm1,pw) - m1/pw*BESSI(m1,pw)
c$$$          d2Im =                   BESSI(mm2,pw)
c$$$     .         -  (2d0*m1-1d0)/ pw    *BESSI(mm1,pw)
c$$$     .         + 1d0*m1*(1d0*m1+1d0)/(pw*pw)*BESSI(m1 ,pw)
c$$$          dImv = BESSI(mm1,pv) - 1d0*m1/pv*BESSI(m1,pv)
c$$$          Km   = BESSELK(1d0*m1 ,pw)                   
c$$$          dKm  =-BESSELK(1d0*mm1,pw) - m1/pw*BESSELK(1d0*m1,pw)
c$$$          d2Km =                   BESSELK(1d0*mm2,pw)
c$$$     .         +  (2*m1-1)/ pw    *BESSELK(1d0*mm1,pw)
c$$$     .         + m1*(m1+1)/(pw*pw)*BESSELK(1d0*m1 ,pw)
c$$$          dKmv =-BESSELK(1d0*mm1,pv) - m1/pv*BESSELK(1d0*m1,pv)
c$$$          Im   = BESSELI(1d0*m1 ,pw)
c$$$          dIm  = BESSELI(1d0*mm1,pw) - m1/pw*BESSELI(1d0*m1,pw)
c$$$          d2Im =                   BESSELI(1d0*mm2,pw)
c$$$     .         -  (2*m1-1)/ pw    *BESSELI(1d0*mm1,pw)
c$$$     .         + m1*(m1+1)/(pw*pw)*BESSELI(1d0*m1 ,pw)
c$$$          dImv = BESSELI(1d0*mm1,pv) - m1/pv*BESSELI(1d0*m1,pv)
c$$$          dum = (dImv*Km-dKmv*Im)/((dImv*dKm-dKmv*dIm)*pw)
c$$$          dum = -1d0*m1/m_h*(dImv*Km-dKmv*Im)/((dImv*dKm-dKmv*dIm)*pw)
c$$$
c$$$          write (*,*) m1,dum
c$$$        enddo
c$$$        stop
      endif

c     Find old-time EXTERNAL magnetic field components outside RW

      !Bth [multiply by -sqrt(-1)*m/m*]
      Brn_m = mm/m_h*phii_m
      Bin_m =-mm/m_h*phir_m
      call ifft(nyg,Bth,theta,nyg,Brn_m,Bin_m,mm,0)

      !Get LOCAL helical (cov) magnetic field components @ r=1 
      do j = 0,ny+1
        jglobal = j + g_def%jlo(1) - 1
        if (jglobal > nyg) jglobal = jglobal - nyg
        if (jglobal < 1  ) jglobal = nyg + jglobal
        bcov_ext(j,:,1) = 0d0 !Not needed
        bcov_ext(j,:,2) = Bth(jglobal)
        bcov_ext(j,:,3) = 0d0 !implies Bz(jglobal) = k_h*bcov_ext(j,:,2)
      enddo

c     Find new-time EXTERNAL magnetic field components outside RW

      if (new_time) then
        !New-time Bth  [multiply by -sqrt(-1)*m/m*]
        Brn_m = mm/m_h*phii_m_np
        Bin_m =-mm/m_h*phir_m_np
        call ifft(nyg,Bth,theta,nyg,Brn_m,Bin_m,mm,0)

        !Get LOCAL helical (cov) magnetic field components @ r=1 
        do j = 0,ny+1
          jglobal = j + g_def%jlo(1) - 1
          if (jglobal > nyg) jglobal = jglobal - nyg
          if (jglobal < 1  ) jglobal = nyg + jglobal
          bcov_ext_np(j,:,1) = 0d0 !Not needed
          bcov_ext_np(j,:,2) = Bth(jglobal)
          bcov_ext_np(j,:,3) = 0d0 !implies Bz(jglobal) = k_h*bcov_ext_np(j,:,2)
        enddo
     
        if (dbg.and.save_tn) then
          Bnl(1,1:ny,1)=(bcov_ext_np(1:ny,1,2)/bcov_ext(1:ny,1,2)-1d0)
     $                 /dt
#if defined(petsc)
          dum = psum(Bnl(1,1:ny,1),mpi_comm=g_def%MPI_COMM_Y)/nyg
#else
          dum = sum(Bnl(1,1:ny,1))/nyg
#endif
          if (debug_io)
     $      write(*,*) "Bt norm. decay rate, m =",nint(m_h)," is",dum
        endif

        !Return new-time EXTERNAL magnetic field components outside RW
        bcov_ext = bcov_ext_np
      endif
       
      deallocate(Bnl,Bnml,Bn,Bnm)
      
      deallocate(theta,mm,phir_m,phii_m,Br,Bth,Bz)

      deallocate(Brn_m,Bin_m,Br_m,Bi_m)

      deallocate(gBr_m,gBi_m)

      deallocate(phir_m_np,phii_m_np)

c$$$      contains
c$$$
c$$$c     filter_br
c$$$c     ####################################################################
c$$$      subroutine filter_br()
c$$$
c$$$        do m=1,nyg
c$$$          m1   = abs(nint(mm(m)))
c$$$          if (m1 > m_max) then
c$$$            Brn_m(m) = 0d0 ; Bin_m(m) = 0d0
c$$$          endif
c$$$        enddo
c$$$
c$$$      end subroutine filter_br
      
      end subroutine vacuum_hel

c     vacuum_cyl
c     ####################################################################
c$$$      subroutine vacuum_cyl(g_def,igr,bcnvn_int,bcnvn_int_m,bcov_ext)
      subroutine vacuum_cyl(g_def,igr,bcnvn_int,bcov_ext)
      
c     -----------------------------------------------------------------
c     Performs RW BC. Inputs normal internal B-component, outputs
c     tangential external B-components
c     -----------------------------------------------------------------

      implicit none

c     Call variables
        
      type(grid_mg_def),pointer :: g_def

      integer :: igr
      
      real(8) :: bcnvn_int(0:,0:),bcov_ext(0:,0:,:)!,bcnvn_int_m(0:,0:)

c     Local variables
      
      integer:: m,n,jglobal,kglobal,i,j,k,m1,n1,mm2,mm1,nx,ny,nz,nyg,nzg
      real(8),dimension(:),pointer :: theta,zz,mm,nn
      real(8),dimension(:,:),pointer :: phir_m,phii_m
     .                                 ,Brn_m,Bin_m,Br_m,Bi_m
     .                                 ,gBr_m,gBi_m,Bth,Bz

      real(8) :: dImv,dIm,Im,d2Im,dKmv,dKm,Km,d2Km,pv,pw,RR,rm
     .          ,efac,gamma,dum
     .          ,bcov_ext_np(0:size(bcov_ext,1)-1
     .                      ,0:size(bcov_ext,2)-1
     .                      ,  size(bcov_ext,3)  )

      real(8),allocatable,dimension(:,:,:) :: Bn,Bnl!,Bnml,Bnm

c     Begin program

      !Cylindrical params
      RR = g_def%params(3)

      nx = g_def%nxv(igr) 
      ny = g_def%nyv(igr)
      nz = g_def%nzv(igr)

      nyg = g_def%nygl(igr)
      nzg = g_def%nzgl(igr)

      allocate(Bnl(1,ny,nz),Bn(1,nyg,nzg))

      allocate(theta(nyg),zz(nzg),mm(nyg),nn(nzg)
     $        ,phir_m(nyg,nzg),phii_m(nyg,nzg)
     .        ,Bth(nyg,nzg),Bz(nyg,nzg))

      allocate(Brn_m(nyg,nzg),Bin_m(nyg,nzg))

c$$$      allocate(Bnml(1,ny,nz),Bnm(1,nyg,nzg))
c$$$      allocate(Br_m(nyg,nzg),Bi_m(nyg,nzg))
c$$$      allocate(gBr_m(nyg,nzg),gBi_m(nyg,nzg))
c$$$      allocate(phir_m_np(nyg,nzg),phii_m_np(nyg,nzg))
      
c     Global angular mesh

      theta = g_def%yg(1:nyg)
      zz    = g_def%zg(1:nzg)

c     FFT Bnormal at plasma wall (transform to cov by dividing by J=r, r=1)

      Bnl(1,:,:) = bcnvn_int(1:ny,1:nz)  

      !Gather Bnormal in parallel along boundary
#if defined(petsc)
      call find_global_nobc(Bnl,Bn,mpi_comm=g_def%MPI_COMM_YZ)
#else
      call find_global_nobc(Bnl,Bn)
#endif

      !Perform 2D FFT
cc>>>TODO      call fft2d(nyg,Bn(1,:,:),theta,nyg,Brn_m,Bin_m,mm,0,0)
      
c$$$c     Process gradient of Bnormal (new time)
c$$$      
c$$$      !Compute gradient of Bn (transform to cov by dividing by J=r, r=1)
c$$$      Bnml(1,:,:) = bcnvn_int_m(1:ny,1:nz)
c$$$
c$$$      !Gather gradient in parallel along boundary
c$$$#if defined(petsc)
c$$$      call find_global_nobc(Bnml,Bnm,mpi_comm=g_def%MPI_COMM_YZ)
c$$$#else
c$$$      call find_global_nobc(Bnml,Bnm)
c$$$#endif
c$$$
c$$$      !Perform 2D FFT
c$$$cc>>>TODO      call fft2d(nyg,Bnm(1,:,1),theta,nyg,gBr_m,gBi_m,mm,0,0)

c     Find VACUUM potential solution Phi outside of resistive wall

      do n=1,nzg
        do m=1,nyg
          n1 = abs(nint(nn(n)))
          m1 = abs(nint(mm(m)))

          if (m1 == 0.and.n1 == 0) then
            phir_m(m,n) = 0d0 ; phii_m(m,n) = 0d0
            cycle
          elseif (n1 == 0) then
            rm = r_v**(2*mm(m))

c$$$            !Evolve normal B component
c$$$            gamma = 1 + mm(m)*(rm+1d0)/(rm-1d0)
c$$$            efac  = exp(-gamma*dt/tau_w)
c$$$            Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
c$$$            Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)

            !Find phi_nm
            dum = (rm+1d0)/((1d0-rm)*mm(m))
c$$$            phir_m_np(m) = Br_m(m)*dum
c$$$            phii_m_np(m) = Bi_m(m)*dum
            phir_m   (m,n) = Brn_m(m,n)*dum
            phii_m   (m,n) = Bin_m(m,n)*dum
            cycle
          endif
          
          mm1  = abs(m1-1)
          mm2  = abs(m1-2)

          pw = 1d0*n1/RR
          pv = pw*r_v

          !Bessel function evals
          Km   = BESSK(m1 ,pw)                   
          dKm  =-BESSK(mm1,pw) - m1/pw*BESSK(m1,pw)
          d2Km =                   BESSK(mm2,pw)
     .         +  (2*m1-1)/ pw    *BESSK(mm1,pw)
     .         + m1*(m1+1)/(pw*pw)*BESSK(m1 ,pw)
          dKmv =-BESSK(mm1,pv) - m1/pv*BESSK(m1,pv)
          Im   = BESSI(m1 ,pw)
          dIm  = BESSI(mm1,pw) - m1/pw*BESSI(m1,pw)
          d2Im =                   BESSI(mm2,pw)
     .         -  (2*m1-1)/ pw    *BESSI(mm1,pw)
     .         + m1*(m1+1)/(pw*pw)*BESSI(m1 ,pw)
          dImv = BESSI(mm1,pv) - m1/pv*BESSI(m1,pv)

c$$$          !Evolve normal B component
c$$$          gamma =-pw*(dImv*d2Km-dKmv*d2Im)/(dImv*dKm -dKmv*dIm)
c$$$          efac  = exp(-gamma*dt/tau_w)
c$$$          Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
c$$$          Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)

          !Find phi_nm
          dum = (dImv*Km-dKmv*Im)/((dImv*dKm-dKmv*dIm)*pw)
c$$$          phir_m_np(m) = Br_m(m)*dum
c$$$          phii_m_np(m) = Bi_m(m)*dum
          phir_m   (m,n) = Brn_m(m,n)*dum
          phii_m   (m,n) = Bin_m(m,n)*dum
        enddo
      enddo

c     Find old-time EXTERNAL magnetic field components outside RW

      !Bth [multiply by -sqrt(-1)*m]
      do n=1,nzg
        Brn_m(:,n) = mm(:)*phii_m(:,n)
        Bin_m(:,n) =-mm(:)*phir_m(:,n)
      enddo
cc>>>TODO      call ifft2D(nyg,Bth,theta,nyg,Brn_m,Bin_m,mm,0)

      !Bz [multiply by -sqrt(-1)*k]
      do m=1,nyg
        Brn_m(m,:) = nn(:)*phii_m(m,:)/RR
        Bin_m(m,:) =-nn(:)*phir_m(m,:)/RR
      enddo
cc>>>TODO      call ifft2D(nyd,Bz,theta,nyd,Br_m,Bi_m,mm,0)

      !Get LOCAL cov magnetic field components @ r=1 
      do k = 0,ny+1
        do j = 0,ny+1
          jglobal = j + g_def%jlo(1) - 1
          if (jglobal > nyg) jglobal = jglobal - nyg
          if (jglobal < 1  ) jglobal = nyg + jglobal
          kglobal = k + g_def%klo(1) - 1
          if (kglobal > nzg) kglobal = kglobal - nzg
          if (kglobal < 1  ) kglobal = nzg + kglobal
          bcov_ext(j,k,1) = 0d0 !Not needed
          bcov_ext(j,k,2) = Bth(jglobal,kglobal)
          bcov_ext(j,k,3) = Bz (jglobal,kglobal)
        enddo
      enddo
      
c$$$c     Find new-time EXTERNAL magnetic field components outside RW
c$$$
c$$$      if (new_time) then
c$$$        !New-time normal component Br
c$$$        call ifft(nyg,Br,theta,nyg,Br_m,Bi_m,mm,0)
c$$$
c$$$        !New-time Bth  [multiply by -sqrt(-1)*m/m*]
c$$$        Brn_m = mm/m_h*phii_m_np
c$$$        Bin_m =-mm/m_h*phir_m_np
c$$$        call ifft(nyg,Bth,theta,nyg,Brn_m,Bin_m,mm,0)
c$$$
c$$$        !Get LOCAL helical (cov) magnetic field components @ r=1 
c$$$        do j = 0,ny+1
c$$$          jglobal = j + g_def%jlo(1) - 1
c$$$          if (jglobal > nyg) jglobal = jglobal - nyg
c$$$          if (jglobal < 1  ) jglobal = nyg + jglobal
c$$$          bcov_ext_np(j,:,1) = 0d0 !Not needed
c$$$          bcov_ext_np(j,:,2) = Bth(jglobal)
c$$$          bcov_ext_np(j,:,3) = 0d0 !implies Bz(jglobal) = k_h*bcov_ext_np(j,:,2)
c$$$        enddo
c$$$
c$$$        !Return new-time EXTERNAL magnetic field components outside RW
c$$$        bcov_ext = bcov_ext_np
c$$$      endif
       
      deallocate(Bnl,Bn)
      
      deallocate(theta,zz,mm,nn,phir_m,phii_m,Bth,Bz)

      deallocate(Brn_m,Bin_m)

c$$$      deallocate(Bnml,Bnm)
c$$$      deallocate(Br_m,Bi_m)
c$$$      deallocate(gBr_m,gBi_m)
c$$$      deallocate(phir_m_np,phii_m_np)
      
      end subroutine vacuum_cyl

c     vacuum_tor_symm
c     ####################################################################
      subroutine vacuum_tor_symm(g_def,igr,bcnvn_int,bcov_ext,dbg)
      
c     -----------------------------------------------------------------
c     Performs RW BC in toroidal geometry. Inputs normal internal
c     B-component, outputs tangential external B-components
c     -----------------------------------------------------------------

      use rw_bim_mod, only: rw_bim_symm_solve

      implicit none

c     Call variables
        
      type(grid_mg_def),pointer :: g_def

      integer :: igr
      
      real(8) :: bcnvn_int(0:,0:),bcov_ext(0:,0:,:)

      logical :: dbg

      real(8),allocatable,dimension(:,:,:) :: Bn,Bng
      real(8),allocatable,dimension(:)     :: phi
      real(8),allocatable,dimension(:)     :: theta

c     Local variables
      
      integer :: jgl,jglp,jglm,i,j,k,nx,ny,nz,nyg

c     Begin program

      nx = g_def%nxv(igr) 
      ny = g_def%nyv(igr)
      nz = g_def%nzv(igr)

      nyg = g_def%nygl(igr)

      allocate(Bn(1,ny,1),phi(nyg),Bng(1,nyg,1),theta(nyg))
      
      theta = g_def%yg(1:nyg)

c     Bnormal at plasma wall (pointing inward)

      Bn(1,:,1) =-bcnvn_int(1:ny,1)

#if defined(petsc)
      call find_global_nobc(Bn,Bng,mpi_comm=g_def%MPI_COMM_Y)
#else
      call find_global_nobc(Bn,Bng)
#endif
      
c     Find VACUUM potential solution Phi outside of resistive wall

      phi = 0d0

      call rw_bim_symm_solve(Bng(1,:,1),phi)

c     Find EXTERNAL magnetic field components outside RW

      !Get LOCAL (cov) magnetic field components @ r=1 
      do j = 0,ny+1
        jgl = j + g_def%jlo(1) - 1
        jglp = jgl + 1
        jglm = jgl - 1
        if (jglp > nyg) jglp = jglp - nyg
        if (jglm < 1  ) jglm = nyg + jglm
        bcov_ext(j,:,1) = 0d0 !Not needed
        bcov_ext(j,:,2) = (phi  (jglp)-phi  (jglm))
     $                   /(theta(jglp)-theta(jglm))
        bcov_ext(j,:,3) = 0d0 !From toroidal symmetry
      enddo

      deallocate(Bn,Bng,phi,theta)

      end subroutine vacuum_tor_symm

c     rw_A_BCs
c     ##################################################################
      subroutine rw_A_BCs(g_def,acov_n,bcnv_n,acov,acnv,bcnd_a,bcnd_b
     .                   ,Erw,order)
c     ------------------------------------------------------------------
c     Defines RW BCs for vector potential
c     ------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      real(8),dimension(0:,0:,0:,:) :: acov_n,acov,bcnv_n,Erw,acnv

      integer :: bcnd_a(6,3),bcnd_b(6,3)

      integer,optional :: order
      
c     Local variables

      integer :: nnx,nny,nnz,igr,it
      real(8) :: xi_rw,mag
      real(8),allocatable,dimension(:,:,:,:) :: v0,vdum,bcnv

c     Begin program

      igr = 1
      
      nnx = g_def%nxv(igr) 
      nny = g_def%nyv(igr)
      nnz = g_def%nzv(igr)

      allocate(v0  (0:nnx+1,0:nny+1,0:nnz+1,3)
     .        ,vdum(0:nnx+1,0:nny+1,0:nnz+1,3)
     .        ,bcnv(0:nnx+1,0:nny+1,0:nnz+1,3))

      v0   = 0d0
      vdum = 0d0

      if (rw_implicit_bc) then
        call rw_A_imp
      elseif (rw_exp_pc_bc) then
        call rw_A_exp_pc
      else !default
        call rw_A_exp_ssp
      endif
      
      deallocate(v0,vdum,bcnv)

      contains

c     rw_A_imp
c     ##################################################################
      subroutine rw_A_imp

      implicit none

      if (save_tn) return
      
      !Second-order BDF2 implicit update of RW BC
      if (itime > inewtime + 2) then
        vdum(:,:,:,1) =-cn /cnp*u_n %array_var(IAX)%array
     $                 -cnm/cnp*u_nm%array_var(IAX)%array
        vdum(:,:,:,2) =-cn /cnp*u_n %array_var(IAY)%array
     $                 -cnm/cnp*u_nm%array_var(IAY)%array
        vdum(:,:,:,3) =-cn /cnp*u_n %array_var(IAZ)%array
     $                 -cnm/cnp*u_nm%array_var(IAZ)%array
      else
        vdum = acov_n
      endif

      call find_B(gv%gparams,igr,B0,acov,bcnv,.false.
     $           ,bcond_a=bcnd_a,bcond_b=bcnd_b)

      !Find RW E-field
      call rw_Ebc(gv%gparams,igr,bcnv_n,Erw,bcnvp=bcnv
     $           ,dbg=(ilevel>0))

      !A BCs
      v0 = 0d0
      call XferBCs(vdum,v0,.true.)
      if (itime > inewtime + 2) then
        v0(nnx+1,:,:,:) = v0(nnx+1,:,:,:)-dt/cnp*Erw(nnx+1,:,:,:)
      else
        v0(nnx+1,:,:,:) = v0(nnx+1,:,:,:)-dt    *Erw(nnx+1,:,:,:)
      endif

      call setBC(gv%gparams,IAX,3,nnx,nny,nnz,acnv,acov,v0
     $          ,bcnd_a,igr,igr,igr,iorder=order,is_cnv=.false.)
              
c$$$              if (.not.save_tn) then
c$$$                !Predict A @ BC by extrapolation and find B, Erw
c$$$                call find_B(gv%gparams,igr,B0,acov,bcnv,.false.
c$$$     $                ,bcond_a=bcnd_a
c$$$     $                ,bcond_b=bcnd_b)
c$$$
c$$$                !Correct A @ BC and find new Erw and iterate
c$$$                itmax = 2
c$$$                do it=1,itmax
c$$$                  v0 = 0d0
c$$$c$$$                  call rw_Ebc(gv%gparams,igr,bcnv,v0,extrp_b=.true.
c$$$                  call rw_Ebc(gv%gparams,igr,bcnv_n,v0,bcnvp=bcnv
c$$$     $                       ,rw_dbg=(it==itmax.and.ilevel>0))
c$$$
c$$$c$$$                  if (it == 1) then
c$$$                    Erw(nnx+1,:,:,:) = v0(nnx+1,:,:,:)
c$$$c$$$                  else
c$$$c$$$                    Erw(nnx+1,:,:,:) = 0.9*v0  (nnx+1,:,:,:)
c$$$c$$$     $                                + 0.1*Erw(nnx+1,:,:,:)
c$$$c$$$                  endif
c$$$                  
c$$$c$$$                  mag = rms(Erw(nnx+1,1:nny,1,3))
c$$$c$$$                  if (my_rank == 0) write (*,*) "Ebc rms it=",it,mag
c$$$                  
c$$$                  !A BCs
c$$$                  v0 = 0d0
c$$$                  call XferBCs(vdum,v0,.true.)
c$$$                  if (itime > inewtime + 2) then
c$$$                    v0(nnx+1,:,:,:) = v0(nnx+1,:,:,:)
c$$$     $                               -dt/cnp*Erw(nnx+1,:,:,:)
c$$$                  else
c$$$                    v0(nnx+1,:,:,:) = v0(nnx+1,:,:,:)
c$$$     $                               -dt*Erw(nnx+1,:,:,:)
c$$$                  endif
c$$$
c$$$                  call setBC(gv%gparams,IAX,3,nnx,nny,nnz,acnv,acov,v0
c$$$     $                    ,bcnd_a,igr,igr,igr,iorder=order,is_cnv=.false.)
c$$$
c$$$                  if (itmax == 1) exit
c$$$                  
c$$$                  !Correct B
c$$$                  call find_B(gv%gparams,igr,B0,acov,bcnv,.false.)
c$$$
c$$$                  v0 = 0d0 ; prs1 = rho*tmp
c$$$                  call J_pinch_bc(gv%gparams,igr,bcnv,prs1,eeta,Erw
c$$$     $                           ,v0)
c$$$                  bcov = XformVector(gv%gparams,igr,bcnv,.false.)
c$$$
c$$$                  call setBC(gv%gparams,IBX,3,nnx,nny,nnz,bcnv,bcov,v0
c$$$     $                ,bcnd_b,igr,igr,igr,iorder=order,is_cnv=.true.)
c$$$                enddo
c$$$  endif
        
      end subroutine rw_A_imp
      
c     rw_A_exp_ssp
c     ##################################################################
      subroutine rw_A_exp_ssp

      implicit none

      !Compute RW E-field using time-n B
      if (save_tn) call rw_Ebc(gv%gparams,igr,bcnv_n,Erw)
      
      !Second-order SSP explicit update of RW BC
      if ((itime > inewtime + 2).and.(.not.rw_fft_tupdate)
     $    .and.rw_exp_ssp_bc) then
        xi_rw = dt/(dt_n+dt_nm)
        vdum(:,:,:,1) = u_n%array_var(IAX)%array
     $                 -xi_rw**2*(u_n  %array_var(IAX)%array
     $                           -u_nm2%array_var(IAX)%array)
        vdum(:,:,:,2) = u_n%array_var(IAY)%array
     $                 -xi_rw**2*(u_n  %array_var(IAY)%array
     $                           -u_nm2%array_var(IAY)%array)
        vdum(:,:,:,3) = u_n%array_var(IAZ)%array
     $                 -xi_rw**2*(u_n  %array_var(IAZ)%array
     $                           -u_nm2%array_var(IAZ)%array)
      else
        xi_rw = 0d0
        vdum = acov_n
      endif

      call XferBCs(vdum,v0,.true.)
      v0(nnx+1,:,:,:) = v0(nnx+1,:,:,:)
     $                 -dt*(1d0+xi_rw)*Erw(nnx+1,:,:,:)

      !Impose BC
      if (.not.init_time_eval) then
        call setBC(gv%gparams,IAX,3,nnx,nny,nnz,acnv,acov,v0
     .        ,bcnd_a,igr,igr,igr,iorder=order,is_cnv=.false.)
      endif
          
      end subroutine rw_A_exp_ssp
      
c     rw_A_exp_pc
c     ##################################################################
      subroutine rw_A_exp_pc

      implicit none

      !Predictor-corrector explicit update on RW BC
      if (save_tn) then
        bcnv = bcnv_n
        do it=1,2
          call rw_Ebc(gv%gparams,igr,bcnv,Erw
     $               ,dbg=(ilevel>0.and.it==2))
      
          call XferBCs(acov_n,v0,.true.)
          v0(nnx+1,:,:,:) = v0(nnx+1,:,:,:)-dt*Erw(nnx+1,:,:,:)

          call setBC(gv%gparams,IAX,3,nnx,nny,nnz,acnv,acov
     .      ,v0,bcnd_a,igr,igr,igr,iorder=order,is_cnv=.false.)

          if (it==1) then
            !Correct B
            call find_B(gv%gparams,igr,B0,acov_n,bcnv,.false.
     $                 ,bcond_b=bcnd_b)

          !WARNING: This makes RW BC unstable for arbitrary dt for small tau_w
c$$$                call find_B(gv%gparams,igr,B0,acov,bcnv,.false.)
c$$$                v0 = 0d0 ; prs1 = rho*tmp
c$$$                call J_pinch_bc(gv%gparams,igr,bcnv,prs1,eeta,Erw
c$$$     $                         ,v0)
c$$$                bcov = XformVector(gv%gparams,igr,bcnv,.false.)
c$$$
c$$$                call setBC(gv%gparams,IBX,3,nnx,nny,nnz,bcnv,bcov
c$$$     $              ,v0,bcnd_b,igr,igr,igr,iorder=order,is_cnv=.true.)
          endif
        enddo
      else
        call XferBCs(acov_n,v0,.true.)
        v0(nnx+1,:,:,:) = v0(nnx+1,:,:,:) - dt*Erw(nnx+1,:,:,:)

        call setBC(gv%gparams,IAX,3,nnx,nny,nnz,acnv,acov,v0
     .            ,bcnd_a,igr,igr,igr,iorder=order,is_cnv=.false.)
      endif
      
      end subroutine rw_A_exp_pc
      
      end subroutine rw_A_BCs

      end module rw_bc
