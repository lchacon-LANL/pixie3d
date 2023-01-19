c module rw_bc
c ####################################################################
      module rw_bc

        use local_BCS_variables
      
        !Resistive wall vars
        logical :: rw_bcs(6)      =.false.
     $            ,rw_fft_tupdate =.false.

        real(8) :: tau_w=1d30,r_v=1d0  !Ideal wall by default

        real(8),pointer,dimension(:,:,:,:) :: E_bc

        logical,private :: debug_io=.false.
     .                    ,rw_dbg  =.true.

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

        integer :: nx,ny,nz
      
c     Begin program

        nx = gv%gparams%nxv(1) 
        ny = gv%gparams%nyv(1)
        nz = gv%gparams%nzv(1)

        allocate(E_bc(0:nx+1,0:ny+1,0:nz+1,3))

        E_bc = 0d0

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
      subroutine rw_Ebc(g_def,igr,bcnv,E_0,extrp_b)

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

      logical,optional :: extrp_b
      
c     Local variables

      integer :: i_bc,j,k,nx,ny,nz,nn,mm,my_rank_y,ig,jg,kg,nxg
     $          ,my_rank_yz

      real(8) :: phinm   (0:size(bcnv,2)-1
     .                   ,0:size(bcnv,3)-1)
     .          ,bcovext (0:size(bcnv,2)-1
     .                   ,0:size(bcnv,3)-1,3)
     .          ,bcnvint (0:size(bcnv,2)-1
     .                   ,0:size(bcnv,3)-1,3)
     .          ,gbcnvint(0:size(bcnv,2)-1
     .                   ,0:size(bcnv,3)-1,3)
     .          ,bcovint (0:size(bcnv,2)-1
     .                   ,0:size(bcnv,3)-1,3)

      real(8) :: dum,dh,bcov_avg(2:3),send_buf(2),rec_buf(2),nfactor
      real(8) :: mag

      logical :: extrap
      
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
      
#if defined(petsc)
      call MPI_Comm_rank(g_def%MPI_COMM_YZ,my_rank_yz,mpierr)
#else
      my_rank_yz = 0
#endif
      
      debug_io = rw_dbg.and.my_rank_yz==0.and.(.not.test)

c     Compute INTERNAL tangential magnetic field components, EXCLUDING m=n=0 mode (CRITICAL)

      if (extrap) then
c$$$        bcnvint (:,:,:) = bcnv(i_bc-1,:,:,:)
        !Extrapolate to face from interior
        call getMGmap(g_def,i_bc,1,1,igr,igr,igr,ig,jg,kg)
        bcnvint (:,:,:) = bcnv(i_bc-2,:,:,:)
     $                  +(bcnv(i_bc-1,:,:,:)-bcnv(i_bc-2,:,:,:))
     $                   /(g_def%xx(ig-1)-g_def%xx(ig-2))
     $                    *(0.5*g_def%xx(ig  )
     $                     +0.5*g_def%xx(ig-1)
     $                     -    g_def%xx(ig-2))
      else
        !Interpolate to face
        bcnvint (:,:,:)= 0.5*bcnv(i_bc  ,:,:,:)
     .                  +0.5*bcnv(i_bc-1,:,:,:)
      endif

      do k=0,nz+1
        do j=0,ny+1
          bcovint(j,k,:)
     .       =XformToCov_ijk_face(g_def,i_bc-1,j,k,igr,bcnvint(j,k,:),1)
        enddo
      enddo

      gbcnvint(:,:,:)=(bcnv(i_bc  ,:,:,:)
     .                -bcnv(i_bc-1,:,:,:))
     .                /(g_def%xg(nxg+1)-g_def%xg(nxg))

      !Perform angular average
      bcov_avg(2) = sum(bcovint(1:ny,1:nz,2))
      bcov_avg(3) = sum(bcovint(1:ny,1:nz,3))

#if defined(petsc)
      send_buf = bcov_avg
      call MPI_Allreduce(send_buf,rec_buf,2,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,g_def%MPI_COMM_YZ,mpierr)
      bcov_avg = rec_buf
#endif

      bcov_avg = bcov_avg/(nyd*nzd)
      
      !Subtract average
      bcovint(:,:,2) = bcovint(:,:,2) - bcov_avg(2)
      bcovint(:,:,3) = bcovint(:,:,3) - bcov_avg(3)

c     Compute VACUUM response

      select case(coords)
      case('hel')
       call vacuum_hel(g_def,igr,bcnvint(:,:,1),gbcnvint(:,:,1),bcovext)
c$$$      case('cyl')
c$$$      case('tor')
      case default
        call critical('rw_bc'
     .       ,'VACUUM response not implemented for '//coords)
      end select

c     Find electric field components @ r=1 (store in (nx+1))

      !Cnv from curl
      E_0(i_bc,:,:,1) = 0d0
      E_0(i_bc,:,:,2) =-(bcovext(:,:,3)-bcovint(:,:,3))/tau_w
      E_0(i_bc,:,:,3) = (bcovext(:,:,2)-bcovint(:,:,2))/tau_w
      
      !Xform E_0 to Cov @ boundary
      do k=0,nz+1
        do j=0,ny+1
          E_0(i_bc,j,k,:)
     .      =XformToCov_ijk_face(g_def,i_bc-1,j,k,igr,E_0(i_bc,j,k,:),1)
        enddo
      enddo

c     Debug IO
      
#if defined(petsc)
      mag = pmax(maxval(abs(E_0(i_bc,1:ny,1:nz,:)))
     $          ,mpi_comm=g_def%MPI_COMM_YZ)
#else
      mag = maxval(abs(E_0(i_bc,1:ny,1:nz,:)))
#endif
      if (debug_io) write (*,*) "max RW Ebc=",mag

c     End program

      end subroutine rw_Ebc

c     vacuum_hel
c     ####################################################################
      subroutine vacuum_hel(g_def,igr,bcnvn_int,gbcnvn_int,bcov_ext)
      
c     -----------------------------------------------------------------
c     Performs RW BC. Inputs normal internal B-component, outputs
c     tangential external B-components
c     -----------------------------------------------------------------

      implicit none

c     Call variables
        
      type(grid_mg_def),pointer :: g_def

      integer :: igr
      
      real(8) :: bcnvn_int(0:,0:),gbcnvn_int(0:,0:),bcov_ext(0:,0:,:)

c     Local variables
      
      integer :: m,n,jglobal,i,j,m1,mm2,mm1,nx,ny,nz,nyg
      real(8),dimension(:),pointer :: theta,mm,phir_m,phii_m
     .                               ,Brn_m,Bin_m,Br_m,Bi_m
     .                               ,gBr_m,gBi_m,Br,Bth,Bz
     .                               ,phir_m_np,phii_m_np

      real(8) :: dImv,dIm,Im,d2Im,dKmv,dKm,Km,d2Km,pv,pw,RR,k_h,rm
     .          ,efac,gamma,dum,m_h
     .          ,bcov_ext_np(0:size(bcov_ext,1)-1
     .                      ,0:size(bcov_ext,2)-1
     .                      ,  size(bcov_ext,3)  )

      real(8),allocatable,dimension(:,:,:) :: Bn,Bnl,gBnl,gBn

c     Begin program

      !Helical pitch modes
      m_h = g_def%params(1)
      k_h = g_def%params(2)

      nx = g_def%nxv(igr) 
      ny = g_def%nyv(igr)
      nz = g_def%nzv(igr)

      nyg = g_def%nygl(igr)

      allocate(gBnl(1,ny,1),gBn(1,nyg,1))
      allocate( Bnl(1,ny,1), Bn(1,nyg,1))

      allocate(theta(nyg),mm(nyg),phir_m(nyg),phii_m(nyg)
     .        ,Br(nyg),Bth(nyg),Bz(nyg))

      allocate(Brn_m(nyg),Bin_m(nyg),Br_m(nyg),Bi_m(nyg))

      allocate(gBr_m(nyg),gBi_m(nyg))

      allocate(phir_m_np(nyg),phii_m_np(nyg))
      
c     Global angular mesh

      theta = g_def%yg(1:nyg)

c     FFT Bnormal at plasma wall (transform to cov by dividing by J=r/m*, r=1)

      Bnl(1,:,1) = m_h*bcnvn_int(1:ny,1)  

      !Gather Bnormal in parallel along boundary
#if defined(petsc)
      call find_global_nobc(Bnl,Bn,mpi_comm=g_def%MPI_COMM_Y)
#else
      call find_global_nobc(Bnl,Bn)
#endif

      !Perform FFT
      call fft(nyg,Bn(1,:,1),theta,nyg,Brn_m,Bin_m,mm,0,0)
      
c     Process gradient of Bnormal (new time)

      if (rw_fft_tupdate) then
        !Compute gradient of Bn (transform to cov by dividing by J=r/m*, r=1)
        gBnl(1,:,1) = m_h*gbcnvn_int(1:ny,1)

        !Gather gradient in parallel along boundary
#if defined(petsc)
        call find_global_nobc(gBnl,gBn,mpi_comm=g_def%MPI_COMM_Y)
#else
        call find_global_nobc(gBnl,gBn)
#endif
        !Perform FFT
        call fft(nyg,gBn(1,:,1),theta,nyg,gBr_m,gBi_m,mm,0,0)
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
          if (rw_fft_tupdate) then
            gamma = 1 + mm(m)*(rm+1d0)/(rm-1d0)
            efac  = exp(-MIN(gamma*dt/tau_w,175d0))
            Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
            Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)
          endif
          
          !Find phi_nm
          dum = (rm+1d0)/((1d0-rm)*mm(m))
          if (rw_fft_tupdate) then
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

          dum = (dImv*Km-dKmv*Im)/((dImv*dKm-dKmv*dIm)*pw)

          !Find phi_nm
          phir_m(m) = Brn_m(m)*dum
          phii_m(m) = Bin_m(m)*dum

          if (rw_fft_tupdate) then
            !Evolve normal B component
            gamma =-pw*(dImv*d2Km-dKmv*d2Im)/(dImv*dKm -dKmv*dIm)
            efac  = exp(-MIN(gamma*dt/tau_w,175d0))
            Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
            Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)
c$$$            efac  = 1d0/(1d0 + gamma*dt/tau_w)
c$$$            Br_m(m) = efac*(Brn_m(m)-dt/tau_w*gBr_m(m))
c$$$            Bi_m(m) = efac*(Bin_m(m)-dt/tau_w*gBi_m(m))

            phir_m_np(m) = Br_m(m)*dum
            phii_m_np(m) = Bi_m(m)*dum
          endif
        enddo
      endif

c     Find old-time EXTERNAL magnetic field components outside RW

c diag ****
c$$$#if defined(plotting)
c$$$        do i=-nyd/2,nyd/2
c$$$          do j=1,nyd
c$$$            if (nint(mm(j))-i == 0) exit
c$$$          enddo
c$$$          write (*,'(a,i0,a,2e14.7)')
c$$$     $    "Br FFT m=",i," ",time,sqrt(Brn_m(j)**2+Bin_m(j)**2)
c$$$        enddo
c$$$#endif
c diag ****

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

      if (rw_fft_tupdate) then
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
     
        if (rw_dbg.and.save_tn) then
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
       
      deallocate(Bnl,gBnl,Bn,gBn)
      
      deallocate(theta,mm,phir_m,phii_m,Br,Bth,Bz)

      deallocate(Brn_m,Bin_m,Br_m,Bi_m)

      deallocate(gBr_m,gBi_m)

      deallocate(phir_m_np,phii_m_np)
      
      end subroutine vacuum_hel

c     vacuum_cyl
c     ####################################################################
c$$$      subroutine vacuum_cyl(g_def,igr,bcnvn_int,gbcnvn_int,bcov_ext)
      subroutine vacuum_cyl(g_def,igr,bcnvn_int,bcov_ext)
      
c     -----------------------------------------------------------------
c     Performs RW BC. Inputs normal internal B-component, outputs
c     tangential external B-components
c     -----------------------------------------------------------------

      implicit none

c     Call variables
        
      type(grid_mg_def),pointer :: g_def

      integer :: igr
      
      real(8) :: bcnvn_int(0:,0:),bcov_ext(0:,0:,:)!,gbcnvn_int(0:,0:)

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

      real(8),allocatable,dimension(:,:,:) :: Bn,Bnl!,gBnl,gBn

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

c$$$      allocate(gBnl(1,ny,nz),gBn(1,nyg,nzg))
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
c$$$      gBnl(1,:,:) = gbcnvn_int(1:ny,1:nz)
c$$$
c$$$      !Gather gradient in parallel along boundary
c$$$#if defined(petsc)
c$$$      call find_global_nobc(gBnl,gBn,mpi_comm=g_def%MPI_COMM_YZ)
c$$$#else
c$$$      call find_global_nobc(gBnl,gBn)
c$$$#endif
c$$$
c$$$      !Perform 2D FFT
c$$$cc>>>TODO      call fft2d(nyg,gBn(1,:,1),theta,nyg,gBr_m,gBi_m,mm,0,0)

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

c$$$      deallocate(gBnl,gBn)
c$$$      deallocate(Br_m,Bi_m)
c$$$      deallocate(gBr_m,gBi_m)
c$$$      deallocate(phir_m_np,phii_m_np)
      
      end subroutine vacuum_cyl

      end module rw_bc
