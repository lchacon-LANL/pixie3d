c module rw_bc
c ####################################################################
      module rw_bc

        use local_BCS_variables
      
        !Resistive wall vars
        logical :: res_wall_bc   =.false.
     .            ,rw_bcs(6)     =.false.

        real(8) :: tau_w=1d0,r_v=1.5d0

        real(8),pointer,dimension(:,:,:,:) :: E_bc

        logical,private :: debug_io=.false.
     .                    ,rw_dbg  =.false.

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

c$$$c     apply_rw_bc_A
c$$$c     ###############################################################
c$$$      subroutine apply_rw_bc_A(g_def,a_covn,a_cov,E_bc)
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$        
c$$$        type(grid_mg_def),pointer :: g_def
c$$$
c$$$        real(8) :: a_cov (0:,0:,0:,:)
c$$$     $            ,a_covn(0:,0:,0:,:)
c$$$     $            ,E_bc  (0:,0:,0:,:)
c$$$
c$$$c     Local variables
c$$$
c$$$        integer :: igr,nnx,nny,nnz,i,j,k
c$$$        
c$$$        real(8),pointer,dimension(:,:,:,:) :: v0,b_covn
c$$$        
c$$$c     Begin program
c$$$
c$$$        if (my_rank == 0.and.ilevel > 0.and.save_tn.and.time/=0d0)
c$$$     $       write (*,*) "Apply RW BC..."
c$$$        
c$$$        if (g_def%ihi(1) /= nxd) return
c$$$      
c$$$        igr = 1
c$$$        
c$$$        nnx = g_def%nxv(igr) 
c$$$        nny = g_def%nyv(igr)
c$$$        nnz = g_def%nzv(igr)
c$$$
c$$$c     Find electric field @ wall (store in (nx+1) ghost cell)
c$$$
c$$$        if (save_tn.and.(time/=0d0)) then
c$$$          allocate(b_covn(0:nnx+1,0:nny+1,0:nnz+1,3))
c$$$          
c$$$          !Find Bcov
c$$$          do k = 0,nnz+1
c$$$            do j = 0,nny+1
c$$$              do i = 0,nnx+1
c$$$                b_covn(i,j,k,:) =
c$$$     .                 XformToCurv_ijk(g_def,i,j,k,igr,B0,.true.)
c$$$     .           + curl(g_def,i,j,k,nnx,nny,nnz,igr,igr,igr,a_covn)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$          b_covn = XformToCov(g_def,igr,b_covn)
c$$$
c$$$          !Find E_bc
c$$$          call rw_Ebc(g_def,b_covn,E_bc)
c$$$
c$$$          deallocate(b_covn)
c$$$        endif
c$$$
c$$$c     Update old-time vector potential
c$$$
c$$$        allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))
c$$$
c$$$        v0 = 0d0
c$$$
c$$$        !Interpolate to face
c$$$        call XferBCs(a_covn,v0,.true.)
c$$$
c$$$        !Update face value
c$$$        a_cov(nnx+1,:,:,:) = v0(nnx+1,:,:,:)-dt*E_bc(nnx+1,:,:,:)
c$$$
c$$$        deallocate(v0)
c$$$       
c$$$      end subroutine apply_rw_bc_A

c     rw_Ebc
c     #################################################################
      subroutine rw_Ebc(g_def,bcov,E_0)

c     -----------------------------------------------------------------
c     Finds tangential cov electric field components according to
c     resistive wall BCs.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      real(8),dimension(0:,0:,0:,:) :: bcov
      real(8),dimension(0:,0:,0:,:) :: E_0

c     Local variables

      integer :: i_bc,j,k,nx,ny,nz,nn,mm,my_rank_y,igr

      real(8) :: phinm   (0:size(bcov,2)-1
     .                   ,0:size(bcov,3)-1)
     .          ,bcovext (0:size(bcov,2)-1
     .                   ,0:size(bcov,3)-1,3)
c$$$     .          ,bcnvext (0:size(bcov,2)-1
c$$$     .                   ,0:size(bcov,3)-1,3)
     .          ,bcovint (0:size(bcov,2)-1
     .                   ,0:size(bcov,3)-1,3)
     .          ,gbcovint(0:size(bcov,2)-1
     .                   ,0:size(bcov,3)-1,3)
c$$$     .          ,bcnvintnp(0:size(bcov,2)-1
c$$$     .                    ,0:size(bcov,3)-1,3)
c$$$     .          ,bcnvextnp(0:size(bcov,2)-1
c$$$     .                    ,0:size(bcov,3)-1,3)
     .          ,bcovextnp(0:size(bcov,2)-1
     .                    ,0:size(bcov,3)-1,3)

      real(8) :: dum,dh,bcov_avg(2:3),send_buf(2),rec_buf(2),nfactor

c     Begin program

      igr = 1
      
      if (g_def%ihi(igr) /= nxd) return
      
      nx = g_def%nxv(igr) 
      ny = g_def%nyv(igr)
      nz = g_def%nzv(igr)

      i_bc = nx + 1
      
#if defined(petsc)
      call MPI_Comm_rank(g_def%MPI_COMM_Y,my_rank_y,mpierr)
#else
      my_rank_y = 0
#endif
      
cc      debug_io = old_time_eval.and.rw_dbg.and.my_rank_y==0
      debug_io = rw_dbg.and.my_rank_y==0

c     Compute INTERNAL tangential magnetic field components, EXCLUDING m=n=0 mode (CRITICAL)

      !Extrapolate to face from interior
      bcovint (:,:,:)= 0.5*bcov(i_bc  ,:,:,:)
     .                +0.5*bcov(i_bc-1,:,:,:)
      gbcovint(:,:,:)=(bcov(i_bc  ,:,:,:)
     .                -bcov(i_bc-1,:,:,:))
     .                /(g_def%xg(nxd+1)-g_def%xg(nxd))

      !Perform angular average
      bcov_avg(2) = sum(bcovint(1:ny,1:nz,2))
      bcov_avg(3) = sum(bcovint(1:ny,1:nz,3))

#if defined(petsc)
      send_buf = bcov_avg(2:3)
      call MPI_Allreduce(send_buf,rec_buf,2,MPI_DOUBLE_PRECISION
     .                  ,MPI_SUM,g_def%MPI_COMM_YZ,mpierr)
      bcov_avg = rec_buf
#endif

      bcov_avg = bcov_avg/(nyd*nzd)
      
      !Subtract average
      bcovint(:,:,2) = bcovint(:,:,2) - bcov_avg(2)
      bcovint(:,:,3) = bcovint(:,:,3) - bcov_avg(3)

c     Compute EXTERNAL tangential magnetic field components from NORMAL internal ones

      select case(coords)
      case('hel')
        call rw_hel(g_def,bcovint(:,:,1),gbcovint(:,:,1),bcovext
     $             ,bcovextnp)
cc      case('cyl')
cc      case('tor')
      case default
        call pstop('rw_bc'
     .       ,'Resistive wall BC not implemented for '//coords)
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

c     End program

      end subroutine rw_Ebc

c     rw_hel
c     ####################################################################
      subroutine rw_hel(g_def,bcovn_int,gbcovn_int,bcov_ext,bcov_ext_np)
      
c     -----------------------------------------------------------------
c     Performs RW BC. Inputs normal internal B-component, outputs
c     tangential external B-components at old and new times
c     -----------------------------------------------------------------

      implicit none

c     Call variables
        
      type(grid_mg_def),pointer :: g_def

      real(8) :: bcovn_int(0:,0:),gbcovn_int(0:,0:)
      real(8) :: bcov_ext(0:,0:,:),bcov_ext_np(0:,0:,:)

c     Local variables
      
      integer :: m,n,jglobal,i,j,m1,mm2,mm1,nx,ny,nz
      real(8) :: theta(nyd)
     .          ,Bn(1,nyd,1),gBn(1,nyd,1)
     .          ,Br_m(nyd),Bi_m(nyd),Brn_m(nyd),Bin_m(nyd)
     .          ,gBr_m(nyd),gBi_m(nyd),mm(nyd)
     .          ,phir_m(nyd),phii_m(nyd)
     .          ,phir_m_np(nyd),phii_m_np(nyd)
     .          ,d2phir_m(nyd),d2phii_m(nyd)
     .          ,dBr(nyd),Br(nyd),Bth(nyd),Bz(nyd)
      real(8) :: dImv,dIm,Im,d2Im,dKmv,dKm,Km,d2Km,pv,pw,RR,kk,rm,rbc
     .          ,efac,gamma,dum,mmm

      real(8),allocatable,dimension(:,:,:) :: Bnl,gBnl

c     Begin program
      
      nx = size(bcov,1)-2
      ny = size(bcov,2)-2
      nz = size(bcov,3)-2

      allocate(Bnl(1,ny,1),gBnl(1,ny,1))

c     Global angular mesh
      
      theta = g_def%yg(1:nyd)

c     Process Bnormal at plasma wall (old time) to determine potential BC

      Bnl(1,:,1) = bcovn_int(1:ny,1)

      !Gather Bnormal in parallel along boundary
#if defined(petsc)
      call find_global_nobc(Bnl ,Bn,mpi_comm=g_def%MPI_COMM_Y)
#else
      call find_global_nobc(Bnl ,Bn)
#endif

      !Perform FFT
      call fft(nyd,Bn(1,:,1),theta,nyd,Brn_m,Bin_m,mm,0,0)

c     Process gradient of Bnormal (new time)
      
      !Compute gradient of Bn
      gBnl(1,:,1) = gbcovn_int(1:ny,1)

      !Gather gradient in parallel along boundary
#if defined(petsc)
      call find_global_nobc(gBnl,gBn,mpi_comm=g_def%MPI_COMM_Y)
#else
      call find_global_nobc(gBnl,gBn)
#endif

      !Perform FFT
      call fft(nyd,gBn(1,:,1),theta,nyd,gBr_m,gBi_m,mm,0,0)

c     Find VACUUM potential solution Phi outside of resistive wall

      mmm = g_def%params(1)
      kk  = g_def%params(2)

      if (kk == 0d0) then
        do m=1,nyd
          m1   = abs(nint(mm(m)))
          if (m1 == 0) then
            phir_m(m) = 0d0 ; phii_m(m) = 0d0
          else
            rm = r_v**(2*mm(m))

            !Evolve normal B component
            if (.not.source_eval) then
              gamma = 1 + mm(m)*(rm+1d0)/(rm-1d0)
              efac  = exp(-gamma*dt/tau_w)

              Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
              Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)

c$$$              if(m1==nint(mmm).and.debug_io) then
c$$$                if (Brn_m(m)*Br_m(m) > 0d0) then  
c$$$                  dum=log(Br_m(m)/Brn_m(m))*tau_w/dt
c$$$                else
c$$$                  dum=0d0
c$$$                endif
c$$$                write (*,*)"Br norm. decay k=0 =",dum
c$$$              endif
            else
              Br_m(m) = Brn_m(m)
              Bi_m(m) = Bin_m(m)
            endif

            !Find phi_nm
            dum = (rm+1d0)/((1d0-rm)*mm(m))
            phir_m_np(m) = Br_m(m)*dum
            phii_m_np(m) = Bi_m(m)*dum
            phir_m   (m) = Brn_m(m)*dum
            phii_m   (m) = Bin_m(m)*dum
          endif
        enddo
      else
        pw = abs(kk)
        pv = abs(kk)*r_v
        do m=1,nyd
          m1   = abs(nint(mm(m)))
          mm1  = abs(m1-1)
          mm2  = abs(m1-2)

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

          !Evolve normal B component
          if (.not.source_eval) then
            gamma =-pw*(dImv*d2Km-dKmv*d2Im)/(dImv*dKm -dKmv*dIm)
            efac  = exp(-gamma*dt/tau_w)
            Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
            Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)

c$$$            if(m1==nint(mmm).and.debug_io) then
c$$$              if (Brn_m(m)*Br_m(m) > 0d0) then
c$$$                dum=log(Br_m(m)/Brn_m(m))*tau_w/dt
c$$$              else
c$$$                dum=0d0
c$$$              endif
c$$$              write (*,*) m,"Br norm. decay k/=0, m =",mmm," is",dum
c$$$            endif
          else
            Br_m(m) = Brn_m(m)
            Bi_m(m) = Bin_m(m)
          endif

          !Find phi_nm
          dum = (dImv*Km-dKmv*Im)/((dImv*dKm-dKmv*dIm)*pw)
          phir_m_np(m) = Br_m(m)*dum
          phii_m_np(m) = Bi_m(m)*dum
          phir_m   (m) = Brn_m(m)*dum
          phii_m   (m) = Bin_m(m)*dum
        enddo
      endif

c     Find old-time EXTERNAL magnetic field components outside RW

      !Bth [multiply by -sqrt(-1)*m; minus sign from convention of iFFT]
      Brn_m = mm*phii_m
      Bin_m =-mm*phir_m
      call ifft(nyd,Bth,theta,nyd,Brn_m,Bin_m,mm,0)

c$$$      !Bz [multiply by -sqrt(-1)*k; minus sign from convention of iFFT]
c$$$      Br_m = kk*phii_m
c$$$      Bi_m =-kk*phir_m
c$$$      call ifft(nyd,Bz,theta,nyd,Br_m,Bi_m,mm,0)

      !Get LOCAL helical (cov) magnetic field components @ r=1 
      rbc = 1d0

      do j = jlom,jhip
        jglobal = j + g_def%jlo(1) - 1
        if (jglobal > nyd) jglobal = jglobal - nyd
        if (jglobal < 1  ) jglobal = nyd + jglobal
        bcov_ext(j,:,1) = bcovn_int(j,:)
        bcov_ext(j,:,2) = rbc/mmm*Bth(jglobal)
        bcov_ext(j,:,3) = 0d0 !Bz(jglobal) - kk*bcov_ext(j,:,2)
      enddo

c     Find new-time EXTERNAL magnetic field components outside RW

      !Br
      call ifft(nyd,Br,theta,nyd,Br_m,Bi_m,mm,0)

      !Bth [multiply by -sqrt(-1)*m; minus sign from convention of iFFT]
      Brn_m = mm*phii_m_np
      Bin_m =-mm*phir_m_np
      call ifft(nyd,Bth,theta,nyd,Brn_m,Bin_m,mm,0)

c$$$      !Bz [multiply by -sqrt(-1)*k; minus sign from convention of iFFT]
c$$$      Br_m = kk*phii_m_np
c$$$      Bi_m =-kk*phir_m_np
c$$$      call ifft(nyd,Bz,theta,nyd,Br_m,Bi_m,mm,0)

      !Get LOCAL helical (cov) magnetic field components @ r=1 
      rbc = 1d0

      do j = jlom,jhip
        jglobal = j + g_def%jlo(1) - 1
        if (jglobal > nyd) jglobal = jglobal - nyd
        if (jglobal < 1  ) jglobal = nyd + jglobal
        bcov_ext_np(j,:,1) = Br(jglobal)
        bcov_ext_np(j,:,2) = rbc/mmm*Bth(jglobal)
        bcov_ext_np(j,:,3) = 0d0 !Bz(jglobal) - kk*bcov_ext_np(j,:,2)
      enddo

      if (rw_dbg) then
        dum = 0d0
        do j = jlo,jhi
          dum = dum + (bcov_ext_np(j,1,2)-bcov_ext(j,1,2))
     $                /bcov_ext(j,1,2)/dt*tau_w
        enddo
        dum = psum((/dum/))/nyd
        
        if (debug_io)
     $    write(*,*) "Bt norm. decay rate k/=0, m =",nint(mmm)," is",dum
      endif
      
      deallocate(Bnl,gBnl)

      end subroutine rw_hel

c$$$c     rw_cyl
c$$$c     ####################################################################
c$$$      subroutine rw_cyl
c$$$
c$$$      integer :: i_bc,m,n,jglobal,i,j,m1,mm2,mm1
c$$$      real(8) :: theta(nyd)
c$$$     .          ,Bn(1,nyd,1),gBn(1,nyd,1)
c$$$     .          ,Bnl(1,ny,1),gBnl(1,ny,1)
c$$$     .          ,Br_m(nyd),Bi_m(nyd)
c$$$     .          ,gBr_m(nyd),gBi_m(nyd),mm(nyd)
c$$$     .          ,phir_m(nyd),phii_m(nyd)
c$$$     .          ,d2phir_m(nyd),d2phii_m(nyd)
c$$$     .          ,dBr(nyd),Bth(nyd),Bz(nyd)
c$$$      real(8) :: dImv,dIm,Im,d2Im,dKmv,dKm,Km,d2Km,pv,pw,RR,kk,rm,rbc
c$$$     .          ,efac,gamma,dum
c$$$
c$$$c     Begin program
c$$$      
c$$$      i_bc = nx + 1
c$$$
c$$$c     Global angular mesh
c$$$      
c$$$      theta = g_def%yg(1:nyd)
c$$$
c$$$c     Process Bnormal (old time)
c$$$
c$$$#if !defined(plotting)
c$$$      if ((old_time_eval.or.source_eval)) then
c$$$#endif
c$$$        if (.not.allocated(Brn_m)) allocate(Brn_m(nyd),Bin_m(nyd))
c$$$        
c$$$        !Extrapolate to r=1 old-time Bn
c$$$        Bnl(1,:,1) = 1.5*bcov_n(i_bc-1,1:ny,1,1)
c$$$     .              -0.5*bcov_n(i_bc-2,1:ny,1,1)
c$$$
c$$$        !Gather Bnormal in parallel along boundary
c$$$#if defined(petsc)
c$$$        call find_global_nobc(Bnl ,Bn ,mpi_comm=g_def%MPI_COMM_Y)
c$$$#else
c$$$        call find_global_nobc(Bnl ,Bn )
c$$$#endif
c$$$
c$$$        !Perform FFT
c$$$        call fft(nyd,Bn (1,:,1),theta,nyd,Brn_m ,Bin_m ,mm,0,0)
c$$$
c$$$#if !defined(plotting)
c$$$      endif
c$$$#endif
c$$$      
c$$$c     Process gradient of Bnormal (new time)
c$$$      
c$$$      !Compute new-time gradient of Bn
c$$$      dh = g_def%xg(nxd)-g_def%xg(nxd-1)
c$$$      gBnl(1,:,1) =(bcov(i_bc-1,1:ny,1,1)
c$$$     .             -bcov(i_bc-2,1:ny,1,1))/dh
c$$$
c$$$      !Gather gradient in parallel along boundary
c$$$#if defined(petsc)
c$$$      call find_global_nobc(gBnl,gBn,mpi_comm=g_def%MPI_COMM_Y)
c$$$#else
c$$$      call find_global_nobc(gBnl,gBn)
c$$$#endif
c$$$
c$$$      !Perform FFT
c$$$      call fft(nyd,gBn(1,:,1),theta,nyd,gBr_m,gBi_m,mm,0,0)
c$$$
c$$$      !Find potential solution Phi outside of resistive wall
c$$$      kk = g_def%params(2)
c$$$
c$$$      if (kk == 0d0) then
c$$$        do m=1,nyd
c$$$          m1   = abs(nint(mm(m)))
c$$$          if (mm(m) == 0) then
c$$$            phir_m(m) = 0d0 ; phii_m(m) = 0d0
c$$$          else
c$$$            rm = r_v**(2*mm(m))
c$$$            
c$$$            !Evolve normal B component
c$$$            if (.not.source_eval) then
c$$$              gamma = 1 + mm(m)*(rm+1d0)/(rm-1d0)
c$$$              efac  = exp(-gamma*dt/tau_w)
c$$$
c$$$              Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
c$$$              Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)
c$$$
c$$$              if(m1==nint(g_def%params(1)).and.old_time_eval.and.rw_dbg)
c$$$     .        then
c$$$                if (Brn_m(m) /= 0d0.and. Br_m(m) < 0d0) then  
c$$$                  dum=log(Br_m(m)/Brn_m(m))*tau_w/dt
c$$$                else
c$$$                  dum=0d0
c$$$                endif
c$$$                write (*,*)"Br norm. decay k=0 =",dum
c$$$              endif
c$$$            endif
c$$$
c$$$            !Find phi_nm
c$$$            phir_m(m) = Br_m(m)/mm(m)*(rm+1d0)/(1d0-rm)
c$$$            phii_m(m) = Bi_m(m)/mm(m)*(rm+1d0)/(1d0-rm)
c$$$          endif
c$$$        enddo
c$$$      else
c$$$        pw = abs(kk)
c$$$        pv = abs(kk)*r_v
c$$$        do m=1,nyd
c$$$          m1   = abs(nint(mm(m)))
c$$$          mm1  = abs(m1-1)
c$$$          mm2  = abs(m1-2)
c$$$
c$$$          !Bessel function evals
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
c$$$
c$$$          !Evolve normal B component
c$$$          if (.not.source_eval) then
c$$$            gamma =-pw*(dImv*d2Km-dKmv*d2Im)/(dImv*dKm -dKmv*dIm)
c$$$            efac  = exp(-gamma*dt/tau_w)
c$$$            Br_m(m) = efac*Brn_m(m)-(1d0-efac)/gamma*gBr_m(m)
c$$$            Bi_m(m) = efac*Bin_m(m)-(1d0-efac)/gamma*gBi_m(m)
c$$$
c$$$            if(m1==nint(g_def%params(1)).and.old_time_eval.and.rw_dbg)
c$$$     .      then
c$$$              if (Brn_m(m) < 0d0.and.Br_m(m) < 0d0) then
c$$$                dum=log(Br_m(m)/Brn_m(m))*tau_w/dt
c$$$              else
c$$$                dum=0d0
c$$$              endif
c$$$              write (*,*)"Br norm. decay k/=0 =",dum
c$$$            endif
c$$$          endif
c$$$          
c$$$          !Find phi_nm
c$$$          phir_m(m)=Br_m(m)/pw*(dImv*Km-dKmv*Im)/(dImv*dKm-dKmv*dIm)
c$$$          phii_m(m)=Bi_m(m)/pw*(dImv*Km-dKmv*Im)/(dImv*dKm-dKmv*dIm)
c$$$c$$$          d2phir_m(m)=-Br_m(m)*gamma
c$$$c$$$          d2phii_m(m)=-Bi_m(m)*gamma
c$$$        enddo
c$$$      endif
c$$$
c$$$c$$$      !Reconstruct GLOBAL cyl magnetic field components from grad(Phi) outside plasma wall
c$$$c$$$      !Bth [multiply by -sqrt(-1)*m; minus sign from convention of iFFT]
c$$$c$$$      call ifft(nyd,dBr,theta,nyd,d2phir_m,d2phii_m,mm,0)
c$$$      
c$$$      !Bth [multiply by -sqrt(-1)*m; minus sign from convention of iFFT]
c$$$      Br_m = mm*phii_m
c$$$      Bi_m =-mm*phir_m
c$$$      call ifft(nyd,Bth,theta,nyd,Br_m,Bi_m,mm,0)
c$$$
c$$$      !Bz [multiply by -sqrt(-1)*k; minus sign from convention of iFFT]
c$$$      m    = g_def%params(1)
c$$$      Br_m = mm*kk/m*phii_m  !mm vector, m scalar; for helical only
c$$$      Bi_m =-mm*kk/m*phir_m  !mm vector, m scalar; for helical only
c$$$      call ifft(nyd,Bz,theta,nyd,Br_m,Bi_m,mm,0)
c$$$
c$$$      !Get LOCAL helical (cov) magnetic field components @ r=1 (store in (nx+1))
c$$$      rbc = 1d0
c$$$
c$$$      do k = klom,khip
c$$$        do j = jlom,jhip
c$$$          jglobal = j + g_def%jlo(1) - 1
c$$$          if (jglobal > nyd) jglobal = jglobal - nyd
c$$$          if (jglobal < 1  ) jglobal = nyd + jglobal
c$$$c$$$          dbcnv1ext(j,k) = rbc*dBr(jglobal)
c$$$          bcovext(j,k,2) = rbc*Bth(jglobal)
c$$$          bcovext(j,k,3) = Bz (jglobal) - rbc*kk/m*Bth(jglobal)
c$$$        enddo
c$$$      enddo
c$$$    
c$$$      end subroutine rw_cyl

      end module rw_bc
