c postProcessSol
c####################################################################
      subroutine postProcessSol(vnm,vn,vnp)

c--------------------------------------------------------------------
c     Postprocess solution (e.g. filtering). In this instance, we
c     do divergence cleaning of the magnetic field.
c--------------------------------------------------------------------

      use nlfunction_setup

      use imposeBCinterface

      use PC_var, only: pc_B_solve,pc_type,pc_divclean

      use xdraw_io

      use B_tools

      use par_int

      use PC_setup, only: pc_as_solver

      use local_BCS_variables, ONLY: order_bc

      use rw_bc, ONLY: res_wall_bc,E_bc,rw_Ebc
      
      implicit none

c Call variables

      type(var_array),pointer :: vnm,vn,vnp

c Local variables

      integer :: i,j,k,igx,igy,igz,nx,ny,nz,ieq,igrid
      real(8) :: tol,acc

c Interfaces

      INTERFACE
        subroutine perturbEquilibrium_ppnch_bc(varray,br_pert,br_phase)
        use var_setup
        type(var_array),pointer :: varray
        real(8) :: br_pert,br_phase
        end subroutine perturbEquilibrium_ppnch_bc

        function tau_est(g_def,vf) result(tol)
        use var_setup
        real(8) :: tol
        type(var_array),pointer :: vf
        type(grid_mg_def),pointer :: g_def
        end function tau_est
      END INTERFACE

c Begin program

      igrid = 1

      !Setup grid info
      igx = igrid
      igy = igrid
      igz = igrid

      nx = gv%gparams%nxv(igx)
      ny = gv%gparams%nyv(igy)
      nz = gv%gparams%nzv(igz)

c Add random noise in EOM (to break symmetry)

      if (noise) then
        if (nc_eom_v) then
          do k=1,nz
            do j=1,ny
              do i=1,nx
                do ieq=IVX,IVZ
                  call random_number(acc) !Between 0 and 1
                  acc = noise_lev*sqrt(dt)*(1-2.*acc)
                  vnp%array_var(ieq)%array =vnp%array_var(ieq)%array
     .                                     +acc
                enddo
              enddo
            enddo
          enddo
        else
          do k=1,nz
            do j=1,ny
              do i=1,nx
                do ieq=IVX,IVZ
                  call random_number(acc) !Between 0 and 1
                  acc = noise_lev*sqrt(dt)*(1-2.*acc)
                  vnp%array_var(ieq)%array=vnp%array_var(ieq)%array
     .                                    +acc*vnp%array_var(IRHO)%array
                enddo
              enddo
            enddo
          enddo
        endif
      endif

c -------------------------------------------------------------------
c SAMRAI does not execute code below
c -------------------------------------------------------------------

#if !defined(samrai)

c$$$c Filter single mode in helical geometry (only on V for now)
c$$$
c$$$      if (coords == 'hel'.and.post_filter_hel) then
c$$$        call post__filter_V
c$$$      endif
      
c Smooth magnetic and velocity fields

      if (post_smooth_B) then
        call post__smooth_B
        call post__smooth_V
      endif
      
c Divergence-clean magnetic field when PC is on

#if !defined(vec_pot)
      if (post_divclean.or.res_wall_bc) call post__divclean
#endif

c Perform Lagrangian temperature step

      if (lagrangian) call post__lagrangian

c Add rotating magnetic perturbation

#if !defined(vec_pot)
      if (br_pert_freq /= 0d0) then

        !Substract old one
        call perturbEquilibrium_ppnch_bc(vnp,-br_pert_bc
     .                                  ,br_pert_phase)

        br_pert_phase = br_pert_phase + br_pert_freq*dt

        !Add new one
        call perturbEquilibrium_ppnch_bc(vnp, br_pert_bc
     .                                  ,br_pert_phase)
      endif
#endif

c Update BCs

      call applyBC(igx,vnp,gv%aux)

c Compute new Newton update tolerance based on tau error estimator

c$$$      tol = tau_est(gv%gparams,vnp)
c$$$
c$$$      if (my_rank == 0) write (*,*) 'New Newton update tolerance=',tol

c Dump parareal convergence info

#if defined(pit)
      allocate(array(0:nx+1,0:ny+1,0:nz+1))

      !Magnetic energy
      array = scalarProduct(igrid,bcnv,bcov)
      Em = 0.5*integral(igx,igy,igz,nx,ny,nz,array,average=.false.)

      !Kinetic energy
      array = scalarProduct(igrid,vcnv,vcov)
      Ek = 0.5*integral(igx,igy,igz,nx,ny,nz,array,average=.false.)

      !Thermal energy
      array = rho*tmp

      if (gamma /= 1d0) then
        Et = 1d0/(gamma-1)
     .      *integral(igx,igy,igz,nx,ny,nz,array,average=.false.)
      else
        Et = 0d0
      endif

      if (my_rank == 0) then
        pit_file = 'pixie3d_pit_conv.txt'

        if (itime==inewtime) then
           open(unit=upit,file=trim(pit_file)
     $         ,status='unknown')
        else
           open(unit=upit,file=trim(pit_file)
     $         ,status='unknown',position='append')
        endif

        if (itime == inewtime) write (upit,*) "3" !Number of diagnostics
        write(upit,*) time,Em,Ek,Et
        close(upit)
      endif

      deallocate(array)

      !Output PARAREAL checkpointing file
      if (    (tmax .gt.0d0.and.time .ge.(tmax-1d-9*dt))
     .    .or.(numtime.ge.0.and.itime.ge.(numtime+inewtime-1))) then
        call writeRecordFile("pixie3d_pit_chkpnt.bin",0,0d0,dt,u_0
     .                      ,init=.true.)
        call writeRecordFile("pixie3d_pit_chkpnt.bin",itime,time
     .                      ,dt,vnp)
      endif

#endif

c End program

#endif

      contains

c$$$c     post__filter_V
c$$$c     ###############################################################
c$$$      subroutine post__filter_V
c$$$
c$$$        implicit none
c$$$
c$$$        integer :: ieq,i,j,k,jg,jlog
c$$$        real(8) :: ac(nx),as(nx),lsum(nx,2),gsum(nx,2),mm
c$$$        real(8),allocatable,dimension(:,:,:,:) :: dv
c$$$
c$$$c     Begin program
c$$$
c$$$        if (my_rank == 0) then
c$$$          if (ilevel > 0) write (*,*)
c$$$          write (*,*) 'Performing filtering of V in postproc...'
c$$$        endif
c$$$
c$$$c     Find dV=V^(n+1)-V^n
c$$$
c$$$        allocate(dv(0:nx+1,0:ny+1,0:nz+1,3))
c$$$
c$$$        dv(:,:,:,1) = vnp%array_var(IVX)%array
c$$$     .               -vn %array_var(IVX)%array
c$$$        dv(:,:,:,2) = vnp%array_var(IVY)%array
c$$$     .               -vn %array_var(IVY)%array
c$$$        dv(:,:,:,3) = vnp%array_var(IVZ)%array
c$$$     .               -vn %array_var(IVZ)%array
c$$$
c$$$c     Filter by Fourier transform
c$$$
c$$$        mm = gv%gparams%params(1)
c$$$
c$$$        jlog = gv%gparams%jlo(igy)
c$$$
c$$$        do ieq=1,3
c$$$
c$$$          k = 1
c$$$
c$$$          ac = 0d0 ; as = 0d0
c$$$          
c$$$          !Perform Fourier transform in theta
c$$$          do j = 1,ny
c$$$            jg = j + jlog - 1   !Global index
c$$$            ac = ac + dv(1:nx,j,k,ieq)*cos(-2*pi*mm*(jg-1)/nyd)
c$$$            as = as + dv(1:nx,j,k,ieq)*sin(-2*pi*mm*(jg-1)/nyd)
c$$$          enddo
c$$$          
c$$$#if defined(petsc)
c$$$          lsum(:,1) = ac
c$$$          lsum(:,2) = as
c$$$          call MPI_Allreduce(lsum,gsum,2*nx,MPI_DOUBLE_PRECISION
c$$$     .                      ,MPI_SUM,gv%gparams%MPI_COMM_Y,mpierr)
c$$$          ac = gsum(:,1)
c$$$          as = gsum(:,2)
c$$$#endif
c$$$
c$$$          !Inverse Fourier transform
c$$$          do j = 1,ny
c$$$            jg = j + jlog - 1   !Global index
c$$$            dv(1:nx,j,k,ieq) =2*(ac*cos(2*pi*mm*(jg-1)/nyd) !Factor of 2 to
c$$$     .                          -as*sin(2*pi*mm*(jg-1)/nyd))/nyd !account for (N-m) contrib
c$$$          enddo
c$$$
c$$$        enddo
c$$$
c$$$c     Find new V^(n+1)
c$$$
c$$$c$$$        write (*,*) "1"
c$$$c$$$     .              ,sqrt(sum((vnp%array_var(IVX)%array(1:nx,1:ny,1:nz)
c$$$c$$$     .                        -vn %array_var(IVX)%array(1:nx,1:ny,1:nz)
c$$$c$$$     .                       -dv(1:nx,1:ny,1:nz,1))**2))
c$$$c$$$     .              /sqrt(sum((dv(1:nx,1:ny,1:nz,1))**2))
c$$$c$$$        write (*,*) "2"
c$$$c$$$     .              ,sqrt(sum((vnp%array_var(IVY)%array(1:nx,1:ny,1:nz)
c$$$c$$$     .                        -vn %array_var(IVY)%array(1:nx,1:ny,1:nz)
c$$$c$$$     .                       -dv(1:nx,1:ny,1:nz,2))**2))
c$$$c$$$     .              /sqrt(sum((dv(1:nx,1:ny,1:nz,2))**2))
c$$$c$$$        write (*,*) "3"
c$$$c$$$     .              ,sqrt(sum((vnp%array_var(IVZ)%array(1:nx,1:ny,1:nz)
c$$$c$$$     .                        -vn %array_var(IVZ)%array(1:nx,1:ny,1:nz)
c$$$c$$$     .                       -dv(1:nx,1:ny,1:nz,3))**2))
c$$$c$$$     .              /sqrt(sum((dv(1:nx,1:ny,1:nz,3))**2))
c$$$c$$$        k = 1
c$$$c$$$        do j = 1,1 !ny
c$$$c$$$          do i = 1,nx
c$$$c$$$            write (*,*) "1"
c$$$c$$$     .                 ,(vnp%array_var(IVX)%array(i,j,k)
c$$$c$$$     .                  -vn %array_var(IVX)%array(i,j,k))
c$$$c$$$     .                 /dv(i,j,k,1)
c$$$c$$$            write (*,*) "2"
c$$$c$$$     .                 ,(vnp%array_var(IVY)%array(i,j,k)
c$$$c$$$     .                  -vn %array_var(IVY)%array(i,j,k))
c$$$c$$$     .                 /dv(i,j,k,2)
c$$$c$$$            write (*,*) "3"
c$$$c$$$     .                 ,(vnp%array_var(IVZ)%array(i,j,k)
c$$$c$$$     .                  -vn %array_var(IVZ)%array(i,j,k))
c$$$c$$$     .                 /dv(i,j,k,3)
c$$$c$$$          enddo
c$$$c$$$        enddo
c$$$c$$$        
c$$$c$$$        stop
c$$$
c$$$        vnp%array_var(IVX)%array = dv(:,:,:,1)
c$$$     .                            +vn%array_var(IVX)%array
c$$$        vnp%array_var(IVY)%array = dv(:,:,:,2)
c$$$     .                            +vn%array_var(IVY)%array
c$$$        vnp%array_var(IVZ)%array = dv(:,:,:,3)
c$$$     .                            +vn%array_var(IVZ)%array
c$$$
c$$$        deallocate(dv)
c$$$        
c$$$      end subroutine post__filter_V

c     post__smooth_B
c     ###############################################################
      subroutine post__smooth_B

        implicit none

        integer :: nn,bcnd(6,3),filter_it,ig,jg,kg,ii,guess,iout
        real(8) :: opt_om,dvol
        real(8),allocatable,dimension(:,:) :: xrhs,x

        external  sm_mtvc_B

c     Begin program
        
        nn = nx*ny*nz

        if (my_rank == 0) then
          if (ilevel > 0) write (*,*)
#if defined(vec_pot)
          write (*,*) 'Performing A smoothing in postproc...'
#else
          write (*,*) 'Performing B smoothing in postproc...'
#endif
        endif

c     Find db=B^(n+1)-B^n

        allocate(x(nn,3),xrhs(nn,3))

#if defined(vec_pot)
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

              dvol = gv%gparams%gmetric%grid(igx)%dlvol(i,j,k)

              x   (ii,1) = vnp%array_var(IAX)%array(i,j,k)
     .                    -vn %array_var(IAX)%array(i,j,k)
              x   (ii,2) = vnp%array_var(IAY)%array(i,j,k)
     .                    -vn %array_var(IAY)%array(i,j,k)
              x   (ii,3) = vnp%array_var(IAZ)%array(i,j,k)
     .                    -vn %array_var(IAZ)%array(i,j,k)

              xrhs(ii,:) = x(ii,:)*dvol
            enddo
          enddo
        enddo

        do ieq=IAX,IAZ
          bcnd(:,ieq-IAX+1) = vnp%array_var(ieq)%bconds(:)
        enddo
#else
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

              dvol = gv%gparams%gmetric%grid(igx)%dlvol(i,j,k)

              x   (ii,1) = vnp%array_var(IBX)%array(i,j,k)
     .                    -vn %array_var(IBX)%array(i,j,k)
              x   (ii,2) = vnp%array_var(IBY)%array(i,j,k)
     .                    -vn %array_var(IBY)%array(i,j,k)
              x   (ii,3) = vnp%array_var(IBZ)%array(i,j,k)
     .                    -vn %array_var(IBZ)%array(i,j,k)

              xrhs(ii,:) = x(ii,:)*dvol
            enddo
          enddo
        enddo

        do ieq=IBX,IBZ
          bcnd(:,ieq-IBX+1) = vnp%array_var(ieq)%bconds(:)
        enddo
#endif

c     Smooth A or B

        filter_it = 4

        dvol = 1d0
        if (nxd > 1) dvol = dvol*(xmax-xmin)/nxd
        if (nyd > 1) dvol = dvol*(ymax-ymin)/nyd
        if (nzd > 1) dvol = dvol*(zmax-zmin)/nzd

        opt_om = dt*di/dvol

        opt_om = (1.+2.*opt_om)/(1.+4.*opt_om)

        iout = ilevel - 1
        guess = 1

        call cSolver(gv%gparams,3,nn,xrhs,x,bcnd
     .              ,igx,iout,guess,sm_mtvc_B,.true.
     .              ,mg_smooth     = 'jb'
     .              ,sm_it         = filter_it
     .              ,sm_omega      = opt_om
     .              ,sm_ncolors    = 2
     .              ,mg_smooth_only= .true.
     .              )

c     Find new A^(n+1) or B^(n+1)

#if defined(vec_pot)
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              vnp%array_var(IAX)%array(i,j,k) = x(ii,1)
     .                 +vn%array_var(IAX)%array(i,j,k)
              vnp%array_var(IAY)%array(i,j,k) = x(ii,2)
     .                 +vn%array_var(IAY)%array(i,j,k)
              vnp%array_var(IAZ)%array(i,j,k) = x(ii,3)
     .                 +vn%array_var(IAZ)%array(i,j,k)
            enddo
          enddo
        enddo
#else
        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              vnp%array_var(IBX)%array(i,j,k) = x(ii,1)
     .                 +vn%array_var(IBX)%array(i,j,k)
              vnp%array_var(IBY)%array(i,j,k) = x(ii,2)
     .                 +vn%array_var(IBY)%array(i,j,k)
              vnp%array_var(IBZ)%array(i,j,k) = x(ii,3)
     .                 +vn%array_var(IBZ)%array(i,j,k)
            enddo
          enddo
        enddo
#endif
        deallocate(x,xrhs)
      
      end subroutine post__smooth_B

c     post__smooth_V
c     ###############################################################
      subroutine post__smooth_V

        implicit none

        integer :: nn,bcnd(6,3),filter_it,ig,jg,kg,ii,guess,iout
        real(8) :: opt_om,dvol
        real(8),allocatable,dimension(:,:) :: xrhs,x

        external  sm_mtvc_V

c     Begin program
        
        nn = nx*ny*nz

        if (my_rank == 0) then
          if (ilevel > 0) write (*,*)
          write (*,*) 'Performing V smoothing in postproc...'
        endif

c     Find db=B^(n+1)-B^n

        allocate(x(nn,3),xrhs(nn,3))

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

              dvol = gv%gparams%gmetric%grid(igx)%dlvol(i,j,k)

              x   (ii,1) = vnp%array_var(IVX)%array(i,j,k)
     .                    -vn %array_var(IVX)%array(i,j,k)
              x   (ii,2) = vnp%array_var(IVY)%array(i,j,k)
     .                    -vn %array_var(IVY)%array(i,j,k)
              x   (ii,3) = vnp%array_var(IVZ)%array(i,j,k)
     .                    -vn %array_var(IVZ)%array(i,j,k)

              xrhs(ii,:) = x(ii,:)*dvol
            enddo
          enddo
        enddo

        do ieq=IVX,IVZ
          bcnd(:,ieq-IVX+1) = vnp%array_var(ieq)%bconds(:)
        enddo

c     Smooth A or B

        filter_it = 4

        dvol = 1d0
        if (nxd > 1) dvol = dvol*(xmax-xmin)/nxd
        if (nyd > 1) dvol = dvol*(ymax-ymin)/nyd
        if (nzd > 1) dvol = dvol*(zmax-zmin)/nzd

        opt_om = dt*di/dvol

        opt_om = (1.+2.*opt_om)/(1.+4.*opt_om)

        iout = ilevel - 1
        guess = 1

        call cSolver(gv%gparams,3,nn,xrhs,x,bcnd
     .              ,igx,iout,guess,sm_mtvc_V,.true.
     .              ,mg_smooth     = 'jb'
     .              ,sm_it         = filter_it
     .              ,sm_omega      = opt_om
     .              ,sm_ncolors    = 2
     .              ,mg_smooth_only= .true.
     .              )

c     Find new A^(n+1) or B^(n+1)

        do k = 1,nz
          do j = 1,ny
            do i = 1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              vnp%array_var(IVX)%array(i,j,k) = x(ii,1)
     .                 +vn%array_var(IVX)%array(i,j,k)
              vnp%array_var(IVY)%array(i,j,k) = x(ii,2)
     .                 +vn%array_var(IVY)%array(i,j,k)
              vnp%array_var(IVZ)%array(i,j,k) = x(ii,3)
     .                 +vn%array_var(IVZ)%array(i,j,k)
            enddo
          enddo
        enddo

        deallocate(x,xrhs)
      
      end subroutine post__smooth_V

c     post__divclean
c     ###############################################################
      subroutine post__divclean

        implicit none

        real(8) :: tol_divcl
        real(8),allocatable,dimension(:,:,:,:) :: db_cnv,da_cov

        if (ilevel > 0 .and. my_rank == 0) then
          write (*,*)
          if (res_wall_bc) then
            write (*,*) 'Applying RW BC in postproc...'
          else
            write (*,*) 'Performing divergence cleaning in postproc...'
          endif
        endif

c       Find db=B^(n+1)-B^n

        allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     $          ,da_cov(0:nx+1,0:ny+1,0:nz+1,3))

        db_cnv(:,:,:,1) = vnp%array_var(IBX)%array
     .                   -vn %array_var(IBX)%array
        db_cnv(:,:,:,2) = vnp%array_var(IBY)%array
     .                   -vn %array_var(IBY)%array
        db_cnv(:,:,:,3) = vnp%array_var(IBZ)%array
     .                   -vn %array_var(IBZ)%array

c       Find vector potential RW boundary condition at face

        da_cov = 0d0
        if (res_wall_bc) call XferBCs(-dt*E_bc,da_cov,.false.)
        
c       Clean db (and impose RW BC, if applicable)

        tol_divcl = rtol
        call vlap_div_clean(gv%gparams,ilevel,nx,ny,nz
     .                     ,igx,igy,igz,db_cnv,da0cov=da_cov
     .                     ,gm_it=it_divcl,rtol=tol_divcl
     .                     )!,diag_plots=.true.)
        
c       Find new B^(n+1)

        vnp%array_var(IBX)%array = db_cnv(:,:,:,1)
     .                            +vn%array_var(IBX)%array
        vnp%array_var(IBY)%array = db_cnv(:,:,:,2)
     .                            +vn%array_var(IBY)%array
        vnp%array_var(IBZ)%array = db_cnv(:,:,:,3)
     .                            +vn%array_var(IBZ)%array

c       Update helical flux and magnetic field

        bcnv(:,:,:,1) = vnp%array_var(IBX)%array
        bcnv(:,:,:,2) = vnp%array_var(IBY)%array
        bcnv(:,:,:,3) = vnp%array_var(IBZ)%array

        deallocate(db_cnv,da_cov)

      end subroutine post__divclean

c     post__lagrangian
c     ###############################################################
      subroutine post__lagrangian

        implicit none

        !Lagrangian
        real(8),allocatable,dimension(:,:,:)   :: tmp_nm,tmp_n,tmp_np
     $                                           ,chipar,array
        logical :: B_weight,do_3d
        real(8) :: G_tol
        integer :: sp_order

        !Parareal
        integer :: upit=4321
        character(1024) :: pit_file
        real(8) :: Et,Ek,Em

        allocate(tmp_nm(0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_n (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_np(0:nx+1,0:ny+1,0:nz+1)
     $          ,chipar(0:nx+1,0:ny+1,0:nz+1))

        chipar = chi_par - chi
        G_tol = 1d-8
        sp_order = 2
        B_weight = .false.
        do_3d    = .true.

        if (maxval(chipar) == 0d0) return

        if (ilevel > 0 .and. my_rank == 0) then
          write (*,*)
          write (*,*) 'Performing divergence cleaning in postproc...'
        endif

        call lag_init(gv%gparams,igrid,bcnv,2d0,ilevel
     .               ,G_tol,sp_order,B_weight,do_3d)

        !Find temp from pressure (automatically takes care of jacobian factor)
        tmp_nm = vnm%array_var(IPRS)%array
     $          /vnm%array_var(IRHO)%array
        tmp_n  = vn %array_var(IPRS)%array
     $          /vn %array_var(IRHO)%array
        tmp_np = vnp%array_var(IPRS)%array
     $          /vnp%array_var(IRHO)%array

        call lag_update_bdf2(gv%gparams,igrid,IPRS
     .                      ,vnp%array_var(IPRS)%bconds
     $                      ,tmp_nm,tmp_n,tmp_np,chipar,ilevel)
c$$$        lag_init = .true.
c$$$        E_diag   = .true.
c$$$        call lagr_update(igrid,bcnv
c$$$     .                  ,u_nm %array_var(IPRS)%array
c$$$     .                  ,vn  %array_var(IPRS)%array
c$$$     .                  ,vnp%array_var(IPRS)%array
c$$$     .                  ,chi_par,chi,lag_init,E_diag,ilevel-3)

        call lag_destroy

        !Return pressure (with jacobian)
        vnp%array_var(IPRS)%array = tmp_np*vnp%array_var(IRHO)%array

        deallocate(tmp_nm,tmp_n,tmp_np,chipar)
        
      end subroutine post__lagrangian

      end subroutine postProcessSol

c preProcessSol
c####################################################################
      subroutine preProcessSol(v1,v2)

c--------------------------------------------------------------------
c     Preprocess solution (place holder)
c--------------------------------------------------------------------

      use variables

      implicit none

c Call variables

      type(var_array),pointer :: v1,v2

c End program

      end subroutine preProcessSol

c tau_est
c ###################################################################
      function tau_est(g_def,vf) result(tol)

c -------------------------------------------------------------------
c     Computes Newton update tolerance based on tau error estimator
c -------------------------------------------------------------------

      use var_setup

      use ts_base

      use mg_xfer

      implicit none

c Call variables

      real(8) :: tol

      type(var_array),pointer :: vf

      type(grid_mg_def),pointer :: g_def

c Local variables

      type(var_array),pointer :: vc => null()

      integer  :: neq,ntotf,ntotc,igf,igc,nxf,nyf,nzf,nxc,nyc,nzc
      real(8),allocatable,dimension(:) :: res_f,res_c,res_c2

      INTERFACE
        subroutine evaluate_NLF(igr,varray,fi)
        use parameters
        use var_setup
        integer :: igr
        real(8) :: fi(ntotd)
        type(var_array),pointer :: varray
        end subroutine evaluate_NLF
      END INTERFACE

c Begin program

      igf = 1
      igc = igf + 1

      nxf = g_def%nxv(igf)
      nyf = g_def%nyv(igf)
      nzf = g_def%nzv(igf)

      nxc = g_def%nxv(igc)
      nyc = g_def%nyv(igc)
      nzc = g_def%nzv(igc)

      neq = vf%nvar

      ntotf = neq*nxf*nyf*nzf
      ntotc = neq*nxc*nyc*nzc

      allocate(res_f(ntotf),res_c(ntotc),res_c2(ntotc))

      call initializeDerivedType(vc,igr=igc)

c Evaluate nonlinear residuals

      !Restrict vf -> vc
      call restrictDerivedType(g_def,igf,vf,igc,vc,0)

      !Evaluate res_c
      call evaluate_NLF(igc,vc,res_c)

      !Evaluate res_f
      call evaluate_NLF(igf,vf,res_f)

      !Restrict res_f -> res_c2
      call restrictVector(g_def,neq,igf,res_f,igc,res_c2,0,.false.)

      !Compute residual difference norm
      tol = 0.1*dt*sqrt(sum((res_c - res_c2)**2))

c Destroy work storage

      deallocate(res_f,res_c,res_c2)
      call deallocateDerivedType(vc)

c End program

      end function tau_est
