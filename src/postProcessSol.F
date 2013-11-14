c postProcessSol
c####################################################################
      subroutine postProcessSol(v_arr)

c--------------------------------------------------------------------
c     Postprocess solution (e.g. filtering). In this instance, we
c     do divergence cleaning of the magnetic field.
c--------------------------------------------------------------------

      use nlfunction_setup

      use imposeBCinterface

      use precond_variables, only: pc_B_solve,pc_type,pc_divclean

      use xdraw_io

      use B_tools

      use par_int

      use precond_setup, only: pc_as_solver

      implicit none

c Call variables

      type(var_array),pointer :: v_arr

c Local variables

      integer :: i,j,k,igx,igy,igz,nx,ny,nz,ieq,igrid
      real(8) :: mag,acc
      real(8),allocatable,dimension(:,:,:,:) :: db_cnv,da_cov

      !Lagrangian
      real(8),allocatable,dimension(:,:,:)   :: tmp_nm,tmp_n,tmp_np
     $                                         ,chipar,array
      logical :: B_weight
      real(8) :: G_tol
      integer :: sp_order

      !Parareal
      integer :: upit=4321
      character(1024) :: pit_file
      real(8) :: Et,Ek,Em

c Interfaces

c Begin program

      igrid = 1

      !Setup grid info
      igx = igrid
      igy = igrid
      igz = igrid

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Add random noise in EOM (to break symmetry)

      if (noise) then
        if (nc_eom_v) then
          do k=1,nz
            do j=1,ny
              do i=1,nx
                do ieq=IVX,IVZ
                  call random_number(acc) !Between 0 and 1
                  acc = noise_lev*sqrt(dt)*(1-2.*acc)
                  v_arr%array_var(ieq)%array=v_arr%array_var(ieq)%array
     .                                      +acc
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
                  v_arr%array_var(ieq)%array=v_arr%array_var(ieq )%array
     .                                  +acc*v_arr%array_var(IRHO)%array
                enddo
              enddo
            enddo
          enddo
        endif
      endif

c Divergence-clean magnetic field when PC is on
#if !defined(samrai)
#if !defined(vec_pot)

      if (.not.((.not.pc_B_solve)
     .    .or.   pc_divclean
     .    .or.  (.not.solenoidal)
     .    .or.  (pc_type == 'id')
     .    .or.  (nxd==1.and.nyd==1)
     .    .or.  (nxd==1.and.nzd==1)
     .    .or.  (nyd==1.and.nzd==1)) ) then

        if (ilevel > 0 .and. my_rank == 0) then
          write (*,*)
          write (*,*) 'Performing divergence cleaning in postproc...'
        endif

c       Find db=B^(n+1)-B^n

        allocate(db_cnv(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,da_cov(0:nx+1,0:ny+1,0:nz+1,3))

        db_cnv(:,:,:,1) = v_arr%array_var(IBX)%array
     .                   -u_n  %array_var(IBX)%array
        db_cnv(:,:,:,2) = v_arr%array_var(IBY)%array
     .                   -u_n  %array_var(IBY)%array
        db_cnv(:,:,:,3) = v_arr%array_var(IBZ)%array
     .                   -u_n  %array_var(IBZ)%array

c       Clean db

cc        call grad_div_clean(db_cnv)

        call vlap_div_clean(ilevel,nx,ny,nz
     .                   ,igx,igy,igz,db_cnv,dacov=da_cov
     .                   ,w_bcs=.true.,gm_it=it_divcl)!,diag_plots=.true.)

cc        call curl_div_clean(nx,ny,nz,igx,igy,igz,db_cnv)

c       Find new B^(n+1)

        v_arr%array_var(IBX)%array = db_cnv(:,:,:,1)
     .                  +u_n%array_var(IBX)%array
        v_arr%array_var(IBY)%array = db_cnv(:,:,:,2)
     .                  +u_n%array_var(IBY)%array
        v_arr%array_var(IBZ)%array = db_cnv(:,:,:,3)
     .                  +u_n%array_var(IBZ)%array

c       Update helical flux and magnetic field

cc        psi = psi + da_cov(:,:,:,3)
        bcnv(:,:,:,1) = v_arr%array_var(IBX)%array
        bcnv(:,:,:,2) = v_arr%array_var(IBY)%array
        bcnv(:,:,:,3) = v_arr%array_var(IBZ)%array

        deallocate(db_cnv,da_cov)

      endif

#endif

c Perform Lagrangian temperature step

      if (lagrangian) then

        allocate(tmp_nm(0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_n (0:nx+1,0:ny+1,0:nz+1)
     $          ,tmp_np(0:nx+1,0:ny+1,0:nz+1)
     $          ,chipar(0:nx+1,0:ny+1,0:nz+1))

        chipar = chi_par - chi
        G_tol = 1d-8
        sp_order = 2
        B_weight = .false.

        if (maxval(chipar) == 0d0) return

        if (ilevel > 0 .and. my_rank == 0)
     .     write (*,*) '  --> Performing parallel transport step...'

        call lag_init(igrid,nx,ny,nz,bcnv,2d0,ilevel-3
     .               ,G_tol,sp_order,B_weight,.false.)

        !Find temp from pressure (automatically takes care of jacobian factor)
        tmp_nm = u_nm %array_var(ITMP)%array
     $          /u_nm %array_var(IRHO)%array
        tmp_n  = u_n  %array_var(ITMP)%array
     $          /u_n  %array_var(IRHO)%array
        tmp_np = v_arr%array_var(ITMP)%array
     $          /v_arr%array_var(IRHO)%array

        call lag_update_bdf2(igrid,ITMP,v_arr%array_var(ITMP)%bconds
     $                      ,tmp_nm,tmp_n,tmp_np,chipar,ilevel-3)
c$$$        lag_init = .true.
c$$$        E_diag   = .true.
c$$$        call lagr_update(igrid,bcnv
c$$$     .                  ,u_nm %array_var(ITMP)%array
c$$$     .                  ,u_n  %array_var(ITMP)%array
c$$$     .                  ,v_arr%array_var(ITMP)%array
c$$$     .                  ,chi_par,chi,lag_init,E_diag,ilevel-3)

        call lag_destroy

        !Update BCs
        call setMGBC(0,1,nx,ny,nz,igrid,tmp_np
     .              ,v_arr%array_var(ITMP)%bconds
     .              ,iorder=order_bc,icomp=(/ITMP/))

        !Return pressure (with jacobian pressure)
        v_arr%array_var(ITMP)%array = tmp_np*v_arr%array_var(IRHO)%array

        deallocate(tmp_nm,tmp_n,tmp_np,chipar)
      endif

cc      call applyBC (v_arr,igx,igy,igz)

c     Dump parareal convergence info

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
     .                      ,dt,v_arr)
      endif

#endif

c End program

#endif

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
