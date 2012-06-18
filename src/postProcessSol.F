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

      implicit none

c Call variables

      type(var_array),pointer :: v_arr

c Local variables

      integer :: i,j,k,igx,igy,igz,nx,ny,nz,ieq,igrid
      real(8) :: mag,acc
      real(8),allocatable,dimension(:,:,:,:) :: db_cnv,da_cov
      logical :: E_diag,lag_init

c Interfaces

#if !defined(vec_pot) && !defined(samrai)
cc      INTERFACE
cc        subroutine applyBC (varray,iigx,iigy,iigz,obc)
cc        use variable_setup
cc        integer    :: iigx,iigy,iigz
cc        integer,optional :: obc
cc        type(var_array),pointer :: varray
cc        end subroutine applyBC
cc      END INTERFACE
#endif

c Begin program

      igrid = 1

      !Setup grid info
      igx = igrid
      igy = igrid
      igz = igrid

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Postprocess adiabatic equation for density update

cc      if (solve_prs .and. (.not.solve_rho)) then
cc        if (my_rank == 0) write (*,*) 'Postprocessing density'
cccc        v_arr%array_var(IRHO)%array(1:nx,1:ny,1:nz) =
cccc     .       u_0  %array_var(IRHO)%array(1:nx,1:ny,1:nz)
cccc     .     *(v_arr%array_var(ITMP)%array(1:nx,1:ny,1:nz)
cccc     .      /u_0  %array_var(ITMP)%array(1:nx,1:ny,1:nz))**(1./gamma)
cc        v_arr%array_var(IRHO)%array =
cc     .       u_0  %array_var(IRHO)%array
cc     .     *(v_arr%array_var(ITMP)%array
cc     .      /u_0  %array_var(ITMP)%array)**(1./gamma)
cc      endif

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
     .                   ,w_bcs=.true.)!,diag_plots=.true.)

cc        call curl_div_clean(db_cnv)

c       Find new B^(n+1)

        v_arr%array_var(IBX)%array = db_cnv(:,:,:,1)
     .                  +u_n%array_var(IBX)%array
        v_arr%array_var(IBY)%array = db_cnv(:,:,:,2)
     .                  +u_n%array_var(IBY)%array
        v_arr%array_var(IBZ)%array = db_cnv(:,:,:,3)
     .                  +u_n%array_var(IBZ)%array

c       Update helical flux

cc        psi = psi + da_cov(:,:,:,3)

        deallocate(db_cnv,da_cov)

      endif

#endif

c Perform Lagrangian temperature step

      if (lagrangian) then
        bcnv(:,:,:,1) = v_arr%array_var(IBX)%array
        bcnv(:,:,:,2) = v_arr%array_var(IBY)%array
        bcnv(:,:,:,3) = v_arr%array_var(IBZ)%array

        lag_init = .true.
        E_diag   = .true.
        call lagr_update(igrid,bcnv
     .                  ,u_nm %array_var(ITMP)%array
     .                  ,u_n  %array_var(ITMP)%array
     .                  ,v_arr%array_var(ITMP)%array
     .                  ,chi_par,chi,lag_init,E_diag,ilevel-3)

        !Update BCs
        call setMGBC(0,1,nx,ny,nz,igrid
     .              ,v_arr%array_var(ITMP)%array
     .              ,v_arr%array_var(ITMP)%bconds
     .              ,iorder=bc_order,icomp=(/ITMP/))
      endif

cc      call applyBC (v_arr,igx,igy,igz)

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
