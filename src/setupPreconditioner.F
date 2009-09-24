c setupPreconditioner
c####################################################################
      subroutine setupPreconditioner(v_arr)

c--------------------------------------------------------------------
c     Finds velocity and magnetic field components in all grids
c     for matrix-free preconditioning.
c--------------------------------------------------------------------

      use grid

      use precond_variables

      use matvec

      use operators

      use imposeBCinterface

      implicit none

c Call variables

      type(var_array),pointer :: v_arr

c Local variables

      integer :: igx,igy,igz,igrid,ii,jj,isig,nx,ny,nz

c Diag

      integer :: nxx,nyy,nzz,nt,icomp
      real(8),allocatable,dimension(:,:) :: v_mat

c Externals

      external   v_mtvc,v_mtvc2,v_mtvc3,b_mtvc,rho_mtvc,tmp_mtvc
     .          ,mstar_mtvc,hall_mtvc,a_mtvc,b_mtvc_hyp,a_mtvc_hyp

c Interfaces

      INTERFACE
        subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz,obc)
        use variable_setup
        integer    :: iigx,iigy,iigz
        integer,optional :: obc
        type(var_array),pointer :: varray
        end subroutine imposeBoundaryConditions
      END INTERFACE

c Begin program

      if (precon == 'id') return

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Setup PC

      call PC_setup

c Unpack vector x

      call equateDerivedType(varray,v_arr)

      call imposeBoundaryConditions(varray,igx,igy,igz)

c Extract arrays and BC's

      rho => varray%array_var(IRHO)%array
      px  => varray%array_var(IVX )%array
      py  => varray%array_var(IVY )%array
      pz  => varray%array_var(IVZ )%array
#if !defined(vec_pot)
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
#else
      ax  => varray%array_var(IAX )%array
      ay  => varray%array_var(IAY )%array
      az  => varray%array_var(IAZ )%array
#endif
      tmp => varray%array_var(ITMP)%array

c Find coefficients for linearized systems

      if (jit == 1) call findPCcoeffs(igx,igy,igz)

c Setup BCs for diagonal formation

      !Set boundary conditions
      do ieq=1,neqd
        bcs(:,ieq) = varray%array_var(ieq)%bconds(:)
      enddo

cc#if !defined(vec_pot)
cc      bcs(:,IJX:IJZ) = bcs(:,IBX:IBZ)
cc      where (bcs(:,IJX:IJZ) == -NEU) bcs(:,IJX:IJZ) = -DIR  !Use covariant components 
cc#else
#if defined(vec_pot)
      where (bcs(:,IAX:IAZ) == -EQU) bcs(:,IAX:IAZ) = -DIR  !Use covariant components 

      !Magnetic field conditions
      bcs(:,IBX:IBZ) = bcs(:,IAX:IAZ)
      where (bcs(:,IBX:IBZ) == -DIR) bcs(:,IBX:IBZ) = -NEU  !Use covariant components 
#endif

      bcs(:,IJX:IJZ) = bcs(:,IBX:IBZ)
      where (bcs(:,IJX:IJZ) == -NEU) bcs(:,IJX:IJZ) = -DIR  !Use covariant components 

      if (asm_PC.and.(np>1)) then
        call setASM_PC_BC(igx,igy,igz,bcs)
      else
        call setASMflag(.false.)
      endif

#if defined(petsc)
      if (ilevel >= 5) then
        write (*,*) 'Proc ',my_rank,': asm     in setupPC',asm
        write (*,*) 'Proc ',my_rank,': par_bc  in setupPC',par_bc
        write (*,*) 'Proc ',my_rank,': asm_dir in setupPC',asm_dir
      endif
#endif

c Find required diagonals in all grids

      if (jit == 1) then

        call find_mf_diag(1,1,ntotdp,rho_mtvc,1,bcs(:,IRHO)
     .                   ,rho_diag,ncolors,vol_wgt,mg_glrkin)

        call find_mf_diag(1,1,ntotdp,tmp_mtvc,1,bcs(:,ITMP)
     .                   ,tmp_diag,ncolors,vol_wgt,mg_glrkin)

        if (.not.gm_smooth) then
#if !defined(vec_pot)
          if (di > 0d0 .and. (.not.hall_new)) then
            bcs_hyp(:,1:3) = bcs(:,IBX:IBZ)
            bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
            call find_mf_diag(6,1,6*ntotdp,b_mtvc_hyp,1,bcs_hyp
     .                       ,b_diag_hyp,ncolors,vol_wgt,mg_glrkin)
          else
            call find_mf_diag(3,1,3*ntotdp,b_mtvc,1,bcs(:,IBX:IBZ)
     .                       ,b_diag,ncolors,vol_wgt,mg_glrkin)
          endif
#else
          if (di > 0d0 .and. (.not.hall_new)) then
            bcs_hyp(:,1:3) = bcs(:,IAX:IAZ)
            bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
            call find_mf_diag(6,1,6*ntotdp,a_mtvc_hyp,1,bcs_hyp
     .                       ,a_diag_hyp,ncolors,vol_wgt,mg_glrkin)
          else
            call find_mf_diag(3,1,3*ntotdp,a_mtvc,1,bcs(:,IAX:IAZ)
     .                       ,a_diag,ncolors,vol_wgt,mg_glrkin)
          endif
#endif
          if (v_solve) then
            call find_mf_diag(3,1,3*ntotdp,v_mtvc2,1,bcs(:,IVX:IVZ)
     .                       ,v_diag2,ncolors,vol_wgt,mg_glrkin)
          endif

          if (precon == 's1' .or. precon == 's2') then
            call find_mf_diag(3,1,3*ntotdp,v_mtvc,1,bcs(:,IVX:IVZ)
     .                       ,v_diag,ncolors,vol_wgt,mg_glrkin)
          elseif (precon == 's3') then
            call find_mf_diag(3,1,3*ntotdp,v_mtvc3,1,bcs(:,IVX:IVZ)
     .                       ,v_diag,ncolors,vol_wgt,mg_glrkin)
          endif
        endif

      endif

c Initialize PC Jacobi it. counter

      si_it = 0

c diag ***** DIAGONAL WRITE UP
cc      isig = 1
cc      do igrid = 1,ngrid
cc
cc        write (*,*) 'Grid level',igrid
cc        write (*,*) 'Grid level pointer:',isig
cc
cc        nxx = grid_params%nxv(igrid)
cc        nyy = grid_params%nyv(igrid)
cc        nzz = grid_params%nzv(igrid)
cc
cc        nt = nxx*nyy*nzz
cc
cc        do k=1,nzz
cc          do j=1,nyy
cc            do i=1,nxx
cc        
cc              ii = i + nxx*(j-1) + nxx*nyy*(k-1)
cc
cc              jj  = (ii-1)*3 + isig - 1
cc
cc              write (*,*) 'Grid pos:',i,j,k,jj+1
cc              write (*,*) v_diag(:,jj+1)
cc              write (*,*) v_diag(:,jj+2)
cc              write (*,*) v_diag(:,jj+3)
cc              write (*,*)
cc            enddo
cc          enddo
cc        enddo
cc
cc        isig = isig + 3*nt
cc      enddo
cc      stop
c diag *****

c diag ***** DIAGONAL GENERATION FOR XDRAW PLOTTING
cc      do k=1,nzz
cc        do j=1,nyy
cc          do i=1,nxx
cc            ii = i + nxx*(j-1)+nxx*nyy*(k-1)
cc            debug(i,j,k) = 1./sqrt(sum(v_diag(:,3*ii-2:3*ii)**2)
cc          enddo
cc        enddo
cc      enddo
cc
cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .    ,status='replace')
cc      call contour(debug(1:nxx,1:nyy,1),nxx,nyy,0d0,xmax,0d0,ymax,0,110)
cc      close(110)
cc      stop
c diag *****

c diag ***** MATRIX GENERATION FOR MATLAB PROCESSING
cc      if (debug) then
cc      write (*,*) 'DIAG -- setupPC: Building matrix...'
cc
cc      igrid = 1
cc
cc      nxx = grid_params%nxv(igrid)
cc      nyy = grid_params%nyv(igrid)
cc      nzz = grid_params%nzv(igrid)
cc
cc      nt = 3*nxx*nyy*nzz
cc
cc      allocate(v_mat(nt,nt))
cc
cc      icomp = IVX
cc      call find_mf_mat(3,nt,v_mtvc,igrid,bcs(:,IVX:IVZ),v_mat)
cc
cccc      open(unit=110,file='debug.bin',form='unformatted'
cccc     .       ,status='replace')
cccc      call contour(v_mat,nt,nt,0d0,1d0,0d0,1d0,0,110)
cccc      close(110)
cc
cc      open(unit=110,file='debug.mat',status='replace')
cc      do i=1,nt
cc        write (110,*) (v_mat(i,j),j=1,nt)
cc      enddo
cc      close(110)
cc
cc      deallocate(v_mat)
cc
cc      write (*,*) 'Done!'
cc
cccc      stop
cc      endif
c diag *****

c End program

      end subroutine setupPreconditioner

c killPreconditioner
c###################################################################
      subroutine killPreconditioner

      use precond_variables

c     Deallocate variables

#if defined(samrai)
      call deallocAuxVariables
#endif

      end subroutine killPreconditioner
