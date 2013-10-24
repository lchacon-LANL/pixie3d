c setup_app_PC
c####################################################################
      subroutine setup_app_PC(igrid,varray,nwt_it)

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

      integer :: igrid,nwt_it
      type (var_array),pointer :: varray

c Local variables

      integer :: igx,igy,igz

c Interfaces

c Begin program

      igx = igrid
      igy = igrid
      igz = igrid

      nwt_iter = nwt_it

c Setup PC parameters

      call PC_setup

c Impose BCs

#if !defined(samrai)
      call applyBC(igx,varray,bc_order=order_bc)
#else
      write (*,*) 'SAMRAI version of PC not ready'
      stop
#endif

c Find coefficients for linearized systems

      if (nwt_it == 1) call findPCcoeffs(igx,igy,igz,varray,gv%aux)

c Setup PC BCs

      call setup_PC_BC

c Find required diagonals in all grids

      if (nwt_it == 1) then
        form_diag = .true.
        call form_diagonals
        form_diag = .false.
      endif

c Initialize PC Jacobi it. counter

      iter_si = 0

c End program

      contains

c     setup_PC_BC
c     ##############################################################
      subroutine setup_PC_BC

c     --------------------------------------------------------------
c     Setup preconditioner BCs
c     --------------------------------------------------------------

        integer :: ieq

c     Set boundary conditions

        do ieq=1,neqd
          bcs(:,ieq) = varray%array_var(ieq)%bconds(:)
        enddo

c     Set BCs for auxiliary variables

#if defined(vec_pot)
        bcs(:,IBX:IBZ) = gv%aux%vec_list(IBCNV)%bconds
#endif

        bcs(:,IJX:IJZ) = gv%aux%vec_list(IJCNV)%bconds

c     Set parallel BCs

        call setASM_PC_BC(igx,igy,igz,bcs,pc_asm)

      end subroutine setup_PC_BC

c     form_diagonals
c     ##############################################################
      subroutine form_diagonals

c     --------------------------------------------------------------
c     Form require diagonals for smoothing in MG
c     --------------------------------------------------------------

c     Diag

        integer :: nxx,nyy,nzz,nt,icomp,neq,igrid,ii,jj,isig,ieq
        real(8),allocatable,dimension(:,:) :: v_mat

c     Externals

        external   v_mtvc,v_mtvc2,v_mtvc_af,b_mtvc,rho_mtvc,tmp_mtvc
     .            ,mstar_mtvc,hall_mtvc,a_mtvc,b_mtvc_hyp,a_mtvc_hyp
     .            ,v_mtvc_hyp

c     Begin program

c     Density

        call find_mf_diag(1,1,ntotdp,rho_mtvc,1,bcs(:,IRHO)
     .                   ,rho_diag,sm_ncolors,vol_wgt,mg_glrkin)

c     Temperature

        call find_mf_diag(1,1,ntotdp,tmp_mtvc,1,bcs(:,ITMP)
     .                   ,tmp_diag,sm_ncolors,vol_wgt,mg_glrkin)

        if (sm_type /= 'gm') then

c       Magnetic field

#if !defined(vec_pot)
          if (di > 0d0 .and. (.not.pc_ion_hall)) then
            bcs_hyp(:,1:3) = bcs(:,IBX:IBZ)
            bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
            call find_mf_diag(6,1,6*ntotdp,b_mtvc_hyp,1,bcs_hyp
     .                       ,b_diag_hyp,sm_ncolors,vol_wgt,mg_glrkin)
          else
            call find_mf_diag(3,1,3*ntotdp,b_mtvc,1,bcs(:,IBX:IBZ)
     .                       ,b_diag,sm_ncolors,vol_wgt,mg_glrkin)
          endif
#else
          if (di > 0d0 .and. (.not.pc_ion_hall)) then
            bcs_hyp(:,1:3) = bcs(:,IAX:IAZ)
            if (fake_ve) then
              bcs_hyp(:,4:6) = gv%aux%vec_list(IVEFCNV)%bconds
            else
              bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
            endif
            call find_mf_diag(6,1,6*ntotdp,a_mtvc_hyp,1,bcs_hyp
     .                       ,a_diag_hyp,sm_ncolors,vol_wgt,mg_glrkin)
          endif

          call find_mf_diag(3,1,3*ntotdp,a_mtvc,1,bcs(:,IAX:IAZ)
     .                       ,a_diag,sm_ncolors,vol_wgt,mg_glrkin)
#endif

c        Velocity

          if (pc_v_solve) then
            call find_mf_diag(3,1,3*ntotdp,v_mtvc2,1,bcs(:,IVX:IVZ)
     .                       ,v_diag2,sm_ncolors,vol_wgt,mg_glrkin)
          endif

          if (     pc_type == 's1'
     .        .or. pc_type == 's2'
     .        .or. pc_type == 's4') then
            if (di == 0d0 .or. (.not.pc_ion_hall)) then
              call find_mf_diag(3,1,3*ntotdp,v_mtvc,1,bcs(:,IVX:IVZ)
     .                         ,v_diag,sm_ncolors,vol_wgt,mg_glrkin)
            else
              bcs_hyp(:,1:3) = bcs(:,IVX:IVZ)
              if (pc_ion_hall.and.fake_ve) then
                bcs_hyp(:,4:6) = gv%aux%vec_list(IVEFCNV)%bconds
              else
                bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
              endif
              call find_mf_diag(6,1,6*ntotdp,v_mtvc_hyp,1,bcs_hyp
     .                         ,v_diag_hyp,sm_ncolors,vol_wgt,mg_glrkin)
            endif
          elseif (pc_type == 's3') then
            call find_mf_diag(3,1,3*ntotdp,v_mtvc_af,1,bcs(:,IVX:IVZ)
     .                       ,v_diag,sm_ncolors,vol_wgt,mg_glrkin)
          endif
        endif

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
cccc      if (debug) then
cc      write (*,*) 'DIAG -- setupPC: Building matrix...'
cc
cc      neq = 6
cc      igrid = 1
cc
cc      nxx = grid_params%nxv(igrid)
cc      nyy = grid_params%nyv(igrid)
cc      nzz = grid_params%nzv(igrid)
cc
cc      nt = neq*nxx*nyy*nzz
cc
cc      allocate(v_mat(nt,nt))
cc
cccc      icomp = IVX
cccc      call find_mf_mat(neq,nt,v_mtvc,igrid,bcs(:,IVX:IVZ),v_mat)
cc      icomp = IAX
cc      bcs_hyp(:,1:3) = bcs(:,IAX:IAZ)
cc      bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
cc      call find_mf_mat(neq,nt,a_mtvc_hyp,igrid,bcs_hyp,v_mat)
cc
cccc      open(unit=110,file='debug.bin',form='unformatted'
cccc     .       ,status='replace')
cccc      call contour(v_mat,nt,nt,0d0,1d0,0d0,1d0,0,110)
cccc      close(110)
cc
cc      open(unit=110,file='debug.mat',status='replace')
cc      do ii=1,nt
cc        write (110,*) (v_mat(ii,jj),jj=1,nt)
cc      enddo
cc      close(110)
cc
cc      deallocate(v_mat)
cc
cc      write (*,*) 'Done!'
cc
cc      stop
cccc      endif
c diag *****

      end subroutine form_diagonals

      end subroutine setup_app_PC

c kill_app_PC
c###################################################################
      subroutine kill_app_PC

      use precond_variables

      end subroutine kill_app_PC
