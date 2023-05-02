c setup_app_PC
c####################################################################
      subroutine setup_app_PC(igrid,varray)

c--------------------------------------------------------------------
c     Finds velocity and magnetic field components in all grids
c     for matrix-free preconditioning.
c--------------------------------------------------------------------

      use grid

      use PC_var

      use matvec

      use operators

      use imposeBCinterface

      implicit none

c Call variables

      integer :: igrid
      type (var_array),pointer :: varray

c Local variables

      integer :: igx,igy,igz

c Interfaces

c Begin program

      igx = igrid
      igy = igrid
      igz = igrid

c Configure PC parameters

      call PC_config

c Allocate PC variables

      if (.not.pc_vars_allcd) call allocPCVariables

c Impose BCs

#if !defined(samrai)
      nwt_iter = nk_current_Newton_iteration()

      call applyBC(igx,varray,gv%aux,bc_order=order_bc)
#else
      write (*,*) 'SAMRAI version of PC not ready'
      stop
#endif

c Find coefficients for linearized systems

      if (nwt_iter == 1) call findPCcoeffs(igrid,varray,gv%aux)

c Find required diagonals in all grids

      if (nwt_iter == 1) call form_diagonals

c Initialize PC Jacobi it. counter

      iter_si = 0

c End program

      contains

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
     .            ,v_mtvc_hyp,prs_mtvc

c     Begin program

c     Density

        call find_mf_diag(gv%gparams,1,1,ntotdp,rho_mtvc,1,bcs(:,IRHO)
     .                   ,rho_diag,sm_ncolors,vol_wgt,mg_glrkin)

c     Temperature

        if (solve_prs) then
         call find_mf_diag(gv%gparams,1,1,ntotdp,prs_mtvc,1,bcs(:,IPRS)
     .                    ,prs_diag,sm_ncolors,vol_wgt,mg_glrkin)
        else
         call find_mf_diag(gv%gparams,1,1,ntotdp,tmp_mtvc,1,bcs(:,ITMPE)
     .                    ,tmp_diag,sm_ncolors,vol_wgt,mg_glrkin)
        endif
        
        if (sm_type == 'gm') return

c      Magnetic field

#if !defined(vec_pot)
        if (pc_e_hall) then
          bcs_hyp(:,1:3) = bcs(:,IBX:IBZ)
          bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
          call find_mf_diag(gv%gparams,6,1,6*ntotdp,b_mtvc_hyp,1
     .        ,bcs_hyp,b_diag_hyp,sm_ncolors,vol_wgt,mg_glrkin)
        else
          call find_mf_diag(gv%gparams,3,1,3*ntotdp,b_mtvc,1,bcs(:
     .        ,IBX:IBZ),b_diag,sm_ncolors,vol_wgt,mg_glrkin)
        endif
#else
        if (pc_e_hall) then
          bcs_hyp(:,1:3) = bcs(:,IAX:IAZ)
          bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
          call find_mf_diag(gv%gparams,6,1,6*ntotdp,a_mtvc_hyp,1
     .        ,bcs_hyp,a_diag_hyp,sm_ncolors,vol_wgt,mg_glrkin)
        endif

        call find_mf_diag(gv%gparams,3,1,3*ntotdp,a_mtvc,1,bcs(:
     .      ,IAX:IAZ),a_diag,sm_ncolors,vol_wgt,mg_glrkin)
#endif

c       Velocity

        if (pc_v_solve) then
          call find_mf_diag(gv%gparams,3,1,3*ntotdp,v_mtvc2,1,bcs(:
     .         ,IVX:IVZ),v_diag2,sm_ncolors,vol_wgt,mg_glrkin)
        endif

        select case(pc_type)
        case('s1','s2','s4')
          if (di == 0d0 .or. (.not.pc_ion_hall)) then
            call find_mf_diag(gv%gparams,3,1,3*ntotdp,v_mtvc,1,bcs(:
     .           ,IVX:IVZ),v_diag,sm_ncolors,vol_wgt,mg_glrkin)
          else
            bcs_hyp(:,1:3) = bcs(:,IVX:IVZ)
            bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
            call find_mf_diag(gv%gparams,6,1,6*ntotdp,v_mtvc_hyp,1
     .           ,bcs_hyp,v_diag_hyp,sm_ncolors,vol_wgt,mg_glrkin)
          endif
        case('s3')
          call find_mf_diag(gv%gparams,3,1,3*ntotdp,v_mtvc_af,1,bcs(:
     .         ,IVX:IVZ),v_diag,sm_ncolors,vol_wgt,mg_glrkin)
        end select

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
c$$$cc      if (debug) then
c$$$      write (*,*) 'DIAG -- setupPC: Building matrix...'
c$$$
c$$$      neq = 6
c$$$      igrid = 1
c$$$
c$$$      nxx = gv%gparams%nxv(igrid)
c$$$      nyy = gv%gparams%nyv(igrid)
c$$$      nzz = gv%gparams%nzv(igrid)
c$$$
c$$$      nt = neq*nxx*nyy*nzz
c$$$
c$$$      allocate(v_mat(nt,nt))
c$$$
c$$$cc      icomp = IVX
c$$$cc      call find_mf_mat(neq,nt,v_mtvc,igrid,bcs(:,IVX:IVZ),v_mat)
c$$$cc      icomp = IAX
c$$$cc      bcs_hyp(:,1:3) = bcs(:,IAX:IAZ)
c$$$cc      bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
c$$$cc      call find_mf_mat(neq,nt,a_mtvc_hyp,igrid,bcs_hyp,v_mat)
c$$$      icomp = IVX
c$$$      bcs_hyp(:,1:3) = bcs(:,IVX:IVZ)
c$$$      bcs_hyp(:,4:6) = bcs(:,IJX:IJZ)
c$$$      call find_mf_mat(gv%gparams,neq,nt,v_mtvc_hyp,igrid,bcs_hyp,v_mat)
c$$$
c$$$cc      open(unit=110,file='debug.bin',form='unformatted'
c$$$cc     .       ,status='replace')
c$$$cc      call contour(v_mat,nt,nt,0d0,1d0,0d0,1d0,0,110)
c$$$cc      close(110)
c$$$
c$$$      open(unit=110,file='debug.mat',status='replace')
c$$$      do ii=1,nt
c$$$        write (110,*) (v_mat(ii,jj),jj=1,nt)
c$$$      enddo
c$$$      close(110)
c$$$
c$$$      deallocate(v_mat)
c$$$
c$$$      write (*,*) 'Done!'
c$$$
c$$$      stop
c$$$cc      endif
c diag *****

      end subroutine form_diagonals

      end subroutine setup_app_PC

c kill_app_PC
c###################################################################
      subroutine kill_app_PC

      use PC_var

c     Deallocate LOCAL auxiliary variables

      call destroyAppBCs

cc      if (pc_vars_allcd) call deallocPCVariables
      
      end subroutine kill_app_PC
