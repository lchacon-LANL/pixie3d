c fixDivB
c ######################################################################
      subroutine fixDivB(patch_var)
      use variable_setup
      use problem_def
      use variables
      use variable_setup
      implicit none

      type(patch), TARGET :: patch_var
      integer :: nbcgrps
      real(8),pointer :: B_cnv(:,:,:,:)
      real(8),pointer :: B_cov(:,:,:,:)

c Get some basic properties

      gv => patch_var
      grid_params => gv%gparams
      gmetric => gv%gparams%gmetric
      call setVectorDimensions

      B_cnv => patch_var%aux%vec_list(IBCNV)%vec
      B_cov => patch_var%aux%vec_list(IBCOV)%vec

      nxbc = grid_params%nxv(1)
      nybc = grid_params%nyv(1)
      nzbc = grid_params%nzv(1)
      iimnm = grid_params%ilo(1)
      iimxp = grid_params%ihi(1)
      jjmnm = grid_params%jlo(1)
      jjmxp = grid_params%jhi(1)
      kkmnm = grid_params%klo(1)
      kkmxp = grid_params%khi(1)
      igxbc = 1
      igybc = 1
      igzbc = 1

c Correct B cnv
      call solen_bc(B_cnv,1,0,.false.)
      call solen_bc(B_cnv,2,0,.false.)
      call solen_bc(B_cnv,3,0,.false.)
      call solen_bc(B_cnv,1,1,.false.)
      call solen_bc(B_cnv,2,1,.false.)
      call solen_bc(B_cnv,3,1,.false.)

c Correct B cov
      call XformVector_BC(igxbc,B_cnv,B_cov,.false.)

      end subroutine fixDivB


