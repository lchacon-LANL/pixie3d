c postProcessSol
c####################################################################
      subroutine postProcessSol(v_arr)

c--------------------------------------------------------------------
c     Postprocess solution (e.g. filtering). In this instance, we
c     do divergence cleaning of the magnetic field.
c--------------------------------------------------------------------

      use nlfunction_setup

      use imposeBCinterface

      implicit none

c Call variables

      type(var_array),pointer :: v_arr

c Local variables

      integer :: i,j,k
      real(8),allocatable,dimension(:,:,:,:) :: b_cnv,a_cov

c Interfaces

      INTERFACE
        subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz)
        use variable_setup
        integer    :: iigx,iigy,iigz
        type(var_array),pointer :: varray
        end subroutine imposeBoundaryConditions
      END INTERFACE

c Begin program

#if !defined(vec_pot) 

      !Setup grid info
      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Apply boundary conditions

      call imposeBoundaryConditions(v_arr,igx,igy,igz)

cc      b_cnv => bcnv

      allocate(a_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,b_cnv(0:nx+1,0:ny+1,0:nz+1,3))

cc      b_cnv = bcnv
      b_cnv(:,:,:,1) = v_arr%array_var(IBX)%array(:,:,:)
      b_cnv(:,:,:,2) = v_arr%array_var(IBY)%array(:,:,:)
      b_cnv(:,:,:,3) = v_arr%array_var(IBZ)%array(:,:,:)

c Find dA* from dB* (w/ ghost cells)

      call curl_inv(nx,ny,nz,igx,igy,igz,b_cnv,a_cov)

cc#if !defined(petsc)
cc      if (isSP(1,1,1,igrid,igrid,igrid)) then
cc        !Impose BCs
cc        bcnd(:,1) = bcond
cc        bcnd(:,2) = bcond
cc        bcnd(:,3) = bcond
cc        where (bcnd /= SP) bcnd = DEF
cccc        where (bcnd == DEF) bcnd = EXT
cc
cc        call setMGBC(0,3,nx,ny,nz,igrid,da_cov,bcnd,icomp=IBX
cc     .              ,is_vec=.true.,is_cnv=.false.,iorder=2)
cccc        call setBC(IBX,3,nx,ny,nz,da_cov,da_cnv,vzeros,bcnd,igrid,igrid,igrid
cccc     .            ,is_cnv=.false.,iorder=2)
cc
cc      endif
cc#endif

c Find clean dB=curl(dA)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            b_cnv(i,j,k,:)=curl(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,a_cov(:,:,:,1)
     .                    ,a_cov(:,:,:,2)
     .                    ,a_cov(:,:,:,3))
            v_arr%array_var(IBX)%array(i,j,k) = b_cnv(i,j,k,1)
            v_arr%array_var(IBY)%array(i,j,k) = b_cnv(i,j,k,2)
            v_arr%array_var(IBZ)%array(i,j,k) = b_cnv(i,j,k,3)
          enddo
        enddo
      enddo

c End program

      deallocate(a_cov,b_cnv)

#endif

      end subroutine postProcessSol
