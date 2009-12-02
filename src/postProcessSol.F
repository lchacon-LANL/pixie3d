c postProcessSol
c####################################################################
      subroutine postProcessSol(v_arr)

c--------------------------------------------------------------------
c     Postprocess solution (e.g. filtering). In this instance, we
c     do divergence cleaning of the magnetic field.
c--------------------------------------------------------------------

      use nlfunction_setup

      use imposeBCinterface

      use precond_variables, only: inv_B,divCln,precon

      implicit none

c Call variables

      type(var_array),pointer :: v_arr

c Local variables

      integer :: i,j,k,igx,igy,igz,nx,ny,nz
      real(8) :: mag
      real(8),allocatable,dimension(:,:,:,:) :: b_cnv,a_cov,dbg

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

      !Setup grid info
      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      if (adiabatic .and. (.not.solve_rho)) then
        if (my_rank == 0) write (*,*) 'Postprocessing density'
        v_arr%array_var(IRHO)%array(1:nx,1:ny,1:nz) =
     .       u_0  %array_var(IRHO)%array(1:nx,1:ny,1:nz)
     .     *(v_arr%array_var(ITMP)%array(1:nx,1:ny,1:nz)
     .      /u_0  %array_var(ITMP)%array(1:nx,1:ny,1:nz))**(1./gamma)
      endif

c Save previous time step information for next time step (for BCs)

cc      bcnv_n = bcnv

c Divergence-clean magnetic field when PC is on

      allocate(dbg(nx,ny,nz,1))
      do k=1,nz
        do j=1,ny
          do i=1,nx
            dbg(i,j,k,1) = div(i,j,k,nx,ny,nz,igx,igy,igz
     $                        ,bcnv(:,:,:,1)
     $                        ,bcnv(:,:,:,2)
     $                        ,bcnv(:,:,:,3))
          enddo
        enddo
      enddo

      mag = sqrt(dot(ntotdp,dbg,dbg))
      if (my_rank == 0) write (*,*) 'L-2 norm of div(dB)',mag

      deallocate(dbg)

#if !defined(vec_pot) 

      if ((.not.(inv_B.and.divCln)).or.(precon == 'id')) return

cc      if (ilevel > 0 .and. my_rank == 0) then
      if (my_rank == 0) then
        write (*,*)
        write (*,*) 'Performing divergence cleaning...'
      endif

c Apply boundary conditions

      call imposeBoundaryConditions(v_arr,igx,igy,igz)

c Find perturbations in magnetic field, dB

      allocate(a_cov(0:nx+1,0:ny+1,0:nz+1,3)
     .        ,b_cnv(0:nx+1,0:ny+1,0:nz+1,3))

      b_cnv(:,:,:,1) = v_arr%array_var(IBX)%array
     .                -u_n  %array_var(IBX)%array
      b_cnv(:,:,:,2) = v_arr%array_var(IBY)%array
     .                -u_n  %array_var(IBY)%array
      b_cnv(:,:,:,3) = v_arr%array_var(IBZ)%array
     .                -u_n  %array_var(IBZ)%array

c Find dA from dB (w/ ghost cells)

      call curl_inv(nx,ny,nz,igx,igy,igz,b_cnv,a_cov)

c Find clean B= B_n + curl(dA)

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx
            b_cnv(i,j,k,:)=curl(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,a_cov(:,:,:,1)
     .                    ,a_cov(:,:,:,2)
     .                    ,a_cov(:,:,:,3))
            v_arr%array_var(IBX)%array(i,j,k) = b_cnv(i,j,k,1)
     .                +u_n  %array_var(IBX)%array(i,j,k)
            v_arr%array_var(IBY)%array(i,j,k) = b_cnv(i,j,k,2)
     .                +u_n  %array_var(IBY)%array(i,j,k)
            v_arr%array_var(IBZ)%array(i,j,k) = b_cnv(i,j,k,3)
     .                +u_n  %array_var(IBZ)%array(i,j,k)
          enddo
        enddo
      enddo

c End program

      deallocate(a_cov,b_cnv)

#endif

      end subroutine postProcessSol
