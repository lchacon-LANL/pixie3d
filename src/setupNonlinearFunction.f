c setupNonlinearFunction
c#################################################################
      subroutine setupNonlinearFunction(varray)
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use parameters

      use variables

      use nlfunction_setup

      use equilibrium

      use constants

      use timeStepping

      implicit none

c Call variables

      type (var_array),target :: varray

c Local variables

      integer(4) :: i,j,k,ig,jg,kg

      integer(4) :: dim,loc,ibc,ieq,bctype

      real(8)    :: dh,vmax

c Externals

      real*8       vis,res
      external     vis,res

c Begin program

c Impose boundary conditions and find auxiliary quantities

      call imposeBoundaryConditions(varray)

c Safeguards

      if (minval(rho) < 0d0) then
        write (*,*) 'Warning: negative densities are occurring'
        write (*,*) rho(:,:,:)
        stop
      endif

      if (minval(tmp) < 0d0) then
        write (*,*) 'Warning: negative temperatures are occurring'
      endif

c Transport parameters

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            nuu  (i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
            eeta (i,j,k) = res(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

c End program

      end subroutine setupNonlinearFunction
