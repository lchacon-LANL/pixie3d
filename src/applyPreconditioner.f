c applyPreconditioner
c###################################################################
      subroutine applyPreconditioner(ntot,y,x,iout)

c-------------------------------------------------------------------
c     This subroutine solves P x = y for the vector x. 
c     The parameter iout controls output level:
c       * iout <= 0 --> No output
c       * iout >  0 --> Level of output design by user.
c-------------------------------------------------------------------

      use parameters

      use grid

      use precond_variables

      use iosetup

      implicit none

c Call variables

      integer ::  ntot,iout
      real*8      x(ntot),y(ntot)

c Local variables

      real*8      xx(ntotdp,neqd)
      real*8      yy(ntotdp,neqd)

      integer*4   i,j,k,ii,iii,ieq

c Externals

c Begin program

c *******************************************************************
c     Identity preconditioner
c *******************************************************************

      if (precon.eq.'id') then

        x = y*dt   !Diagonal scaling

      else

c     Set up size of operator split vectors

        do ieq=1,neqd
          do k = 1,nzd
            do j = 1,nyd
              do i = 1,nxd
                ii  = i + nxd*(j-1) + nxd*nyd*(k-1)
                iii = ieq + neqd*(i-1 + nxd*(j-1) + nxd*nyd*(k-1))

                yy(ii,ieq) = y(iii)
                xx(ii,ieq) = 0d0
              enddo
            enddo
          enddo
        enddo

c     Process vectors

c     Store solution in "x"

        do k = 1,nzd
          do j = 1,nyd
            do i = 1,nxd
              do ieq=1,neqd
                ii  = i + nxd*(j-1) + nxd*nyd*(k-1)
                iii = ieq + neqd*(i-1 + nxd*(j-1) + nxd*nyd*(k-1))
                x(iii) = xx(ii,ieq)
              enddo
            enddo
          enddo
        enddo

      endif

c End program

      end subroutine applyPreconditioner
