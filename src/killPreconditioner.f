c killPreconditioner
c###################################################################
      subroutine killPreconditioner

      use precond_variables

c     Deallocate variables
cc
cc      deallocate (bxx,byy,vxx,vyy,vxe,vye)
cc      deallocate (pp)
cc      deallocate (diag_mu)
cc      deallocate (idiagp,bcs)
cc      deallocate (ue)
cc      deallocate (d_sc,d_bc,d_pc,d_bh,d_alf)

      end subroutine killPreconditioner
