c killPreconditioner
c###################################################################
      subroutine killPreconditioner

      use precond_variables

c     Deallocate variables

      deallocate (bcnv,vcnv)
      deallocate (divrgV,divV)
      deallocate (bcs)
      deallocate (rho_diag,tmp_diag,b_diag,v_diag)
cc      deallocate (pp)
cc      deallocate (diag_mu)
cc      deallocate (idiagp,bcs)
cc      deallocate (ue)
cc      deallocate (d_sc,d_bc,d_pc,d_bh,d_alf)

      end subroutine killPreconditioner
