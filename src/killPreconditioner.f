c killPreconditioner
c###################################################################
      subroutine killPreconditioner

      use precond_variables

      use newton_gmres

c     Deallocate variables

      if (precon == 'id') return

      deallocate (mgj0cnv)
      deallocate (mgdivV0)
      deallocate (mgadvdiffV0)
      deallocate (bcs)

      call deallocateMGArray(gp0)
      call deallocateMGArray(grho0)
      call deallocateMGArray(gb0)
      call deallocateMGArray(gv0)

      call deallocateDerivedType(varray)

      end subroutine killPreconditioner
