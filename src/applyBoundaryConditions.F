c defineBoundaryConditions
c####################################################################
      subroutine defineBoundaryConditions (neq,bbcs)
c--------------------------------------------------------------------
c     Defines boundary conditions of physical quantities.
c     On input:
c       * neq -> number of equations
c     On output:
c       * bbcs -> real array of size (6,neq) containing BC setup:
c           + bbcs(1) ---> at x0
c           + bbcs(2) ---> at x1
c           + bbcs(3) ---> at y0
c           + bbcs(4) ---> at y1
c           + bbcs(5) ---> at z0
c           + bbcs(6) ---> at z1
c     Definition of BC identifiers is given in "grid_mod.f". In vectors,
c     a negative BC identifier means that BCs are to be imposed on
c     covariant components instead of on (default) contravariant comps.
c--------------------------------------------------------------------

      use icond

      use grid

      use equilibrium

      use variables

      use auxiliaryVariables

      use local_BCS_variables

      implicit none

c Call variables

      integer    :: neq,bbcs(6,neq)

c Local variables

      integer    :: ieq,bcsq(6)

c Begin program

c Default boundary conditions

      bcsq = bbcs(:,IRHO)
      where (bcsq == DEF) bcsq = NEU
      bbcs(:,IRHO) = bcsq

      bcsq = bbcs(:,IVX)
      where (bcsq == DEF) bcsq = DIR
      bbcs(:,IVX) = bcsq

      bcsq = bbcs(:,IVY)
      if (     coords == 'cyl'
     .    .or. coords == 'hel'
     .    .or. coords == 'tor'
     .    .or. coords == 'car') then
        where (bcsq == DEF) bcsq =  NEU !On contravariant components
      else
        where (bcsq == DEF) bcsq = -NEU !On covariant components
      endif
      bbcs(:,IVY) = bcsq

      bcsq = bbcs(:,IVZ)
      if (     coords == 'cyl'
     .    .or. coords == 'hel'
     .    .or. coords == 'tor'
     .    .or. coords == 'car') then
        where (bcsq == DEF) bcsq =  NEU !On contravariant components
      else
        where (bcsq == DEF) bcsq = -NEU !On covariant components
      endif
      bbcs(:,IVZ) = bcsq

#if !defined(vec_pot)
      bcsq = bbcs(:,IBX)
      where (bcsq == DEF) bcsq = DIR
      bbcs(:,IBX) = bcsq

      bcsq = bbcs(:,IBY)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IBY) = bcsq

      bcsq = bbcs(:,IBZ)
      where (bcsq == DEF) bcsq = -NEU  !On covariant components
      bbcs(:,IBZ) = bcsq
#else
      bcsq = bbcs(:,IAX)
      where (bcsq == DEF) bcsq = DIR   !On contravariant components
cc      where (bcsq == DEF) bcsq = EXT   !On contravariant components
cc      where (bcsq == DEF) bcsq = -EQU  !On contravariant components
      bbcs(:,IAX) = bcsq

      bcsq = bbcs(:,IAY)
      where (bcsq == DEF) bcsq = -EQU  !On covariant components
      bbcs(:,IAY) = bcsq

      bcsq = bbcs(:,IAZ)
      where (bcsq == DEF) bcsq = -EQU  !On covariant components
      bbcs(:,IAZ) = bcsq
#endif

      bcsq = bbcs(:,ITMP)
      where (bcsq == DEF) bcsq = NEU !To allow isothermal case
      bbcs(:,ITMP) = bcsq

c Exceptions for specific equilibria

      select case (equil)

      case ('rfp1')

cc        bbcs(2,IBY) = EQU  !Imposed by equilibrium
cc
cc        bbcs(2,IBZ) = EQU  !Imposed by equilibrium

      case ('tmcyz','3dkh')

        where (bbcs(:,IVX) == DIR) bbcs(:,IVX) = -NEU
        where (bbcs(:,IVY) == NEU) bbcs(:,IVY) =  DIR

        where (bbcs(:,IBX) == DIR) bbcs(:,IBX) = -NEU
        where (bbcs(:,IBY) == NEU) bbcs(:,IBY) =  DIR

      case ('ppnch','ppn3d','ppnst')   !No stress BCs

cc        where (bbcs(:,IRHO) == NEU) bbcs(:,IRHO) = EQU

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX ) = EQU

      case ('ppnsl')                   !No slip BCs

cc        where (bbcs(:,IRHO) == NEU) bbcs(:,IRHO) = EQU

        where (bbcs(:,IVX ) == DIR) bbcs(:,IVX ) = EQU

        !Dirichlet boundary conditions for tangential components
        where (bbcs(:,IVY ) == NEU) bbcs(:,IVY ) = DIR
        where (bbcs(:,IVZ ) == NEU) bbcs(:,IVZ ) = DIR

c$$$      case ('ic')
c$$$
c$$$        bbcs(2,:) = EQU

      case ('vmec','tor')

        !Dirichlet boundary conditions for tangential components
        where (abs(bbcs(:,IVY )) == NEU) bbcs(:,IVY) = DIR
        where (abs(bbcs(:,IVZ )) == NEU) bbcs(:,IVZ) = DIR
cc        where (abs(bbcs(:,ITMP)) == NEU) bbcs(:,ITMP) = DIR

      end select

c BCs for auxiliary variables

      call defineAuxBCs

c Define BC sequence

      call setBCsched

c End subroutine

      end subroutine defineBoundaryConditions

# if 1

c imposeBoundaryConditions
c####################################################################
      subroutine imposeBoundaryConditions(varray,iigx,iigy,iigz)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz

      type(var_array),pointer :: varray

c Local variables

      integer    :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz,order

      real(8)    :: bbcov(3),bbcnv(3),b2,J0(3),E00(3)

      real(8),allocatable,dimension(:,:,:,:) :: v0

c Begin program

      order = 2   !Order of extrapolation in BCs

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Allocate auxiliary variables in local domain

      allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))

c Density BC

      call setBC(IRHO,nnx,nny,nnz
     .               ,varray%array_var(IRHO)%array
     .               ,u_0   %array_var(IRHO)%array
     .               ,varray%array_var(IRHO)%bconds
     .               ,iigx,iigy,iigz,iorder=order)

#if !defined(vec_pot)

c Magnetic field BC

c     Resistivity (needed for BCs; may depend on T)

      do k=0,nnz+1
        do j=0,nny+1
          do i=0,nnx+1
            eeta (i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
          enddo
        enddo
      enddo

c     BC setup

      bcnd(:,1) = varray%array_var(IBX)%bconds
      bcnd(:,2) = varray%array_var(IBY)%bconds
      bcnd(:,3) = varray%array_var(IBZ)%bconds

      bcnv(:,:,:,1) = varray%array_var(IBX)%array
      bcnv(:,:,:,2) = varray%array_var(IBY)%array
      bcnv(:,:,:,3) = varray%array_var(IBZ)%array

      v0(:,:,:,1) = u_0%array_var(IBX)%array
      v0(:,:,:,2) = u_0%array_var(IBY)%array
      v0(:,:,:,3) = u_0%array_var(IBZ)%array

c     Fill ghost nodes

      call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

      varray%array_var(IBX)%array = bcnv(:,:,:,1)
      varray%array_var(IBY)%array = bcnv(:,:,:,2)
      varray%array_var(IBZ)%array = bcnv(:,:,:,3)

#else

c Vector Pot. BC

c     Resistivity (needed for BCs; may depend on T)

      do k=0,nnz+1
        do j=0,nny+1
          do i=0,nnx+1
            eeta (i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
          enddo
        enddo
      enddo

c     BC setup

      bcnd(:,1) = varray%array_var(IAX)%bconds
      bcnd(:,2) = varray%array_var(IAY)%bconds
      bcnd(:,3) = varray%array_var(IAZ)%bconds

      acov(:,:,:,1) = varray%array_var(IAX)%array
      acov(:,:,:,2) = varray%array_var(IAY)%array
      acov(:,:,:,3) = varray%array_var(IAZ)%array

      v0(:,:,:,1) = u_0%array_var(IAX)%array
      v0(:,:,:,2) = u_0%array_var(IAY)%array
      v0(:,:,:,3) = u_0%array_var(IAZ)%array

c     Fill ghost nodes

      call setBC(IAX,3,nnx,nny,nnz,acnv,acov,v0,bcnd,iigx,iigy,iigz
     .          ,is_cnv=.false.,iorder=order)

c     Postprocessing

      varray%array_var(IAX)%array = acov(:,:,:,1)
      varray%array_var(IAY)%array = acov(:,:,:,2)
      varray%array_var(IAZ)%array = acov(:,:,:,3)

c Magnetic field BC

c     BC setup  (only topological constraints)

      bcsq = bcond
      where (bcsq == DEF) bcsq = DIR
      bcnd(:,1) = bcsq

cc      if (equil == 'spnch' .or. equil == 'rfp2') then
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cc        bcnd(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = -EXT !On covariant components
cc        bcnd(:,3) = bcsq
cc      else
        bcsq = bcond
        where (bcsq == DEF) bcsq = -NEU !On covariant components
        bcnd(:,2) = bcsq

        bcsq = bcond
        where (bcsq == DEF) bcsq = -NEU !On covariant components
        bcnd(:,3) = bcsq
cc      endif

c     B=curl(A)

      do k = 0,nnz+1
        do j = 0,nny+1
          do i = 0,nnx+1
            do icomp=1,3
              !Here, v_cov=a per previous step
              bcnv(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                                ,acov(:,:,:,1)
     .                                ,acov(:,:,:,2)
     .                                ,acov(:,:,:,3),icomp)
            enddo
          enddo
        enddo
      enddo

      v0 = bcnv

c     Fill ghost nodes

      call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

#endif

c Current BC

c     BC setup

#if !defined(vec_pot)
      bcnd(:,1) = varray%array_var(IBX)%bconds
      bcnd(:,2) = varray%array_var(IBY)%bconds
      bcnd(:,3) = varray%array_var(IBZ)%bconds
      where (bcnd == -NEU)
cc        bcnd = -DIR  !Use covariant components for tangential BCs
        bcnd = -EQU  !Use covariant components for tangential BCs
      end where
#else
cc      if (equil == 'spnch' .or. equil == 'rfp2') then
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = DIR
cc        bcnd(:,1) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = EXT
cccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc        bcnd(:,2) = bcsq
cc
cc        bcsq = bcond
cc        where (bcsq == DEF) bcsq = EXT
cccc        where (bcsq == DEF) bcsq = -EQU  !On covariant components
cc        bcnd(:,3) = bcsq
cc      else
        bcnd(:,1) = varray%array_var(IAX)%bconds
        bcnd(:,2) = varray%array_var(IAY)%bconds
        bcnd(:,3) = varray%array_var(IAZ)%bconds
cc      endif
#endif

      do k = 0,nnz+1
        do j = 0,nny+1
          do i = 0,nnx+1
            do icomp=1,3
              jcnv(i,j,k,icomp)=curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                                ,bcov(:,:,:,1)
     .                                ,bcov(:,:,:,2)
     .                                ,bcov(:,:,:,3),icomp)
            enddo
          enddo
        enddo
      enddo

      if (    equil == 'ppnch' .or. equil == 'ppn3d'
     .    .or.equil == 'ppnsl' .or. equil == 'ppnst' ) then       !Jt = alpha*Bt in cov representation
cc      if (sum(abs(E0)) /= 0d0) then  !Jt = alpha*Bt in cov representation

        i = nnx

        do k=1,nnz
          do j=1,nny
#if !defined(vec_pot)
            bbcnv = 0.5*(/ varray%array_var(IBX)%array(i  ,j,k)
     .                    ,varray%array_var(IBY)%array(i  ,j,k)
     .                    ,varray%array_var(IBZ)%array(i  ,j,k) /)
     .             +0.5*(/ varray%array_var(IBX)%array(i+1,j,k)
     .                    ,varray%array_var(IBY)%array(i+1,j,k)
     .                    ,varray%array_var(IBZ)%array(i+1,j,k) /)
#else
            bbcnv = 0.5*(/ bcnv(i  ,j,k,1)
     .                    ,bcnv(i  ,j,k,2)
     .                    ,bcnv(i  ,j,k,3) /)
     .             +0.5*(/ bcnv(i+1,j,k,1)
     .                    ,bcnv(i+1,j,k,2)
     .                    ,bcnv(i+1,j,k,3) /)
#endif
            bbcov = matmul(0.5*(gmetric%grid(iigx)%gsub(i  ,j,k,:,:)
     .                         +gmetric%grid(iigx)%gsub(i+1,j,k,:,:))
     .                    ,bbcnv)

            b2   = 2*dot_product(bbcnv,bbcov)
     .            /(gmetric%grid(iigx)%jac(i  ,j,k)
     .             +gmetric%grid(iigx)%jac(i+1,j,k))

            J0 = E0*2./(eeta(i,j,k)+eeta(i+1,j,k))

            v0(i+1,j,k,:) = scalarProduct(i,j,k,iigx,iigy,iigz
     .                                 ,J0   (1),J0   (2),J0   (3)
     .                                 ,bbcnv(1),bbcnv(2),bbcnv(3))/b2
     .                     *bbcov(:)
          enddo
        enddo

      else

        v0 = 0d0

      endif

c     Fill ghost nodes

      call setBC(IJX,3,nnx,nny,nnz,jcnv,jcov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c Velocity BC

c     BC setup

      bcnd(:,1) = varray%array_var(IVX)%bconds
      bcnd(:,2) = varray%array_var(IVY)%bconds
      bcnd(:,3) = varray%array_var(IVZ)%bconds

      if (.not.nc_eom_v) then
        where (varray%array_var(IRHO)%array /= 0d0)
          vcnv(:,:,:,1) = varray%array_var(IVX )%array
     .                    /varray%array_var(IRHO)%array
          vcnv(:,:,:,2) = varray%array_var(IVY )%array
     .                    /varray%array_var(IRHO)%array
          vcnv(:,:,:,3) = varray%array_var(IVZ )%array
     .                    /varray%array_var(IRHO)%array
        end where
      else
        vcnv(:,:,:,1) = varray%array_var(IVX)%array
        vcnv(:,:,:,2) = varray%array_var(IVY)%array
        vcnv(:,:,:,3) = varray%array_var(IVZ)%array
      endif

      if (    equil == 'ppnch' .or. equil == 'ppn3d'
     .    .or.equil == 'ppnsl' .or. equil == 'ppnst' ) then

        i = nnx

        do k=1,nnz
          do j=1,nny
#if !defined(vec_pot)
            bbcnv = 0.5*(/ varray%array_var(IBX)%array(i  ,j,k)
     .                    ,varray%array_var(IBY)%array(i  ,j,k)
     .                    ,varray%array_var(IBZ)%array(i  ,j,k) /)
     .             +0.5*(/ varray%array_var(IBX)%array(i+1,j,k)
     .                    ,varray%array_var(IBY)%array(i+1,j,k)
     .                    ,varray%array_var(IBZ)%array(i+1,j,k) /)
#else
            bbcnv = 0.5*(/ bcnv(i  ,j,k,1)
     .                    ,bcnv(i  ,j,k,2)
     .                    ,bcnv(i  ,j,k,3) /)
     .             +0.5*(/ bcnv(i+1,j,k,1)
     .                    ,bcnv(i+1,j,k,2)
     .                    ,bcnv(i+1,j,k,3) /)
#endif
            bbcov = matmul(0.5*(gmetric%grid(iigx)%gsub(i  ,j,k,:,:)
     .                         +gmetric%grid(iigx)%gsub(i+1,j,k,:,:))
     .                    ,bbcnv)

            b2   = 2*dot_product(bbcnv,bbcov)
     .            /(gmetric%grid(iigx)%jac(i  ,j,k)
     .             +gmetric%grid(iigx)%jac(i+1,j,k))

            if (b2 > 0d0) then
cc              E00(1) = E0(1) - 0.5*(eeta(i  ,j,k)*jcov(i  ,j,k,1)
cc     .                             +eeta(i+1,j,k)*jcov(i+1,j,k,1))
              E00(2) = E0(2) - 0.5*(eeta(i  ,j,k)*jcov(i  ,j,k,2)
     .                             +eeta(i+1,j,k)*jcov(i+1,j,k,2))
              E00(3) = E0(3) - 0.5*(eeta(i  ,j,k)*jcov(i  ,j,k,3)
     .                             +eeta(i+1,j,k)*jcov(i+1,j,k,3))

              v0(i+1,j,k,1) = -(bbcov(2)*E00(3)-bbcov(3)*E00(2))/b2 !ExB velocity
            else
              v0(i+1,j,k,1) = 0d0
            endif
          enddo
        enddo

      else

        v0(:,:,:,1) = u_0%array_var(IVX)%array

      endif

      v0(:,:,:,2) = u_0%array_var(IVY)%array
      v0(:,:,:,3) = u_0%array_var(IVZ)%array

c     Fill ghost nodes

      call setBC(IVX,3,nnx,nny,nnz,vcnv,vcov,v0,bcnd,iigx,iigy,iigz
     .          ,iorder=order)

c     Postprocessing

      if (.not.nc_eom_v) then
        varray%array_var(IVX)%array = vcnv(:,:,:,1)
     .                               *varray%array_var(IRHO)%array
        varray%array_var(IVY)%array = vcnv(:,:,:,2)
     .                               *varray%array_var(IRHO)%array
        varray%array_var(IVZ)%array = vcnv(:,:,:,3)
     .                               *varray%array_var(IRHO)%array
      else
        varray%array_var(IVX)%array = vcnv(:,:,:,1)
        varray%array_var(IVY)%array = vcnv(:,:,:,2)
        varray%array_var(IVZ)%array = vcnv(:,:,:,3)
      endif

c Temperature BCs

      call setBC(ITMP,nnx,nny,nnz
     .               ,varray%array_var(ITMP)%array
     .               ,u_0   %array_var(ITMP)%array
     .               ,varray%array_var(ITMP)%bconds
     .               ,iigx,iigy,iigz,iorder=order)

c Deallocate variables 

      deallocate(v0)

c End

      end subroutine imposeBoundaryConditions

#else

c imposeBoundaryConditions
c####################################################################
      subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz

      type(var_array),pointer :: varray

c Local variables

      integer :: i

c Interfaces

      INTERFACE
        subroutine applyDepVarBC(bc_inf,varray,vaux,iigx,iigy,iigz)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine applyDepVarBC
      END INTERFACE

      INTERFACE
        subroutine applyAuxVarBC(bc_inf,varray,vaux,iigx,iigy,iigz)
          use variable_setup
          integer    :: iigx,iigy,iigz,bc_inf(2)
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine applyAuxVarBC
      END INTERFACE

c Begin program

      call setBCsched

#if !defined(samrai)
      do i = 1,gv%nbc_seq
        if (gv%bc_seq(i,1)>0) then  !Process dep. variable
          call applyDepVarBC(gv%bc_seq(i,:),varray,gv%aux
     .                      ,iigx,iigy,iigz)
        else                          !Process aux. variable
          call applyAuxVarBC(gv%bc_seq(i,:),varray,gv%aux
     .                      ,iigx,iigy,iigz)
        endif
      enddo
#endif

c End program

      end subroutine imposeBoundaryConditions

c applyDepVarBC
c####################################################################
      subroutine applyDepVarBC(bc_inf,varray,vaux,iigx,iigy,iigz)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

c Local variables

      integer    :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz,order

      real(8)    :: bbcov(3),bbcnv(3),b2,J0(3),E00(3)

      real(8),allocatable,dimension(:,:,:,:) :: v0
      real(8),pointer,dimension(:,:,:,:) :: b_cnv,b_cov,j_cnv,j_cov

c Begin program

      order = 2   !Order of extrapolation in BCs

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Allocate auxiliary variables in local domain

      allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))

      j_cnv => vaux%vec_list(IJCNV)%vec
      j_cov => vaux%vec_list(IJCOV)%vec

#if defined(vec_pot)
      b_cnv => vaux%vec_list(IBCNV)%vec
      b_cov => vaux%vec_list(IBCOV)%vec
#endif

c Select variable

      select case(abs(bc_inf(1)))
      case(IRHO)

c     Density BCs

        if (bc_inf(2) /= 0) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be scalar (0)')
        endif

        call setBC(IRHO,nnx,nny,nnz
     .               ,varray%array_var(IRHO)%array
     .               ,u_0   %array_var(IRHO)%array
     .               ,varray%array_var(IRHO)%bconds
     .               ,iigx,iigy,iigz,iorder=order)

#if defined(vec_pot)
      case(IAX)

c     Vector Pot. BC

        if (bc_inf(2) /= 1) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be vector (1)')
        endif

c       Resistivity (needed for BCs; may depend on T)

        do k=0,nnz+1
          do j=0,nny+1
            do i=0,nnx+1
              eeta (i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
            enddo
          enddo
        enddo

c       BC setup

        bcnd(:,1) = varray%array_var(IAX)%bconds
        bcnd(:,2) = varray%array_var(IAY)%bconds
        bcnd(:,3) = varray%array_var(IAZ)%bconds

        acov(:,:,:,1) = varray%array_var(IAX)%array
        acov(:,:,:,2) = varray%array_var(IAY)%array
        acov(:,:,:,3) = varray%array_var(IAZ)%array

        v0(:,:,:,1) = u_0%array_var(IAX)%array
        v0(:,:,:,2) = u_0%array_var(IAY)%array
        v0(:,:,:,3) = u_0%array_var(IAZ)%array

c       Fill ghost nodes

        call setBC(IAX,3,nnx,nny,nnz,acnv,acov,v0,bcnd,iigx,iigy,iigz
     .            ,is_cnv=.false.,iorder=order)

        varray%array_var(IAX)%array = acov(:,:,:,1)
        varray%array_var(IAY)%array = acov(:,:,:,2)
        varray%array_var(IAZ)%array = acov(:,:,:,3)

c       Postprocessing: B=curl(A)

        do k = 0,nnz+1
          do j = 0,nny+1
            do i = 0,nnx+1
              do icomp=1,3
                b_cnv(i,j,k,icomp)
     .                 = curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                        ,acov(:,:,:,1)
     .                        ,acov(:,:,:,2)
     .                        ,acov(:,:,:,3),icomp)

                call transformFromCurvToCurv(i,j,k,iigx,iigy,iigz
     .            ,b_cov(i,j,k,1),b_cov(i,j,k,2),b_cov(i,j,k,3)
     .            ,b_cnv(i,j,k,1),b_cnv(i,j,k,2),b_cnv(i,j,k,3)
     .            ,.false.)
              enddo
            enddo
          enddo
        enddo
#else
      case(IBX)

c     Magnetic field BC

        if (bc_inf(2) /= 1) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be vector (1)')
        endif

c       Resistivity (needed for BCs; may depend on T)

        do k=0,nnz+1
          do j=0,nny+1
            do i=0,nnx+1
              eeta (i,j,k) = res(i,j,k,nnx,nny,nnz,iigx,iigy,iigz)
            enddo
          enddo
        enddo

c       BC setup

        bcnd(:,1) = varray%array_var(IBX)%bconds
        bcnd(:,2) = varray%array_var(IBY)%bconds
        bcnd(:,3) = varray%array_var(IBZ)%bconds

        bcnv(:,:,:,1) = varray%array_var(IBX)%array
        bcnv(:,:,:,2) = varray%array_var(IBY)%array
        bcnv(:,:,:,3) = varray%array_var(IBZ)%array

        v0(:,:,:,1) = u_0%array_var(IBX)%array
        v0(:,:,:,2) = u_0%array_var(IBY)%array
        v0(:,:,:,3) = u_0%array_var(IBZ)%array

c       Fill ghost nodes

        call setBC(IBX,3,nnx,nny,nnz,bcnv,bcov,v0,bcnd,iigx,iigy,iigz
     .            ,iorder=order)

        varray%array_var(IBX)%array = bcnv(:,:,:,1)
        varray%array_var(IBY)%array = bcnv(:,:,:,2)
        varray%array_var(IBZ)%array = bcnv(:,:,:,3)

c       Postprocessing: J=curl(B)

        do k = 0,nnz+1
          do j = 0,nny+1
            do i = 0,nnx+1
              do icomp=1,3
                j_cnv(i,j,k,icomp)=
     .                     curl2(i,j,k,nnx,nny,nnz,iigx,iigy,iigz
     .                          ,bcov(:,:,:,1)
     .                          ,bcov(:,:,:,2)
     .                          ,bcov(:,:,:,3),icomp)

                call transformFromCurvToCurv(i,j,k,iigx,iigy,iigz
     .            ,j_cov(i,j,k,1),j_cov(i,j,k,2),j_cov(i,j,k,3)
     .            ,j_cnv(i,j,k,1),j_cnv(i,j,k,2),j_cnv(i,j,k,3)
     .            ,.false.)
              enddo
            enddo
          enddo
        enddo
#endif

      case(IVX)

c     Velocity BC

        if (bc_inf(2) /= 1) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be vector (1)')
        endif

c       BC setup

        bcnd(:,1) = varray%array_var(IVX)%bconds
        bcnd(:,2) = varray%array_var(IVY)%bconds
        bcnd(:,3) = varray%array_var(IVZ)%bconds

        if (.not.nc_eom_v) then
          where (varray%array_var(IRHO)%array /= 0d0)
            vcnv(:,:,:,1) = varray%array_var(IVX )%array
     .                      /varray%array_var(IRHO)%array
            vcnv(:,:,:,2) = varray%array_var(IVY )%array
     .                      /varray%array_var(IRHO)%array
            vcnv(:,:,:,3) = varray%array_var(IVZ )%array
     .                      /varray%array_var(IRHO)%array
          end where
        else
          vcnv(:,:,:,1) = varray%array_var(IVX)%array
          vcnv(:,:,:,2) = varray%array_var(IVY)%array
          vcnv(:,:,:,3) = varray%array_var(IVZ)%array
        endif

        if (    equil == 'ppnch' .or. equil == 'ppn3d'
     .      .or.equil == 'ppnsl' .or. equil == 'ppnst' ) then

          i = nnx

          do k=1,nnz
            do j=1,nny
#if !defined(vec_pot)
              bbcnv = 0.5*(/ varray%array_var(IBX)%array(i  ,j,k)
     .                      ,varray%array_var(IBY)%array(i  ,j,k)
     .                      ,varray%array_var(IBZ)%array(i  ,j,k) /)
     .               +0.5*(/ varray%array_var(IBX)%array(i+1,j,k)
     .                      ,varray%array_var(IBY)%array(i+1,j,k)
     .                      ,varray%array_var(IBZ)%array(i+1,j,k) /)
#else
              bbcnv = 0.5*(/ b_cnv(i  ,j,k,1)
     .                      ,b_cnv(i  ,j,k,2)
     .                      ,b_cnv(i  ,j,k,3) /)
     .               +0.5*(/ b_cnv(i+1,j,k,1)
     .                      ,b_cnv(i+1,j,k,2)
     .                      ,b_cnv(i+1,j,k,3) /)
#endif
              bbcov = matmul(0.5*(gmetric%grid(iigx)%gsub(i  ,j,k,:,:)
     .                           +gmetric%grid(iigx)%gsub(i+1,j,k,:,:))
     .                      ,bbcnv)

              b2   = 2*dot_product(bbcnv,bbcov)
     .              /(gmetric%grid(iigx)%jac(i  ,j,k)
     .               +gmetric%grid(iigx)%jac(i+1,j,k))

              if (b2 > 0d0) then
cc                E00(1) = E0(1) - 0.5*(eeta(i  ,j,k)*jcov(i  ,j,k,1)
cc       .                             +eeta(i+1,j,k)*jcov(i+1,j,k,1))
                E00(2) = E0(2) - 0.5*(eeta(i  ,j,k)*jcov(i  ,j,k,2)
     .                               +eeta(i+1,j,k)*jcov(i+1,j,k,2))
                E00(3) = E0(3) - 0.5*(eeta(i  ,j,k)*jcov(i  ,j,k,3)
     .                               +eeta(i+1,j,k)*jcov(i+1,j,k,3))

                v0(i+1,j,k,1) = -(bbcov(2)*E00(3)-bbcov(3)*E00(2))/b2 !ExB velocity
              else
                v0(i+1,j,k,1) = 0d0
              endif
            enddo
          enddo

        else

          v0(:,:,:,1) = u_0%array_var(IVX)%array

        endif

        v0(:,:,:,2) = u_0%array_var(IVY)%array
        v0(:,:,:,3) = u_0%array_var(IVZ)%array

c       Fill ghost nodes

        call setBC(IVX,3,nnx,nny,nnz,vcnv,vcov,v0,bcnd,iigx,iigy,iigz
     .            ,iorder=order)

c       Postprocessing

        if (.not.nc_eom_v) then
          varray%array_var(IVX)%array = vcnv(:,:,:,1)
     .                                 *varray%array_var(IRHO)%array
          varray%array_var(IVY)%array = vcnv(:,:,:,2)
     .                                 *varray%array_var(IRHO)%array
          varray%array_var(IVZ)%array = vcnv(:,:,:,3)
     .                                 *varray%array_var(IRHO)%array
        else
          varray%array_var(IVX)%array = vcnv(:,:,:,1)
          varray%array_var(IVY)%array = vcnv(:,:,:,2)
          varray%array_var(IVZ)%array = vcnv(:,:,:,3)
        endif

      case(ITMP)

c     Temperature BCs

        if (bc_inf(2) /= 0) then
          call pstop('applyDepVarBC'
     .              ,'Inconsistent qty def: should be scalar (0)')
        endif

        call setBC(ITMP,nnx,nny,nnz
     .               ,varray%array_var(ITMP)%array
     .               ,u_0   %array_var(ITMP)%array
     .               ,varray%array_var(ITMP)%bconds
     .               ,iigx,iigy,iigz,iorder=order)

      end select

c Deallocate variables 

      deallocate(v0)
      nullify(j_cnv,j_cov,b_cnv,b_cov)

c End

      end subroutine applyDepVarBC

c applyAuxVarBC
c####################################################################
      subroutine applyAuxVarBC(bc_inf,varray,vaux,iigx,iigy,iigz)
c--------------------------------------------------------------------
c     Sets adequate boundary conditions on array structure varray.
c--------------------------------------------------------------------

      use local_BCS_variables

      use imposeBCinterface

      use auxiliaryVariables

      use variable_setup

      implicit none

c Call variables

      integer    :: iigx,iigy,iigz,bc_inf(2)

      type(var_array),pointer :: varray
      type(aux_array),pointer :: vaux

c Local variables

      integer    :: i,j,k,icomp,bcnd(6,3),bcsq(6),nnx,nny,nnz,order

      real(8)    :: bbcov(3),bbcnv(3),b2,J0(3),E00(3)

      real(8),allocatable,dimension(:,:,:,:) :: v0
      real(8),pointer,dimension(:,:,:,:) :: b_cnv,b_cov,j_cnv,j_cov

c Begin program

      order = 2   !Order of extrapolation in BCs

c Local grid sizes

      nnx = grid_params%nxv(iigx) 
      nny = grid_params%nyv(iigy)
      nnz = grid_params%nzv(iigz)

c Allocate auxiliary variables in local domain

      allocate(v0(0:nnx+1,0:nny+1,0:nnz+1,3))

      j_cnv => vaux%vec_list(IJCNV)%vec
      j_cov => vaux%vec_list(IJCOV)%vec

#if defined(vec_pot)
      b_cnv => vaux%vec_list(IBCNV)%vec
      b_cov => vaux%vec_list(IBCOV)%vec
#endif

c Set aux BCs

      call defineAuxBCs

c Select variable

      select case(abs(bc_inf(1)))
      case(IJX)

c     Current BC

        bcnd = vaux%vec_list(IJCNV)%bconds

        if (    equil == 'ppnch' .or. equil == 'ppn3d'
     .      .or.equil == 'ppnsl' .or. equil == 'ppnst' ) then       !Jt = alpha*Bt in cov representation

          i = nnx

          do k=1,nnz
            do j=1,nny
#if !defined(vec_pot)
              bbcnv = 0.5*(/ varray%array_var(IBX)%array(i  ,j,k)
     .                      ,varray%array_var(IBY)%array(i  ,j,k)
     .                      ,varray%array_var(IBZ)%array(i  ,j,k) /)
     .               +0.5*(/ varray%array_var(IBX)%array(i+1,j,k)
     .                      ,varray%array_var(IBY)%array(i+1,j,k)
     .                      ,varray%array_var(IBZ)%array(i+1,j,k) /)
#else
              bbcnv = 0.5*(/ b_cnv(i  ,j,k,1)
     .                      ,b_cnv(i  ,j,k,2)
     .                      ,b_cnv(i  ,j,k,3) /)
     .               +0.5*(/ b_cnv(i+1,j,k,1)
     .                      ,b_cnv(i+1,j,k,2)
     .                      ,b_cnv(i+1,j,k,3) /)
#endif
              bbcov = matmul(0.5*(gmetric%grid(iigx)%gsub(i  ,j,k,:,:)
     .                           +gmetric%grid(iigx)%gsub(i+1,j,k,:,:))
     .                      ,bbcnv)

              b2   = 2*dot_product(bbcnv,bbcov)
     .              /(gmetric%grid(iigx)%jac(i  ,j,k)
     .               +gmetric%grid(iigx)%jac(i+1,j,k))

              J0 = E0*2./(eeta(i,j,k)+eeta(i+1,j,k))

              v0(i+1,j,k,:) = scalarProduct(i,j,k,iigx,iigy,iigz
     .                                   ,J0   (1),J0   (2),J0   (3)
     .                                   ,bbcnv(1),bbcnv(2),bbcnv(3))/b2
     .                       *bbcov(:)
            enddo
          enddo

        else

          v0 = 0d0

        endif

c       Fill ghost nodes

        call setBC(IJX,3,nnx,nny,nnz,j_cnv,j_cov,v0,bcnd,iigx,iigy,iigz
     .            ,iorder=order)

#if defined(vec_pot)
      case(IBX)

c     Magnetic field BC

        bcnd = vaux%vec_list(IBCNV)%bconds

        v0 = b_cnv

        call setBC(IBX,3,nnx,nny,nnz,b_cnv,b_cov,v0,bcnd,iigx,iigy,iigz
     .            ,iorder=order)
#endif
      end select

c Deallocate variables 

      deallocate(v0)
      nullify(j_cnv,j_cov,b_cnv,b_cov)

c End

      end subroutine applyAuxVarBC

#endif
