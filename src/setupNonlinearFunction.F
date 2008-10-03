c setupNonlinearFunction
c#################################################################
      subroutine setupNonlinearFunction(varray)
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use parameters

      use variable_setup

      use nlfunction_setup

      use equilibrium

      use constants

      use timeStepping

      use error

      use imposeBCinterface

      use newtongm

      implicit none

c Call variables

      type(var_array),pointer :: varray

c Local variables

      integer    :: i,j,k,ig,jg,kg

      integer    :: dim,loc,ibc,ieq,bctype,bcnd(6,3)

      real(8)    :: dh,vmax,cnv(3),cov(3),minrho,mintmp,lrmin,ltmin
     .             ,pp(ilom:ihip,jlom:jhip,klom:khip)

c Interfaces

      INTERFACE
        subroutine imposeBoundaryConditions (varray,iigx,iigy,iigz)
        use variable_setup
        integer    :: iigx,iigy,iigz
        type(var_array),pointer :: varray
        end subroutine imposeBoundaryConditions
      END INTERFACE

c Begin program

      igx = 1
      igy = 1
      igz = 1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Impose boundary conditions and find auxiliary quantities

      call imposeBoundaryConditions(varray,igx,igy,igz)

c Set aliases

      rho => varray%array_var(IRHO)%array
      px  => varray%array_var(IVX )%array
      py  => varray%array_var(IVY )%array
      pz  => varray%array_var(IVZ )%array
#if !defined(vec_pot)
      bx  => varray%array_var(IBX )%array
      by  => varray%array_var(IBY )%array
      bz  => varray%array_var(IBZ )%array
#else
      ax  => varray%array_var(IAX )%array
      ay  => varray%array_var(IAY )%array
      az  => varray%array_var(IAZ )%array
#endif
      tmp => varray%array_var(ITMP)%array

c Safeguards

      !Mass density positivity
      minrho = minval(rho(1:nx,1:ny,1:nz))

#if defined(petsc)
      lrmin = minrho
      call MPI_Allreduce(lrmin,minrho,1
     .                  ,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif

      if (minrho < 0d0) then
        close (urecord)  !Close record file before exiting
cc        if (my_rank == 0) write (*,*) 'rho=',minrho
        write(*,*)'proc=',my_rank,'rho=',lrmin,'loc='
     $            ,minloc(rho(1:nx,1:ny,1:nz))
        call pstop('setupNonlinearFunction',
     .             'Warning: negative densities are occurring')
      endif

      !Temperature positivity (enforce)
      mintmp = minval(tmp(1:nx,1:ny,1:nz))

#if defined(petsc)
      ltmin = mintmp
      call MPI_Allreduce(ltmin,mintmp,1
     .                  ,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif

      if (mintmp < 0d0) then
        close (urecord) !Close record file before exiting
        write(*,*)'proc=',my_rank,'tmp=',ltmin,'loc='
     $            ,minloc(tmp(1:nx,1:ny,1:nz))
        call pstop('setupNonlinearFunction',
     .             'Warning: negative temperatures are occurring')
      endif

c Viscosity (resistivity found in BC routine)

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            nuu(i,j,k) = vis(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

c Resistive term

cc      if (source .and. jit == -2) jcov_0 = jcov

#if !defined(vec_pot)
      if (source) then
        E_ni(:,:,:,1) = -eeta*(jcov  (:,:,:,1)
     .                        -jcov_0(:,:,:,1))
        E_ni(:,:,:,2) = -eeta*(jcov  (:,:,:,2)
     .                        -jcov_0(:,:,:,2))
        E_ni(:,:,:,3) = -eeta*(jcov  (:,:,:,3)
     .                        -jcov_0(:,:,:,3))
      else
        E_ni(:,:,:,1) = -eeta*jcov  (:,:,:,1)
        E_ni(:,:,:,2) = -eeta*jcov  (:,:,:,2)
        E_ni(:,:,:,3) = -eeta*jcov  (:,:,:,3)
      endif
#else
      E_ni = 0d0
#endif

c Electron forces (covariant vector)

      !Resistive MHD
      if (di == 0d0) then

        vecnv = vcnv

      !Hall MHD
      else

        !Standard Hall formulation: E_ni = di*(heta*lap(j) + grad(pe)/rho)
        if (.not.hall_new) then
 
          !Electron velocity
#if defined(vec_pot)
          vecnv(:,:,:,1) = vcnv(:,:,:,1) - di*jcnv(:,:,:,1)/rho
          vecnv(:,:,:,2) = vcnv(:,:,:,2) - di*jcnv(:,:,:,2)/rho
          vecnv(:,:,:,3) = vcnv(:,:,:,3) - di*jcnv(:,:,:,3)/rho
#else
          vecnv = vcnv
#endif

          do k=0,nz+1
            do j=0,ny+1
              do i=0,nx+1

              !Electron pressure (only tmp*grad(rho)/rho survives curl)

                call grad(i,j,k,nx,ny,nz,igx,igy,igz,rho
     .                   ,cov(1),cov(2),cov(3))

                E_ni(i,j,k,:) = E_ni(i,j,k,:)
     $                        + di*tmp(i,j,k)*cov/rho(i,j,k)

#if defined(vec_pot)
                call grad(i,j,k,nx,ny,nz,igx,igy,igz,tmp
     .                   ,cov(1),cov(2),cov(3))

                E_ni(i,j,k,:) = E_ni(i,j,k,:) + di*cov
#endif
              !Electron hall (normally included directly via ve in vexB tem)

cc                cov(1) = (vecnv(i,j,k,2)*bcnv(i,j,k,3)
cc     .                   -vecnv(i,j,k,3)*bcnv(i,j,k,2))
cc     .                   /gmetric%grid(igx)%jac(i,j,k)
cc                cov(2) = (vecnv(i,j,k,3)*bcnv(i,j,k,1)
cc     .                   -vecnv(i,j,k,1)*bcnv(i,j,k,3))
cc     .                   /gmetric%grid(igx)%jac(i,j,k)
cc                cov(3) = (vecnv(i,j,k,1)*bcnv(i,j,k,2)
cc     .                   -vecnv(i,j,k,2)*bcnv(i,j,k,1))
cc     .                   /gmetric%grid(igx)%jac(i,j,k)
cc
cc                E_ni(i,j,k,:) = E_ni(i,j,k,:) + cov
#if !defined(vec_pot)
                cov(1) = (jcnv(i,j,k,2)*bcnv(i,j,k,3)
     .                   -jcnv(i,j,k,3)*bcnv(i,j,k,2))
     .                   /gmetric%grid(igx)%jac(i,j,k)
                cov(2) = (jcnv(i,j,k,3)*bcnv(i,j,k,1)
     .                   -jcnv(i,j,k,1)*bcnv(i,j,k,3))
     .                   /gmetric%grid(igx)%jac(i,j,k)
                cov(3) = (jcnv(i,j,k,1)*bcnv(i,j,k,2)
     .                   -jcnv(i,j,k,2)*bcnv(i,j,k,1))
     .                   /gmetric%grid(igx)%jac(i,j,k)

                E_ni(i,j,k,:) = E_ni(i,j,k,:) - di*cov/rho(i,j,k)
#endif

              enddo
            enddo
          enddo

          !Impose topological BCs
          bcnd(:,1) = bcond
          bcnd(:,2) = bcond
          bcnd(:,3) = bcond

          call setBC(IVX,3,nx,ny,nz,vdummy,E_ni,vzeros,bcnd
     .              ,igx,igy,igz,iorder=2,is_cnv=.false.)

        !Alternate Hall formulation
        else

          !Form E_ni = di*(-d_t v + v x curl(v) + nu*lap(v) + heta*lap(j) - grad(pi)/rho)
          vecnv = vcnv

          if (jit < 0) vcov_n = vcov  !Store time n info

          do k=0,nz+1
            do j=0,ny+1
              do i=0,nx+1

              !Temporal term

                E_ni(i,j,k,:) = E_ni(i,j,k,:)
     $                         -di*(vcov(i,j,k,:)-vcov_n(i,j,k,:))
     $                            *cnp/dt/alpha

              !Ion pressure term (only T*grad(rho)/rho survives curl)

                if (a_p - 1d0 > 0d0) then
                  call grad(i,j,k,nx,ny,nz,igx,igy,igz,rho
     .                   ,cov(1),cov(2),cov(3))

                  E_ni(i,j,k,:) = E_ni(i,j,k,:)
     $                          - di*(a_p-1d0)*tmp(i,j,k)*cov/rho(i,j,k)

#if defined(vec_pot)
                  call grad(i,j,k,nx,ny,nz,igx,igy,igz,tmp
     .                   ,cov(1),cov(2),cov(3))

                  E_ni(i,j,k,:) = E_ni(i,j,k,:) - di*(a_p-1d0)*cov
#endif
                endif

cc              !Add v.grad(v)
cc                nabla_v = fnabla_v_cov(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                                ,vcov(:,:,:,1)
cc     .                                ,vcov(:,:,:,2)
cc     .                                ,vcov(:,:,:,3),0)
cc
cc                do ieq=1,3
cc                  cov(ieq) =( vcov(i,j,k,1)*nabla_v(1,ieq)
cc     .                       +vcov(i,j,k,2)*nabla_v(2,ieq)
cc     .                       +vcov(i,j,k,3)*nabla_v(3,ieq))
cc     .                      /gmetric%grid(igx)%jac(i,j,k)
cc                enddo
cc
cc                E_ni(i,j,k,:) = E_ni(i,j,k,:) - di*cov

              !Add v x curl(v)
 
                do ieq=1,3
                  cnv(ieq) = curl2(i,j,k,nx,ny,nz,igx,igy,igz
     .                            ,vcov(:,:,:,1)
     .                            ,vcov(:,:,:,2)
     .                            ,vcov(:,:,:,3),ieq)
                enddo

                !Find v x curl(v)
                cov(1) = (vcnv(i,j,k,2)*cnv(3)
     .                   -vcnv(i,j,k,3)*cnv(2))
                cov(2) = (vcnv(i,j,k,3)*cnv(1)
     .                   -vcnv(i,j,k,1)*cnv(3))
                cov(3) = (vcnv(i,j,k,1)*cnv(2)
     .                   -vcnv(i,j,k,2)*cnv(1))

                cov = cov/gmetric%grid(igx)%jac(i,j,k)

                E_ni(i,j,k,:) = E_ni(i,j,k,:) + di*cov

cc              !Add div(nu*grad(v))/rho
cc
cc                if (nu > 0d0) then
cc                  cov = veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,vcov,nuu
cc     .                            ,vol=.false.)/rho(i,j,k)
cc
cc                  E_ni(i,j,k,:) = E_ni(i,j,k,:) + di*cov
cc                endif

cc              !Add hyperresistive term (heta*lap(j))
cc
cc                if (heta > 0d0) then
cc                  cov = veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,jcov
cc     .                            ,vol=.false.)
cc
cc                  E_ni(i,j,k,:) = E_ni(i,j,k,:) + di*heta*cov
cc                endif

              enddo
            enddo
          enddo

cc          !Impose BCs on dissipation (perfect conductor E x n = 0)
cc          bcnd(:,1) = varray%array_var(IBX)%bconds
cc          bcnd(:,2) = varray%array_var(IBY)%bconds
cc          bcnd(:,3) = varray%array_var(IBZ)%bconds
cc          where (bcnd == -NEU) bcnd = -EQU  !Use covariant components for tangential BCs
cc
cc          call setBC(IJX,3,nx,ny,nz,vdummy,E_ni,vzeros,bcnd
cc     .              ,igx,igy,igz,iorder=2,is_cnv=.false.)

          !Impose topological BCs
          bcnd(:,1) = bcond
          bcnd(:,2) = bcond
          bcnd(:,3) = bcond

          call setBC(IVX,3,nx,ny,nz,vdummy,E_ni,vzeros,bcnd
     .              ,igx,igy,igz,iorder=2,is_cnv=.false.)

          !Add viscosity (with BCs)
          do k=1,nz
            do j=1,ny
              do i=1,nx

              !Add div(nu*grad(v))/rho

                if (nu > 0d0) then
                  vlap(i,j,k,:) =
     .                  veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,vcov,nuu
     .                            ,vol=.false.)/rho(i,j,k)
                endif

cc              !Add hyperresistive term (heta*lap(j))
cc
cc                if (heta > 0d0) then
cc                  cov = veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,jcov
cc     .                            ,vol=.false.)
cc
cc                  vlap(i,j,k,:) = vlap(i,j,k,:) + heta*cov
cccc                  vlap(i,j,k,:) = heta*cov
cc                endif

              enddo
            enddo
          enddo

          !Impose BCs on dissipation (extrapolation)
          bcnd(:,1) = bcond
          bcnd(:,2) = bcond
          bcnd(:,3) = bcond
          where (bcnd == DEF) bcnd =-EXT !Extrapolation BCs on covariant components

          call setBC(IVX,3,nx,ny,nz,vdummy,vlap,vzeros,bcnd
     .              ,igx,igy,igz,iorder=2,is_cnv=.false.)

          E_ni = E_ni + di*vlap

        endif

        !Add hyperresistive term
        if (di*heta > 0d0) then

#if !defined(vec_pot)
          if (solenoidal) then

            !Solenoidal: E_ni = E_ni + di*nu_e*lap(j-j0)
            vdummy = jcov - jcov_0
            do k=1,nz
              do j=1,ny
                do i=1,nx
                   vlap(i,j,k,:)
     .                 = veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,vdummy
     .                              ,vol=.false.)
                enddo
              enddo
            enddo

            !Impose BCs on dissipation (extrapolation)
            bcnd(:,1) = bcond
            bcnd(:,2) = bcond
            bcnd(:,3) = bcond
            where (bcnd == DEF) bcnd =-EXT !Extrapolation BCs on covariant components

            call setBC(IVX,3,nx,ny,nz,vdummy,vlap,vzeros,bcnd
     .                ,igx,igy,igz,iorder=2,is_cnv=.false.)

            E_ni = E_ni + di*heta*vlap

          else

            !Not solenoidal, but compact di*nu_e*lap^2(B)
            do k=1,nz
              do j=1,ny
                do i=1,nx
                  vlap(i,j,k,:)
     .                =veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,bcnv
     .                             ,.false.,vol=.false.)
                enddo
              enddo
            enddo

cc            bcnd(:,1) = bcond
cc            bcnd(:,2) = bcond
cc            bcnd(:,3) = bcond
cc            where (bcnd == DEF) bcnd = EXT !Extrapolation BCs on contravariant components
            bcnd(:,1) = varray%array_var(IBX)%bconds
            bcnd(:,2) = varray%array_var(IBY)%bconds
            bcnd(:,3) = varray%array_var(IBZ)%bconds
            where (bcnd == -NEU) bcnd = -EQU !Use covariant components for tangential BCs

            call setBC(IJX,3,nx,ny,nz,vdummy,vlap,vzeros,bcnd
     .                ,igx,igy,igz,iorder=2,is_cnv=.true.)
 
          endif

#else

          vlap_a = -(jcov - jcov_0)
cc          do k=1,nz
cc            do j=1,ny
cc              do i=1,nx
cccc              vlap_a(i,j,k,:) =
cccc       .           veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,acov,vol=.false.)
cc                cnv=-curlcurl(i,j,k,nx,ny,nz,igx,igy,igz,acov
cc     .                       ,vol=.false.)
cc                call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .             ,vlap_a(i,j,k,1),vlap_a(i,j,k,2),vlap_a(i,j,k,3)
cc     .             ,cnv(1),cnv(2),cnv(3),.false.)
cc              enddo
cc            enddo
cc          enddo
cc
cc          bcnd(:,1) = varray%array_var(IAX)%bconds
cc          bcnd(:,2) = varray%array_var(IAY)%bconds
cc          bcnd(:,3) = varray%array_var(IAZ)%bconds
cc
cc          call setBC(IJX,3,nx,ny,nz,vdummy,vlap_a,vzeros,bcnd
cc     .              ,igx,igy,igz,iorder=2,is_cnv=.false.)

#endif

       endif

      endif

c Other auxiliary variables

      pcnv(:,:,:,1) = px
      pcnv(:,:,:,2) = py
      pcnv(:,:,:,3) = pz

c SI operator

      if (k_si > 0d0) then
        b_n(:,:,:,1) = u_n%array_var(IBX)%array
        b_n(:,:,:,2) = u_n%array_var(IBY)%array
        b_n(:,:,:,3) = u_n%array_var(IBZ)%array

        p_n(:,:,:,1) = 2*u_n%array_var(IRHO)%array
     .                  *u_n%array_var(ITMP)%array

        dv_dt(:,:,:,1) = (vcnv(:,:,:,1)-u_n%array_var(IVX)%array)/dt
        dv_dt(:,:,:,2) = (vcnv(:,:,:,2)-u_n%array_var(IVY)%array)/dt
        dv_dt(:,:,:,3) = (vcnv(:,:,:,3)-u_n%array_var(IVZ)%array)/dt

        max_dv_dt = maxval(dv_dt)
      endif

c End program

      end subroutine setupNonlinearFunction

c killNonlinearFunction
c#################################################################
      subroutine killNonlinearFunction
c------------------------------------------------------------------
c     This function calculates auxiliary quantities for the
c     Jacobian-free product
c------------------------------------------------------------------

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

c Begin program

cc#if !defined(vec_pot)
cc      nullify(rho,bx,by,bz,px,py,pz,tmp)
cc#else
cc      nullify(rho,ax,ay,az,px,py,pz,tmp)
cc#endif

c End program

      end subroutine killNonlinearFunction
