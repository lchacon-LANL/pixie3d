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

      use nk_mod

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
        if (my_rank == 0) write (*,*) 'rho=',minrho
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

      E_ni = 0d0

      if (source)  then           !Source processing
        if (jit == -2) jcov_0 = jcov
      else
        jcov_0 = 0d0
      endif

#if !defined(vec_pot)
      E_ni(:,:,:,1) = E_ni(:,:,:,1) - eeta*(jcov  (:,:,:,1)
     .                                     -jcov_0(:,:,:,1))
      E_ni(:,:,:,2) = E_ni(:,:,:,2) - eeta*(jcov  (:,:,:,2)
     .                                     -jcov_0(:,:,:,2))
      E_ni(:,:,:,3) = E_ni(:,:,:,3) - eeta*(jcov  (:,:,:,3)
     .                                     -jcov_0(:,:,:,3))
#endif

c Electron forces (covariant vector)

      if (di > 0d0) then

        !Standard Hall formulation: E_ni = di*(heta*lap(j) + grad(pe)/rho)
        if (.not.hall_new) then
 
          !Electron velocity
cc          vecnv(:,:,:,1) = vcnv(:,:,:,1) - di*jcnv(:,:,:,1)/rho
cc          vecnv(:,:,:,2) = vcnv(:,:,:,2) - di*jcnv(:,:,:,2)/rho
cc          vecnv(:,:,:,3) = vcnv(:,:,:,3) - di*jcnv(:,:,:,3)/rho

          vecnv = vcnv

          !Electron pressure
          pp = rho*tmp

          do k=0,nz+1
            do j=0,ny+1
              do i=0,nx+1
              !Electron pressure

                call grad(i,j,k,nx,ny,nz,igx,igy,igz,pp
     .                   ,cov(1),cov(2),cov(3))

              !Electron hall (normally included directly via ve in vexB tem)

                cov(1) = cov(1) - (jcnv(i,j,k,2)*bcnv(i,j,k,3)
     .                            -jcnv(i,j,k,3)*bcnv(i,j,k,2))
     .                            /gmetric%grid(igx)%jac(i,j,k)
                cov(2) = cov(2) - (jcnv(i,j,k,3)*bcnv(i,j,k,1)
     .                            -jcnv(i,j,k,1)*bcnv(i,j,k,3))
     .                            /gmetric%grid(igx)%jac(i,j,k)
                cov(3) = cov(3) - (jcnv(i,j,k,1)*bcnv(i,j,k,2)
     .                            -jcnv(i,j,k,2)*bcnv(i,j,k,1))
     .                            /gmetric%grid(igx)%jac(i,j,k)

                E_ni(i,j,k,:) = E_ni(i,j,k,:) + di*cov/rho(i,j,k)

cc                cov(1) = cov(1) + (vecnv(i,j,k,2)*bcnv(i,j,k,3)
cc     .                            -vecnv(i,j,k,3)*bcnv(i,j,k,2))
cc     .                            /gmetric%grid(igx)%jac(i,j,k)
cc                cov(2) = cov(2) + (vecnv(i,j,k,3)*bcnv(i,j,k,1)
cc     .                            -vecnv(i,j,k,1)*bcnv(i,j,k,3))
cc     .                            /gmetric%grid(igx)%jac(i,j,k)
cc                cov(3) = cov(3) + (vecnv(i,j,k,1)*bcnv(i,j,k,2)
cc     .                            -vecnv(i,j,k,2)*bcnv(i,j,k,1))
cc     .                            /gmetric%grid(igx)%jac(i,j,k)
cc
cc                E_ni(i,j,k,:) = E_ni(i,j,k,:) + cov

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

          pp = (a_p - 1d0)*rho*tmp  !Ion pressure

          if (jit < 0) vcov_n = vcov  !Store time n info

cc          do k=1,nz
cc            do j=1,ny
cc              do i=1,nx
          do k=0,nz+1
            do j=0,ny+1
              do i=0,nx+1

              !Temporal term

                 E_ni(i,j,k,:) = E_ni(i,j,k,:)
     $                          -di*(vcov(i,j,k,:)-vcov_n(i,j,k,:))
     $                          *cnp/dt/alpha

              !Add -grad(pi)/rho

                if (a_p - 1d0 > 0d0) then
                  call grad(i,j,k,nx,ny,nz,igx,igy,igz,pp
     .                   ,cov(1),cov(2),cov(3))

                  E_ni(i,j,k,:) = E_ni(i,j,k,:) - di*cov/rho(i,j,k)
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
 
                !Find curl(v)
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
        if (heta > 0d0) then

#if !defined(vec_pot)
          if (solenoidal) then

            !Solenoidal: de*nu_e*curl(lap(j))
            do k=1,nz
              do j=1,ny
                do i=1,nx
                   vlap(i,j,k,:)
     .                 = veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,jcov
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

cc        vlap_a = -jcov
          do k=1,nz
            do j=1,ny
              do i=1,nx
cc              vlap_a(i,j,k,:) =
cc       .           veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,acov,vol=.false.)
                cnv=-curlcurl(i,j,k,nx,ny,nz,igx,igy,igz,acov,vol=.false.)
                call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .             ,vlap_a(i,j,k,1),vlap_a(i,j,k,2),vlap_a(i,j,k,3)
     .             ,cnv(1),cnv(2),cnv(3),.false.)
              enddo
            enddo
          enddo

          bcnd(:,1) = varray%array_var(IAX)%bconds
          bcnd(:,2) = varray%array_var(IAY)%bconds
          bcnd(:,3) = varray%array_var(IAZ)%bconds

          call setBC(IJX,3,nx,ny,nz,vdummy,vlap_a,vzeros,bcnd,igx,igy,igz
     .              ,iorder=2,is_cnv=.false.)

#endif

       endif

      !Resistive MHD
      else

        E_ni = 0d0

        vecnv = vcnv

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
