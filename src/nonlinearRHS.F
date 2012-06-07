#if defined(flux_rhs)

c nonlinearRHS
c#################################################################
      subroutine nonlinearRHS(i,j,k,igx,igy,igz,varray,ff)

c-----------------------------------------------------------------
c     This function computes the nonlinear rhs for the 
c     2-D Hall MHD problem. Nonlinear functions Fi should be
c     coded with the following convention:
c        dt Ui + Fi(Uj) = 0
c     where Ui are the dependent variables, given in varray.
c     The result of Fi(Uj) at node (i,j) for all equations neqd
c     is returned in the vector ff of dimension neqd.
c-----------------------------------------------------------------

      use variable_setup

      use timeStepping

      use nlfunction_setup

      use equilibrium

      use transport_params
      
      use operators

      implicit none

c Call variables

      real(8) :: ff(neqd)

      integer :: i,j,k,igx,igy,igz

      type(var_array),pointer :: varray

c Local variables

      integer :: ieq,nx,ny,nz,ig,jg,kg
      real(8) :: dxx,dyy,dzz

c Begin program

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = grid_params%dxh(ig)
      dyy = grid_params%dyh(jg)
      dzz = grid_params%dzh(kg)

c Assemble fluxes

      ff = (gv%flux_x(i+1,j,k,:)-gv%flux_x(i,j,k,:))/dxx
     .    +(gv%flux_y(i,j+1,k,:)-gv%flux_y(i,j,k,:))/dyy
     .    +(gv%flux_z(i,j,k+1,:)-gv%flux_z(i,j,k,:))/dzz
     .    + gv%src   (i,j,k,:)

c End

      end subroutine nonlinearRHS

#else

c nonlinearRHS
c#################################################################
      subroutine nonlinearRHS(i,j,k,igx,igy,igz,varray,ff)

c-----------------------------------------------------------------
c     This function computes the nonlinear rhs for the 
c     2-D Hall MHD problem. Nonlinear functions Fi should be
c     coded with the following convention:
c        dt Ui + Fi(Uj) = 0
c     where Ui are the dependent variables, given in varray.
c     The result of Fi(Uj) at node (i,j) for all equations neqd
c     is returned in the vector ff of dimension neqd.
c-----------------------------------------------------------------

      use variable_setup

      use timeStepping

      use nlfunction_setup

      use equilibrium

      use transport_params
      
      use operators

      implicit none

c Call variables

      real(8) :: ff(neqd)

      integer :: i,j,k,igx,igy,igz

      type(var_array),pointer :: varray

c Local variables

      integer    :: ieq,nx,ny,nz

      real(8)    :: cnv(3),cov(3)

      ! rho equation
      real(8)    :: advc,diffus

      ! tmp equation
      real(8)    :: joule,viscous,heat_flx,heat_src,gm,sgn

c Begin program

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      if (solve_rho) then
        advc = flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,rho,advect)
cc        advc = flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                  ,vcnv,rho,advect,sp=bcSP(),zip_vel=.true.)  !sp is for conservation

        if (dd /= 0d0) then
          diffus = dd*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,rho)
        else
          diffus = 0d0
        endif

        ff(IRHO) = advc - diffus
      else
        ff(IRHO) = 0d0
      endif

c Faraday's law <<<<<<<<<<<<<<<<<<<<<<

#if !defined(vec_pot)

      !Magnetic field
      ff(IBX:IBZ)= curl(i,j,k,nx,ny,nz,igx,igy,igz,E_ni)

#else

      !Vector potential
      ff(IAX:IAZ) = E_ni(i,j,k,:)

#endif

c Temperature <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      if ((gamma == 1d0) .and. (.not.solve_prs)) then

        ff(ITMP) = 0d0

      else

c       Heat flux

        heat_flx = 0d0

        if (chi /= 0d0) then
          heat_flx =-chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,tmp
cc     .                            ,rho,vol=.false.)  !chi_perp=chi*rho
     .                            )                  !chi_perp=chi
        endif

        if (chi_par /= chi .and.(.not.lagrangian)) then
          heat_flx = heat_flx
     .              -(chi_par-chi)
     .              *par_diff(i,j,k,nx,ny,nz,igx,igy,igz,tmp
     .                   ,bhat(:,:,:,1),bhat(:,:,:,2),bhat(:,:,:,3),2
     .                   ,symmetric=.true.,limited=.false.,vol=.false.)
        endif

        if (.not.solve_prs) heat_flx = heat_flx/rho(i,j,k)

c       Advection

        if (solve_prs) then
          gm = gamma
          sgn =-1d0
        else
          gm = gamma - 1d0
          sgn = 1d0
        endif

        advc=     gm    *flx_advec   (i,j,k,nx,ny,nz,igx,igy,igz
     .                               ,    vscnv,tmp,advect)
     .      -sgn*(gm-1.)*flx_nc_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                               ,sgn*vscnv,tmp,1)

c       Sources

        if (.not.adiabatic) then
          !Joule heating
          joule = eeta(i,j,k)
cc     .           *dot_product(jcnv(i,j,k,:),jcov(i,j,k,:))/jac
     .           *scalarProduct_ijk(i,j,k,igx,jcnv(i,j,k,:)
     .                                       ,jcov(i,j,k,:))

          !Ion viscosity
          heat_src = joule
     .              +viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz
     .                               ,nuu,vcnv,vcnv)

          !Electron viscosity
          !Comment out for old TM Hall MHD performance results
          if (di > 0d0) then
            if (fake_ve) then
              heat_src = heat_src
     .               +viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz
     .                                ,h_eta,vefcnv,vecnv)
            else
              heat_src = heat_src
     .               +viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz
     .                                ,h_eta,vecnv ,vecnv)
            endif
          endif

          if (.not.solve_prs) heat_src = heat_src/rho(i,j,k)

          if (heat_src < 0d0) then
            write (*,*) 'Negative heat source'
            stop
          endif
        else
          heat_src = 0d0
        endif

        !Add contributions
        ff(ITMP) = advc + (gamma-1d0)*(heat_flx-heat_src)/a_p

      endif

c EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<

c     JxB volume force

      if (nc_eom_jxb) then
        cnv = EOM_jxb(i,j,k,nx,ny,nz,igx,igy,igz
cc     .               ,bcnv,jcnv,.true.
     .               ,bcov,jcov,.false.
     .               ,nc_eom_jxb)
      else
        cnv = EOM_jxb(i,j,k,nx,ny,nz,igx,igy,igz
     .               ,bcnv,bcnv,.true.
     .               ,nc_eom_jxb)
      endif

      ff(IVX:IVZ) = -cnv

c     grad(p) volume force

      cnv = EOM_gp(i,j,k,nx,ny,nz,igx,igy,igz,rho,tmp,nc_eom_gp)

      ff(IVX:IVZ) = ff(IVX:IVZ) + cnv

c     Stress tensor: div(Pi)

      if (slava .or. di == 0d0) then
        cnv = div_pi(i,j,k,:) !Match Slava's EOM
      else
        cnv = div_pi(i,j,k,:) + div_pe(i,j,k,:)
      endif

      ff(IVX:IVZ) = ff(IVX:IVZ) + cnv

c     Advection

      if (nc_eom_v) then
        cnv = tensor_nc_advc(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,vcnv
     .                      ,v_advect)

        ff(IVX:IVZ) = ff(IVX:IVZ)/rho(i,j,k) + cnv  !1/rho factor because
                                                    !we evolve dv/dt
      else
        cnv = tensor_advc(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,pcnv,v_advect)

        ff(IVX:IVZ) = ff(IVX:IVZ) + cnv
      endif

c End

      end subroutine nonlinearRHS

#endif

c defineTSParameters
c####################################################################
      subroutine defineTSParameters
c--------------------------------------------------------------------
c     Calculates nonlinear residuals. 
c--------------------------------------------------------------------

      use parameters

      use equilibrium

      use timeStepping

      use transport_params

      use nlfunction_setup

      implicit none

c Call variables

c Local variables

c Begin program

      if (bdf2) then
        bdfp  = cnp
        bdfn  = cn
        bdfnm = cnm
        cnf   = 0d0
      else
        bdfp  = 1d0
        bdfn  =-1d0
        bdfnm = 0d0
        cnf   = cnfactor
      endif

      if (dt == 0d0) call pstop('defineTSparams','dt=0!!!')

      one_over_dt = 1d0/dt

      !Setup enslavement of rho
      if (.not.solve_rho) then
        one_over_dt(IRHO) = 0d0
        cnf(IRHO) = 0d0
      endif

c End program

      end subroutine defineTSParameters

c set_nlf_flux_x
c#################################################################
      function set_nlf_flux_x(igr,nx,ny,nz) result(flx_x)

c-----------------------------------------------------------------
c     This function computes the X-flux of the nonlinear function
c     evaluation, with convention:
c            dt U + div(flux) + src = 0
c     It returns the flux value at X-faces for all variables
c     over the whole mesh.
c-----------------------------------------------------------------

      use variable_setup

      use timeStepping

      use nlfunction_setup

      use equilibrium

      use transport_params
      
      use operators

      implicit none

c Call variables

      integer :: igr,nx,ny,nz

      real(8) :: flx_x(nx+1,ny,nz,neqd)

c Local variables

      integer    :: i,j,k,ig,jg,kg!,ip,im,jp,jm,kp,km,hex,hey,hez
      integer    :: ieq

      real(8)    :: cnv(3),cov(3),jac

      ! rho equation
      real(8)    :: advc,diffus

      ! tmp equation
      real(8)    :: joule,viscous,heat_flx,heat_src,gm,sgn

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 0,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (solve_rho) then
              advc = scl_advc_flx_x(i,j,k,nx,ny,nz,igr,vcnv,rho,advect)

              if (dd /= 0d0) then
                diffus = dd*lap_flx_x_ijk(i,j,k,igr,rho)
              else
                diffus = 0d0
              endif

              flx_x(i+1,j,k,IRHO) = advc - diffus
            else
              flx_x(i+1,j,k,IRHO) = 0d0
            endif

            !Faraday's law <<<<<<<<<<<<<<<<<<<<<<
#if !defined(vec_pot)
            !Magnetic field
            flx_x(i+1,j,k,IBX:IBZ) = 0d0
#else
            !Vector potential
            flx_x(i+1,j,k,IAX:IAZ) = 0d0
#endif

            !Temperature <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if ((gamma == 1d0) .and. (.not.solve_prs)) then
              flx_x(i+1,j,k,ITMP) = 0d0
            else
              if (dd /= 0d0.and.solve_prs) then
                heat_flx =-chi*lap_flx_x_ijk(i,j,k,igr,tmp)
              else
                heat_flx = 0d0
              endif

              if (solve_prs) then
                gm = gamma
                sgn =-1d0
              else
                gm = gamma - 1d0
                sgn = 1d0
              endif

              advc = gm*scl_advc_flx_x(i,j,k,nx,ny,nz,igr,vscnv,tmp
     .                                ,advect)

              flx_x(i+1,j,k,ITMP) = advc + (gamma-1d0)*heat_flx/a_p
            endif

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (.not.nc_eom_v) then
              !JxB
              if (.not.nc_eom_jxb) then
                vec1 => bcnv
                vec2 => bcnv
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_jxb_x)
                nullify(vec1,vec2)
              else
                cnv = 0d0
              endif

              flx_x(i+1,j,k,IVX:IVZ) = -cnv

              !grad(P)
              if (.not.nc_eom_gp) then
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_grad_p_x)
                nullify(sc1,sc2)
              else
                cnv = 0d0
              endif

              flx_x(i+1,j,k,IVX:IVZ) = flx_x(i+1,j,k,IVX:IVZ) + cnv

              !Advection
              adv_meth = v_advect
              vec1 => vcnv
              vec2 => pcnv
              cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                           ,tnsr_advc_x)
              nullify(vec1,vec2)

              flx_x(i+1,j,k,IVX:IVZ) = flx_x(i+1,j,k,IVX:IVZ) + cnv
            else
              flx_x(i+1,j,k,IVX:IVZ) = 0d0
            endif

          enddo
        enddo
      enddo

c End program

      end function set_nlf_flux_x

c set_nlf_flux_y
c#################################################################
      function set_nlf_flux_y(igr,nx,ny,nz) result(flx_y)

c-----------------------------------------------------------------
c     This function computes the X-flux of the nonlinear function
c     evaluation, with convention:
c            dt U + div(flux) + src = 0
c     It returns the flux value at X-faces for all variables
c     over the whole mesh.
c-----------------------------------------------------------------

      use variable_setup

      use timeStepping

      use nlfunction_setup

      use equilibrium

      use transport_params
      
      use operators

      implicit none

c Call variables

      integer :: igr,nx,ny,nz

      real(8) :: flx_y(nx,ny+1,nz,neqd)

c Local variables

      integer    :: i,j,k,ig,jg,kg!,ip,im,jp,jm,kp,km,hex,hey,hez
      integer    :: ieq

      real(8)    :: cnv(3),cov(3),jac

      ! rho equation
      real(8)    :: advc,diffus

      ! tmp equation
      real(8)    :: joule,viscous,heat_flx,heat_src,gm,sgn

c Begin program

      do k = 1,nz
        do j = 0,ny
          do i = 1,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (solve_rho) then
              advc = scl_advc_flx_y(i,j,k,nx,ny,nz,igr,vcnv,rho,advect)

              if (dd /= 0d0) then
                diffus = dd*lap_flx_y_ijk(i,j,k,igr,rho)
              else
                diffus = 0d0
              endif

              flx_y(i,j+1,k,IRHO) = advc - diffus
            else
              flx_y(i,j+1,k,IRHO) = 0d0
            endif

            !Faraday's law <<<<<<<<<<<<<<<<<<<<<<
#if !defined(vec_pot)
            !Magnetic field
            flx_y(i,j+1,k,IBX:IBZ) = 0d0
#else
            !Vector potential
            flx_y(i,j+1,k,IAX:IAZ) = 0d0
#endif

            !Temperature <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if ((gamma == 1d0) .and. (.not.solve_prs)) then
              flx_y(i,j+1,k,ITMP) = 0d0
            else
              if (dd /= 0d0.and.solve_prs) then
                heat_flx =-chi*lap_flx_y_ijk(i,j,k,igr,tmp)
              else
                heat_flx = 0d0
              endif

              if (solve_prs) then
                gm = gamma
                sgn =-1d0
              else
                gm = gamma - 1d0
                sgn = 1d0
              endif

              advc = gm*scl_advc_flx_y(i,j,k,nx,ny,nz,igr,vscnv,tmp
     .                                ,advect)

              flx_y(i,j+1,k,ITMP) = advc + (gamma-1d0)*heat_flx/a_p
            endif

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (.not.nc_eom_v) then
              !JxB
              if (.not.nc_eom_jxb) then
                vec1 => bcnv
                vec2 => bcnv
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_jxb_y)
                nullify(vec1,vec2)
              else
                cnv = 0d0
              endif

              flx_y(i,j+1,k,IVX:IVZ) = -cnv

              !grad(P)
              if (.not.nc_eom_gp) then
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_grad_p_y)
                nullify(sc1,sc2)
              else
                cnv = 0d0
              endif

              flx_y(i,j+1,k,IVX:IVZ) = flx_y(i,j+1,k,IVX:IVZ) + cnv

              !Advection
              adv_meth = v_advect
              vec1 => vcnv
              vec2 => pcnv
              cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                           ,tnsr_advc_y)
              nullify(vec1,vec2)

              flx_y(i,j+1,k,IVX:IVZ) = flx_y(i,j+1,k,IVX:IVZ) + cnv
            else
              flx_y(i,j+1,k,IVX:IVZ) = 0d0
            endif
          enddo
        enddo
      enddo

c End program

      end function set_nlf_flux_y

c set_nlf_flux_z
c#################################################################
      function set_nlf_flux_z(igr,nx,ny,nz) result(flx_z)

c-----------------------------------------------------------------
c     This function computes the X-flux of the nonlinear function
c     evaluation, with convention:
c            dt U + div(flux) + src = 0
c     It returns the flux value at X-faces for all variables
c     over the whole mesh.
c-----------------------------------------------------------------

      use variable_setup

      use timeStepping

      use nlfunction_setup

      use equilibrium

      use transport_params
      
      use operators

      implicit none

c Call variables

      integer :: igr,nx,ny,nz

      real(8) :: flx_z(nx,ny,nz+1,neqd)

c Local variables

      integer    :: i,j,k,ig,jg,kg!,ip,im,jp,jm,kp,km,hex,hey,hez
      integer    :: ieq

      real(8)    :: cnv(3),cov(3),jac

      ! rho equation
      real(8)    :: advc,diffus

      ! tmp equation
      real(8)    :: joule,viscous,heat_flx,heat_src,gm,sgn

c Begin program

      do k = 0,nz
        do j = 1,ny
          do i = 1,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (solve_rho) then
              advc = scl_advc_flx_z(i,j,k,nx,ny,nz,igr,vcnv,rho,advect)

              if (dd /= 0d0) then
                diffus = dd*lap_flx_z_ijk(i,j,k,igr,rho)
              else
                diffus = 0d0
              endif

              flx_z(i,j,k+1,IRHO) = advc - diffus
            else
              flx_z(i,j,k+1,IRHO) = 0d0
            endif

            !Faraday's law <<<<<<<<<<<<<<<<<<<<<<
#if !defined(vec_pot)
            !Magnetic field
            flx_z(i,j,k+1,IBX:IBZ) = 0d0
#else
            !Vector potential
            flx_z(i,j,k+1,IAX:IAZ) = 0d0
#endif

            !Temperature <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if ((gamma == 1d0) .and. (.not.solve_prs)) then
              flx_z(i,j,k+1,ITMP) = 0d0
            else
              if (dd /= 0d0.and.solve_prs) then
                heat_flx =-chi*lap_flx_z_ijk(i,j,k,igr,tmp)
              else
                heat_flx = 0d0
              endif

              if (solve_prs) then
                gm = gamma
                sgn =-1d0
              else
                gm = gamma - 1d0
                sgn = 1d0
              endif

              advc = gm*scl_advc_flx_z(i,j,k,nx,ny,nz,igr,vscnv,tmp
     .                                ,advect)

              flx_z(i,j,k+1,ITMP) = advc + (gamma-1d0)*heat_flx/a_p
            endif

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (.not.nc_eom_v) then
              !JxB
              if (.not.nc_eom_jxb) then
                vec1 => bcnv
                vec2 => bcnv
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_jxb_z)
                nullify(vec1,vec2)
              else
                cnv = 0d0
              endif

              flx_z(i,j,k+1,IVX:IVZ) = -cnv

              !grad(P)
              if (.not.nc_eom_gp) then
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_grad_p_z)
                nullify(sc1,sc2)
              else
                cnv = 0d0
              endif

              flx_z(i,j,k+1,IVX:IVZ) = flx_z(i,j,k+1,IVX:IVZ) + cnv

              !Advection
              adv_meth = v_advect
              vec1 => vcnv
              vec2 => pcnv
              cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                           ,tnsr_advc_z)
              nullify(vec1,vec2)

              flx_z(i,j,k+1,IVX:IVZ) = flx_z(i,j,k+1,IVX:IVZ) + cnv
            else
              flx_z(i,j,k+1,IVX:IVZ) = 0d0
            endif

          enddo
        enddo
      enddo

c End program

      end function set_nlf_flux_z

c set_nlf_src
c#################################################################
      function set_nlf_src(igr,nx,ny,nz) result(src)

c-----------------------------------------------------------------
c     This function computes the X-flux of the nonlinear function
c     evaluation, with convention:
c            dt U + div(flux) + src = 0
c     It returns the flux value at X-faces for all variables
c     over the whole mesh.
c-----------------------------------------------------------------

      use variable_setup

      use timeStepping

      use nlfunction_setup

      use equilibrium

      use transport_params
      
      use operators

      implicit none

c Call variables

      integer :: igr,nx,ny,nz

      real(8) :: src(nx,ny,nz,neqd)

c Local variables

      integer    :: i,j,k,ig,jg,kg
      integer    :: ieq

      real(8)    :: cnv(3),cov(3),jac

      ! rho equation
      real(8)    :: advc,diffus

      ! tmp equation
      real(8)    :: joule,viscous,heat_src,heat_flx,gm,sgn

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            src(i,j,k,IRHO) = 0d0

            !Faraday's law <<<<<<<<<<<<<<<<<<<<<<
#if !defined(vec_pot)
            !Magnetic field
            src(i,j,k,IBX:IBZ) = curl(i,j,k,nx,ny,nz,igr,igr,igr,E_ni)
#else
            !Vector potential
            src(i,j,k,IAX:IAZ) = E_ni(i,j,k,:)
#endif

            !Temperature <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if ((gamma == 1d0) .and. (.not.solve_prs)) then
              src(i,j,k,ITMP) = 0d0
            else
              !Heat Flux
              !(parallel piece; should be made conservative)
              if (chi_par /= chi .and.(.not.lagrangian)) then
                heat_flx = -(chi_par-chi)
     .              *par_diff(i,j,k,nx,ny,nz,igr,igr,igr,tmp
     .                   ,bhat(:,:,:,1),bhat(:,:,:,2),bhat(:,:,:,3),2
     .                   ,symmetric=.true.,limited=.false.,vol=.false.)
              endif

              if (.not.solve_prs) then
                if (chi /= 0d0) then
                  heat_flx=heat_flx
     .                    -chi*laplacian(i,j,k,nx,ny,nz,igr,igr,igr,tmp
cc     .                            ,rho,vol=.false.)  !chi_perp=chi*rho
     .                            )                  !chi_perp=chi
                endif

                heat_flx = heat_flx/rho(i,j,k)
              endif

              !Advection (non-conservative piece)
              if (solve_prs) then
                gm = gamma
                sgn =-1d0
              else
                gm = gamma - 1d0
                sgn = 1d0
              endif

              advc =-sgn*(gm-1.)*flx_nc_advec(i,j,k,nx,ny,nz,igr,igr,igr
     .                                       ,sgn*vscnv,tmp,1)

              !Sources
              if (.not.adiabatic) then
                !Joule heating
                joule = eeta(i,j,k)
     .                 *scalarProduct_ijk(i,j,k,igr,jcnv(i,j,k,:)
     .                                             ,jcov(i,j,k,:))

                !Ion viscosity
                heat_src = joule
     .                    +viscous_heat_src(i,j,k,nx,ny,nz,igr,igr,igr
     .                                     ,nuu,vcnv,vcnv)

                !Electron viscosity
                !Comment out for old TM Hall MHD performance results
                if (di > 0d0) then
                  if (fake_ve) then
                    heat_src = heat_src
     .                     +viscous_heat_src(i,j,k,nx,ny,nz,igr,igr,igr
     .                                      ,h_eta,vefcnv,vecnv)
                  else
                    heat_src = heat_src
     .                     +viscous_heat_src(i,j,k,nx,ny,nz,igr,igr,igr
     .                                      ,h_eta,vecnv ,vecnv)
                  endif
                endif

                if (.not.solve_prs) heat_src = heat_src/rho(i,j,k)

                if (heat_src < 0d0) then
                  write (*,*) 'Negative heat source'
                  stop
                endif
              else
                heat_src = 0d0
              endif

              src(i,j,k,ITMP) = advc
     .                        + (gamma-1d0)*(heat_flx-heat_src)/a_p
            endif

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (nc_eom_v) then
              !JxB force
              if (nc_eom_jxb) then
                cnv = EOM_jxb(i,j,k,nx,ny,nz,igr,igr,igr
cc             .               ,bcnv,jcnv,.true.
     .                       ,bcov,jcov,.false.
     .                       ,nc_eom_jxb)
              else
                cnv = EOM_jxb(i,j,k,nx,ny,nz,igr,igr,igr
     .                       ,bcnv,bcnv,.true.
     .                       ,nc_eom_jxb)
              endif

              src(i,j,k,IVX:IVZ) = -cnv

              !grad(p) volume force
              cnv = EOM_gp(i,j,k,nx,ny,nz,igr,igr,igr,rho,tmp,nc_eom_gp)

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) + cnv

              !Stress tensor: div(Pi)
              if (slava .or. di == 0d0) then
                cnv = div_pi(i,j,k,:) !Match Slava's EOM
              else
                cnv = div_pi(i,j,k,:) + div_pe(i,j,k,:)
              endif

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) + cnv

              !Advection
              cnv = tensor_nc_advc(i,j,k,nx,ny,nz,igr,igr,igr,vcnv,vcnv
     .                            ,v_advect)

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ)/rho(i,j,k) + cnv !1/rho factor because
                                                                       !we evolve dv/dt
            else
              !JxB
              if (nc_eom_jxb) then
                cnv = EOM_jxb(i,j,k,nx,ny,nz,igr,igr,igr
cc     .                       ,bcnv,jcnv,.true.
     .                       ,bcov,jcov,.false.
     .                       ,nc_eom_jxb)
              else
                vec1 => bcnv
                vec2 => bcnv
                cnv =  dtnsr_src_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                              ,eom_jxb_x,eom_jxb_y,eom_jxb_z)
                nullify(vec1,vec2)
              endif

              src(i,j,k,IVX:IVZ) = -cnv

              !Grad(p)
              if (nc_eom_gp) then
                cnv=EOM_gp(i,j,k,nx,ny,nz,igr,igr,igr,rho,tmp,nc_eom_gp)
              else
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_src_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                          ,eom_grad_p_x,eom_grad_p_y,eom_grad_p_z)
                nullify(sc1,sc2)
              endif

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) + cnv

              !Stress tensor: div(Pi) (should be made conservative)
              if (slava .or. di == 0d0) then
                cnv = div_pi(i,j,k,:) !Match Slava's EOM
              else
                cnv = div_pi(i,j,k,:) + div_pe(i,j,k,:)
              endif

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) + cnv

              !Advection
              adv_meth = v_advect
              vec1 => vcnv
              vec2 => pcnv
              cnv = dtnsr_src_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                           ,tnsr_advc_x,tnsr_advc_y,tnsr_advc_z)
              nullify(vec1,vec2)

              src(i,j,k,IVX:IVZ)=src(i,j,k,IVX:IVZ) + cnv
            endif
          enddo
        enddo
      enddo

c End program

      end function set_nlf_src
