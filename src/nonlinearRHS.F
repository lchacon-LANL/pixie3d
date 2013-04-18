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

      real(8)    :: cnv(3),cov(3),jac

      ! rho equation
      real(8)    :: advc,diffus

      ! tmp equation
      real(8)    :: joule,viscous,heat_flx,heat_src,gm,sgn

c Begin program

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      jac = gmetric%grid(igx)%jac(i,j,k)

c Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      if (solve_rho) then
        advc = flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,rho,advect)
cc        advc = flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                  ,vcnv,rho,advect,sp=bcSP())  !sp is for conservation

        if (dd /= 0d0) then
          diffus = dd*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,rho)
        else
          diffus = 0d0
        endif

        ff(IRHO) = (advc - diffus)*jac
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

c Total kinetic pressure equation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
cc
cc      if ((gamma == 1d0) .and. (.not.use_p_eom)) then
cc
cc        ff(ITMP) = 0d0
cc
cc      else
cc
ccc       Heat flux
cc
cc        heat_flx = 0d0
cc
cc        if (chi /= 0d0) then
cc          heat_flx =-chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,tmp
cccc     .                            ,rho,vol=.false.)  !chi_perp=chi*rho
cc     .                            )                  !chi_perp=chi
cc        endif
cc
cc        if (chi_par /= chi .and.(.not.lagrangian)) then
cc          heat_flx = heat_flx
cc     .              -(chi_par-chi)
cc     .              *par_diff(i,j,k,nx,ny,nz,igx,igy,igz,tmp
cc     .                   ,bhat(:,:,:,1),bhat(:,:,:,2),bhat(:,:,:,3),2
cc     .                   ,symmetric=.true.,limited=.false.,vol=.false.)
cc        endif
cc
cc        if (.not.use_p_eom) heat_flx = heat_flx/rho(i,j,k)
cc
ccc       Advection
cc
cc        if (use_p_eom) then
cc          gm = gamma
cc          sgn =-1d0
cc        else
cc          gm = gamma - 1d0
cc          sgn = 1d0
cc        endif
cc
cc        advc=     gm    *flx_advec   (i,j,k,nx,ny,nz,igx,igy,igz
cc     .                               ,    vscnv,tmp,advect)
cc     .      -sgn*(gm-1.)*flx_nc_advec(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                               ,sgn*vscnv,tmp,1)
cc
ccc       Sources
cc
cc        if (.not.adiabatic) then
cc          !Joule heating
cc          joule = eeta(i,j,k)
cccc     .           *dot_product(jcnv(i,j,k,:),jcov(i,j,k,:))/jac
cc     .           *scalarProduct_ijk(i,j,k,igx,jcnv(i,j,k,:)
cc     .                                       ,jcov(i,j,k,:))
cc
cc          !Ion viscosity
cc          heat_src = joule
cc     .              +viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                               ,nuu,vcnv,vcnv)
cc
cc          !Electron viscosity
cc          !Comment out for old TM Hall MHD performance results
cc          if (di > 0d0) then
cc            if (fake_ve) then
cc              heat_src = heat_src
cc     .               +viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                                ,h_eta,vefcnv,vecnv)
cc            else
cc              heat_src = heat_src
cc     .               +viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                                ,h_eta,vecnv ,vecnv)
cc            endif
cc          endif
cc
cc          if (.not.use_p_eom) heat_src = heat_src/rho(i,j,k)
cc
cc          if (heat_src < 0d0) then
cc            write (*,*) 'Negative heat source'
cc            stop
cc          endif
cc        else
cc          heat_src = 0d0
cc        endif
cc
cc        !Add contributions
cc        ff(ITMP) = (advc + (gamma-1d0)*(heat_flx-heat_src)/a_p)*jac
cc
cc      endif

c       Heat flux

        heat_flx = 0d0

        if (chi /= 0d0) then
          if (use_p_diff) then
            heat_flx =-chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,prs      !Total pressure
cc     .                            ,rho)  !chi_perp=chi*rho
     .                            )      !chi_perp=chi
          else
            heat_flx =-chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,tmp  !Total temperature
cc     .                            ,rho)  !chi_perp=chi*rho
     .                            )      !chi_perp=chi
          endif
        endif

        if (chi_par /= chi .and.(.not.lagrangian)) then
          heat_flx = heat_flx
     .              -(chi_par-chi)
     .              *par_diff(i,j,k,nx,ny,nz,igx,igy,igz,tmpe
     .                   ,bhat(:,:,:,1),bhat(:,:,:,2),bhat(:,:,:,3),2
     .                   ,symmetric=.true.,limited=.false.,vol=.false.)
        endif

c       Advection

        gm = gamma
        sgn =-1d0

        advc=     gm    *flx_advec   (i,j,k,nx,ny,nz,igx,igy,igz
     .                               ,vscnv,prs,advect)
cc     .                               ,vscnv,prs,advect,sp=bcSP())
     .          +(gm-1.)*flx_nc_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                               ,vscnv,prs,1,reverse=.true.)
cc     .      -sgn*(gm-1.)*flx_nc_advec(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                               ,sgn*vscnv,prs,1)

c       Sources

        if (.not.adiabatic) then
          !Joule heating
          joule = eeta(i,j,k)
     .           *scalarProduct_ijk(i,j,k,igx,jcnv(i,j,k,:)
     .                                       ,jcov(i,j,k,:))

          !Ion viscosity
          heat_src = joule
     .              +viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz
     .                               ,nuu,vcnv,vcnv)

          !Electron viscosity
          !Comment out for old TM Hall MHD performance results
          if (di > 0d0.and.(.not.slava)) then
c$$$            if (fake_ve) then
c$$$              heat_src = heat_src
c$$$     .               +viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz
c$$$     .                                ,h_eta,vefcnv,vefcnv)
c$$$            else
              heat_src = heat_src
     .               +viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz
     .                                ,h_eta,vecnv ,vecnv)
c$$$            endif
          endif

          if (heat_src < 0d0) then
            write (*,*) 'Negative heat source'
            stop
          endif
        else
          heat_src = 0d0
        endif

        !Add contributions
        ff(ITMP) = (advc + (gamma-1d0)*(heat_flx-heat_src))*jac

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

c     Stress tensor: div(Pi+Pe)

c$$$      if (slava) then
        ff(IVX:IVZ) = ff(IVX:IVZ) + div_pi(i,j,k,:)
c$$$      else
c$$$        ff(IVX:IVZ) = ff(IVX:IVZ) + div_pi(i,j,k,:) + div_pe(i,j,k,:)
c$$$      endif

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

      if (alt__eom()) ff(IVY) = ff(IVY)*jac

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

      integer :: i,j,k,ig,jg,kg

      real(8) :: cnv(3),cov(3),jac

      ! rho equation
      real(8) :: advc,diffus

      ! tmp equation
      real(8) :: joule,viscous,heat_flx,heat_src,gm,sgn

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 0,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (solve_rho) then
              advc = scl_advc_flx_x(i,j,k,nx,ny,nz,igr,vcnv(:,:,:,1)
     .                             ,rho,advect)!,sp=bcSP())

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
            if (chi /= 0d0) then
              if (use_p_diff) then
                heat_flx =-chi*lap_flx_x_ijk(i,j,k,igr,prs)
              else
                heat_flx =-chi*lap_flx_x_ijk(i,j,k,igr,tmp)
              endif
            else
              heat_flx = 0d0
            endif

            advc = scl_advc_flx_x(i,j,k,nx,ny,nz,igr,vscnv(:,:,:,1)
     .                           ,prs,advect)!,sp=bcSP())

            flx_x(i+1,j,k,ITMP) = gamma*advc + (gamma-1d0)*heat_flx

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (nc_eom_v) then
              flx_x(i+1,j,k,IVX:IVZ) = 0d0
            else
              !JxB
              if (nc_eom_jxb) then
                cnv = 0d0
              else
                vec1 => bcnv
                vec2 => bcnv
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_jxb_x)
                nullify(vec1,vec2)
              endif

              flx_x(i+1,j,k,IVX:IVZ) = cnv

              !grad(P)
              if (nc_eom_gp) then
                cnv = 0d0
              else
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_grad_p_x)
                nullify(sc1,sc2)
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

      integer :: i,j,k,ig,jg,kg

      real(8) :: cnv(3),cov(3),jac

      ! rho equation
      real(8) :: advc,diffus

      ! tmp equation
      real(8) :: joule,viscous,heat_flx,heat_src,gm,sgn

c Begin program

      do k = 1,nz
        do j = 0,ny
          do i = 1,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (solve_rho) then
              advc = scl_advc_flx_y(i,j,k,nx,ny,nz,igr,vcnv(:,:,:,2)
     .                             ,rho,advect)

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
            if (chi /= 0d0) then
              if (use_p_diff) then
                heat_flx =-chi*lap_flx_y_ijk(i,j,k,igr,prs)
              else
                heat_flx =-chi*lap_flx_y_ijk(i,j,k,igr,tmp)
              endif
            else
              heat_flx = 0d0
            endif

            advc = scl_advc_flx_y(i,j,k,nx,ny,nz,igr,vscnv(:,:,:,2)
     .                           ,prs,advect)

            flx_y(i,j+1,k,ITMP) = gamma*advc + (gamma-1d0)*heat_flx

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (nc_eom_v) then
              flx_y(i,j+1,k,IVX:IVZ) = 0d0
            else
              !JxB
              if (nc_eom_jxb) then
                cnv = 0d0
              else
                vec1 => bcnv
                vec2 => bcnv
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_jxb_y)
                nullify(vec1,vec2)
              endif

              flx_y(i,j+1,k,IVX:IVZ) = cnv

              !grad(P)
              if (nc_eom_gp) then
                cnv = 0d0
              else
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_grad_p_y)
                nullify(sc1,sc2)
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

      integer :: i,j,k,ig,jg,kg

      real(8) :: cnv(3),cov(3),jac

      ! rho equation
      real(8) :: advc,diffus

      ! tmp equation
      real(8) :: joule,viscous,heat_flx,heat_src,gm,sgn

c Begin program

      do k = 0,nz
        do j = 1,ny
          do i = 1,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (solve_rho) then
              advc = scl_advc_flx_z(i,j,k,nx,ny,nz,igr,vcnv(:,:,:,3)
     .                             ,rho,advect)

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
            if (chi /= 0d0) then
              if (use_p_diff) then
                heat_flx =-chi*lap_flx_z_ijk(i,j,k,igr,prs)
              else
                heat_flx =-chi*lap_flx_z_ijk(i,j,k,igr,tmp)
              endif
            else
              heat_flx = 0d0
            endif

            advc = scl_advc_flx_z(i,j,k,nx,ny,nz,igr,vscnv(:,:,:,3)
     .                           ,prs,advect)

            flx_z(i,j,k+1,ITMP) = gamma*advc + (gamma-1d0)*heat_flx

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (nc_eom_v) then
              flx_z(i,j,k+1,IVX:IVZ) = 0d0
            else
              !JxB
              if (nc_eom_jxb) then
                cnv = 0d0
              else
                vec1 => bcnv
                vec2 => bcnv
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_jxb_z)
                nullify(vec1,vec2)
              endif

              flx_z(i,j,k+1,IVX:IVZ) = cnv

              !grad(P)
              if (nc_eom_gp) then
                cnv = 0d0
              else
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_flx_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                             ,eom_grad_p_z)
                nullify(sc1,sc2)
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

      integer :: i,j,k,ig,jg,kg

      real(8) :: cnv(3),cov(3),jac

      ! rho equation
      real(8) :: advc,diffus

      ! tmp equation
      real(8) :: joule,viscous,heat_src,heat_flx,gm,sgn

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx

            jac = gmetric%grid(igr)%jac(i,j,k)

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

            !Heat Flux
            !(parallel piece; should be made conservative)
            if (chi_par /= chi .and.(.not.lagrangian)) then
              heat_flx = -(chi_par-chi)
     .            *par_diff(i,j,k,nx,ny,nz,igr,igr,igr,tmpe
     .                 ,bhat(:,:,:,1),bhat(:,:,:,2),bhat(:,:,:,3),2
     .                 ,symmetric=.true.,limited=.false.,vol=.false.)
            else
              heat_flx = 0d0
            endif

            !Advection (non-conservative piece, with reversed sign for velocity)
            advc = flx_nc_advec(i,j,k,nx,ny,nz,igr,igr,igr,vscnv,prs,1
     .                         ,reverse=.true.)

            !Sources
            if (.not.adiabatic) then
              !Joule heating
              joule = eeta(i,j,k)
     .               *scalarProduct_ijk(i,j,k,igr,jcnv(i,j,k,:)
     .                                           ,jcov(i,j,k,:))

              !Ion viscosity
              heat_src = joule
     .                  +viscous_heat_src(i,j,k,nx,ny,nz,igr,igr,igr
     .                                   ,nuu,vcnv,vcnv)

              !Electron viscosity
              !Comment out for old TM Hall MHD performance results
              if (di > 0d0.and.(.not.slava)) then
c$$$                if (fake_ve) then
c$$$                  heat_src = heat_src
c$$$     .                   +viscous_heat_src(i,j,k,nx,ny,nz,igr,igr,igr
c$$$     .                                    ,h_eta,vefcnv,vefcnv)
c$$$                else
                  heat_src = heat_src
     .                   +viscous_heat_src(i,j,k,nx,ny,nz,igr,igr,igr
     .                                    ,h_eta,vecnv ,vecnv)
c$$$                endif
              endif

              if (heat_src < 0d0) then
                write (*,*) 'Negative heat source'
                stop
              endif
            else
              heat_src = 0d0
            endif

            src(i,j,k,ITMP) = (gamma-1d0)*(advc+(heat_flx-heat_src))*jac

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

              !Stress tensor: div(Pi+Pe)
              if (slava) then
                src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ)
     .                             + div_pi(i,j,k,:)
              else
                src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ)
     .                             + div_pi(i,j,k,:) + div_pe(i,j,k,:)
              endif

              !Advection
              cnv = tensor_nc_advc(i,j,k,nx,ny,nz,igr,igr,igr,vcnv,vcnv
     .                            ,v_advect)

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ)/rho(i,j,k) + cnv !1/rho factor because
                                                                       !we evolve dv/dt

              if (alt__eom()) src(i,j,k,IVY) = src(i,j,k,IVY)*jac

            else
              !JxB
              if (nc_eom_jxb) then
                cnv = EOM_jxb(i,j,k,nx,ny,nz,igr,igr,igr
cc     .                       ,bcnv,jcnv,.true.
     .                       ,bcov,jcov,.false.
     .                       ,nc_eom_jxb)

                if (alt__eom()) cnv(2) = cnv(2)*jac
              else
                vec1 => bcnv
                vec2 => bcnv
                cnv = -dtnsr_src_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                              ,eom_jxb_x,eom_jxb_y,eom_jxb_z)
                nullify(vec1,vec2)
              endif

              src(i,j,k,IVX:IVZ) = -cnv

              !Grad(p)
              if (nc_eom_gp) then
                cnv=EOM_gp(i,j,k,nx,ny,nz,igr,igr,igr,rho,tmp,nc_eom_gp)

                if (alt__eom()) cnv(2) = cnv(2)*jac
              else
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_src_ijk(i,j,k,nx,ny,nz,igr,alt__eom()
     .                          ,eom_grad_p_x,eom_grad_p_y,eom_grad_p_z)
                nullify(sc1,sc2)
              endif

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) + cnv

              !Stress tensor: div(Pi+Pe) (should be made conservative)
              if (slava) then
                cnv = div_pi(i,j,k,:)
              else
                cnv = div_pi(i,j,k,:) + div_pe(i,j,k,:)
              endif
              if (alt__eom()) cnv(2) = cnv(2)*jac

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
