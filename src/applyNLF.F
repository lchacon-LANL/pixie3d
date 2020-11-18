#if defined(flux_rhs)

c nonlinearRHS
c#################################################################
      subroutine nonlinearRHS(i,j,k,igx,igy,igz,varray,ff)

c-----------------------------------------------------------------
c     This function computes the nonlinear rhs for the 
c     3D Hall MHD problem. Nonlinear functions Fi should be
c     coded with the following convention:
c        dt Ui + Fi(Uj) = 0
c     where Ui are the dependent variables, given in varray.
c     The result of Fi(Uj) at node (i,j) for all equations neqd
c     is returned in the vector ff of dimension neqd.
c-----------------------------------------------------------------

      use var_setup

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

      nx = gv%gparams%nxv(igx)
      ny = gv%gparams%nyv(igy)
      nz = gv%gparams%nzv(igz)

      call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = gv%gparams%dxh(ig)
      dyy = gv%gparams%dyh(jg)
      dzz = gv%gparams%dzh(kg)

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
c     3D Hall MHD problem. Nonlinear functions Fi should be
c     coded with the following convention:
c        dt Ui + Fi(Uj) = 0
c     where Ui are the dependent variables, given in varray.
c     The result of Fi(Uj) at node (i,j) for all equations neqd
c     is returned in the vector ff of dimension neqd.
c-----------------------------------------------------------------

      use var_setup

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

      integer :: ieq,nx,ny,nz

      real(8) :: cnv(3),cov(3),jac

      ! rho equation
      real(8) :: advc,diffus

      ! tmp equation
      real(8) :: joule,viscous,heat_flx,heat_srce

      ! EOM equation
      real(8) :: MS(3)

c Begin program

      nx = gv%gparams%nxv(igx)
      ny = gv%gparams%nyv(igy)
      nz = gv%gparams%nzv(igz)

      jac = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

c Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

      if (solve_rho) then
        advc = rho_advc(i,j,k)

        if ((.not.source_eval).and.(dd /= 0d0)) then
          diffus =lap(gv%gparams,i,j,k,nx,ny,nz,igx,igy,igz,rho,dff=dff)
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
      ff(IBX:IBZ) = curl(gv%gparams,i,j,k,nx,ny,nz,igx,igy,igz,E_ni)

#else

      !Vector potential
      ff(IAX:IAZ) = E_ni(i,j,k,:)

#endif

c Total kinetic pressure equation <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c     Advection

      advc = prs_advc(i,j,k)

c     Heat flux

      if (source_eval.and.(.not.adiabatic)) then
        heat_flx = 0d0
      else
        heat_flx = heat_xport(i,j,k)
      endif
      
c     Sources

      if ((.not.source_eval).and.(.not.adiabatic)) then
        heat_srce = heat_src(i,j,k)
      else
        heat_srce = 0d0
      endif

c     Add contributions

      ff(IPRS) = (advc + (gamma-1d0)*(heat_flx-heat_srce))*jac

c EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<

c     JxB volume force

      if (nc_eom_jxb) then
        cnv = crossProduct(gv%gparams,i,j,k,igx
     .                    ,jcov(i,j,k,:),bcov(i,j,k,:),.false.)
      else
        cnv = EOM_jxb(gv%gparams,i,j,k,nx,ny,nz,igx,bcnv,bcnv)
      endif

      ff(IVX:IVZ) = -cnv

c     grad(p) volume force

      if (use_p_eom) then
        cnv = EOM_gp (gv%gparams,i,j,k,nx,ny,nz,igx,prs    ,nc_eom_gp)
      else
        cnv = EOM_gnT(gv%gparams,i,j,k,nx,ny,nz,igx,rho,tmp,nc_eom_gp)
      endif

      ff(IVX:IVZ) = ff(IVX:IVZ) + cnv

c     Stress tensor: div(Pi+Pe)

      if (no_eom_divpe) then
        ff(IVX:IVZ) = ff(IVX:IVZ) + div_pi(i,j,k,:)
      else
        ff(IVX:IVZ) = ff(IVX:IVZ) + div_pi(i,j,k,:) + div_pe(i,j,k,:)
      endif

c     Advection

      if (nc_eom_v) then
        ff(IVX:IVZ) = ff(IVX:IVZ)/rho(i,j,k)   !1/rho factor because
     .              + v_advc(i,j,k,:)          !we evolve dv/dt
      else
        ff(IVX:IVZ) = ff(IVX:IVZ) + v_advc(i,j,k,:)
      endif

c     Sources

      if (.not.source_eval) then
        MS = ext_mom_src(gv%gparams,i,j,k,igx)
      else
        MS = 0d0
      endif

      if (nc_eom_v) then
        ff(IVX:IVZ) = ff(IVX:IVZ) - MS/rho(i,j,k)
      else
        ff(IVX:IVZ) = ff(IVX:IVZ) - MS
      endif

      if (alt__eom()) ff(IVY) = ff(IVY)*jac

c End

      end subroutine nonlinearRHS

#endif

c app_TS_Parameters
c####################################################################
      subroutine app_TS_Parameters(dt,one_over_dt,cnf)
c--------------------------------------------------------------------
c     Calculates application-dependent time step quantities for
c     nonlinear residuals.
c--------------------------------------------------------------------

      use equilibrium

      implicit none

c Call variables

      real(8) :: dt,one_over_dt(neqd),cnf(neqd)

c Local variables

c Begin program

      one_over_dt = 1d0/dt

      !Setup enslavement of rho
      if (.not.solve_rho) then
        one_over_dt(IRHO) = 0d0
        cnf(IRHO) = 0d0
      endif

c End program

      end subroutine app_TS_Parameters

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

      use var_setup

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
      real(8) :: joule,viscous,heat_flx

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 0,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (solve_rho) then
              advc = scl_advc_flx_x(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,vcnv(:,:,:,1),rho,advect)!,sp=bcSP())

              if (dd /= 0d0) then
                diffus = dd*lap_flx_x_ijk(gv%gparams,i,j,k,igr,rho)
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
            if ((chi /= 0d0).and.(gamma > 1d0).and.(chi_par==0d0)) then
              if (use_p_diff) then
                heat_flx = -lap_flx_x_ijk(gv%gparams,i,j,k,igr,prs
     .                                   ,diff=cchi)
              else
                heat_flx = -lap_flx_x_ijk(gv%gparams,i,j,k,igr,tmp
     .                                   ,diff=cchi)
              endif
            else
              heat_flx = 0d0
            endif

            advc = scl_advc_flx_x(gv%gparams,i,j,k,nx,ny,nz,igr
     .                           ,vscnv(:,:,:,1),prs,advect)!,sp=bcSP())

            flx_x(i+1,j,k,IPRS) = gamma*advc + (gamma-1d0)*heat_flx

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (nc_eom_v) then
              !Advection
              if (v_advect > 2) then
                adv_meth = v_advect
                vec1 => vcnv
                vec2 => vcnv
                cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,alt__eom(),tnsr_advc_x)
                nullify(vec1,vec2)
              else
                cnv = 0d0
              endif

              flx_x(i+1,j,k,IVX:IVZ) = cnv
            else
              !JxB
              if (nc_eom_jxb) then
                cnv = 0d0
              else
                vec1 => bcnv
                vec2 => bcnv
                cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,alt__eom(),eom_jxb_x)
                nullify(vec1,vec2)
              endif

              flx_x(i+1,j,k,IVX:IVZ) = cnv

              !grad(P)
              if (nc_eom_gp) then
                cnv = 0d0
              else
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,alt__eom(),eom_grad_p_x)
                nullify(sc1,sc2)
              endif

              flx_x(i+1,j,k,IVX:IVZ) = flx_x(i+1,j,k,IVX:IVZ) + cnv

              !Advection
              adv_meth = v_advect
              vec1 => vcnv
              vec2 => pcnv
              cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                           ,alt__eom(),tnsr_advc_x)
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

      use var_setup

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
      real(8) :: joule,viscous,heat_flx

c Begin program

      do k = 1,nz
        do j = 0,ny
          do i = 1,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (solve_rho) then
              advc = scl_advc_flx_y(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,vcnv(:,:,:,2),rho,advect)

              if (dd /= 0d0) then
                diffus = dd*lap_flx_y_ijk(gv%gparams,i,j,k,igr,rho)
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
            if ((chi /= 0d0).and.(gamma > 1d0).and.(chi_par==0d0)) then
              if (use_p_diff) then
                heat_flx = -lap_flx_y_ijk(gv%gparams,i,j,k,igr,prs
     .                                   ,diff=cchi)
              else
                heat_flx = -lap_flx_y_ijk(gv%gparams,i,j,k,igr,tmp
     .                                   ,diff=cchi)
              endif
            else
              heat_flx = 0d0
            endif

            advc = scl_advc_flx_y(gv%gparams,i,j,k,nx,ny,nz,igr
     .                           ,vscnv(:,:,:,2),prs,advect)

            flx_y(i,j+1,k,IPRS) = gamma*advc + (gamma-1d0)*heat_flx

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (nc_eom_v) then
              if (v_advect > 2) then
                adv_meth = v_advect
                vec1 => vcnv
                vec2 => vcnv
                cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,alt__eom(),tnsr_advc_y)
                nullify(vec1,vec2)
              else
                cnv = 0d0
              endif

              flx_y(i,j+1,k,IVX:IVZ) = cnv
            else
              !JxB
              if (nc_eom_jxb) then
                cnv = 0d0
              else
                vec1 => bcnv
                vec2 => bcnv
                cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,alt__eom(),eom_jxb_y)
                nullify(vec1,vec2)
              endif

              flx_y(i,j+1,k,IVX:IVZ) = cnv

              !grad(P)
              if (nc_eom_gp) then
                cnv = 0d0
              else
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,alt__eom(),eom_grad_p_y)
                nullify(sc1,sc2)
              endif

              flx_y(i,j+1,k,IVX:IVZ) = flx_y(i,j+1,k,IVX:IVZ) + cnv

              !Advection
              adv_meth = v_advect
              vec1 => vcnv
              vec2 => pcnv
              cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                           ,alt__eom(),tnsr_advc_y)
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

      use var_setup

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
      real(8) :: joule,viscous,heat_flx

c Begin program

      do k = 0,nz
        do j = 1,ny
          do i = 1,nx

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (solve_rho) then
              advc = scl_advc_flx_z(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,vcnv(:,:,:,3),rho,advect)

              if (dd /= 0d0) then
                diffus = dd*lap_flx_z_ijk(gv%gparams,i,j,k,igr,rho)
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
            if ((chi /= 0d0).and.(gamma > 1d0).and.(chi_par==0d0)) then
              if (use_p_diff) then
                heat_flx = -lap_flx_z_ijk(gv%gparams,i,j,k,igr,prs
     .                                   ,diff=cchi)
              else
                heat_flx = -lap_flx_z_ijk(gv%gparams,i,j,k,igr,tmp
     .                                   ,diff=cchi)
              endif
            else
              heat_flx = 0d0
            endif

            advc = scl_advc_flx_z(gv%gparams,i,j,k,nx,ny,nz,igr
     .                           ,vscnv(:,:,:,3),prs,advect)

            flx_z(i,j,k+1,IPRS) = gamma*advc + (gamma-1d0)*heat_flx

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (nc_eom_v) then
              if (v_advect > 2) then
                adv_meth = v_advect
                vec1 => vcnv
                vec2 => vcnv
                cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,alt__eom(),tnsr_advc_z)
                nullify(vec1,vec2)
              else
                cnv = 0d0
              endif

              flx_z(i,j,k+1,IVX:IVZ) = cnv
            else
              !JxB
              if (nc_eom_jxb) then
                cnv = 0d0
              else
                vec1 => bcnv
                vec2 => bcnv
                cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,alt__eom(),eom_jxb_z)
                nullify(vec1,vec2)
              endif

              flx_z(i,j,k+1,IVX:IVZ) = cnv

              !grad(P)
              if (nc_eom_gp) then
                cnv = 0d0
              else
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                             ,alt__eom(),eom_grad_p_z)
                nullify(sc1,sc2)
              endif

              flx_z(i,j,k+1,IVX:IVZ) = flx_z(i,j,k+1,IVX:IVZ) + cnv

              !Advection
              adv_meth = v_advect
              vec1 => vcnv
              vec2 => pcnv
              cnv = dtnsr_flx_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                           ,alt__eom(),tnsr_advc_z)
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

      use var_setup

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
      real(8) :: joule,viscous,heat_srce,heat_flx

      ! EOM equation
      real(8) :: MS(3)

c Begin program

      do k = 1,nz
        do j = 1,ny
          do i = 1,nx

            jac = gv%gparams%gmetric%grid(igr)%jac(i,j,k)

            !Rho <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            src(i,j,k,IRHO) = 0d0

            !Faraday's law <<<<<<<<<<<<<<<<<<<<<<
#if !defined(vec_pot)
            !Magnetic field
            src(i,j,k,IBX:IBZ)
     .           = curl(gv%gparams,i,j,k,nx,ny,nz,igr,igr,igr,E_ni)
#else
            !Vector potential
            src(i,j,k,IAX:IAZ) = E_ni(i,j,k,:)
#endif

            !Temperature <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

            !Heat Flux (should be made completely conservative)
            heat_flx = heat_xport(i,j,k)

            !Heat sources
            if (gamma > 1d0) then
              !Advection (non-conservative piece, with reversed sign for velocity)
              advc = flx_nc_advec(gv%gparams,i,j,k,nx,ny,nz,igr,igr,igr
     .                           ,vscnv,prs,1,reverse=.true.)

              !Sources
              heat_srce = heat_src(i,j,k)

              src(i,j,k,IPRS)=(gamma-1d0)*(advc+heat_flx-heat_srce)*jac

            else

              src(i,j,k,IPRS) = 0d0

            endif

            !EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<
            if (nc_eom_v) then
              !JxB force
              if (nc_eom_jxb) then
                cnv = crossProduct(gv%gparams,i,j,k,igr
     .                    ,jcov(i,j,k,:),bcov(i,j,k,:),.false.)
              else
                cnv = EOM_jxb(gv%gparams,i,j,k,nx,ny,nz,igr,bcnv,bcnv)
              endif

              src(i,j,k,IVX:IVZ) = -cnv

              !grad(p) volume force
              cnv = EOM_gnT(gv%gparams,i,j,k,nx,ny,nz,igr,rho,tmp
     .                     ,nc_eom_gp)

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) + cnv

              !Stress tensor: div(Pi+Pe)
              if (no_eom_divpe) then
                src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ)
     .                             + div_pi(i,j,k,:)
              else
                src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ)
     .                             + div_pi(i,j,k,:) + div_pe(i,j,k,:)
              endif


              !Advection
              if (v_advect <= 2) then
                cnv= tensor_nc_advc(gv%gparams,i,j,k,nx,ny,nz
     .                             ,igr,igr,igr,vcnv,vcnv,v_advect)
              else
                cnv= -vcnv(i,j,k,:)
     .               *div(gv%gparams,i,j,k,nx,ny,nz,igr,igr,igr,vcnv)
              endif

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ)/rho(i,j,k) + cnv !1/rho factor because
                                                                       !we evolve dv/dt

              !Momentum source
              MS = ext_mom_src(gv%gparams,i,j,k,igr)

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) - MS/rho(i,j,k)

              if (alt__eom()) src(i,j,k,IVY) = src(i,j,k,IVY)*jac

            else
              !JxB
              if (nc_eom_jxb) then
                cnv = crossProduct(gv%gparams,i,j,k,igr
     .                    ,jcov(i,j,k,:),bcov(i,j,k,:),.false.)
                if (alt__eom()) cnv(2) = cnv(2)*jac
              else
                vec1 => bcnv
                vec2 => bcnv
                cnv = -dtnsr_src_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                       ,alt__eom(),eom_jxb_x,eom_jxb_y,eom_jxb_z)
                nullify(vec1,vec2)
              endif

              src(i,j,k,IVX:IVZ) = -cnv

              !Grad(p)
              if (nc_eom_gp) then
                cnv=EOM_gnT(gv%gparams,i,j,k,nx,ny,nz,igr,rho,tmp
     .                     ,nc_eom_gp)

                if (alt__eom()) cnv(2) = cnv(2)*jac
              else
                sc1 => rho
                sc2 => tmp
                cnv = dtnsr_src_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .               ,alt__eom()
     .               ,eom_grad_nT_x,eom_grad_nT_y,eom_grad_nT_z)
                nullify(sc1,sc2)
              endif

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) + cnv

              !Stress tensor: div(Pi+Pe) (should be made conservative)
              if (no_eom_divpe) then
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
              cnv = dtnsr_src_ijk(gv%gparams,i,j,k,nx,ny,nz,igr
     .                  ,alt__eom(),tnsr_advc_x,tnsr_advc_y,tnsr_advc_z)
              nullify(vec1,vec2)

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) + cnv

              !Momentum source
              MS = ext_mom_src(gv%gparams,i,j,k,igr)

              src(i,j,k,IVX:IVZ) = src(i,j,k,IVX:IVZ) - MS

            endif
          enddo
        enddo
      enddo

c End program

      end function set_nlf_src
