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

      integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,hex,hey,hez
      integer    :: ieq,nx,ny,nz

      real(8)    :: cnv(3),cov(3),jac

      ! rho equation
      real(8)    :: advc,diffus

      ! tmp equation
      real(8)    :: joule,viscous,heat_flx,heat_src

c Begin program

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      alt_eom = alt__eom()

c Grid parameters

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      jac = gmetric%grid(igx)%jac (i,j,k)

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

c       Advection + sources

        if (solve_prs) then  !Pressure equation

cc          advc=flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                  ,vscnv,tmp,advect,vol=.false.)
cc     .        +(gamma-1.)*tmp(i,j,k)
cc     .                   *div(i,j,k,nx,ny,nz,igx,igy,igz,vscnv)

          advc=   gamma *flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                            ,vscnv,tmp,advect)
cc     .       +(gamma-1.)*flx_nc_advec(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                       ,-vscnv,tmp,1,vol=.false.)
     .       +(1.-gamma)*dot_product(vscnv(i,j,k,:)
     .                  ,grad(i,j,k,nx,ny,nz,igx,igy,igz,tmp
     .                       ,hex=int(sign(1d0,vscnv(i,j,k,1)))
     .                       ,hey=int(sign(1d0,vscnv(i,j,k,2)))
     .                       ,hez=int(sign(1d0,vscnv(i,j,k,3)))))
     .                  /gmetric%grid(igx)%jac(i,j,k)

          ff(ITMP) = advc + (gamma-1d0)*heat_flx/a_p

        else                 !Temperature equation

cc          advc=flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                  ,vscnv,tmp,advect,vol=.false.)
cc     .        +(gamma-2.)*tmp(i,j,k)
cc     .                   *div(i,j,k,nx,ny,nz,igx,igy,igz,vscnv)

          advc=(gamma-1.)*flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                             ,vscnv,tmp,advect)
     .        -(gamma-2.)*flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                             ,vscnv,tmp,advect
     .                             ,conserv=.false.
     .                             ,zip_vel=.false.)
cc     .        +(2.-gamma)*dot_product(vscnv(i,j,k,:)
cc     .                  ,grad(i,j,k,nx,ny,nz,igx,igy,igz,tmp
cc     .                       ,hex=int(sign(1d0,-vscnv(i,j,k,1)))
cc     .                       ,hey=int(sign(1d0,-vscnv(i,j,k,2)))
cc     .                       ,hez=int(sign(1d0,-vscnv(i,j,k,3)))))
cc     .                  /gmetric%grid(igx)%jac(i,j,k)

          ff(ITMP) = advc + (gamma-1.)*heat_flx/a_p/rho(i,j,k)
        endif

        if (.not.adiabatic) then
          !Joule heating
          joule = eeta(i,j,k)
     .           *dot_product(jcnv(i,j,k,:),jcov(i,j,k,:))/jac

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

          if (heat_src < 0d0) then
            write (*,*) 'Negative heat source'
            stop
          endif

          !Add contributions
          if (solve_prs) then
            ff(ITMP) = ff(ITMP) - (gamma-1.)*heat_src/a_p
          else
            ff(ITMP) = ff(ITMP) - (gamma-1.)*heat_src/rho(i,j,k)/a_p
          endif
        endif

      endif

c EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<

c     JxB volume force

      if (nc_eom_jxb) then
        cnv = EOM_jxb(i,j,k,nx,ny,nz,igx,igy,igz
cc     .               ,bcnv,jcnv,.true.
     .               ,bcov,jcov,.false.
     .               ,nc_eom_jxb,alt_eom)
      else
        cnv = EOM_jxb(i,j,k,nx,ny,nz,igx,igy,igz,bcnv,bcnv,.true.
     .               ,nc_eom_jxb,alt_eom)
      endif

      ff(IVX:IVZ) = -cnv

c     grad(p) volume force

      cnv = EOM_gp(i,j,k,nx,ny,nz,igx,igy,igz,rho,tmp
     .            ,nc_eom_gp,alt_eom)

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
        cnv = tensor_advc(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,vcnv
     .                   ,v_advect,cons=.not.nc_eom_v)

        ff(IVX:IVZ) = ff(IVX:IVZ)/rho(i,j,k) + cnv  !1/rho factor because
                                                    !we evolve dv/dt
      else
        cnv = tensor_advc(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,pcnv
     .                   ,v_advect,cons=.not.nc_eom_v)

        ff(IVX:IVZ) = ff(IVX:IVZ) + cnv
      endif

c End

      end subroutine nonlinearRHS

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
