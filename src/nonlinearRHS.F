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

      !Marder B cleaning
cc      real(8)    :: dS1,dS2,dS3,dxx,dyy,dzz,kappa
cc      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

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
        advc = flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                  ,vcnv,rho,advect,sp=bcSP(),vol=.false.)

        if (dd /= 0d0) then
          diffus = dd*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,rho
     .                         ,vol=.false.)
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

      if (gamma == 1d0) then

        ff(ITMP) = 0d0

      else

c       Heat flux

        heat_flx = 0d0

        if (chi /= 0d0) then
          heat_flx =-a_p  !Factor of a_p needed because variable should be total T=Ti+Te=a_p*Te
     .              *chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,tmp
cc     .                            ,rho,vol=.false.)  !chi_perp=chi*rho
     .                            ,vol=.false.)  !chi_perp=chi
        endif

        if (chi_par /= chi) then
          heat_flx = heat_flx
     .              -a_p  !Factor of a_p needed because variable should be total T=Ti+Te=a_p*Te
     .              *(chi_par-chi)
     .              *par_diff(i,j,k,nx,ny,nz,igx,igy,igz,tmp
     .                   ,bhat(:,:,:,1),bhat(:,:,:,2),bhat(:,:,:,3),2
     .                   ,symmetric=.true.,limited=.false.,vol=.false.)
        endif

c       Advection + sources

        if (adiabatic) then  !Pressure equation

          advc=flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                  ,vscnv,tmp,advect,vol=.false.)
     .        +(gamma-1.)*tmp(i,j,k)
     .                   *div(i,j,k,nx,ny,nz,igx,igy,igz,vscnv)

          ff(ITMP) = advc + (gamma-1d0)*heat_flx/a_p

        else                 !Temperature equation

          advc=flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                  ,vscnv,tmp,advect,vol=.false.)
     .        +(gamma-2.)*tmp(i,j,k)
     .                   *div(i,j,k,nx,ny,nz,igx,igy,igz,vscnv)

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
          ff(ITMP) = advc
     .             +(gamma-1d0)*(heat_flx - heat_src)/rho(i,j,k)/a_p

        endif

      endif

c EOM <<<<<<<<<<<<<<<<<<<<<<<<<<<

c     Volume force: (jxB-grad(p))

      if (nc_eom_jxb) then
        cnv = EOM_f(i,j,k,nx,ny,nz,igx,igy,igz,rho,tmp
     .             ,bcov,jcov,.false.,nc_eom_jxb,nc_eom_gp,alt_eom
cc     .             ,bcnv,jcnv,.true.,nc_eom_jxb,nc_eom_gp,alt_eom
     $             ,vol=.false.)
      else
        cnv = EOM_f(i,j,k,nx,ny,nz,igx,igy,igz,rho,tmp
     .             ,bcnv,bcnv,.true. ,nc_eom_jxb,nc_eom_gp,alt_eom
     $             ,vol=.false.)
      endif

      ff(IVX:IVZ) = -cnv

c     Stress tensor: div(Pi)

      if (slava) then
        cnv = div_pi(i,j,k,:) !Match Slava's EOM
      else
        cnv = div_pi(i,j,k,:) + div_pe(i,j,k,:)
      endif

      ff(IVX:IVZ) = ff(IVX:IVZ) + cnv

c     Advection

      if (nc_eom_v) then
        cnv = tensor_advc(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,vcnv
     .                   ,advect,cons=.not.nc_eom_v)

        ff(IVX:IVZ) = ff(IVX:IVZ)/rho(i,j,k) + cnv  !1/rho factor because
                                                    !we evolve dv/dt
      else
        cnv = tensor_advc(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,pcnv
cc     .                   ,2,cons=.not.nc_eom_v)
     .                   ,advect,cons=.not.nc_eom_v)
cc     .                   ,min(advect,2),cons=.not.nc_eom_v)

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

      one_over_dt = 1d0/dt

      !Setup enslavement of rho
      if (.not.solve_rho) then
        one_over_dt(IRHO) = 0d0
        cnf(IRHO) = 0d0
      endif

c End program

      end subroutine defineTSParameters
