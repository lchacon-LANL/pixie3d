c nonlinearRHS
c#################################################################
      subroutine nonlinearRHS(i,j,k,varray,ff)

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

#if defined(petsc) && !defined(petsc_c)
      use nk_mod
#endif

      implicit none

c Call variables

      real(8) :: ff(neqd)

      integer :: i,j,k

      type(var_array),pointer :: varray

c Local variables

      integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,hex,hey,hez

      !Marder B cleaning
      real(8)    :: dS1,dS2,dS3,dxx,dyy,dzz
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      ! rho equation
      real(8)    :: advc,diffus

      ! tmp equation
      real(8)    :: heat_flx,heat_src,joule,viscous

      integer    :: ieq
      real(8)    :: cnv(3),cov(3),kappa,vlap1(3),vlap2(3)

c Begin program

      ip = i+1
      im = i-1
      jp = j+1
      jm = j-1
      kp = k+1
      km = k-1

      alt_eom = alt__eom()
cc      alt_eom = .false.

c     Grid parameters

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      dxx = dxh(ig)
      dyy = dyh(jg)
      dzz = dzh(kg)

      dS1 = dyy*dzz
      dS2 = dxx*dzz
      dS3 = dxx*dyy

      gsub   = gmetric%grid(igx)%gsub(i,j,k,:,:)
      gsuper = gmetric%grid(igx)%gsup(i,j,k,:,:)
      jac    = gmetric%grid(igx)%jac (i,j,k)

      jacip  = gmetric%grid(igx)%jac(ip,j,k)
      jacim  = gmetric%grid(igx)%jac(im,j,k)
      jacjp  = gmetric%grid(igx)%jac(i,jp,k)
      jacjm  = gmetric%grid(igx)%jac(i,jm,k)
      jackp  = gmetric%grid(igx)%jac(i,j,kp)
      jackm  = gmetric%grid(igx)%jac(i,j,km)

      if (.not.adiabatic .or. nc_eom_v)
     .   nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                  ,vcnv(:,:,:,1)
     .                  ,vcnv(:,:,:,2)
     .                  ,vcnv(:,:,:,3),0)

c     Rho

      if (solve_rho) then
        advc =flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,vcnv(:,:,:,1)
     .                 ,vcnv(:,:,:,2)
     .                 ,vcnv(:,:,:,3),rho,advect,sp=bcSP(),vol=.false.)

        if (dd /= 0d0) then
          diffus = dd*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,rho
     .                         ,vol=.false.)
        else
          diffus = 0d0
        endif

        ff(IRHO) = advc - diffus
      else
cc        !If this is added, need to add stuff in defineTSparameters below
cc        if (adiabatic) then
cc          ff(IRHO) = rho(i,j,k)
cc     .             - u_0%array_var(IRHO)%array(i,j,k)
cc     .             *(tmp(i,j,k)
cc     .              /u_0%array_var(ITMP)%array(i,j,k))**(1./gamma)
cc        else
          ff(IRHO) = 0d0
cc        endif
      endif

c     Faraday's law

#if !defined(vec_pot)
      vec1 => vecnv
      vec2 => bcnv

      ff(IBX:IBZ)= div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .                       ,btensor_x,btensor_y,btensor_z
     .                       ,vol=.false.)
     .            -curl(i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,E_ni(:,:,:,1),E_ni(:,:,:,2),E_ni(:,:,:,3))

      if (heta > 0d0 .and. (.not.solenoidal) ) then
        vlap2=veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,vlap
     .                     ,.false.,vol=.false.)

        ff(IBX:IBZ) = ff(IBX:IBZ) + di*heta*vlap2
      endif

cc      if (.not.solenoidal) then
cc      ff(IBX:IBZ)= div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
cc     .                         ,btensor_x,btensor_y,btensor_z
cc     .                         ,vol=.false.)
cc     .             - eeta(i,j,k)*veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                                       ,bcnv,.false.,vol=.false.)
cc     .             - curl(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                   ,E_ni(:,:,:,1)
cc     .                   ,E_ni(:,:,:,2)
cc     .                   ,E_ni(:,:,:,3))
cccc        !Marder divergence cleaning
cccc
cccccc        flxip = div(i ,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=1)
cccccc        flxim = div(im,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=1)
cccccc
cccccc        flxjp = div(i,j ,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=2)
cccccc        flxjm = div(i,jm,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=2)
cccccc
cccccc        flxkp = div(i,j,k ,nx,ny,nz,igx,igy,igz,bx,by,bz,he=3) 
cccccc        flxkm = div(i,j,km,nx,ny,nz,igx,igy,igz,bx,by,bz,he=3) 
cccc
cccc        flxip = 0.5*(div(ip,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
cccc     .              +div(i ,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
cccc        flxim = 0.5*(div(im,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
cccc     .              +div(i ,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
cccc
cccc        flxjp = 0.5*(div(i,jp,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
cccc     .              +div(i,j ,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
cccc        flxjm = 0.5*(div(i,jm,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
cccc     .              +div(i,j ,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
cccc
cccc        flxkp = 0.5*(div(i,j,kp,nx,ny,nz,igx,igy,igz,bx,by,bz) 
cccc     .              +div(i,j,k ,nx,ny,nz,igx,igy,igz,bx,by,bz))
cccc        flxkm = 0.5*(div(i,j,km,nx,ny,nz,igx,igy,igz,bx,by,bz) 
cccc     .              +div(i,j,k ,nx,ny,nz,igx,igy,igz,bx,by,bz))
cccc
cccc        !Take homogeneous dirichlet BCs at boundaries
cccc        if (i == nx) flxip = 0d0
cccc        if (i == 1 ) flxim = 0d0
cccc        if (j == ny) flxjp = 0d0
cccc        if (j == 1 ) flxjm = 0d0
cccc        if (k == nz) flxkp = 0d0
cccc        if (k == 1 ) flxkm = 0d0
cccc
cccc        kappa = 10*eeta(i,j,k)
cccc
cccc        ff(IBX:IBZ) = ff(IBX:IBZ)
cccc     .              - kappa*gsuper(:,1)*jac*dS1*(flxip-flxim)
cccc     .              - kappa*gsuper(:,2)*jac*dS2*(flxjp-flxjm)
cccc     .              - kappa*gsuper(:,3)*jac*dS3*(flxkp-flxkm)
cc      endif


      nullify(vec1,vec2)

#else

c     Vector potential

      cov(1) = (vecnv(i,j,k,2)*bcnv(i,j,k,3)
     .         -vecnv(i,j,k,3)*bcnv(i,j,k,2))/jac + E_ni(i,j,k,1)

      cov(2) = (vecnv(i,j,k,3)*bcnv(i,j,k,1)
     .         -vecnv(i,j,k,1)*bcnv(i,j,k,3))/jac + E_ni(i,j,k,2)

      cov(3) = (vecnv(i,j,k,1)*bcnv(i,j,k,2)
     .         -vecnv(i,j,k,2)*bcnv(i,j,k,1))/jac + E_ni(i,j,k,3)

      vlap1  = veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,acov,vol=.false.)
cc      vlap1  = vlap_a(i,j,k,:)
cc      cnv   = veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,acnv
cc     .                     ,.false.,vol=.false.)
cc
cc      call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,vlap1(1),vlap1(2),vlap1(3)
cc     .                        ,cnv (1),cnv (2),cnv (3)
cc     .                        ,.false.)

      if (heta > 0d0) then
        vlap2=veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,vlap_a,vol=.false.)
      else
        vlap2=0d0
      endif

      ff(IAX:IAZ) = -cov - eeta(i,j,k)*vlap1 + di*heta*vlap2

#endif

c     Temperature

      if (gamma == 1d0) then

        ff(ITMP) = 0d0

      else

c       Heat flux

        heat_flx = 0d0

        if (chi /= 0d0) then
          heat_flx =-chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,tmp
     .                            ,vol=.false.)
        endif

        if (chi_par /= 0d0) then
          heat_flx = heat_flx
     .              -(chi_par-chi)
     .              *par_flux_2nd(i,j,k,nx,ny,nz,igx,igy,igz,tmp
     .                           ,bhat(:,:,:,1)
     .                           ,bhat(:,:,:,2)
     .                           ,bhat(:,:,:,3),vol=.false.)
        endif

c       Advection + sources

        if (adiabatic) then  !Pressure equation
          advc=flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                  ,vcnv(:,:,:,1)
     .                  ,vcnv(:,:,:,2)
     .                  ,vcnv(:,:,:,3),tmp,advect,vol=.false.)
     .        +(gamma-1.)*tmp(i,j,k)
     .                   *div(i,j,k,nx,ny,nz,igx,igy,igz
     .                       ,vcnv(:,:,:,1)
     .                       ,vcnv(:,:,:,2)
     .                       ,vcnv(:,:,:,3))

          ff(ITMP) = advc + (gamma-1.)*heat_flx
        else                 !Temperature equation
          advc=flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                  ,vcnv(:,:,:,1)
     .                  ,vcnv(:,:,:,2)
     .                  ,vcnv(:,:,:,3),tmp,advect,vol=.false.)
     .        +(gamma-2.)*tmp(i,j,k)
     .                   *div(i,j,k,nx,ny,nz,igx,igy,igz
     .                       ,vcnv(:,:,:,1)
     .                       ,vcnv(:,:,:,2)
     .                       ,vcnv(:,:,:,3))

          !Joule heating
          joule = eeta(i,j,k)*( jcnv(i,j,k,1)*(jcov  (i,j,k,1)
     .                                        -jcov_0(i,j,k,1))
     .                         +jcnv(i,j,k,2)*(jcov  (i,j,k,2)
     .                                        -jcov_0(i,j,k,2))
     .                         +jcnv(i,j,k,3)*(jcov  (i,j,k,3) 
     .                                        -jcov_0(i,j,k,3)))/jac

          !Viscous heating
          cov_tnsr = matmul(nabla_v,gsub)

          if (sym_st) cov_tnsr = cov_tnsr + transpose(cov_tnsr)

          cnv_tnsr = matmul(gsuper,nabla_v)

          viscous = nuu(i,j,k)*sum(cov_tnsr*cnv_tnsr)/jac**2

          heat_src = joule + viscous

          if (viscous < 0d0) then
            call pstop('nonlinearRHS','Viscous heat source is negative')
          endif

          !Add contributions
          ff(ITMP) = advc
     .             +(gamma-1d0)*(heat_flx - heat_src)/rho(i,j,k)/a_p

        endif

      endif

c     EOM

      !SI contribution
      if (k_si > 0d0 .and. max_dv_dt /= 0d0) then
        ff(IVX:IVZ) = ff(IVX:IVZ)
     .              + si_op(i,j,k,nx,ny,nz,igx,igy,igz)
      else
        ff(IVX:IVZ) = 0d0
      endif

      !Volume force
      if (nc_eom_jxb) then
        cnv = eom_f(i,j,k,nx,ny,nz,igx,igy,igz,rho,tmp
     .           ,bcov,jcov,.false.,nc_eom_jxb,nc_eom_gp,alt_eom
cc        cnv = eom_f(i,j,k,nx,ny,nz,igx,igy,igz,rho,tmp
cc     .           ,bcnv,jcnv,.true.,nc_eom_jxb,nc_eom_gp,alt_eom
     $           ,vol=.false.)
      else
        cnv = eom_f(i,j,k,nx,ny,nz,igx,igy,igz,rho,tmp
     .           ,bcnv,bcnv,.true. ,nc_eom_jxb,nc_eom_gp,alt_eom
     $           ,vol=.false.)
      endif

      ff(IVX:IVZ) = ff(IVX:IVZ) - cnv

      !Advection
      if (nc_eom_v) then
        if (advect == 1) then !Upwinding
          hex = 1
          hey = 1
          hez = 1
          if (vcnv(i,j,k,1) > 0d0) hex = -1
          if (vcnv(i,j,k,2) > 0d0) hey = -1
          if (vcnv(i,j,k,3) > 0d0) hez = -1

          nabla_v = fnabla_v_upwd(i,j,k,nx,ny,nz,igx,igy,igz
     .                           ,vcnv(:,:,:,1)
     .                           ,vcnv(:,:,:,2)
     .                           ,vcnv(:,:,:,3),hex,hey,hez)
        endif

        do ieq=1,3
          cnv(ieq) =( vcnv(i,j,k,1)*nabla_v(1,ieq)
     .               +vcnv(i,j,k,2)*nabla_v(2,ieq)
     .               +vcnv(i,j,k,3)*nabla_v(3,ieq))/jac
        enddo

        ff(IVX:IVZ) = ff(IVX:IVZ)/rho(i,j,k) + cnv  !1/rho factor because
                                                    !we evolve dv/dt
cc        vec1 => vcnv
cc        vec2 => vcnv
cc        ff(IVX:IVZ) = ff(IVX:IVZ)/rho(i,j,k)
cc     .              + div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                          ,eom_advc_x,eom_advc_y,eom_advc_z
cc     .                          ,vol=.false.)
cc     .              - vcnv(i,j,k,:)*div(i,j,k,nx,ny,nz,igx,igy,igz
cc     .                                 ,vcnv(:,:,:,1)
cc     .                                 ,vcnv(:,:,:,2)
cc     .                                 ,vcnv(:,:,:,3))
cc        nullify(vec1,vec2)

      else
        vec1 => vcnv
        vec2 => pcnv
        ff(IVX:IVZ) = ff(IVX:IVZ)
     .              + div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                          ,eom_advc_x,eom_advc_y,eom_advc_z
     .                          ,vol=.false.)
        nullify(vec1,vec2)
      endif

      !Diffusion
      if (nu /= 0d0) then
        vec1 => vcnv
        coef => nuu
        cnv = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                  ,eom_diff_x,eom_diff_y,eom_diff_z
     .                  ,vol=.false.)
        nullify(vec1,coef)

#if !defined(RFX)
        if (nc_eom_v) cnv = cnv/rho(i,j,k)
#endif
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

cc      if (adiabatic) then
cc        one_over_dt(IRHO) = 0d0
cc        cnf(IRHO) = 0d0
cc      endif

c End program

      end subroutine defineTSParameters
