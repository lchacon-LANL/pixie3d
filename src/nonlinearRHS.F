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

      implicit none

c Call variables

      real(8) :: ff(neqd)

      integer :: i,j,k

      type(var_array),pointer :: varray

c Local variables

      integer    :: ig,jg,kg,ip,im,jp,jm,kp,km

      !Marder B cleaning
      real(8)    :: dS1,dS2,dS3,dxx,dyy,dzz
      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

      ! rho equation
      real(8)    :: advc,diffus

      ! tmp equation
      real(8)    :: heat_flx,heat_src,joule,viscous

      integer    :: ieq
      real(8)    :: cnv(3),cov(3),kappa,vlap(3),vlap2(3)

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

c     Rho

      if (solve_rho) then
        !Various flavors of advection (old)
cc        advc =c_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,rho,sp=bcSP())
cc        advc =c_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,rho)
cc        advc =c_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,rho,sp=bcSP()
cc     .               ,upwind=.true.)

        !General advection (w/ SP treatment for particle conservation)
        advc =flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,vcnv(:,:,:,1)
     .                 ,vcnv(:,:,:,2)
     .                 ,vcnv(:,:,:,3),rho,advect,sp=bcSP(),vol=.false.)
        !Hardwired upwinding (generalizes ZIP)
cc        advc =flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,rho,1
cc     .                 ,sp=bcSP(),vol=.false.)
        !Hardwired SMART (generalizes ZIP)
cc        advc =flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,rho,4
c                        ,vol=.false.)

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

c     Faraday's law

#if !defined(vec_pot)
      vec1 => vecnv
      vec2 => bcnv

      if (solenoidal) then
        ff(IBX:IBZ)= div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .                         ,btensor_x,btensor_y,btensor_z
     .                         ,vol=.false.)
     .             + curl(i,j,k,nx,ny,nz,igx,igy,igz,grd(:,:,:,1)
     .                   ,grd(:,:,:,2),grd(:,:,:,3))
      else
        ff(IBX:IBZ)= div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .                         ,btensor_x,btensor_y,btensor_z
     .                         ,vol=.false.)
     .             - eeta(i,j,k)*veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                                       ,bcnv,.false.,vol=.false.)
     .             - curl(i,j,k,nx,ny,nz,igx,igy,igz,grd(:,:,:,1)
     .                   ,grd(:,:,:,2),grd(:,:,:,3))
cc        !Marder divergence cleaning
cc
cccc        flxip = div(i ,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=1)
cccc        flxim = div(im,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=1)
cccc
cccc        flxjp = div(i,j ,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=2)
cccc        flxjm = div(i,jm,k,nx,ny,nz,igx,igy,igz,bx,by,bz,he=2)
cccc
cccc        flxkp = div(i,j,k ,nx,ny,nz,igx,igy,igz,bx,by,bz,he=3) 
cccc        flxkm = div(i,j,km,nx,ny,nz,igx,igy,igz,bx,by,bz,he=3) 
cc
cc        flxip = 0.5*(div(ip,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
cc     .              +div(i ,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
cc        flxim = 0.5*(div(im,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
cc     .              +div(i ,j,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
cc
cc        flxjp = 0.5*(div(i,jp,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
cc     .              +div(i,j ,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
cc        flxjm = 0.5*(div(i,jm,k,nx,ny,nz,igx,igy,igz,bx,by,bz)
cc     .              +div(i,j ,k,nx,ny,nz,igx,igy,igz,bx,by,bz))
cc
cc        flxkp = 0.5*(div(i,j,kp,nx,ny,nz,igx,igy,igz,bx,by,bz) 
cc     .              +div(i,j,k ,nx,ny,nz,igx,igy,igz,bx,by,bz))
cc        flxkm = 0.5*(div(i,j,km,nx,ny,nz,igx,igy,igz,bx,by,bz) 
cc     .              +div(i,j,k ,nx,ny,nz,igx,igy,igz,bx,by,bz))
cc
cc        !Take homogeneous dirichlet BCs at boundaries
cc        if (i == nx) flxip = 0d0
cc        if (i == 1 ) flxim = 0d0
cc        if (j == ny) flxjp = 0d0
cc        if (j == 1 ) flxjm = 0d0
cc        if (k == nz) flxkp = 0d0
cc        if (k == 1 ) flxkm = 0d0
cc
cc        kappa = 10*eeta(i,j,k)
cc
cc        ff(IBX:IBZ) = ff(IBX:IBZ)
cc     .              - kappa*gsuper(:,1)*jac*dS1*(flxip-flxim)
cc     .              - kappa*gsuper(:,2)*jac*dS2*(flxjp-flxjm)
cc     .              - kappa*gsuper(:,3)*jac*dS3*(flxkp-flxkm)
      endif

      nullify(vec1,vec2)

#else

c     Vector potential

      cov(1) = (vecnv(i,j,k,2)*bcnv(i,j,k,3)
     .         -vecnv(i,j,k,3)*bcnv(i,j,k,2))/jac + grd(i,j,k,1)

      cov(2) = (vecnv(i,j,k,3)*bcnv(i,j,k,1)
     .         -vecnv(i,j,k,1)*bcnv(i,j,k,3))/jac + grd(i,j,k,2)

      cov(3) = (vecnv(i,j,k,1)*bcnv(i,j,k,2)
     .         -vecnv(i,j,k,2)*bcnv(i,j,k,1))/jac + grd(i,j,k,3)

      vlap  = veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,acov,vol=.false.)
cc      vlap  = vlap_a(i,j,k,:)
cc      cnv   = veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,acnv
cc     .                     ,.false.,vol=.false.)
cc
cc      call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,vlap(1),vlap(2),vlap(3)
cc     .                        ,cnv (1),cnv (2),cnv (3)
cc     .                        ,.false.)

      if (heta > 0d0) then
        vlap2=veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz,vlap_a,vol=.false.)
      else
        vlap2=0d0
      endif

      ff(IAX:IAZ) = -(cov + eeta(i,j,k)*vlap - heta*vlap2)

#endif

c     Temperature

      if (gamma == 1d0) then

        ff(ITMP) = 0d0

      else

c       Advection

        advc= flx_advec(i,j,k,nx,ny,nz,igx,igy,igz
     .                 ,vcnv(:,:,:,1)
     .                 ,vcnv(:,:,:,2)
     .                 ,vcnv(:,:,:,3),tmp,advect,vol=.false.)
     .      +(gamma-2.)*tmp(i,j,k)
     .                 *div(i,j,k,nx,ny,nz,igx,igy,igz
     .                     ,vcnv(:,:,:,1)
     .                     ,vcnv(:,:,:,2)
     .                     ,vcnv(:,:,:,3))

c       Heat flux

        if (chi /= 0d0) then
          heat_flx =-chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,tmp
     .                            ,vol=.false.)
        else
          heat_flx = 0d0
        endif

c       Heat sources

        if (.not.adiabatic) then
          !Joule heating
          joule = eeta(i,j,k)*( jcnv(i,j,k,1)*jcov(i,j,k,1)
     .                         +jcnv(i,j,k,2)*jcov(i,j,k,2)
     .                         +jcnv(i,j,k,3)*jcov(i,j,k,3) )/jac

          !Viscous heating
          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,vcnv(:,:,:,1)
     .                      ,vcnv(:,:,:,2)
     .                      ,vcnv(:,:,:,3),0)

          cov_tnsr = matmul(nabla_v,gsub)

          if (sym_st) cov_tnsr = cov_tnsr + transpose(cov_tnsr)

          cnv_tnsr = matmul(gsuper,nabla_v)

          viscous = nuu(i,j,k)*sum(cov_tnsr*cnv_tnsr)/jac**2

          heat_src = joule + viscous

          if (heat_src < 0d0) then
            call pstop('nonlinearRHS','Heat source is negative')
          endif
        else
          heat_src = 0d0
        endif

c       Add contributions
        ff(ITMP) = advc+(gamma-1d0)*(heat_flx - heat_src)/rho(i,j,k)/a_p

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
     $           ,vol=.false.)
      else
        cnv = eom_f(i,j,k,nx,ny,nz,igx,igy,igz,rho,tmp
     .           ,bcnv,jcnv,.true. ,nc_eom_jxb,nc_eom_gp,alt_eom
     $           ,vol=.false.)
      endif

      ff(IVX:IVZ) = ff(IVX:IVZ) - cnv

      !Advection
      if (nc_eom_v) then
        do ieq=1,3
          cnv(ieq) =( vcnv(i,j,k,1)*nabla_v(1,ieq)
     .               +vcnv(i,j,k,2)*nabla_v(2,ieq)
     .               +vcnv(i,j,k,3)*nabla_v(3,ieq))/jac
        enddo

        cnv = cnv
     .      - veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,vcnv
     .                    ,alt_eom,nuu,vol=.false.)
     .       /rho(i,j,k)

        ff(IVX:IVZ) = ff(IVX:IVZ)/rho(i,j,k) + cnv  !1/rho factor because
                                                    !we evolve dv/dt
      else
        vec1 => vcnv
        vec2 => pcnv
        ff(IVX:IVZ) = ff(IVX:IVZ)
     .              + div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                          ,eom_advc_x,eom_advc_y,eom_advc_z
     .                          ,vol=.false.)
        nullify(vec1,vec2)
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

c End program

      end subroutine defineTSParameters
