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

      use variables

      use timeStepping

      use nlfunction_setup

      use equilibrium

      implicit none

c Call variables

      real(8) :: ff(neqd)

      integer :: i,j,k

      type (var_array) :: varray

c Local variables

      integer(4) :: ig,jg,kg,ip,im,jp,jm,kp,km

      real(8)    :: dvol,ivol,dS1,dS2,dS3,dxx,dyy,dzz,xhm,yhm,zhm

      real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,dummy

      ! rho equation
      real(8)    :: advc,diffus

      ! Faraday's law
      real(8)    :: Ez_jp,Ez_jm,Ey_kp,Ey_km,Ex_kp,Ex_km
     .             ,Ez_ip,Ez_im,Ey_ip,Ey_im,Ex_jp,Ex_jm

      ! EOM
      real(8)    :: t11p,t12p,t13p,t11m,t12m,t13m,t11o,t12o,t13o
     .             ,t21p,t22p,t23p,t21m,t22m,t23m,t21o,t22o,t23o
     .             ,t31p,t32p,t33p,t31m,t32m,t33m,t31o,t32o,t33o

      ! tmp equation
      real(8)    :: heat_flx,heat_src,joule,viscous

      integer(4) :: ieq
      real(8)    :: cnv(3),cov(3),kappa

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

      dvol = dxx*dyy*dzz

      gsub   = gmetric%grid(igx)%gsub(i,j,k,:,:)
      gsuper = gmetric%grid(igx)%gsup(i,j,k,:,:)
      jac    = gmetric%grid(igx)%jac (i,j,k)

      jacip  = gmetric%grid(igx)%jac(ip,j,k)
      jacim  = gmetric%grid(igx)%jac(im,j,k)
      jacjp  = gmetric%grid(igx)%jac(i,jp,k)
      jacjm  = gmetric%grid(igx)%jac(i,jm,k)
      jackp  = gmetric%grid(igx)%jac(i,j,kp)
      jackm  = gmetric%grid(igx)%jac(i,j,km)

      ivol = 1d0/volume(i,j,k,igx,igy,igz)

c     Rho

      if (solve_rho) then
cc        advc =c_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,rho,sp=bcSP())
cc        advc =c_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,rho)
cc        advc =c_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,rho,sp=bcSP()
cc     .               ,upwind=.true.)

        advc =flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,rho,4)
cc     .                 ,sp=bcSP())

        if (dd /= 0d0) then
          diffus = dd*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,rho)
        else
          diffus = 0d0
        endif

        ff(IRHO) = advc - diffus
      else
        ff(IRHO) = 0d0
      endif

c     Faraday's law

      vec1 => vecnv
      vec2 => bcnv

      if (solenoidal) then
        ff(IBX:IBZ)= div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .                         ,btensor_x,btensor_y,btensor_z)
     .           + jac*dvol*curl(i,j,k,nx,ny,nz,igx,igy,igz,ejx,ejy,ejz)
      else
        ff(IBX:IBZ)= div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .                         ,btensor_x,btensor_y,btensor_z)
     .             - eeta(i,j,k)*veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz
     .                                       ,bcnv,.false.)
     .           + jac*dvol*curl(i,j,k,nx,ny,nz,igx,igy,igz,ejx,ejy,ejz)

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

c     Temperature

      if (gamma == 1d0) then
        ff(ITMP) = 0d0
      else
        advc =
     .        flx_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,tmp,4)
cc     .        c_advec(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,tmp)
     .        +(gamma-2.)/ivol*tmp(i,j,k)
     .                        *div(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz)

        !Heat flux
        if (chi /= 0d0) then
          heat_flx =-chi*laplacian(i,j,k,nx,ny,nz,igx,igy,igz,tmp)
        else
          heat_flx = 0d0
        endif

        !Heat sources
cc        if (gamma > 1d0) then
cc
cc          !!!Joule heating
cc          joule = dvol*eeta(i,j,k)*( jx(i,j,k)*jx_cov(i,j,k)
cc       .                            +jy(i,j,k)*jy_cov(i,j,k)
cc       .                            +jz(i,j,k)*jz_cov(i,j,k) )
cc
cc          !!!Viscous heating
cc          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)
cc
cc          cov_tnsr = matmul(nabla_v,gsub   )
cc          cnv_tnsr = matmul(gsuper ,nabla_v)
cc
cc          viscous = dvol*nuu(i,j,k)*sum(cov_tnsr*cnv_tnsr)/jac
cc
cc          heat_src = joule + viscous
cc
cc          if (heat_src < 0d0) then
cc            write (*,*) 'Heat source is negative'
cc            write (*,*) 'Aborting...'
cc            stop
cc          endif
cc
cc        else

          heat_src = 0d0

cc        endif

        !Construct temperature equation
        ff(ITMP) = advc+(gamma-1d0)*(heat_flx - heat_src)/rho(i,j,k)/a_p
      endif

c     EOM

      if (k_si > 0d0 .and. max_dv_dt /= 0d0) then
        ff(IVX:IVZ) = ff(IVX:IVZ)
     .              + si_op(i,j,k,nx,ny,nz,igx,igy,igz)/ivol
      else
        ff(IVX:IVZ) = 0d0
      endif

      !Lorentz force
      if (nc_eom_jxb) then
        cnv(1) = jy_cov(i,j,k)*bz_cov(i,j,k)
     .         - jz_cov(i,j,k)*by_cov(i,j,k)

        cnv(2) = jz_cov(i,j,k)*bx_cov(i,j,k)
     .         - jx_cov(i,j,k)*bz_cov(i,j,k)

        cnv(3) = jx_cov(i,j,k)*by_cov(i,j,k)
     .         - jy_cov(i,j,k)*bx_cov(i,j,k)

cc        cov(1) =(jy(i,j,k)*bz(i,j,k)
cc     .          -jz(i,j,k)*by(i,j,k))/jac
cc
cc        cov(2) =(jz(i,j,k)*bx(i,j,k)
cc     .          -jx(i,j,k)*bz(i,j,k))/jac
cc
cc        cov(3) =(jx(i,j,k)*by(i,j,k)
cc     .          -jy(i,j,k)*bx(i,j,k))/jac
cc
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                              ,cov(1),cov(2),cov(3)
cc     .                              ,cnv(1),cnv(2),cnv(3),.true.)

        ff(IVX:IVZ) = ff(IVX:IVZ) - cnv/ivol
      else
        ff(IVX:IVZ) = ff(IVX:IVZ)
     .               +div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                          ,eom_jxb_x,eom_jxb_y,eom_jxb_z)
      endif

      !Pressure force
      if (nc_eom_gp) then
        cov(1) = 0.5*( rho(i,j,k)*(tmp(ip,j,k)-tmp(im,j,k))/dxx
     .                +tmp(i,j,k)*(rho(ip,j,k)-rho(im,j,k))/dxx)

        cov(2) = 0.5*( rho(i,j,k)*(tmp(i,jp,k)-tmp(i,jm,k))/dyy
     .                +tmp(i,j,k)*(rho(i,jp,k)-rho(i,jm,k))/dyy)
        
        cov(3) = 0.5*( rho(i,j,k)*(tmp(i,j,kp)-tmp(i,j,km))/dzz
     .                +tmp(i,j,k)*(rho(i,j,kp)-rho(i,j,km))/dzz)

cc        cov(1)=0.5*(rho(ip,j,k)*tmp(ip,j,k)-rho(im,j,k)*tmp(im,j,k))/dxx
cc        cov(2)=0.5*(rho(i,jp,k)*tmp(i,jp,k)-rho(i,jm,k)*tmp(i,jm,k))/dyy
cc        cov(3)=0.5*(rho(i,j,kp)*tmp(i,j,kp)-rho(i,j,km)*tmp(i,j,km))/dzz

        cov = cov*a_p    !Multiply by alpha_p=1 + Ti/Te

        call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                              ,cov(1),cov(2),cov(3)
     .                              ,cnv(1),cnv(2),cnv(3),.true.)

        ff(IVX:IVZ) = ff(IVX:IVZ) + cnv/ivol
      else
        ff(IVX:IVZ) = ff(IVX:IVZ)
     .               +div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                          ,eom_grad_p_x,eom_grad_p_y,eom_grad_p_z)
      endif

      !Advective part
      if (nc_eom_v) then
        nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz,vx,vy,vz,0)

        do ieq=1,3
          cnv(ieq) =( vx(i,j,k)*nabla_v(1,ieq)
     .               +vy(i,j,k)*nabla_v(2,ieq)
     .               +vz(i,j,k)*nabla_v(3,ieq))/jac
        enddo

        cnv = cnv - nuu(i,j,k)
     .             *veclaplacian(i,j,k,nx,ny,nz,igx,igy,igz,vcnv
     .                          ,alt_eom,vol=.false.)

        ff(IVX:IVZ) = ff(IVX:IVZ)/rho(i,j,k) + cnv/ivol
      else
        ff(IVX:IVZ) = ff(IVX:IVZ)
     .               +div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                          ,eom_advc_x,eom_advc_y,eom_advc_z)
      endif

c Divide by cell volume factor

      ff = ff*ivol

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
