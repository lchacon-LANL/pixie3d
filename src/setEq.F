c setEquilibrium
c####################################################################
      subroutine setEquilibrium(iigx,varray)

c--------------------------------------------------------------------
c     Set initial conditions for initial value calculation, define
c     boundary conditions, and label physical quantities.
c
c     Boundary conditions specification in bcs is defined in routine
c     applyBoundaryConditions.f. 
c
c     Equilibrium quantities are:
c       * ieq=IRHO: density
c       * ieq=IVX:IVZ: momentum if conservative, velocity otherwise
c       * ieq=IBX:IBZ: magnetic field
c       * ieq=IAX:IAZ: vector potential
c       * ieq=ITMP: total pressure is found as p=n*(Ti+Te)=a_p*n*Te.
c--------------------------------------------------------------------

      use parameters

      use equilibrium

      use grid

      use timeStepping

      use nlfunction_setup

      use error

      use oned_int

      use local_BCS_variables

      use app_iosetup

      use B_tools, ONLY:vlap_div_clean

      implicit none

c Call variables

      integer :: iigx

      type(var_array),pointer :: varray

c Local variables

      integer :: i,j,k,ig,jg,kg,ieq,nmax,iglobal,jglobal
     .          ,ierr,sh,it,icomp,bcnd(6,3),brank
     .          ,igx,igy,igz,nx,ny,nz,igr,nxx,nyy,nzz

      real(8) :: ldaflow,x1,y1,z1,cx,cy,cz,x00

      logical :: covariant,to_cartsn,cartsn,a_bc

      real(8) :: r,jac1,bnorm,aaa,bbb,ccc,qq,qqp,q0,RR,ff,nn,bb
     .          ,aspect_ratio,mm,kk,Iz,bz0,btheta,vr,rrho,r1
     .          ,r0(1),aa,dummy1,dummy2,dummy3,eta_in,err,bz_avg,dr
     .          ,bzh,bth,rh,e0_over_eta,m0_over_nu,theta_ppnch,p0,eps
     .          ,phi,rho_norm,Rtor,sigma,gsub(3,3),bb0(3)
     .          ,hlx,y_char,pw

      real(8) :: rint(0:nxdp),dh,offset(1),bth0(0:nxdp)
     .          ,bzz0(0:nxdp),vr0(0:nxdp),rho0(0:nxdp),vz0(0:nxdp)
     .          ,eth0(0:nxdp),ez0(0:nxdp),jth0(0:nxdp)
     .          ,jz0 (0:nxdp),bth0_old(0:nxdp),p00(0:nxdp)
     .          ,bzz0_old(0:nxdp),dummy(0:nxdp),bmag2(0:nxdp)
     .          ,q00(0:nxdp),dum (ilom:ihip,2)
     .          ,ldum(ilom:ihip,2)
     .          ,r2  (ilom:ihip,jlom:jhip)
     .          ,br02d (0:nxdp,0:nydp),bth02d(0:nxdp,0:nydp)
     .          ,bzz02d(0:nxdp,0:nydp)
     .          ,vth02d(0:nxdp,0:nydp),vr02d(0:nxdp,0:nydp)
     .          ,psi0(0:nxdp)

      real(8) :: pp(ilom:ihip,jlom:jhip)

      real(8) :: a (ilom:ihip,jlom:jhip,klom:khip,3)

      real(8) :: var(ilom:ihip,jlom:jhip,klom:khip,neqd)

      character(20),pointer,dimension(:) :: label

      integer :: bcs(6,neqd)

#if defined(vec_pot)
      real(8),allocatable,dimension(:,:,:,:) :: ag,bg
#endif
      real(8),allocatable,dimension(:,:,:,:) :: vmec_car

      INTERFACE
        subroutine setAppBCs(varray,vaux)
          use var_setup
          type(var_array),pointer :: varray
          type(aux_array),pointer :: vaux
        end subroutine setAppBCs
      END INTERFACE

c Functions

      real(8) :: qprof,pprof
      external   qprof,pprof

c Begin program

      ldaflow = dlambda/rshear

      E0 = 0d0 ; B0 = 0d0

      label => varray%array_var%descr

      do ieq=1,neqd
        var(:,:,:,ieq) = varray%array_var(ieq)%array
cc        var(:,:,:,ieq:ieq) => varray%array_var(ieq:ieq)%array
      enddo

c Initialize required local grid information

      igx = iigx
      igy = iigx
      igz = iigx

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

c Label variables
     
      label(IRHO) = 'J*Rho'
#if !defined(vec_pot)
      label(IBX)  = 'B^1'
      label(IBY)  = 'B^2'
      label(IBZ)  = 'B^3'
#else
      label(IAX)  = 'A_1'
      label(IAY)  = 'A_2'
      label(IAZ)  = 'A_3'
#endif
      label(IVX)  = 'P^1'
      if (alt__eom()) then
        label(IVY)  = 'J*P^2'
      else
        label(IVY)  = 'P^2'
      endif
      label(IVZ)  = 'P^3'
      label(ITMP) = 'J*(pi+pe)'

c Define boundary conditions for dependent variables

      call defineDepVariables(varray)

c Set initial guess

      select case (trim(equil))

c     --------------------------------------------
c     ---------------  WAVE TESTS ----------------
c     --------------------------------------------

      case ('msw') !Magnetosonic wave

        gamma = 1d0

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

#if !defined(vec_pot)
        var(:,:,:,IBX)  = 0d0
        var(:,:,:,IBY)  = 0d0
        var(:,:,:,IBZ)  = 1d0
#else
        B0(3) = 1d0

        var(:,:,:,IAX:IAZ) = 0d0
#endif

      case ('3wv') !3-wave MHD test

        gamma = 1.6667

        phi   = grid_params%params(1)   !Angle of B with x-axis, in degrees
        beta  = grid_params%params(2)   !plasma beta

        phi   = phi/180.*pi             !Phi in radians

#if !defined(PER_BC_SYNC)
        if (di > 0d0) then
          call pstop('setEquilibrium','Compile with PER_BC_SYNC=t')
        endif
#endif

c     Check input

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (test
     .     .and.(   (adiabatic)
     .          .or.(a_p /= 1d0)
     .          .or.(.not.solve_rho)
cc     .          .or.(.not.use_p_eom)
     .          )
     .     ) then
          messg = 'Wrong input for this equilibrium. '
cc     .          //'Check adiab, Ti/Te, solve_rho, use_p_eom=F,0.0,T,T'
     .          //'Check adiab, Ti/Te, solve_rho=F,0.0,T'
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 0.5*beta/a_p  !Electron temperature

        B0(1) = cos(phi)
        B0(2) = sin(phi)
        B0(3) = 0d0

#if !defined(vec_pot)
        var(:,:,:,IBX)  = B0(1)
        var(:,:,:,IBY)  = B0(2)
        var(:,:,:,IBZ)  = B0(3)
#else
        var(:,:,:,IAX:IAZ) = 0d0
#endif

      case ('mswsn') !Magnetosonic wave in sinusoidal geometry

        gamma = 1d0

c     Check coordinates

        if (coords /= 'sin') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0
        var(:,:,:,IRHO) = 1d0
        var(:,:,:,ITMP) = 1d0

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              jac1 = gmetric%grid(igx)%jac(i,j,k)

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = 0d0
              var(i,j,k,IBZ)  = jac1*1d0
            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')
cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc
cc              jac1 = gmetric%grid(igx)%jac(i,j,k)  !Not defined in the global domain
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = 0d0
cc              b(i,j,k,3) = jac1*1d0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('whslr') !Whistler wave

        gamma = 1d0

c     Check input

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (test
     .     .and.(   (a_p /= 2d0)
     .          .or.(.not.solve_rho))) then
          messg = 'Wrong input. Check Ti/Te, solve_rho =F,1.0,T'
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

#if !defined(vec_pot)
        var(:,:,:,IBX)  = 1d0
        var(:,:,:,IBY)  = 0d0
        var(:,:,:,IBZ)  = 0d0
#else
        B0(1) = 1d0

        var(:,:,:,IAX:IAZ) = 0d0
#endif

      case ('kaw') !KAW wave

        gamma = 1d0

        beta = 1d4

        di = (xmax-xmin)/2/pi/nh1*sqrt(beta/(beta-2))  !di k_par = 1

c     Check input

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (test
     .     .and.(   (adiabatic)
     .          .or.(a_p /= 1d0)
     .          .or.(.not.solve_rho)
cc     .          .or.(use_p_eom)
     .          )
     .     ) then
cc        if ((adiabatic).or.(a_p /= 1d0).or. (.not.solve_rho)) then
          messg = 'Wrong input for this equilibrium. '
cc     .          //'Check adiab,Ti/Te,solve_rho,use_p_eom=F,0.,T,F'
     .          //'Check adiab,Ti/Te,solve_rho=F,0.,T'
          call pstop('setEquilibrium',messg)
        endif

c     Uniform medium with constant magnetic field in cartesian coordinates

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 0.5*beta/a_p

#if !defined(vec_pot)
        var(:,:,:,IBX)  = 1d0
        var(:,:,:,IBY)  = 0d0
        var(:,:,:,IBZ)  = 0d0
#else
        B0(1) = 1d0

        var(:,:,:,IAX:IAZ) = 0d0
#endif

c     --------------------------------------------
c     ------------  ADVECTION TESTS --------------
c     --------------------------------------------

      case ('adv1') !Advection test (Germaschewski)

        gamma = 2d0

        eps   = grid_params%params(1)
        sigma = grid_params%params(2)
        x00   = grid_params%params(3)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Opposite Gaussian perturbations in rho, T, uniform V, no B

        var(:,:,:,IVX)  = 1d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

c$$$              var(i,j,k,IRHO) = 1d0 + eps*exp(-0.5*(x1-0.75)**2/sigma**2)
c$$$              var(i,j,k,ITMP) = 1d0 - eps*exp(-0.5*(x1-0.75)**2/sigma**2)
              var(i,j,k,IRHO) = 1d0 +eps*exp(-0.5*(x1-x00)**2/sigma**2)
              var(i,j,k,ITMP) = 1d0 -eps*exp(-0.5*(x1-x00)**2/sigma**2)

            enddo
          enddo
        enddo

#if !defined(vec_pot)
        var(:,:,:,IBX:IBZ) = 0d0
#else
        var(:,:,:,IAX:IAZ) = 0d0
#endif

      case ('adv2') !Advection test (Jasak et al, IJNMF 31, 1999; Lenard, CMAME 88, 1991)

        gamma = 2d0

        eps   = grid_params%params(1)
        sigma = grid_params%params(2)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Opposite Gaussian perturbations in rho, T, uniform V, no B

        var(:,:,:,IVX)  = 0.866
        var(:,:,:,IVY)  = 0.5
        var(:,:,:,IVZ)  = 0d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              y_char = 1./6.+0.5/0.866*x1
              var(i,j,k,IRHO) = 1d0+eps*0.5*(1d0+sign(1d0,y1-y_char))
              var(i,j,k,ITMP) = 1d0-eps*0.5*(1d0+sign(1d0,y1-y_char))

            enddo
          enddo
        enddo

#if !defined(vec_pot)
        var(:,:,:,IBX:IBZ) = 0d0
#else
        var(:,:,:,IAX:IAZ) = 0d0
#endif

      case ('sod') !Hydrodynamic Sod problem

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Opposite Gaussian perturbations in rho, T, uniform V, no B

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)

              var(i,j,k,IRHO) = 1d0-0.5*(1d0+sign(1d0,x1-0.5))*0.875
              var(i,j,k,ITMP) = 1d0-0.5*(1d0+sign(1d0,x1-0.5))*0.9  !Pressure

            enddo
          enddo
        enddo

        var(:,:,:,ITMP) = var(:,:,:,ITMP)/var(:,:,:,IRHO)/a_p !Temperature

#if !defined(vec_pot)
        var(:,:,:,IBX:IBZ) = 0d0
#else
        var(:,:,:,IAX:IAZ) = 0d0
#endif

c     --------------------------------------------------
c     ---------------  KELVIN-HELMHOLTZ ----------------
c     --------------------------------------------------

      case ('khcar')

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (.not.solve_rho)then
          messg = 'Wrong input for this equilibrium. '
     .          //'Check solve_rho=T'
          call pstop('setEquilibrium',messg)
        endif

c     Kelvin-Helmholtz with constant magnetic field and pressure in cartesian coordinates

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)
              var(i,j,k,IVY)=vperflow
     .                      *tanh((x1-0.5d0*(xmax+xmin))/dlambda)
            enddo
          enddo
        enddo

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVZ)  = 0d0

#if !defined(vec_pot)
        var(:,:,:,IBX)  = 0d0
        var(:,:,:,IBY)  = 0d0
        var(:,:,:,IBZ)  = 1d0
#else
        call fillVectorPotential(igx,var(:,:,:,IAX:IAZ))
#endif

      case ('3dkh')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(igx,a)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

cc        if ((.not.adiabatic).or.(a_p /= 2d0).or. (.not.solve_rho)) then
cc          messg = 'Wrong input for this equilibrium. '
cc     .          //'Check solve_rho=T'
cc          call pstop('setEquilibrium',messg)
cc        endif

c     3D KHI w/ differential rotation (Knoll and Brackbill, PoP (2002))

        beta = 1d0

        var(:,:,:,ITMP) = 0.5*beta/a_p  !Electron temperature

#if !defined(vec_pot)
        var(:,:,:,IBX:IBZ) = curl(igx,a)
#else
        var(:,:,:,IAX:IAZ) = a
#endif

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)
#if !defined(vec_pot)
             bb0 = var(i,j,k,IBX:IBZ)
#else
             bb0 = curl(i,j,k,nx,ny,nz,igx,igy,igz,a)
#endif
             var(i,j,k,IRHO) = 2 - bb0(3)**2

             var(i,j,k,IVX) = vperflow*tanh((y1-.5*(ymax+ymin))/ldaflow)
     .                                *sin (2*pi*z1/(zmax-zmin))
             var(i,j,k,IVY) = 0d0
             var(i,j,k,IVZ) = 0d0

            enddo
          enddo
        enddo

      case ('khihk')

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (.not.solve_rho)then
          messg = 'Wrong input for this equilibrium. '
     .          //'Check solve_rho=T'
          call pstop('setEquilibrium',messg)
        endif

c     Kelvin-Helmholtz with constant magnetic field and pressure in cartesian coordinates

        var(:,:,:,IRHO) = 1d0

        beta = 0.3d0
        var(:,:,:,ITMP) = 0.5*beta/a_p

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                    ,x1,y1,z1)
              var(i,j,k,IVY)=vperflow
     .                      *tanh((x1-0.5d0*(xmax+xmin))/dlambda)
            enddo
          enddo
        enddo

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        B0(1) = 0d0
        B0(2) = 0.04998d0
        B0(3) = 0.998d0

#if !defined(vec_pot)
        var(:,:,:,IBX)  = B0(1)
        var(:,:,:,IBY)  = B0(2)
        var(:,:,:,IBZ)  = B0(3)
#else
        call fillVectorPotential(igx,var(:,:,:,IAX:IAZ))
#endif
c     ----------------------------------------------
c     ---------------  TEARING MODE ----------------
c     ----------------------------------------------

      case ('tmcar')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(igx,a)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

cc        if (test
cc     .     .and.(   (.not.adiabatic)
cc     .          .or.(a_p /= 2d0)
cc     .          .or.(.not.solve_rho)
cc     .          .or.(.not.use_p_eom))) then
cc          write (*,*) 'Ti/Te,solve_rho,adiabatic,use_p_eom='
cc     .               ,temp_ratio,solve_rho,adiabatic,use_p_eom
cc          messg = 'Wrong input. '
cc     .     //'Check Ti/Te,solve_rho,adiabatic,use_p_eom = 1.0,T,T,T'
cc          call pstop('setEquilibrium',messg)
cc        endif

c     Tearing mode in cartesian coordinates

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

        bz0 = 1d0

#if !defined(vec_pot)
        var(:,:,:,IBX:IBZ) = curl(igx,a)
        var(:,:,:,IBZ)=sqrt(bz0**2 - var(:,:,:,IBY)**2)
#else
        var(:,:,:,IAX:IAZ) = a
#endif

      case ('tmsin')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(igx,a)

c     Check coordinates

        if (coords /= 'sin') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Tearing mode in sinusoidal coordinates

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d0

        bz0 = 1d0

#if !defined(vec_pot)
        var(:,:,:,IBX:IBZ) = curl(igx,a)

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              !X-Y equilibrium
             gsub = gmetric%grid(igx)%gsub(i,j,k,:,:)
             jac1 = gmetric%grid(igx)%jac (i,j,k)

             bnorm= vectorNorm(i,j,k,igx,var(i,j,k,IBX:IBZ),.false.)

             ccc = jac1*(bz0**2 - bnorm)
             bbb = gsub(3,2)*var(i,j,k,IBY) + gsub(3,1)*var(i,j,k,IBX)
             aaa = gsub(3,3)

             var(i,j,k,IBZ)  = (-bbb+sqrt(bbb**2+4*aaa*ccc))/2./aaa

             if (bbb**2+4*aaa*ccc < 0d0) then
               write (*,*) var(i,j,k,IBZ),bnorm,aaa,bbb**2+4*aaa*ccc
             endif
            enddo
          enddo
        enddo
#else
        var(:,:,:,IAX:IAZ) = a
#endif

cc      case ('gem')
cc
ccc     Define vector potential (in curvilinear coordinates) for initialization
cc
cc        call fillVectorPotential(igx,a)
cc
ccc     Check coordinates
cc
cc        if (coords /= 'car' .and. coords /= 'scl') then
cc          messg = 'Wrong coordinates for equilibrium '//equil
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
cc        if ((a_p /= 6d0).or. (.not.solve_rho).or.(adiabatic)) then
cc          messg = 'Wrong input for this equilibrium. '
cc     .          //'Check Ti/Te, solve_rho, adiabatic = 5.0, T, F'
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
ccc     GEM challenge
cc
cc        var(:,:,:,IVX)  = 0d0
cc        var(:,:,:,IVY)  = 0d0
cc        var(:,:,:,IVZ)  = 0d0
cc
cc        beta = 1d0
cc
cc        var(:,:,:,ITMP) = 0.5*beta/a_p !Electron temperature required for force balance
cc
cc        do k = klo,khi
cc          do j = jlo,jhi
cc            do i = ilo,ihi
cc             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                   ,x1,y1,z1)
cc
cc             var(i,j,k,IRHO)=1d0/cosh((x1-0.5d0*(xmax+xmin))/dlambda)**2
cc     .                      +0.2
cc
cccc             var(i,j,k,IBX) = 0d0
cccc             var(i,j,k,IBY) = tanh((x1-0.5*(xmax+xmin))/dlambda)
cccc             var(i,j,k,IBZ) = 0d0
cc            enddo
cc          enddo
cc        enddo
cc
cc#if !defined(vec_pot)
cc        var(:,:,:,IBX:IBZ) = curl(igx,a)
cc#else
cc        var(:,:,:,IAX:IAZ) = a
cc#endif

      case ('gem')

c     Check coordinates

        if (coords /= 'car' .and. coords /= 'lpk') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (test
     .     .and.(   (adiabatic)
     .          .or.(a_p /= 6d0)
     .          .or.(.not.solve_rho)
cc     .          .or.(.not.sym_st)
     .          )
     .     ) then
          if (my_rank == 0) then
            write (*,*)
     .           'Check Ti/Te,solve_rho,adiabatic,sym_st=5.0, T, F, T'
            write (*,*) 'Currently:',a_p - 1,solve_rho,adiabatic,sym_st
          endif
          messg = 'Wrong input for this equilibrium. '
          call pstop('setEquilibrium',messg)
        endif

c     GEM challenge

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        beta = 1d0

        var(:,:,:,ITMP) = 0.5*beta/a_p !Electron temperature required for force balance

        if (bcond(1) /= FSYM .and. bcond(2) == FSYM) then
          hlx = xmax
        elseif (bcond(1) == FSYM .and. bcond(2) /= FSYM) then
          hlx = xmin
        else
          hlx = 0.5d0*(xmax+xmin)
        endif

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,IRHO)=1d0/cosh((x1-hlx)/dlambda)**2 + 0.2
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        call fillVectorPotential(igx,a)

        var(:,:,:,IBX:IBZ) = curl(igx,a)
#else
        call fillVectorPotential(igx,var(:,:,:,IAX:IAZ))
#endif

      case ('kai1')

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

cc        if ((.not.adiabatic).or.(a_p /= 2d0).or.(.not.solve_rho)) then
cc          messg = 'Wrong input for this equilibrium. '
cc     .          //'Check adiab, Ti/Te, solve_rho=T,1.0,T'
cc          call pstop('setEquilibrium',messg)
cc        endif

c     Kai TM test

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,ITMP) = 1d0 !Electron pressure

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,IRHO)=0.5/cosh((x1-0.5d0*(xmax+xmin))/dlambda)**2
     .                      +1d0
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        call fillVectorPotential(igx,a)

        var(:,:,:,IBX:IBZ) = curl(igx,a)
#else
        call fillVectorPotential(igx,var(:,:,:,IAX:IAZ))
#endif

      case ('dtm')

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(igx,a)

c     Check coordinates

        if (coords /= 'car' .and. coords /= 'lpk') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

cc        if ((a_p /= 2d0).or. (.not.solve_rho) .or. (adiabatic)) then
cc          messg = 'Wrong input for this equilibrium. '
cc     .          //'Check Ti/Te, solve_rho, adiabatic = 5.0, T, F'
cc          call pstop('setEquilibrium',messg)
cc        endif

c     Tearing mode in cartesian coordinates

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,ITMP) = 0.5d0  !Electron temperature

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,IRHO)=0.5/cosh((x1-xmin)/dlambda)**2
     .                      +0.5/cosh((x1-xmax)/dlambda)**2
     .                      +1d0
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        var(:,:,:,IBX:IBZ) = curl(igx,a)
#else
        var(:,:,:,IAX:IAZ) = a
#endif

c     -------------------------------------------------
c     ---------------  ISLAND COALESCENCE -------------
c     -------------------------------------------------

      case ('ic')

        p0     = grid_params%params(1)
        eps    = grid_params%params(2)
        dlambda = (ymax-ymin)/4./pi

c     Define vector potential (in curvilinear coordinates) for initialization

        call fillVectorPotential(igx,a)

c     Check coordinates

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Island coalescence in cartesian coordinates

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,ITMP)=(p0 + 0.5*(1-eps**2)
     $                       /(    cosh(x1/dlambda)
     $                        +eps*cos (y1/dlambda))**2)
     $                       /a_p/var(i,j,k,IRHO)
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        var(:,:,:,IBX:IBZ) = curl(igx,a)
#else
        var(:,:,:,IAX:IAZ) = a
#endif

      case ('ic-hk')  !Homa Karimabadi's equilibrium (use symmetry BCs)

        beta   = eq_params(1)
        eps    = eq_params(2)

cc        dlambda = 0.5/pi !(Ly/4/pi) with Ly=2 in full domain
        if (bcond(3) == FSYM .or. bcond(4) == FSYM) then
          dlambda = (ymax-ymin)/2./pi !(2*Ly/4/pi)
        else
          dlambda = (ymax-ymin)/4./pi !(Ly/4/pi)
        endif

c     Check coordinates

        if (coords /= 'car')then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Island coalescence in cartesian coordinates

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,ITMP) = 0.5*beta/a_p

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                   ,x1,y1,z1)

             var(i,j,k,IRHO)=(0.2
     .                      +(1-eps**2)/(     cosh(x1/dlambda)
     .                                   +eps*cos (y1/dlambda))**2 )
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        call fillVectorPotential(igx,a)
        var(:,:,:,IBX:IBZ) = curl(igx,a)

cc        do k = klo,khi
cc          do j = jlo,jhi
cc            do i = ilo,ihi
cc             call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                   ,x1,y1,z1)
cc
cc             var(i,j,k,IBX)=eps
cc     .                     *sin (y1/dlambda)/(     cosh(x1/dlambda)
cc     .                                        +eps*cos (y1/dlambda))
cc             var(i,j,k,IBY)=sinh(x1/dlambda)/(     cosh(x1/dlambda)
cc     .                                        +eps*cos (y1/dlambda))
cc             var(:,:,:,IBZ) = 0d0
cc            enddo
cc          enddo
cc        enddo
cc        call XformVector(igx,var(:,:,:,IBX:IBZ),'car','cnv')
#else
        call fillVectorPotential(igx,var(:,:,:,IAX:IAZ))
#endif

c     ------------------------------------------------
c     ---------------  3D CYLINDRICAL ----------------
c     ------------------------------------------------

      case ('3dkai')

cc        mm = grid_params%params(1)
cc        kk = grid_params%params(2)
        RR = grid_params%params(3)

c     Check coordinates

        if (coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (abs(RR-zmax/2./pi) > 1d-5) then
          messg = 'Major radius and zmax do not agree'
          call pstop('setEquilibrium',messg)
        endif

        if ((a_p /= 2d0).or. (.not.solve_rho)) then
          messg = 'Wrong input for this equilibrium. '
     .          //'Check Ti/Te, solve_rho=1.0,T'
          call pstop('setEquilibrium',messg)
        endif

c     Tokamak equilibria (analytical)

        Iz = 2.22
        bb = 0.5*Iz

cc        nh2 = mm       !To set the right perturbation wavelength
cc        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2      !For perturbations to be of the right order
                        !(nh2,nh3 set by input file)

        if (nh2 == 0 .and. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*x1/(1+(x1/dlambda)**2)
              bz0    = RR

#if !defined(vec_pot)
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta
              var(i,j,k,IBZ)  = x1*bz0
#endif
              var(i,j,k,ITMP) = 1d-5
     .                     + 1/16.*(Iz*dlambda/(1+(x1/dlambda)**2))**2

            enddo
          enddo
        enddo

#if defined(vec_pot)
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              btheta = bb*x1/(1+(x1/dlambda)**2)
cc              bz0    = RR
cc
cc              b(i,j,k,1)  = 0d0                  
cc              b(i,j,k,2)  = btheta
cc              b(i,j,k,3)  = x1*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

cc      case ('3dpar')  !EQ for parallel transport test
cc
cccc        mm = grid_params%params(1)
cccc        kk = grid_params%params(2)
cc        RR = grid_params%params(3)
cc        q0 = grid_params%params(4)
cc
ccc     Check coordinates
cc
cc        if (coords /= 'cyl') then
cc          messg = 'Wrong coordinates for equilibrium '//equil
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
cc        if (abs(RR-zmax/2./pi) > 1d-5) then
cc          messg = 'Major radius and zmax do not agree'
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
ccc     Simple screw pinch equilibria
cc
cccc        nh2 = mm       !To set the right perturbation wavelength
cccc        nh3 = kk*RR    !To set the right perturbation wavelength
cc
cc        nh1 = nh2      !For perturbations to be of the right order
cc                       !(nh2,nh3 set by input file)
cc
cc        if (nh2 == 0 .or. nh3 == 0) then
cc          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
cc     .          //' equilibrium not set'
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
cc        !Build equilibrium
cc
cc        var(:,:,:,IVX)  = 0d0
cc        var(:,:,:,IVY)  = 0d0
cc        var(:,:,:,IVZ)  = 0d0
cc
cc        var(:,:,:,IRHO) = 1d0
cc
cc        var(:,:,:,ITMP) = 1d-5
cc
cc#if !defined(vec_pot)
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)
cc
cccc              btheta = x1  !Green's function orbit test
cccc              bz0 = 1d0
cccc
cccc              btheta = 2*x1  !CosT orbit test
cccc              bz0 = 1d0
cc
cc              bz0 = 1d0
cc              bb  = q0*bz0*dlambda/RR
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc
cc              var(i,j,k,IBX) = 0d0
cc              var(i,j,k,IBY) = btheta
cc              var(i,j,k,IBZ) = x1*bz0
cc
cc            enddo
cc          enddo
cc        enddo
cc#else
cc        call pstop('setEquilibrium','Equilibrium not working')
cc
cccc        do k = 0,nzdp
cccc          do j = 0,nydp
cccc            do i = 0,nxdp
cccc              x1 = grid_params%xg(i)  !Global coordinate
cccc
cccc              if (coords == 'hl2') x1 = sign(sqrt(2*abs(x1)),x1)
cccc
cccc              bb     = (dlambda**2+aa)
cccc     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
cccc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cccc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cccc
cccc              if (coords == 'hel') then
cccc                b(i,j,k,1)  = 0d0
cccc                b(i,j,k,2)  = btheta + kk*x1/mm*bz0
cccc                b(i,j,k,3)  = x1*bz0
cccc              else
cccc                b(i,j,k,1)  = 0d0
cccc                b(i,j,k,2)  = btheta/x1 + kk/mm*bz0
cccc                b(i,j,k,3)  = bz0
cccc              endif
cccc            enddo
cccc          enddo
cccc        enddo
cccc
cccc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cccc     .                          ,b
cccc     .                          ,var(:,:,:,IAX:IAZ))
cc#endif

      case ('3drfp')

cc        mm = grid_params%params(1)
cc        kk = grid_params%params(2)
        RR = grid_params%params(3)
        aa = grid_params%params(4)

c     Check coordinates

        if (coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (abs(RR-zmax/2./pi) > 1d-5) then
          messg = 'Major radius and zmax do not agree'
          call pstop('setEquilibrium',messg)
        endif

c     RFP equilibria (analytical)

cc        nh2 = mm       !To set the right perturbation wavelength
cc        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2      !For radial perturbations to be of the right order
                       !(nh2,nh3 set by input file)

        if (nh2 == 0 .and. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d-5

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta
              var(i,j,k,IBZ)  = x1*bz0
            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              bb     = (dlambda**2+aa)
cc     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = btheta
cc              b(i,j,k,3) = x1*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

c     ----------------------------------------------------
c     ---------------  HELICAL EQUILIBRIA ----------------
c     ----------------------------------------------------

      case ('spnch')  !SCREW PINCH EQ.

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     RFP equilibria (analytical)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        aa = grid_params%params(4)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d-5

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              if (coords == 'hel') then
                var(i,j,k,IBX)  = 0d0
                var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
                var(i,j,k,IBZ)  = x1*bz0
              endif
            enddo
          enddo
        enddo
#else
        allocate(bg(0:nxdp,0:ny+1,0:nz+1,3))

        do k = 0,nz+1
          do j = 0,ny+1
            do i = 0,nxdp

              x1 = grid_params%xg(i)  !Global coordinate

              bb     = (dlambda**2+aa)
     .                 /sqrt((dlambda**2+aa)**2-dlambda**4)
              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              bg(i,j,k,1) = 0d0
              bg(i,j,k,2) = btheta + kk*x1/mm*bz0
              bg(i,j,k,3) = x1*bz0
            enddo
          enddo
        enddo

        call findVectorPotential(nxd,ny,nz,var(:,:,:,IVX:IVZ)
     .                          ,bg,var(:,:,:,IAX:IAZ))

        deallocate(bg)
#endif

      case ('2dtok')

c     Cylindrical Tokamak analytical equilibrium (Furth et al., 1973)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        q0 = grid_params%params(4)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        bb = dlambda/RR/q0

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,ITMP) = 1d-3

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))

              !X-Y equilibrium
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
              var(i,j,k,IBZ)  = x1*bz0
            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = btheta + kk*x1/mm*bz0
cc              b(i,j,k,3) = x1*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('2dkai')

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (test
     .     .and.(    (a_p /= 2d0)
     .          .or. (.not.adiabatic)
     .          .or. (.not.solve_rho)
cc     .          .or. (use_p_eom)
     .          )
     .     ) then
          messg = 'Wrong input. '
cc     .          //'Check adiab,Ti/Te,solve_rho,use_p_eom=T,1.0,T,F'
     .          //'Check adiab,Ti/Te,solve_rho=T,1.0,T'
          call pstop('setEquilibrium',messg)
        endif

c     Tokamak equilibria (analytical)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)

        Iz           = 2.22
        bb           = 0.5*Iz

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              btheta = bb*x1/(1+(x1/dlambda)**2)
              bz0    = RR

#if !defined(vec_pot)
              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = btheta + kk*x1/mm*bz0
              var(i,j,k,IBZ)  = x1*bz0
#endif
              var(i,j,k,ITMP) = 1d0
     .                    + 1./16.*(Iz*dlambda/(1+(x1/dlambda)**2))**2
            enddo
          enddo
        enddo

#if defined(vec_pot)
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              btheta = bb*x1/(1+(x1/dlambda)**2)
cc              bz0    = RR
cc
cc              b(i,j,k,1)  = 0d0                  
cc              b(i,j,k,2)  = btheta + kk*x1/mm*bz0
cc              b(i,j,k,3)  = x1*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('rfp83')

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     RFP EQUILIBRIA (Caramana et al, PoP, 1983)

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        !Find pressure and q-profiles
        bz0  = 1d0
        beta = 0d0

        do i=0,nxdp
          x1  = grid_params%xg(i)
          p00(i) = 2d-5
          q00(i) = RR*qprof(x1)
        enddo

        !Initialize toroidal field
        bzz0 = bz0  

        !Solve 1D Grad-Shafranov equation
        call GS_equil_1d(q00,p00,bth0,bzz0)

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

              var(i,j,k,IBX)  = 0d0
              var(i,j,k,IBY)  = bth0(iglobal) + kk*x1/mm*bzz0(iglobal)
              var(i,j,k,IBZ)  = x1*bzz0(iglobal)

              var(i,j,k,ITMP) = p00(iglobal)/a_p/var(i,j,k,IRHO)

            enddo
          enddo
        enddo
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = bth0(i) + kk*x1/mm*bzz0(i)
cc              b(i,j,k,3) = x1*bzz0(i)
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case ('ntm2d')

        bootstrap= .true.

        mm   = grid_params%params(1) !Poloidal mode number
        kk   = grid_params%params(2) !Helical paramter (-nn/RR, nn-> toroidal mode number)
        RR   = grid_params%params(3) !Major radius (determines aspect ratio)
        bz0  = grid_params%params(4) !Toroidal field at magnetic axis
        beta = grid_params%params(5) !Beta toroidal

        nh2 = mm  !To set the right perturbation wavelength
        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Equilibria (numerical)

        !Find pressure and q-profiles
        do i=0,nxdp
          x1  = grid_params%xg(i)
          p00(i) = 0.5*beta*bz0**2*pprof(x1,0)/pprof(0d0,0)
          q00(i) = RR*qprof(x1)
        enddo

        !Initialize toroidal field
        bzz0 = bz0  

        !Solve 1D Grad-Shafranov equation
        call GS_equil_1d(q00,p00,bth0,bzz0)

        !Build equilibrium

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

#if !defined(vec_pot)
              var(i,j,k,IBX) = 0d0
              var(i,j,k,IBY) = bth0(iglobal) + kk*x1/nh2*bzz0(iglobal)
              var(i,j,k,IBZ) = x1*bzz0(iglobal)
#endif

              var(i,j,k,ITMP) = 0.5*p00(iglobal)/var(i,j,k,IRHO)

            enddo
          enddo
        enddo

#if defined(vec_pot)
cc        allocate(bg(0:nxdp,0:nydp,0:nzdp,3))
cc
cc        do k = 0,nzdp
cc          do j = 0,nydp
cc            do i = 0,nxdp
cc              x1 = grid_params%xg(i)  !Global coordinate
cc
cc              bg(i,j,k,1) = 0d0
cc              bg(i,j,k,2) = bth0(i) + kk*x1/nh2*bzz0(i)
cc              bg(i,j,k,3) = x1*bzz0(i)
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,bg
cc     .                          ,var(:,:,:,IAX:IAZ))
cc
cc        deallocate(bg)
#endif

cc      case ('3dtok')
cc
ccc     Cylindrical Tokamak analytical equilibrium (Furth et al., 1973)
cc
cc        RR = grid_params%params(1)
cccc        mm = grid_params%params(4)
cccc        nn = grid_params%params(5)
cc        q0 = grid_params%params(6)
cc
cc        kk = 0d0   !To avoid affecting q-profile
cc
cc        bb = dlambda/RR/q0
cc
ccc     Check coordinates
cc
cc        if (coords /= 'tor') then
cc          messg = 'Need "tor" coordinates for equilibrium '//equil
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
ccc     Failsafes
cc
cccc        nh2 =  mm       !To set the right perturbation wavelength
cccc        nh3 = -nn       !To set the right perturbation wavelength
cc
cc        nh1 = nh2       !For perturbations to be of the right order
cc                        !(nh2, nh3 set in input file)
cc
cc        if (nh2 == 0 .or. nh3 == 0) then
cc          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
cc     .          //' equilibrium not set'
cc          call pstop('setEquilibrium',messg)
cc        endif
cc
ccc     Build equilibrium
cc
cc        var(:,:,:,IVX)  = 0d0
cc        var(:,:,:,IVY)  = 0d0
cc        var(:,:,:,IVZ)  = 0d0
cc
cc        var(:,:,:,IRHO) = 1d0
cc
cc        var(:,:,:,ITMP) = 1d-3
cc
cc#if !defined(vec_pot)
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)
cc
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,x1,y1,z1)
cc
cc              !X-Y equilibrium
cc              var(i,j,k,IBX)  = 0d0
cc              var(i,j,k,IBY)  = btheta
cc              var(i,j,k,IBZ)  = gmetric%grid(igx)%gsup(i,j,k,3,3)*bz0
cc            enddo
cc          enddo
cc        enddo
cc#else
cc        call pstop('setEquilibrium','Equilibrium not working')
cc
cccc        do k = klom,khip
cccc          do j = jlom,jhip
cccc            do i = ilom,ihip
cccc
cccc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cccc     .                                      ,x1,y1,z1)
cccc
cccc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cccc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cccc
cccc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cccc     .                                    ,x1,y1,z1)
cccc
cccc              !X-Y equilibrium
cccc              b(i,j,k,1) = 0d0
cccc              b(i,j,k,2) = btheta
cccc              b(i,j,k,3) = gmetric%grid(igx)%gsup(i,j,k,3,3)*bz0
cccc            enddo
cccc          enddo
cccc        enddo
cccc
cccc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cccc     .                          ,b
cccc     .                          ,var(:,:,:,IAX:IAZ))
cc#endif

#if defined(vmec)

c     --------------------------------------------------
c     ---------------  VMEC EQUILIBRIUM ----------------
c     --------------------------------------------------

      case ('vmec')    

        nh1 = nh2       !For perturbations to be of the right order in radius
                        !(nh2,nh3 set by input file)

        if (nh2 == 0 .and. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

c     Read VMEC coordinates and fill grid hierarchy

        if (coords /= 'ext') then
          messg = 'Need "ext" coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        call vmec_map(.false.)

c     Read VMEC equilibrium

#if !defined(vec_pot)
        call vmec_equ(igx,nx,ny,nz,var(:,:,:,IBX)
     .                            ,var(:,:,:,IBY)
     .                            ,var(:,:,:,IBZ)
     .                            ,var(:,:,:,ITMP)
     .                            ,var(:,:,:,IRHO),gamma)
#else
        call pstop('setEquilibrium','Equilibrium not working')

cc        call vmec_equ(igx,nx,ny,nz,b(:,:,:,1),b(:,:,:,2)
cc     .                            ,b(:,:,:,3),var(:,:,:,ITMP))
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

c     Other equilibrium quantities

        var(:,:,:,IVX:IVZ)  = 0d0

cc        var(:,:,:,IRHO) = 1d0

        !Transform to e temperature
        where (var(:,:,:,IRHO) /= 0d0)
          var(:,:,:,ITMP) = var(:,:,:,ITMP)/var(:,:,:,IRHO)/a_p
        elsewhere
          var(:,:,:,ITMP) = 0d0
        end where
#endif

c     --------------------------------------------------
c     ---------------  OHMIC EQUILIBRIA ----------------
c     --------------------------------------------------

      case ('ppnch','ppnsl','ppnst')

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        theta_ppnch = grid_params%params(4)
        e0_over_eta = grid_params%params(5)
        beta        = grid_params%params(6)
        m0_over_nu  = grid_params%params(7)
        pw = grid_params%params(8)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Equilibria (numerical)

        if (e0_over_eta == 0d0) then     !External file equilibrium

          if (my_rank == 0)
     .      write (*,*) 'Reading numerical equilibrium ',equil,'...'

          !Read file
          open (unit=u_equf,file=trim(equ_file),status='old')

          read(1000,'(a)')
          do i=0,nxdp
            read(1000,'(9f12.8)') dummy1,bth0(i),bzz0(i),jth0(i),jz0(i)
     .                  ,vr0(i),ez0(i),eth0(i),rho0(i)
cc            write (*,*)
cc            write (*,*) dummy1,bth0(i),bzz0(i),jth0(i),jz0(i)
cc     .                  ,vr0(i),ez0(i),eth0(i),rho0(i)
          enddo

          read(1000,'(a)')
          read(1000,'(a)')
          read(1000,'(f7.1)') eta_in

          close(1000)

          if (my_rank ==0) write (*,*) 'Done!'

          !Find electric field at face r=1d0
          E0(2) = 0.5*(eth0(nxd)+eth0(nxdp))
          E0(3) = 0.5*(ez0 (nxd)+ez0 (nxdp))
          E0 = E0*eta/eta_in

        else   !Calculate equilibrium

          if (my_rank == 0)
     .       write (*,*) 'Calculating equilibrium ',equil,'...'

          bz0 = 1d0

          !Find pressure and q-profiles
          do i=0,nxdp
            x1  = grid_params%xg(i)
            p00(i) = 0.5*beta*bz0**2*pprof(x1,0)
          enddo

	  !Initialize toroidal field
          bzz0 = bz0

! Daniele, 11/10/2013
          if (gamma > 1d0) then
             call Ohm_equil_1d_xport(-kk,e0_over_eta,theta_ppnch,pw
     .                              ,p00,bth0,bzz0,vr0,psi0)
          else
             call Ohm_equil_1d      (-kk,e0_over_eta,theta_ppnch
     .                              ,p00,bth0,bzz0,vr0,psi0)
          endif

	  !Find electric field
          E0(1:2) = 0d0
          E0(3)   = e0_over_eta*eta

!Daniele, 11/10/2013
c$$$          if (spitzer) then
c$$$            E0(3) = e0_over_eta*eta*res_spitzer(p00(1)/a_p)
c$$$          else
c$$$            E0(3) = e0_over_eta*eta
c$$$          endif

          if (spitzer) te0_ref = p00(1)/a_p !Reference Te for eta computation

          eta_in = eta

          !Adjust radial flow to actual resistivity
          vr0 = eta*vr0

          !Find axial flow from constant momentum source (Dirichlet at wall)
          vz0 = m0_over_nu*0.25*(1-grid_params%xg**2)

          !Find momentum source (cylindrical)
          M0(1:2) = 0d0
          M0(3) = m0_over_nu*nu

	  !Find density
          rho0 = 1d0

        endif

#if defined(RFX)
        if (my_rank == 0.and.(.not.test)) then
          write (*,*) 'r,bth0,bzz0,q0,vr0,p0'
          do i=0,nxdp
            x1 = grid_params%xg(i) !Global coordinate
            write(*,10) x1,bth0(i),bzz0(i)
     .             ,(x1*bzz0(i))/(RR*bth0(i)),vr0(i),p00(i)
          enddo
        endif
#endif

        !Convert current to covariant at r=1d0
        E0(3) = -kk/mm*E0(2) + E0(3)

        !Build equilibrium
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

#if !defined(vec_pot)
              var(i,j,k,IBX) = 0d0
              var(i,j,k,IBY) = bth0(iglobal) + kk*x1/mm*bzz0(iglobal)
              var(i,j,k,IBZ) = x1*bzz0(iglobal)
#endif
              var(i,j,k,IVX) = x1*vr0(iglobal)*eta/eta_in
              var(i,j,k,IVY) = kk*x1/mm*vz0(iglobal)
              var(i,j,k,IVZ) = x1*vz0(iglobal)

              if (solve_rho) then
                var(i,j,k,IRHO) = rho0(iglobal)
              else
                var(i,j,k,IRHO) = 1d0
              endif

              var(i,j,k,ITMP) = p00(iglobal)/a_p/var(i,j,k,IRHO)

cc              psi(i,j,k) = psi0(iglobal)  !Store flux function

            enddo
          enddo
        enddo

#if defined(vec_pot)
        allocate(bg(0:nxdp,0:ny+1,0:nz+1,3))

        do k = 0,nz+1
          do j = 0,ny+1
            do i = 0,nxdp

              x1 = grid_params%xg(i)  !Global coordinate

              bg(i,j,k,1) = 0d0
              bg(i,j,k,2) = bth0(i) + kk*x1/mm*bzz0(i)
              bg(i,j,k,3) = x1*bzz0(i)
            enddo
          enddo
        enddo

        call findVectorPotential(nxd,ny,nz,var(:,:,:,IVX:IVZ)
     .                          ,bg,var(:,:,:,IAX:IAZ))

        deallocate(bg)
#endif

      case ('2dhel')

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        theta_ppnch = grid_params%params(4)
        e0_over_eta = grid_params%params(5)
        beta = grid_params%params(6)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Initial guess: 1D equilibrium

        var(:,:,:,IBX) = 0d0

        bzz0 = 1d0
        p00  = 0d0

        if (gamma > 1d0) then
           call Ohm_equil_1d_xport(-kk,e0_over_eta,theta_ppnch,0d0
     .                            ,p00,bth0,bzz0,vr0,psi0)
        else
           call Ohm_equil_1d      (-kk,e0_over_eta,theta_ppnch
     .                            ,p00,bth0,bzz0,vr0,psi0)
        endif

        do i=0,nx+1
          iglobal = i + grid_params%ilo(igx) - 1
          var(i,:,:,IBY) = bth0(iglobal)
          var(i,:,:,IBZ) = bzz0(iglobal)
        enddo

cc        !Find poloidal flow
cc        var(:,:,:,IVX) = -(E0(3)*var(:,:,:,IBY)
cc     .                    -E0(2)*var(:,:,:,IBZ))
cc     .                    /(var(:,:,:,IBX)**2
cc     .                     +var(:,:,:,IBY)**2
cc     .                     +var(:,:,:,IBZ)**2)*eta/eta_in
cc        var(:,:,:,IVY) =  (E0(3)*var(:,:,:,IBX)
cc     .                    -E0(1)*var(:,:,:,IBZ))
cc     .                    /(var(:,:,:,IBX)**2
cc     .                     +var(:,:,:,IBY)**2
cc     .                     +var(:,:,:,IBZ)**2)*eta/eta_in
cc
cc        !Build equilibrium
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)
cc
cc              !X-Y equilibrium
cc#if !defined(vec_pot)
cc              var(i,j,k,IBX) = x1*var(i,j,k,IBX)
cc              var(i,j,k,IBY) = var(i,j,k,IBY) + kk*x1/mm*var(i,j,k,IBZ)
cc              var(i,j,k,IBZ) = x1*var(i,j,k,IBZ)
cc#endif
cc              var(i,j,k,IVX) = x1*var(i,j,k,IVX)*eta/eta_in
cc              var(i,j,k,IVZ) = 0d0
cc              var(i,j,k,IVY) = var(i,j,k,IVY) + kk*x1/mm*var(i,j,k,IVZ)
cc
cc            enddo
cc          enddo
cc        enddo

c     Helical equilibria

	!Solve for helical eq (returns B and V in cnv representation)
        call Ohm_hel_equil(nx,ny,nz,-kk/mm,nh2,e0_over_eta,theta_ppnch
     .                    ,psi0,var(:,:,:,IBX:IBZ),var(:,:,:,IVX:IVY))

        var(:,:,:,IVX:IVY) = var(:,:,:,IVX:IVY)*eta
        var(:,:,:,IVZ) = 0d0

	!Find electric field
        E0(1:2) = 0d0
        E0(3) = e0_over_eta*eta

        !Convert current to covariant at r=1d0
        E0(3) = -kk/mm*E0(2) + E0(3)

c     Density and temperature

        var(:,:,:,IRHO) = 1d0
        var(:,:,:,ITMP) = 1d-5

#if defined(vec_pot)
        allocate(bg(0:nxdp,0:ny+1,0:nz+1,3))

        do k = 0,nz+1
          do j = 0,ny+1
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              bg(i,j,k,1) = 0d0
              bg(i,j,k,2) = bth0(i) + kk*x1/mm*bzz0(i)
              bg(i,j,k,3) = x1*bzz0(i)
            enddo
          enddo
        enddo

        call findVectorPotential(nxd,ny,nz,var(:,:,:,IVX:IVZ)
     .                          ,bg,var(:,:,:,IAX:IAZ))

        deallocate(bg)
#endif

      case ('ppn3d','p3nsl')

        mm = grid_params%params(1)
        kk = grid_params%params(2)
        RR = grid_params%params(3)
        theta_ppnch = grid_params%params(4)
        e0_over_eta = grid_params%params(5)
        beta        = grid_params%params(6)
        m0_over_nu  = grid_params%params(7)
        pw = grid_params%params(8)

        nh2 = mm       !To set the right perturbation wavelength
        nh3 = kk*RR    !To set the right perturbation wavelength

        nh1 = nh2      !For perturbations to be of the right order

        if (nh2 == 0 .and. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

c     Check coordinates

        if (coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (abs(RR-zmax/2./pi) > 1d-5) then
          messg = 'Major radius and zmax do not agree'
          call pstop('setEquilibrium',messg)
        endif

c     Equilibria (numerical)

        if (e0_over_eta == 0d0) then     !External file equilibrium

          if (my_rank == 0)
     .      write (*,*) 'Reading numerical equilibrium ',equil,'...'

          !Read file
          open (unit=1000,file=trim(equ_file),status='old')

          read(1000,'(a)')
          do i=0,nxdp
            read(1000,'(9f12.8)') dummy1,bth0(i),bzz0(i),jth0(i),jz0(i)
     .                  ,vr0(i),ez0(i),eth0(i),rho0(i)
          enddo

          read(1000,'(a)')
          read(1000,'(a)')
          read(1000,'(f7.1)') eta_in

          close(1000)

          if (my_rank ==0) write (*,*) 'Done!'

          !Find electric field at face r=1d0
          E0(2) = 0.5*(eth0(nxd)+eth0(nxdp))
          E0(3) = 0.5*(ez0 (nxd)+ez0 (nxdp))
          E0 = E0*eta/eta_in

        else   !Calculate equilibrium

          if (my_rank == 0)
     .       write (*,*) 'Calculating equilibrium ',equil,'...'

          bz0 = 1d0

          !Find pressure and q-profiles
          do i=0,nxdp
            x1  = grid_params%xg(i)
            p00(i) = 0.5*beta*bz0**2*pprof(x1,0)
          enddo

	  !Initialize toroidal field
          bzz0 = bz0

          if (gamma > 1d0) then
             call Ohm_equil_1d_xport(-kk,e0_over_eta,theta_ppnch,pw
     .                              ,p00,bth0,bzz0,vr0,psi0)
          else
             call Ohm_equil_1d      (-kk,e0_over_eta,theta_ppnch
     .                              ,p00,bth0,bzz0,vr0,psi0)
          endif

          !Find electric field
          E0(1:2) = 0d0
          E0(3)   = e0_over_eta*eta
          eta_in  = eta

          !Adjust radial flow to actual resistivity
          vr0 = eta*vr0

          !Find axial flow from constant momentum source (Dirichlet at wall)
          vz0 = m0_over_nu*0.25*(1-grid_params%xg**2)

          !Find momentum source (cylindrical)
          M0(1:2) = 0d0
          M0(3) = m0_over_nu*nu

          !Find density
          rho0 = 1d0

        endif

#if defined(RFX)
        if (my_rank == 0) then
          write (*,*) 'r,bth0,bzz0,q0,vr0,vz0'
          do i=0,nxdp
            x1 = grid_params%xg(i) !Global coordinate
            write(*,10) x1,bth0(i),bzz0(i)
     .             ,(x1*bzz0(i))/(RR*bth0(i)),vr0(i),vz0(i)
          enddo
        endif
#else
cc        write (*,*) 'r,bth0,bzz0,jth0,jz0,vr0'
cc        do i=0,nxdp
cc          call getCurvilinearCoordinates(i,1,1,igx,igy,igz,ig,jg,kg
cc     .                                  ,x1,y1,z1)
cc          write(*,*) x1,bth0(i),bzz0(i),vr0(i)
cc        enddo
cc        stop
#endif

        !Build equilibrium
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
     .                                      ,x1,y1,z1)

              iglobal = i + grid_params%ilo(igx) - 1

#if !defined(vec_pot)
              var(i,j,k,IBX) = 0d0
              var(i,j,k,IBY) = bth0(iglobal)
              var(i,j,k,IBZ) = x1*bzz0(iglobal)
#endif
              var(i,j,k,IVX) = x1*vr0(iglobal)*eta/eta_in
              var(i,j,k,IVY) = 0d0
              var(i,j,k,IVZ) = x1*vz0(iglobal)

              if (solve_rho) then
                var(i,j,k,IRHO) = rho0(iglobal)
              else
                var(i,j,k,IRHO) = 1d0
              endif

              if (equil == 'ppn3d') then
                var(i,j,k,ITMP) = 1d-5
              else
                var(i,j,k,ITMP) = p00(iglobal)/a_p/var(i,j,k,IRHO)
              endif

cc              psi(i,j,k) = psi0(iglobal)  !Store flux function

            enddo
          enddo
        enddo

#if defined(vec_pot)
        allocate(bg(0:nxdp,0:ny+1,0:nz+1,3))

        do k = 0,nz+1
          do j = 0,ny+1
            do i = 0,nxdp
              x1 = grid_params%xg(i)  !Global coordinate

              bg(i,j,k,1) = 0d0
              bg(i,j,k,2) = bth0(i)
              bg(i,j,k,3) = x1*bzz0(i)
            enddo
          enddo
        enddo

        call findVectorPotential(nxd,ny,nz,var(:,:,:,IVX:IVZ)
     .                          ,bg,var(:,:,:,IAX:IAZ))

        deallocate(bg)
#endif

      case ('ohtor')

        !Parameters 2 and 3 define elliptical shape; used in grid_anal_map
        RR = grid_params%params(1)              !Major radius
cc        mm = grid_params%params(4)              !m mode
cc        nn = grid_params%params(5)              !n mode
        e0_over_eta = grid_params%params(6)*RR  !Covariant component is constant, needs RR factor
        beta = grid_params%params(7)     !Toroidal beta

        kk = 0d0   !To avoid affecting q-profile

        !Now these must be initalized directly in input file.
cc        nh2 =  mm       !To set the right perturbation wavelength
cc        nh3 = -nn       !To set the right perturbation wavelength

        nh1 = nh2       !For perturbations to be of the right order
                        !(nh2, nh3 set in input file)

        if (nh2 == 0 .and. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

c     Check coordinates

        if (coords /= 'tor') then
          messg = 'Need "tor" coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Ohmic GS equilibrium

        !Build equilibrium
        call Ohm_GS_equil_2d(RR,e0_over_eta,beta,nx,ny,nz
     .                   ,var(:,:,:,IBX),var(:,:,:,IBY),var(:,:,:,IBZ)
     .                   ,var(:,:,:,IVX),var(:,:,:,IVY),var(:,:,:,ITMP))

	!Find electric field
        E0(1:2) = 0d0
        E0(3)   = e0_over_eta*eta
        eta_in  = eta

        !Adjust pinch flow to current resistivity
        var(:,:,:,IVX) = var(:,:,:,IVX)*eta
        var(:,:,:,IVY) = var(:,:,:,IVY)*eta
        var(:,:,:,IVZ) = 0d0

        !Fill other arrays
        var(:,:,:,IRHO) = 1d0
        var(:,:,:,ITMP) = var(:,:,:,ITMP)/a_p/var(:,:,:,IRHO)

#if defined(vec_pot)
        call pstop('setEquilibrium','Equilibrium not working')

cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)
cc
cc              btheta = bb*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-bb**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              call getCartesianCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,x1,y1,z1)
cc
cc              !X-Y equilibrium
cc              b(i,j,k,1) = 0d0
cc              b(i,j,k,2) = btheta
cc              b(i,j,k,3) = gmetric%grid(igx)%gsup(i,j,k,3,3)*bz0
cc            enddo
cc          enddo
cc        enddo
cc
cc        call findVectorPotential(var(:,:,:,IVX:IVZ)
cc     .                          ,b
cc     .                          ,var(:,:,:,IAX:IAZ))
#endif

      case default

        messg = 'Equilibrium '//trim(equil)//' undefined'
        call pstop('setEquilibrium',messg)

      end select

c Find total pressure

      var(:,:,:,ITMP) = a_p*var(:,:,:,ITMP)*var(:,:,:,IRHO)

c Find momentum components

      if (.not.nc_eom_v) then
        var(:,:,:,IVX) = var(:,:,:,IRHO)*var(:,:,:,IVX)
        var(:,:,:,IVY) = var(:,:,:,IRHO)*var(:,:,:,IVY)
        var(:,:,:,IVZ) = var(:,:,:,IRHO)*var(:,:,:,IVZ)
      endif

c Apply Jacobian factors

      var(:,:,:,IRHO) = var(:,:,:,IRHO)*gmetric%grid(iigx)%jac
      var(:,:,:,ITMP) = var(:,:,:,ITMP)*gmetric%grid(iigx)%jac

      if (alt__eom()) then
        var(:,:,:,IVY) = var(:,:,:,IVY)*gmetric%grid(iigx)%jac
      endif

c Transfer equilibrium to varray

      do ieq=1,neqd
        varray%array_var(ieq)%array = var(:,:,:,ieq) 
      enddo

c Transfer BC fields info to postprocessor via input file

      call write_bc_field(E0,'E0')
      call write_bc_field(B0,'B0')

c End program

10    format (4f7.3,2e12.3)

      contains

c     write_bc_field
c     ###############################################################
      subroutine write_bc_field(var,desc)

c     ---------------------------------------------------------------
c     Writes boundary electric and magnetic fields to pixie3d.in
c     for later reading.
c     ---------------------------------------------------------------

c     Call variables

      character(*) :: desc
      real(8) :: var(3)

c     Local variables

      character(200):: command,line

c     Begin program

c     Transfer electric field info to postprocessor via input file

      if (sum(var) /= 0d0 .and. my_rank == 0) then
        command = 'grep '//trim(desc)//' '//trim(inputfile)
     .           // ' > /dev/null'
        ierr = shellcmd(trim(command))

        if (ierr == 0) ierr = delete_line(inputfile,desc)

        write (line,*) '   '//trim(desc)//'=',real(var(1),4)
     .                                   ,',',real(var(2),4)
     .                                   ,',',real(var(3),4)

        command = 'grep datin ' // trim(inputfile) // ' > /dev/null'
        ierr = shellcmd(trim(command))

        if (ierr == 0) then
          ierr = add_line_after_txt(inputfile,'datin',line)
        else
          ierr = add_line_after_txt(inputfile,'phys',line)
        endif

        if (ierr /= 0) then
          command = 'Could not write '//trim(desc)//' to input file'
          call pstop('setEq',command)
        endif
      endif

      end subroutine write_bc_field

      end subroutine setEquilibrium

c     GS_equil_1d
c     ###############################################################
      subroutine GS_equil_1d(Rxq0,p0,bth0,bzz0)

c     ---------------------------------------------------------------
c     Solves 1D Grad-Shafranov equation to find cylindrical components
c     of magnetic field B_theta, B_z from specified q-profile and
c     pressure profile. On input:
c       * RR: major radius (aspect ratio -> 1/RR)
c       * Rxq0,p0: normalized q-profile (RR x q) and pressure profile
c       * bth0,bzz0: output magnetic field components (bzz0 contains
c                    initial guess, and bzz0(0) contains value at axis)
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use grid

        use app_iosetup

        implicit none

c     Call variables

        real(8)    :: bth0(0:nxdp),bzz0(0:nxdp)
     .               ,p0(0:nxdp),Rxq0(0:nxdp)

c     Local variables

        integer    :: i,j,k,ig,jg,kg,it
        real(8)    :: bz_avg,x1,bzh,bth,gs_err,rh,dh

        real(8)    :: bz0,rint(0:nxdp),dummy(0:nxdp)
     .               ,p0h(0:nxdp),q0h(0:nxdp)
     .               ,xh(0:nxd),dr(0:nxd)

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '************************************'
          write (*,*) '1D Grad-Shafranov equilibrium solver'
          write (*,*) '************************************'
          write (*,*)
        endif

c       Initialize iteration

        bz0 = bzz0(0)

        !Find half-mesh positions
        do i=0,nxd
          xh(i) = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
          dr(i) =      grid_params%xg(i+1)-grid_params%xg(i)
        enddo

        !Average pressure and q-profiles to half mesh
        do i=0,nxd
          p0h(i) = 0.5*(p0(i)  +p0  (i+1)) !Average pressure  to half mesh
          q0h(i) = 0.5*(Rxq0(i)+Rxq0(i+1)) !Average q-profile to half mesh
        enddo

c       Iteration

        do it=1,100

cc          !Perform integral of B_z/qprof on half-mesh
cc          dummy(0) = 0d0
cc          bzz0 (0) = bz0
cc          do i=1,nxd
cc            rh  = grid_params%xg(i)
cc
cc            bzh = 0.5*(bzz0(i)+bzz0(i-1))
cc
cc            dummy(i) = dummy(i-1) + 2*(bzh/Rxq0(i))**2*rh*dr(i)
cc          enddo

          !Perform integral of B_th^2/r on half-mesh
          dummy(0) = 0d0
          do i=1,nxd
            rh  = grid_params%xg(i)
            bth = 0.5*(bth0(i)+bth0(i-1))
            dummy(i) = dummy(i-1) + 2*bth**2/rh*dr(i)
          enddo

          !Safeguard for solvability
          where (dummy>bz0**2+2*(p0h(0)-p0h))
     .           dummy=bz0**2+2*(p0h(0)-p0h)

          !Calculate Bz, Btheta on half-mesh
          do i=0,nxd
            x1  = xh(i)
            bth0(i) = x1*bzz0(i)/q0h(i)
            bzz0(i) = sign(1d0,q0h(i))
     .               *sqrt((bz0**2+2*(p0h(0)-p0h(i))-dummy(i))
     .                    /((x1/q0h(i))**2+1))
          enddo

          !Check convergence (force balance)
          gs_err = 0d0
          do i=1,nxd
            rh  = grid_params%xg(i)

            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))
            gs_err = gs_err
     .                +(bzh*(bzz0(i)-bzz0(i-1))/dr(i)
     .                     +(p0h (i)-p0h (i-1))/dr(i)
cc     $                 +bth/rh*(bth0(i)*i-bth0(i-1)*(i-1)))**2
     $                 +0.5*(bth0(i)**2-bth0(i-1)**2)/dr(i)
     .                 +bth**2/rh)**2
          enddo
          gs_err = sqrt(gs_err/nxd)

          if (ilevel > 0 .and. my_rank == 0) then
            write (*,*) '1D Grad-Shafranov equilibrium iter =',it
     .                 ,' ; Force balance rror=',gs_err
          endif

          dh = xmax/nxd
          if (gs_err < 0.01*dh**2) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
           write (*,'(a,i3,a,1pe10.2)')
     .          ' 1D equilibrium converged in '
     .          ,it,' iterations with force-balance error =',gs_err
        endif

cc        !Check force balance
cc        do i=1,nxd
cc          rh  = grid_params%xg(i)
cc
cc          bzh = 0.5*(bzz0(i)+bzz0(i-1))
cc          bth = 0.5*(bth0(i)+bth0(i-1))
cc          dummy(i) =  bzh*(bzz0(i)-bzz0(i-1))/dr(i)
cc     .                   +(p0h (i)-p0h (i-1))/dr(i)
cc     $              +bth/rh*(bth0(i)*i-bth0(i-1)*(i-1))
cc        enddo
cc
cc        if (my_rank == 0) then
cc          write(*,*)
cc          write(*,'(a,1pe14.7)') ' Force-balance error ='
cc     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
cccc          write (*,*) dummy(1:nxd)
cc        endif

        !Calculate equilibrium properties
        bz_avg = 0d0
        do i=1,nxd
          x1  = xh(i)
          bz_avg = bz_avg + 2*bzz0(i)*x1*dr(i)
        enddo

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' 1D Grad-Shafranov equilibrium properties:'
          write (*,999) '   Theta  =',bth0(nxd)/bz_avg
          write (*,999) '   Bz(r=0)=',bzz0(0)
          write (*,999) '   I0     =',2*pi*xmax*bth0(nxd)
          write (*,999) '   beta_t =',2*p0h(0)/(bzz0(0)**2)
          write (*,999) '   beta_p =',2*p0h(0)/(bth0(nxd)**2)
 999      format (a,f10.5)
        endif

        !Save some equilibrium properties (to calculate bootstrap current)
        eq_params(1) = 1d0/grid_params%params(3)       !Aspect ratio
        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
        eq_params(3) = q0h(nxd)/grid_params%params(3)  !q(r=a)

        !Transfer to integer mesh
        call IntDriver1d(nxdp,xh(0:nxd),bth0(0:nxd)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bth0 = dummy

        call IntDriver1d(nxdp,xh(0:nxd),bzz0(0:nxd)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bzz0 = dummy

        !Extrapolate linearly in radius @ r=0
        bth0(0) = -bth0(1)
        bzz0(0) = bzz0(1)

        !Dump equilibrium for eigen code

        if (my_rank == 0 .and. equil == 'ntm2d') then
          open(unit=123,file='ntm2d.txt',status='unknown')
          write (123,*) nxdp+1
          write (123,*) grid_params%xg(0:nxdp)
          write (123,*) bth0
          write (123,*) bzz0
          write (123,*) 0.5*p0
          write (123,*) grid_params%params(3),eq_params(2),eq_params(3)
          close(123)
cc          stop
        endif

      end subroutine GS_equil_1d

c     Ohm_equil_1d
c     ###############################################################
      subroutine Ohm_equil_1d(kk,Eoeta,theta,p0,bth0,bzz0,vr0,psi0)

c     ---------------------------------------------------------------
c     Solves 1D Ohmic eq. equations to find cylindrical components
c     of magnetic field B_theta, B_z and radial velocity vr from 
c     imposed electric field, theta parameter, and pressure profile.
c     On input:
c       * kk: helical pitch (>0 by convention)
c       * Eoeta: electric field over eta
c       * theta: RFP theta parameter
c       * p0: input presure profile
c     On output:
c       * bth0,bzz0 (output): magnetic field components (bzz0 contains
c                    initial guess, and bzz0(0) contains value at axis)
c       * vr0  (output): radial pinch flow
c       * psi0 (output): helical poloidal flux
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use nlfunction_setup, ONLY: res_rfx,res_spitzer,eta,spitzer,a_p

        use grid

        use app_iosetup

        implicit none

c     Call variables

        real(8),intent(IN) :: kk,Eoeta,theta,p0(0:nxdp)

        real(8),intent(OUT):: bth0(0:nxdp),bzz0(0:nxdp)
     .                       ,vr0 (0:nxdp),psi0(0:nxdp)

c     Local variables

        integer :: i,j,k,ig,jg,kg,it
        real(8) :: bz0,bz_avg,x1,bzh,bth,err,rh,dh,p00,etah

        real(8) :: rint(0:nxdp),dummy(0:nxdp),bth0_old(0:nxdp)
     .            ,bzz0_old(0:nxdp),xh(0:nxd),dr(0:nxd)

        real(8) :: pprof
        external   pprof

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '************************************'
          write (*,*) '1D Ohmic equilibrium solver'
          write (*,*) '************************************'
          write (*,*)
        endif

c       Initialize iteration

        bz0 = bzz0(0)
        p00 = maxval(p0)

        !Find half-mesh positions
        do i=0,nxd
          xh(i) = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
          dr(i) =      grid_params%xg(i+1)-grid_params%xg(i)
        enddo

        !Initial guess for B_theta (defined in half mesh)
        bth0 = 0d0
        do i=0,nxd
          bth0(i) = 0.5*xh(i)*Eoeta
        enddo

c       Iteration

        do it=1,100

          bth0_old = bth0
          bzz0_old = bzz0

          !Perform integral of B_theta
          dummy(0) = 0d0
          bth0 (0) = 0d0
          do i=1,nxd
            x1 = xh(i)
            rh = grid_params%xg(i)

            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))

            if (spitzer) then
cc              etah = res_spitzer(p0(i)/a_p)
              etah = res_spitzer(p0(i)/p0(1))
            else
              etah = res_rfx(rh)
            endif

            dummy(i)=dummy(i-1)
     .              +(Eoeta/etah*bzh**2
     .               -bth*(p0(i+1)-p0(i-1))/2./dr(i))
     .                /(bzh**2+bth**2)*rh*dr(i)

            bth0(i) = dummy(i)/x1
          enddo

          !Perform integral of Bz
          dummy = 0d0
          do i=1,nxd
            x1 = xh(i)
            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))
            rh  = grid_params%xg(i)

            if (spitzer) then
cc              etah = res_spitzer(p0(i)/a_p)
              etah = res_spitzer(p0(i)/p0(1))
            else
              etah = res_rfx(rh)
            endif

            dummy(i) = dummy(i-1)
     .                +(Eoeta/etah*bth+(p0(i+1)-p0(i-1))/2./dr(i))
     .                /(bzh**2+bth**2)*dr(i)
          enddo

          bzz0 = exp(-dummy)

          if (theta /= 0d0) then
            bz_avg = 0d0
            do i=1,nxd
               bz_avg = bz_avg + (bzz0(i)*i + bzz0(i-1)*(i-1))*dr(i)**2
            enddo
            bz_avg = bz_avg/xmax**2
            bzz0 = bth0(nxd)/theta/bz_avg*bzz0
          elseif (bz0 /= 0d0) then
            bzz0 = bz0*bzz0
          else
            call pstop('setEquilibrium'
     $                ,'Cannot find PPNCH equilibrium:'
     $                 //' Both Bz0 and theta are zero')
          endif

          !Pressure solve

          !Check convergence
          err = sqrt(0.5*(sum((bth0-bth0_old)**2)
     .                   +sum((bzz0-bzz0_old)**2))/nxdp)

          if (ilevel > 0 .and. my_rank == 0) then
            write (*,*) '1D Ohmic equilibrium iter =',it
     .                 ,' ; Error=',err
          endif

          dh = xmax/nxd
          if (err < 0.1*dh**2) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
           write (*,'(a,i3,a,1pe10.2)')
     .          ' 1D equilibrium converged in '
     .          ,it,' iterations with error =',err
        endif

        !Check force balance
        do i=1,nxd
          rh  = grid_params%xg(i)

          bzh = 0.5*(bzz0(i)+bzz0(i-1))
          bth = 0.5*(bth0(i)+bth0(i-1))
          dummy(i) = bzh*(bzz0(i)-bzz0(i-1))/dr(i)
     $              +bth/rh*(bth0(i)*i-bth0(i-1)*(i-1))
     .              +(p0(i+1)-p0(i-1))/2./dr(i)
        enddo

        if (my_rank == 0) then
          write(*,*)
          write(*,'(a,1pe14.7)') ' Force-balance error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
        endif

        !Check Ohm`s law
        do i=1,nxd
          bzh = 0.5*(bzz0(i)+bzz0(i-1))
          bth = 0.5*(bth0(i)+bth0(i-1))
          rh  = grid_params%xg(i)

          if (spitzer) then
cc            etah = res_spitzer(p0(i)/a_p)
            etah = res_spitzer(p0(i)/p0(1))
          else
            etah = res_rfx(rh)
          endif

          dummy(i) = Eoeta*bzh/etah
     $              -bzh/rh*(bth0(i)*i-bth0(i-1)*(i-1))
     $              +bth*(bzz0(i)-bzz0(i-1))/dr(i)
        enddo

        if (my_rank == 0) then
          write(*,'(a,1pe14.7)') ' Ohms law      error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
        endif

        !Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)
        bz_avg = 0d0
        do i=1,nxd
          bz_avg = bz_avg + (bzz0(i)*i + bzz0(i-1)*(i-1))*dr(i)**2
        enddo

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' Pinch equilibrium features:'
          write (*,999) '   Theta  =',bth0(nxd)/bz_avg
          write (*,999) '   Bz(r=0)=',bzz0(0)
          write (*,999) '   I0     =',2*pi*xmax*bth0(nxd)
          write (*,999) '   E0/eta =',Eoeta
          write (*,999) '   T_flux =',pi*xmax**2*bz_avg
 999      format (a,f10.7)
        endif

        !Save some equilibrium properties (to calculate bootstrap current)
cc        eq_params(1) = 1d0/grid_params%params(3)       !Aspect ratio
cc        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
cc        eq_params(3) = q0h(nxd)/grid_params%params(3)  !q(r=a)

        !Find helical flux (integer mesh)
        psi0(0) = 0d0
        do i=1,nxd+1
          psi0(i) = psi0(i-1)+dr(i-1)*(bth0(i-1)-kk*xh(i-1)*bzz0(i-1))
        enddo
        psi0 = psi0 - 0.5*(psi0(nxd)+psi0(nxdp)) !Homogeneous Dirichlet BC

        !Transfer B to integer mesh
        call IntDriver1d(nxdp,xh(0:nxd),bth0(0:nxd)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bth0 = dummy

        call IntDriver1d(nxdp,xh(0:nxd),bzz0(0:nxd)
     .                  ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bzz0 = dummy

        !Extrapolate linearly in radius @ r=0
        bth0(0) = -bth0(1)
        bzz0(0) = bzz0(1)

        !Find pinch flow (normalized to nominal resistivity)
        vr0 = 0d0
        do i=1,nxd
          rh  = grid_params%xg(i)
          if (spitzer) then
cc            etah = res_spitzer(p0(i)/a_p)
            etah = res_spitzer(p0(i)/p0(1))
          else
            etah = res_rfx(rh)
          endif

          vr0(i) = -(Eoeta*bth0(i)+etah*(p0(i+1)-p0(i-1))/2./dr(i))
     .             /(bzz0(i)**2+bth0(i)**2)
        enddo

        !Find pressure
cc        do i=1,nxd
cc          rh  = grid_params%xg(i)
cc          p0(i) = p00*pprof(rh,0)
cc        enddo

        !Dump equilibrium for eigen code
cc        if (my_rank == 0) then
cc          open(unit=123,file='ntm2d.txt',status='unknown')
cc          write (123,*) nxdp+1
cc          write (123,*) grid_params%xg(0:nxdp)
cc          write (123,*) bth0
cc          write (123,*) bzz0
cc          write (123,*) 0.5*p0
cc          write (123,*) grid_params%params(3),eq_params(2),eq_params(3)
cc          close(123)
cccc          stop
cc        endif

      end subroutine Ohm_equil_1d

c     Ohm_equil_1d_xport
c     ###############################################################
      subroutine Ohm_equil_1d_xport(kk,Eoeta,theta,pw,p0,bth0,bzz0,vr0
     .                             ,psi0)

c     ---------------------------------------------------------------
c     Solves 1D Ohmic eq. equations coupled with transport to find
c     cylindrical components  of magnetic field B_theta, B_z ,
c     and pressure profile, and radial velocity vr from imposed 
c     electric field, theta parameter, and beta parameter.
c     On input:
c       * kk (input): helical pitch (>0 by convention)
c       * Eoeta (input): electric field over eta
c       * theta (input): RFP theta parameter
c       * pw (input): pressure at wall boundary (r=1)
c     On output:
c       * p0 (output): presure profile (satisfies p(r=1)=pw)
c       * bth0,bzz0 (output): magnetic field components (bzz0 contains
c                    initial guess, and bzz0(0) contains value at axis)
c       * vr0 (output): radial pinch flow
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use nlfunction_setup, ONLY: res_rfx,eta,gamma,dfvty_rfx
     .                             ,res_spitzer,spitzer,a_p

        use grid

        use app_iosetup

        implicit none

c     Call variables

        real(8) :: kk,bth0(0:nxdp),bzz0(0:nxdp),vr0(0:nxdp)
     .            ,p0(0:nxdp),Eoeta,theta,b_tor,psi0(0:nxdp)
     .            ,pw

c     Local variables

        integer :: i,j,k,ig,jg,kg,it,info
        real(8) :: bz_avg,x1,bzh,bth,err,rh,dh,p00,rhm,rhp
! Daniele, 02/07/2014
        real(8) :: chirh,chirhm,chirhp

        real(8) :: bz0,rint(0:nxdp),dummy(0:nxdp),bth0_old(0:nxdp)
     .            ,bzz0_old(0:nxdp),xh(0:nxd),dr(0:nxd)
     .            ,bzf(0:nxdp),btf(0:nxdp),p0_old(0:nxdp)
     .            ,jzf(nxd),jtf(nxd),dl(nxd),d(nxd),du(nxd)
     .            ,etah(0:nxdp)

        real(8) :: pprof
        external   pprof

! Daniele, 11/10/2013
        real(8) :: alpha,etah0
        real(8) :: chih(0:nxd)

! Daniele, 1/8/2014
        real(8) :: pa

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '**************************************'
          write (*,*) '1D Ohmic equilibrium solver with Xport'
          write (*,*) '**************************************'
          write (*,*)
        endif

! Daniele, 11/10/2013
!        alpha = 0d0 !No underrelaxation
!        alpha = 0.5d0 !underrelaxation constant
! Daniele, 02/07/2014
        alpha = 0.8d0 !underrelaxation constant
        etah0 = 1d0 ! needed if not spitzer

! Daniele, 1/8/2014; Luis 1/23/2014
        pa = pw

c       Initialize iteration

        bz0 = bzz0(0)

        !Find half-mesh positions
        do i=0,nxd
          dr(i) =      grid_params%xg(i+1)-grid_params%xg(i)
          xh(i) = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
! Daniele, 02/07/2014
c$$$!          chih(i)=chi*(1.-(1.-1d-2)*sin(pi*xh(i)/1.75))
c$$$          chih(i)=dfvty_rfx(xh(i))
! Luis 04/15/2014, to mimic what is done in PIXIE3D
          chih(i)=0.5*(dfvty_rfx(grid_params%xg(i+1))
     $                +dfvty_rfx(grid_params%xg(i  )))
        enddo

        !Initial guess for B_theta (defined in half mesh)
        bth0 = 0d0
        do i=0,nxd
          bth0(i) = 0.5*xh(i)*Eoeta
        enddo

! Daniele, 11/10/2013
        p0_old = p0

c       Iteration

        do it=1,100

          !Daniele, 11/10/2013
          !underrelaxation of p
          p0 = (1d0-alpha)*p0 + alpha*p0_old

          !Luis, 11/15/2013, find resistivity profile
          if (spitzer) then
            do i=1,nxdp
c$$$              etah(i) = res_spitzer(p0(i)/a_p)
              etah(i) = res_spitzer(p0(i)/p0(1))
            enddo
c$$$            etah0 = etah(1)
c$$$            etah  = etah/etah0
          else
            do i=1,nxdp
              rh = grid_params%xg(i)
              etah(i) = res_rfx(rh)
            enddo
          endif

          bth0_old = bth0
          bzz0_old = bzz0
          p0_old   = p0

          !Perform integral of B_theta
          dummy(0) = 0d0
          bth0 (0) = 0d0
          do i=1,nxd
            x1 = xh(i)
            rh = grid_params%xg(i)

            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))

            dummy(i)=dummy(i-1)
     .              +(Eoeta/etah(i)*bzh**2
     .               -bth*(p0(i+1)-p0(i-1))/2./dr(i))
     .                /(bzh**2+bth**2)*rh*dr(i)

            bth0(i) = dummy(i)/x1
          enddo

          !Perform integral of Bz
          dummy = 0d0
          do i=1,nxd
            x1 = xh(i)
            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))
            rh  = grid_params%xg(i)

            dummy(i) = dummy(i-1)
     .                +(Eoeta/etah(i)*bth+(p0(i+1)-p0(i-1))/2./dr(i))
     .                /(bzh**2+bth**2)*dr(i)
          enddo

          bzz0 = exp(-dummy)

          if (theta /= 0d0) then
            bz_avg = 0d0
            do i=1,nxd
               bz_avg = bz_avg + (bzz0(i)*i + bzz0(i-1)*(i-1))*dr(i)**2
            enddo
            bz_avg = bz_avg/xmax**2
            bzz0 = bth0(nxd)/theta/bz_avg*bzz0
          elseif (bz0 /= 0d0) then
            bzz0 = bz0*bzz0
          else
            call pstop('setEquilibrium'
     $                ,'Cannot find PPNCH equilibrium:'
     $                 //' Both Bz0 and theta are zero')
          endif

! Daniele, 12/7/2011
! f stands for full mesh
          !Transfer to integer mesh
          call IntDriver1d(nxdp,xh(0:nxd),bth0(0:nxd)
     .                    ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                    ,2,0)
          btf = dummy

          call IntDriver1d(nxdp,xh(0:nxd),bzz0(0:nxd)
     .                    ,nxdp,grid_params%xg(1:nxdp),dummy(1:nxdp)
     .                    ,2,0)
          bzf = dummy

          !Extrapolate linearly in radius @ r=0
          btf(0) = -btf(1)
          bzf(0) =  bzf(1)

          !Find pinch flow (NOT normalized to nominal resistivity)
          vr0 = 0d0
          do i=1,nxd
            vr0(i) = -(Eoeta*btf(i)
     .             +etah(i)*(p0(i+1)-p0(i-1))/2./dr(i))
     .             /(bzf(i)**2+btf(i)**2)
          enddo
! Daniele, 12/7/2011
! dr(nxdp)=0. so I use dr(nxd)
          vr0(nxdp)= -(Eoeta*btf(nxdp)
c$$$     .         +etah(nxdp)
c$$$     $         *(1.5*p0(nxdp)-2.*p0(nxd)+0.5*p0(nxd-1))/dr(nxd))
! Daniele, 1/8/2014
     .         +etah(nxdp)*(p0(nxdp)-p0(nxd))/dr(nxd))
     .         /(bzf(nxdp)**2+btf(nxdp)**2)

          vr0 = etah0*vr0

          vr0(0) = -vr0(1)

          !Find current j on integer mesh
          do i=1,nxd
            rh = grid_params%xg(i)
            jzf(i)=(xh(i)*bth0(i)-xh(i-1)*bth0(i-1))/(rh*dr(i))
            jtf(i)=-(bzz0(i)-bzz0(i-1))/dr(i)
          enddo

          !Pressure solve

          !Find heat source on integer mesh
          do i=1,nxd
            rh = grid_params%xg(i)
            p0(i)=rh*etah(i)*(jtf(i)*jtf(i)+jzf(i)*jzf(i))
          enddo
          p0 = etah0*p0

! Daniele, 1/8/2014. BCs with finite p(a)
          rh = grid_params%xg(nxd)
          p0(nxd) = p0(nxd) + pa/dr(nxd)*
     .         (-gamma/(gamma-1.)*xh(nxd)*vr0(nxd)
     .          +vr0(nxd)*rh
     .          +2.*chih(nxd)/eta*xh(nxd)/dr(nxd))

          !Build the tridiagonal matrix to be inverted
          do i=1,nxd
            rh  = grid_params%xg(i)

            dl(i)= -gamma/(gamma-1.)*xh(i-1)*vr0(i)/(2.*dr(i))
     .            +vr0(i)*rh/(2.*dr(i))
!     .            -chi/eta*xh(i-1)/dr(i)**2.
! Daniele, 11/10/2013
     .            -chih(i-1)/eta*xh(i-1)/dr(i)**2.

            d(i)=  gamma/(gamma-1.)
     .            *(xh(i)*vr0(i+1)-xh(i-1)*vr0(i-1))/(2.*dr(i))
     .            +0.
!     .            +chi/eta*(xh(i)+xh(i-1))/dr(i)**2.
! Daniele, 11/10/2013
     .            +1./eta*(chih(i)*xh(i)+chih(i-1)*xh(i-1))/dr(i)**2.
           
            du(i)= gamma/(gamma-1.)
     .            *(xh(i)*vr0(i))/(2.*dr(i))
     .            -vr0(i)*rh/(2.*dr(i))
!     .            -chi/eta*xh(i)/dr(i)**2.
! Daniele, 11/10/2013
     .            -chih(i)/eta*xh(i)/dr(i)**2.
          enddo

          !BCs at i=1
          rh  = grid_params%xg(1)
          d(1)=d(1)+vr0(1)*rh/(2.*dr(1))

          !BSc at i=nxd
          rh  = grid_params%xg(nxd)
          d(nxd)=d(nxd)
     .          -gamma/(gamma-1.)*xh(nxd)*vr0(nxd)/(2.*dr(nxd))
     .          +vr0(nxd)*rh/(2.*dr(nxd))
!     .         +chi/eta*xh(nxd)/dr(nxd)**2.
! Daniele, 11/10/2013
     .         +chih(nxd)/eta*xh(nxd)/dr(nxd)**2.

          !Invert tridiagonal matrix
          call DGTSV(nxd,1,dl(2:nxd),d,du(1:nxd-1),p0(1:nxd),nxd,info)

          !BCs
          p0(0   ) =  p0(1)
!          p0(nxdp) = -p0(nxd)
! Daniele, 1/8/2014. BCs with finite p(a)
          p0(nxdp) = -p0(nxd) + 2.*pa

          !Check convergence
          err = sqrt(0.3*(sum((bth0-bth0_old)**2)
     .                   +sum((bzz0-bzz0_old)**2)
     .                   +sum((p0  -p0_old  )**2))/nxdp)

          if (ilevel > 0 .and. my_rank == 0) then
            write (*,*) '1D Ohmic equilibrium iter =',it
     .                  ,' ; Error=',err
          endif

          dh = xmax/nxd
          if (err < 0.1*dh**2) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
          write (*,'(a,i3,a,1pe10.2)')
     .          ' 1D equilibrium converged in '
     .          ,it,' iterations with error =',err
        endif

        !Recompute resistivity (Luis, 11/15/2013)
        if (spitzer) then
          do i=1,nxd
c$$$            etah(i) = res_spitzer(p0(i)/a_p)
            etah(i) = res_spitzer(p0(i)/p0(1))
          enddo
c$$$          etah0 = etah(1)
c$$$          etah  = etah/etah0
        else
          do i=1,nxd
            rh = grid_params%xg(i)
            etah(i) = res_rfx(rh)
          enddo
        endif

        !Check force balance
        do i=1,nxd
          rh  = grid_params%xg(i)

          bzh = 0.5*(bzz0(i)+bzz0(i-1))
          bth = 0.5*(bth0(i)+bth0(i-1))
          dummy(i) = bzh*(bzz0(i)-bzz0(i-1))/dr(i)
     $               +bth/rh*(bth0(i)*i-bth0(i-1)*(i-1))
     .               +(p0(i+1)-p0(i-1))/2./dr(i)
        enddo

        if (my_rank == 0) then
          write(*,*)
          write(*,'(a,1pe14.7)') ' Force-balance error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
        endif

        !Check Ohm's law
        do i=1,nxd
          bzh = 0.5*(bzz0(i)+bzz0(i-1))
          bth = 0.5*(bth0(i)+bth0(i-1))
          rh  = grid_params%xg(i)

          dummy(i) = Eoeta*bzh/etah(i)
     $              -bzh/rh*(bth0(i)*i-bth0(i-1)*(i-1))
     $              +bth*(bzz0(i)-bzz0(i-1))/dr(i)
        enddo

        if (my_rank == 0) then
          write(*,'(a,1pe14.7)') ' Ohms law      error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
        endif

        !Check energy equation
        etah = etah*etah0
        do i=1,nxd
          rh  = grid_params%xg(i)
          rhm = grid_params%xg(i-1)
          rhp = grid_params%xg(i+1)

! Daniele, 02/07/2014
          chirh=dfvty_rfx(rh)
          chirhm=dfvty_rfx(rhm)
          chirhp=dfvty_rfx(rhp)

! Daniele, 11/10/2013
          dummy(i) =
     .          gamma/(gamma-1.)
     .           *(rhp*p0(i+1)*vr0(i+1)-rhm*p0(i-1)*vr0(i-1))/(2.*dr(i))
     .          -vr0(i)*rh*(p0(i+1)-p0(i-1))/(2.*dr(i))
!     .          -chi/eta*(0.5*(rhp+rh)*(p0(i+1)-p0(i  ))/(rhp-rh)
!     .                   +0.5*(rhm+rh)*(p0(i  )-p0(i-1))/(rhm-rh))/dr(i)
! Daniele, 02/07/2014
     .          -1./eta*
     .          (0.5*(chirhp+chirh)
     .          *0.5*(rhp+rh)*(p0(i+1)-p0(i  ))/(rhp-rh)
     .          +0.5*(chirhm+chirh)
     .          *0.5*(rhm+rh)*(p0(i  )-p0(i-1))/(rhm-rh))
     .          /dr(i)
     .          -rh*etah(i)*(jtf(i)*jtf(i)+jzf(i)*jzf(i))
        enddo
! Daniele, 02/07/2014
        dummy = dummy*eta

        if (my_rank == 0) then
          write(*,'(a,1pe14.7)') ' Energy equation error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
        endif

        !Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)
        bz_avg = 0d0
        do i=1,nxd
          bz_avg = bz_avg + (bzz0(i)*i + bzz0(i-1)*(i-1))*dr(i)**2
        enddo

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' Pinch equilibrium features:'
          write (*,999) '   Theta  =',bth0(nxd)/bz_avg
          write (*,999) '   Bz(r=0)=',bzz0(0)
          write (*,999) '   I0     =',2*pi*xmax*bth0(nxd)
          write (*,999) '   E0/eta =',Eoeta
          write (*,999) '   T_flux =',pi*xmax**2*bz_avg
 999      format (a,f10.7)
        endif

        !Save some equilibrium properties (to calculate bootstrap current)
cc        eq_params(1) = 1d0/grid_params%params(3)       !Aspect ratio
cc        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
cc        eq_params(3) = q0h(nxd)/grid_params%params(3)  !q(r=a)

        !Find helical flux
        psi0(0) = 0d0
        do i=1,nxd+1
          psi0(i) = psi0(i-1)+dr(i-1)*(bth0(i-1)-kk*xh(i-1)*bzz0(i-1))
        enddo
        psi0 = psi0 - 0.5*(psi0(nxd)+psi0(nxdp)) !Homogeneous Dirichlet BC

        !Transfer B to integer mesh
        bth0 = btf
        bzz0 = bzf

101     format (5f7.3,2e12.3)

      end subroutine Ohm_equil_1d_xport

c     Ohm_hel_equil
c     ###############################################################
      subroutine Ohm_hel_equil(nx,ny,nz,kk,mm,Eoeta,theta,psig,bb,vv)

c     ---------------------------------------------------------------
c     Solves 2D helical Ohmic equil. equations to find cylindrical
c     components of magnetic field. On input:
c       * kk: helical parameter (n/m/R) (defined as theta-kk*z,
c             i.e, kk>0 for resonance)
c       * mm: poloidal perturbation number
c       * Eoeta: electric field over eta
c       * theta: RFP theta parameter
c       * psig: global 1D helical flux function (input)
c       * bb: on input, 1D equilibrium;
c             on output, helical eq. magnetic field components at
c             cell centers.
c       * vv: output poloidal flow at cell centers
c     In the calculation below, Psi and II are both at cell centers.
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use grid

        use app_iosetup

        use mg_solver

        use imposeBCinterface

        use local_BCS_variables, ONLY:order_bc

        use ts_setup, ONLY: vol_wgt

        use transport_params

        use orbit

#if defined(FPA)
        use fixed_point_accelerator
#endif
        implicit none

c     Call variables

        integer :: mm,nx,ny,nz
        real(8) :: bb(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,vv(0:nx+1,0:ny+1,0:nz+1,2)
     .            ,Eoeta,theta,kk

        real(8) :: psig(0:nxdp)

c     Local variables

        integer :: i,j,k,ii,jj,it,guess,bcnd(6)
     .            ,ip,im,jp,jm,kp,km,igrid,nn
        real(8) :: alpha,eps

        real(8) :: psi    (0:nx+1,0:ny+1,0:nz+1)
     .            ,III    (0:nx+1,0:ny+1,0:nz+1)
     .            ,dpsi   (0:nx+1,0:ny+1,0:nz+1)
     .            ,dIII   (0:nx+1,0:ny+1,0:nz+1)
     .            ,zeros  (0:nx+1,0:ny+1,0:nz+1)
     .            ,jpar_B (0:nx+1,0:ny+1,0:nz+1)
     .            ,ohm_err(0:nx+1,0:ny+1,0:nz+1)
     .            ,gs_err (0:nx+1,0:ny+1,0:nz+1)

        real(8) :: x(nx*ny*nz),r(nx*ny*nz)

#if defined(FPA)
        type (fpa_state) :: fpav
#endif
        logical :: phaseout_psibc=.false.

c     Begin program

        igrid = 1

        nn = nx*ny*nz

cc        have_jparB = .true.
cc
cc        allocate(jpar_B(0:nx+1,0:ny+1,0:nz+1))

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '***********************************'
          write (*,*) '2D Ohmic HELICAL equilibrium solver'
          write (*,*) '***********************************'
          write (*,*)
        endif

        if (Eoeta == 0d0) then
          messg = 'E0/eta = 0!!'
          call pstop('Ohm_hel_equil',messg)
        endif

        zeros = 0d0

        alpha = 0d0 !No underrelaxation
cc        alpha = 0.1d0 !underrelaxation constant

c     Initialize psi and B_3=III

        call init_hel_eq(psig,bb,psi,III)

c     Perturb Psi and magnetic field components with helical perturbation

        call perturb_1d_eq(bb,dpsi,dIII)

        psi = psi + dpsi
        III = III + dIII

c     Store BC information

        call allocateMGArray(1,gpsi0)
        call allocateMGArray(1,gbz0)

        gpsi0%grid(igrid)%array(:,:,:,1) = dpsi
        gbz0 %grid(igrid)%array(:,:,:,1) = dIII

        bcnd = bcond
        where (bcnd == DEF) bcnd = EXT

        call restrictMGArray(IRHO,1,gpsi0,bcnd,igrid,order_bc)
        call restrictMGArray(IRHO,1,gbz0 ,bcnd,igrid,order_bc)

c     Gather solution

        x = reshape(psi(1:nx,1:ny,1:nz),shape(x))

c     Iteration

#if defined(FPA)
        !Create fixed-point accelerator
        call fpa_create(fpav,x,maxv=5)
#endif

        do it=1,200

          call scatter_vec(x,psi)

          !Find new I=B_3
          III = solve_par_Ohm(psi)

          !Compute cnv magnetic field components
          bb = compute_B(psi,III)

          !Perform flux average
          jpar_B = flux_average(bb)

          !Dump plots (diag)
cc          call dumpPlots

          !Find new psi
          psi = solve_hel_GS(psi,III)
cc          psi = solve_hel_GS(zeros,III)

          !Compute residual vector
          r = x - reshape(psi(1:nx,1:ny,1:nz),shape(x))

#if defined(FPA)
          !Apply fixed-point accelerator
          call fpa_correction(fpav,r,dp=dot2)
cc          alpha = 0d0
#endif
          !Update solution
          x = x - (1d0-alpha)*r

          !Check convergence
          if (chk_conv(r)) exit

        enddo

#if defined(FPA)
        !Destroy fixed-point accelerator
        call fpa_destroy(fpav)
#endif

c     Compute converged magnetic field components

        call scatter_vec(x,psi)

        III = solve_par_Ohm(psi)

        bb = compute_B(psi,III)

c     Check force balance: j3 = -B^3 dI/d(psi) = B^3 j||/B

        call chk_hel_GS

c     Check PARALLEL Ohm's law:  j.B = dI/d(psi) * (B)^2

        call chk_par_Ohm

c     Dump final plots

        call dumpPlots

c     Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)

        call pinch_eq_features

c     Find poloidal flow (in cnv representation)

        vv = compute_pinch_flow(bb)

c     Deallocate memory

        call deallocateMGArray(gpsi0)
        call deallocateMGArray(gbz0)

      contains

c     init_hel_eq
c     ##################################################################
      subroutine init_hel_eq(psig,bb,psi,III)

        implicit none

c     ------------------------------------------------------------------
c     Initialize psi and III from 1D equilibrium quantities
c     ------------------------------------------------------------------

c     Call variables

        real(8) :: bb (0:nx+1,0:ny+1,0:nz+1,3)
     .            ,psi(0:nx+1,0:ny+1,0:nz+1)
     .            ,III(0:nx+1,0:ny+1,0:nz+1),psig(0:nxdp)

c     Local variables

        integer :: iglobal,ig,jg,kg
        real(8) :: x1,y1,z1

c     Begin program

c     Initialize I(psi)

        k=1 ; j=1
        do i=0,nx+1
          call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
     .                                  ,ig,jg,kg,x1,y1,z1)
          III(i,:,:) = bb(i,:,:,3)  + kk*x1*bb(i,:,:,2)
        enddo

c     Initialize Psi from 1D global psi

        do i=0,nx+1
          iglobal = i + grid_params%ilo(igrid) - 1
          psi(i,:,:) = psig(iglobal)
        enddo

c     End program

      end subroutine init_hel_eq

c     perturb_1d_eq
c     ##################################################################
      subroutine perturb_1d_eq(bb,dpsi,dIII)

        implicit none

c     ------------------------------------------------------------------
c     Initialize psi and III from 1D equilibrium quantities
c     ------------------------------------------------------------------

c     Call variables

        real(8) :: bb  (0:nx+1,0:ny+1,0:nz+1,3)
     .            ,dpsi(0:nx+1,0:ny+1,0:nz+1)
     .            ,dIII(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer :: ig,jg,kg!,bcnd(6)
        real(8) :: x1,y1,z1,eps

c     Begin program

        jpar_B = Eoeta*bb(:,:,:,3)/(bb(:,:,:,1)**2
     .                             +bb(:,:,:,2)**2
     .                             +bb(:,:,:,3)**2)

        eps = 1e-2

        k=1
        do j=0,ny+1
          do i=0,nx+1
            call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
     .                                    ,ig,jg,kg,x1,y1,z1)

cc            dpsi =-eps/mm*sin(pi*x1/grid_params%xg(nxd+1))*sin(y1*mm)  !Produces finite Br perturbation at
cc                                                                       !face, but cancels at ghost cells to
cc                                                                       !avoid field lines from getting out
cc                                                                       !of domain
            dpsi(i,j,:) = - eps/mm*x1*sin(y1*mm)   !Produces finite Br perturbation at wall
            dIII(i,j,:) = - jpar_B(i,j,:)*dpsi(i,j,:)
          enddo
        enddo

c     End program

      end subroutine perturb_1d_eq

c     scatter_vec
c     ##################################################################
      subroutine scatter_vec(x,psi)

        implicit none

c     ------------------------------------------------------------------
c     Map psi from interation vector
c     ------------------------------------------------------------------

c     Call variables

        real(8) :: x  (nx*ny*nz)
     .            ,psi(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

        integer :: bcnd(6),igr,ig,jg,kg
        real(8) :: lI0,I0,x1,y1,z1

        real(8) :: aa=0.3,eps=1e-1

c     Begin program

        psi(1:nx,1:ny,1:nz) = reshape(x(1:nn)
     .                               ,shape(psi(1:nx,1:ny,1:nz)))

        bcnd = bcond
        where (bcnd == DEF) bcnd = EQU
        call setBC(IRHO,nx,ny,nz,psi
     .            ,gpsi0%grid(igrid)%array(:,:,:,1)
     .            ,bcnd,igrid,igrid,igrid,iorder=order_bc)

c     Phase out psi BC

        if (phaseout_psibc) then
          aa = aa*0.5
          if (aa > (xmin-xmax)/nxd) then

            if (my_rank ==0) write (*,*) 'Phasing out psi BC...'

            do i=1,nx
              call getCurvilinearCoordinates(i,1,1,igrid,igrid,igrid
     .                                    ,ig,jg,kg,x1,y1,z1)
           
              psi(i,:,:) = psi(i,:,:)*tanh((1-x1)/eps+aa)
            enddo

            do igr=1,ngrid
              gpsi0%grid(igr)%array = aa*gpsi0%grid(igr)%array
            enddo
          else
            do igr=1,ngrid
              gpsi0%grid(igr)%array = 0d0
            enddo
          endif
        endif

c     End program

      end subroutine scatter_vec

c     flux_average
c     ##############################################################
      function flux_average(bb) result(ff)

        implicit none

c     --------------------------------------------------------------
c     Performs flux averages according to psi and flux_avg routine.
c     --------------------------------------------------------------

c     Call variables

        real(8) :: bb(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,ff(0:nx+1,0:ny+1,0:nz+1)

        INTERFACE
          function hel_flux_avg(fwd,xo,x,yo,y,zo,z,so,ss)
          logical :: hel_flux_avg,fwd
          real(8) :: xo,x,yo,y,zo,z,so,ss
          end function hel_flux_avg
        END INTERFACE

c     Local variables

        integer :: ig,jg,kg,ierr,bcnd(6)
        real(8) :: rh,th
        real(8),allocatable,dimension(:,:,:)   :: etag,jac,b2
        real(8),allocatable,dimension(:,:,:,:) :: xmap,bbg,bcar

c     Begin program

c     Initialize orbit integrator work space

        allocate(bbg  (0:nxdp,0:nydp,0:nzdp,3)
     .          ,bcar (0:nxdp,0:nydp,0:nzdp,3)
     .          ,jac  (0:nxdp,0:nydp,0:nzdp)
     .          ,xmap (0:nxdp,0:nydp,0:nzdp,3)
     .          ,etag (0:nxdp,0:nydp,0:nzdp))

c     Get GLOBAL contravariant magnetic field components

        call find_global(bb,bbg)

c     Orbit integrator setup: cnv components (with BCs)

        !Periodic BCs
        bbg(:,0    ,:,:)=bbg(:,nyd,:,:)
        bbg(:,nyd+1,:,:)=bbg(:,1  ,:,:)

        bbg(:,:,0    ,:)=bbg(:,:,nzd,:)
        bbg(:,:,nzd+1,:)=bbg(:,:,1  ,:)

        !SP BCs
        do k=0,nzd+1
          do j=0,nyd+1
            jj = mod(j+nyd/2,nyd)
            if (jj == 0) jj = nyd

            bbg(0,j,k,1) = bbg(1,jj,k,1)
            bbg(0,j,k,2) =-bbg(1,jj,k,2)
            bbg(0,j,k,3) =-bbg(1,jj,k,3)
          enddo
        enddo

        !Cartesian components, maps
        do k=0,nzd+1
          do j=0,nyd+1
            do i=0,nxd+1
              rh = grid_params%xg(i)

              !Cartesian components
              th = grid_params%yg(j)

              bcar(i,j,k,1) = bbg(i,j,k,1)/rh*cos(th)
     .                      -(bbg(i,j,k,2)+kk*bbg(i,j,k,3))*sin(th)
              bcar(i,j,k,2) = bbg(i,j,k,1)/rh*sin(th)
     .                      +(bbg(i,j,k,2)+kk*bbg(i,j,k,3))*cos(th)
              bcar(i,j,k,3) = bbg(i,j,k,3)/rh

              !Map (cylindrical)
              jac (i,j,k)  = rh
              xmap(i,j,k,1)= rh*cos(th)
              xmap(i,j,k,2)= rh*sin(th)
              xmap(i,j,k,3)= grid_params%zg(k)
            enddo
          enddo
        enddo

c     Orbit setup (input Az directly)

cc        bxg = 0d0
cc        byg = 0d0
cc        bzg = -psig
        call orbit_setup(nxd+2,nyd+2,nzd+2
     .                  ,grid_params%xg
     .                  ,grid_params%yg
     .                  ,grid_params%zg
     .                  ,bbg (:,:,:,1),bbg (:,:,:,2),bbg (:,:,:,3)
     .                  ,bcar(:,:,:,1),bcar(:,:,:,2),bcar(:,:,:,3)
     .                  ,jac,xmap,bcond
     .                  ,dtime=1d-2
     .                  ,solen=.true.
cc     $                  ,lsode_dt_lim=.true.
cc     .                  ,B_input_is_A=.true.
     .                  ,r_min=1d-3)

        !Orbit average (at cell centers)
        if (my_rank == 0.and.(ilevel > 1)) then
          write (*,FMT="(a)",ADVANCE="NO") 'Orbit averaging: '
        endif

c     Spline GLOBAL resistivity (normalized to nominal resistivity)

        do i=0,nxd+1
          etag(i,:,:) = res_rfx(grid_params%xg(i))
        enddo

        call splineFld(etag)

c     Orbit integration (local processor)

        ff = 0d0

        do k=1,nz
          do j=1,ny
            do i=1,nx
              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)
              x0 = grid_params%xx(ig)
              y0 = grid_params%yy(jg)
              z0 = grid_params%zz(kg)

              bzz_avg   = 0d0
              etab2_avg = 0d0
              ssmax     = 0d0

              call orbit_find(ilevel-3,x0,y0,z0,.true.
     .                       ,line_int=hel_flux_avg,ierror=ierr)

              select case(ierr)
              case(ORB_OK)
                if ((my_rank == 0).and.(ilevel > 1))
     .            write (*,FMT="(a)",ADVANCE="NO") '+'

                ff(i,j,k) = Eoeta*bzz_avg/etab2_avg
              case(ORB_OUT_DOM)
                if ((my_rank == 0).and.(ilevel > 1))
     .            write (*,FMT="(a)",ADVANCE="NO") 'X'

                ff(i,j,k) = quad_int(grid_params%xx(ig-4)
     .                              ,grid_params%xx(ig-3)
     .                              ,grid_params%xx(ig-2)
     .                              ,grid_params%xx(ig-1)
     .                              ,ff(i-4,j,k)
     .                              ,ff(i-3,j,k)
     .                              ,ff(i-2,j,k)
     .                              ,ff(i-1,j,k)
     .                              ,x0,2)
cc                ff(i,j,k) = quad_int(psi(i-4,j,k)
cc     .                              ,psi(i-3,j,k)
cc     .                              ,psi(i-2,j,k)
cc     .                              ,psi(i-1,j,k)
cc     .                              ,ff(i-4,j,k)
cc     .                              ,ff(i-3,j,k)
cc     .                              ,ff(i-2,j,k)
cc     .                              ,ff(i-1,j,k)
cc     .                              ,psi(i,j,k),2)
              case default
                call pstop('flux_average'
     $                    ,'Orbit error '//int2char(ierr)//' at '
     $                     //int2char(i)//','//int2char(j))
              end select

            enddo
          enddo
        enddo

c     Impose BCs

        bcnd = bcond
        where (bcnd == DEF) bcnd = EXT
        call setBC(IRHO,nx,ny,nz,ff,zeros,bcnd,igrid,igrid,igrid
     .            ,iorder=order_bc)
cc     .            ,iorder=2)

c     Free up orbit integral work space

        if ((my_rank == 0).and.(ilevel > 1)) write (*,*)

        call orbit_destroy

        deallocate(bbg,bcar,jac,xmap)

      end function flux_average

c     solve_hel_GS
c     #################################################################
      function solve_hel_GS(psi,III) result(psi_new)

        implicit none

c     -----------------------------------------------------------------
c     Solves helical Grad-Shafranov equation
c     -----------------------------------------------------------------

c       Call variables

        real(8) :: psi    (0:nx+1,0:ny+1,0:nz+1)
     .            ,psi_new(0:nx+1,0:ny+1,0:nz+1)
     .            ,III    (0:nx+1,0:ny+1,0:nz+1)
     .            ,psiv(nx*ny*nz)

c       Local variables

        integer :: ig,jg,kg
        real(8) :: x1,y1,z1

        real(8) :: rhs (nn)
        integer :: bcnd(6)

        external del_hel

c       Begin program

        do k=1,nz
          do j=1,ny
            do i=1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
     .                                      ,ig,jg,kg,x1,y1,z1)

              rhs(ii)=-III(i,j,k)/(1+(kk*x1)**2)
     .                *(2*kk/(1+(kk*x1)**2)-jpar_B(i,j,k))

              if (vol_wgt) rhs(ii)=rhs(ii)
     .                            *gmetric%grid(igrid)%dvol(i,j,k)
            enddo
          enddo
        enddo

c       Recalibrate psi BC

c       Solve Helical Grad-Shafranov equation (find psi at cell centers)

        bcnd = bcond
        where (bcnd == DEF) bcnd = EQU

        psiv = reshape(psi(1:nx,1:ny,1:nz),shape(psiv))

        guess = 1
        call cSolver(1,nn,rhs,psiv,bcnd
     .              ,igrid,ilevel-1,guess,del_hel,vol_wgt
     .              ,tol           =1d-7
     .              ,gm_driver     =.true.
     .              ,ks_it         =20
     .              ,mg_vcyc       = 1
     .              ,mg_order_res  = 0
     .              ,mg_order_prol = 2
     .              ,mg_smooth     = 'jb'
     .              ,mg_gm_coarse_solve  =.false.
     .              ,sm_omega      = 0.7d0
     .              ,sm_ncolors    = 2  !Star stencil
     .              ,sm_zebra_relax=.false.
     .              )

        !Impose BCs
        psi_new(1:nx,1:ny,1:nz)= reshape(psiv
     .                                  ,shape(psi_new(1:nx,1:ny,1:nz)))

        bcnd = bcond
        where (bcnd == DEF) bcnd = EQU
        call setBC(IRHO,nx,ny,nz,psi_new
     .            ,gpsi0%grid(igrid)%array(:,:,:,1)
     .            ,bcnd,igrid,igrid,igrid,iorder=order_bc)

      end function solve_hel_GS

c     solve_par_Ohm
c     #####################################################################
      function solve_par_Ohm(psi) result(III)

        implicit none

c     -----------------------------------------------------------------
c       Solves parallel Ohm's law jpar/B=f(psi)
c     -----------------------------------------------------------------

c       Call variables

        real(8) :: psi(0:nx+1,0:ny+1,0:nz+1)
     .            ,III(0:nx+1,0:ny+1,0:nz+1)

c       Local variables

        integer :: bcnd(6)
        real(8) :: rhs(nn),IIv(nx*ny*nz)
        real(8) :: lI0,I0

        external lap_mtvc2

c       Begin program

c       Find Laplacian RHS (at cell centers, for computation of I=B_3)

        rhs = 0d0
        do k=1,nz
          do j=1,ny
            do i=1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)
              rhs(ii)=-lap(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                          ,psi,dff=jpar_B,vol=vol_wgt)
            enddo
          enddo
        enddo

c       Solve for B_3=I (at cell centers)

        bcnd = bcond
        where (bcnd == DEF) bcnd = EQU
cc        where (bcnd == DEF) bcnd = DIR

        IIv = 0d0

        guess = 0
        call cSolver(1,nn,rhs,IIv,bcnd
     .              ,igrid,ilevel-1,guess,lap_mtvc2,vol_wgt
     .              ,tol           = 1d-7
     .              ,gm_driver     = .true.
     .              ,ks_it         = 20
     .              ,mg_order_res  = 0
     .              ,mg_order_prol = 2
     .              ,mg_smooth     = 'jb'
     .              ,mg_gm_coarse_solve =.true.
     .              ,mg_vcyc       = 1
     .              ,sm_omega      = 0.7d0
     .              ,sm_ncolors    = 4  !9-pt stencil
     .              ,sm_zebra_relax=.false.
     .              )

        !Set BCs
        III(1:nx,1:ny,1:nz) = reshape(IIv,shape(III(1:nx,1:ny,1:nz)))

        call setBC(IRHO,nx,ny,nz,III,gbz0%grid(igrid)%array(:,:,:,1)
     .            ,bcnd,igrid,igrid,igrid,iorder=order_bc)

c       Renormalize III so that Bz(r=0)=1d0

        if (isBdry(1,igrid,1)) then
          I0 = sum(III(1,1:ny,1))/nyd
        else
          I0 = 0d0
        endif

#if defined(petsc)
        lI0=I0
        call MPI_Allreduce(lI0,I0,1,MPI_DOUBLE_PRECISION
     .       ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

        III = III - I0 + 1d0

      end function solve_par_Ohm

c     compute_B
c     ##################################################################
      function compute_B(psi,III) result(bb)

      implicit none

c     ------------------------------------------------------------------
c     Finds B in cnv representation
c     ------------------------------------------------------------------

c     Call variables

        real(8) :: psi(0:nx+1,0:ny+1,0:nz+1)
     .            ,III(0:nx+1,0:ny+1,0:nz+1)
     .            ,bb (0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

        integer :: ig,jg,kg,bcnd(6,3)
        real(8) :: x1,y1,z1,gr(3)

c     Begin program

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
     .                                      ,ig,jg,kg,x1,y1,z1)

              gr = grad_ijk(i,j,k,nx,ny,nz,igrid,igrid,igrid,psi)

cc              !Cylindrical
cc              bb(i,j,k,1) = -gr(2)/x1
cc              bb(i,j,k,2) = (gr(1)+kk*x1*III(i,j,k))/(1+(kk*x1)**2)
cc              bb(i,j,k,3) = III(i,j,k) - kk*x1*bb(i,j,k,2)
              !Contravariant
              bb(i,j,k,1) = -gr(2)
              bb(i,j,k,2) =  gr(1)
              bb(i,j,k,3) = x1*(III(i,j,k)-kk*x1*gr(1))/(1+(kk*x1)**2)
            enddo
          enddo
        enddo

        !Impose topological BCs
        bcnd(:,1) = bcond
        bcnd(:,2) = bcond
        bcnd(:,3) = bcond
        where (bcnd == DEF) bcnd = EXT
        call setMGBC(0,3,nx,ny,nz,igrid,bb,bcnd
     .              ,icomp=(/IBX/),is_cnv=.true.,is_vec=.true.
     $              ,iorder=order_bc)

      end function compute_B

c     chk_conv
c     ##################################################################
      function chk_conv(dx)

        implicit none

c     ------------------------------------------------------------------
c     Finds B in cnv representation
c     ------------------------------------------------------------------

c     Call variables

        logical :: chk_conv
        real(8) :: dx(:)

c     Local variables

        real(8) :: dh2,err,tol

c     Begin program

        dh2 = ((xmax-xmin)/nxd)**2
cc        tol = 0.0001*dh2
        tol = 0.01*dh2

        err = sqrt(dot2(dx,dx))

        if (ilevel > 0 .and. my_rank == 0) then
          write (*,111) '2D_Ohmic_equilibrium_iter=',it
     .                 ,'   Error=',err
 111      format (a,i5,a,1pe10.3)
        endif

        chk_conv = (err < tol)

cc        phaseout_psibc = (err < 100*tol).and.(err > tol)

        if (chk_conv) then
          if (ilevel == 0 .and. my_rank == 0.and.chk_conv) then
            write (*,'(a,i3,a,1pe10.2)')
     .          ' 2D equilibrium converged in '
     .          ,it,' iterations with error =',err
          endif
        endif

      end function chk_conv

c     chk_hel_GS
c     ##################################################################
      subroutine chk_hel_GS

      implicit none

c     ------------------------------------------------------------------
c     Checks residual of helical GS equation from magnetic field
c     components
c     ------------------------------------------------------------------

c     Call variables

c     Local variables

        integer :: ig,jg,kg
        real(8) :: x1,y1,z1,j3,err

c     Begin program

        gs_err = 0d0

        do k=1,nz
          do j=1,ny-1
            do i=1,nx
              call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
     .                                    ,ig,jg,kg,x1,y1,z1)

              j3=(grid_params%xx(ig+1)*(bb(i+1,j,k,2)+kk*bb(i+1,j,k,3))
     .           -grid_params%xx(ig-1)*(bb(i-1,j,k,2)+kk*bb(i-1,j,k,3)))
     .           /2./grid_params%dxh(ig)  ! d(B_2)/dx1
     .          -(bb(i,j+1,k,1)-bb(i,j-1,k,1))/x1/2./grid_params%dyh(jg)
                                          !-d(B_1)/dx2

              gs_err(i,j,k) = j3-bb(i,j,k,3)*jpar_B(i,j,k)
            enddo
          enddo
        enddo

        err = integral(igrid,igrid,igrid,nx,ny,nz,gs_err**2)
        err = sqrt(err)

        if (my_rank == 0) then
          write(*,*)
          write(*,'(a,1pe14.7)') ' Force-balance error =',err
        endif

      end subroutine chk_hel_GS

c     chk_par_Ohm
c     ##################################################################
      subroutine chk_par_Ohm

      implicit none

c     ------------------------------------------------------------------
c     Checks residual of parallel Ohm's law
c     ------------------------------------------------------------------

c     Call variables

c     Local variables

        integer :: ig,jg,kg
        real(8) :: x1,y1,z1,j1,j2,j3,b1,b2,b3,br,bth,bz,err,bnorm2

c     Begin program

        ohm_err = 0d0

        do k=1,nz
          do j=1,ny-1
            do i=1,nx
              call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
     .                                      ,ig,jg,kg,x1,y1,z1)

              !Cyl B
              br = bb(i,j,k,1)/x1
              bth= bb(i,j,k,2)+kk*bb(i,j,k,3)
              bz = bb(i,j,k,3)/x1

              bnorm2 = br**2+bth**2+bz**2

              !Cnv J
              j1= (III(i,j+1,k)-III(i,j-1,k))/2./grid_params%dyh(jg) ! d(B3)/dx2
              j2=-(III(i+1,j,k)-III(i-1,j,k))/2./grid_params%dxh(ig) !-d(B3)/dx1
              j3=(grid_params%xx(ig+1)*(bb(i+1,j,k,2)+kk*bb(i+1,j,k,3))
     .           -grid_params%xx(ig-1)*(bb(i-1,j,k,2)+kk*bb(i-1,j,k,3)))
     .           /2./grid_params%dxh(ig)  ! d(B_2)/dx1
     .          -(bb(i,j+1,k,1)-bb(i,j-1,k,1))/x1/2./grid_params%dyh(jg)
     .                                    !-d(B_1)/dx2

              !Cov B
              b1 = br
              b2 = x1*bth
              b3 = III(i,j,k)

              ohm_err(i,j,k)=(j1*b1+j2*b2+j3*b3)/x1-jpar_B(i,j,k)*bnorm2
            enddo
          enddo
        enddo

        err = integral(igrid,igrid,igrid,nx,ny,nz,ohm_err**2)
        err = sqrt(err)

        if (my_rank == 0) then
          write(*,'(a,1pe14.7)') ' Parallel Ohms law error =',err
        endif

      end subroutine chk_par_Ohm

c     pinch_eq_features
c     ##################################################################
      subroutine pinch_eq_features

      implicit none

c     ------------------------------------------------------------------
c     Compute's helical PINCH equilibrium features
c     ------------------------------------------------------------------

c     Call variables

c     Local variables

        integer :: ig,jg,kg
        real(8) :: x1,y1,z1,bz_avg,bth_avg,lpsi0
     .            ,dummy(0:nx+1,0:ny+1,0:nz+1)

c     Begin program

        !Average Bz in domain
        do i=1,nx
          call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
     .                                  ,ig,jg,kg,x1,y1,z1)
          dummy(i,:,:) = bb(i,:,:,3)/x1
        enddo

        bz_avg = integral(igrid,igrid,igrid,nx,ny,nz,dummy
     .                   ,average=.true.)

        !Average Bth @ boundary
        if (isBdry(nx,igrid,2)) then
          bth_avg = sum(bb(nx,1:ny,1,2)+kk*bb(nx,1:ny,1,3))/nyd
        else
          bth_avg = 0d0
        endif

#if defined (petsc)
        lpsi0 = bth_avg
        call MPI_Allreduce(lpsi0,bth_avg,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' Pinch equilibrium features:'
          write (*,999) '   Theta  =',bth_avg/bz_avg
          write (*,999) '   Bz(r=0)=',dummy(1,1,1)
          write (*,999) '   I0     =',2*pi*xmax*bth_avg
          write (*,999) '   E0/eta =',Eoeta
          write (*,999) '   T_flux =',pi*xmax**2*bz_avg
 999      format (a,f10.7)
        endif

      end subroutine pinch_eq_features

c     dumpPlots
c     ##################################################################
      subroutine dumpPlots

      implicit none

c     ------------------------------------------------------------------
c     Dump plots
c     ------------------------------------------------------------------

c     Call variables

c     Local variables

      real(8) :: dum(0:nx+1,0:ny+1,0:nz+1)
     $          ,dummy(0:nxd+1,0:nyd+1,0:nzd+1)

c     Begin program

      call createDrawInCfile(6,'hel_eq.bin','Helical Eq.','t','x','y'
     $        ,(/'A_3    ','B_3    ','j||/B  '
     $          ,'Resid  ','GS err ','Ohm err'/)
     $        ,'-c -X0 -L57','drawhel_eq.in')

      if (my_rank == 0) then
        open(unit=110,file='hel_eq.bin',form='unformatted'
     .      ,status='replace')
      endif

      !Psi
      call find_global(psi,dummy)
      if (my_rank == 0)
     $  call contour(dummy(0:nxd+1,0:nyd,1),nxd+2,nyd+1
     .              ,0d0,xmax,0d0,ymax,0,110)

      !III
      call find_global(III,dummy)
      if (my_rank ==0)
     $  call contour(dummy(0:nxd+1,0:nyd,1),nxd+2,nyd+1
     .              ,0d0,xmax,0d0,ymax,1,110)

cc      dIII=0d0
cc      dIII=jpar_B*sqrt(vectorNorm(igrid,bb,.false.))
cc      call contour(dIII (0:nxd+1,0:nyd,1),nxd+2,nyd+1
      !j||/B
      call find_global(jpar_B,dummy)
      if (my_rank == 0)
     $  call contour(dummy(0:nxd+1,0:nyd,1),nxd+2,nyd+1
     .              ,0d0,xmax,0d0,ymax,1,110)

      !Iteration residual
      dum=0d0
      dum(1:nx,1:ny,1:nz) = reshape(r,shape(psi(1:nx,1:ny,1:nz)))
      call find_global(dum,dummy)
      if (my_rank == 0)
     $  call contour(dummy(0:nxd+1,0:nyd,1),nxd+2,nyd+1
     .              ,0d0,xmax,0d0,ymax,1,110)

      !GS error
      call find_global(gs_err,dummy)
      if (my_rank == 0)
     $  call contour(dummy(0:nxd+1,0:nyd,1),nxd+2,nyd+1
     .              ,0d0,xmax,0d0,ymax,1,110)

      !Ohm's law error
      call find_global(ohm_err,dummy)
      if (my_rank == 0)
     $  call contour(dummy(0:nxd+1,0:nyd,1),nxd+2,nyd+1
     .              ,0d0,xmax,0d0,ymax,1,110)

      if (my_rank == 0) close(110)

      end subroutine dumpPlots

c     compute_pinch_flow
c     ##################################################################
      function compute_pinch_flow(bb) result(vv)

        implicit none

c     ------------------------------------------------------------------
c     Compute's helical PINCH equilibrium features
c     ------------------------------------------------------------------

c     Call variables

        real(8) :: vv(0:nx+1,0:ny+1,0:nz+1,2)
     .            ,bb(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

        integer :: ig,jg,kg
        real(8) :: x1,y1,z1,br,bth,bz,bnorm2,b1,b2

c     Begin program

        do k=1,nz
          do j=1,ny
            do i=1,nx
              call getCurvilinearCoordinates(i,j,k,igrid,igrid,igrid
     .                                      ,ig,jg,kg,x1,y1,z1)

              !Cyl B
              br = bb(i,j,k,1)/x1
              bth= bb(i,j,k,2)+kk*bb(i,j,k,3)
              bz = bb(i,j,k,3)/x1

              bnorm2 = br**2+bth**2+bz**2

              !Cov B
              b1 = br
              b2 = grid_params%xx(ig)*bth

              vv(i,j,k,1) = -(Eoeta*b2)/bnorm2
              vv(i,j,k,2) =  (Eoeta*b1)/bnorm2

            enddo
          enddo
        enddo

      end function compute_pinch_flow

      end subroutine Ohm_hel_equil

c     hel_flux_avg
c     ###################################################################
      function hel_flux_avg(fwd,xo,x,yo,y,zo,z,so,ss) result(exloop)

c     -------------------------------------------------------------------
c     Performs flux average of various quantities for 2D helical
c     equilibrium
c     -------------------------------------------------------------------

        use equilibrium

        use par_int

        implicit none

c     Call variables

        logical :: exloop,fwd
        real(8) :: xo,x,yo,y,zo,z,so,ss

c     Local variables

        integer :: ierr
        real(8) :: xh,yh,zh,bb,b2,bx,by,bz
     .            ,xs,ys,zs,ts,ds,eta,jac

c     Begin program

        ierr = ORB_OK

        if (orb_init) call per_orbit_init(xo,x,yo,y,zo,z,order=1)

c     Check whether orbit should terminate

        call per_orbit_chk(so,ss,xo,x,yo,y,zo,z,ts,xs,ys,zs) !--> Returns ts <> 0 if intersects

        exloop = (ts /= 0d0)

c     Update current position

        if (exloop) then
          x = xs ; y = ys ; z = zs
        else
          ts = ss
        endif

c     Find middle points

        xh = 0.5*(x +xo)
        yh = 0.5*(y +yo)
        zh = 0.5*(z +zo)

c     Find magnetic field Cartesian components

        call getB(xh,yh,zh,bx,by,bz,.false.,.true.,ierr)

        if (ierr /= ORB_OK) return

        bb = sqrt(bx**2+by**2+bz**2)

        eta = evalFld(xh,yh,zh,ierr)

c     Perform flux averages

        ds = (ts-so)  !Arc length

        bzz_avg   = bzz_avg    + bz*ds/bb
        etab2_avg = etab2_avg  + eta*bb*ds
        ssmax     = ssmax      + ds/bb
cc        b2_avg  = b2_avg  + bb*ds
cc        ssmax   = ssmax   + ds/bb*b2  !Projected poloidal arc length

      end function hel_flux_avg

c     Ohm_GS_equil_2d
c     ###############################################################
      subroutine Ohm_GS_equil_2d(R0,Eoeta,b_tor,nx,ny,nz,b1,b2,b3
     .                          ,v1,v2,pp)

c     ---------------------------------------------------------------
c     Solves 2D Ohmic Grad-Shafranov equilibrium equations to find 
c     contravariant components of magnetic field. On input:
c       * R0: Major radius
c       * Eoeta: electric field over eta
c       * b_tor: toroidal beta
c       * nx,ny,nz: local (processor) dimensions
c       * b1,b2,b3: output magnetic field contravariant components 
c       * v1,v2: output velocity contravariant components
c
c     In the calculation below, psi is at cell centers, and B_3=I(psi)
c     is defined at faces.
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use grid

        use grid_debug

        use app_iosetup

        use PC_var, ONLY: csolver

        use mg_xfer, ONLY: mapMGVectorToArray

        use imposeBCinterface

        use ts_setup, ONLY: vol_wgt

        use operators, ONLY: lap

        use transport_params

        use nlfunction_setup, ONLY: res_rfx,eta

        use local_BCS_variables, ONLY: order_bc

        implicit none

c     Call variables

        integer :: nx,ny,nz
        real(8) :: b1(0:nx+1,0:ny+1,0:nz+1)
     .            ,b2(0:nx+1,0:ny+1,0:nz+1)
     .            ,b3(0:nx+1,0:ny+1,0:nz+1)
     .            ,v1(0:nx+1,0:ny+1,0:nz+1)
     .            ,v2(0:nx+1,0:ny+1,0:nz+1)
     .            ,pp(0:nx+1,0:ny+1,0:nz+1)
     .            ,Eoeta,R0,b_tor

c     Local variables

        integer :: i,j,k,ig,jg,kg,ii,it,bcnd(6),igrid,guess
     .            ,iglobal,jglobal,pos(2),nn
     .            ,np_RAD,my_rank_RAD,np_Z,my_rank_Z

        real(8) :: bz_avg,bth_avg,gserr,oherr,rh,dh,dt,I0,bnorm2
     .            ,psir,psith,x1,y1,z1,IIh,spsi,psi0,psia,gpsi2
     .            ,cc,p00,etah,lpsi0,dsdpsi,dpdpsi,spsip,spsim
     .            ,dpdr,dpdth,f_1,f_2,f_3,j1,j2,j3,jac,e1,e2,e3
     .            ,psip,psim,q0,qa,tor_flx,bz2_avg,p_avg,b_pol
     .            ,cov(3),e4,e5,e6

        real(8) :: xh(0:nxd),dr(0:nxd),dth(0:nyd)

        real(8) :: iR2  (0:nx+1,0:ny+1,0:nz+1)
     .            ,zeros(0:nx+1,0:ny+1,0:nz+1)
     .            ,qprof(0:nx+1,0:ny+1)
     .            ,dummy(0:nx+1,0:ny+1)
     .            ,III  (0:nx+1,0:ny+1)
     .            ,RR   (0:nx+1,0:ny+1)
     .            ,psi(0:nx+1,0:ny+1,0:nz+1)

        real(8),allocatable,dimension(:) :: psi_rhs,psiv,errv
     .                                     ,mpivec,mpivec_in

        real(8) :: pprof
        external del_star,pprof

c     Begin program

        igrid = 1

        allocate(psi_rhs(nx*ny*nz),psiv(nx*ny*nz),errv(nx*ny*nz))

        psi_rhs = 0d0 ; psiv = 0d0 ; errv = 0d0

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '************************************'
          write (*,*) '2D Ohmic GS equilibrium solver'
          write (*,*) '************************************'
          write (*,*)
        endif

        if (Eoeta == 0d0) then
          messg = 'E0/eta = 0!!'
          call pstop('Ohm_GS_equil_2d',messg)
        endif

        zeros = 0d0 ; iR2 = 0d0 ; dummy = 0d0
        qprof = 0d0 ; III = 0d0 !; psig = 0d0

c       Initialize iteration

        !Find half-mesh positions
        do i=0,nxd
          xh(i) = 0.5*(grid_params%xg(i+1)+grid_params%xg(i))
          dr(i) =      grid_params%xg(i+1)-grid_params%xg(i)
        enddo

        do j=0,nyd
          dth(j)=      grid_params%yg(j+1)-grid_params%yg(j)
        enddo

        !Initial guess for psi
        psi = 0d0

        !Initalize toroidal radius (R=R0+r*cos(theta)) and 1/R^2
        do j=0,ny+1
          do i=0,nx+1
            iglobal = i + grid_params%ilo(igrid) - 1
            jglobal = j + grid_params%jlo(igrid) - 1
            RR(i,j) = R0 + grid_params%xg(iglobal)
     .                    *cos(grid_params%yg(jglobal))
          enddo
        enddo

        if (minval(RR) < 0d0) then
          call pstop('Ohm_GS_equil_2d','r > R in torus')
        endif

        do k=0,nz+1
          iR2(:,:,k) = RR**(-2)
        enddo

        !Initial guess for III=B_phi*R
        III = RR   !B_phi=1 uniform

        !Maximum pressure
        p00 = 0.5*b_tor

c       Iteration

        psi0 =-1d0
        psia = 0d0
        oherr= 0d0

        do it=1,100

          !Find GS equation RHS (at cell centers)
          do k=1,nz
            do j=1,ny
              do i=1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              rh = grid_params%xx (ig)
              dh = grid_params%dxh(ig)
              dt = grid_params%dyh(jg)

              IIh = 0.5*(III(i,j)+III(i-1,j))

              psir   = 0.5*(psi(i+1,j,k)-psi(i-1,j,k))/dh
              psith  = 0.5*(psi(i,j+1,k)-psi(i,j-1,k))/dt/rh
              gpsi2  = psir**2 + psith**2
              bnorm2 = IIh**2 + gpsi2  !This is R^2*B^2

              etah =res_rfx(rh)

              spsip = sqrt(1d0-psi(i+1,j,k)/psi0)
              spsim = sqrt(1d0-psi(i-1,j,k)/psi0)
              dpdr  = 0.5*p00*(pprof(spsip,0)-pprof(spsim,0))/dh

              spsip = sqrt(1d0-psi(i,j+1,k)/psi0)
              spsim = sqrt(1d0-psi(i,j-1,k)/psi0)
              dpdth = 0.5*p00*(pprof(spsip,0)-pprof(spsim,0))/dt/rh

              psi_rhs(ii)= (Eoeta/etah*IIh**2*iR2(i,j,k)
     .                     -psir*dpdr-psith*dpdth)/bnorm2

              !Multiply by R to ensure proper limit R->inf
              errv(ii) = (lap(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                             ,psi,dff=iR2,vol=.false.)
     .                   -psi_rhs(ii))*RR(i,j)

              if (vol_wgt) psi_rhs(ii)=psi_rhs(ii)
     .                                *gmetric%grid(igrid)%dvol(i,j,k)

              oherr = oherr
     .               +(Eoeta/etah*IIh*iR2(i,j,k)*psir
     .                +IIh*dpdr
     .                +bnorm2*iR2(i,j,k)  !R^2 factor needed to normalize bnorm2
     .                *(III(i,j)-III(i-1,j))/dh)**2
              enddo
            enddo
          enddo

          !Calculate L2 norm of error
          gserr = sqrt(dot2(errv,errv)/nxd/nyd/nzd)

#if defined(petsc)
          e1 = oherr
          call MPI_Allreduce(e1,oherr,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
          oherr = sqrt(oherr/nxd/nyd/nzd)

          !Solve GS Poisson eq (find psi at cell centers)
          !Define BCs
          bcnd = bcond
          where (bcnd == DEF) bcnd = DIR

          guess = 1  !To reuse previous iteration

          call cSolver(1,nx*ny*nz,psi_rhs,psiv,bcnd
     .              ,1,ilevel,guess,del_star,vol_wgt
     .              ,tol           = 5d-1
     .              ,gm_driver     =.true.
     .              ,ks_it         = 10
     .              ,mg_order_res  = 0
     .              ,mg_order_prol = 1
     .              ,mg_gm_coarse_solve=.false.
cc     .              ,mg_gm_coarse_solve=.true.
     .              ,mg_vcyc       = 1
     .              ,mg_smooth     = 'jb'
     .              ,sm_ncolors    = 2  !Star stencil
     .              ,sm_it         = 5
     .              )

          !Map to array (w BCs)
          call mapMGVectorToArray(0,1,psiv,nx,ny,nz,psi,igrid,.false.)

          call setBC(IRHO,nx,ny,nz,psi,zeros,bcnd,igrid,igrid,igrid
     .              ,iorder=order_bc)

          !Perform integral of III (at half-mesh on global grid in radius)

          psi0 = minval(psi)
#if defined (petsc)
          lpsi0 = psi0
          call MPI_Allreduce(lpsi0,psi0,1,MPI_DOUBLE_PRECISION
     .                      ,MPI_MIN,MPI_COMM_WORLD,mpierr)
#endif

          k = 1
          do j=1,ny
            dummy(0,j) = 0d0
            do i=1,nx
              IIh = 0.5*(III(i,j)+III(i-1,j))

              call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

              rh = grid_params%xx(ig)
              dh = grid_params%dxh(ig)
              dt = grid_params%dyh(jg)

              psir  =0.5*(psi(i+1,j  ,k)
     .                   -psi(i-1,j  ,k))/dh
              psith =0.5*(psi(i  ,j+1,k)
     .                   -psi(i  ,j-1,k))/dt/rh
              bnorm2= IIh**2 + psir**2 + psith**2

              psip  = 0.5*(psi(i+1,j,k)
     .                    +psi(i  ,j,k))
              psim  = 0.5*(psi(i-1,j,k)
     .                    +psi(i  ,j,k))
              spsip = sqrt(1d0-psip/psi0)
              spsim = sqrt(1d0-psim/psi0)
              dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/dh

              etah  = res_rfx(rh)

              dummy(i,j) = dummy(i-1,j)
     .                         - RR(i,j)**2/bnorm2
     .                          *(Eoeta*psir/etah*iR2(i,j,k)
     .                           +dpdr)*grid_params%dx(ig)!!*dr(iglobal)
            enddo
          enddo

#if defined(petsc)
          !Broadcast poloidal integral along Z
          call initMPI(MPI_COMM_Z,np_Z,my_rank_Z)

          nn = size(dummy)
          call MPI_Bcast(dummy,nn,MPI_DOUBLE_PRECISION
     .                  ,0,MPI_COMM_Z,mpierr)

          !Communicate radial integration constants
          if (npx > 1) then
            allocate(mpivec(npx*ny))
            call MPI_Allgather(dummy(nx,1:ny)
     .                             ,ny,MPI_DOUBLE_PRECISION
     .                      ,mpivec,ny,MPI_DOUBLE_PRECISION
     .                      ,MPI_COMM_RAD,mpierr)

            !Sum integration constants
            do j=2,npx
              mpivec(1+ny*(j-1):ny*j) = mpivec(1+ny*(j-1):ny* j   )
     .                                 +mpivec(1+ny*(j-2):ny*(j-1))
            enddo

            !Add integration constants
            call initMPI(MPI_COMM_RAD,np_rad,my_rank_rad)

            j = my_rank_rad
            if (j > 0) then
              do i=0,nx+1
                dummy(i,1:ny) = dummy(i,1:ny) + mpivec(1+ny*(j-1):ny*j)
              enddo
            endif

            deallocate(mpivec)
          endif
#endif

          !Find III (local array)
          III = R0*exp(dummy)

          !Check convergence
          if (ilevel > 0 .and. my_rank == 0) then
            write (*,'(a,i3,a,1pe10.2,a,e10.2)')
     .                '2D Ohmic equilibrium iter =',it
     .                 ,' ; GS error=',gserr,' Ohm error=',oherr
            if (ilevel > 1) write (*,*)
          endif

          dh = xmax/nxd
          if (gserr < 0.001*dh**2) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
           write (*,'(a,i3,a,1pe10.2)')
     .          ' 2D equilibrium converged in '
     .          ,it,' iterations with error =',gserr
        endif

c       Find vectors and scalars

        k = 1
        do j=1,ny
          do i=1,nx
            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            rh = grid_params%xx(ig)
            dh = grid_params%dxh(ig)
            dt = grid_params%dyh(jg)

            IIh    = 0.5*(III(i,j)+III(i-1,j))
            psir   = 0.5*(psi(i+1,j,k)-psi(i-1,j,k))/dh
            psith  = 0.5*(psi(i,j+1,k)-psi(i,j-1,k))/dt/rh
            bnorm2 = (IIh**2+psir**2+psith**2)*iR2(i,j,k)

           !Find contravariant magnetic field components in integer mesh
            b1(i,j,:) =-psith*rh
            b2(i,j,:) = psir
            b3(i,j,:) = IIh*rh/RR(i,j)

            !Find contravariant velocity components
            etah =res_rfx(rh)

            spsip = sqrt(1d0-psi(i+1,j,k)/psi0)
            spsim = sqrt(1d0-psi(i-1,j,k)/psi0)
            dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dh

            spsip = sqrt(1d0-psi(i,j+1,k)/psi0)
            spsim = sqrt(1d0-psi(i,j-1,k)/psi0)
            dpdth = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dt/rh

            v1(i,j,:) =-(Eoeta*psir *iR2(i,j,k)+etah*dpdr )
     .                  *RR(i,j)*rh/bnorm2
            v2(i,j,:) =-(Eoeta*psith*iR2(i,j,k)+etah*dpdth)
     .                  *RR(i,j)/bnorm2

            !Find pressure
            spsi = sqrt(1d0-psi(i,j,k)/psi0)
            pp(i,j,:) = p00*pprof(spsi,0)
          enddo
        enddo

c       Calculate GS convergence properties

        e1 = 0d0
        e2 = 0d0
        e3 = 0d0
        e4 = 0d0
        e5 = 0d0
        e6 = 0d0

        k = 1
        do j=2,ny-1
          do i=2,nx-1
cc        do j=1,ny
cc          do i=1,nx
            iglobal = i + grid_params%ilo(igrid) - 1
            jglobal = j + grid_params%jlo(igrid) - 1

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            rh = grid_params%xx(ig)
            dh = grid_params%dxh(ig)
            dt = grid_params%dyh(jg)

            IIh    = 0.5*(III(i,j)+III(i-1,j))
            psir   = (psi(i+1,j,k)-psi(i-1,j,k))/2./dh
            psith  = (psi(i,j+1,k)-psi(i,j-1,k))/2./dt/rh
            bnorm2 = (IIh**2+psir**2+psith**2)*iR2(i,j,k)

            jac = rh*RR(i,j)

            !Find pressure derivatives
            spsip = sqrt(1d0-psi(i+1,j,k)/psi0)
            spsim = sqrt(1d0-psi(i-1,j,k)/psi0)
            dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dh

            spsip = sqrt(1d0-psi(i,j+1,k)/psi0)
            spsim = sqrt(1d0-psi(i,j-1,k)/psi0)
            dpdth = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dt

            !Find contravariant current components
            j1 = (III(i,j+1)+III(i-1,j+1)
     .           -III(i,j-1)-III(i-1,j-1))/4./dt
            j2 =-(III(i,j)  -III(i-1,j))/dh
            j3 = jac*lap(i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                        ,psi,dff=iR2,vol=.false.)

            !Find covariant EOM force components
            f_1 = (j2*b3(i,j,k)-j3*b2(i,j,k))/jac - dpdr
            f_2 = (j3*b1(i,j,k)-j1*b3(i,j,k))/jac - dpdth
            f_3 = (j1*b2(i,j,k)-j2*b1(i,j,k))/jac

            !Accumulate errors
            e1 = e1 + f_1**2*gmetric%grid(igrid)%dvol(i,j,k)
            e2 = e2 + f_2**2*gmetric%grid(igrid)%dvol(i,j,k)
            e3 = e3 + f_3**2*gmetric%grid(igrid)%dvol(i,j,k)

            !Find covariant Ohm's law components
            etah =res_rfx(rh)

            cov = matmul(gmetric%grid(1)%gsub(i,j,k,:,:),(/ j1,j2,j3/))

            f_1 =      (v2(i,j,k)*b3(i,j,k)                    )/jac
            f_2 =      (                   -v1(i,j,k)*b3(i,j,k))/jac
            f_3 =Eoeta+(v1(i,j,k)*b2(i,j,k)-v2(i,j,k)*b1(i,j,k))/jac

            f_1 =f_1-etah*cov(1)
            f_2 =f_2-etah*cov(2)
            f_3 =f_3-etah*cov(3)

            !Accumulate errors
            e4 = e4 + f_1**2*gmetric%grid(igrid)%dvol(i,j,k)
            e5 = e5 + f_2**2*gmetric%grid(igrid)%dvol(i,j,k)
            e6 = e6 + f_3**2*gmetric%grid(igrid)%dvol(i,j,k)
          enddo
        enddo

#if defined (petsc)
        lpsi0 = e1
        call MPI_Allreduce(lpsi0,e1,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)

        lpsi0 = e2
        call MPI_Allreduce(lpsi0,e2,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)

        lpsi0 = e3
        call MPI_Allreduce(lpsi0,e3,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)

        lpsi0 = e4
        call MPI_Allreduce(lpsi0,e4,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)

        lpsi0 = e5
        call MPI_Allreduce(lpsi0,e5,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)

        lpsi0 = e6
        call MPI_Allreduce(lpsi0,e6,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)
#endif

        gserr = e1 + e2 + e3
        gserr = sqrt(gserr/pi/xmax**2)

        e1 = sqrt(e1/pi/xmax**2)
        e2 = sqrt(e2/pi/xmax**2)
        e3 = sqrt(e3/pi/xmax**2)

        oherr = e4 + e5 + e6
        oherr = sqrt(oherr/pi/xmax**2)

        e4 = sqrt(e4/pi/xmax**2)
        e5 = sqrt(e5/pi/xmax**2)
        e6 = sqrt(e6/pi/xmax**2)

        if (my_rank == 0) then
          write (*,*) 'EOM force residual comp 1=',e1
          write (*,*) 'EOM force residual comp 2=',e2
          write (*,*) 'EOM force residual comp 3=',e3
          write (*,*) 'EOM force residual total =',gserr
          write (*,*)
          write (*,*) 'Ohms law residual comp 1=',e4
          write (*,*) 'Ohms law residual comp 2=',e5
          write (*,*) 'Ohms law residual comp 3=',e6
          write (*,*) 'Ohms law residual total =',oherr
        endif

c       Calculate GS equilibrium properties

        tor_flx = 0d0
        p_avg   = 0d0
        bz2_avg = 0d0
        bth_avg = 0d0
        I0      = 0d0

        k = 1
        do j=1,ny
          do i=1,nx
            iglobal = i + grid_params%ilo(igrid) - 1
cc            jglobal = j + grid_params%jlo(igrid) - 1

            call getMGmap(i,j,k,igrid,igrid,igrid,ig,jg,kg)

            rh = grid_params%xx(ig)
            dh = grid_params%dxh(ig)
            dt = grid_params%dyh(jg)

            !Q-profile
            qprof(i,j) = b3(i,j,k)/b2(i,j,k)

            !Toroidal flux
            tor_flx = tor_flx + b3(i,j,k)*dh*dt

            !Bz^2 average
            bz2_avg = bz2_avg+(b3(i,j,k))**2/rh*dh*dt

            !Pressure average
            p_avg = p_avg + pp(i,j,k)*rh*dh*dt

            !Toroidal current
            if (iglobal == nxd) then
              I0 = I0 + (psi(i+1,j,k)-psi(i,j,k))/dh
     .                  *2./(RR(i,j)+RR(i+1,j))*dt
            endif

          enddo
        enddo

#if defined (petsc)
        lpsi0 = tor_flx
        call MPI_Allreduce(lpsi0,tor_flx,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)

        lpsi0 = bz2_avg
        call MPI_Allreduce(lpsi0,bz2_avg,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)

        lpsi0 = p_avg
        call MPI_Allreduce(lpsi0,p_avg,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)

        lpsi0 = I0
        call MPI_Allreduce(lpsi0,I0,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_XY,mpierr)
#endif

        !Q-profile diagnostics (not working)
cc        q0 = mval(qprof)
cc#if defined (petsc)
cc        lpsi0 = q0
cc        call MPI_Allreduce(lpsi0,q0,1,MPI_DOUBLE_PRECISION
cc     .                    ,MPI_MIN,MPI_COMM_WORLD,mpierr)
cc#endif
cc        qa = maxval(qprof)
cc#if defined (petsc)
cc        lpsi0 = qa
cc        call MPI_Allreduce(lpsi0,qa,1,MPI_DOUBLE_PRECISION
cc     .                    ,MPI_MAX,MPI_COMM_WORLD,mpierr)
cc#endif

        bth_avg = I0/2/pi   !Average Btheta

        b_pol = 2*p00*pprof(0d0,0)/bth_avg**2  !Poloidal beta

        b_tor = 2*p_avg/bz2_avg !Toroidal beta

        bz_avg = 0.5*tor_flx/pi/xmax**2  !Average Bz

        if (my_rank == 0) then
          write (*,*) 
          write (*,999) ' Toroidal Ohmic equilibrium features:'
          write (*,999) '   Theta   =',bth_avg/bz_avg
          write (*,999) '   Bz(r=0) =',1d0
          write (*,999) '   I0      =',I0
          write (*,999) '   E0/eta  =',Eoeta/R0
          write (*,999) '   Tor_flx =',tor_flx
cc          write (*,999) '   q(psi=0)=',q0
cc          write (*,999) '   q(psi=1)=',qa
          write (*,999) '   beta_pol [2 p0/Bth(a)^2]=',b_pol
          write (*,999) '   beta_tor [2 <p>/<Bz^2>] =',b_tor
 999      format (a,f10.7)
        endif

c     Change sign of psi (to agree with A_3 sign convention)

        psi = -psi

        !Save some equilibrium properties (to calculate bootstrap current)
cc        eq_params(1) = 1d0/grid_params%params(3)       !Aspect ratio
cc        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
cc        eq_params(3) = q0h(nxd)/grid_params%params(3)  !q(r=a)

        !Dump equilibrium for eigen code
cc        if (my_rank == 0) then
cc          open(unit=123,file='ntm2d.txt',status='unknown')
cc          write (123,*) nxdp+1
cc          write (123,*) grid_params%xg(0:nxdp)
cc          write (123,*) bth0
cc          write (123,*) bzz0
cc          write (123,*) 0.5*p0
cc          write (123,*) grid_params%params(3),eq_params(2),eq_params(3)
cc          close(123)
cccc          stop
cc        endif

c diag ***** 
cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .    ,status='replace')
cccc      call contour(psi(:,:,1),nx+2,ny+2,0d0,1d0,0d0,1d0,0,110)
cc      call contour(errv            ,nx,ny,0d0,1d0,0d0,1d0,0,110)
cc      call contour(b1 (1:nx,1:ny)  ,nx,ny,0d0,1d0,0d0,1d0,1,110)
cc      call contour(b2 (1:nx,1:ny)  ,nx,ny,0d0,1d0,0d0,1d0,1,110)
cc      call contour(b3 (1:nx,1:ny)  ,nx,ny,0d0,1d0,0d0,1d0,1,110)
cc      call contour(v1 (1:nx,1:ny)  ,nx,ny,0d0,1d0,0d0,1d0,1,110)
cc      call contour(v2 (1:nx,1:ny)  ,nx,ny,0d0,1d0,0d0,1d0,1,110)
cc      close (110)
cc      stop
c diag *****

        deallocate(psi_rhs,psiv,errv)

      end subroutine Ohm_GS_equil_2d

c     qprof
c     ###############################################################
      function qprof(r)

c     ---------------------------------------------------------------
c     Finds cylindrical components of magnetic field B_theta, B_z
c     from q profile.
c     ---------------------------------------------------------------

        use equilibrium

        use app_iosetup

        implicit none

c     Call variables

        real(8)    :: qprof,r

c     Local variables

        integer :: p
        real(8) :: q0,qa

c     Begin program

        select case(equil)
        case('ntm2d')

          q0 = 1.02
          qa = 4.
          p  = 1

          qprof = q0*(1+r**(2*p)*( (qa/q0)**p - 1 ))**(1./p)

        case('rfp83')

          qprof = 0.6125*(1 - 1.8748*r**2 + 0.8323*r**4)

        case default

          call pstop('qprof','Undefined q profile')

        end select

c     End program

      end function qprof

c     pprof
c     ###############################################################
      function pprof(s,der)

c     ---------------------------------------------------------------
c     Specifies pressure profile (der = 0) and its derivative (der=1)
c     as a function of poloidal flux.
c     ---------------------------------------------------------------

        use equilibrium

        implicit none

c     Call variables

        integer    :: der
        real(8)    :: pprof,s

c     Local variables

c     Begin program

        select case(equil)
        case('ntm2d','ppnsl')

          if (der == 0) then
            pprof = (1-s**2)**4! + 1d-3
          else
            pprof =-8*s*(1-s**2)**3
          endif

        case('ohtor')

          if (der == 0) then
            pprof = (1-s**2)**2! + 1d-3
          else
            pprof =-4*s*(1-s**2)
          endif

        case('p3nsl')

          if (der == 0) then
            pprof = (1-s**2)*(1.4-s**2)! + 1d-3
          else
            pprof =-4*s*(1-s**2)
          endif

        case default

          pprof = 0d0

        end select

c     End program

      end function pprof

c fillVectorPotential
c####################################################################
      subroutine fillVectorPotential(igr,a)

c--------------------------------------------------------------------
c     Defines COVARIANT vector potential for initialization of 
c     equilibrium.
c--------------------------------------------------------------------

      use parameters

      use grid

      use equilibrium

      use imposeBCinterface

      implicit none

c Call variables

      integer :: igr
      real(8) :: a(ilom:ihip,jlom:jhip,klom:khip,3)

c Local variables

      integer :: ig,jg,kg,i,j,k,bcnd(6,3),nx,ny,nz
      real(8) :: x1,y1,z1,eps,car(3),hlx
      real(8) :: a_cnv(ilom:ihip,jlom:jhip,klom:khip,3)
      logical :: cartsn

c Begin program

      nx = grid_params%nxv(igr)
      ny = grid_params%nyv(igr)
      nz = grid_params%nzv(igr)

      select case (trim(equil))

      case ('khcar')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igr,igr,igr,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = x1
              a(i,j,k,3) = 0d0
            enddo
          enddo
        enddo

      case ('gem','kai1')

        if (bcond(1) /= FSYM .and. bcond(2) == FSYM) then
          hlx = xmax
        elseif (bcond(1) == FSYM .and. bcond(2) /= FSYM) then
          hlx = xmin
        else
          hlx = 0.5d0*(xmax+xmin)
        endif

cc        write (*,*) 'DIAG -- fillVectorPot'

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igr,igr,igr,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = -dlambda*dlog(cosh((x1-hlx)/dlambda)) 
cc              a(i,j,k,3) = 0d0  !<<<DIAG

            enddo
          enddo
        enddo

      case ('tmcar','tmsin')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igr,igr,igr,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 2*dlambda
     .               *atan(exp((x1-0.5d0*(xmax+xmin))/dlambda))
              a(i,j,k,3) = -dlambda
     .               *log(cosh((x1-0.5d0*(xmax+xmin))/dlambda)) 
            enddo
          enddo
        enddo

      case ('dtm')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igr,igr,igr,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = -dlambda
     .                 *( log(cosh((x1-xmin)/dlambda))
     .                   -log(cosh((x1-xmax)/dlambda))
     .                   -x1/dlambda)
            enddo
          enddo
        enddo

      case ('ic')

        eps = grid_params%params(2)
        dlambda = (ymax-ymin)/2./pi/2.

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igr,igr,igr,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = -dlambda*dlog(    cosh(x1/dlambda)
     $                                  +eps*cos (y1/dlambda)) 
            enddo
          enddo
        enddo

      case ('ic-hk')

        eps = eq_params(2)
        dlambda = (ymax-ymin)/pi !(Ly/4/pi)

cc        if (bcond(1) /= FSYM .and. bcond(2) == FSYM) then
cc          hlx = xmax
cc        elseif (bcond(1) == FSYM .and. bcond(2) /= FSYM) then
cc          hlx = xmin
cc        else
cc          hlx = 0.5d0*(xmax+xmin)
cc        endif

        if (bcond(3) == FSYM .or. bcond(4) == FSYM) then
          dlambda = (ymax-ymin)/2./pi !(2*Ly/4/pi)
        else
          dlambda = (ymax-ymin)/4./pi !(Ly/4/pi)
        endif

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igr,igr,igr,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = -dlambda*dlog(    cosh(x1/dlambda)
     $                                  +eps*cos (y1/dlambda)) 
            enddo
          enddo
        enddo

      case ('tmcyz','3dkh')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(i,j,k,igr,igr,igr,ig,jg,kg
     .                                    ,x1,y1,z1)
              a(i,j,k,1) = dlambda
     .             *dlog(cosh((y1-0.5d0*(ymax+ymin))/dlambda)) 
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = 0d0
            enddo
          enddo
        enddo

      end select

c Transform to Covariant

      call XformVector_mesh2(igr,a,'car','cov')

c End program

      end subroutine fillVectorPotential

c findVectorPotential
c####################################################################
      subroutine findVectorPotential(nxg,nyg,nzg,v,b,a_cov)

c--------------------------------------------------------------------
c     Finds COVARIANT vector potential from CONTRAVARIANT magnetic
c     field components with gauge:
c                grad(phi + eta div(a)) = 0
c
c     On input:
c       * v -> cnv velocity
c       * b -> cnv magnetic field
c     On output:
c       * a -> cov vector potential
c--------------------------------------------------------------------

      use parameters

      use equilibrium

      use operators

      use local_BCS_variables, ONLY: order_bc

      use imposeBCinterface

      use timeStepping

      use B_tools

      implicit none

c Call variables

      integer    :: nxg,nyg,nzg
      real(8)    :: a_cov(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,v(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,b(0:nxg+1,0:nyg+1,0:nzg+1,3)

c Local variables

      integer    :: i,j,k,igrid,ii,bcnd(6,3),nx,ny,nz
     .             ,igm,igp,jgm,jgp,kgm,kgp
      real(8)    :: cm,c0,cp,dxx,dyy,dvol,cov(3)
     .             ,cai_rhs(ntotdp),caiv(ntotdp)
     .             ,phi(ilom:ihip,jlom:jhip,klom:khip)
     .             ,cai(ilom:ihip,jlom:jhip,klom:khip)
     .             ,a_cnv(ilom:ihip,jlom:jhip,klom:khip,3)
     .             ,a(0:nxg+1,0:nyg+1,0:nzg+1,3)
     .             ,vzeros(ilom:ihip,jlom:jhip,klom:khip,3)

      logical    :: spoint

c Begin program

      igrid = 1

      nx = grid_params%nxv(igrid)
      ny = grid_params%nyv(igrid)
      nz = grid_params%nzv(igrid)

      vzeros = 0d0

c Find vector potential with gauge Ax = 0 in GLOBAL domain

      call curl_inv(nxg,nyg,nzg,igrid,igrid,igrid,b,a,global=.true.)

c Restrict solution to current domain

      if (nxg == nxd) then
        igm = ilomg ; igp = ihipg
      else
        igm = ilom  ; igp = ihip
      endif

      if (nyg == nyd) then
        jgm = jlomg ; jgp = jhipg
      else
        jgm = jlom  ; jgp = jhip
      endif

      if (nzg == nzd) then
        kgm = klomg ; kgp = khipg
      else
        kgm = klom  ; kgp = khip
      endif

cc      a_cov(ilom:ihip,jlom:jhip,klom:khip,:)
cc     .     = a(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)
      a_cov(ilom:ihip,jlom:jhip,klom:khip,:)
     .     = a(igm:igp,jgm:jgp,kgm:kgp,:)

c Set topological BCs

      bcnd(:,1) = bcond
      bcnd(:,2) = bcond
      bcnd(:,3) = bcond
cc      where (bcnd == DEF) bcnd = -IFC

      a_cnv = XformVector(igrid,a_cov,.true.)
      call setBC(IAX,3,nx,ny,nz,a_cnv,a_cov,vzeros,bcnd
     .          ,igrid,igrid,igrid,is_cnv=.false.
     .          ,iorder=order_bc)

c Find equilibrium (d_t A = 0) scalar potential (vxb electric field)

cc      call findScalarPotential(v
cc     .                        ,b(ilomg:ihipg,jlomg:jhipg,klomg:khipg,:)
cc     .                        ,phi)

ccc Impose **topological constraints** and find a_cnv
cc
cc      bcnd(:,1) = bcond
cc      bcnd(:,2) = bcond
cc      bcnd(:,3) = bcond
cccc      where (bcnd == DEF) bcnd = EXT
cc
cc      a_cnv = XformVector(igrid,a_cov,.true.)
cc      call setBC(IAX,3,nx,ny,nz,a_cnv,a_cov,vzeros,bcnd
cc     .          ,igrid,igrid,igrid,is_cnv=.false.,iorder=2)
cc
ccc Transform A to gauge div(A)=0
cc
cc      !Find gauge transformation potential cai
cc      do k=klo,khi
cc        do j=jlo,jhi
cc          do i=ilo,ihi
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc            if (vol_wgt) then
cc              dvol = gmetric%grid(igrid)%dvol(i,j,k)
cc            else
cc              dvol = 1d0
cc            endif
cc
cc            cai_rhs(ii)=-dvol*(div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                            ,a_cnv(:,:,:,1)
cc     .                            ,a_cnv(:,:,:,2)
cc     .                            ,a_cnv(:,:,:,3)))
cc          enddo
cc        enddo
cc      enddo
cc
cc      !Transform A (if rhs /= 0d0)
cc      if (sqrt(sum(cai_rhs*cai_rhs)) > 0d0) then
cc
cc        !Solve Poisson problem to find Gauge potential cai
cc#if defined(petsc)
cc        call pstop('findVectorPotential'
cc     .          ,'Poisson solver not implemented in parallel')
cc#else
cc        write (*,*) 'Solving for Gauge potential'
cc
cc        !Define BCs
cc        bcnd(:,1) = bcond
cc        where (bcnd == DEF) bcnd = DIR
cc
cc        !Solve Poisson eq
cc        call poisson(ntotdp,cai_rhs,caiv,30,1d-6,bcnd(:,1),ilevel,0)
cc
cc        !Map to array (w BCs)
cc        call mapMGVectorToArray(0,1,caiv,nx,ny,nz,cai,igrid,.false.)
cc
cc        call setBC(IRHO,nx,ny,nz,cai,zeros,bcnd(:,1),igrid,igrid,igrid
cc     .            ,iorder=2)
cc#endif
cc
cc        !Gauge transformation A' = A + grad(cai)
cc        do k=klom,khip
cc          do j=jlom,jhip
cc            do i=ilom,ihip
cc              call grad(i,j,k,nx,ny,nz,igrid,igrid,igrid,cai
cc     .                 ,cov(1),cov(2),cov(3))
cc              a_cov(i,j,k,:) = a_cov(i,j,k,:) + cov
cc            enddo
cc          enddo
cc        enddo
cc
cc      endif

c End program

      end subroutine findVectorPotential

ccc findScalarPotential
ccc####################################################################
cc      subroutine findScalarPotential(v,b,phi)
cc
ccc--------------------------------------------------------------------
ccc     Finds scalar equilibrium potential from equilibrium velocity and 
ccc     magnetic fields, obeying the relation:
ccc     
ccc         lap(phi_0) = div(v0 x B0)
ccc
ccc     This automatically takes into account ExB flows at the boundary.
ccc--------------------------------------------------------------------
cc
cc      use parameters
cc
cc      use grid
cc
cc      use equilibrium
cc
cc      use PC_var
cc
cc      implicit none
cc
ccc Call variables
cc
cc      real(8)    :: v(ilom:ihip,jlom:jhip,klom:khip,3)
cc     .             ,b(ilom:ihip,jlom:jhip,klom:khip,3)
cc     .             ,phi(ilom:ihip,jlom:jhip,klom:khip)
cc
ccc Local variables
cc
cc      integer    :: igrid,ii,bcnd(6,1)
cc      real(8)    :: phi_rhs(ntotdp),phiv(ntotdp),cov(3)
cc     .             ,dum(ilom:ihip,jlom:jhip,klom:khip,3),dvol
cc
ccc Begin program
cc
cc      igrid = 1
cc
cc      nx = grid_params%nxv(igrid)
cc      ny = grid_params%nyv(igrid)
cc      nz = grid_params%nzv(igrid)
cc
ccc Find rhs
cc
cc      do k=klom,khip
cc        do j=jlom,jhip
cc          do i=ilom,ihip
cc            jac = gmetric%grid(igrid)%jac(i,j,k)
cc
cc            cov(1) = v(i,j,k,2)*b(i,j,k,3)
cc     .             - v(i,j,k,3)*b(i,j,k,2)
cc
cc            cov(2) = v(i,j,k,3)*b(i,j,k,1)
cc     .             - v(i,j,k,1)*b(i,j,k,3)
cc
cc            cov(3) = v(i,j,k,1)*b(i,j,k,2)
cc     .             - v(i,j,k,2)*b(i,j,k,1)
cc
cc            cov = cov/jac
cc
cc            call transformFromCurvToCurv(i,j,k,igrid,igrid,igrid
cc     .                                  ,cov(1),cov(2),cov(3)
cc     .                                  ,dum(i,j,k,1)
cc     .                                  ,dum(i,j,k,2)
cc     .                                  ,dum(i,j,k,3),.true.)
cc
cc          enddo
cc        enddo
cc      enddo
cc
cc      do k=klo,khi
cc        do j=jlo,jhi
cc          do i=ilo,ihi
cc            ii  = i + nx*(j-1) + nx*ny*(k-1)
cc
cc            if (vol_wgt) then
cc              dvol = volume(i,j,k,igrid,igrid,igrid)
cc            else
cc              dvol = 1d0
cc            endif
cc
cc            phi_rhs(ii) = dvol*div(i,j,k,nx,ny,nz,igrid,igrid,igrid
cc     .                        ,dum(:,:,:,1),dum(:,:,:,2),dum(:,:,:,3))
cc          enddo
cc        enddo
cc      enddo
cc
ccc Solve for potential
cc
cc      if (sqrt(sum(phi_rhs*phi_rhs))>0d0) then
cc
cc#if defined(petsc)
cc
cc        call pstop('findScalarPotential'
cc     .          ,'Poisson solver not implemented in parallel')
cc
cc#else
cc
cc        !Define BCs
cc        bcnd(:,1) = bcond
cc        where (bcnd == DEF) bcnd = DIR
cc
cc        !Solve Poisson eq
cc        call poisson(ntotdp,phi_rhs,phiv,30,1d-6,bcnd,ilevel,0)
cc
cc        !Map to array (w/o BCs)
cc        call mapMGVectorToArray(0,1,phiv,nx,ny,nz,phi,igrid,.false.)
cc#endif
cc
cc      else
cc
cc        phi = 0d0
cc
cc      endif
cc
ccc End program
cc
cc      end subroutine findScalarPotential
