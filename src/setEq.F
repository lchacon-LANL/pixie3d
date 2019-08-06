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
c       * ieq=IPRS: total pressure is found as p=n*(Ti+Te)=a_p*n*Te.
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

      use B_tools, ONLY: B_div_clean,findVectorPotential

      implicit none

c Call variables

      integer :: iigx

      type(var_array),pointer :: varray

c Local variables

      integer :: i,j,k,ig,jg,kg,ieq,nmax,iglobal,jglobal
     .          ,ierr,sh,it,icomp,bcnd(6,3),brank
     .          ,igx,igy,igz,nx,ny,nz,igr,nxx,nyy,nzz,mmm

      real(8) :: ldaflow,x1,y1,z1,cx,cy,cz,x00

      logical :: covariant,to_cartsn,cartsn,a_bc,find_a,div_clean

      real(8) :: r,jac1,bnorm,aaa,bbb,ccc,qq,qqp,q0,RR,ff,nn,b_0
     .          ,aspect_ratio,mm,kk,Iz,bz0,btheta,rrho,r1
     .          ,r0(1),a1,b1,dummy1,dummy2,dummy3,eta_in,err,bz_avg,dr
     .          ,bzh,rh,e0_over_eta,m0_over_nu,theta_ppnch,p0,prsw,eps
     .          ,phi,rho_norm,Rtor,sigma,gsub(3,3),bb0(3)
     .          ,hlx,y_char,besselIp,besselI,br,bth,bz,kn,pv,pw
     .          ,Impv,Kmpv,Kmp,Imp,Km,Im,vr,vth,vz
     .          ,Er,Eth,Ez,btho,bzo,bthi,bzi

      real(8) :: dh,bth0(0:nxdp)
      real(8) :: bzz0(0:nxdp),vr0(0:nxdp),rho0(0:nxdp),vz0(0:nxdp)
      real(8) :: eth0(0:nxdp),ez0(0:nxdp),jth0(0:nxdp)
      real(8) :: jz0 (0:nxdp),bth0_old(0:nxdp),p00(0:nxdp)
      real(8) :: bzz0_old(0:nxdp),dummy(0:nxdp),bmag2(0:nxdp)
      real(8) :: q00(0:nxdp),dum (ilom:ihip,2)
      real(8) :: ldum(ilom:ihip,2)
      real(8) :: r2  (ilom:ihip,jlom:jhip)
      real(8) :: br02d (0:nxdp,0:nydp),bth02d(0:nxdp,0:nydp)
      real(8) :: bzz02d(0:nxdp,0:nydp)
      real(8) :: vth02d(0:nxdp,0:nydp),vr02d(0:nxdp,0:nydp)
      real(8) :: psi0(0:nxdp),rg(0:nxdp)

      real(8) :: pp(ilom:ihip,jlom:jhip)

      real(8) :: aa(ilom:ihip,jlom:jhip,klom:khip,3)
     .          ,bb(ilom:ihip,jlom:jhip,klom:khip,3)

      real(8) :: var(ilom:ihip,jlom:jhip,klom:khip,neqd)

      character(20),pointer,dimension(:) :: label

      integer :: bcs(6,neqd)

c Functions

      logical :: fillVectorPotential
      external   fillVectorPotential

      real(8) :: qprof,pprof
      external   qprof,pprof

c Begin program

      find_a = .true.

      ldaflow = dlambda/rshear

      E0 = 0d0 ; B0 = 0d0 ; bth0 = 0d0

      label => varray%array_var%descr

      do ieq=1,neqd
        var(:,:,:,ieq) = varray%array_var(ieq)%array
      enddo

c Initialize required local grid information

      igx = iigx
      igy = iigx
      igz = iigx

      igr = iigx
      
      nx = gv%gparams%nxv(igx)
      ny = gv%gparams%nyv(igy)
      nz = gv%gparams%nzv(igz)

c Label variables
     
      label(IRHO) = 'J.Rho'
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
      label(IPRS) = 'J.(pi+pe)'

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

        var(:,:,:,IPRS) = 1d0

#if !defined(vec_pot)
        bb(:,:,:,1)  = 0d0
        bb(:,:,:,2)  = 0d0
        bb(:,:,:,3)  = 1d0
#else
        B0(3) = 1d0

        aa = 0d0 ; find_a = .false.
#endif

      case ('3wv') !3-wave MHD test

        gamma = 1.6667

        phi   = gv%gparams%params(1)   !Angle of B with x-axis, in degrees
        beta  = gv%gparams%params(2)   !plasma beta

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

        var(:,:,:,IPRS) = 0.5*beta/a_p  !Electron temperature

        B0(1) = cos(phi)
        B0(2) = sin(phi)
        B0(3) = 0d0

#if !defined(vec_pot)
        bb(:,:,:,1)  = B0(1)
        bb(:,:,:,2)  = B0(2)
        bb(:,:,:,3)  = B0(3)
#else
        aa = 0d0 ; find_a = .false.
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
        var(:,:,:,IPRS) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

              bb(i,j,k,1)  = 0d0
              bb(i,j,k,2)  = 0d0
              bb(i,j,k,3)  = jac1*1d0
            enddo
          enddo
        enddo

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

        var(:,:,:,IPRS) = 1d0

#if !defined(vec_pot)
        bb(:,:,:,1)  = 1d0
        bb(:,:,:,2)  = 0d0
        bb(:,:,:,3)  = 0d0
#else
        B0(1) = 1d0

        aa = 0d0 ; find_a = .false.
#endif

      case ('kaw') !KAW wave

        beta = 1d4

        di = (xmax-xmin)/2/pi/nh1*sqrt(beta/(beta-2))  !di k_par = 1

        call write_real_field_to_inputf(di,"di")

c     Check input

        if (coords /= 'car') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (test
     .     .and.(   (adiabatic)
     .          .or.(a_p /= 1d0)
     .          .or.(.not.solve_rho)
     .          .or.(gamma /= 1d0)
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

        var(:,:,:,IPRS) = 0.5*beta/a_p

#if !defined(vec_pot)
        bb(:,:,:,1)  = 1d0
        bb(:,:,:,2)  = 0d0
        bb(:,:,:,3)  = 0d0
#else
        B0(1) = 1d0

        aa = 0d0 ; find_a = .false.
#endif

c     --------------------------------------------
c     ------------  ADVECTION TESTS --------------
c     --------------------------------------------

      case ('adv1') !Advection test (Germaschewski)

        gamma = 2d0

        eps   = gv%gparams%params(1)
        sigma = gv%gparams%params(2)
        x00   = gv%gparams%params(3)

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
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .             ,ig,jg,kg,x1,y1,z1)

c$$$              var(i,j,k,IRHO) = 1d0 + eps*exp(-0.5*(x1-0.75)**2/sigma**2)
c$$$              var(i,j,k,IPRS) = 1d0 - eps*exp(-0.5*(x1-0.75)**2/sigma**2)
              var(i,j,k,IRHO) = 1d0 +eps*exp(-0.5*(x1-x00)**2/sigma**2)
              var(i,j,k,IPRS) = 1d0 -eps*exp(-0.5*(x1-x00)**2/sigma**2)

            enddo
          enddo
        enddo

#if !defined(vec_pot)
        bb = 0d0
#else
        aa = 0d0 ; find_a = .false.
#endif

      case ('adv2') !Advection test (Jasak et al, IJNMF 31, 1999; Lenard, CMAME 88, 1991)

        gamma = 2d0

        eps   = gv%gparams%params(1)
        sigma = gv%gparams%params(2)

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
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .             ,ig,jg,kg,x1,y1,z1)

              y_char = 1./6.+0.5/0.866*x1
              var(i,j,k,IRHO) = 1d0+eps*0.5*(1d0+sign(1d0,y1-y_char))
              var(i,j,k,IPRS) = 1d0-eps*0.5*(1d0+sign(1d0,y1-y_char))

            enddo
          enddo
        enddo

#if !defined(vec_pot)
        bb = 0d0
#else
        aa = 0d0 ; find_a = .false.
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
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .             ,ig,jg,kg,x1,y1,z1)

              var(i,j,k,IRHO) = 1d0-0.5*(1d0+sign(1d0,x1-0.5))*0.875
              var(i,j,k,IPRS) = 1d0-0.5*(1d0+sign(1d0,x1-0.5))*0.9  !Pressure

            enddo
          enddo
        enddo

        var(:,:,:,IPRS) = var(:,:,:,IPRS)/var(:,:,:,IRHO)/a_p !Temperature

#if !defined(vec_pot)
        bb = 0d0
#else
        aa = 0d0 ; find_a = .false.
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

        var(:,:,:,IPRS) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .                                    ,ig,jg,kg,x1,y1,z1)
              var(i,j,k,IVY)=vperflow
     .                      *tanh((x1-0.5d0*(xmax+xmin))/dlambda)
            enddo
          enddo
        enddo

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVZ)  = 0d0

#if !defined(vec_pot)
        bb(:,:,:,1)  = 0d0
        bb(:,:,:,2)  = 0d0
        bb(:,:,:,3)  = 1d0
#else
        find_a = fillVectorPotential(igx,aa)
#endif

      case ('3dkh')

c     Define vector potential (in curvilinear coordinates) for initialization

        find_a = fillVectorPotential(igx,aa)

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

        var(:,:,:,IPRS) = 0.5*beta/a_p  !Electron temperature

#if !defined(vec_pot)
        bb = curl(gv%gparams,igx,aa)
#endif

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .                                    ,ig,jg,kg,x1,y1,z1)
#if !defined(vec_pot)
              bb0 = bb(i,j,k,:)
#else
              bb0 = curl(gv%gparams,i,j,k,nx,ny,nz,igx,igy,igz,aa)
#endif
              var(i,j,k,IRHO) = 2 - bb0(3)**2

              var(i,j,k,IVX) =vperflow*tanh((y1-.5*(ymax+ymin))/ldaflow)
     .                                *sin (2*pi*z1/(zmax-zmin))
              var(i,j,k,IVY) = 0d0
              var(i,j,k,IVZ) = 0d0

            enddo
          enddo
        enddo

      case ('khihk')  !Homa Karimabadi's KHI

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
        var(:,:,:,IPRS) = 0.5*beta/a_p

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .             ,ig,jg,kg,x1,y1,z1)
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
        bb(:,:,:,1)  = B0(1)
        bb(:,:,:,2)  = B0(2)
        bb(:,:,:,3)  = B0(3)
#else
        aa = 0d0 ; find_a = .false.
#endif

      case ('khidp')

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

c     Double periodic Kelvin-Helmholtz with constant magnetic field and pressure
c     in cartesian coordinates
        
        var(:,:,:,IRHO) = 1d0

        beta = 0.3d0
        var(:,:,:,IPRS) = 0.5*beta/a_p

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .                                    ,ig,jg,kg,x1,y1,z1)

              var(i,j,k,IVY)=vperflow
     .                      *(tanh((x1-(xmin+(xmax-xmin)/4d0))/dlambda)
     .                   -tanh((x1-(xmin+3d0*(xmax-xmin)/4d0))/dlambda)
     .                    -1d0)
            enddo
          enddo
        enddo

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        B0(1) = 0d0
        B0(2) = 0.04998d0
        B0(3) = 0.998d0

#if !defined(vec_pot)
        bb(:,:,:,1)  = B0(1)
        bb(:,:,:,2)  = B0(2)
        bb(:,:,:,3)  = B0(3)
#else
        aa = 0d0 ; find_a = .false.
#endif

c     ----------------------------------------------
c     ---------------  TEARING MODE ----------------
c     ----------------------------------------------

      case ('tmcar')

c     Define vector potential (in curvilinear coordinates) for initialization

        find_a = fillVectorPotential(igx,aa)

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

        var(:,:,:,IPRS) = 1d0

        bz0 = 1d0

#if !defined(vec_pot)
        bb = curl(gv%gparams,igx,aa)
        bb(:,:,:,3)=sqrt(bz0**2 - bb(:,:,:,2)**2)
#endif

      case ('tmsin')

c     Define vector potential (in curvilinear coordinates) for initialization

        find_a = fillVectorPotential(igx,aa)

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

        var(:,:,:,IPRS) = 1d0

        bz0 = 1d0

#if !defined(vec_pot)
        bb = curl(gv%gparams,igx,aa)

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

             !X-Y equilibrium
             gsub = gv%gparams%gmetric%grid(igx)%gsub(i,j,k,:,:)
             jac1 = gv%gparams%gmetric%grid(igx)%jac (i,j,k)

             bnorm =vectorNorm(gv%gparams,i,j,k,igx,bb(i,j,k,:),.false.)

             ccc = jac1*(bz0**2 - bnorm)
             bbb = gsub(3,2)*bb(i,j,k,2) + gsub(3,1)*bb(i,j,k,1)
             aaa = gsub(3,3)

             bb(i,j,k,3)  = (-bbb+sqrt(bbb**2+4*aaa*ccc))/2./aaa

             if (bbb**2+4*aaa*ccc < 0d0) then
               write (*,*) bb(i,j,k,3),bnorm,aaa,bbb**2+4*aaa*ccc
             endif
            enddo
          enddo
        enddo
#endif

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

        var(:,:,:,IPRS) = 0.5*beta/a_p !Electron temperature required for force balance

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
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .                                    ,ig,jg,kg,x1,y1,z1)

             var(i,j,k,IRHO)=1d0/cosh((x1-hlx)/dlambda)**2 + 0.2
            enddo
          enddo
        enddo

        find_a = fillVectorPotential(igx,aa)

#if !defined(vec_pot)
        bb = curl(gv%gparams,igx,aa)
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

        var(:,:,:,IPRS) = 1d0 !Electron pressure

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .             ,ig,jg,kg,x1,y1,z1)

             var(i,j,k,IRHO)=0.5/cosh((x1-0.5d0*(xmax+xmin))/dlambda)**2
     .                      +1d0
            enddo
          enddo
        enddo

        find_a = fillVectorPotential(igx,aa)

#if !defined(vec_pot)
        bb = curl(gv%gparams,igx,aa)
#endif

      case ('dtm')

c     Define vector potential (in curvilinear coordinates) for initialization

        find_a = fillVectorPotential(igx,aa)

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

        var(:,:,:,IPRS) = 0.5d0  !Electron temperature

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .             ,ig,jg,kg,x1,y1,z1)

             var(i,j,k,IRHO)=0.5/cosh((x1-xmin)/dlambda)**2
     .                      +0.5/cosh((x1-xmax)/dlambda)**2
     .                      +1d0
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        bb = curl(gv%gparams,igx,aa)
#endif

c     --------------------------------------------------------------------
c     ----------  FORCE FREE HARRIS WITH TANH PRESSURE PROFILE -----------
c     --------------------------------------------------------------------

       case ('diamh')

c     Check coordinates                                                                                                                                              
         if (coords /= 'car' .and. coords /= 'lpk') then
           messg = 'Wrong coordinates for equilibrium '//equil
           call pstop('setEquilibrium',messg)
         endif

         var(:,:,:,IVX)  = 0d0
         var(:,:,:,IVY)  = 0d0
         var(:,:,:,IVZ)  = 0d0

         beta=1.0d0 ! For force balance in Hesse equilibrium
         var(:,:,:,IPRS) = 1.0*beta/a_p !Electron temperature required for force balance
 
         if (bcond(1) /= FSYM .and. bcond(2) == FSYM) then
           hlx = xmax
         elseif (bcond(1) == FSYM .and. bcond(2) /= FSYM) then
           hlx = xmin
         else
           hlx = 0.5d0*(xmax+xmin)
         endif

         var(:,:,:,IPRS) = beta/a_p !Electron temperature required for force balance

         do k = klom,khip
           do j = jlom,jhip
             do i = ilom,ihip
               call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .                                     ,ig,jg,kg,x1,y1,z1)

               var(i,j,k,IRHO)=1.0d0 + 0.5d0*tanh((x1-hlx)/rshear)
                 
             enddo
           enddo
         enddo

         if (my_rank == 0) write(*,*) 'About to fill B'

         do k = klom,khip
           do j = jlom,jhip
             do i = ilom,ihip
               call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .              ,ig,jg,kg,x1,y1,z1)

               bb(i,j,k,1)  = 0d0                                                                       
               bb(i,j,k,2)  = tanh((x1-hlx)/dlambda)
               bb(i,j,k,3)  = sqrt(4.0d0-(tanh((x1-hlx)/dlambda))**2
     .                      - tanh((x1-hlx)/rshear))
             enddo
           enddo
         enddo

         if (my_rank == 0) write(*,*) 'B filled'

c     -------------------------------------------------
c     ---------------  ISLAND COALESCENCE -------------
c     -------------------------------------------------

      case ('ic')

        p0  = eq_params(1)
        eps = eq_params(2)

        if (bcond(3) == FSYM .or. bcond(4) == FSYM) then
          dlambda = (ymax-ymin)/2./pi !(2*Ly/4/pi)
        else
          dlambda = (ymax-ymin)/4./pi !(Ly/4/pi)
        endif

c     Define vector potential (in curvilinear coordinates) for initialization

        find_a = fillVectorPotential(igx,aa)

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
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .             ,ig,jg,kg,x1,y1,z1)

             var(i,j,k,IPRS)=(p0 + 0.5*(1-eps**2)
     $                       /(    cosh(x1/dlambda)
     $                        +eps*cos (y1/dlambda))**2)
     $                       /a_p/var(i,j,k,IRHO)
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        bb = curl(gv%gparams,igx,aa)
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

        var(:,:,:,IPRS) = 0.5*beta/a_p

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz
     .             ,ig,jg,kg,x1,y1,z1)

             var(i,j,k,IRHO)=(0.2
     .                      +(1-eps**2)/(     cosh(x1/dlambda)
     .                                   +eps*cos (y1/dlambda))**2 )
            enddo
          enddo
        enddo

        find_a = fillVectorPotential(igx,aa)

#if !defined(vec_pot)
        bb = curl(gv%gparams,igx,aa)
#endif

c     ------------------------------------------------
c     ---------------  3D CYLINDRICAL ----------------
c     ------------------------------------------------

      case ('3dkai')

cc        mm = gv%gparams%params(1)
cc        kk = gv%gparams%params(2)
        RR = gv%gparams%params(3)

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

        Iz  = 2.22
        b_0 = 0.5*Iz

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

              call getCurvilinearCoordinates(gv%gparams,i,j,k,igx,igy
     .             ,igz,ig,jg,kg,x1,y1,z1)

              btheta = b_0*x1/(1+(x1/dlambda)**2)
              bz0    = RR

              bb(i,j,k,1)  = 0d0
              bb(i,j,k,2)  = btheta
              bb(i,j,k,3)  = x1*bz0

              var(i,j,k,IPRS) = 1d-5
     .                     + 1/16.*(Iz*dlambda/(1+(x1/dlambda)**2))**2

            enddo
          enddo
        enddo

cc      case ('3dpar')  !EQ for parallel transport test
cc
cccc        mm = gv%gparams%params(1)
cccc        kk = gv%gparams%params(2)
cc        RR = gv%gparams%params(3)
cc        q0 = gv%gparams%params(4)
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
cc        var(:,:,:,IPRS) = 1d-5
cc
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
cc              b_0  = q0*bz0*dlambda/RR
cc              btheta = b_0*(x1/dlambda)/(1+(x1/dlambda)**2)
cc
cc              bb(i,j,k,1) = 0d0
cc              bb(i,j,k,2) = btheta
cc              bb(i,j,k,3) = x1*bz0
cc
cc            enddo
cc          enddo
cc        enddo

      case ('3drfp')

cc        mm = gv%gparams%params(1)
cc        kk = gv%gparams%params(2)
        RR = gv%gparams%params(3)
        a1 = gv%gparams%params(4)

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

        var(:,:,:,IPRS) = 1d-5

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call getCurvilinearCoordinates(gv%gparams,i,j,k,igx,igy
     .                                      ,igz,ig,jg,kg,x1,y1,z1)

              b_0    = (dlambda**2+a1)
     .                 /sqrt((dlambda**2+a1)**2-dlambda**4)
              btheta = b_0*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-b_0**2*(1-1./(1+(x1/dlambda)**2)**2))

              bb(i,j,k,1)  = 0d0
              bb(i,j,k,2)  = btheta
              bb(i,j,k,3)  = x1*bz0
            enddo
          enddo
        enddo

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

        mm = gv%gparams%params(1)  !Poloidal perturbation
        kk = gv%gparams%params(2)  !Axial perturbation (-n)
        RR = gv%gparams%params(3)  !Toroidal radius
        a1 = gv%gparams%params(4)  !a^2 in screw pinch eq.

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        !Build equilibrium

        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,IPRS) = 1d-5

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call find_polar_coords(gv%gparams,igx,i,j,k,x1,y1,z1
     .                              ,coords)

              jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

              b_0    = (dlambda**2+a1)
     .                 /sqrt((dlambda**2+a1)**2-dlambda**4)
              btheta = b_0*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-b_0**2*(1-1./(1+(x1/dlambda)**2)**2))

              bb(i,j,k,1)  = 0d0
              bb(i,j,k,2)  = jac1*(btheta/x1 + kk/mm*bz0)
              bb(i,j,k,3)  = jac1*bz0
            enddo
          enddo
        enddo

      case ('bpnch')  !BENNET SCREW PINCH BALANCED BY DENSITY

c     Check coordinates

        if (coords /= 'hel' .and. coords /= 'cyl' ) then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     RFP equilibria (analytical)

        mm = gv%gparams%params(1) ! Azimuthal mode number.
        kk = gv%gparams%params(2) ! -n where n is the axial mode number.
        RR = gv%gparams%params(3) ! Major radius R=Lz/2pi.
        a1 = gv%gparams%params(4) ! Radius^2 of q=1 surface.

        if (coords == 'hel') then
           nh2 = mm             ! To set the right perturbation wavelength
           nh1 = nh2            ! For perturbations to be of the right order
        end if

        bz0 = RR/(1d0 + a1)/dlambda

c$$$        if (my_rank == 0) then
c$$$          write (*,*) 'm,k,R0,rc,lam,bz0:',mm,kk,RR,a1,dlambda,bz0
c$$$        endif

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        beta = 1d0
        var(:,:,:,IPRS) = 0.5*beta/a_p !Electron temperature required for force balance

        do k = klom,khip
          do j = jlom,jhip
             do i = ilom,ihip
              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy
     .                                      ,igz,ig,jg,kg,x1,y1,z1)
              
              r1 = sqrt(x1**2+y1**2)

              var(i,j,k,IRHO) = 1d0/(1d0+(r1/dlambda)**2)**2
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call find_polar_coords(gv%gparams,igr,i,j,k,x1,y1,z1
     .                              ,coords)

              btheta = (r1/dlambda)/(1+(r1/dlambda)**2)

              ! Calculate in Cartesian, then call XformVector_mesh2
              ! to calculate curvilinear and include the effects of packing.
              bb(i,j,k,1) =-btheta*sin(y1)
              bb(i,j,k,2) = btheta*cos(y1)
              bb(i,j,k,3) = bz0
            enddo
          enddo
        enddo

        call XformVector_mesh2(gv%gparams,igx,bb,'car','cnv')
#else
        find_a = fillVectorPotential(igx,aa)
#endif

      case ('2dtok')

c     Cylindrical Tokamak analytical equilibrium (Furth et al., 1973)

        mm = gv%gparams%params(1)
        kk = gv%gparams%params(2)
        RR = gv%gparams%params(3)
        q0 = gv%gparams%params(4)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        b_0 = dlambda/RR/q0

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        var(:,:,:,IPRS) = 1d-3

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call find_polar_coords(gv%gparams,igx,i,j,k,x1,y1,z1
     .                              ,coords)

              btheta = b_0*(x1/dlambda)/(1+(x1/dlambda)**2)
              bz0    = sqrt(1.-b_0**2*(1-1./(1+(x1/dlambda)**2)**2))

              !X-Y equilibrium
              bb(i,j,k,1)  = 0d0
              bb(i,j,k,2)  = btheta + kk*x1/mm*bz0
              bb(i,j,k,3)  = x1*bz0
            enddo
          enddo
        enddo

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

        mm = gv%gparams%params(1)
        kk = gv%gparams%params(2)
        RR = gv%gparams%params(3)

        nh2 = mm  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

        call find_global_radial_mesh(rg)

        Iz  = 2.22
        b_0 = 0.5*Iz
        do i=0,nxdp
          x1  = rg(i)
          bth0(i) = b_0*x1/(1+(x1/dlambda)**2)
          bzz0(i) = RR
          p00(i)  = a_p*(1d0+1./16.*(Iz*dlambda/(1+(x1/dlambda)**2))**2)
        enddo

        !Dump equilibrium for eigen code
        call dump_1d_eq
        
        !Build equilibrium
        var(:,:,:,IVX)  = 0d0
        var(:,:,:,IVY)  = 0d0
        var(:,:,:,IVZ)  = 0d0

        var(:,:,:,IRHO) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call find_polar_coords(gv%gparams,igr,i,j,k,x1,y1,z1
     .                                ,coords)

              jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

              iglobal = i + gv%gparams%ilo(igr) - 1

              bb(i,j,k,1)=0d0
              bb(i,j,k,2)=jac1*(bth0(iglobal)/x1+kk/mm*bzz0(iglobal))
              bb(i,j,k,3)=jac1*bzz0(iglobal)

              var(i,j,k,IPRS) = p00(iglobal)/a_p/var(i,j,k,IRHO)

            enddo
          enddo
        enddo

c$$$        do k = klom,khip
c$$$          do j = jlom,jhip
c$$$            do i = ilom,ihip
c$$$              call find_polar_coords(gv%gparams,igx,i,j,k,x1,y1,z1
c$$$     .                              ,coords)
c$$$
c$$$              jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)
c$$$
c$$$              btheta = b_0*x1/(1+(x1/dlambda)**2)
c$$$              bz0    = RR
c$$$
c$$$              bb(i,j,k,1)  = 0d0
c$$$              bb(i,j,k,2)  = jac1*(btheta/x1 + kk/mm*bz0)
c$$$              bb(i,j,k,3)  = jac1*bz0
c$$$
c$$$              var(i,j,k,IPRS) = 1d0
c$$$     .                    + 1./16.*(Iz*dlambda/(1+(x1/dlambda)**2))**2
c$$$            enddo
c$$$          enddo
c$$$        enddo

      case ('rfp83','vela1','vela2','tbbr1')

c     Check coordinates

        if (coords /= 'hel'.and.coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        mm   = gv%gparams%params(1) !Poloidal mode number
        kk   = gv%gparams%params(2) !Helical paramter (-nn/RR, nn-> toroidal mode number)
        RR   = gv%gparams%params(3) !Major radius (determines aspect ratio)
        bz0  = gv%gparams%params(4) !Toroidal field at magnetic axis
        beta = gv%gparams%params(5) !Beta toroidal

        if (nh2 == 0) nh2 = nint(mm)       !To set the right perturbation wavelength
        nh1 = nh2 !For perturbations to be of the right order

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

c$$$     RFP EQUILIBRIA (Caramana et al, PoP, 1983)

c$$$        bz0  = 1d0
c$$$        beta = 4d-5

        call find_global_radial_mesh(rg)

        do i=0,nxdp
          x1  = rg(i)
          p00(i) = 0.5*beta*bz0**2*pprof(x1,0)/pprof(0d0,0)
          q00(i) = RR*qprof(x1)
        enddo

        !Initialize toroidal field
        bzz0 = bz0  

        !Solve 1D Grad-Shafranov equation
        call GS_equil_1d(rg,q00,p00,bth0,bzz0)

        !Dump equilibrium for eigen code
        call dump_1d_eq

        !Build equilibrium
        var(:,:,:,IVX:IVZ)  = 0d0

        if (equil(1:4) == 'vela') then
          var(:,:,:,IRHO) = 1d0/(gamma-1d0)
        else
          var(:,:,:,IRHO) = 1d0
        endif

        select case(coords)
        case('hel')
          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                call find_polar_coords(gv%gparams,igr,i,j,k,x1,y1,z1
     .                                ,coords)

                jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

                iglobal = i + gv%gparams%ilo(igr) - 1

                bb(i,j,k,1)=0d0
                bb(i,j,k,2)=jac1*(bth0(iglobal)/x1+kk/mm*bzz0(iglobal))
                bb(i,j,k,3)=jac1*bzz0(iglobal)

                var(i,j,k,IPRS) = p00(iglobal)/a_p/var(i,j,k,IRHO)

              enddo
            enddo
          enddo
        case('cyl')
          if (nh3 == 0) nh3 = nint(kk*RR) !To set the right perturbation wavelength

          if (nh2 == 0 .and. nh3 == 0) then
            messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
            call pstop('setEquilibrium',messg)
          else
            messg=' Perturbation harmonics nh1, nh2, nh3 for '
     .          //trim(equil)//' equilibrium are'
            if (my_rank==0) then
              write (*,*)
              write (*,'(a,3i3)') trim(messg),nh1,nh2,nh3
            endif
          endif

          do k = klom,khip
            do j = jlom,jhip
              do i = ilom,ihip
                call find_polar_coords(gv%gparams,igr,i,j,k,x1,y1,z1
     .                                ,coords)

                jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

                iglobal = i + gv%gparams%ilo(igr) - 1

                bb(i,j,k,1)=0d0
                bb(i,j,k,2)=jac1*bth0(iglobal)/x1
                bb(i,j,k,3)=jac1*bzz0(iglobal)

                var(i,j,k,IPRS) = p00(iglobal)/a_p/var(i,j,k,IRHO)

              enddo
            enddo
          enddo
        end select
        
      case ('ntm2d')

        bootstrap= .true.

        mm   = gv%gparams%params(1) !Poloidal mode number
        kk   = gv%gparams%params(2) !Helical paramter (-nn/RR, nn-> toroidal mode number)
        RR   = gv%gparams%params(3) !Major radius (determines aspect ratio)
        bz0  = gv%gparams%params(4) !Toroidal field at magnetic axis
        beta = gv%gparams%params(5) !Beta toroidal

        nh2 = mm  !To set the right perturbation wavelength
        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Equilibria (numerical)

        call find_global_radial_mesh(rg)
        
        !Find pressure and q-profiles
        do i=0,nxdp
          x1  = rg(i)
          p00(i) = 0.5*beta*bz0**2*pprof(x1,0)/pprof(0d0,0)
          q00(i) = RR*qprof(x1)
        enddo

        !Initialize toroidal field
        bzz0 = bz0  

        !Solve 1D Grad-Shafranov equation
        call GS_equil_1d(rg,q00,p00,bth0,bzz0)

        !Build equilibrium

        var(:,:,:,IVX:IVZ) = 0d0

        var(:,:,:,IRHO) = 1d0

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call find_polar_coords(gv%gparams,igr,i,j,k,x1,y1,z1
     .                              ,coords)

              iglobal = i + gv%gparams%ilo(igx) - 1

              bb(i,j,k,1) = 0d0
              bb(i,j,k,2) = bth0(iglobal) + kk*x1/mm*bzz0(iglobal)
              bb(i,j,k,3) = x1*bzz0(iglobal)

              var(i,j,k,IPRS) = p00(iglobal)/a_p/var(i,j,k,IRHO)

            enddo
          enddo
        enddo

cc      case ('3dtok')
cc
ccc     Cylindrical Tokamak analytical equilibrium (Furth et al., 1973)
cc
cc        RR = gv%gparams%params(1)
cccc        mm = gv%gparams%params(4)
cccc        nn = gv%gparams%params(5)
cc        q0 = gv%gparams%params(6)
cc
cc        kk = 0d0   !To avoid affecting q-profile
cc
cc        b_0 = dlambda/RR/q0
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
cc        var(:,:,:,IPRS) = 1d-3
cc
cc        do k = klom,khip
cc          do j = jlom,jhip
cc            do i = ilom,ihip
cc
cc              call getCurvilinearCoordinates(i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                      ,x1,y1,z1)
cc
cc              btheta = b_0*(x1/dlambda)/(1+(x1/dlambda)**2)
cc              bz0    = sqrt(1.-b_0**2*(1-1./(1+(x1/dlambda)**2)**2))
cc
cc              call getCartesianCoordinates(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg
cc     .                                    ,x1,y1,z1)
cc
cc              !X-Y equilibrium
cc              bb(i,j,k,1)  = 0d0
cc              bb(i,j,k,2)  = btheta
cc              bb(i,j,k,3)  = gv%gparams%gmetric%grid(igx)%gsup(i,j,k,3,3)*bz0
cc            enddo
cc          enddo
cc        enddo

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

        if (coords /= 'ext') then
          messg = 'Need "ext" coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        if (.not.is_file(equ_file)) then
          call pstop('setEq','VMEC EQU file not found')
        endif

c     Read VMEC coordinates and fill grid hierarchy

cc        call vmec_map(equ_file,gv%gparams) !Moved to createNLF

c     Read VMEC equilibrium

        call vmec_equ(ilevel,igx,nx,ny,nz,bb,var(:,:,:,IPRS)
     .                                      ,var(:,:,:,IRHO)
     .               ,gamma,equ_file,dcon,div_clean)

        find_a = .false.

        !Divergence clean field if not Tokamak
#if !defined(vec_pot)
        if (div_clean) then
#endif
c$$$          call B_div_clean(gv%gparams,ilevel,igx,bb,"vlap"
c$$$     &                    ,rtol=1d-7,dacov=aa)
#if !defined(vec_pot)
        endif
#endif

c     Other equilibrium quantities

        var(:,:,:,IVX:IVZ)  = 0d0

        !Transform to e temperature
        where (var(:,:,:,IRHO) /= 0d0)
          var(:,:,:,IPRS) = var(:,:,:,IPRS)/var(:,:,:,IRHO)/a_p
        elsewhere
          var(:,:,:,IPRS) = 0d0
        end where
#endif

c     --------------------------------------------------
c     ---------------  EFIT EQUILIBRIUM ----------------
c     --------------------------------------------------

      case ('efit')    

        nh1 = nh2       !For perturbations to be of the right order in radius
                        !(nh2,nh3 set by input file)

        if (nh2 == 0 .and. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

        if (coords/='tor'.and.coords/='tsq'.and.coords/='sh3') then
          messg = 'Need tor or tsq coordinates for EFIT equilibrium '
     .            //equil
          call pstop('setEquilibrium',messg)
        endif

        if (.not.is_file(equ_file)) then
          call pstop('setEq','EFIT EQU file not found')
        endif

c     Read EFIT equilibrium

        call efit_equ(ilevel,igx,nx,ny,nz,bb,var(:,:,:,IPRS)
     .                                      ,var(:,:,:,IRHO)
     .               ,gamma,equ_file)

        find_a = .false.

        !Divergence clean field if not Tokamak
#if !defined(vec_pot)
        if (div_clean) then
#endif
c$$$          call B_div_clean(gv%gparams,ilevel,igx,bb,"vlap"
c$$$     &                    ,rtol=1d-7,dacov=aa)
#if !defined(vec_pot)
        endif
#endif

c     Other equilibrium quantities

        var(:,:,:,IVX:IVZ)  = 0d0

        !Transform to e temperature
        where (var(:,:,:,IRHO) /= 0d0)
          var(:,:,:,IPRS) = var(:,:,:,IPRS)/var(:,:,:,IRHO)/a_p
        elsewhere
          var(:,:,:,IPRS) = 0d0
        end where

c     ------------------------------------------
c     ---------------  RWM TEST ----------------
c     ------------------------------------------

      case ('rwhel')

        mm   = gv%gparams%params(1)
        kk   = gv%gparams%params(2)
        RR   = gv%gparams%params(3)
        beta = gv%gparams%params(4)
        prsw = gv%gparams%params(5)

        if (nh2 == 0) nh2 = nint(mm)  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel'.and.coords /= 'cyl') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

        !Amplitude of external B at bdry
        kn = abs(kk)
        if (kn /= 0d0) then
          pw = kn
          pv = kn*r_v
          Km   = BESSK(nh2  ,pw)                   
          Kmp  =-BESSK(nh2-1,pw) - nh2/pw*BESSK(nh2,pw)
          Kmpv =-BESSK(nh2-1,pv) - nh2/pv*BESSK(nh2,pv)
          Im   = BESSI(nh2  ,pw)
          Imp  = BESSI(nh2-1,pw) - nh2/pw*BESSI(nh2,pw)
          Impv = BESSI(nh2-1,pv) - nh2/pv*BESSI(nh2,pv)
          btho =-mm*Imp*(Km *Impv-Kmpv*Im )
     .                 /(Kmp*Impv-Kmpv*Imp)
          bzo  =-kk*Imp*(Km *Impv-Kmpv*Im )
     .                 /(Kmp*Impv-Kmpv*Imp)
        
          !Amplitude of internal B at bdry
          bthi =-mm*Im
          bzi  =-kk*Im
        else
          btho =-mm*(r_v**(2*mm)+1d0)/(1d0-r_v**(2*mm))
          bzo  = 0d0
        
          !Amplitude of internal B at bdry
          bthi =-mm
          bzi  = 0d0
        endif
c$$$        !Build vacuum potential
c$$$        do k = klom,khip
c$$$          do j = jlom,jhip
c$$$            do i = ilom,ihip
c$$$              call find_polar_coords(gv%gparams,igx,i,j,k,x1,y1,z1
c$$$     .                              ,coords)
c$$$              r = abs(x1)
c$$$              besselI = BESSI(nh2,kn*r)
c$$$              aa(i,j,k,1) = besselI*cos(y1)
c$$$            enddo
c$$$          enddo
c$$$        enddo

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call find_polar_coords(gv%gparams,igx,i,j,k,x1,y1,z1
     .                              ,coords)

              jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

              !Cnv helical B field
              r = abs(x1)
              if (kn /= 0d0) then
                besselIp = BESSI(nh2-1,kn*r)-nh2/(kn*r)*BESSI(nh2,kn*r)
                besselI  = BESSI(nh2  ,kn*r)
                br  = kn*besselIp*cos(y1)
                bth =-mm*besselI *sin(y1)/x1
                bz  =-kk*besselI *sin(y1)
              else
                br  = mm*x1**(mm-1)*cos(y1)
                bth =-mm*x1**(mm-1)*sin(y1)
                bz  = 0d0
              endif
              
              bb(i,j,k,1)=x1*br
              bb(i,j,k,2)=jac1*(bth/x1+kk/mm*bz)
              bb(i,j,k,3)=jac1*bz

c$$$              bb(i,j,k,:) = grad_ijk(gv%gparams,i,j,k,nx,ny,nz,igx,igx
c$$$     .                              ,igx,aa(:,:,:,1))
c$$$
c$$$              bb(i,j,k,:) = XformToCnv(gv%gparams,i,j,k,1,bb(i,j,k,:))

              !Boundary cyl electric field @ r=1
              Er  = 0d0
              Eth =-(bzo -bzi )/tau_w*sin(y1)
              Ez  = (btho-bthi)/tau_w*sin(y1)

              !Cyl pinch flow
              vr  = (Eth*bz-Ez*bth)/(br**2+bz**2+bth**2)
              vth = ( Ez*br       )/(br**2+bz**2+bth**2)
              vz  = (      -Eth*br)/(br**2+bz**2+bth**2)

              if (i==nx+1.and.k==1.and.(j>jlom.and.j<jhip)) then
c$$$                !Inside cyl B
c$$$                write (*,*) i,j,br,bth,bz
c$$$                !================
c$$$                if (j==1) write (*,*) "Inside helical cov B @ setEq"
c$$$                bth = x1*bth
c$$$                bz  = bz - kk/mm*bth
c$$$                write (*,*) i,j,bth,bz
c$$$                !================
c$$$                !Outer helical cov B
c$$$                if (j==1) write (*,*) "Outside helical cov B"
c$$$                write (*,*) i,j,btho*sin(y1),(bzo- kk/mm*btho)*sin(y1)
c$$$                !================
c$$$                !Helical cov E
c$$$                if (j==1) write (*,*) "Helical cov E"
c$$$                Ez  = Ez - kk/mm*Eth
c$$$                write (*,*) i,j,Eth,Ez
              endif

              var(i,j,k,IVX)=x1*vr
              var(i,j,k,IVY)=jac1*(vth/x1+kk/mm*vz)
              var(i,j,k,IVZ)=jac1*vz
cc              var(i,j,k,IVX:IVZ)=0d0

              var(i,j,k,IRHO) = 1d0

              var(i,j,k,IPRS) = prsw
            enddo
          enddo
        enddo

#if !defined(vec_pot)
        call B_div_clean(gv%gparams,ilevel,igx,bb,"vlap")
#endif

c     --------------------------------------------------
c     ---------------  OHMIC EQUILIBRIA ----------------
c     --------------------------------------------------

      case ('ppnch','ppnsl','ppnst')

        mm          = gv%gparams%params(1)
        kk          = gv%gparams%params(2)
        RR          = gv%gparams%params(3)
        theta_ppnch = gv%gparams%params(4)
        e0_over_eta = gv%gparams%params(5)
        beta        = gv%gparams%params(6)
        m0_over_nu  = gv%gparams%params(7)
        prsw        = gv%gparams%params(8)

        if (nh2 == 0) nh2 = nint(mm)  !To set the right perturbation wavelength

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

          call find_global_radial_mesh(rg)

          !Find pressure profile
          do i=0,nxdp
            x1  = rg(i)
            p00(i) = 0.5*beta*bz0**2*pprof(x1,0)
          enddo

	  !Initialize toroidal field
          bzz0 = bz0

! Daniele, 11/10/2013
          if (gamma > 1d0) then
             call Ohm_equil_1d_xport(-kk,RR,e0_over_eta,theta_ppnch,prsw
     .                              ,p00,bth0,bzz0,vr0,psi0)
          else
             call Ohm_equil_1d      (-kk,RR,e0_over_eta,theta_ppnch
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
          vz0 = m0_over_nu*0.25*(1-rg**2)

          !Find momentum source (cylindrical)
          M0(1:2) = 0d0
          M0(3) = m0_over_nu*nu

	  !Find density
          rho0 = 1d0

        endif

        !Convert electric field to covariant at r=1d0
        E0(3) = -kk/mm*E0(2) + E0(3)

        !Build equilibrium
        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call find_polar_coords(gv%gparams,igx,i,j,k,x1,y1,z1
     .                              ,coords)

              jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

              iglobal = i + gv%gparams%ilo(igx) - 1

              bb(i,j,k,1)=0d0
              bb(i,j,k,2)=jac1*(bth0(iglobal)/x1+kk/mm*bzz0(iglobal))
              bb(i,j,k,3)=jac1*bzz0(iglobal)

              var(i,j,k,IVX)=x1*vr0(iglobal)*eta/eta_in
              var(i,j,k,IVY)=kk*jac1/mm*vz0(iglobal)
              var(i,j,k,IVZ)=jac1*vz0(iglobal)

              if (solve_rho) then
                var(i,j,k,IRHO) = rho0(iglobal)
              else
                var(i,j,k,IRHO) = 1d0
              endif

              var(i,j,k,IPRS) = p00(iglobal)/a_p/var(i,j,k,IRHO)

cc              psi(i,j,k) = psi0(iglobal)  !Store flux function

            enddo
          enddo
        enddo

      case ('2dhel')

        mm = gv%gparams%params(1)
        kk = gv%gparams%params(2)
        RR = gv%gparams%params(3)
        theta_ppnch = gv%gparams%params(4)
        e0_over_eta = gv%gparams%params(5)
        beta = gv%gparams%params(6)

        if (nh2 == 0) nh2 = nint(mm)  !To set the right perturbation wavelength

        nh1 = nh2 !For perturbations to be of the right order

c     Check coordinates

        if (coords /= 'hel') then
          messg = 'Wrong coordinates for equilibrium '//equil
          call pstop('setEquilibrium',messg)
        endif

c     Initial guess: 1D equilibrium

        bb(:,:,:,1) = 0d0

        bzz0 = 1d0
        p00  = 0d0

        if (gamma > 1d0) then
           call Ohm_equil_1d_xport(-kk,RR,e0_over_eta,theta_ppnch,0d0
     .                            ,p00,bth0,bzz0,vr0,psi0)
        else
           call Ohm_equil_1d      (-kk,RR,e0_over_eta,theta_ppnch
     .                            ,p00,bth0,bzz0,vr0,psi0)
        endif

        do i=0,nx+1
          iglobal = i + gv%gparams%ilo(igx) - 1
          bb(i,:,:,2) = bth0(iglobal)
          bb(i,:,:,3) = bzz0(iglobal)
        enddo

c     Helical equilibria

	!Solve for helical eq (returns B and V in cnv representation)
        call Ohm_hel_equil(nx,ny,nz,-kk/mm,nh2,e0_over_eta,theta_ppnch
     .                    ,psi0,bb,var(:,:,:,IVX:IVY))

        var(:,:,:,IVX:IVY) = var(:,:,:,IVX:IVY)*eta
        var(:,:,:,IVZ) = 0d0

	!Find electric field
        E0(1:2) = 0d0
        E0(3) = e0_over_eta*eta

        !Convert current to covariant at r=1d0
        E0(3) = -kk/mm*E0(2) + E0(3)

c     Density and temperature

        var(:,:,:,IRHO) = 1d0
        var(:,:,:,IPRS) = 1d-5

c     Magnetic field

        do k = 0,nz+1
          do j = 0,ny+1
            do i = 0,nx+1
              call find_polar_coords(gv%gparams,igx,i,j,k,x1,y1,z1
     .                              ,'cyl')

              jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

              iglobal = i + gv%gparams%ilo(igx) - 1

              bb(i,j,k,1) = 0d0
              bb(i,j,k,2) = jac1*(bth0(iglobal)/x1+kk/mm*bzz0(iglobal))
              bb(i,j,k,3) = jac1*bzz0(iglobal)
            enddo
          enddo
        enddo

      case ('ppn3d','p3nsl')

        mm = gv%gparams%params(1)
        kk = gv%gparams%params(2)
        RR = gv%gparams%params(3)
        theta_ppnch = gv%gparams%params(4)
        e0_over_eta = gv%gparams%params(5)
        beta        = gv%gparams%params(6)
        m0_over_nu  = gv%gparams%params(7)
        prsw = gv%gparams%params(8)

        if (nh2 == 0) nh2 = nint(mm)       !To set the right perturbation wavelength
        if (nh3 == 0) nh3 = nint(kk*RR)    !To set the right perturbation wavelength

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

        if (abs(RR-zmax/(2.*pi)) > 1d-5.and.nzd > 1) then
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

          call find_global_radial_mesh(rg)

          !Find pressure and q-profiles
          do i=0,nxdp
            x1  = rg(i)
            p00(i) = 0.5*beta*bz0**2*pprof(x1,0)
          enddo

	  !Initialize toroidal field
          bzz0 = bz0

          if (gamma > 1d0) then
             call Ohm_equil_1d_xport(-kk,RR,e0_over_eta,theta_ppnch,prsw
     .                              ,p00,bth0,bzz0,vr0,psi0)
          else
             call Ohm_equil_1d      (-kk,RR,e0_over_eta,theta_ppnch
     .                              ,p00,bth0,bzz0,vr0,psi0)
          endif

          !Find electric field
          E0(1:2) = 0d0
          E0(3)   = e0_over_eta*eta
          eta_in  = eta

          !Adjust radial flow to actual resistivity
          vr0 = eta*vr0

          !Find axial flow from constant momentum source (Dirichlet at wall)
          vz0 = m0_over_nu*0.25*(1-rg**2)

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
            x1 = rg(i) !Global coordinate
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
              call find_polar_coords(gv%gparams,igr,i,j,k,x1,y1,z1
     .                              ,coords)

              jac1 = gv%gparams%gmetric%grid(igx)%jac(i,j,k)

              iglobal = i + gv%gparams%ilo(igx) - 1

              bb(i,j,k,1) = 0d0
              bb(i,j,k,2) = jac1*bth0(iglobal)/x1
              bb(i,j,k,3) = jac1*bzz0(iglobal)

              var(i,j,k,IVX) = x1*vr0(iglobal)*eta/eta_in
              var(i,j,k,IVY) = 0d0
              var(i,j,k,IVZ) = jac1*vz0(iglobal)

              if (solve_rho) then
                var(i,j,k,IRHO) = rho0(iglobal)
              else
                var(i,j,k,IRHO) = 1d0
              endif

              if (equil == 'ppn3d') then
                var(i,j,k,IPRS) = 1d-5
              else
                var(i,j,k,IPRS) = p00(iglobal)/a_p/var(i,j,k,IRHO)
              endif

cc              psi(i,j,k) = psi0(iglobal)  !Store flux function

            enddo
          enddo
        enddo

      case ('ohtor')

        !Parameters 2 and 3 define elliptical shape; used in grid_anal_map
        RR  = gv%gparams%params(1)            !Major radius
        a1  = gv%gparams%params(2)            !Horizontal elliptical radius
        b1  = gv%gparams%params(3)            !Vertical   elliptical radius
        mm  = gv%gparams%params(4)            !m mode
        nn  = gv%gparams%params(5)            !n mode
        e0_over_eta = gv%gparams%params(6)*RR !Covariant component is constant, needs RR factor
        beta= gv%gparams%params(7)            !Toroidal beta
        prsw= gv%gparams%params(8)            !Pressure at wall
        
        kk = 0d0   !To avoid affecting q-profile

        !Now these must be initialized directly in input file.
        if (nh2 == 0) nh2 =  nint(mm)       !To set the right perturbation wavelength
        if (nh3 == 0) nh3 = -nint(nn)       !To set the right perturbation wavelength

        nh1 = nh2       !For perturbations to be of the right order

        if (nh2 == 0 .and. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

        if (eta == 0d0) then
          messg='Need finite resistivity for this equilibrium'
          call pstop('setEquilibrium',messg)
        endif

c     Check coordinates

        if (      coords /= 'tor'
     .      .and. coords /= 'sha' 
     .      .and. coords /= 'sh2') then
          messg = 
     .         'Need "tor", "sha" or "sh2" coordinates for equilibrium '
     .           //equil
          call pstop('setEquilibrium',messg)
        endif

c     Ohmic GS equilibrium

        !Build equilibrium
        call Ohm_GS_equil_2d(RR,e0_over_eta,beta,nx,ny,nz
#if defined(vec_pot)
     .                   ,aa(:,:,:,1),aa(:,:,:,2),aa(:,:,:,3)
#else
     .                   ,bb(:,:,:,1),bb(:,:,:,2),bb(:,:,:,3)
#endif
     .                   ,var(:,:,:,IVX),var(:,:,:,IVY),var(:,:,:,IPRS))

        find_a = .false.

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
        var(:,:,:,IPRS) = var(:,:,:,IPRS) + prsw
        var(:,:,:,IPRS) = var(:,:,:,IPRS)/a_p/var(:,:,:,IRHO)

      case ('stell')

        !Parameters 2 and 3 define elliptical shape; used in grid_anal_map
        RR  = gv%gparams%params(1)            !Major radius
        a1  = gv%gparams%params(2)            !Horizontal elliptical radius
cc        bb  = gv%gparams%params(3)            !Vertical   elliptical radius
        mm  = gv%gparams%params(4)            !m mode
        nn  = gv%gparams%params(5)            !n mode
        e0_over_eta = gv%gparams%params(6)*RR !Covariant component is constant, needs RR factor
        beta= gv%gparams%params(7)            !Toroidal beta
        prsw= gv%gparams%params(8)            !Pressure at wall
        
        kk = 0d0   !To avoid affecting q-profile

        !Now these must be initalized directly in input file.
        if (nh2 == 0) nh2 =  nint(mm)       !To set the right perturbation wavelength
        if (nh3 == 0) nh3 = -nint(nn)       !To set the right perturbation wavelength

        nh1 = nh2       !For perturbations to be of the right order

        if (nh2 == 0 .and. nh3 == 0) then
          messg='Perturbation harmonics nh2, nh3 for '//trim(equil)
     .          //' equilibrium not set'
          call pstop('setEquilibrium',messg)
        endif

c     Check coordinates

        if (coords /= 'stl') then
          messg = 'Need "stl" coordinates for equilibrium '
     .           //equil
          call pstop('setEquilibrium',messg)
        endif

c$$$c     Ohmic GS equilibrium
c$$$
c$$$        !Build equilibrium
c$$$        call Ohm_GS_equil_2d(RR,e0_over_eta,beta,nx,ny,nz
c$$$#if defined(vec_pot)
c$$$     .                   ,var(:,:,:,IAX),var(:,:,:,IAY),var(:,:,:,IAZ)
c$$$#else
c$$$     .                   ,bb(:,:,:,1),bb(:,:,:,2),bb(:,:,:,3)
c$$$#endif
c$$$     .                   ,var(:,:,:,IVX),var(:,:,:,IVY),var(:,:,:,IPRS))
c$$$        
c$$$	!Find electric field
c$$$        E0(1:2) = 0d0
c$$$        E0(3)   = e0_over_eta*eta
c$$$        eta_in  = eta
c$$$
c$$$        !Adjust pinch flow to current resistivity
c$$$        var(:,:,:,IVX) = var(:,:,:,IVX)*eta
c$$$        var(:,:,:,IVY) = var(:,:,:,IVY)*eta
c$$$        var(:,:,:,IVZ) = 0d0
c$$$
c$$$        !Fill other arrays
c$$$        var(:,:,:,IRHO) = 1d0
c$$$        var(:,:,:,IPRS) = var(:,:,:,IPRS) + prsw
c$$$        var(:,:,:,IPRS) = var(:,:,:,IPRS)/a_p/var(:,:,:,IRHO)

        var(:,:,:,IRHO) = 1d0
        var(:,:,:,IPRS) = 0d0
        bb = 0d0
        var(:,:,:,IVX) = 0d0
        var(:,:,:,IVY) = 0d0
        var(:,:,:,IVZ) = 0d0

      case default

        messg = 'Equilibrium '//trim(equil)//' undefined'
        call pstop('setEquilibrium',messg)

      end select

c Transfer magnetic field/vector potential

#if defined(vec_pot)
      if (find_a) then
        call findVectorPotential(gv%gparams,ilevel,igx,bb,aa,'vlap')
      endif
      var(:,:,:,IAX:IAZ) = aa
#else
      var(:,:,:,IBX:IBZ) = bb
#endif

c Find total pressure

      var(:,:,:,IPRS) = a_p*var(:,:,:,IPRS)*var(:,:,:,IRHO)

c Find momentum components

      if (.not.nc_eom_v) then
        var(:,:,:,IVX) = var(:,:,:,IRHO)*var(:,:,:,IVX)
        var(:,:,:,IVY) = var(:,:,:,IRHO)*var(:,:,:,IVY)
        var(:,:,:,IVZ) = var(:,:,:,IRHO)*var(:,:,:,IVZ)
      endif

c Apply Jacobian factors

      var(:,:,:,IRHO) =var(:,:,:,IRHO)*gv%gparams%gmetric%grid(iigx)%jac
      var(:,:,:,IPRS) =var(:,:,:,IPRS)*gv%gparams%gmetric%grid(iigx)%jac

      if (alt__eom()) then
        var(:,:,:,IVY) =var(:,:,:,IVY)*gv%gparams%gmetric%grid(iigx)%jac
      endif

c Transfer equilibrium to varray

      do ieq=1,neqd
        varray%array_var(ieq)%array = var(:,:,:,ieq) 
      enddo

c Transfer BC fields info to postprocessor via input file

      call write_bc_field_to_inputf(E0,'E0')
      call write_bc_field_to_inputf(B0,'B0')
      call write_bc_field_to_inputf(M0,'M0')

c End program
      
#if defined(RFX)
10    format (4f7.3,2e12.3)
#endif

      contains

c     dump_1d_eq
c     ###############################################################
      subroutine dump_1d_eq

c     ---------------------------------------------------------------
c     Dump 1d equilibrium for eigen analysis
c     ---------------------------------------------------------------

      implicit none

      if (my_rank == 0) then
        open(unit=123,file=trim(equil)//'-eq.txt',status='unknown')
        write (123,*) nxdp+1
        do i=0,nxdp
          write (123,*) rg(i),bth0(i),bzz0(i),p00(i)/a_p
        enddo
        write (123,*) gv%gparams%params(3),eq_params(2),eq_params(3)
        close(123)
      endif

      end subroutine dump_1d_eq
      
      end subroutine setEquilibrium

c     find_global_radial_mesh
c     ###############################################################
      subroutine find_global_radial_mesh(rg)

c     ---------------------------------------------------------------
c     Find radial global-mesh 
c     ---------------------------------------------------------------

      use parameters

      use equilibrium

      implicit none

c     Call variables
      
      real(8) :: rg(0:nxdp)

c     Local variables

      integer :: i
      real(8) :: rl(nxl),th,zz

c     Begin program

      do i=1,nxl
        call find_polar_coords(gv%gparams,1,i,1,1,rl(i),th,zz,coords)
      enddo

#if defined(petsc)
      call MPI_Allgather(rl(1:nxl),nxl,MPI_DOUBLE_PRECISION
     .                  ,rg(1:nxd),nxl,MPI_DOUBLE_PRECISION
     .                  ,gv%gparams%MPI_COMM_X,mpierr)
#else
      rg(1:nxd) = rl
#endif

      rg(0)    = -rg(1)         !@ SP
      rg(nxdp) = 2*xmax - rg(nxd) !@ radial boundary

c     End program

      end subroutine find_global_radial_mesh
      
c     GS_equil_1d
c     ###############################################################
      subroutine GS_equil_1d(rg,Rxq0,p0,bth0,bzz0)

c     ---------------------------------------------------------------
c     Solves 1D Grad-Shafranov equation to find cylindrical components
c     of magnetic field B_theta, B_z from specified q-profile and
c     pressure profile. On input:
c       * rg: global (possibly non-uniform) radial mesh
c       * Rxq0,p0: normalized q-profile (RR x q) and pressure profile
c       * bth0,bzz0: output magnetic field components (bzz0 contains
c                    initial guess, and bzz0(0) contains value at axis)
c     ---------------------------------------------------------------

        use parameters

        use equilibrium

        use grid

        use app_iosetup

        use nlfunction_setup, ONLY: a_p

        implicit none

c     Call variables

        real(8)    :: bth0(0:nxdp),bzz0(0:nxdp)
     .               ,p0(0:nxdp),Rxq0(0:nxdp),rg(0:nxdp)

c     Local variables

        integer    :: i,j,k,ig,jg,kg,it
        real(8)    :: bz_avg,bz0,x1,bzh,bth,gs_err,dh,rr

        real(8)    :: bth_f(0:nxdp)
     .               ,xg(0:nxdp),dr(0:nxd),dx(0:nxd)
     .               ,p0h(0:nxd),q0h(0:nxd),rh(0:nxd)
     .               ,xh(0:nxd),jach(0:nxd),dummy(0:nxdp)

        real(8),parameter :: cm=-1d0/12d0       !Upwinded
     .                      ,c0=2d0/3d0 
     .                      ,cp=5d0/12d0

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '************************************'
          write (*,*) '1D Grad-Shafranov equilibrium solver'
          write (*,*) '************************************'
          write (*,*)
        endif

c     Initialize mesh qtys

        xg = gv%gparams%xg
        
        !Find half-mesh positions
        do i=0,nxd
          rh(i) = 0.5*(rg(i+1)+rg(i))
          dr(i) =      rg(i+1)-rg(i)
          xh(i) = 0.5*(xg(i+1)+xg(i))
          dx(i) =      xg(i+1)-xg(i)
        enddo

c     Initialize iteration

        !Average pressure and q-profiles to half mesh
c$$$        do i=0,nxd
c$$$          p0h(i) = 0.5*(p0(i)  +p0  (i+1)) !Average pressure  to half mesh
c$$$          q0h(i) = 0.5*(Rxq0(i)+Rxq0(i+1)) !Average q-profile to half mesh
c$$$        enddo
        call IntDriver1d(nxd+2,xg,p0  ,nxd+1,xh(0:nxd),p0h,3,0)
        call IntDriver1d(nxd+2,xg,Rxq0,nxd+1,xh(0:nxd),q0h,3,0)
        call IntDriver1d(nxd+2,xg,rg  ,nxd+1,xh(0:nxd),jach,3,1)

        bz0 = bzz0(0)

c     Iteration

        do it=1,100

cc          call IntDriver1d(nxd+1,rh(0:nxd),bth0 (0:nxd)
cc     .                    ,nxd  ,rg(1:nxd),bth_f(1:nxd),3,0)

          !Perform integral of B_th^2/r on half-mesh (3rd order)
          dummy(0) = 0d0
          do i=1,nxd
cc            dummy(i) = dummy(i-1) + 2*bth_f(i)**2/rg(i)*jach(i)*dx(i)
cc            dummy(i) = dummy(i-1) + 2*bth0(i-1)*bth0(i)/rg(i)*dr(i)
            if (i == 1) then
              dummy(i) = dummy(i-1)
     .                 + 0d0
     .                 + 2*c0*bth0(i  )**2/rh(i  )*jach(i  )*dx(i  )
     .                 + 2*cm*bth0(i+1)**2/rh(i+1)*jach(i+1)*dx(i+1)
            elseif (i == 2) then
              dummy(i) = dummy(i-1)
     .                 + 0d0
     .                 + 2*c0*bth0(i-1)**2/rh(i-1)*jach(i-1)*dx(i-1)
     .                 + 2*cp*bth0(i  )**2/rh(i  )*jach(i  )*dx(i  )
            else
              dummy(i) = dummy(i-1)
     .                 + 2*cm*bth0(i-2)**2/rh(i-2)*jach(i-2)*dx(i-2)
     .                 + 2*c0*bth0(i-1)**2/rh(i-1)*jach(i-1)*dx(i-1)
     .                 + 2*cp*bth0(i  )**2/rh(i  )*jach(i  )*dx(i  )
            endif
          enddo

          !Safeguard for solvability
          where (dummy(0:nxd)>bz0**2+2*(p0h(0)-p0h))
     .           dummy(0:nxd)=bz0**2+2*(p0h(0)-p0h)

          !Calculate Bz, Btheta on half-mesh
          do i=0,nxd
            x1  = rh(i)   !Half mesh
            bzz0(i) = sign(1d0,q0h(i))
     .               *sqrt((bz0**2+2*(p0h(0)-p0h(i))-dummy(i))
     .                    /((x1/q0h(i))**2+1))
          enddo

          !Check convergence
          gs_err = sqrt(sum(bth0(0:nxd)-rh*bzz0(0:nxd)/q0h(0:nxd))**2)  !RMS error
          
          if (ilevel > 0 .and. my_rank == 0) then
            write (*,*) '1D Grad-Shafranov equilibrium iter =',it
     .                 ,' ; Force balance error=',gs_err
          endif

          !Updated poloidal B-field component
          bth0(0:nxd) = rh*bzz0(0:nxd)/q0h

          if (gs_err < 1d-7) exit

        enddo

        if (ilevel == 0 .and. my_rank == 0) then
           write (*,'(a,i3,a,1pe10.2)')
     .          ' 1D equilibrium converged in '
     .          ,it,' iterations with force-balance error =',gs_err
        endif

        !Check force balance
        do i=1,nxd
          rr  = rg(i)

          bzh = 0.5*(bzz0(i)+bzz0(i-1))
          bth = 0.5*(bth0(i)+bth0(i-1))
          gs_err = gs_err
     .              +(bzh*(bzz0(i)-bzz0(i-1))/dr(i)
     .                   +(p0h (i)-p0h (i-1))/dr(i)
     .               +0.5*(bth0(i)**2-bth0(i-1)**2)/dr(i)
     .               +bth**2/rr)**2
        enddo
        gs_err = sqrt(gs_err/nxd)

        if (my_rank == 0) then
          write(*,*)
          write(*,'(a,1pe14.7)') ' Force-balance error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
cc          write (*,*) dummy(1:nxd)
        endif

        !Calculate equilibrium properties
        bz_avg = 0d0
        do i=1,nxd
          x1  = rh(i)
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
        eq_params(1) = 1d0/gv%gparams%params(3)       !Aspect ratio
        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)        !Poloidal beta
        eq_params(3) = q0h(nxd)/gv%gparams%params(3)  !q(r=a)

        !Transfer to integer mesh
        call IntDriver1d(nxdp,rh(0:nxd),bth0(0:nxd)
     .                  ,nxdp,rg(1:nxdp),dummy(1:nxdp),3,0)
        bth0 = dummy

        call IntDriver1d(nxdp,rh(0:nxd),bzz0(0:nxd)
     .                  ,nxdp,rg(1:nxdp),dummy(1:nxdp),3,0)
        bzz0 = dummy

        !Extrapolate linearly in radius @ r=0
        bth0(0) =-bth0(1)
        bzz0(0) = bzz0(1)

      end subroutine GS_equil_1d

c     Ohm_equil_1d
c     ###############################################################
      subroutine Ohm_equil_1d(kk,RR,Eoeta,theta,p0,bth0,bzz0,vr0,psi0)

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

        real(8),intent(IN) :: kk,RR,Eoeta,theta,p0(0:nxdp)

        real(8),intent(OUT):: bth0(0:nxdp),bzz0(0:nxdp)
     .                       ,vr0 (0:nxdp),psi0(0:nxdp)

c     Local variables

        integer :: i,j,k,ig,jg,kg,it
        real(8) :: bz0,bz_avg,x1,bzh,bth,err,rh,dh,p00,th,zz
     .            ,dx1,dx2,gp,dvol

        real(8) :: dummy(0:nxdp),bth0_old(0:nxdp)
     .            ,bzz0_old(0:nxdp),xh(0:nxd),dr(0:nxd),rg(0:nxdp)
     .            ,rl(nxl),drh(0:nxd),eta0(0:nxdp)

        real(8) :: pprof
        external   pprof

c     Begin program

        if (my_rank == 0) then
          write (*,*)
          write (*,*) '***************************'
          write (*,*) '1D Ohmic equilibrium solver'
          write (*,*) '***************************'
          write (*,*)
        endif

c       Initialize iteration

        bz0 = bzz0(0)
        p00 = maxval(p0)

c       Find global-mesh radial positions and increments

        !Integer mesh positions
        do i=1,nxl
          call find_polar_coords(gv%gparams,1,i,1,1,rl(i),th,zz
     .                          ,coords)
        enddo

#if defined(petsc)
        call MPI_Allgather(rl(1:nxl),nxl,MPI_DOUBLE_PRECISION
     .                    ,rg(1:nxd),nxl,MPI_DOUBLE_PRECISION
     .                    ,gv%gparams%MPI_COMM_X,mpierr)
#else
        rg(1:nxd) = rl
#endif
        rg(0)    = -rg(1)           !@ SP
        rg(nxdp) = 2*xmax - rg(nxd) !@ radial boundary

        !Half-mesh positions and spacings
        do i=0,nxd
          xh (i) = 0.5*(rg(i+1)+rg(i))
          drh(i) =      rg(i+1)-rg(i)
        enddo

        !Integer mesh spacings
        do i=1,nxd
          dr(i) =      xh(i)-xh(i-1)
        enddo
        dr(0) = dr(1)

        !Initial guess for B_theta (defined in half mesh)
        bth0 = 0d0
        do i=0,nxd
          bth0(i) = 0.5*xh(i)*Eoeta
        enddo

        !Resistivity profile
        do i=1,nxdp
          if (spitzer) then
            eta0(i) = res_spitzer(p0(i)/p0(1))
          else
            eta0(i) = res_rfx(rg(i))
          endif
        enddo

c       Iteration

        do it=1,100

          bth0_old = bth0
          bzz0_old = bzz0

          !Perform integral of B_theta
          dummy(0) = 0d0
          bth0 (0) = 0d0
          do i=1,nxd
            x1 = xh(i)  !Half mesh
            rh = rg(i)  !Full mesh

            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))

            !Second order pressure gradient
            dx1=drh(i-1)
            dx2=drh(i)

            gp = ( p0(i+1)*dx1/dx2/(dx1+dx2)
     .            +p0(i  )*(dx2-dx1)/dx2/dx1
     .            -p0(i-1)*dx2/dx1/(dx1+dx2))

            !Accumulate integral on half mesh
            dvol = 0.5*(xh(i)**2-xh(i-1)**2)
            dummy(i)=dummy(i-1)
     .              +(Eoeta/eta0(i)*bzh**2-bth*gp)/(bzh**2+bth**2)*dvol

            !Transform to Bth on half-mesh
            bth0(i) = dummy(i)/xh(i)
          enddo

          !Perform integral of Bz
          dummy = 0d0
          do i=1,nxd
            x1 = xh(i)  !Half mesh
            rh = rg(i)  !Full mesh

            bzh = 0.5*(bzz0(i)+bzz0(i-1))
            bth = 0.5*(bth0(i)+bth0(i-1))

            !Second order pressure gradient
            dx1 = drh(i-1)
            dx2 = drh(i)

            gp = ( p0(i+1)*dx1/dx2/(dx1+dx2)
     .            +p0(i  )*(dx2-dx1)/dx2/dx1
     .            -p0(i-1)*dx2/dx1/(dx1+dx2))

            !Accumulate integral
            dummy(i) = dummy(i-1)
     .               +(Eoeta/eta0(i)*bth+gp)/(bzh**2+bth**2)*dr(i)
          enddo

          !Transform to Bz
          bzz0 = exp(-dummy)

          if (theta /= 0d0) then
            bz_avg = 0d0
            do i=1,nxd
               bz_avg = bz_avg + (bzz0(i)*xh(i)+bzz0(i-1)*xh(i-1))*dr(i)
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
          rh  = rg(i)

          !Second order pressure gradient
          dx1 = drh(i-1)
          dx2 = drh(i)

          gp = ( p0(i+1)*dx1/dx2/(dx1+dx2)
     .          +p0(i  )*(dx2-dx1)/dx2/dx1
     .          -p0(i-1)*dx2/dx1/(dx1+dx2))

          bzh = 0.5*(bzz0(i)+bzz0(i-1))
          bth = 0.5*(bth0(i)+bth0(i-1))
          dummy(i) = bzh*(bzz0(i)-bzz0(i-1))/dr(i)
     $              +bth/rh*(bth0(i)*xh(i)-bth0(i-1)*xh(i-1))/dr(i)
     .              +gp
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
          rh  = rg(i)

          dummy(i) = Eoeta*bzh/eta0(i)
     $              -bzh/rh*(bth0(i)*xh(i)-bth0(i-1)*xh(i-1))/dr(i)
     $              +bth*(bzz0(i)-bzz0(i-1))/dr(i)
        enddo

        if (my_rank == 0) then
          write(*,'(a,1pe14.7)') ' Ohms law      error ='
     .                        ,sqrt(sum(dummy(1:nxd)**2)/nxd)
        endif

        !Calculate pinch equilibrium features (theta parameter, Bz(r=0), I0)
        bz_avg = 0d0
        do i=1,nxd
          bz_avg = bz_avg + (bzz0(i)*xh(i)+bzz0(i-1)*xh(i-1))*dr(i)
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
cc        eq_params(1) = 1d0/gv%gparams%params(3)       !Aspect ratio
cc        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
cc        eq_params(3) = q0h(nxd)/gv%gparams%params(3)  !q(r=a)

        !Find helical flux (integer mesh)
        psi0(0) = 0d0
        do i=1,nxd+1
          psi0(i) = psi0(i-1)+dr(i-1)*(bth0(i-1)-kk*xh(i-1)*bzz0(i-1))
        enddo
        psi0 = psi0 - 0.5*(psi0(nxd)+psi0(nxdp)) !Homogeneous Dirichlet BC

        !Transfer B and V to integer mesh
        call IntDriver1d(nxdp,xh(0:nxd),bth0(0:nxd)
     .                  ,nxdp,rg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bth0 = dummy

        call IntDriver1d(nxdp,xh(0:nxd),bzz0(0:nxd)
     .                  ,nxdp,rg(1:nxdp),dummy(1:nxdp)
     .                  ,2,0)
        bzz0 = dummy

        !Extrapolate linearly in radius @ r=0
        bth0(0) = -bth0(1)
        bzz0(0) = bzz0(1)

        !Find pinch flow at cell centers (normalized to nominal resistivity)
        vr0 = 0d0

        do i=1,nxdp
          rh  = rg(i)

          dx1 = drh(i-1)

          if (i < nxdp) then
            !Second order pressure gradient
            dx2 = drh(i)

            gp = ( p0(i+1)*dx1/dx2/(dx1+dx2)
     .            +p0(i  )*(dx2-dx1)/dx2/dx1
     .            -p0(i-1)*dx2/dx1/(dx1+dx2))
          else
            gp = (p0(i)-p0(i-1))/dx1
          endif

          vr0(i) = -(Eoeta*bth0(i)+eta0(i)*gp)/(bzz0(i)**2+bth0(i)**2)
        enddo

        !i=0
        vr0(0) = -vr0(1)

        !Find pressure
cc        do i=1,nxd
cc          rh  = rg(i)
cc          p0(i) = p00*pprof(rh,0)
cc        enddo

        !Dump equilibrium for eigen code
cc        if (my_rank == 0) then
cc          open(unit=123,file='ntm2d.txt',status='unknown')
cc          write (123,*) nxdp+1
cc          write (123,*) rg(0:nxdp)
cc          write (123,*) bth0
cc          write (123,*) bzz0
cc          write (123,*) 0.5*p0
cc          write (123,*) gv%gparams%params(3),eq_params(2),eq_params(3)
cc          close(123)
cccc          stop
cc        endif

#if defined(RFX)
        if (my_rank == 0.and.(.not.test)) then
          write (*,*) 'r,bth0,bzz0,q0,vr0,p0'
          do i=0,nxdp
            write(*,10) rg(i),bth0(i),bzz0(i)
     .             ,(rg(i)*bzz0(i))/(RR*bth0(i)),vr0(i),p0(i)
          enddo
        endif

 10     format (4f7.3,2e12.3)
#endif

      end subroutine Ohm_equil_1d

c     Ohm_equil_1d_xport
c     ###############################################################
      subroutine Ohm_equil_1d_xport(kk,RR,Eoeta,theta,pw,p0,bth0,bzz0
     .                             ,vr0,psi0)

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
! Daniele, 12 Nov. 2014
     .                             ,ext_heat_src_prof

        use grid

        use app_iosetup

        implicit none

c     Call variables

        real(8) :: kk,bth0(0:nxdp),bzz0(0:nxdp),vr0(0:nxdp)
     .            ,p0(0:nxdp),Eoeta,theta,b_tor,psi0(0:nxdp)
     .            ,pw,RR

c     Local variables

        integer :: i,j,k,ig,jg,kg,it,info
        real(8) :: bz_avg,x1,bzh,bth,err,rh,dh,p00,rhm,rhp
! Daniele, 02/07/2014
        real(8) :: chirh,chirhm,chirhp

        real(8) :: bz0,dummy(0:nxdp),bth0_old(0:nxdp)
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
          dr(i) =      gv%gparams%xg(i+1)-gv%gparams%xg(i)
          xh(i) = 0.5*(gv%gparams%xg(i+1)+gv%gparams%xg(i))
! Daniele, 02/07/2014
c$$$!          chih(i)=chi*(1.-(1.-1d-2)*sin(pi*xh(i)/1.75))
c$$$          chih(i)=dfvty_rfx(xh(i))
! Luis 04/15/2014, to mimic what is done in PIXIE3D
          chih(i)=0.5*(dfvty_rfx(gv%gparams%xg(i+1))
     $                +dfvty_rfx(gv%gparams%xg(i  )))
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
              rh = gv%gparams%xg(i)
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
            rh = gv%gparams%xg(i)

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
            rh  = gv%gparams%xg(i)

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
     .                    ,nxdp,gv%gparams%xg(1:nxdp),dummy(1:nxdp)
     .                    ,2,0)
          btf = dummy

          call IntDriver1d(nxdp,xh(0:nxd),bzz0(0:nxd)
     .                    ,nxdp,gv%gparams%xg(1:nxdp),dummy(1:nxdp)
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
            rh = gv%gparams%xg(i)
            jzf(i)=(xh(i)*bth0(i)-xh(i-1)*bth0(i-1))/(rh*dr(i))
            jtf(i)=-(bzz0(i)-bzz0(i-1))/dr(i)
          enddo

          !Pressure solve

          !Find heat source on integer mesh
          do i=1,nxd
            rh = gv%gparams%xg(i)
            p0(i)=rh*etah(i)*(jtf(i)*jtf(i)+jzf(i)*jzf(i))
! Daniele, 12 Nov. 2014
     .           +rh*ext_heat_src_prof(rh)/eta
          enddo
          p0 = etah0*p0

! Daniele, 1/8/2014. BCs with finite p(a)
          rh = gv%gparams%xg(nxd)
          p0(nxd) = p0(nxd) + pa/dr(nxd)*
     .         (-gamma/(gamma-1.)*xh(nxd)*vr0(nxd)
     .          +vr0(nxd)*rh
     .          +2.*chih(nxd)/eta*xh(nxd)/dr(nxd))

          !Build the tridiagonal matrix to be inverted
          do i=1,nxd
            rh  = gv%gparams%xg(i)

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
          rh  = gv%gparams%xg(1)
          d(1)=d(1)+vr0(1)*rh/(2.*dr(1))

          !BSc at i=nxd
          rh  = gv%gparams%xg(nxd)
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
            rh = gv%gparams%xg(i)
            etah(i) = res_rfx(rh)
          enddo
        endif

        !Check force balance
        do i=1,nxd
          rh  = gv%gparams%xg(i)

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
          rh  = gv%gparams%xg(i)

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
          rh  = gv%gparams%xg(i)
          rhm = gv%gparams%xg(i-1)
          rhp = gv%gparams%xg(i+1)

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
cc        eq_params(1) = 1d0/gv%gparams%params(3)       !Aspect ratio
cc        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
cc        eq_params(3) = q0h(nxd)/gv%gparams%params(3)  !q(r=a)

        !Find helical flux
        psi0(0) = 0d0
        do i=1,nxd+1
          psi0(i) = psi0(i-1)+dr(i-1)*(bth0(i-1)-kk*xh(i-1)*bzz0(i-1))
        enddo
        psi0 = psi0 - 0.5*(psi0(nxd)+psi0(nxdp)) !Homogeneous Dirichlet BC

        !Transfer B to integer mesh
        bth0 = btf
        bzz0 = bzf

#if defined(RFX)
        if (my_rank == 0.and.(.not.test)) then
          write (*,*) 'r,bth0,bzz0,q0,vr0,p0'
          do i=0,nxdp
            rh  = gv%gparams%xg(i)
            write(*,10) rh,bth0(i),bzz0(i)
     .             ,(rh*bzz0(i))/(RR*bth0(i)),vr0(i),p0(i)
          enddo
        endif
 10     format (4f7.3,2e12.3)
#endif

 101    format (5f7.3,2e12.3)

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

        use fixed_point_accelerator

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

        type (fpa_state) :: fpav

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

        call allocateMGArray(gv%gparams,1,gpsi0)
        call allocateMGArray(gv%gparams,1,gbz0)

        gpsi0%grid(igrid)%array(:,:,:,1) = dpsi
        gbz0 %grid(igrid)%array(:,:,:,1) = dIII

        bcnd = bcond
        where (bcnd == DEF) bcnd = EXT

        call restrictMGArray(IRHO,1,gpsi0,bcnd,igrid,order_bc)
        call restrictMGArray(IRHO,1,gbz0 ,bcnd,igrid,order_bc)

c     Gather solution

        x = reshape(psi(1:nx,1:ny,1:nz),shape(x))

c     Iteration

        !Create fixed-point accelerator
        call fpa_create(fpav,x,maxv=5)

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

          !Apply fixed-point accelerator
          call fpa_correction(fpav,r,dp=dot2)
cc          alpha = 0d0

          !Update solution
          x = x - (1d0-alpha)*r

          !Check convergence
          if (chk_conv(r)) exit

        enddo

        !Destroy fixed-point accelerator
        call fpa_destroy(fpav)

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
          call getCurvilinearCoordinates(gv%gparams,i,j,k,igrid,igrid
     .         ,igrid,ig,jg,kg,x1,y1,z1)
          III(i,:,:) = bb(i,:,:,3)  + kk*x1*bb(i,:,:,2)
        enddo

c     Initialize Psi from 1D global psi

        do i=0,nx+1
          iglobal = i + gv%gparams%ilo(igrid) - 1
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
            call getCurvilinearCoordinates(gv%gparams,i,j,k,igrid,igrid
     .           ,igrid,ig,jg,kg,x1,y1,z1)

cc            dpsi =-eps/mm*sin(pi*x1/gv%gparams%xg(nxd+1))*sin(y1*mm)  !Produces finite Br perturbation at
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
        call setBC(gv%gparams,IRHO,nx,ny,nz,psi
     .            ,gpsi0%grid(igrid)%array(:,:,:,1)
     .            ,bcnd,igrid,igrid,igrid,iorder=order_bc)

c     Phase out psi BC

        if (phaseout_psibc) then
          aa = aa*0.5
          if (aa > (xmin-xmax)/nxd) then

            if (my_rank ==0) write (*,*) 'Phasing out psi BC...'

            do i=1,nx
              call getCurvilinearCoordinates(gv%gparams,i,1,1,igrid
     .             ,igrid,igrid,ig,jg,kg,x1,y1,z1)
           
              psi(i,:,:) = psi(i,:,:)*tanh((1-x1)/eps+aa)
            enddo

            do igr=1,gv%gparams%ngrid
              gpsi0%grid(igr)%array = aa*gpsi0%grid(igr)%array
            enddo
          else
            do igr=1,gv%gparams%ngrid
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
              rh = gv%gparams%xg(i)

              !Cartesian components
              th = gv%gparams%yg(j)

              bcar(i,j,k,1) = bbg(i,j,k,1)/rh*cos(th)
     .                      -(bbg(i,j,k,2)+kk*bbg(i,j,k,3))*sin(th)
              bcar(i,j,k,2) = bbg(i,j,k,1)/rh*sin(th)
     .                      +(bbg(i,j,k,2)+kk*bbg(i,j,k,3))*cos(th)
              bcar(i,j,k,3) = bbg(i,j,k,3)/rh

              !Map (cylindrical)
              jac (i,j,k)  = rh
              xmap(i,j,k,1)= rh*cos(th)
              xmap(i,j,k,2)= rh*sin(th)
              xmap(i,j,k,3)= gv%gparams%zg(k)
            enddo
          enddo
        enddo

c     Orbit setup (input Az directly)

cc        bxg = 0d0
cc        byg = 0d0
cc        bzg = -psig
        call orbit_setup(gv%gparams
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
          etag(i,:,:) = res_rfx(gv%gparams%xg(i))
        enddo

        call splineFld(etag)

c     Orbit integration (local processor)

        ff = 0d0

        do k=1,nz
          do j=1,ny
            do i=1,nx
              call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)
              x0 = gv%gparams%xx(ig)
              y0 = gv%gparams%yy(jg)
              z0 = gv%gparams%zz(kg)

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

                ff(i,j,k) = quad_int(gv%gparams%xx(ig-4)
     .                              ,gv%gparams%xx(ig-3)
     .                              ,gv%gparams%xx(ig-2)
     .                              ,gv%gparams%xx(ig-1)
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
        call setBC(gv%gparams,IRHO,nx,ny,nz,ff,zeros,bcnd
     .            ,igrid,igrid,igrid,iorder=order_bc)
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

              call getCurvilinearCoordinates(gv%gparams,i,j,k,igrid
     .             ,igrid,igrid,ig,jg,kg,x1,y1,z1)

              rhs(ii)=-III(i,j,k)/(1+(kk*x1)**2)
     .                *(2*kk/(1+(kk*x1)**2)-jpar_B(i,j,k))

              if (vol_wgt) rhs(ii)=rhs(ii)
     .                       *gv%gparams%gmetric%grid(igrid)%dvol(i,j,k)
            enddo
          enddo
        enddo

c       Recalibrate psi BC

c       Solve Helical Grad-Shafranov equation (find psi at cell centers)

        bcnd = bcond
        where (bcnd == DEF) bcnd = EQU

        psiv = reshape(psi(1:nx,1:ny,1:nz),shape(psiv))

        guess = 1
        call cSolver(gv%gparams,1,nn,rhs,psiv,bcnd
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
        call setBC(gv%gparams,IRHO,nx,ny,nz,psi_new
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
              rhs(ii)=-lap(gv%gparams,i,j,k,nx,ny,nz,igrid,igrid,igrid
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
        call cSolver(gv%gparams,1,nn,rhs,IIv,bcnd
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

        call setBC(gv%gparams,IRHO,nx,ny,nz,III
     .            ,gbz0%grid(igrid)%array(:,:,:,1)
     .            ,bcnd,igrid,igrid,igrid,iorder=order_bc)

c       Renormalize III so that Bz(r=0)=1d0

        if (isBdry(gv%gparams,1,igrid,1)) then
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
              call getCurvilinearCoordinates(gv%gparams,i,j,k,igrid
     .             ,igrid,igrid,ig,jg,kg,x1,y1,z1)

              gr = grad_ijk(gv%gparams,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                     ,psi)

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
        call setMGBC(gv%gparams,0,3,nx,ny,nz,igrid,bb,bcnd
     .              ,icomp=(/IBX/),is_cnv=.true.,is_vec=.true.
     $              ,iorder=order_bc)

      end function compute_B

c     chk_conv
c     ##################################################################
      function chk_conv(dx)

        implicit none

c     ------------------------------------------------------------------
c     Checks iteration convergence
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

        err = rms(dx)

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
              call getCurvilinearCoordinates(gv%gparams,i,j,k,igrid
     .             ,igrid,igrid,ig,jg,kg,x1,y1,z1)

              j3=(gv%gparams%xx(ig+1)*(bb(i+1,j,k,2)+kk*bb(i+1,j,k,3))
     .           -gv%gparams%xx(ig-1)*(bb(i-1,j,k,2)+kk*bb(i-1,j,k,3)))
     .           /2./gv%gparams%dxh(ig)  ! d(B_2)/dx1
     .          -(bb(i,j+1,k,1)-bb(i,j-1,k,1))/x1/2./gv%gparams%dyh(jg)
                                          !-d(B_1)/dx2

              gs_err(i,j,k) = j3-bb(i,j,k,3)*jpar_B(i,j,k)
            enddo
          enddo
        enddo

        err = integral(gv%gparams,igrid,gs_err**2)
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
              call getCurvilinearCoordinates(gv%gparams,i,j,k,igrid
     .             ,igrid,igrid,ig,jg,kg,x1,y1,z1)

              !Cyl B
              br = bb(i,j,k,1)/x1
              bth= bb(i,j,k,2)+kk*bb(i,j,k,3)
              bz = bb(i,j,k,3)/x1

              bnorm2 = br**2+bth**2+bz**2

              !Cnv J
              j1= (III(i,j+1,k)-III(i,j-1,k))/2./gv%gparams%dyh(jg) ! d(B3)/dx2
              j2=-(III(i+1,j,k)-III(i-1,j,k))/2./gv%gparams%dxh(ig) !-d(B3)/dx1
              j3=(gv%gparams%xx(ig+1)*(bb(i+1,j,k,2)+kk*bb(i+1,j,k,3))
     .           -gv%gparams%xx(ig-1)*(bb(i-1,j,k,2)+kk*bb(i-1,j,k,3)))
     .           /2./gv%gparams%dxh(ig)  ! d(B_2)/dx1
     .          -(bb(i,j+1,k,1)-bb(i,j-1,k,1))/x1/2./gv%gparams%dyh(jg)
     .                                    !-d(B_1)/dx2

              !Cov B
              b1 = br
              b2 = x1*bth
              b3 = III(i,j,k)

              ohm_err(i,j,k)=(j1*b1+j2*b2+j3*b3)/x1-jpar_B(i,j,k)*bnorm2
            enddo
          enddo
        enddo

        err = integral(gv%gparams,igrid,ohm_err**2)
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
          call getCurvilinearCoordinates(gv%gparams,i,j,k,igrid,igrid
     .         ,igrid,ig,jg,kg,x1,y1,z1)
          dummy(i,:,:) = bb(i,:,:,3)/x1
        enddo

        bz_avg = integral(gv%gparams,igrid,dummy,average=.true.)

        !Average Bth @ boundary
        if (isBdry(gv%gparams,nx,igrid,2)) then
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
              call getCurvilinearCoordinates(gv%gparams,i,j,k,igrid
     .             ,igrid,igrid,ig,jg,kg,x1,y1,z1)

              !Cyl B
              br = bb(i,j,k,1)/x1
              bth= bb(i,j,k,2)+kk*bb(i,j,k,3)
              bz = bb(i,j,k,3)/x1

              bnorm2 = br**2+bth**2+bz**2

              !Cov B
              b1 = br
              b2 = gv%gparams%xx(ig)*bth

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

        use mg_xfer, ONLY: fillArray

        use imposeBCinterface

        use ts_setup, ONLY: vol_wgt

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
     .            ,Eoeta,b_tor,R0

c     Local variables

        integer :: i,j,k,ig,jg,kg,ii,it,bcnd(6,1),igrid,guess
     .            ,iglobal,jglobal,pos(2),nn
     .            ,np_RAD,my_rank_RAD,np_Z,my_rank_Z

        real(8) :: bz_avg,bth_avg,gserr,oherr,oherr_n,I0,bnorm2
     .            ,psir,psith,x1,y1,z1,IIh,spsi,psi0,psia
     .            ,cc,p00,etah,lpsi0,dsdpsi,dpdpsi,spsip,spsim,cnv(3)
     .            ,dpdr,dpdth,f_1,f_2,f_3,j1,j2,j3,jac,e1,e2,e3,II2
     .            ,psip,psim,q0,qa,tor_flx,bz2_avg,p_avg,b_pol,gpsi(3)
     .            ,cov(3),e4,e5,e6,bcnv(3),gp(3),gp_cnv(3),dp_dpsi,gpsi2
     .            ,bcov(3),th,zz,dh,dt,rh

        real(8) :: xh(0:nxd),drh(0:nxd),dr(0:nxd)

        real(8) :: iR2  (0:nx+1,0:ny+1,0:nz+1)
     .            ,zeros(0:nx+1,0:ny+1,0:nz+1)
     .            ,qprof(0:nx+1,0:ny+1)
     .            ,dummy(0:nx+1,0:ny+1)
     .            ,III  (0:nx+1,0:ny+1)
     .            ,RR   (0:nx+1,0:ny+1)
     .            ,psi(0:nx+1,0:ny+1,0:nz+1,1)

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

        !Find half-mesh positions and spacings
        xh (0) = 0.5*gv%gparams%xg(1)
        drh(0) =     gv%gparams%xg(1)
        do i=1,nxd
          xh (i) = 0.5*(gv%gparams%xg(i+1)+gv%gparams%xg(i))
          drh(i) =      gv%gparams%xg(i+1)-gv%gparams%xg(i)
        enddo

        !Integer mesh spacings
        do i=1,nxd
          dr(i) = xh(i)-xh(i-1)
        enddo
        dr(0) = dr(1)

        !Initial guess for psi
        psi = 0d0

        !Initalize toroidal radius (R=R0+r*cos(theta)) and 1/R^2
        do j=0,ny+1
          do i=0,nx+1
            call getCartesianCoordinates(gv%gparams,i,j,1
     .           ,igrid,igrid,igrid,ig,jg,kg,x1,y1,z1)
            RR(i,j) = sqrt(x1**2+y1**2)
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
          oherr_n = oherr

          do k=1,nz
            do j=1,ny
              do i=1,nx
              ii  = i + nx*(j-1) + nx*ny*(k-1)

              call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)

              call find_polar_coords(gv%gparams,igrid,i,j,k,rh,th,zz
     .                              ,coords)
c$$$              rh = gv%gparams%xx (ig)

              dh = gv%gparams%dxh(ig)
              dt = gv%gparams%dyh(jg)

              jac = gv%gparams%gmetric%grid(igrid)%jac(i,j,k)
              !Compute (R*|B|)^2
              IIh = 0.5*(III(i,j)+III(i-1,j))

              gpsi = grad(gv%gparams,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                   ,psi(:,:,:,1))
              gpsi(3) = 0d0

              gpsi2 = vectorNorm_ijk(gv%gparams,i,j,k,igrid,gpsi,.true.)

              bcnv(1) = -gpsi(2)
              bcnv(2) = gpsi(1)
              bcnv(3) = jac*IIh*iR2(i,j,k)
              bnorm2=vectorNorm_ijk(gv%gparams,i,j,k,igrid,bcnv,.false.)
              bnorm2=bnorm2*RR(i,j)*RR(i,j)  !This is R^2*B^2

              !Compute grad(p)
              spsip = sqrt(1d0-psi(i+1,j,k,1)/psi0)
              spsim = sqrt(1d0-psi(i-1,j,k,1)/psi0)
              gp(1) = 0.5*p00*(pprof(spsip,0)-pprof(spsim,0))/dh

              spsip = sqrt(1d0-psi(i,j+1,k,1)/psi0)
              spsim = sqrt(1d0-psi(i,j-1,k,1)/psi0)
              gp(2) = 0.5*p00*(pprof(spsip,0)-pprof(spsim,0))/dt

              gp(3) = 0d0

              !Compute dp/dpsi=grad_psi.grad_p/|grad_psi|^2
              gp_cnv = XformToCnv_ijk(gv%gparams,i,j,k,igrid,gp)

              dp_dpsi = scalarProduct_ijk(gv%gparams,i,j,k,igrid
     .                                   ,gpsi,gp_cnv)
c$$$              if (gpsi2 > 0d0) then
c$$$                dp_dpsi = dp_dpsi/gpsi2
c$$$              else
c$$$                dp_dpsi = 0d0
c$$$              endif

              !Compute GS RHS
              etah =res_rfx(rh)
c$$$              II2 = IIh*IIh/bnorm2
c$$$              psi_rhs(ii)= (Eoeta/etah*iR2(i,j,k)*II2-dp_dpsi*(1d0-II2))
              II2 = IIh*IIh
              psi_rhs(ii)= (Eoeta/etah*iR2(i,j,k)*II2-dp_dpsi)/bnorm2

              !Multiply by R to ensure proper limit R->inf
              errv(ii) =(lap(gv%gparams,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                      ,psi(:,:,:,1),dff=iR2,vol=.false.)
     .                  -psi_rhs(ii))*RR(i,j)

              if (vol_wgt) psi_rhs(ii)=psi_rhs(ii)
     .                       *gv%gparams%gmetric%grid(igrid)%dvol(i,j,k)

              !Check convergence of Ohm's law
              psir  = gpsi(1)
              dpdr  = gp(1)
              oherr = oherr
     .               +(Eoeta/etah*IIh*iR2(i,j,k)*psir
     .                +IIh*dpdr
     .                +bnorm2*iR2(i,j,k)  !R^2 factor needed to normalize bnorm2
     .                *(III(i,j)-III(i-1,j))/dh)**2
              enddo
            enddo
          enddo

          !Calculate L2 norm of error
          gserr = rms(errv)!sqrt(dot2(errv,errv)/nxd/nyd/nzd)

#if defined(petsc)
          e1 = oherr
          call MPI_Allreduce(e1,oherr,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,MPI_COMM_WORLD,mpierr)
#endif
          oherr = sqrt(oherr/nxd/nyd/nzd)

          !Solve GS Poisson eq (find psi at cell centers)
          !Define BCs
          bcnd(:,1) = bcond
          where (bcnd == DEF) bcnd = DIR

          guess = 1  !To reuse previous iteration

          call cSolver(gv%gparams,1,nx*ny*nz,psi_rhs,psiv,bcnd
     .              ,1,ilevel,guess,del_star,vol_wgt
     .              ,tol           = 1d-1
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

          !Fill array (w BCs)
          call fillArray(gv%gparams,igrid,1,psiv,psi,bcnd
     .                  ,iorder=order_bc,icomp=(/IRHO/))

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
              call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)

              call find_polar_coords(gv%gparams,igrid,i,j,k,rh,th,zz
     .                              ,coords)
c$$$              rh = gv%gparams%xx(ig)

              dh = gv%gparams%dxh(ig)
              dt = gv%gparams%dyh(jg)

              jac = gv%gparams%gmetric%grid(igrid)%jac(i,j,k)

              !Compute |B|^2
              gpsi = grad(gv%gparams,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                   ,psi(:,:,:,1))
              gpsi(3) = 0d0

              psir   = gpsi(1)

              IIh = 0.5*(III(i,j)+III(i-1,j))

              bcnv(1) = -gpsi(2)
              bcnv(2) = gpsi(1)
              bcnv(3) = jac*IIh*iR2(i,j,k)
              bnorm2=vectorNorm_ijk(gv%gparams,i,j,k,igrid,bcnv,.false.)

              !Compute grad(p)
              psip  = 0.5*(psi(i+1,j,k,1)
     .                    +psi(i  ,j,k,1))
              psim  = 0.5*(psi(i-1,j,k,1)
     .                    +psi(i  ,j,k,1))
              spsip = sqrt(1d0-psip/psi0)
              spsim = sqrt(1d0-psim/psi0)
              dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/dh

              etah  = res_rfx(rh)

              dummy(i,j) = dummy(i-1,j)
     .                    -(Eoeta*psir/etah*iR2(i,j,k)+dpdr)/bnorm2
     .                    *gv%gparams%dx(ig)
            enddo
          enddo

#if defined(petsc)
          !Broadcast poloidal integral along Z
          call initMPI(gv%gparams%MPI_COMM_Z,np_Z,my_rank_Z)

          nn = size(dummy)
          call MPI_Bcast(dummy,nn,MPI_DOUBLE_PRECISION
     .                  ,0,gv%gparams%MPI_COMM_Z,mpierr)

          !Communicate radial integration constants
          if (npx > 1) then
            allocate(mpivec(npx*ny))
            call MPI_Allgather(dummy(nx,1:ny)
     .                             ,ny,MPI_DOUBLE_PRECISION
     .                      ,mpivec,ny,MPI_DOUBLE_PRECISION
     .                      ,gv%gparams%MPI_COMM_RAD,mpierr)

            !Sum integration constants
            do j=2,npx
              mpivec(1+ny*(j-1):ny*j) = mpivec(1+ny*(j-1):ny* j   )
     .                                 +mpivec(1+ny*(j-2):ny*(j-1))
            enddo

            !Add integration constants
            call initMPI(gv%gparams%MPI_COMM_RAD,np_rad,my_rank_rad)

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
          if (gserr<0.001*dh**2 .or. abs((1d0-oherr_n/oherr))<1d-3) exit

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
            call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            call find_polar_coords(gv%gparams,igrid,i,j,k,rh,th,zz
     .                            ,coords)
c$$$            rh = gv%gparams%xx(ig)
            dh = gv%gparams%dxh(ig)
            dt = gv%gparams%dyh(jg)

            jac = gv%gparams%gmetric%grid(igrid)%jac(i,j,k)

            !Find contravariant magnetic field components in integer mesh
            IIh    = 0.5*(III(i,j)+III(i-1,j))

            gpsi = grad(gv%gparams,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                 ,psi(:,:,:,1))
            gpsi(3) = 0d0

            gpsi2 = vectorNorm_ijk(gv%gparams,i,j,k,igrid,gpsi,.true.)

            bcnv(1) = -gpsi(2)
            bcnv(2) = gpsi(1)
            bcnv(3) = jac*IIh*iR2(i,j,k)
            bnorm2=vectorNorm_ijk(gv%gparams,i,j,k,igrid,bcnv,.false.)

            b1(i,j,:) = bcnv(1)
            b2(i,j,:) = bcnv(2)
            b3(i,j,:) = bcnv(3)

            bcov = XformToCov_ijk(gv%gparams,i,j,k,igrid,bcnv)

            !Find contravariant velocity components
            etah =res_rfx(rh)

            spsip = sqrt(1d0-psi(i+1,j,k,1)/psi0)
            spsim = sqrt(1d0-psi(i-1,j,k,1)/psi0)
            gp(1) = 0.5*p00*(pprof(spsip,0)-pprof(spsim,0))/dh

            spsip = sqrt(1d0-psi(i,j+1,k,1)/psi0)
            spsim = sqrt(1d0-psi(i,j-1,k,1)/psi0)
            gp(2) = 0.5*p00*(pprof(spsip,0)-pprof(spsim,0))/dt

            gp(3) = 0d0

            gp_cnv = XformToCnv_ijk(gv%gparams,i,j,k,igrid,gp)

            v1(i,j,:) =-(Eoeta*bcov(2)+etah*gp_cnv(1))/bnorm2
            v2(i,j,:) = (Eoeta*bcov(1)-etah*gp_cnv(2))/bnorm2

            !Find pressure
            spsi = sqrt(1d0-psi(i,j,k,1)/psi0)
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
            call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)

            call find_polar_coords(gv%gparams,igrid,i,j,k,rh,th,zz
     .                            ,coords)

c$$$            rh = gv%gparams%xx(ig)
            dh = gv%gparams%dxh(ig)
            dt = gv%gparams%dyh(jg)

            jac = gv%gparams%gmetric%grid(igrid)%jac(i,j,k)

            !Compute |B|^2
            IIh    = 0.5*(III(i,j)+III(i-1,j))

            gpsi = grad(gv%gparams,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                 ,psi(:,:,:,1))
            gpsi(3) = 0d0

            gpsi2 = vectorNorm_ijk(gv%gparams,i,j,k,igrid,gpsi,.true.)

            bcnv(1) = -gpsi(2)
            bcnv(2) = gpsi(1)
            bcnv(3) = jac*IIh*iR2(i,j,k)
            bnorm2=vectorNorm_ijk(gv%gparams,i,j,k,igrid,bcnv,.false.)

            !Find pressure derivatives
            spsip = sqrt(1d0-psi(i+1,j,k,1)/psi0)
            spsim = sqrt(1d0-psi(i-1,j,k,1)/psi0)
            dpdr  = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dh

            spsip = sqrt(1d0-psi(i,j+1,k,1)/psi0)
            spsim = sqrt(1d0-psi(i,j-1,k,1)/psi0)
            dpdth = p00*(pprof(spsip,0)-pprof(spsim,0))/2./dt

            !Find contravariant current components
            j1 = (III(i,j+1)+III(i-1,j+1)
     .           -III(i,j-1)-III(i-1,j-1))/4./dt
            j2 =-(III(i,j)  -III(i-1,j))/dh
            j3 = jac*lap(gv%gparams,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                  ,psi(:,:,:,1),dff=iR2,vol=.false.)

            !Find covariant EOM force components
            f_1 = (j2*b3(i,j,k)-j3*b2(i,j,k))/jac - dpdr
            f_2 = (j3*b1(i,j,k)-j1*b3(i,j,k))/jac - dpdth
            f_3 = (j1*b2(i,j,k)-j2*b1(i,j,k))/jac

            !Accumulate errors
            e1 = e1 + f_1**2*gv%gparams%gmetric%grid(igrid)%dvol(i,j,k)
            e2 = e2 + f_2**2*gv%gparams%gmetric%grid(igrid)%dvol(i,j,k)
            e3 = e3 + f_3**2*gv%gparams%gmetric%grid(igrid)%dvol(i,j,k)

            !Find covariant Ohm's law components
            etah =res_rfx(rh)

            cov = matmul(gv%gparams%gmetric%grid(1)%gsub(i,j,k,:,:)
     .                  ,(/ j1,j2,j3/))

            f_1 =      (v2(i,j,k)*b3(i,j,k)                    )/jac
            f_2 =      (                   -v1(i,j,k)*b3(i,j,k))/jac
            f_3 =Eoeta+(v1(i,j,k)*b2(i,j,k)-v2(i,j,k)*b1(i,j,k))/jac

            f_1 =f_1-etah*cov(1)
            f_2 =f_2-etah*cov(2)
            f_3 =f_3-etah*cov(3)

            !Accumulate errors
            e4 = e4 + f_1**2*gv%gparams%gmetric%grid(igrid)%dvol(i,j,k)
            e5 = e5 + f_2**2*gv%gparams%gmetric%grid(igrid)%dvol(i,j,k)
            e6 = e6 + f_3**2*gv%gparams%gmetric%grid(igrid)%dvol(i,j,k)
          enddo
        enddo

#if defined (petsc)
        lpsi0 = e1
        call MPI_Allreduce(lpsi0,e1,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)

        lpsi0 = e2
        call MPI_Allreduce(lpsi0,e2,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)

        lpsi0 = e3
        call MPI_Allreduce(lpsi0,e3,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)

        lpsi0 = e4
        call MPI_Allreduce(lpsi0,e4,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)

        lpsi0 = e5
        call MPI_Allreduce(lpsi0,e5,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)

        lpsi0 = e6
        call MPI_Allreduce(lpsi0,e6,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)
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
          write (*,*)
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
            call getMGmap(gv%gparams,i,j,k,igrid,igrid,igrid,ig,jg,kg)

c$$$            rh = gv%gparams%xx(ig)
c$$$            dh = gv%gparams%dxh(ig)
c$$$            dt = gv%gparams%dyh(jg)

            !Q-profile
            qprof(i,j) = b3(i,j,k)/b2(i,j,k)

            !Toroidal flux
            tor_flx = tor_flx + b3(i,j,k)
     .               *gv%gparams%gmetric%grid(igrid)%dlvol(ig,jg,kg) !dh*dt

            !Bz^2 average
            bz2_avg = bz2_avg+(b3(i,j,k))**2
     .           *gv%gparams%gmetric%grid(igrid)%dvol(ig,jg,kg)
     .           /gv%gparams%gmetric%grid(igrid)%jac (ig,jg,kg)!/rh*dh*dt

            !Pressure average
            p_avg = p_avg + pp(i,j,k)
     .            *gv%gparams%gmetric%grid(igrid)%dvol(ig,jg,kg) !*rh*dh*dt

            !Toroidal current
            iglobal = i + gv%gparams%ilo(igrid) - 1
            if (iglobal == nxd) then
              I0 = I0 + (psi(i+1,j,k,1)-psi(i,j,k,1))/dh
     .                  *2./(RR(i,j)+RR(i+1,j))*dt
            endif

          enddo
        enddo

#if defined (petsc)
        lpsi0 = tor_flx
        call MPI_Allreduce(lpsi0,tor_flx,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)

        lpsi0 = bz2_avg
        call MPI_Allreduce(lpsi0,bz2_avg,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)

        lpsi0 = p_avg
        call MPI_Allreduce(lpsi0,p_avg,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)

        lpsi0 = I0
        call MPI_Allreduce(lpsi0,I0,1,MPI_DOUBLE_PRECISION
     .                    ,MPI_SUM,gv%gparams%MPI_COMM_XY,mpierr)
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

c     Transform to vector potential
        
#if defined(vec_pot)

cc        !Find B^3 at faces (III is defined at radial faces)
cccc        III(0,:) = 0d0
cc        do i=0,nx
cc          iglobal = i + gv%gparams%ilo(igrid) - 1
cc          
cc          III(i,:)= 2.*III(i,:)/(RR(i+1,:)+RR(i,:))*xh(iglobal)
cccc          write (*,*) xh(iglobal),III(i,1),(RR(i+1,1)+RR(i,1))/2.
cc        enddo
cc
cc        !Find A_2 at cell centers by radial integral of b^3 at cell faces
cc        dummy(0,:) = 0d0
cc        do i=1,nx+1
cc          iglobal = i + gv%gparams%ilo(igrid) - 1
cc
cc          dummy(i,:) = dummy(i-1,:) +  III(i-1,:)*dr(iglobal-1)
cc        enddo

        !Find A_2 at cell centers by staggered radial integral of b^3 at cell centers
        dummy = 0d0
        if (isSP2(gv%gparams,1,igrid)) then
          write (*,*) my_rank
          !Estimate B^3 at r=0
          iglobal = gv%gparams%ilo(igrid) - 1
          dummy(1,:)= 2.*III(0,:)/(RR(1,:)+RR(0,:))*xh(iglobal)
          !Find A_2 at first cell
          dummy(1,:) = dummy(1,:)*drh(iglobal)
          !Extrapolate to ghost cell
          dummy(0,:) = dummy(1,:)
        endif
        
        !Integrate staggered
        do i=1,nx,2
          call getMGmap(gv%gparams,i,1,1,igrid,igrid,igrid,ig,jg,kg)

          dummy(i+1,:) = dummy(i-1,:) + b3(i,:,1)*2*gv%gparams%dxh(ig)
        enddo

#if defined(petsc)
        !Communicate radial integration constants
        if (npx > 1) then
          nn = ny+2
          allocate(mpivec(npx*nn))
          call MPI_Allgather(dummy(nx,0:ny+1)
     .                           ,nn,MPI_DOUBLE_PRECISION
     .                    ,mpivec,nn,MPI_DOUBLE_PRECISION
     .                    ,gv%gparams%MPI_COMM_RAD,mpierr)

          !Sum integration constants
          do j=2,npx
            mpivec(1+nn*(j-1):nn*j) = mpivec(1+nn*(j-1):nn* j   )
     .                               +mpivec(1+nn*(j-2):nn*(j-1))
          enddo

          !Add integration constants
          call initMPI(gv%gparams%MPI_COMM_RAD,np_rad,my_rank_rad)

          j = my_rank_rad
          if (j > 0) then
            do i=0,nx+1,2
              dummy(i,:) = dummy(i,:) + mpivec(1+nn*(j-1):nn*j)
            enddo
          endif

          deallocate(mpivec)
        endif
#endif

        do i=2,nx,2
          call getMGmap(gv%gparams,i,1,1,igrid,igrid,igrid,ig,jg,kg)

          dummy(i+1,:) = dummy(i-1,:) + b3(i,:,1)*2*gv%gparams%dxh(ig)
          write (*,*) my_rank,i+1
        enddo

#if defined(petsc)
        !Communicate radial integration constants
        if (npx > 1) then
          nn = ny+2
          allocate(mpivec(npx*nn))
          call MPI_Allgather(dummy(nx+1,0:ny+1)
     .                           ,nn,MPI_DOUBLE_PRECISION
     .                    ,mpivec,nn,MPI_DOUBLE_PRECISION
     .                    ,gv%gparams%MPI_COMM_RAD,mpierr)

          !Sum integration constants
          do j=2,npx
            mpivec(1+nn*(j-1):nn*j) = mpivec(1+nn*(j-1):nn* j   )
     .                               +mpivec(1+nn*(j-2):nn*(j-1))
          enddo

          !Add integration constants
          call initMPI(gv%gparams%MPI_COMM_RAD,np_rad,my_rank_rad)

          j = my_rank_rad
          if (j > 0) then
            do i=1,nx+1,2
              dummy(i,:) = dummy(i,:) + mpivec(1+nn*(j-1):nn*j)
            enddo
          endif

          deallocate(mpivec)
        endif
#endif

        !Return vector potential components
        b1 = 0d0
        b3 = -psi(:,:,:,1)

        do k=0,nz+1
          b2(:,:,k) = dummy
        enddo

#endif
        
        !Save some equilibrium properties (to calculate bootstrap current)
cc        eq_params(1) = 1d0/gv%gparams%params(3)       !Aspect ratio
cc        eq_params(2) = 2*p0h(0)/(bth0(nxd)**2)         !Poloidal beta
cc        eq_params(3) = q0h(nxd)/gv%gparams%params(3)  !q(r=a)

        !Dump equilibrium for eigen code
cc        if (my_rank == 0) then
cc          open(unit=123,file='ntm2d.txt',status='unknown')
cc          write (123,*) nxdp+1
cc          write (123,*) gv%gparams%xg(0:nxdp)
cc          write (123,*) bth0
cc          write (123,*) bzz0
cc          write (123,*) 0.5*p0
cc          write (123,*) gv%gparams%params(3),eq_params(2),eq_params(3)
cc          close(123)
cccc          stop
cc        endif

c diag ***** 
cc      open(unit=110,file='debug.bin',form='unformatted'
cc     .    ,status='replace')
cccc      call contour(psi(:,:,1),nx+2,ny+2,0d0,1d0,0d0,1d0,0,110)
cc      call contour(dummy,nx+2,ny+2,0d0,1d0,0d0,1d0,0,110)
cc      call contour(III,nx+2,ny+2,0d0,1d0,0d0,1d0,0,110)
cc      call contour(b1 (1:nx,1:ny)  ,nx,ny,0d0,1d0,0d0,1d0,1,110)
cc      call contour(b2 (1:nx,1:ny)  ,nx,ny,0d0,1d0,0d0,1d0,1,110)
cc      call contour(b3 (1:nx,1:ny,1)  ,nx,ny,0d0,1d0,0d0,1d0,0,110)
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

        real(8) :: q0,qa,qstar,RR,aa,jz1,p

c     Begin program

        select case(equil)
        case('ntm2d')

          q0 = 1.02
          qa = 4.
          p  = 1.

          qprof = q0*(1+r**(2*p)*( (qa/q0)**p - 1 ))**(1./p)

        case('rfp83')

          RR   = gv%gparams%params(3)
          
          qprof = 0.6125*(1 - 1.8748*r**2 + 0.8323*r**4)/RR

        case('vela1')

          qstar = 1.5

          RR   = gv%gparams%params(3)
          
          aa = 1./(qstar*RR)
          
          qprof = r**2/(aa*RR)

        case('vela2')

          qprof = log(1+r)+exp(r)

        case('tbbr1')

          RR   = gv%gparams%params(3)
          
          q0   = gv%gparams%params(6)  !q0*R
          jz1  = gv%gparams%params(7)
          aa   = gv%gparams%params(8)
          p    = gv%gparams%params(9)

          qa = r**2/(0.5*jz1*r**2+aa**2*(1-0.5*jz1))

          qprof = q0/RR*(1d0+qa**p)**(1./p)
          
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

        real(8) :: qstar,aa,RR

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
            pprof = (1-s**2)**2
          else
            pprof =-4*s*(1-s**2)
          endif

        case('p3nsl')

          if (der == 0) then
            pprof = (1-s**2)*(1.4-s**2)! + 1d-3
          else
            pprof =-4*s*(1-s**2)
          endif

        case('rfp83','tbbr1')

          if (der == 0) then
            pprof = 1d0
          else
            pprof = 0d0
          endif

        case('vela1')

          qstar = 1.5

          RR   = gv%gparams%params(3)
          
          aa = 1./(qstar*RR)
          
          pprof = 0.5*(1-s*s)**2*(1+2*s**2-s**4+2./3.*aa**2*(5-2*s*s))

        case('vela2')

          pprof = 1.1 + 0.1*s**2 - 1.1*s**4

        case default

          pprof = 0d0

        end select

c     End program

      end function pprof

c fillVectorPotential
c####################################################################
      function fillVectorPotential(igr,a) result(find_a)

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

      logical :: find_a
      integer :: igr
      real(8) :: a(ilom:ihip,jlom:jhip,klom:khip,3)

c Local variables

      integer :: ig,jg,kg,i,j,k,bcnd(6,3),nx,ny,nz
      real(8) :: x1,y1,z1,r1,eps,car(3),hlx,RR,aa,bz0
      real(8) :: a_cnv(ilom:ihip,jlom:jhip,klom:khip,3)
      logical :: cartsn

c Begin program

      find_a = .false.

      nx = gv%gparams%nxv(igr)
      ny = gv%gparams%nyv(igr)
      nz = gv%gparams%nzv(igr)

      select case (trim(equil))

      case ('khcar')

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(gv%gparams,i,j,k,igr,igr,igr
     .                                    ,ig,jg,kg,x1,y1,z1)
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

              call getCartesianCoordinates(gv%gparams,i,j,k,igr,igr,igr
     .                                    ,ig,jg,kg,x1,y1,z1)
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

              call getCartesianCoordinates(gv%gparams,i,j,k,igr,igr,igr
     .                                    ,ig,jg,kg,x1,y1,z1)
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

              call getCartesianCoordinates(gv%gparams,i,j,k,igr,igr,igr
     .                                    ,ig,jg,kg,x1,y1,z1)
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

        eps = eq_params(2)

        if (bcond(3) == FSYM .or. bcond(4) == FSYM) then
          dlambda = (ymax-ymin)/2./pi !(2*Ly/4/pi)
        else
          dlambda = (ymax-ymin)/4./pi !(Ly/4/pi)
        endif

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(gv%gparams,i,j,k,igr,igr,igr
     .                                    ,ig,jg,kg,x1,y1,z1)
              a(i,j,k,1) = 0d0
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = -dlambda*dlog(    cosh(x1/dlambda)
     $                                  +eps*cos (y1/dlambda)) 
            enddo
          enddo
        enddo

      case ('ic-hk')

        eps = eq_params(2)

        if (bcond(3) == FSYM .or. bcond(4) == FSYM) then
          dlambda = (ymax-ymin)/2./pi !(2*Ly/4/pi)
        else
          dlambda = (ymax-ymin)/4./pi !(Ly/4/pi)
        endif

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip

              call getCartesianCoordinates(gv%gparams,i,j,k,igr,igr,igr
     .                                    ,ig,jg,kg,x1,y1,z1)
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

              call getCartesianCoordinates(gv%gparams,i,j,k,igr,igr,igr
     .                                    ,ig,jg,kg,x1,y1,z1)
              a(i,j,k,1) = dlambda
     .                    *dlog(cosh((y1-0.5d0*(ymax+ymin))/dlambda)) 
              a(i,j,k,2) = 0d0
              a(i,j,k,3) = 0d0
            enddo
          enddo
        enddo

      case ('bpnch')

        RR = gv%gparams%params(3) ! Major radius R=Lz/2pi.
        aa = gv%gparams%params(4) ! Radius^2 of q=1 surface.

        bz0 = RR/(1d0 + aa)/dlambda

        do k = klom,khip
          do j = jlom,jhip
            do i = ilom,ihip
              call find_polar_coords(gv%gparams,igr,i,j,k,x1,y1,z1
     .                              ,coords)
c$$$              call getCartesianCoordinates(gv%gparams,i,j,k,igr,igr
c$$$     .                                    ,igr,ig,jg,kg,x1,y1,z1)
c$$$              
c$$$              r1 = sqrt(x1**2+y1**2)
c$$$
c$$$              call getCurvilinearCoordinates(gv%gparams,i,j,k,igr,igr
c$$$     .                                      ,igr,ig,jg,kg,x1,y1,z1)
c$$$
c$$$              r1 = r1*sign(1d0,x1)

              a(i,j,k,1) = -0.5*Bz0*r1*sin(y1)
              a(i,j,k,2) =  0.5*Bz0*r1*cos(y1)
              a(i,j,k,3) = -0.5*dlambda*log(1+(r1/dlambda)**2)
          enddo
        enddo
      enddo

      case default

        call pstop("fillVectorPotential","Equilibrium def not found")

      end select

c Transform to Covariant

      call XformVector_mesh2(gv%gparams,igr,a,'car','cov')

c End program

      end function fillVectorPotential
