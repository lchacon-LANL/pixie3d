c module app_iosetup
c ######################################################################
      module app_iosetup

        use iosetup

        integer :: it_divcl=0

        logical :: dcon=.false.

        logical :: poincare=.false.,fourier=.false.,test=.false.

      end module app_iosetup

c module problem_def
c ######################################################################
      module problem_def

        use variables

        implicit none

#if !defined(vec_pot)
        integer, parameter :: IRHO=1,IVX=2,IVY=3,IVZ=4,IBX=5,IBY=6,IBZ=7
     .                 ,ITMP=8,IJX=9,IJY=10,IJZ=11,IAX=12,IAY=13,IAZ=14

        integer, parameter :: IJCNV    =1
     .                       ,IJCNV_0  =2
cc     .                       ,IJCNV_N  =23
     .                       ,IJCOV    =3
     .                       ,IJCOV_0  =4

        integer, parameter :: IBCNV    =5
cc     .                       ,IBCNV_N  =20
     .                       ,IBCNV_0  =6
     .                       ,IBCOV    =7
     .                       ,IBCOV_0  =8
cc     .                       ,IBHAT    =9

        integer, parameter :: IVCNV    =9
     .                       ,IVCNV_N  =10
     .                       ,IVCNV_0  =11
     .                       ,IVCOV    =12
     .                       ,IVCOV_N  =13

        integer, parameter :: IVECNV   =14
     .                       ,IVECOV   =15
cc     .                       ,IVECNV_0 =12
     .                       ,IVECOV_N =16

        integer, parameter :: IVEFCNV  =17
cc     .                       ,IVSCNV   =27

        integer, parameter :: IENI     =18
     .                       ,IENI_0   =19
     .                       ,IEH      =20

        integer, parameter :: IDIVPI   =21
     .                       ,IDIVPE   =22
cc     .                       ,IPCNV    =17
#else
        integer, parameter :: IRHO=1,IVX=2,IVY=3,IVZ=4,IAX=5,IAY=6,IAZ=7
     .                  ,ITMP=8,IBX=9,IBY=10,IBZ=11,IJX=12,IJY=13,IJZ=14

        integer, parameter :: IJCNV   =1
     .                       ,IJCNV_0 =2
cc     .                       ,IJCNV_N =26
     .                       ,IJCOV   =3
     .                       ,IJCOV_0 =4

        integer, parameter :: IBCNV   =5
cc     .                       ,IBCNV_N =22
     .                       ,IBCNV_0 =6
     .                       ,IBCOV   =7
     .                       ,IBCOV_0 =8
cc     .                       ,IBHAT   =9

        integer, parameter :: IACNV   =9
     .                       ,IACOV   =10
     .                       ,IACOV_0 =11

        integer, parameter :: IVCNV   =12
     .                       ,IVCNV_0 =13
     .                       ,IVCNV_N =14
     .                       ,IVCOV   =15
     .                       ,IVCOV_N =16

        integer, parameter :: IVECNV  =17
cc     .                       ,IVECNV_0=19
     .                       ,IVECOV  =18
     .                       ,IVECOV_N=19

        integer, parameter :: IVEFCNV =20
cc     .                       ,IVSCNV  =30

        integer, parameter :: IENI    =21
     .                       ,IEH     =22
     .                       ,IENI_0  =23

        integer, parameter :: IDIVPI  =24
     .                       ,IDIVPE  =25
cc     .                       ,IPCNV   =21
#endif

        integer, parameter :: IETA=1
     .                       ,IETA_0=2
     .                       ,INU=3
cc     .                       ,IPSI=3
     .                       ,IHETA=4
     .                       ,IRHO_0=5
     .                       ,ITMP_0=6

        character(5)  :: equil
        integer :: u_equf=1000,u_prtf=2000,u_mapf=3000
        character(20) :: equ_file,prt_file,map_file

        logical :: adiabatic=.false.   !Whether to use pressure eq with adiabatic EOS
     .            ,ion_hall=.false.    !Whether to use ion EOM for Hall E
     .            ,pinch_flow=.false.  !Whether equilibrium has pinch flow
     .            ,slava     =.false.  !Whether to use Slava Lukin's Hall MHD implementation

      contains

c     setupProblemDims
c     ################################################################
      subroutine setupProblemDims

c     ----------------------------------------------------------------
c     Define problem dimensions
c     ----------------------------------------------------------------

      neqd = 8  !Number of equations

#if !defined(vec_pot)
      NAUXV=22  !Number of auxiliary vectors
#else
      NAUXV=25  !Number of auxiliary vectors
#endif

      NAUXS=6   !Number of auxiliary scalars

      end subroutine setupProblemDims

      end module problem_def

c module transport_params
c ######################################################################
      module transport_params

        use grid

        use problem_def

        real(8) :: nu,eta,dd,chi,chi_par,gamma,di,de,a_p,temp_ratio
     .            ,heta=0d0,c_hyper

        real(8) :: aa_eta,bb_eta,cc_eta
     .            ,aa_nu ,bb_nu ,cc_nu

      end module transport_params

c module auxiliaryVariables
c ######################################################################
      module auxiliaryVariables

        use variables

        use problem_def

        use transport_params

#if defined(vec_pot)
        real(8),pointer,dimension(:,:,:,:) :: acov
     .                                       ,acov_0
     .                                       ,acnv
#endif

        real(8),pointer,dimension(:,:,:)   :: nuu
     .                                       ,eeta
     .                                       ,eeta_0
     .                                       ,h_eta
     .                                       ,psi
     .                                       ,rho_0
     .                                       ,tmp_0

        real(8),pointer,dimension(:,:,:,:) :: bcnv
     .                                       ,bcnv_0
     .                                       ,bcov
     .                                       ,bcov_0
cc     .                                       ,bcnv_n

        real(8),pointer,dimension(:,:,:,:) :: vcnv
     .                                       ,vcnv_0
     .                                       ,vcnv_n
     .                                       ,vcov
     .                                       ,vcov_n
     .                                       ,pcnv

        real(8),pointer,dimension(:,:,:,:) :: bhat
     .                                       ,E_ni
     .                                       ,E_ni_0
     .                                       ,E_h
     .                                       ,div_pi
     .                                       ,div_pe
     .                                       ,jcnv
cc     .                                       ,jcnv_n
     .                                       ,jcnv_0
     .                                       ,jcov
     .                                       ,jcov_0

        real(8),pointer,dimension(:,:,:,:) :: vecnv
cc     .                                       ,vecnv_0
     .                                       ,vecov
     .                                       ,vecov_n
     .                                       ,vscnv
     .                                       ,vefcnv

        !Grid aliases
cc        real(8),pointer,dimension(:) :: xx,yy,zz,dxh,dyh,dzh,dx,dy,dz

      contains

c     allocAuxVariables
c     ###################################################################
      subroutine allocAuxVariables(vaux)

c     -------------------------------------------------------------------
c     Allocates auxiliary variable storage, sets up pointers, defines
c     variable dependencies for BC scheduler.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(aux_array),pointer :: vaux

c     Local variables

        integer :: nx,ny,nz

c     Begin program

c     Allocate storage

        if (.not.associated(vaux)) call allocateAuxStruct(vaux)

c     Setup pointers

        !Equilibrium variables
        vaux%var_list(IETA_0)%descr = 'eta_0'
        eeta_0=>vaux%var_list(IETA_0)%array

        vaux%var_list(IRHO_0)%descr = 'rho_0'
        rho_0=> vaux%var_list(IRHO_0)%array

        vaux%var_list(ITMP_0)%descr = 'tmp_0'
        tmp_0=> vaux%var_list(ITMP_0)%array

#if defined(vec_pot)
        vaux%vec_list(IACOV_0)%cnv   = .false.
        vaux%vec_list(IACOV_0)%descr = 'A_0 cov'
        acov_0 => vaux%vec_list(IACOV_0)%vec
#endif

        vaux%vec_list(IBCNV_0)%cnv   = .true.
        vaux%vec_list(IBCNV_0)%descr = 'B_0 cnv'
        bcnv_0 => vaux%vec_list(IBCNV_0)%vec

        vaux%vec_list(IBCOV_0)%cnv = .false.
        vaux%vec_list(IBCOV_0)%descr = 'B_0 cov'
        bcov_0 => vaux%vec_list(IBCOV_0)%vec

        vaux%vec_list(IJCNV_0)%cnv = .true.
        vaux%vec_list(IJCNV_0)%descr = 'J_0 cnv'
        jcnv_0 => vaux%vec_list(IJCNV_0)%vec

        vaux%vec_list(IJCOV_0)%cnv = .false.
        vaux%vec_list(IJCOV_0)%descr = 'J_0 cov'
        jcov_0 => vaux%vec_list(IJCOV_0)%vec

        vaux%vec_list(IVCNV_0)%cnv = .true.
        vaux%vec_list(IVCNV_0)%descr = 'V_0 cnv'
        vcnv_0 => vaux%vec_list(IVCNV_0)%vec

        vaux%vec_list(IENI_0)%cnv   = .false.
        vaux%vec_list(IENI_0)%descr = 'E_0 cov'
        E_ni_0 => vaux%vec_list(IENI_0)%vec

cc        vaux%vec_list(IVECNV_0)%cnv   = .true.
cc        vaux%vec_list(IVECNV_0)%descr = 'Ve_0 cnv'
cc        vecnv_0 => vaux%vec_list(IVECNV_0)%vec

        !Old-time variables
        vaux%vec_list(IVCOV_N)%cnv = .false.
        vaux%vec_list(IVCOV_N)%descr = 'V_n cov'
        vcov_n => vaux%vec_list(IVCOV_N)%vec

        vaux%vec_list(IVCNV_N)%cnv = .true.
        vaux%vec_list(IVCNV_N)%descr = 'V_n cnv'
        vcnv_n => vaux%vec_list(IVCNV_N)%vec

        vaux%vec_list(IVECOV_N)%cnv   = .false.
        vaux%vec_list(IVECOV_N)%descr = 'Ve_n cov'
        vecov_n => vaux%vec_list(IVECOV_N)%vec

cc        vaux%vec_list(IBCNV_N)%cnv   = .true.
cc        vaux%vec_list(IBCNV_N)%descr = 'B_n cnv'
cc        bcnv_n => vaux%vec_list(IBCNV_N)%vec

cc        vaux%vec_list(IJCNV_N)%cnv = .true.
cc        vaux%vec_list(IJCNV_N)%descr = 'J_n cnv'
cc        jcnv_n => vaux%vec_list(IJCNV_N)%vec

        !New-time variables
        vaux%vec_list(IBCNV)%cnv   = .true.
        vaux%vec_list(IBCNV)%descr = 'B cnv'
        bcnv => vaux%vec_list(IBCNV)%vec

        vaux%vec_list(IBCOV)%cnv   = .false.
        vaux%vec_list(IBCOV)%descr = 'B cov'
        bcov => vaux%vec_list(IBCOV)%vec

        vaux%vec_list(IJCNV)%cnv = .true.
        vaux%vec_list(IJCNV)%descr = 'J cnv'
        jcnv => vaux%vec_list(IJCNV)%vec

        vaux%vec_list(IJCOV)%cnv = .false.
        vaux%vec_list(IJCOV)%descr = 'J cov'
        jcov => vaux%vec_list(IJCOV)%vec

        vaux%vec_list(IVCNV)%cnv   = .true.
        vaux%vec_list(IVCNV)%descr = 'V cnv'
        vcnv => vaux%vec_list(IVCNV)%vec

        vaux%vec_list(IVCOV)%cnv = .false.
        vaux%vec_list(IVCOV)%descr = 'V cov'
        vcov => vaux%vec_list(IVCOV)%vec

#if defined(vec_pot)
        vaux%vec_list(IACNV)%cnv   = .true.
        vaux%vec_list(IACNV)%descr = 'A cnv'
        acnv => vaux%vec_list(IACNV)%vec

        vaux%vec_list(IACOV)%cnv   = .false.
        vaux%vec_list(IACOV)%descr = 'A cov'
        acov => vaux%vec_list(IACOV)%vec
#endif

        vaux%vec_list(IVECNV)%cnv   = .true.
        vaux%vec_list(IVECNV)%descr = 'Ve cnv'
        vecnv => vaux%vec_list(IVECNV)%vec

        vaux%vec_list(IVECOV)%cnv   = .false.
        vaux%vec_list(IVECOV)%descr = 'Ve cov'
        vecov => vaux%vec_list(IVECOV)%vec

        vaux%vec_list(IVEFCNV)%cnv   = .true.
        vaux%vec_list(IVEFCNV)%descr = 'Fake Ve cnv'
        vefcnv => vaux%vec_list(IVEFCNV)%vec

        vaux%vec_list(IENI)%cnv   = .false.
        vaux%vec_list(IENI)%descr = 'E cov'
        E_ni => vaux%vec_list(IENI)%vec

        vaux%vec_list(IEH)%cnv   = .false.
        vaux%vec_list(IEH)%descr = 'Hall E cov'
        E_h => vaux%vec_list(IEH)%vec

        vaux%vec_list(IDIVPI)%cnv   = .true.
        vaux%vec_list(IDIVPI)%descr = 'div(Pi_i)'
        div_pi => vaux%vec_list(IDIVPI)%vec

        vaux%vec_list(IDIVPE)%cnv   = .true.
        vaux%vec_list(IDIVPE)%descr = 'div(Pi_e)'
        div_pe => vaux%vec_list(IDIVPE)%vec

        vaux%var_list(IETA)%descr = 'eta'
        eeta => vaux%var_list(IETA)%array

        vaux%var_list(IHETA)%descr = 'heta'
        h_eta => vaux%var_list(IHETA)%array

        vaux%var_list(INU )%descr = 'nu'
        nuu  => vaux%var_list(INU )%array

cc        vaux%vec_list(IVSCNV)%cnv   = .true.
cc        vaux%vec_list(IVSCNV)%descr = 'Vstar cnv'
cc        vscnv => vaux%vec_list(IVSCNV)%vec

cc        vaux%vec_list(IPCNV)%cnv   = .true.
cc        vaux%vec_list(IPCNV)%descr = 'P=rho*v cnv'
cc        pcnv => vaux%vec_list(IPCNV)%vec

cc        vaux%var_list(IPSI)%descr = 'Psi'
cc        psi  => vaux%var_list(IPSI)%array

cc        vaux%vec_list(IBHAT)%cnv   = .true.
cc        vaux%vec_list(IBHAT)%descr = 'B normalized'
cc        bhat => vaux%vec_list(IBHAT)%vec

c     Define BCs

        call defineAuxBCs(vaux)

c     Setup variable dependencies at boundaries

        call defineAuxBCDeps(vaux)

c     Setup variable dependencies in inner domain

        call defineAuxDomainDeps(vaux)

c     Allocate LOCAL auxiliary variables (not SAMRAI-aware)

        nx = size(tmp_0,1)-2
        ny = size(tmp_0,2)-2
        nz = size(tmp_0,3)-2

        if (.not.associated(pcnv)) then
          allocate(pcnv (0:nx+1,0:ny+1,0:nz+1,3)
     .            ,vscnv(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,bhat (0:nx+1,0:ny+1,0:nz+1,3)
     .            ,psi  (0:nx+1,0:ny+1,0:nz+1))
        endif

c     End program

      end subroutine allocAuxVariables

c     fillLocalAuxVars
c     ###################################################################
      subroutine fillLocalAuxVars(igx,igy,igz,varray,vaux)

      implicit none
c     -------------------------------------------------------------------
c     Fills local auxiliary variables (SAMRAI unaware)
c     -------------------------------------------------------------------

c     Call variables

      integer :: igx,igy,igz

      type(var_array),pointer :: varray

      type(aux_array),pointer :: vaux

c     Local variables

      integer :: i,j,k,ig,jg,kg,nx,ny,nz

      real(8) :: bnorm

      real(8),pointer,dimension(:,:,:) :: rho

c     Begin program

      nx = grid_params%nxv(igx)
      ny = grid_params%nyv(igy)
      nz = grid_params%nzv(igz)

      !Store B/|B| unit vector (for parallel heat conduction and/or boostrap current)
      if (chi_par /= chi) then
        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

              bnorm = vectorNorm(i,j,k,igx,bcnv(i,j,k,:),.false.)

              bhat(i,j,k,:) = bcnv(i,j,k,:)/sqrt(bnorm)
            enddo
          enddo
        enddo
      endif

      !Momentum
      rho => varray%array_var(IRHO)%array

      pcnv(:,:,:,1) = rho*vcnv(:,:,:,1)
      pcnv(:,:,:,2) = rho*vcnv(:,:,:,2)
      pcnv(:,:,:,3) = rho*vcnv(:,:,:,3)

      !Temperature transport velocity, v_star = v - di*j/rho/a_p
      if (slava) then
        vscnv = vcnv            !Also to recover old TM Hall performance results
      else
        vscnv(:,:,:,1) = vcnv(:,:,:,1)-di*jcnv(:,:,:,1)/rho/a_p
        vscnv(:,:,:,2) = vcnv(:,:,:,2)-di*jcnv(:,:,:,2)/rho/a_p
        vscnv(:,:,:,3) = vcnv(:,:,:,3)-di*jcnv(:,:,:,3)/rho/a_p
      endif

      !Psi
      psi = 0d0

cc      !Previous time step B
cc      if (save_tn) bcnv_n = bcnv !Save previous time step
cc
cc      !Previous time step J
cc      if (save_tn .and. ion_hall) jcnv_n = jcnv
cc
cc      !Equilibrium electron velocity
cc      vecnv_0(:,:,:,1) = vcnv_0(:,:,:,1)-di*jcnv_0(:,:,:,1)/rho_0
cc      vecnv_0(:,:,:,2) = vcnv_0(:,:,:,2)-di*jcnv_0(:,:,:,2)/rho_0
cc      vecnv_0(:,:,:,3) = vcnv_0(:,:,:,3)-di*jcnv_0(:,:,:,3)/rho_0

c     End program

      end subroutine fillLocalAuxVars

c     defineAuxBCs
c     ###################################################################
      subroutine defineAuxBCs(vaux)

c     -------------------------------------------------------------------
c     Defines BCs for auxiliary quantities
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(aux_array),pointer :: vaux

c     Local variables

      integer :: ieq,bcs(6,neqd)

c     Begin program

c     Read dependent variable BC setup

      do ieq=1,neqd
        bcs(:,ieq) = gv%u_0%array_var(ieq)%bconds
      enddo

c     Set auxiliary variable BC

      !Equilibrium density and temperature
      vaux%var_list(IRHO_0)%bconds = bcs(:,IRHO)
      vaux%var_list(ITMP_0)%bconds = bcs(:,ITMP)

      !Flow
      vaux%vec_list(IVCNV)%bconds = bcs(:,IVX:IVZ)

      vaux%vec_list(IVCOV)%bconds = vaux%vec_list(IVCNV)%bconds

      vaux%vec_list(IVCNV_0)%bconds = vaux%vec_list(IVCNV)%bconds

#if defined(vec_pot)
      !Vector potential
      vaux%vec_list(IACNV)%bconds = bcs(:,IAX:IAZ)

      vaux%vec_list(IACOV)%bconds = vaux%vec_list(IACNV)%bconds

      !Magnetic field
      vaux%vec_list(IBCNV)%bconds = vaux%vec_list(IACNV)%bconds
      where (vaux%vec_list(IBCNV)%bconds == -EQU)
        vaux%vec_list(IBCNV)%bconds = -NEU
      end where

      vaux%vec_list(IBCOV)%bconds = vaux%vec_list(IBCNV)%bconds
#else
      !Magnetic field
      vaux%vec_list(IBCNV)%bconds = bcs(:,IBX:IBZ)

      !For debugging
cc      where (vaux%vec_list(IBCNV)%bconds == -NEU)
cc        vaux%vec_list(IBCNV)%bconds =-EXT
cc      end where

      vaux%vec_list(IBCOV)%bconds = vaux%vec_list(IBCNV)%bconds
#endif

      !Equilibrium magnetic field
      vaux%vec_list(IBCNV_0)%bconds = vaux%vec_list(IBCNV)%bconds

      where (vaux%vec_list(IBCNV_0)%bconds == -NEU)
        vaux%vec_list(IBCNV_0)%bconds =-EXT
      end where

      vaux%vec_list(IBCOV_0)%bconds = vaux%vec_list(IBCNV_0)%bconds

      !Current
      vaux%vec_list(IJCNV)%bconds = vaux%vec_list(IBCNV)%bconds

      where (bcond == FSYM)
        vaux%vec_list(IJCNV)%bconds(:,1) = FSYM !Recover symmetry BCs
        vaux%vec_list(IJCNV)%bconds(:,2) = FSYM !Recover symmetry BCs
        vaux%vec_list(IJCNV)%bconds(:,3) = FSYM !Recover symmetry BCs
      end where

      where (vaux%vec_list(IJCNV)%bconds == -NEU)
        vaux%vec_list(IJCNV)%bconds = EQU !Use contravariant components for EQU BCs
      end where

      call set_vec_symm_bc(vaux%vec_list(IJCNV)%bconds,-1)  !Same as velocity

      vaux%vec_list(IJCOV)%bconds = vaux%vec_list(IJCNV)%bconds

      !Equilibrium Current (extrapolate covariant components)
      vaux%vec_list(IJCNV_0)%bconds = vaux%vec_list(IJCNV)%bconds
      where (abs(vaux%vec_list(IJCNV_0)%bconds) == EQU)
        vaux%vec_list(IJCNV_0)%bconds =-EXT
      end where
      vaux%vec_list(IJCOV_0)%bconds = vaux%vec_list(IJCNV_0)%bconds

#if defined(vec_pot)
      vaux%vec_list(IACOV_0)%bconds = vaux%vec_list(IJCNV_0)%bconds
#endif

      !Pressure tensors
cc      vaux%vec_list(IDIVPI)%bconds(:,1) = bcond
cc      vaux%vec_list(IDIVPI)%bconds(:,2) = bcond
cc      vaux%vec_list(IDIVPI)%bconds(:,3) = bcond
cc      where (vaux%vec_list(IDIVPI)%bconds == DEF)
cc        vaux%vec_list(IDIVPI)%bconds = EXT
cc      end where
cc
cc      call set_vec_symm_bc(vaux%vec_list(IDIVPI)%bconds,-1)  !Same as velocity

      vaux%vec_list(IDIVPI)%bconds = vaux%vec_list(IVCNV )%bconds

      vaux%vec_list(IDIVPE)%bconds = vaux%vec_list(IDIVPI)%bconds

      !Electric field (Hall piece)
      vaux%vec_list(IEH   )%bconds = vaux%vec_list(IVCOV )%bconds

cc      !Equilibrium electric field (non-ideal piece)
cc      vaux%vec_list(IENI_0)%bconds = vaux%vec_list(IEH   )%bconds

      !Fake electron velocity
      vaux%vec_list(IVEFCNV)%bconds= vaux%vec_list(IJCNV )%bconds

c     End program

      end subroutine defineAuxBCs

c     defineAuxBCDeps
c     ###################################################################
      subroutine defineAuxBCDeps(vaux)

c     -------------------------------------------------------------------
c     Defines intervariable dependencies for auxiliary quantities
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(aux_array),pointer :: vaux

c     Begin program

c     Equilibrium quantities

      !Rho
      vaux%var_list(IRHO_0)%bc_dep_list(1,1) = IRHO_0*AUX
      vaux%var_list(IRHO_0)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IRHO_0)%bc_dep_list(1,3) = 0 !trivial
      vaux%var_list(IRHO_0)%time_0 = .true.

      !Temp
      vaux%var_list(ITMP_0)%bc_dep_list(1,1) = ITMP_0*AUX
      vaux%var_list(ITMP_0)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(ITMP_0)%bc_dep_list(1,3) = 0 !trivial
      vaux%var_list(ITMP_0)%time_0 = .true.

      !Eta
      vaux%var_list(IETA_0)%bc_dep_list(1,1) = ITMP_0*AUX
      vaux%var_list(IETA_0)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IETA_0)%bc_dep_list(1,3) = 0 !trivial
      vaux%var_list(IETA_0)%time_0 = .true.

#if defined(vec_pot) 
      !A cov
      vaux%vec_list(IACOV_0)%bc_dep_list(1,1) = IACOV_0*AUX
      vaux%vec_list(IACOV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IACOV_0)%bc_dep_list(1,3) = 0 !trivial
      vaux%vec_list(IACOV_0)%time_0 = .true.

      !B cnv
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,1) = IACOV_0*AUX
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,3) = 1 !non-trivial
      vaux%vec_list(IBCNV_0)%time_0 = .true.
#else
      !B cnv
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,1) = IBCNV_0*AUX
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_0)%bc_dep_list(1,3) = 0 !trivial
      vaux%vec_list(IBCNV_0)%time_0 = .true.
#endif

      !B cov
      vaux%vec_list(IBCOV_0)%bc_dep_list(1,1) = IBCNV_0*AUX
      vaux%vec_list(IBCOV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCOV_0)%bc_dep_list(1,3) = 0 !trivial
      vaux%vec_list(IBCOV_0)%time_0 = .true.

      !V cnv
      vaux%vec_list(IVCNV_0)%bc_dep_list(1,1) = IRHO_0*AUX
      vaux%vec_list(IVCNV_0)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IVCNV_0)%bc_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IVCNV_0)%bc_dep_list(2,1) = ITMP_0*AUX
      vaux%vec_list(IVCNV_0)%bc_dep_list(2,2) = SCALAR
      vaux%vec_list(IVCNV_0)%bc_dep_list(2,3) = 1 !non-trivial

      vaux%vec_list(IVCNV_0)%bc_dep_list(3,1) = IETA_0*AUX
      vaux%vec_list(IVCNV_0)%bc_dep_list(3,2) = SCALAR
      vaux%vec_list(IVCNV_0)%bc_dep_list(3,3) = 0 !trivial

      vaux%vec_list(IVCNV_0)%bc_dep_list(4,1) = IBCNV_0*AUX
      vaux%vec_list(IVCNV_0)%bc_dep_list(4,2) = VECTOR
      vaux%vec_list(IVCNV_0)%bc_dep_list(4,3) = 0 !trivial

      vaux%vec_list(IVCNV_0)%bc_dep_list(5,1) = IVCNV_0*AUX
      vaux%vec_list(IVCNV_0)%bc_dep_list(5,2) = VECTOR
      vaux%vec_list(IVCNV_0)%bc_dep_list(5,3) = 0 !trivial

      vaux%vec_list(IVCNV_0)%time_0 = .true.

      !J cnv
      vaux%vec_list(IJCNV_0)%bc_dep_list(1,1) = IETA_0*AUX
      vaux%vec_list(IJCNV_0)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IJCNV_0)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IJCNV_0)%bc_dep_list(2,1) = IBCOV_0*AUX
      vaux%vec_list(IJCNV_0)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IJCNV_0)%bc_dep_list(2,3) = 1 !non-trivial

      vaux%vec_list(IJCNV_0)%bc_dep_list(3,1) = IJCNV_0*AUX
      vaux%vec_list(IJCNV_0)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IJCNV_0)%bc_dep_list(3,3) = 0 !trivial

      vaux%vec_list(IJCNV_0)%time_0 = .true.

      !J cov
      vaux%vec_list(IJCOV_0)%bc_dep_list(1,1) = IJCNV_0*AUX
      vaux%vec_list(IJCOV_0)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCOV_0)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IJCOV_0)%time_0 = .true.

      !E-field cov
      vaux%vec_list(IENI_0)%bc_dep_list(1,1) = IETA_0*AUX
      vaux%vec_list(IENI_0)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IENI_0)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IENI_0)%bc_dep_list(2,1) = IBCNV_0*AUX
      vaux%vec_list(IENI_0)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IENI_0)%bc_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IENI_0)%bc_dep_list(3,1) = IVCNV_0*AUX
      vaux%vec_list(IENI_0)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IENI_0)%bc_dep_list(3,3) = 0 !trivial

#if defined(vec_pot)
      vaux%vec_list(IENI_0)%bc_dep_list(4,1) = IACOV_0*AUX
      vaux%vec_list(IENI_0)%bc_dep_list(4,2) = VECTOR
      vaux%vec_list(IENI_0)%bc_dep_list(4,3) = 1 !non-trivial
#else
      vaux%vec_list(IENI_0)%bc_dep_list(4,1) = IJCOV_0*AUX
      vaux%vec_list(IENI_0)%bc_dep_list(4,2) = VECTOR
      vaux%vec_list(IENI_0)%bc_dep_list(4,3) = 0 !trivial
#endif
      vaux%vec_list(IENI_0)%time_0 = .true.

c     Current variables

      !Eta
      vaux%var_list(IETA)%bc_dep_list(1,1) = ITMP
      vaux%var_list(IETA)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IETA)%bc_dep_list(1,3) = 0 !trivial

#if defined(vec_pot)
      !A cov
      vaux%vec_list(IACOV)%bc_dep_list(1,1) = IACOV_0*AUX
      vaux%vec_list(IACOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IACOV)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IACOV)%bc_dep_list(2,1) = IACOV*AUX
      vaux%vec_list(IACOV)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IACOV)%bc_dep_list(2,3) = 0 !trivial

      !A cnv
      vaux%vec_list(IACNV)%bc_dep_list(1,1) = IACOV*AUX
      vaux%vec_list(IACNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IACNV)%bc_dep_list(1,3) = 0 !trivial

      !B cnv
      vaux%vec_list(IBCNV)%bc_dep_list(1,1) = IACOV*AUX
      vaux%vec_list(IBCNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV)%bc_dep_list(1,3) = 1 !non-trivial
#else
      !B cnv
      vaux%vec_list(IBCNV)%bc_dep_list(1,1) = IBCNV*AUX
      vaux%vec_list(IBCNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV)%bc_dep_list(1,3) = 0 !trivial

      if (.not.pinch_flow) then
        vaux%vec_list(IBCNV)%bc_dep_list(2,1) = IJCOV_0*AUX
        vaux%vec_list(IBCNV)%bc_dep_list(2,2) = VECTOR
        vaux%vec_list(IBCNV)%bc_dep_list(2,3) = 0 !trivial
      else
        vaux%vec_list(IBCNV)%bc_dep_list(2,1) = IRHO
        vaux%vec_list(IBCNV)%bc_dep_list(2,2) = SCALAR
        vaux%vec_list(IBCNV)%bc_dep_list(2,3) = 1 !non-trivial

        vaux%vec_list(IBCNV)%bc_dep_list(3,1) = ITMP
        vaux%vec_list(IBCNV)%bc_dep_list(3,2) = SCALAR
        vaux%vec_list(IBCNV)%bc_dep_list(3,3) = 1 !non-trivial

        vaux%vec_list(IBCNV)%bc_dep_list(4,1) = IETA*AUX
        vaux%vec_list(IBCNV)%bc_dep_list(4,2) = SCALAR
        vaux%vec_list(IBCNV)%bc_dep_list(4,3) = 0 !trivial
      endif
#endif
      !B cov
      vaux%vec_list(IBCOV)%bc_dep_list(1,1) = IBCNV*AUX
      vaux%vec_list(IBCOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCOV)%bc_dep_list(1,3) = 0 !trivial

      !Nu=rho*nu_i(T)
      vaux%var_list(INU)%bc_dep_list(1,1) = IRHO
      vaux%var_list(INU)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(INU)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(INU)%bc_dep_list(2,1) = ITMP
      vaux%var_list(INU)%bc_dep_list(2,2) = SCALAR
      vaux%var_list(INU)%bc_dep_list(2,3) = 0 !trivial

      !Heta-rho*nu_e(T)
      vaux%var_list(IHETA)%bc_dep_list(1,1) = IRHO_0*AUX
      vaux%var_list(IHETA)%bc_dep_list(1,2) = SCALAR
      vaux%var_list(IHETA)%bc_dep_list(1,3) = 0 !trivial

      vaux%var_list(IHETA)%bc_dep_list(2,1) = IBCNV_0*AUX
      vaux%var_list(IHETA)%bc_dep_list(2,2) = VECTOR
cc      vaux%var_list(IHETA)%bc_dep_list(2,3) = 0 !trivial
      vaux%var_list(IHETA)%bc_dep_list(2,3) = 1 !trivial

      !V cnv
      vaux%vec_list(IVCNV)%bc_dep_list(1,1) = IVCNV*AUX
      vaux%vec_list(IVCNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCNV)%bc_dep_list(1,3) = 0 !trivial

      if (pinch_flow) then
        vaux%vec_list(IVCNV)%bc_dep_list(2,1) = IRHO
        vaux%vec_list(IVCNV)%bc_dep_list(2,2) = SCALAR
        vaux%vec_list(IVCNV)%bc_dep_list(2,3) = 1 !non-trivial
                                          
        vaux%vec_list(IVCNV)%bc_dep_list(3,1) = ITMP
        vaux%vec_list(IVCNV)%bc_dep_list(3,2) = SCALAR
        vaux%vec_list(IVCNV)%bc_dep_list(3,3) = 1 !non-trivial

        vaux%vec_list(IVCNV)%bc_dep_list(4,1) = IETA*AUX
        vaux%vec_list(IVCNV)%bc_dep_list(4,2) = SCALAR
        vaux%vec_list(IVCNV)%bc_dep_list(4,3) = 0 !trivial

        vaux%vec_list(IVCNV)%bc_dep_list(5,1) = IBCNV*AUX
        vaux%vec_list(IVCNV)%bc_dep_list(5,2) = VECTOR
        vaux%vec_list(IVCNV)%bc_dep_list(5,3) = 0 !trivial
      else
        vaux%vec_list(IVCNV)%bc_dep_list(2,1) = IVCNV_0*AUX
        vaux%vec_list(IVCNV)%bc_dep_list(2,2) = VECTOR
        vaux%vec_list(IVCNV)%bc_dep_list(2,3) = 0 !trivial
      endif

      !V cov
      vaux%vec_list(IVCOV)%bc_dep_list(1,1) = IVCNV*AUX
      vaux%vec_list(IVCOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCOV)%bc_dep_list(1,3) = 0 !trivial

      !J cnv
      vaux%vec_list(IJCNV)%bc_dep_list(1,1) = IJCNV_0*AUX
      vaux%vec_list(IJCNV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCNV)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IJCNV)%bc_dep_list(2,1) = IJCNV*AUX
      vaux%vec_list(IJCNV)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IJCNV)%bc_dep_list(2,3) = 0 !trivial

      if (pinch_flow) then
        vaux%vec_list(IJCNV)%bc_dep_list(3,1) = IRHO
        vaux%vec_list(IJCNV)%bc_dep_list(3,2) = SCALAR
        vaux%vec_list(IJCNV)%bc_dep_list(3,3) = 1 !non-trivial
                                          
        vaux%vec_list(IJCNV)%bc_dep_list(4,1) = ITMP
        vaux%vec_list(IJCNV)%bc_dep_list(4,2) = SCALAR
        vaux%vec_list(IJCNV)%bc_dep_list(4,3) = 1 !non-trivial
                                          
        vaux%vec_list(IJCNV)%bc_dep_list(5,1) = IETA*AUX
        vaux%vec_list(IJCNV)%bc_dep_list(5,2) = SCALAR
        vaux%vec_list(IJCNV)%bc_dep_list(5,3) = 0 !trivial

        vaux%vec_list(IJCNV)%bc_dep_list(6,1) = IBCOV*AUX
        vaux%vec_list(IJCNV)%bc_dep_list(6,2) = VECTOR
        vaux%vec_list(IJCNV)%bc_dep_list(6,3) = 1 !non-trivial
      endif

      !J cov
      vaux%vec_list(IJCOV)%bc_dep_list(1,1) = IJCNV*AUX
      vaux%vec_list(IJCOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCOV)%bc_dep_list(1,3) = 0 !trivial

      !Ve cnv = v - di*j/rho
      vaux%vec_list(IVECNV)%bc_dep_list(1,1) = IRHO
      vaux%vec_list(IVECNV)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IVECNV)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVECNV)%bc_dep_list(2,1) = IVCNV*AUX
      vaux%vec_list(IVECNV)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IVECNV)%bc_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IVECNV)%bc_dep_list(3,1) = IJCNV*AUX
      vaux%vec_list(IVECNV)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IVECNV)%bc_dep_list(3,3) = 0 !trivial

      !Ve cov
      vaux%vec_list(IVECOV)%bc_dep_list(1,1) = IVECNV*AUX
      vaux%vec_list(IVECOV)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IVECOV)%bc_dep_list(1,3) = 0 !trivial

      !Div(pi)
      vaux%vec_list(IDIVPI)%bc_dep_list(1,1) = INU*AUX
      vaux%vec_list(IDIVPI)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IDIVPI)%bc_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IDIVPI)%bc_dep_list(2,1) = IVCNV*AUX
      vaux%vec_list(IDIVPI)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IDIVPI)%bc_dep_list(2,3) = 1 !non-trivial

      vaux%vec_list(IDIVPI)%bc_dep_list(3,1) = IDIVPI*AUX
      vaux%vec_list(IDIVPI)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IDIVPI)%bc_dep_list(3,3) = 0 !trivial

      !Div(pe)
      vaux%vec_list(IDIVPE)%bc_dep_list(1,1) = IHETA*AUX
      vaux%vec_list(IDIVPE)%bc_dep_list(1,2) = SCALAR
      vaux%vec_list(IDIVPE)%bc_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IDIVPE)%bc_dep_list(2,1) = IVECNV*AUX
      vaux%vec_list(IDIVPE)%bc_dep_list(2,2) = VECTOR
      vaux%vec_list(IDIVPE)%bc_dep_list(2,3) = 1 !non-trivial

      vaux%vec_list(IDIVPE)%bc_dep_list(3,1) = IDIVPE*AUX
      vaux%vec_list(IDIVPE)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IDIVPE)%bc_dep_list(3,3) = 0 !trivial

      !E Hall
      if (di > 0d0) then
#if !defined(vec_pot)
        vaux%vec_list(IEH)%bc_dep_list(1,1) = IEH*AUX
        vaux%vec_list(IEH)%bc_dep_list(1,2) = VECTOR
        vaux%vec_list(IEH)%bc_dep_list(1,3) = 0 !trivial
#endif
cc        if (ion_hall) then
cc          vaux%vec_list(IEH)%bc_dep_list(1,1) = IVCOV_N*AUX
cc          vaux%vec_list(IEH)%bc_dep_list(1,2) = VECTOR
cc          vaux%vec_list(IEH)%bc_dep_list(1,3) = 0 !trivial
cc
cc          vaux%vec_list(IEH)%bc_dep_list(2,1) = IVCOV*AUX
cc          vaux%vec_list(IEH)%bc_dep_list(2,2) = VECTOR
cc          vaux%vec_list(IEH)%bc_dep_list(2,3) = 1 !non-trivial
cc
cc          vaux%vec_list(IEH)%bc_dep_list(3,1) = IRHO
cc          vaux%vec_list(IEH)%bc_dep_list(3,2) = SCALAR
cc          vaux%vec_list(IEH)%bc_dep_list(3,3) = 1 !non-trivial
cc
cc          vaux%vec_list(IEH)%bc_dep_list(4,1) = ITMP
cc          vaux%vec_list(IEH)%bc_dep_list(4,2) = SCALAR
cc          vaux%vec_list(IEH)%bc_dep_list(4,3) = 1 !non-trivial
cc
cc          vaux%vec_list(IEH)%bc_dep_list(5,1) = IDIVPI*AUX
cc          vaux%vec_list(IEH)%bc_dep_list(5,2) = VECTOR
cc          vaux%vec_list(IEH)%bc_dep_list(5,3) = 0 !trivial
cc        else
cc          vaux%vec_list(IEH)%bc_dep_list(1,1) = IRHO
cc          vaux%vec_list(IEH)%bc_dep_list(1,2) = SCALAR
cc          vaux%vec_list(IEH)%bc_dep_list(1,3) = 1 !non-trivial
cc
cc          vaux%vec_list(IEH)%bc_dep_list(2,1) = ITMP
cc          vaux%vec_list(IEH)%bc_dep_list(2,2) = SCALAR
cc          vaux%vec_list(IEH)%bc_dep_list(2,3) = 1 !non-trivial
cc
cc          vaux%vec_list(IEH)%bc_dep_list(3,1) = IBCNV*AUX
cc          vaux%vec_list(IEH)%bc_dep_list(3,2) = VECTOR
cc          vaux%vec_list(IEH)%bc_dep_list(3,3) = 0 !trivial
cc
cc          vaux%vec_list(IEH)%bc_dep_list(4,1) = IJCNV*AUX
cc          vaux%vec_list(IEH)%bc_dep_list(4,2) = VECTOR
cc          vaux%vec_list(IEH)%bc_dep_list(4,3) = 0 !trivial
cc
cc          vaux%vec_list(IEH)%bc_dep_list(5,1) = IDIVPE*AUX
cc          vaux%vec_list(IEH)%bc_dep_list(5,2) = VECTOR
cc          vaux%vec_list(IEH)%bc_dep_list(5,3) = 0 !trivial
cc
cc          if (de > 0d0) then
cc            vaux%vec_list(IEH)%bc_dep_list(6,1) = IVECOV_N*AUX
cc            vaux%vec_list(IEH)%bc_dep_list(6,2) = VECTOR
cc            vaux%vec_list(IEH)%bc_dep_list(6,3) = 0 !trivial
cc
cc            vaux%vec_list(IEH)%bc_dep_list(7,1) = IVECOV*AUX
cc            vaux%vec_list(IEH)%bc_dep_list(7,2) = VECTOR
cc            vaux%vec_list(IEH)%bc_dep_list(7,3) = 0 !trivial
cc          endif
cc        endif
      endif

      !E total
#if !defined(vec_pot)
      vaux%vec_list(IENI)%bc_dep_list(1,1) = IENI_0*AUX
      vaux%vec_list(IENI)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IENI)%bc_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IENI)%bc_dep_list(2,1) = IETA*AUX
      vaux%vec_list(IENI)%bc_dep_list(2,2) = SCALAR
      vaux%vec_list(IENI)%bc_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IENI)%bc_dep_list(3,1) = IBCNV*AUX
      vaux%vec_list(IENI)%bc_dep_list(3,2) = VECTOR
      vaux%vec_list(IENI)%bc_dep_list(3,3) = 0 !trivial

      vaux%vec_list(IENI)%bc_dep_list(4,1) = IVCNV*AUX
      vaux%vec_list(IENI)%bc_dep_list(4,2) = VECTOR
      vaux%vec_list(IENI)%bc_dep_list(4,3) = 0 !trivial

      vaux%vec_list(IENI)%bc_dep_list(5,1) = IJCOV*AUX
      vaux%vec_list(IENI)%bc_dep_list(5,2) = VECTOR
      vaux%vec_list(IENI)%bc_dep_list(5,3) = 0 !trivial

      if (di > 0d0) then
        vaux%vec_list(IENI)%bc_dep_list(6,1) = IEH*AUX
        vaux%vec_list(IENI)%bc_dep_list(6,2) = VECTOR
        vaux%vec_list(IENI)%bc_dep_list(6,3) = 0 !trivial
      endif
#else
      vaux%vec_list(IENI)%bc_dep_list(1,1) = IENI*AUX
      vaux%vec_list(IENI)%bc_dep_list(1,2) = VECTOR
      vaux%vec_list(IENI)%bc_dep_list(1,3) = 0 !trivial
#endif

c     End program

      end subroutine defineAuxBCDeps

c     defineAuxDomainDeps
c     ###################################################################
      subroutine defineAuxDomainDeps(vaux)

c     -------------------------------------------------------------------
c     Defines intervariable dependencies for auxiliary quantities
c     -------------------------------------------------------------------

      implicit none

c     Call variables

      type(aux_array),pointer :: vaux

c     Begin program

c     Equilibrium quantities

      !Rho
      vaux%var_list(IRHO_0)%dom_dep_list(1,1) = IRHO_0*AUX
      vaux%var_list(IRHO_0)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(IRHO_0)%dom_dep_list(1,3) = 0 !trivial
      vaux%var_list(IRHO_0)%time_0 = .true.

      !Temp
      vaux%var_list(ITMP_0)%dom_dep_list(1,1) = ITMP_0*AUX
      vaux%var_list(ITMP_0)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(ITMP_0)%dom_dep_list(1,3) = 0 !trivial
      vaux%var_list(ITMP_0)%time_0 = .true.

      !Eta
      vaux%var_list(IETA_0)%dom_dep_list(1,1) = ITMP_0*AUX
      vaux%var_list(IETA_0)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(IETA_0)%dom_dep_list(1,3) = 0 !trivial
      vaux%var_list(IETA_0)%time_0 = .true.

#if defined(vec_pot)
      !A cov
      vaux%vec_list(IACOV_0)%dom_dep_list(1,1) = IACOV_0*AUX
      vaux%vec_list(IACOV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IACOV_0)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IACOV_0)%time_0 = .true.

      !B cnv
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,1) = IACOV_0*AUX
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IBCNV_0)%time_0 = .true.

#else
      !B cnv
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,1) = IBCNV_0*AUX
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV_0)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IBCNV_0)%time_0 = .true.
#endif
      !B cov
      vaux%vec_list(IBCOV_0)%dom_dep_list(1,1) = IBCNV_0*AUX
      vaux%vec_list(IBCOV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCOV_0)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IBCOV_0)%time_0 = .true.

      !V cnv
      vaux%vec_list(IVCNV_0)%dom_dep_list(1,1) = IVCNV_0*AUX
      vaux%vec_list(IVCNV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCNV_0)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVCNV_0)%time_0 = .true.

      !J cnv
      vaux%vec_list(IJCNV_0)%dom_dep_list(1,1) = IBCOV_0*AUX
      vaux%vec_list(IJCNV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCNV_0)%dom_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IJCNV_0)%time_0 = .true.

      !J cov
      vaux%vec_list(IJCOV_0)%dom_dep_list(1,1) = IJCNV_0*AUX
      vaux%vec_list(IJCOV_0)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCOV_0)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IJCOV_0)%time_0 = .true.

      !E-field cov
      vaux%vec_list(IENI_0)%dom_dep_list(1,1) = IETA_0*AUX
      vaux%vec_list(IENI_0)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IENI_0)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IENI_0)%dom_dep_list(2,1) = IBCNV_0*AUX
      vaux%vec_list(IENI_0)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IENI_0)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IENI_0)%dom_dep_list(3,1) = IVCNV_0*AUX
      vaux%vec_list(IENI_0)%dom_dep_list(3,2) = VECTOR
      vaux%vec_list(IENI_0)%dom_dep_list(3,3) = 0 !trivial

      vaux%vec_list(IENI_0)%dom_dep_list(4,1) = IJCOV_0*AUX
      vaux%vec_list(IENI_0)%dom_dep_list(4,2) = VECTOR
      vaux%vec_list(IENI_0)%dom_dep_list(4,3) = 0 !trivial

      vaux%vec_list(IENI_0)%time_0 = .true.

c     Current variables

      !Eta
      vaux%var_list(IETA)%dom_dep_list(1,1) = ITMP
      vaux%var_list(IETA)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(IETA)%dom_dep_list(1,3) = 0 !trivial

#if defined(vec_pot)
      !A cov
      vaux%vec_list(IACOV)%dom_dep_list(1,1) = IAX
      vaux%vec_list(IACOV)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IACOV)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IACOV)%dom_dep_list(2,1) = IAY
      vaux%vec_list(IACOV)%dom_dep_list(2,2) = SCALAR
      vaux%vec_list(IACOV)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IACOV)%dom_dep_list(3,1) = IAZ
      vaux%vec_list(IACOV)%dom_dep_list(3,2) = SCALAR
      vaux%vec_list(IACOV)%dom_dep_list(3,3) = 0 !trivial

      !A cnv
      vaux%vec_list(IACNV)%dom_dep_list(1,1) = IACOV*AUX
      vaux%vec_list(IACNV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IACNV)%dom_dep_list(1,3) = 0 !trivial

      !B cnv
      vaux%vec_list(IBCNV)%dom_dep_list(1,1) = IACOV*AUX
      vaux%vec_list(IBCNV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCNV)%dom_dep_list(1,3) = 1 !non-trivial
#else
      !B cnv
      vaux%vec_list(IBCNV)%dom_dep_list(1,1) = IBX
      vaux%vec_list(IBCNV)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IBCNV)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IBCNV)%dom_dep_list(2,1) = IBY
      vaux%vec_list(IBCNV)%dom_dep_list(2,2) = SCALAR
      vaux%vec_list(IBCNV)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IBCNV)%dom_dep_list(3,1) = IBZ
      vaux%vec_list(IBCNV)%dom_dep_list(3,2) = SCALAR
      vaux%vec_list(IBCNV)%dom_dep_list(3,3) = 0 !trivial
#endif
      !B cov
      vaux%vec_list(IBCOV)%dom_dep_list(1,1) = IBCNV*AUX
      vaux%vec_list(IBCOV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IBCOV)%dom_dep_list(1,3) = 0 !trivial

      !Nu=rho*nu_i(T)
      vaux%var_list(INU)%dom_dep_list(1,1) = IRHO
      vaux%var_list(INU)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(INU)%dom_dep_list(1,3) = 0 !trivial

      vaux%var_list(INU)%dom_dep_list(2,1) = ITMP
      vaux%var_list(INU)%dom_dep_list(2,2) = SCALAR
      vaux%var_list(INU)%dom_dep_list(2,3) = 0 !trivial

      !Heta-rho*nu_e(T)
      vaux%var_list(IHETA)%dom_dep_list(1,1) = IRHO_0*AUX
      vaux%var_list(IHETA)%dom_dep_list(1,2) = SCALAR
      vaux%var_list(IHETA)%dom_dep_list(1,3) = 0 !trivial

      vaux%var_list(IHETA)%dom_dep_list(2,1) = IBCNV_0*AUX
      vaux%var_list(IHETA)%dom_dep_list(2,2) = VECTOR
      vaux%var_list(IHETA)%dom_dep_list(2,3) = 0 !trivial

      !V cnv
      vaux%vec_list(IVCNV)%dom_dep_list(1,1) = IVX
      vaux%vec_list(IVCNV)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IVCNV)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVCNV)%dom_dep_list(2,1) = IVY
      vaux%vec_list(IVCNV)%dom_dep_list(2,2) = SCALAR
      vaux%vec_list(IVCNV)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IVCNV)%dom_dep_list(3,1) = IVZ
      vaux%vec_list(IVCNV)%dom_dep_list(3,2) = SCALAR
      vaux%vec_list(IVCNV)%dom_dep_list(3,3) = 0 !trivial

      !V cov
      vaux%vec_list(IVCOV)%dom_dep_list(1,1) = IVCNV*AUX
      vaux%vec_list(IVCOV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IVCOV)%dom_dep_list(1,3) = 0 !trivial

      !J cnv
      vaux%vec_list(IJCNV)%dom_dep_list(1,1) = IBCOV*AUX
      vaux%vec_list(IJCNV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCNV)%dom_dep_list(1,3) = 1 !non-trivial

      !J cov
      vaux%vec_list(IJCOV)%dom_dep_list(1,1) = IJCNV*AUX
      vaux%vec_list(IJCOV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IJCOV)%dom_dep_list(1,3) = 0 !trivial

      !Ve cnv = v - di*j/rho
      vaux%vec_list(IVECNV)%dom_dep_list(1,1) = IRHO
      vaux%vec_list(IVECNV)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IVECNV)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IVECNV)%dom_dep_list(2,1) = IVCNV*AUX
      vaux%vec_list(IVECNV)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IVECNV)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IVECNV)%dom_dep_list(3,1) = IJCNV*AUX
      vaux%vec_list(IVECNV)%dom_dep_list(3,2) = VECTOR
      vaux%vec_list(IVECNV)%dom_dep_list(3,3) = 0 !trivial

      !Ve cov
      vaux%vec_list(IVECOV)%dom_dep_list(1,1) = IVECNV*AUX
      vaux%vec_list(IVECOV)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IVECOV)%dom_dep_list(1,3) = 0 !trivial

      !Div(pi)
      vaux%vec_list(IDIVPI)%dom_dep_list(1,1) = INU*AUX
      vaux%vec_list(IDIVPI)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IDIVPI)%dom_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IDIVPI)%dom_dep_list(2,1) = IVCNV*AUX
      vaux%vec_list(IDIVPI)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IDIVPI)%dom_dep_list(2,3) = 1 !non-trivial

      !Div(pe)
      vaux%vec_list(IDIVPE)%dom_dep_list(1,1) = IHETA*AUX
      vaux%vec_list(IDIVPE)%dom_dep_list(1,2) = SCALAR
      vaux%vec_list(IDIVPE)%dom_dep_list(1,3) = 1 !non-trivial

      vaux%vec_list(IDIVPE)%dom_dep_list(2,1) = IVECNV*AUX
      vaux%vec_list(IDIVPE)%dom_dep_list(2,2) = VECTOR
      vaux%vec_list(IDIVPE)%dom_dep_list(2,3) = 1 !non-trivial

      !E Hall
      if (di > 0d0) then
        if (ion_hall) then
          vaux%vec_list(IEH)%dom_dep_list(1,1) = IVCOV_N*AUX
          vaux%vec_list(IEH)%dom_dep_list(1,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(1,3) = 0 !trivial

          vaux%vec_list(IEH)%dom_dep_list(2,1) = IVCOV*AUX
          vaux%vec_list(IEH)%dom_dep_list(2,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(2,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(3,1) = IRHO
          vaux%vec_list(IEH)%dom_dep_list(3,2) = SCALAR
          vaux%vec_list(IEH)%dom_dep_list(3,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(4,1) = ITMP
          vaux%vec_list(IEH)%dom_dep_list(4,2) = SCALAR
          vaux%vec_list(IEH)%dom_dep_list(4,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(5,1) = IDIVPI*AUX
          vaux%vec_list(IEH)%dom_dep_list(5,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(5,3) = 0 !trivial
        else
          vaux%vec_list(IEH)%dom_dep_list(1,1) = IRHO
          vaux%vec_list(IEH)%dom_dep_list(1,2) = SCALAR
          vaux%vec_list(IEH)%dom_dep_list(1,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(2,1) = ITMP
          vaux%vec_list(IEH)%dom_dep_list(2,2) = SCALAR
          vaux%vec_list(IEH)%dom_dep_list(2,3) = 1 !non-trivial

          vaux%vec_list(IEH)%dom_dep_list(3,1) = IBCNV*AUX
          vaux%vec_list(IEH)%dom_dep_list(3,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(3,3) = 0 !trivial

          vaux%vec_list(IEH)%dom_dep_list(4,1) = IJCNV*AUX
          vaux%vec_list(IEH)%dom_dep_list(4,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(4,3) = 0 !trivial

          vaux%vec_list(IEH)%dom_dep_list(5,1) = IDIVPE*AUX
          vaux%vec_list(IEH)%dom_dep_list(5,2) = VECTOR
          vaux%vec_list(IEH)%dom_dep_list(5,3) = 0 !trivial

          if (de > 0d0) then
            vaux%vec_list(IEH)%dom_dep_list(6,1) = IVECOV_N*AUX
            vaux%vec_list(IEH)%dom_dep_list(6,2) = VECTOR
            vaux%vec_list(IEH)%dom_dep_list(6,3) = 0 !trivial

            vaux%vec_list(IEH)%dom_dep_list(7,1) = IVECOV*AUX
            vaux%vec_list(IEH)%dom_dep_list(7,2) = VECTOR
            vaux%vec_list(IEH)%dom_dep_list(7,3) = 0 !trivial
          endif
        endif
      endif

      !E total
      vaux%vec_list(IENI)%dom_dep_list(1,1) = IENI_0*AUX
      vaux%vec_list(IENI)%dom_dep_list(1,2) = VECTOR
      vaux%vec_list(IENI)%dom_dep_list(1,3) = 0 !trivial

      vaux%vec_list(IENI)%dom_dep_list(2,1) = IETA*AUX
      vaux%vec_list(IENI)%dom_dep_list(2,2) = SCALAR
      vaux%vec_list(IENI)%dom_dep_list(2,3) = 0 !trivial

      vaux%vec_list(IENI)%dom_dep_list(3,1) = IBCNV*AUX
      vaux%vec_list(IENI)%dom_dep_list(3,2) = VECTOR
      vaux%vec_list(IENI)%dom_dep_list(3,3) = 0 !trivial

      vaux%vec_list(IENI)%dom_dep_list(4,1) = IVCNV*AUX
      vaux%vec_list(IENI)%dom_dep_list(4,2) = VECTOR
      vaux%vec_list(IENI)%dom_dep_list(4,3) = 0 !trivial

      vaux%vec_list(IENI)%dom_dep_list(5,1) = IJCOV*AUX
      vaux%vec_list(IENI)%dom_dep_list(5,2) = VECTOR
      vaux%vec_list(IENI)%dom_dep_list(5,3) = 0 !trivial

      if (di > 0d0) then
        vaux%vec_list(IENI)%dom_dep_list(6,1) = IEH*AUX
        vaux%vec_list(IENI)%dom_dep_list(6,2) = VECTOR
        vaux%vec_list(IENI)%dom_dep_list(6,3) = 0 !trivial
      endif

c     End program

      end subroutine defineAuxDomainDeps

c     deallocLocalAuxVars
c     ###################################################################
      subroutine deallocLocalAuxVars

c     -------------------------------------------------------------------
c     Deallocates LOCAL auxiliary variable storage.
c     -------------------------------------------------------------------

      integer :: istat

c     Begin program

      deallocate(pcnv,vscnv,psi,bhat,stat=istat)

      nullify(bcnv
     .       ,bcnv_0
     .       ,bcov
     .       ,bcov_0
     .       ,jcnv
     .       ,jcnv_0
     .       ,jcov
     .       ,jcov_0
     .       ,vcnv
     .       ,vcnv_0
     .       ,vcnv_n
     .       ,vcov
     .       ,vcov_n
#if defined(vec_pot)
     .       ,acnv
     .       ,acov
     .       ,acov_0
#endif
     .       ,vecnv
     .       ,vecov
     .       ,vecov_n
     .       ,vefcnv
     .       ,E_ni
     .       ,E_ni_0
     .       ,E_h
     .       ,div_pi
     .       ,div_pe
     .       ,eeta
     .       ,eeta_0
     .       ,h_eta
     .       ,nuu
     .       ,rho_0
     .       ,tmp_0
cc     .       ,bcnv_n
cc     .       ,jcnv_n
cc     .       ,vecnv_0
cc     .       ,vscnv
cc     .       ,pcnv
cc     .       ,psi
cc     .       ,bhat
     .       )

c     End program

      end subroutine deallocLocalAuxVars

      end module auxiliaryVariables

c module auxPlotVariables
c ######################################################################
      module auxPlotVariables

        use variables

#if defined(vec_pot)
        real(8),target,allocatable,dimension(:,:,:) ::
     .          ax_car,ay_car,az_car
     .         ,ax_cnv,ay_cnv,az_cnv
     .         ,bx,by,bz
#endif

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_car,by_car,bz_car
     .         ,jx_car,jy_car,jz_car
     .         ,vx_car,vy_car,vz_car
     .         ,divrgB,divrgJ,divrgV,Pflux
     .         ,qfactor,lambda,p_tot
     .         ,vpar,vpsi,vperp,jpar,jpsi,jperp

        real(8),target,allocatable,dimension(:,:,:) ::
     .          bx_cov,by_cov,bz_cov
     .         ,jx,jy,jz,jx_cov,jy_cov,jz_cov
     .         ,vx,vy,vz,vx_cov,vy_cov,vz_cov
     .         ,vex,vey,vez

      contains

c     allocAuxPlotVar
c     ##################################################################
      subroutine allocAuxPlotVar

      implicit none

#if defined(vec_pot)
      allocate (bx    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,by    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,ax_cnv(ilom:ihip,jlom:jhip,klom:khip)
     .         ,ay_cnv(ilom:ihip,jlom:jhip,klom:khip)
     .         ,az_cnv(ilom:ihip,jlom:jhip,klom:khip))

      allocate (ax_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,ay_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,az_car(ilom:ihip,jlom:jhip,klom:khip))
#endif

      allocate (bx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .         ,by_cov(ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz_cov(ilom:ihip,jlom:jhip,klom:khip))

      allocate (jx    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jy    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jz    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .         ,jy_cov(ilom:ihip,jlom:jhip,klom:khip)
     .         ,jz_cov(ilom:ihip,jlom:jhip,klom:khip))

      allocate (vx    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz    (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vx_cov(ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy_cov(ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz_cov(ilom:ihip,jlom:jhip,klom:khip)
     .         ,vex   (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vey   (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vez   (ilom:ihip,jlom:jhip,klom:khip))

      allocate (bx_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,by_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,bz_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,jx_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,jy_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,jz_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,vx_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,vy_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,vz_car(ilom:ihip,jlom:jhip,klom:khip)
     .         ,vpar  (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vpsi  (ilom:ihip,jlom:jhip,klom:khip)
     .         ,vperp (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jpar  (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jpsi  (ilom:ihip,jlom:jhip,klom:khip)
     .         ,jperp (ilom:ihip,jlom:jhip,klom:khip)
     .         ,divrgB(ilom:ihip,jlom:jhip,klom:khip)
     .         ,divrgJ(ilom:ihip,jlom:jhip,klom:khip)
     .         ,divrgV(ilom:ihip,jlom:jhip,klom:khip)
     .         ,Pflux (ilom:ihip,jlom:jhip,klom:khip)
     .         ,p_tot (ilom:ihip,jlom:jhip,klom:khip))

      allocate (qfactor(ilom:ihip,jlom:jhip,klom:khip)
     .         ,lambda (ilom:ihip,jlom:jhip,klom:khip))

      end subroutine allocAuxPlotVar

c     deallocAuxPlotVar
c     ##################################################################
      subroutine deallocAuxPlotVar

      implicit none

#if defined(vec_pot)
      deallocate (bx,by,bz,ax_cnv,ay_cnv,az_cnv)
      deallocate (ax_car,ay_car,az_car)
#endif

      deallocate (bx_cov,by_cov,bz_cov)

      deallocate (jx,jy,jz,jx_cov,jy_cov,jz_cov)

      deallocate (vx,vy,vz,vx_cov,vy_cov,vz_cov,vex,vey,vez)

      deallocate (bx_car,by_car,bz_car,jx_car,jy_car,jz_car
     .           ,vx_car,vy_car,vz_car,vpar,vpsi,vperp 
     .           ,divrgB,divrgJ,divrgV,Pflux,p_tot
     .           ,jpar,jpsi,jperp)

      deallocate (qfactor,lambda)

      end subroutine deallocAuxPlotVar

      end module auxPlotVariables

c module operators
c ######################################################################
      module operators

        use error

        use transport_params

      contains

c     curl_bxv
c     ###################################################################
      function curl_bxv(i,j,k,nx,ny,nz,igx,igy,igz,vv,bb,half_elem)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of curl(bb x vv) at the
c     grid node (i,j,k). One sided derivatives are employed when half_elem=1
c     (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nx,ny,nz,half_elem,igx,igy,igz
        real(8)    :: cnv(3)
        real(8)    :: vv(0:nx+1,0:ny+1,0:nz+1,3)
     $               ,bb(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq
        integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

        real(8)    :: idhx,idhy,idhz,a(3)
        real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .               ,jacp,jacm,jacph,jacmh,jach,jac0
     .               ,ijacip,ijacim,ijacjp,ijacjm,ijackp,ijackm
     .               ,ijacp,ijacm,ijacph,ijacmh,ijach,ijac0

        real(8)    :: vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm

        real(8)    :: bxip,bxim,bxjp,bxjm,bxkp,bxkm
     .               ,byip,byim,byjp,byjm,bykp,bykm
     .               ,bzip,bzim,bzjp,bzjm,bzkp,bzkm

        real(8)    :: one=1d0,half=0.5d0

c     Begin program

c     Defaults

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        idhx = half/grid_params%dxh(ig)
        idhy = half/grid_params%dyh(jg)
        idhz = half/grid_params%dzh(kg)

cc        jac  = gmetric%grid(igx)%jac(i,j,k)

cc        if (coords == 'car') then
cc          ijacip  = one
cc          ijacim  = one
cc          ijacjp  = one
cc          ijacjm  = one
cc          ijackp  = one
cc          ijackm  = one
cc        endif

c     Exceptions

        select case(half_elem)
        case (1)
          idhx = one/grid_params%dx(ig)
          im = i

          if (coords /= 'car') then
            jacip  = gmetric%grid(igx)%jac(ip,j,k)
            jacim  = gmetric%grid(igx)%jac(i ,j,k)
            jacjp  = half*(gmetric%grid(igx)%jac(ip,jp,k)
     .                    +gmetric%grid(igx)%jac(i ,jp,k))
            jacjm  = half*(gmetric%grid(igx)%jac(ip,jm,k)
     .                    +gmetric%grid(igx)%jac(i ,jm,k))
            jackp  = half*(gmetric%grid(igx)%jac(ip,j,kp)
     .                    +gmetric%grid(igx)%jac(i ,j,kp))
            jackm  = half*(gmetric%grid(igx)%jac(ip,j,km)
     .                    +gmetric%grid(igx)%jac(i ,j,km))

            if (isSP2(i+1,igx,dim=1)) then
              jacjp = SP_flsv
              jacjm = SP_flsv
              jackp = SP_flsv
              jackm = SP_flsv
cc              write (*,*) 'DIAG -- curl_bxv'
cc              jacjp = one!SP_flsv
cc              jacjm = one!SP_flsv
cc              jackp = one!SP_flsv
cc              jackm = one!SP_flsv
            endif
          endif

          vxip = vv(ip,j,k,1)
          vxim = vv(i ,j,k,1)
          vyip = vv(ip,j,k,2)
          vyim = vv(i ,j,k,2)
          vzip = vv(ip,j,k,3)
          vzim = vv(i ,j,k,3)

          vxjp = half*(vv(ip,jp,k,1)+vv(i,jp,k,1))
          vxjm = half*(vv(ip,jm,k,1)+vv(i,jm,k,1))
          vyjp = half*(vv(ip,jp,k,2)+vv(i,jp,k,2))
          vyjm = half*(vv(ip,jm,k,2)+vv(i,jm,k,2))
          vzjp = half*(vv(ip,jp,k,3)+vv(i,jp,k,3))
          vzjm = half*(vv(ip,jm,k,3)+vv(i,jm,k,3))

          vxkp = half*(vv(ip,j,kp,1)+vv(i,j,kp,1))
          vxkm = half*(vv(ip,j,km,1)+vv(i,j,km,1))
          vykp = half*(vv(ip,j,kp,2)+vv(i,j,kp,2))
          vykm = half*(vv(ip,j,km,2)+vv(i,j,km,2))
          vzkp = half*(vv(ip,j,kp,3)+vv(i,j,kp,3))
          vzkm = half*(vv(ip,j,km,3)+vv(i,j,km,3))

          bxip = bb(ip,j,k,1)
          bxim = bb(i ,j,k,1)
          byip = bb(ip,j,k,2)
          byim = bb(i ,j,k,2)
          bzip = bb(ip,j,k,3)
          bzim = bb(i ,j,k,3)

          bxjp = half*(bb(ip,jp,k,1)+bb(i,jp,k,1))
          bxjm = half*(bb(ip,jm,k,1)+bb(i,jm,k,1))
          byjp = half*(bb(ip,jp,k,2)+bb(i,jp,k,2))
          byjm = half*(bb(ip,jm,k,2)+bb(i,jm,k,2))
          bzjp = half*(bb(ip,jp,k,3)+bb(i,jp,k,3))
          bzjm = half*(bb(ip,jm,k,3)+bb(i,jm,k,3))

          bxkp = half*(bb(ip,j,kp,1)+bb(i,j,kp,1))
          bxkm = half*(bb(ip,j,km,1)+bb(i,j,km,1))
          bykp = half*(bb(ip,j,kp,2)+bb(i,j,kp,2))
          bykm = half*(bb(ip,j,km,2)+bb(i,j,km,2))
          bzkp = half*(bb(ip,j,kp,3)+bb(i,j,kp,3))
          bzkm = half*(bb(ip,j,km,3)+bb(i,j,km,3))

        case (2)

          idhy = one/grid_params%dy(jg)
          jm = j

          if (coords /= 'car') then
            jacip  = half*(gmetric%grid(igx)%jac(ip,jp,k)
     .                    +gmetric%grid(igx)%jac(ip,j ,k))
            jacim  = half*(gmetric%grid(igx)%jac(im,jp,k)
     .                    +gmetric%grid(igx)%jac(im,j ,k))
            jacjp  = gmetric%grid(igx)%jac(i,jp,k)
            jacjm  = gmetric%grid(igx)%jac(i,j ,k)
            jackp  = half*(gmetric%grid(igx)%jac(i,jp,kp)
     .                    +gmetric%grid(igx)%jac(i,j ,kp))
            jackm  = half*(gmetric%grid(igx)%jac(i,jp,km)
     .                    +gmetric%grid(igx)%jac(i,j ,km))
          endif

          vxip = half*(vv(ip,j,k,1)+vv(ip,jp,k,1))
          vxim = half*(vv(im,j,k,1)+vv(im,jp,k,1))
          vyip = half*(vv(ip,j,k,2)+vv(ip,jp,k,2))
          vyim = half*(vv(im,j,k,2)+vv(im,jp,k,2))
          vzip = half*(vv(ip,j,k,3)+vv(ip,jp,k,3))
          vzim = half*(vv(im,j,k,3)+vv(im,jp,k,3))
                      
          vxjp = vv(i,jp,k,1)
          vxjm = vv(i,j ,k,1)
          vyjp = vv(i,jp,k,2)
          vyjm = vv(i,j ,k,2)
          vzjp = vv(i,jp,k,3)
          vzjm = vv(i,j ,k,3)
                      
          vxkp = half*(vv(i,j,kp,1)+vv(i,jp,kp,1))
          vxkm = half*(vv(i,j,km,1)+vv(i,jp,km,1))
          vykp = half*(vv(i,j,kp,2)+vv(i,jp,kp,2))
          vykm = half*(vv(i,j,km,2)+vv(i,jp,km,2))
          vzkp = half*(vv(i,j,kp,3)+vv(i,jp,kp,3))
          vzkm = half*(vv(i,j,km,3)+vv(i,jp,km,3))
                      
          bxip = half*(bb(ip,j,k,1)+bb(ip,jp,k,1))
          bxim = half*(bb(im,j,k,1)+bb(im,jp,k,1))
          byip = half*(bb(ip,j,k,2)+bb(ip,jp,k,2))
          byim = half*(bb(im,j,k,2)+bb(im,jp,k,2))
          bzip = half*(bb(ip,j,k,3)+bb(ip,jp,k,3))
          bzim = half*(bb(im,j,k,3)+bb(im,jp,k,3))

          bxjp = bb(i,jp,k,1)
          bxjm = bb(i,j ,k,1)
          byjp = bb(i,jp,k,2)
          byjm = bb(i,j ,k,2)
          bzjp = bb(i,jp,k,3)
          bzjm = bb(i,j ,k,3)

          bxkp = half*(bb(i,j,kp,1)+bb(i,jp,kp,1))
          bxkm = half*(bb(i,j,km,1)+bb(i,jp,km,1))
          bykp = half*(bb(i,j,kp,2)+bb(i,jp,kp,2))
          bykm = half*(bb(i,j,km,2)+bb(i,jp,km,2))
          bzkp = half*(bb(i,j,kp,3)+bb(i,jp,kp,3))
          bzkm = half*(bb(i,j,km,3)+bb(i,jp,km,3))

        case (3)
          idhz = one/grid_params%dz(kg)
          km = k

          if (coords /= 'car') then
            jacip  = half*(gmetric%grid(igx)%jac(ip,j,kp)
     .                    +gmetric%grid(igx)%jac(ip,j,k ))
            jacim  = half*(gmetric%grid(igx)%jac(im,j,kp)
     .                    +gmetric%grid(igx)%jac(im,j,k ))
            jacjp  = half*(gmetric%grid(igx)%jac(i,jp,kp)
     .                    +gmetric%grid(igx)%jac(i,jp,k ))
            jacjm  = half*(gmetric%grid(igx)%jac(i,jm,kp)
     .                    +gmetric%grid(igx)%jac(i,jm,k ))
            jackp  = gmetric%grid(igx)%jac(i,j,kp)
            jackm  = gmetric%grid(igx)%jac(i,j,km)
          endif

          vxip = half*(vv(ip,j,k,1)+vv(ip,j,kp,1))
          vxim = half*(vv(im,j,k,1)+vv(im,j,kp,1))
          vyip = half*(vv(ip,j,k,2)+vv(ip,j,kp,2))
          vyim = half*(vv(im,j,k,2)+vv(im,j,kp,2))
          vzip = half*(vv(ip,j,k,3)+vv(ip,j,kp,3))
          vzim = half*(vv(im,j,k,3)+vv(im,j,kp,3))
                      
          vxjp = half*(vv(i,jp,k,1)+vv(i,jp,kp,1))
          vxjm = half*(vv(i,jm,k,1)+vv(i,jm,kp,1))
          vyjp = half*(vv(i,jp,k,2)+vv(i,jp,kp,2))
          vyjm = half*(vv(i,jm,k,2)+vv(i,jm,kp,2))
          vzjp = half*(vv(i,jp,k,3)+vv(i,jp,kp,3))
          vzjm = half*(vv(i,jm,k,3)+vv(i,jm,kp,3))
                      
          vxkp = vv(i,j,kp,1)
          vxkm = vv(i,j,k ,1)
          vykp = vv(i,j,kp,2)
          vykm = vv(i,j,k ,2)
          vzkp = vv(i,j,kp,3)
          vzkm = vv(i,j,k ,3)
                      
          bxip = half*(bb(ip,j,k,1)+bb(ip,j,kp,1))
          bxim = half*(bb(im,j,k,1)+bb(im,j,kp,1))
          byip = half*(bb(ip,j,k,2)+bb(ip,j,kp,2))
          byim = half*(bb(im,j,k,2)+bb(im,j,kp,2))
          bzip = half*(bb(ip,j,k,3)+bb(ip,j,kp,3))
          bzim = half*(bb(im,j,k,3)+bb(im,j,kp,3))
                      
          bxjp = half*(bb(i,jp,k,1)+bb(i,jp,kp,1))
          bxjm = half*(bb(i,jm,k,1)+bb(i,jm,kp,1))
          byjp = half*(bb(i,jp,k,2)+bb(i,jp,kp,2))
          byjm = half*(bb(i,jm,k,2)+bb(i,jm,kp,2))
          bzjp = half*(bb(i,jp,k,3)+bb(i,jp,kp,3))
          bzjm = half*(bb(i,jm,k,3)+bb(i,jm,kp,3))

          bxkp = bb(i,j,kp,1)
          bxkm = bb(i,j,k ,1)
          bykp = bb(i,j,kp,2)
          bykm = bb(i,j,k ,2)
          bzkp = bb(i,j,kp,3)
          bzkm = bb(i,j,k ,3)

        case default

          if (coords /= 'car') then
            jacip  = gmetric%grid(igx)%jac(ip,j,k)
            jacim  = gmetric%grid(igx)%jac(im,j,k)
            jacjp  = gmetric%grid(igx)%jac(i,jp,k)
            jacjm  = gmetric%grid(igx)%jac(i,jm,k)
            jackp  = gmetric%grid(igx)%jac(i,j,kp)
            jackm  = gmetric%grid(igx)%jac(i,j,km)
          endif

          vxip = vv(ip,j,k,1)
          vxim = vv(im,j,k,1)
          vyip = vv(ip,j,k,2)
          vyim = vv(im,j,k,2)
          vzip = vv(ip,j,k,3)
          vzim = vv(im,j,k,3)

          vxjp = vv(i,jp,k,1)
          vxjm = vv(i,jm,k,1)
          vyjp = vv(i,jp,k,2)
          vyjm = vv(i,jm,k,2)
          vzjp = vv(i,jp,k,3)
          vzjm = vv(i,jm,k,3)

          vxkp = vv(i,j,kp,1)
          vxkm = vv(i,j,km,1)
          vykp = vv(i,j,kp,2)
          vykm = vv(i,j,km,2)
          vzkp = vv(i,j,kp,3)
          vzkm = vv(i,j,km,3)

          bxip = bb(ip,j,k,1)
          bxim = bb(im,j,k,1)
          byip = bb(ip,j,k,2)
          byim = bb(im,j,k,2)
          bzip = bb(ip,j,k,3)
          bzim = bb(im,j,k,3)

          bxjp = bb(i,jp,k,1)
          bxjm = bb(i,jm,k,1)
          byjp = bb(i,jp,k,2)
          byjm = bb(i,jm,k,2)
          bzjp = bb(i,jp,k,3)
          bzjm = bb(i,jm,k,3)

          bxkp = bb(i,j,kp,1)
          bxkm = bb(i,j,km,1)
          bykp = bb(i,j,kp,2)
          bykm = bb(i,j,km,2)
          bzkp = bb(i,j,kp,3)
          bzkm = bb(i,j,km,3)

        end select

c     Components

        if (coords /= 'car') then
          ijacip  = one/jacip
          ijacim  = one/jacim
          ijacjp  = one/jacjp
          ijacjm  = one/jacjm
          ijackp  = one/jackp
          ijackm  = one/jackm

          !component 1

          flxjp = ( vyjp*bxjp-vxjp*byjp )*ijacjp
          flxjm = ( vyjm*bxjm-vxjm*byjm )*ijacjm

          flxkp = ( vzkp*bxkp-vxkp*bzkp )*ijackp
          flxkm = ( vzkm*bxkm-vxkm*bzkm )*ijackm

          cnv(1) =  (flxjp-flxjm)*idhy
     .             +(flxkp-flxkm)*idhz

          !component 2

          flxip = ( vxip*byip-vyip*bxip )*ijacip
          flxim = ( vxim*byim-vyim*bxim )*ijacim

          flxkp = ( vzkp*bykp-vykp*bzkp )*ijackp
          flxkm = ( vzkm*bykm-vykm*bzkm )*ijackm

          cnv(2) =  (flxip-flxim)*idhx
     .             +(flxkp-flxkm)*idhz

          !component 3

          flxip = ( vxip*bzip-vzip*bxip )*ijacip
          flxim = ( vxim*bzim-vzim*bxim )*ijacim

          flxjp = ( vyjp*bzjp-vzjp*byjp )*ijacjp
          flxjm = ( vyjm*bzjm-vzjm*byjm )*ijacjm

          cnv(3) =  (flxip-flxim)*idhx
     .             +(flxjp-flxjm)*idhy

        else

          !component 1

          flxjp = ( vyjp*bxjp-vxjp*byjp )
          flxjm = ( vyjm*bxjm-vxjm*byjm )

          flxkp = ( vzkp*bxkp-vxkp*bzkp )
          flxkm = ( vzkm*bxkm-vxkm*bzkm )

          cnv(1) =  (flxjp-flxjm)*idhy
     .             +(flxkp-flxkm)*idhz

          !component 2

          flxip = ( vxip*byip-vyip*bxip )
          flxim = ( vxim*byim-vyim*bxim )

          flxkp = ( vzkp*bykp-vykp*bzkp )
          flxkm = ( vzkm*bykm-vykm*bzkm )

          cnv(2) =  (flxip-flxim)*idhx
     .             +(flxkp-flxkm)*idhz

          !component 3

          flxip = ( vxip*bzip-vzip*bxip )
          flxim = ( vxim*bzim-vzim*bxim )

          flxjp = ( vyjp*bzjp-vzjp*byjp )
          flxjm = ( vyjm*bzjm-vzjm*byjm )

          cnv(3) =  (flxip-flxim)*idhx
     .             +(flxjp-flxjm)*idhy

        endif

      end function curl_bxv

ccc     curl_bxv_vrtx
ccc     ###################################################################
cc      function curl_bxv_vrtx(i,j,k,nx,ny,nz,igx,igy,igz,vv,bb)
cc     .         result(cnv)
cc
ccc     -------------------------------------------------------------------
ccc     Finds contravariant components (a1,a2,a3) of curl(bb x vv) at the
ccc     grid node (i,j,k). One sided derivatives are employed when hex=1
ccc     (i+1/1), hey=1 (j+1/2), and hez=1 (k+1/2).
ccc     -------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,nx,ny,nz,hex,hey,hez,igx,igy,igz
cc        real(8)    :: cnv(3)
cc        real(8)    :: vv(0:nx+1,0:ny+1,0:nz+1,3)
cc     $               ,bb(0:nx+1,0:ny+1,0:nz+1,3)
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq
cc        integer    :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg
cc
cc        real(8)    :: idhx,idhy,idhz,a(3)
cc        real(8)    :: flxip,flxim,flxjp,flxjm,flxkp,flxkm
cc
cc        real(8)    :: jacip,jacim,jacjp,jacjm,jackp,jackm
cc     .               ,jacp,jacm,jacph,jacmh,jach,jac0
cc     .               ,ijacip,ijacim,ijacjp,ijacjm,ijackp,ijackm
cc     .               ,ijacp,ijacm,ijacph,ijacmh,ijach,ijac0
cc
cc        real(8)    :: vxip,vxim,vxjp,vxjm,vxkp,vxkm
cc     .               ,vyip,vyim,vyjp,vyjm,vykp,vykm
cc     .               ,vzip,vzim,vzjp,vzjm,vzkp,vzkm
cc
cc        real(8)    :: bxip,bxim,bxjp,bxjm,bxkp,bxkm
cc     .               ,byip,byim,byjp,byjm,bykp,bykm
cc     .               ,bzip,bzim,bzjp,bzjm,bzkp,bzkm
cc
ccc     Begin program
cc
ccc     Defaults
cc
cc        ip = i+1
cc        im = i-1
cc        jp = j+1
cc        jm = j-1
cc        kp = k+1
cc        km = k-1
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        idhx = 0.5/grid_params%dxh(ig)
cc        idhy = 0.5/grid_params%dyh(jg)
cc        idhz = 0.5/grid_params%dzh(kg)
cc
cc        jac  = gmetric%grid(igx)%jac(i,j,k)
cc
cccc        sing_point = isSP(i,j,k,igx,igy,igz)
cc
ccc     Exceptions
cc
cc        !X
cccc        igp = ig+(ip-i)
cccc        igm = ig-(i-im)
cccc
cccc        dh(1)= (grid_params%xx(igp)-grid_params%xx(igm))
cc
cc
cc        if (hex == 1) then
cc          idhx = 1./dx(ig)
cc          im = i
cc        endif
cc
cc        jacip = face_add(ip,j,k,gmetric%grid(igx)%jac,1)
cc        jacim = face_add(im,j,k,gmetric%grid(igx)%jac,1)
cc
cc        jacjp = face_add(i,jp
cc          jacip  = gmetric%grid(igx)%jac(ip,j,k)
cc          jacim  = gmetric%grid(igx)%jac(i ,j,k)
cc          jacjp  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
cc     .                 +gmetric%grid(igx)%jac(i ,jp,k))
cc          jacjm  = 0.5*(gmetric%grid(igx)%jac(ip,jm,k)
cc     .                 +gmetric%grid(igx)%jac(i ,jm,k))
cc          jackp  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
cc     .                 +gmetric%grid(igx)%jac(i ,j,kp))
cc          jackm  = 0.5*(gmetric%grid(igx)%jac(ip,j,km)
cc     .                 +gmetric%grid(igx)%jac(i ,j,km))
cc
cc          if (isSP(i+1,j,k,igx,igy,igz)) then
cc            jacjp = SP_flsv
cc            jacjm = SP_flsv
cc            jackp = SP_flsv
cc            jackm = SP_flsv
cccc            write (*,*) 'DIAG -- curl_bxv'
cccc            jacjp = 1d0!SP_flsv
cccc            jacjm = 1d0!SP_flsv
cccc            jackp = 1d0!SP_flsv
cccc            jackm = 1d0!SP_flsv
cc          endif
cc
cc          vxip = vv(ip,j,k,1)
cc          vxim = vv(i ,j,k,1)
cc          vyip = vv(ip,j,k,2)
cc          vyim = vv(i ,j,k,2)
cc          vzip = vv(ip,j,k,3)
cc          vzim = vv(i ,j,k,3)
cc
cc          vxjp = 0.5*(vv(ip,jp,k,1)+vv(i,jp,k,1))
cc          vxjm = 0.5*(vv(ip,jm,k,1)+vv(i,jm,k,1))
cc          vyjp = 0.5*(vv(ip,jp,k,2)+vv(i,jp,k,2))
cc          vyjm = 0.5*(vv(ip,jm,k,2)+vv(i,jm,k,2))
cc          vzjp = 0.5*(vv(ip,jp,k,3)+vv(i,jp,k,3))
cc          vzjm = 0.5*(vv(ip,jm,k,3)+vv(i,jm,k,3))
cc
cc          vxkp = 0.5*(vv(ip,j,kp,1)+vv(i,j,kp,1))
cc          vxkm = 0.5*(vv(ip,j,km,1)+vv(i,j,km,1))
cc          vykp = 0.5*(vv(ip,j,kp,2)+vv(i,j,kp,2))
cc          vykm = 0.5*(vv(ip,j,km,2)+vv(i,j,km,2))
cc          vzkp = 0.5*(vv(ip,j,kp,3)+vv(i,j,kp,3))
cc          vzkm = 0.5*(vv(ip,j,km,3)+vv(i,j,km,3))
cc
cc          bxip = bb(ip,j,k,1)
cc          bxim = bb(i ,j,k,1)
cc          byip = bb(ip,j,k,2)
cc          byim = bb(i ,j,k,2)
cc          bzip = bb(ip,j,k,3)
cc          bzim = bb(i ,j,k,3)
cc
cc          bxjp = 0.5*(bb(ip,jp,k,1)+bb(i,jp,k,1))
cc          bxjm = 0.5*(bb(ip,jm,k,1)+bb(i,jm,k,1))
cc          byjp = 0.5*(bb(ip,jp,k,2)+bb(i,jp,k,2))
cc          byjm = 0.5*(bb(ip,jm,k,2)+bb(i,jm,k,2))
cc          bzjp = 0.5*(bb(ip,jp,k,3)+bb(i,jp,k,3))
cc          bzjm = 0.5*(bb(ip,jm,k,3)+bb(i,jm,k,3))
cc
cc          bxkp = 0.5*(bb(ip,j,kp,1)+bb(i,j,kp,1))
cc          bxkm = 0.5*(bb(ip,j,km,1)+bb(i,j,km,1))
cc          bykp = 0.5*(bb(ip,j,kp,2)+bb(i,j,kp,2))
cc          bykm = 0.5*(bb(ip,j,km,2)+bb(i,j,km,2))
cc          bzkp = 0.5*(bb(ip,j,kp,3)+bb(i,j,kp,3))
cc          bzkm = 0.5*(bb(ip,j,km,3)+bb(i,j,km,3))
cc
cc        case (2)
cc
cc          idhy = 1./dy(jg)
cc          jm = j
cc
cc          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,jp,k)
cc     .                 +gmetric%grid(igx)%jac(ip,j ,k))
cc          jacim  = 0.5*(gmetric%grid(igx)%jac(im,jp,k)
cc     .                 +gmetric%grid(igx)%jac(im,j ,k))
cc          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
cc          jacjm  = gmetric%grid(igx)%jac(i,j ,k)
cc          jackp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
cc     .                 +gmetric%grid(igx)%jac(i,j ,kp))
cc          jackm  = 0.5*(gmetric%grid(igx)%jac(i,jp,km)
cc     .                 +gmetric%grid(igx)%jac(i,j ,km))
cc
cc          vxip = (vv(ip,j,k,1)+vv(ip,jp,k,1))*0.5
cc          vxim = (vv(im,j,k,1)+vv(im,jp,k,1))*0.5
cc          vyip = (vv(ip,j,k,2)+vv(ip,jp,k,2))*0.5
cc          vyim = (vv(im,j,k,2)+vv(im,jp,k,2))*0.5
cc          vzip = (vv(ip,j,k,3)+vv(ip,jp,k,3))*0.5
cc          vzim = (vv(im,j,k,3)+vv(im,jp,k,3))*0.5
cc
cc          vxjp = vv(i,jp,k,1)
cc          vxjm = vv(i,j ,k,1)
cc          vyjp = vv(i,jp,k,2)
cc          vyjm = vv(i,j ,k,2)
cc          vzjp = vv(i,jp,k,3)
cc          vzjm = vv(i,j ,k,3)
cc
cc          vxkp = (vv(i,j,kp,1)+vv(i,jp,kp,1))*0.5
cc          vxkm = (vv(i,j,km,1)+vv(i,jp,km,1))*0.5
cc          vykp = (vv(i,j,kp,2)+vv(i,jp,kp,2))*0.5
cc          vykm = (vv(i,j,km,2)+vv(i,jp,km,2))*0.5
cc          vzkp = (vv(i,j,kp,3)+vv(i,jp,kp,3))*0.5
cc          vzkm = (vv(i,j,km,3)+vv(i,jp,km,3))*0.5
cc
cc          bxip = (bb(ip,j,k,1)+bb(ip,jp,k,1))*0.5
cc          bxim = (bb(im,j,k,1)+bb(im,jp,k,1))*0.5
cc          byip = (bb(ip,j,k,2)+bb(ip,jp,k,2))*0.5
cc          byim = (bb(im,j,k,2)+bb(im,jp,k,2))*0.5
cc          bzip = (bb(ip,j,k,3)+bb(ip,jp,k,3))*0.5
cc          bzim = (bb(im,j,k,3)+bb(im,jp,k,3))*0.5
cc
cc          bxjp = bb(i,jp,k,1)
cc          bxjm = bb(i,j ,k,1)
cc          byjp = bb(i,jp,k,2)
cc          byjm = bb(i,j ,k,2)
cc          bzjp = bb(i,jp,k,3)
cc          bzjm = bb(i,j ,k,3)
cc
cc          bxkp = (bb(i,j,kp,1)+bb(i,jp,kp,1))*0.5
cc          bxkm = (bb(i,j,km,1)+bb(i,jp,km,1))*0.5
cc          bykp = (bb(i,j,kp,2)+bb(i,jp,kp,2))*0.5
cc          bykm = (bb(i,j,km,2)+bb(i,jp,km,2))*0.5
cc          bzkp = (bb(i,j,kp,3)+bb(i,jp,kp,3))*0.5
cc          bzkm = (bb(i,j,km,3)+bb(i,jp,km,3))*0.5
cc
cc        case (3)
cc          idhz = 1./dz(kg)
cc          km = k
cc
cc          jacip  = 0.5*(gmetric%grid(igx)%jac(ip,j,kp)
cc     .                 +gmetric%grid(igx)%jac(ip,j,k ))
cc          jacim  = 0.5*(gmetric%grid(igx)%jac(im,j,kp)
cc     .                 +gmetric%grid(igx)%jac(im,j,k ))
cc          jacjp  = 0.5*(gmetric%grid(igx)%jac(i,jp,kp)
cc     .                 +gmetric%grid(igx)%jac(i,jp,k ))
cc          jacjm  = 0.5*(gmetric%grid(igx)%jac(i,jm,kp)
cc     .                 +gmetric%grid(igx)%jac(i,jm,k ))
cc          jackp  = gmetric%grid(igx)%jac(i,j,kp)
cc          jackm  = gmetric%grid(igx)%jac(i,j,km)
cc
cc          vxip = (vv(ip,j,k,1)+vv(ip,j,kp,1))*0.5
cc          vxim = (vv(im,j,k,1)+vv(im,j,kp,1))*0.5
cc          vyip = (vv(ip,j,k,2)+vv(ip,j,kp,2))*0.5
cc          vyim = (vv(im,j,k,2)+vv(im,j,kp,2))*0.5
cc          vzip = (vv(ip,j,k,3)+vv(ip,j,kp,3))*0.5
cc          vzim = (vv(im,j,k,3)+vv(im,j,kp,3))*0.5
cc
cc          vxjp = (vv(i,jp,k,1)+vv(i,jp,kp,1))*0.5
cc          vxjm = (vv(i,jm,k,1)+vv(i,jm,kp,1))*0.5
cc          vyjp = (vv(i,jp,k,2)+vv(i,jp,kp,2))*0.5
cc          vyjm = (vv(i,jm,k,2)+vv(i,jm,kp,2))*0.5
cc          vzjp = (vv(i,jp,k,3)+vv(i,jp,kp,3))*0.5
cc          vzjm = (vv(i,jm,k,3)+vv(i,jm,kp,3))*0.5
cc
cc          vxkp = vv(i,j,kp,1)
cc          vxkm = vv(i,j,k ,1)
cc          vykp = vv(i,j,kp,2)
cc          vykm = vv(i,j,k ,2)
cc          vzkp = vv(i,j,kp,3)
cc          vzkm = vv(i,j,k ,3)
cc
cc          bxip = (bb(ip,j,k,1)+bb(ip,j,kp,1))*0.5
cc          bxim = (bb(im,j,k,1)+bb(im,j,kp,1))*0.5
cc          byip = (bb(ip,j,k,2)+bb(ip,j,kp,2))*0.5
cc          byim = (bb(im,j,k,2)+bb(im,j,kp,2))*0.5
cc          bzip = (bb(ip,j,k,3)+bb(ip,j,kp,3))*0.5
cc          bzim = (bb(im,j,k,3)+bb(im,j,kp,3))*0.5
cc
cc          bxjp = (bb(i,jp,k,1)+bb(i,jp,kp,1))*0.5
cc          bxjm = (bb(i,jm,k,1)+bb(i,jm,kp,1))*0.5
cc          byjp = (bb(i,jp,k,2)+bb(i,jp,kp,2))*0.5
cc          byjm = (bb(i,jm,k,2)+bb(i,jm,kp,2))*0.5
cc          bzjp = (bb(i,jp,k,3)+bb(i,jp,kp,3))*0.5
cc          bzjm = (bb(i,jm,k,3)+bb(i,jm,kp,3))*0.5
cc
cc          bxkp = bb(i,j,kp,1)
cc          bxkm = bb(i,j,k ,1)
cc          bykp = bb(i,j,kp,2)
cc          bykm = bb(i,j,k ,2)
cc          bzkp = bb(i,j,kp,3)
cc          bzkm = bb(i,j,k ,3)
cc
cc        case default
cc          
cc          jacip  = gmetric%grid(igx)%jac(ip,j,k)
cc          jacim  = gmetric%grid(igx)%jac(im,j,k)
cc          jacjp  = gmetric%grid(igx)%jac(i,jp,k)
cc          jacjm  = gmetric%grid(igx)%jac(i,jm,k)
cc          jackp  = gmetric%grid(igx)%jac(i,j,kp)
cc          jackm  = gmetric%grid(igx)%jac(i,j,km)
cc
cc          vxip = vv(ip,j,k,1)
cc          vxim = vv(im,j,k,1)
cc          vyip = vv(ip,j,k,2)
cc          vyim = vv(im,j,k,2)
cc          vzip = vv(ip,j,k,3)
cc          vzim = vv(im,j,k,3)
cc
cc          vxjp = vv(i,jp,k,1)
cc          vxjm = vv(i,jm,k,1)
cc          vyjp = vv(i,jp,k,2)
cc          vyjm = vv(i,jm,k,2)
cc          vzjp = vv(i,jp,k,3)
cc          vzjm = vv(i,jm,k,3)
cc
cc          vxkp = vv(i,j,kp,1)
cc          vxkm = vv(i,j,km,1)
cc          vykp = vv(i,j,kp,2)
cc          vykm = vv(i,j,km,2)
cc          vzkp = vv(i,j,kp,3)
cc          vzkm = vv(i,j,km,3)
cc
cc          bxip = bb(ip,j,k,1)
cc          bxim = bb(im,j,k,1)
cc          byip = bb(ip,j,k,2)
cc          byim = bb(im,j,k,2)
cc          bzip = bb(ip,j,k,3)
cc          bzim = bb(im,j,k,3)
cc
cc          bxjp = bb(i,jp,k,1)
cc          bxjm = bb(i,jm,k,1)
cc          byjp = bb(i,jp,k,2)
cc          byjm = bb(i,jm,k,2)
cc          bzjp = bb(i,jp,k,3)
cc          bzjm = bb(i,jm,k,3)
cc
cc          bxkp = bb(i,j,kp,1)
cc          bxkm = bb(i,j,km,1)
cc          bykp = bb(i,j,kp,2)
cc          bykm = bb(i,j,km,2)
cc          bzkp = bb(i,j,kp,3)
cc          bzkm = bb(i,j,km,3)
cc
cc        end select
cc
ccc     Components
cc
cc        ijacip  = 1d0/jacip
cc        ijacim  = 1d0/jacim
cc        ijacjp  = 1d0/jacjp
cc        ijacjm  = 1d0/jacjm
cc        ijackp  = 1d0/jackp
cc        ijackm  = 1d0/jackm
cc
cccc        write (*,*) 'DIAG -- curl_bxv'
cccc        ijacip  = 1d0
cccc        ijacim  = 1d0
cccc        ijacjp  = 1d0
cccc        ijacjm  = 1d0
cccc        ijackp  = 1d0
cccc        ijackm  = 1d0
cc
cc        !component 1
cc
cc        flxjp = ( vyjp*bxjp-vxjp*byjp )*ijacjp
cc        flxjm = ( vyjm*bxjm-vxjm*byjm )*ijacjm
cc
cc        flxkp = ( vzkp*bxkp-vxkp*bzkp )*ijackp
cc        flxkm = ( vzkm*bxkm-vxkm*bzkm )*ijackm
cc
cc        cnv(1) =  (flxjp-flxjm)*idhy
cc     .           +(flxkp-flxkm)*idhz
cc
cc        !component 2
cc
cc        flxip = ( vxip*byip-vyip*bxip )*ijacip
cc        flxim = ( vxim*byim-vyim*bxim )*ijacim
cc
cc        flxkp = ( vzkp*bykp-vykp*bzkp )*ijackp
cc        flxkm = ( vzkm*bykm-vykm*bzkm )*ijackm
cc
cc        cnv(2) =  (flxip-flxim)*idhx
cc     .           +(flxkp-flxkm)*idhz
cc
cc        !component 3
cc
cc        flxip = ( vxip*bzip-vzip*bxip )*ijacip
cc        flxim = ( vxim*bzim-vzim*bxim )*ijacim
cc
cc        flxjp = ( vyjp*bzjp-vzjp*byjp )*ijacjp
cc        flxjm = ( vyjm*bzjm-vzjm*byjm )*ijacjm
cc
cc        cnv(3) =  (flxip-flxim)*idhx
cc     .           +(flxjp-flxjm)*idhy
cc
cc      contains
cc
cc      function face_add(i,j,k,array,dir) result (sum)
ccc     -----------------------------------------------------------------
ccc     This function adds contributions of 'array' at a given face from 
ccc     different points in the 27-point stencil, depending on the values 
ccc     of hex, hey, hez. The face is determined by the direction 'dir'
ccc     (=1,2,3) and the corresponding index (i,j,k).
ccc     -----------------------------------------------------------------
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,dir
cc        real(8)    :: sum,array(0:nx+1,0:ny+1,0:nz+1)
cc
ccc     Local variables
cc
ccc     Begin program
cc
cc        sum = array(i,j,k)
cc
cc        select case(dir)
cc        case(1)
cc
cc          if (hey == 1) then
cc            sum = sum + array(i,jp,k)
cc            coeff = 0.5
cc          endif
cc
cc          if (hez == 1) then
cc            sum = sum + array(i,j,kp)
cc            coeff = 0.5
cc          endif
cc
cc          if (hey == 1 .and. hez == 1) then
cc            sum = sum + array(i,jp,kp)
cc            coeff = 0.25
cc          endif
cc
cc          sum = coeff*sum
cc
cc        case(2)
cc
cc          if (hex == 1) then
cc            sum = sum + array(ip,j,k)
cc            coeff = 0.5
cc          endif
cc
cc          if (hez == 1) then
cc            sum = sum + array(i,j,kp)
cc            coeff = 0.5
cc          endif
cc
cc          if (hex == 1 .and. hez == 1) then
cc            sum = sum + array(ip,j,kp)
cc            coeff = 0.25
cc          endif
cc
cc          sum = coeff*sum
cc
cc        case(3)
cc
cc          if (hex == 1) then
cc            sum = sum + array(ip,j,k)
cc            coeff = 0.5
cc          endif
cc
cc          if (hey == 1) then
cc            sum = sum + array(i,jp,k)
cc            coeff = 0.5
cc          endif
cc
cc          if (hex == 1 .and. hey == 1) then
cc            sum = sum + array(ip,jp,k)
cc            coeff = 0.25
cc          endif
cc
cc          sum = coeff*sum
cc
cc        end select
cc
cc      end function face_add
cc
cc      end function curl_bxv2

ccc     btensor_x
ccc     #############################################################
cc      subroutine btensor_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t11,t12,t13,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t11-t13 for EOM
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t11,t12,t13
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,ip
cc        real(8)    :: x,y,z
cc        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp
cc
ccc     Begin program
cc
cc        ip = i+1
cc        if (flag == 0) ip = i
cccc        if (flag == 0) then
cccc          t11 = 0d0
cccc          t12 = 0d0
cccc          t13 = 0d0
cccc          return
cccc        endif
cc
cc        jacp = gmetric%grid(igx)%jac(ip,j,k)
cc        jac0 = gmetric%grid(igx)%jac(i ,j,k)
cc        jac  = 0.5*(jacp+jac0)
cc
cccc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cc
cc        ijac  = 1d0/jac
cc        ijac0 = 1d0/jac0
cc        ijacp = 1d0/jacp
cc
cccc        if (isSP(i+1,j,k,igx,igy,igz)) ijac = 1d0/SP_flsv
cccc        if (isSP(i+1,j,k,igx,igy,igz)) jac = 0d0
cc
cc        t11 = 0d0
cc
cc        t12 = 0.5*( vec1(i ,j,k,1)*vec2(i ,j,k,2)*ijac0
cc     .             +vec1(ip,j,k,1)*vec2(ip,j,k,2)*ijacp)*jac
cc     .       -0.5*( vec1(i ,j,k,2)*vec2(i ,j,k,1)*ijac0
cc     .             +vec1(ip,j,k,2)*vec2(ip,j,k,1)*ijacp)*jac
cc
cc        t13 = 0.5*( vec1(i ,j,k,1)*vec2(i ,j,k,3)*ijac0
cc     .             +vec1(ip,j,k,1)*vec2(ip,j,k,3)*ijacp)*jac
cc     .       -0.5*( vec1(i ,j,k,3)*vec2(i ,j,k,1)*ijac0
cc     .             +vec1(ip,j,k,3)*vec2(ip,j,k,1)*ijacp)*jac
cc
cc        if (flag /= 0) then
cc          t11 = t11*ijac
cc          if (.not.alt_eom) t12 = t12*ijac
cc          t13 = t13*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine btensor_x
cc
ccc     btensor_y
ccc     #############################################################
cc      subroutine btensor_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t21,t22,t23,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t11-t13 for EOM
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t21,t22,t23
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,jp
cc        real(8)    :: x,y,z
cc        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp
cc
ccc     Begin program
cc
cc        jp = j+1
cc        if (flag == 0) jp = j
cccc        if (flag == 0) then
cccc          t21 = 0d0
cccc          t22 = 0d0
cccc          t23 = 0d0
cccc          return
cccc        endif
cc
cc        jacp = gmetric%grid(igx)%jac(i,jp,k)
cc        jac0 = gmetric%grid(igx)%jac(i,j ,k)
cc        jac  = 0.5*(jacp+jac0)
cc
cc        ijac  = 1d0/jac
cc        ijacp = 1d0/jacp
cc        ijac0 = 1d0/jac0
cc
cc        t22 = 0d0
cc
cc        t21 = 0.5*( vec1(i,j ,k,2)*vec2(i,j ,k,1)*ijac0
cc     .             +vec1(i,jp,k,2)*vec2(i,jp,k,1)*ijacp)*jac
cc     .       -0.5*( vec1(i,j ,k,1)*vec2(i,j ,k,2)*ijac0
cc     .             +vec1(i,jp,k,1)*vec2(i,jp,k,2)*ijacp)*jac
cc
cc        t23 = 0.5*( vec1(i,j ,k,2)*vec2(i,j ,k,3)*ijac0
cc     .             +vec1(i,jp,k,2)*vec2(i,jp,k,3)*ijacp)*jac
cc     .       -0.5*( vec1(i,j ,k,3)*vec2(i,j ,k,2)*ijac0
cc     .             +vec1(i,jp,k,3)*vec2(i,jp,k,2)*ijacp)*jac
cc
cc        if (flag /= 0) then
cc          t21 = t21*ijac
cc          if (.not.alt_eom) t22 = t22*ijac
cc          t23 = t23*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine btensor_y
cc
ccc     btensor_z
ccc     #############################################################
cc      subroutine btensor_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t31,t32,t33,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t11-t13 for EOM
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t31,t32,t33
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,kp
cc        real(8)    :: x,y,z
cc        real(8)    :: jac,jac0,jacp,ijac,ijacp,ijac0
cc
ccc     Begin program
cc
cc        kp = k+1
cc        if (flag == 0) kp = k
cccc        if (flag == 0) then
cccc          t31 = 0d0
cccc          t32 = 0d0
cccc          t33 = 0d0
cccc          return
cccc        endif
cc
cc        jacp = gmetric%grid(igx)%jac(i,j,kp)
cc        jac0 = gmetric%grid(igx)%jac(i,j,k )
cc        jac  = 0.5*(jacp+jac0)
cc
cc        ijac  = 1d0/jac
cc        ijacp = 1d0/jacp
cc        ijac0 = 1d0/jac0
cc
cc        t33 = 0d0
cc
cc        t31 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,k ,1)*ijac0
cc     .             +vec1(i,j,kp,3)*vec2(i,j,kp,1)*ijacp)*jac
cc     .       -0.5*( vec1(i,j,k ,1)*vec2(i,j,k ,3)*ijac0
cc     .             +vec1(i,j,kp,1)*vec2(i,j,kp,3)*ijacp)*jac
cc
cc        t32 = 0.5*( vec1(i,j,k ,3)*vec2(i,j,k ,2)*ijac0
cc     .             +vec1(i,j,kp,3)*vec2(i,j,kp,2)*ijacp)*jac
cc     .       -0.5*( vec1(i,j,k ,2)*vec2(i,j,k ,3)*ijac0
cc     .             +vec1(i,j,kp,2)*vec2(i,j,kp,3)*ijacp)*jac
cc
cc        if (flag /= 0) then
cc          t31 = t31*ijac
cc          if (.not.alt_eom) t32 = t32*ijac
cc          t33 = t33*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine btensor_z

ccc     lf_x
ccc     #############################################################
cc      subroutine lf_x(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .               ,t11,t12,t13,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t11-t13 for linearized Lorentz
ccc     force term in conservative form.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
cc        real(8)    :: t11,t12,t13
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ip,igrid
cc        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,gsuper(3,3)
cc     .               ,scalar_prod,acnv(3),acnvp(3),bcov(3)
cc     .               ,bcovp(3),bcnv(3),bcnvp(3)
cc
ccc     Begin program
cc
cc        igrid = igx
cc
cc        ip = i+1
cc        if (flag == 0) ip = i
cc
cc        jac    = 0.5*(gmetric%grid(igrid)%jac (ip,j,k)
cc     .               +gmetric%grid(igrid)%jac (i ,j,k))
cc        gsuper = 0.5*(gmetric%grid(igrid)%gsup(ip,j,k,:,:)
cc     .               +gmetric%grid(igrid)%gsup(i ,j,k,:,:))
cc
cc        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
cc     .      .and. bcSP()
cc     .      .and. flag /= 0           ) then
cc          jacp = gmetric%grid(igrid)%jac(ip,j,k)
cc          jac0 = gmetric%grid(igrid)%jac(i ,j,k)
cc        else
cc          jacp = jac
cc          jac0 = jac
cc        endif
cc
cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cc        ijac  = 1d0/jac
cc        ijacp = 1d0/jacp
cc        ijac0 = 1d0/jac0
cc
cccc        acnv(1) = 0.5*(vec1(ip,j,k,1)*ijacp+vec1(i,j,k,1)*ijac0)*jac
cccc        acnv(2) = 0.5*(vec1(ip,j,k,2)     +vec1(i,j,k,2))
cccc        acnv(3) = 0.5*(vec1(ip,j,k,3)*ijacp+vec1(i,j,k,3)*ijac0)*jac
cccc
cccc        bcnv(1) = 0.5*(vec2(ip,j,k,1)*ijacp+vec2(i,j,k,1)*ijac0)*jac
cccc        bcnv(2) = 0.5*(vec2(ip,j,k,2)     +vec2(i,j,k,2))
cccc        bcnv(3) = 0.5*(vec2(ip,j,k,3)*ijacp+vec2(i,j,k,3)*ijac0)*jac
cccc
cccc        if (flag /= 0) then
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=1)
cccc        else
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=0)
cccc        endif
cccc
cccc        scalar_prod = dot_product(acnv,bcov)
cccc
cccc        t11 =( acnv(1)*bcnv(1)
cccc     .        +acnv(1)*bcnv(1)
cccc     .        -gsuper(1,1)*scalar_prod )
cccc
cccc        t12 =( acnv(1)*bcnv(2)
cccc     .        +acnv(2)*bcnv(1)
cccc     .        -gsuper(1,2)*scalar_prod )
cccc
cccc        t13 =( acnv(1)*bcnv(3)
cccc     .        +acnv(3)*bcnv(1)
cccc     .        -gsuper(1,3)*scalar_prod )
cc
cc        acnv(1) = vec1(i,j,k,1)*ijac0*jac
cc        acnv(2) = vec1(i,j,k,2)
cc        acnv(3) = vec1(i,j,k,3)*ijac0*jac
cc
cc        acnvp(1) = vec1(ip,j,k,1)*ijacp*jac
cc        acnvp(2) = vec1(ip,j,k,2)
cc        acnvp(3) = vec1(ip,j,k,3)*ijacp*jac
cc
cccc        bcnv(1) = vec2(i,j,k,1)*ijac0*jac
cccc        bcnv(2) = vec2(i,j,k,2)
cccc        bcnv(3) = vec2(i,j,k,3)*ijac0*jac
cccc
cccc        bcnvp(1) = vec2(ip,j,k,1)*ijacp*jac
cccc        bcnvp(2) = vec2(ip,j,k,2)
cccc        bcnvp(3) = vec2(ip,j,k,3)*ijacp*jac
cc
cc        bcnv(1) = vec2(i,j,k,1)
cc        bcnv(2) = vec2(i,j,k,2)
cc        bcnv(3) = vec2(i,j,k,3)
cc
cc        bcnvp(1) = vec2(ip,j,k,1)
cc        bcnvp(2) = vec2(ip,j,k,2)
cc        bcnvp(3) = vec2(ip,j,k,3)
cc
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        call transformFromCurvToCurv(ip,j,k,igx,igy,igz
cc     .                        ,bcovp(1),bcovp(2),bcovp(3)
cc     .                        ,bcnvp(1),bcnvp(2),bcnvp(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        bcnv(1) = bcnv(1)*ijac0*jac
cc        bcnv(3) = bcnv(3)*ijac0*jac
cc
cc        bcnvp(1) = bcnvp(1)*ijacp*jac
cc        bcnvp(3) = bcnvp(3)*ijacp*jac
cc
cc        bcov (2) = bcov (2)*ijac0*jac
cc        bcovp(2) = bcovp(2)*ijacp*jac
cc
cc        scalar_prod = dot_product(acnvp,bcov)
cc     .               +dot_product(acnv,bcovp)
cc
cc        t11 =0.5*( 2*acnvp(1)*bcnv(1)
cc     .            +2*acnv(1)*bcnvp(1)
cc     .            -gsuper(1,1)*scalar_prod )
cc
cc        t12 =0.5*( acnvp(1)*bcnv(2) + acnv(1)*bcnvp(2)
cc     .            +acnvp(2)*bcnv(1) + acnv(2)*bcnvp(1)
cc     .            -gsuper(1,2)*scalar_prod )
cc
cc        t13 =0.5*( acnvp(1)*bcnv(3) + acnv(1)*bcnvp(3)
cc     .            +acnvp(3)*bcnv(1) + acnv(3)*bcnvp(1)
cc     .            -gsuper(1,3)*scalar_prod )
cc
cccc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cccc        ijac  = 1d0/jac
cc
cc        if (flag /= 0) then
cc          t11 = t11*ijac
cc          if (.not.alt_eom) t12 = t12*ijac
cc          t13 = t13*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine lf_x
cc
ccc     lf_y
ccc     #############################################################
cc      subroutine lf_y(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .               ,t21,t22,t23,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t21-t23 for linearized Lorentz
ccc     force term in conservative form.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
cc        real(8)    :: t21,t22,t23
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: jp,igrid
cc        real(8)    :: jac,ijac,gsuper(3,3),scalar_prod
cc     .               ,acnv(3),acnvp(3),bcov(3),bcovp(3),bcnv(3),bcnvp(3)
cc
ccc     Begin program
cc
cc        igrid = igx
cc
cc        jp = j+1
cc        if (flag == 0) jp = j
cc
cc        jac    = 0.5*(gmetric%grid(igrid)%jac (i,jp,k)
cc     .               +gmetric%grid(igrid)%jac (i,j ,k))
cc        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,jp,k,:,:)
cc     .               +gmetric%grid(igrid)%gsup(i,j ,k,:,:))
cc
cc        ijac = 1d0/jac
cc
cccc        acnv(1) = 0.5*(vec1(i,jp,k,1)+vec1(i,j,k,1))
cccc        acnv(2) = 0.5*(vec1(i,jp,k,2)+vec1(i,j,k,2))
cccc        acnv(3) = 0.5*(vec1(i,jp,k,3)+vec1(i,j,k,3))
cccc                                
cccc        bcnv(1) = 0.5*(vec2(i,jp,k,1)+vec2(i,j,k,1))
cccc        bcnv(2) = 0.5*(vec2(i,jp,k,2)+vec2(i,j,k,2))
cccc        bcnv(3) = 0.5*(vec2(i,jp,k,3)+vec2(i,j,k,3))
cccc
cccc        if (flag /= 0) then
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=2)
cccc        else
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=0)
cccc        endif
cccc
cccc
cccc        scalar_prod = dot_product(acnv,bcov)
cccc
cccc        t21 =( acnv(2)*bcnv(1)
cccc     .        +acnv(1)*bcnv(2)
cccc     .        -gsuper(2,1)*scalar_prod )
cccc
cccc        t22 =( acnv(2)*bcnv(2)
cccc     .        +acnv(2)*bcnv(2)
cccc     .        -gsuper(2,2)*scalar_prod )
cccc
cccc        t23 =( acnv(2)*bcnv(3)
cccc     .        +acnv(3)*bcnv(2)
cccc     .        -gsuper(2,3)*scalar_prod )
cc
cc        acnv(1) = vec1(i,j,k,1)
cc        acnv(2) = vec1(i,j,k,2)
cc        acnv(3) = vec1(i,j,k,3)
cc
cc        acnvp(1) = vec1(i,jp,k,1)
cc        acnvp(2) = vec1(i,jp,k,2)
cc        acnvp(3) = vec1(i,jp,k,3)
cc
cc        bcnv(1) = vec2(i,j,k,1)
cc        bcnv(2) = vec2(i,j,k,2)
cc        bcnv(3) = vec2(i,j,k,3)
cc
cc        bcnvp(1) = vec2(i,jp,k,1)
cc        bcnvp(2) = vec2(i,jp,k,2)
cc        bcnvp(3) = vec2(i,jp,k,3)
cc
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        call transformFromCurvToCurv(i,jp,k,igx,igy,igz
cc     .                        ,bcovp(1),bcovp(2),bcovp(3)
cc     .                        ,bcnvp(1),bcnvp(2),bcnvp(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        scalar_prod = dot_product(acnvp,bcov)
cc     .               +dot_product(acnv,bcovp)
cc
cc        t21 =0.5*( acnvp(2)*bcnv(1) + acnv(2)*bcnvp(1)
cc     .            +acnvp(1)*bcnv(2) + acnv(1)*bcnvp(2)
cc     .            -gsuper(2,1)*scalar_prod )
cc
cc        t22 =0.5*( 2*acnvp(2)*bcnv(2) + 2*acnv(2)*bcnvp(2)
cc     .            -gsuper(2,2)*scalar_prod )
cc
cc        t23 =0.5*( acnvp(2)*bcnv(3) + acnv(2)*bcnvp(3)
cc     .            +acnvp(3)*bcnv(2) + acnv(3)*bcnvp(2)
cc     .            -gsuper(2,3)*scalar_prod )
cc
cc        if (flag /= 0) then
cc          t21 = t21*ijac
cc          if (.not.alt_eom) t22 = t22*ijac
cc          t23 = t23*ijac
cc        endif
cc
cc
ccc     End program
cc
cc      end subroutine lf_y
cc
ccc     lf_z
ccc     #############################################################
cc      subroutine lf_z(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
cc     .                   ,t31,t32,t33,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t31-t33 for linearized Lorentz
ccc     force term in conservative form.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,nxx,nyy,nzz,igx,igy,igz,flag
cc        real(8)    :: t31,t32,t33
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: kp,igrid
cc        real(8)    :: jac,ijac,gsuper(3,3),scalar_prod
cc     .               ,acnv(3),acnvp(3),bcov(3),bcovp(3),bcnv(3),bcnvp(3)
cc
ccc     Begin program
cc
cc        igrid = igx
cc
cc        kp=k+1
cc        if (flag == 0) kp = k
cc
cc        jac    = 0.5*(gmetric%grid(igrid)%jac (i,j,kp)
cc     .               +gmetric%grid(igrid)%jac (i,j,k ))
cc        gsuper = 0.5*(gmetric%grid(igrid)%gsup(i,j,kp,:,:)
cc     .               +gmetric%grid(igrid)%gsup(i,j,k ,:,:))
cc
cc        ijac = 1d0/jac
cc
cccc        acnv(1) = 0.5*(vec1(i,j,kp,1)+vec1(i,j,k,1))
cccc        acnv(2) = 0.5*(vec1(i,j,kp,2)+vec1(i,j,k,2))
cccc        acnv(3) = 0.5*(vec1(i,j,kp,3)+vec1(i,j,k,3))
cccc                                  
cccc        bcnv(1) = 0.5*(vec2(i,j,kp,1)+vec2(i,j,k,1))
cccc        bcnv(2) = 0.5*(vec2(i,j,kp,2)+vec2(i,j,k,2))
cccc        bcnv(3) = 0.5*(vec2(i,j,kp,3)+vec2(i,j,k,3))
cccc
cccc        if (flag /= 0) then
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=3)
cccc        else
cccc          call transformFromCurvToCurv(i,j,k,igx,igy,igz
cccc     .                        ,bcov(1),bcov(2),bcov(3)
cccc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cccc     .                        ,.false.,half_elem=0)
cccc        endif
cccc
cccc        scalar_prod = dot_product(acnv,bcov)
cccc
cccc        t31 =( acnv(3)*bcnv(1)
cccc     .        +acnv(1)*bcnv(3)
cccc     .        -gsuper(3,1)*scalar_prod )
cccc
cccc        t32 =( acnv(3)*bcnv(2)
cccc     .        +acnv(2)*bcnv(3)
cccc     .        -gsuper(3,2)*scalar_prod )
cccc
cccc        t33 =( acnv(3)*bcnv(3)
cccc     .        +acnv(3)*bcnv(3)
cccc     .        -gsuper(3,3)*scalar_prod )
cc
cc        acnv(1) = vec1(i,j,k,1)
cc        acnv(2) = vec1(i,j,k,2)
cc        acnv(3) = vec1(i,j,k,3)
cc
cc        acnvp(1) = vec1(i,j,kp,1)
cc        acnvp(2) = vec1(i,j,kp,2)
cc        acnvp(3) = vec1(i,j,kp,3)
cc
cc        bcnv(1) = vec2(i,j,k,1)
cc        bcnv(2) = vec2(i,j,k,2)
cc        bcnv(3) = vec2(i,j,k,3)
cc
cc        bcnvp(1) = vec2(i,j,kp,1)
cc        bcnvp(2) = vec2(i,j,kp,2)
cc        bcnvp(3) = vec2(i,j,kp,3)
cc
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                        ,bcov(1),bcov(2),bcov(3)
cc     .                        ,bcnv(1),bcnv(2),bcnv(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        call transformFromCurvToCurv(i,j,kp,igx,igy,igz
cc     .                        ,bcovp(1),bcovp(2),bcovp(3)
cc     .                        ,bcnvp(1),bcnvp(2),bcnvp(3)
cc     .                        ,.false.,half_elem=0)
cc
cc        scalar_prod = dot_product(acnvp,bcov)
cc     .               +dot_product(acnv,bcovp)
cc
cc        t31 =0.5*( acnvp(3)*bcnv(1) + acnv(3)*bcnvp(1)
cc     .            +acnvp(1)*bcnv(3) + acnv(1)*bcnvp(3)
cc     .            -gsuper(3,1)*scalar_prod )
cc
cc        t32 =0.5*( acnvp(3)*bcnv(2) + acnv(3)*bcnvp(2)
cc     .            +acnvp(2)*bcnv(3) + acnv(2)*bcnvp(3)
cc     .            -gsuper(3,2)*scalar_prod )
cc
cc        t33 =0.5*( 2*acnvp(3)*bcnv(3)
cc     .           + 2*acnv(3)*bcnvp(3)
cc     .            -gsuper(3,3)*scalar_prod )
cc
cc        if (flag /= 0) then
cc          t31 = t31*ijac
cc          if (.not.alt_eom) t32 = t32*ijac
cc          t33 = t33*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine lf_z

ccc     par_diff_2nd
ccc     ################################################################
cc      function par_diff_2nd(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
cc     .                     ,limited,vol) result(par2)
ccc     ----------------------------------------------------------------
ccc     Calculates div.(BB.grad) in 2D with 2nd order accuracy. Vector
ccc     B *must* be unitary.
ccc     ----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer  :: i,j,k,nx,ny,nz,igx,igy,igz
cc      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1),par2
cc      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
cc     .           ,By(0:nx+1,0:ny+1,0:nz+1)
cc     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)
cc
cc      logical,optional :: limited,vol
cc
ccc     Local variables
cc
cc      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
cc     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
cc     .           ,Fluxz(i-1:i,j-1:j,k-1:k)
cc
cc      real(8)  :: bx2,by2,bxby,dtx,dty,a1,a2,sgn,bb,flx0,flxp
cc      integer  :: ii,jj,kk,ig,jg,kg
cc      logical  :: limit,symm,vol_wgt
cc
cc      real(8),dimension(:,:,:),pointer :: jac
cc
ccc     Begin program
cc
cc      Fluxx=0
cc      Fluxy=0
cc      Fluxz=0
cc
cc      if (PRESENT(limited)) then
cc        limit = limited
cc      else
cc        limit = .true.
cc      endif
cc
cc      vol_wgt = .true.
cc      if (PRESENT(vol)) vol_wgt = vol
cc
cccc      if (PRESENT(symmetric)) then
cccc        symm = symmetric
cccc      else
cccc        symm = .false.
cccc      endif
cc
cc      jac => gmetric%grid(igx)%jac
cc
ccc     Symmetric discretization
cc
cccc      if (symm) then
cc
cc        do kk = k-1,k
cc          do jj = j-1,j
cc            do ii = i-1,i
cc              call getMGmap(ii,jj,kk,igx,igy,igz,ig,jg,kg)
cc
cc              bx2 = (.25*(bx(ii  ,jj  ,kk)/jac(ii  ,jj  ,kk)
cc     .                   +bx(ii+1,jj  ,kk)/jac(ii+1,jj  ,kk)
cc     .                   +bx(ii  ,jj+1,kk)/jac(ii  ,jj+1,kk)
cc     .                   +bx(ii+1,jj+1,kk)/jac(ii+1,jj+1,kk)))
cc     .             *(.25*(bx(ii  ,jj  ,kk)
cc     .                   +bx(ii+1,jj  ,kk)
cc     .                   +bx(ii  ,jj+1,kk)
cc     .                   +bx(ii+1,jj+1,kk)))
cc              bxby= (.25*(bx(ii  ,jj  ,kk)/jac(ii  ,jj  ,kk)
cc     .                   +bx(ii+1,jj  ,kk)/jac(ii+1,jj  ,kk)
cc     .                   +bx(ii  ,jj+1,kk)/jac(ii  ,jj+1,kk)
cc     .                   +bx(ii+1,jj+1,kk)/jac(ii+1,jj+1,kk)))
cc     .             *(.25*(by(ii  ,jj  ,kk)
cc     .                   +by(ii+1,jj  ,kk)
cc     .                   +by(ii  ,jj+1,kk)
cc     .                   +by(ii+1,jj+1,kk)))
cc              by2 = (.25*(by(ii  ,jj  ,kk)/jac(ii  ,jj  ,kk)
cc     .                   +by(ii+1,jj  ,kk)/jac(ii+1,jj  ,kk)
cc     .                   +by(ii  ,jj+1,kk)/jac(ii  ,jj+1,kk)
cc     .                   +by(ii+1,jj+1,kk)/jac(ii+1,jj+1,kk)))
cc     .             *(.25*(by(ii  ,jj  ,kk)
cc     .                   +by(ii+1,jj  ,kk)
cc     .                   +by(ii  ,jj+1,kk)
cc     .                   +by(ii+1,jj+1,kk)))
cc
cc              dtx = .5*(tmp(ii+1,jj  ,kk)-tmp(ii,jj  ,kk)
cc     .                 +tmp(ii+1,jj+1,kk)-tmp(ii,jj+1,kk))/dx(ig)
cc              dty = .5*(tmp(ii  ,jj+1,kk)-tmp(ii  ,jj,kk)
cc     .                 +tmp(ii+1,jj+1,kk)-tmp(ii+1,jj,kk))/dy(jg)
cc
cc              Fluxx(ii,jj,kk)= bx2*dtx + bxby*dty
cc              Fluxy(ii,jj,kk)= by2*dty + bxby*dtx
cc              Fluxz(ii,jj,kk) = 0d0
cc            enddo
cc          enddo
cc        enddo
cc
cc        !Average fluxes from vertices to faces
cc        Fluxx(:,j,k) = 0.25*(Fluxx(:,j,k  )+Fluxx(:,j-1,k  )
cc     .                      +Fluxx(:,j,k-1)+Fluxx(:,j-1,k-1))
cc
cc        Fluxy(i,:,k) = 0.25*(Fluxy(i,:,k  )+Fluxy(i-1,:,k  )
cc     .                      +Fluxy(i,:,k-1)+Fluxy(i-1,:,k-1))
cc
cc        Fluxz(i,j,:) = 0.25*(Fluxz(i,j  ,:)+Fluxz(i-1,j  ,:)
cc     .                      +Fluxz(i,j-1,:)+Fluxz(i-1,j-1,:))
cc
cc        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        par2 = (Fluxx(i,j,k)-Fluxx(i-1,j,k))/dxh(ig)
cc     .        +(Fluxy(i,j,k)-Fluxy(i,j-1,k))/dyh(jg)
cc     .        +(Fluxz(i,j,k)-Fluxz(i,j,k-1))/dzh(kg)
cc
cc        if (vol_wgt) par2=par2*gmetric%grid(igx)%dvol(i,j,k)
cc
ccccc     Uncomment this for dding straight co-directed derivatives
cccc
cccc        !Add diagonally dominant part
cccc        do i = 0,nx
cccc          call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cccc
cccc          Fluxx(i,1:ny,1:nz) = Fluxx(i,1:ny,1:nz)
cccc     .                +(0.5*(bx (i+1,1:ny,1:nz)+bx (i,1:ny,1:nz)))**2
cccc     .                     *(tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
cccc        enddo
cccc
cccc        do j = 0,ny
cccc          call getMGmap(1,j,1,igx,igy,igz,ig,jg,kg)
cccc
cccc          Fluxy(1:nx,j,1:nz) = Fluxy(1:nx,j,1:nz)
cccc     .                +(0.5*(by (1:nx,j+1,1:nz)+by (1:nx,j,1:nz)))**2
cccc     .                     *(tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
cccc        enddo
cccc
cccc        do k = 0,nz
cccc          call getMGmap(1,1,k,igx,igy,igz,ig,jg,kg)
cccc
cccc          Fluxz(1:nx,1:ny,k) = Fluxz(1:nx,1:ny,k)
cccc     .                +(0.5*(bz (1:nx,1:ny,k+1)+bz (1:nx,1:ny,k)))**2
cccc     .                     *(tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
cccc        enddo
cc
ccccc     Asymmetric discretization
cccc
cccc      else
cccc
cccc        if (limit) then !Limited differences (Sharma, Hammett, JCP 227 (2007))
cccc
cccc          !X-component cross flux
cccc          do k = 1,nz
cccc            do j = 1,ny
cccc              do i = 0,nx
cccc                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cccc
cccc                bxby = 0.5*(bx(i  ,j,k)*by(i  ,j,k)
cccc     .                     +bx(i+1,j,k)*by(i+1,j,k))
cccc
cccc                flx0 = mc(tmp(i,j  ,k)-tmp(i,j-1,k)
cccc     .                   ,tmp(i,j+1,k)-tmp(i,j  ,k))
cccc                flxp = mc(tmp(i+1,j  ,k)-tmp(i+1,j-1,k)
cccc     .                   ,tmp(i+1,j+1,k)-tmp(i+1,j  ,k))
cccc
cccc                Fluxx(i,j,k)= bxby/dyh(jg)*mc(flx0,flxp)
cccc              enddo
cccc            enddo
cccc          enddo
cccc
cccc          !Y-component cross flux 
cccc          do k = 1,nz
cccc            do j = 0,ny
cccc              do i = 1,nx
cccc
cccc                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cccc
cccc                bxby = 0.5*(bx(i,j  ,k)*by(i,j  ,k)
cccc     .                     +bx(i,j+1,k)*by(i,j+1,k))
cccc
cccc                flx0 = mc(tmp(i  ,j,k)-tmp(i-1,j,k)
cccc     .                   ,tmp(i+1,j,k)-tmp(i  ,j,k))
cccc                flxp = mc(tmp(i  ,j+1,k)-tmp(i-1,j+1,k)
cccc     .                   ,tmp(i+1,j+1,k)-tmp(i  ,j+1,k))
cccc
cccc                Fluxy(i,j,k)= bxby/dxh(ig)*mc(flx0,flxp)
cccc
cccc              enddo
cccc            enddo
cccc          enddo
cccc
cccc        else !Standard central average
cccc
cccc          !X-component cross flux
cccc          do k = 1,nz
cccc            do j = 1,ny
cccc              do i = 0,nx
cccc                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cccc
cccc                Fluxx(i,j,k)= 0.125/dyh(jg)
cccc     .             *(bx(i  ,j,k)*by(i  ,j,k)
cccc     .              +bx(i+1,j,k)*by(i+1,j,k))
cccc     .             *(tmp(i  ,j+1,k)-tmp(i  ,j-1,k)
cccc     .              +tmp(i+1,j+1,k)-tmp(i+1,j-1,k))
cccccc                Fluxx(i,j,k)= 0.25/dyh(jg)
cccccc     .             *(bx(i  ,j,k)*by(i  ,j,k)
cccccc     .                 *(tmp(i  ,j+1,k)-tmp(i  ,j-1,k))
cccccc     .              +bx(i+1,j,k)*by(i+1,j,k)
cccccc     .                 *(tmp(i+1,j+1,k)-tmp(i+1,j-1,k)))
cccc              enddo
cccc            enddo
cccc          enddo
cccc
cccc          !Y-component cross flux 
cccc          do k = 1,nz
cccc            do j = 0,ny
cccc              do i = 1,nx
cccc                call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cccc
cccc                Fluxy(i,j,k)= 0.125/dxh(ig)
cccc     .             *(bx(i,j  ,k)*by(i,j  ,k)
cccc     .              +bx(i,j+1,k)*by(i,j+1,k))
cccc     .             *(tmp(i+1,j  ,k)-tmp(i-1,j  ,k)
cccc     .              +tmp(i+1,j+1,k)-tmp(i-1,j+1,k))
cccccc                Fluxy(i,j,k)= 0.25/dxh(ig)
cccccc     .             *(bx(i,j  ,k)*by(i,j  ,k)
cccccc     .                 *(tmp(i+1,j  ,k)-tmp(i-1,j  ,k))
cccccc     .              +bx(i,j+1,k)*by(i,j+1,k)
cccccc     .                 *(tmp(i+1,j+1,k)-tmp(i-1,j+1,k)))
cccc              enddo
cccc            enddo
cccc          enddo
cccc
cccc        endif
cccc
cccc        !Add co-derivative terms: standard centered
cccc        do i = 0,nx
cccc          call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cccc
cccc          Fluxx(i,1:ny,1:nz) = Fluxx(i,1:ny,1:nz)
cccc     .                +(0.5*(bx (i+1,1:ny,1:nz)+bx (i,1:ny,1:nz)))**2
cccc     .                     *(tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
cccc        enddo
cccc
cccc        do j = 0,ny
cccc          call getMGmap(1,j,1,igx,igy,igz,ig,jg,kg)
cccc
cccc          Fluxy(1:nx,j,1:nz) = Fluxy(1:nx,j,1:nz)
cccc     .                +(0.5*(by (1:nx,j+1,1:nz)+by (1:nx,j,1:nz)))**2
cccc     .                     *(tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
cccc        enddo
cccc
cccc        do k = 0,nz
cccc          call getMGmap(1,1,k,igx,igy,igz,ig,jg,kg)
cccc
cccc          Fluxz(1:nx,1:ny,k) = Fluxz(1:nx,1:ny,k)
cccc     .                +(0.5*(bz (1:nx,1:ny,k+1)+bz (1:nx,1:ny,k)))**2
cccc     .                     *(tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
cccc        enddo
cccc
cccccc      !Add co-derivative terms: ZIP
cccccc      do i = 0,nx
cccccc        call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cccccc
cccccc        Fluxx(i,1:ny,1:nz) = Fluxx(i,1:ny,1:nz) 
cccccc     .              +bx(i+1,1:ny,1:nz)*bx(i,1:ny,1:nz)
cccccc     .                   *(tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
cccccc      enddo
cccccc
cccccc      do j = 0,ny
cccccc        call getMGmap(1,j,1,igx,igy,igz,ig,jg,kg)
cccccc
cccccc        Fluxy(1:nx,j,1:nz) = Fluxy(1:nx,j,1:nz) 
cccccc     .              +by(1:nx,j+1,1:nz)*by(1:nx,j,1:nz)
cccccc     .                   *(tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
cccccc      enddo
cccccc
cccccc      do k = 0,nz
cccccc        call getMGmap(1,1,k,igx,igy,igz,ig,jg,kg)
cccccc
cccccc        Fluxz(1:nx,1:ny,k) = Fluxz(1:nx,1:ny,k)
cccccc     .              +bz(1:nx,1:ny,k+1)*bz(1:nx,1:ny,k)
cccccc     .                   *(tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
cccccc      enddo
cccc
cccc      endif
cc
ccc     End program
cc
cc      contains
cc
ccc     minmod
ccc     #################################################################
cc      function minmod(a,b)
cc
cc      real(8) :: a,b,minmod
cc
cc      if (a*b <= 0d0) then
cc        minmod = 0d0
cc        return
cc      endif
cc
cc      if (a > 0d0) then
cc        minmod = min(a,b)
cc      else
cc        minmod = max(a,b)
cc      endif
cc
cc      end function minmod
cc
ccc     mc
ccc     #################################################################
cc      function mc(a,b)
cc
cc      real(8) :: a,b,mc
cc
cc      mc = minmod(2*minmod(a,b),0.5*(a+b))
cc
cc      end function mc
cc
cc      end function par_diff_2nd

c     par_diff
c     ################################################################
      function par_diff(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
     .                 ,par_order,symmetric,limited,vol) result(par2)
c     ----------------------------------------------------------------
c     Calculates div.(BB.grad) in 2D with 2nd order accuracy. Vector
c     B *must* be unitary.
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      integer  :: i,j,k,nx,ny,nz,igx,igy,igz,par_order
      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1),par2
      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
     .           ,By(0:nx+1,0:ny+1,0:nz+1)
     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)

      logical,optional :: symmetric,limited,vol

c     Local variables

      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
     .           ,Fluxz(i-1:i,j-1:j,k-1:k)

      real(8)  :: bx2,by2,bxby,dtx,dty,a1,a2,sgn,bb,flx0,flxp
      integer  :: ii,jj,kk,ig,jg,kg
      logical  :: limit,symm,vol_wgt

      real(8),dimension(:,:,:),pointer :: jac

c     Begin program

      Fluxx=0
      Fluxy=0
      Fluxz=0

      if (PRESENT(limited)) then
        limit = limited
      else
        limit = .true.
      endif

      vol_wgt = .true.
      if (PRESENT(vol)) vol_wgt = vol

      if (PRESENT(symmetric)) then
        symm = symmetric
      else
        symm = .false.
      endif

      jac => gmetric%grid(igx)%jac

c     Find fluxes

      select case(par_order)
      case(2) !second order
        call par_flux_2nd(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
     .                 ,Fluxx,Fluxy,Fluxz,symmetric=symm,limited=limit)
      case(4) !fourth order
        call par_flux_4th(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
     .                 ,Fluxx,Fluxy,Fluxz,limited=limit)
      case default
        write (*,*) 'Order not implemented'
        stop
      end select

c     Assemble fluxes

      call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

      par2 = (Fluxx(i,j,k)-Fluxx(i-1,j,k))/grid_params%dxh(ig)
     .      +(Fluxy(i,j,k)-Fluxy(i,j-1,k))/grid_params%dyh(jg)
     .      +(Fluxz(i,j,k)-Fluxz(i,j,k-1))/grid_params%dzh(kg)

      if (vol_wgt) par2=par2*gmetric%grid(igx)%dvol(i,j,k)

      end function par_diff

c     par_flux_2nd
c     ################################################################
      subroutine par_flux_2nd(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
     .                       ,Fluxx,Fluxy,Fluxz,symmetric,limited)
c     ----------------------------------------------------------------
c     Calculates div.(BB.grad) in 3D with 2nd order accuracy
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      integer  :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
     .           ,By(0:nx+1,0:ny+1,0:nz+1)
     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
     .           ,Fluxz(i-1:i,j-1:j,k-1:k)

      logical,optional :: limited,symmetric

c     Local variables

      real(8)  :: bx2,by2,bxby,dtx,dty,a1,a2,sgn,bb,flx0,flxp
      integer  :: ig,jg,kg,ii,jj,kk
      logical  :: limit,symm

      real(8),dimension(:,:,:),pointer :: jac

c     Begin program

      Fluxx=0
      Fluxy=0
      Fluxz=0

      if (PRESENT(limited)) then
        limit = limited
      else
        limit = .true.
      endif

      if (PRESENT(symmetric)) then
        symm = symmetric
      else
        symm = .false.
      endif

      jac => gmetric%grid(igx)%jac

c     Symmetric discretization

      if (symm) then

        do kk = k-1,k
          do jj = j-1,j
            do ii = i-1,i
              call getMGmap(ii,jj,kk,igx,igy,igz,ig,jg,kg)

              bx2 = (.25*(bx(ii  ,jj  ,kk)/jac(ii  ,jj  ,kk)
     .                   +bx(ii+1,jj  ,kk)/jac(ii+1,jj  ,kk)
     .                   +bx(ii  ,jj+1,kk)/jac(ii  ,jj+1,kk)
     .                   +bx(ii+1,jj+1,kk)/jac(ii+1,jj+1,kk)))
     .             *(.25*(bx(ii  ,jj  ,kk)
     .                   +bx(ii+1,jj  ,kk)
     .                   +bx(ii  ,jj+1,kk)
     .                   +bx(ii+1,jj+1,kk)))
              bxby= (.25*(bx(ii  ,jj  ,kk)/jac(ii  ,jj  ,kk)
     .                   +bx(ii+1,jj  ,kk)/jac(ii+1,jj  ,kk)
     .                   +bx(ii  ,jj+1,kk)/jac(ii  ,jj+1,kk)
     .                   +bx(ii+1,jj+1,kk)/jac(ii+1,jj+1,kk)))
     .             *(.25*(by(ii  ,jj  ,kk)
     .                   +by(ii+1,jj  ,kk)
     .                   +by(ii  ,jj+1,kk)
     .                   +by(ii+1,jj+1,kk)))
              by2 = (.25*(by(ii  ,jj  ,kk)/jac(ii  ,jj  ,kk)
     .                   +by(ii+1,jj  ,kk)/jac(ii+1,jj  ,kk)
     .                   +by(ii  ,jj+1,kk)/jac(ii  ,jj+1,kk)
     .                   +by(ii+1,jj+1,kk)/jac(ii+1,jj+1,kk)))
     .             *(.25*(by(ii  ,jj  ,kk)
     .                   +by(ii+1,jj  ,kk)
     .                   +by(ii  ,jj+1,kk)
     .                   +by(ii+1,jj+1,kk)))

              dtx = .5*(tmp(ii+1,jj  ,kk)-tmp(ii,jj  ,kk)
     .                 +tmp(ii+1,jj+1,kk)-tmp(ii,jj+1,kk))
     .                 /grid_params%dx(ig)
              dty = .5*(tmp(ii  ,jj+1,kk)-tmp(ii  ,jj,kk)
     .                 +tmp(ii+1,jj+1,kk)-tmp(ii+1,jj,kk))
     .                 /grid_params%dy(jg)

              Fluxx(ii,jj,kk)= bx2*dtx + bxby*dty
              Fluxy(ii,jj,kk)= by2*dty + bxby*dtx
              Fluxz(ii,jj,kk) = 0d0
            enddo
          enddo
        enddo

        !Average fluxes from vertices to faces
        Fluxx(:,j,k) = 0.25*(Fluxx(:,j,k  )+Fluxx(:,j-1,k  )
     .                      +Fluxx(:,j,k-1)+Fluxx(:,j-1,k-1))

        Fluxy(i,:,k) = 0.25*(Fluxy(i,:,k  )+Fluxy(i-1,:,k  )
     .                      +Fluxy(i,:,k-1)+Fluxy(i-1,:,k-1))

        Fluxz(i,j,:) = 0.25*(Fluxz(i,j  ,:)+Fluxz(i-1,j  ,:)
     .                      +Fluxz(i,j-1,:)+Fluxz(i-1,j-1,:))

c     Asymmetric discretization

      else

        if (limit) then !Limited differences (Sharma, Hammett, JCP 227 (2007))

          !X-component cross flux
          do ii = i-1,i
            call getMGmap(ii,j,k,igx,igy,igz,ig,jg,kg)

            bxby = 0.5*(bx(ii  ,j,k)/jac(ii  ,j,k)
     .                 *by(ii  ,j,k)
     .                 +bx(ii+1,j,k)/jac(ii+1,j,k)
     .                 *by(ii+1,j,k))

            flx0 = mc(tmp(ii  ,j  ,k)-tmp(ii  ,j-1,k)
     .               ,tmp(ii  ,j+1,k)-tmp(ii  ,j  ,k))
            flxp = mc(tmp(ii+1,j  ,k)-tmp(ii+1,j-1,k)
     .               ,tmp(ii+1,j+1,k)-tmp(ii+1,j  ,k))

            Fluxx(ii,j,k)= bxby/grid_params%dyh(jg)*mc(flx0,flxp)
          enddo

          !Y-component cross flux 
          do jj = j-1,j
            call getMGmap(i,jj,k,igx,igy,igz,ig,jg,kg)

            bxby = 0.5*(bx(i,jj  ,k)/jac(i,jj  ,k)
     .                 *by(i,jj  ,k)
     .                 +bx(i,jj+1,k)/jac(i,jj+1,k)
     .                 *by(i,jj+1,k))

            flx0 = mc(tmp(i  ,jj  ,k)-tmp(i-1,jj  ,k)
     .               ,tmp(i+1,jj  ,k)-tmp(i  ,jj  ,k))
            flxp = mc(tmp(i  ,jj+1,k)-tmp(i-1,jj+1,k)
     .               ,tmp(i+1,jj+1,k)-tmp(i  ,jj+1,k))

            Fluxy(i,jj,k)= bxby/grid_params%dxh(ig)*mc(flx0,flxp)
          enddo

        else !Standard central average

          !X-component cross flux
          do ii = i-1,i
            call getMGmap(ii,j,k,igx,igy,igz,ig,jg,kg)

            bxby = 0.5*(bx(ii  ,j,k)/jac(ii  ,j,k)*by(ii  ,j,k)
     .                 +bx(ii+1,j,k)/jac(ii+1,j,k)*by(ii+1,j,k))

            Fluxx(ii,j,k)= 0.25/grid_params%dyh(jg)*bxby
     .             *(tmp(ii  ,j+1,k)-tmp(ii  ,j-1,k)
     .              +tmp(ii+1,j+1,k)-tmp(ii+1,j-1,k))
          enddo

          !Y-component cross flux 
          do jj = j-1,j
            call getMGmap(i,jj,k,igx,igy,igz,ig,jg,kg)

            bxby = 0.5*(bx(i,jj  ,k)/jac(i,jj  ,k)*by(i,jj  ,k)
     .                 +bx(i,jj+1,k)/jac(i,jj+1,k)*by(i,jj+1,k))

            Fluxy(i,jj,k)= 0.25/grid_params%dxh(ig)*bxby
     .                   *(tmp(i+1,jj  ,k)-tmp(i-1,jj  ,k)
     .                    +tmp(i+1,jj+1,k)-tmp(i-1,jj+1,k))
          enddo

        endif

        !Add co-derivative terms: standard centered
        do ii = i-1,i
          call getMGmap(ii,j,k,igx,igy,igz,ig,jg,kg)

          bx2 = (0.5*(bx(ii+1,j,k)/jac(ii+1,j,k)
     .               +bx(ii  ,j,k)/jac(ii  ,j,k)))**2
          Fluxx(ii,j,k) = Fluxx(ii,j,k)
     .                  + bx2*(tmp(ii+1,j,k)-tmp(ii,j,k))
     .                        /grid_params%dx(ig)
        enddo

        do jj = j-1,j
          call getMGmap(i,jj,k,igx,igy,igz,ig,jg,kg)
 
          by2 = (0.5*(by(i,jj+1,k)/jac(i,jj+1,k)
     .               +by(i,jj  ,k)/jac(i,jj  ,k)))**2

          Fluxy(i,jj,k) = Fluxy(i,jj,k)
     .                  + by2*(tmp(i,jj+1,k)-tmp(i,jj,k))
     .                        /grid_params%dy(jg)
        enddo

        do kk = k-1,k
          call getMGmap(i,j,kk,igx,igy,igz,ig,jg,kg)

          Fluxz(i,j,kk) = Fluxz(i,j,kk)
     .                +(0.5*(bz(i,j,kk+1)/jac(i,j,kk+1)
     .                      +bz(i,j,kk  )/jac(i,j,kk  )))**2
     .                     *(tmp(i,j,kk+1)-tmp(i,j,kk))
     .                      /grid_params%dz(kg)
        enddo

cc      !Add co-derivative terms: ZIP
cc      do i = 0,nx
cc        call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cc
cc        Fluxx(i,1:ny,1:nz) = Fluxx(i,1:ny,1:nz) 
cc     .              +bx(i+1,1:ny,1:nz)*bx(i,1:ny,1:nz)
cc     .                   *(tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
cc      enddo
cc
cc      do j = 0,ny
cc        call getMGmap(1,j,1,igx,igy,igz,ig,jg,kg)
cc
cc        Fluxy(1:nx,j,1:nz) = Fluxy(1:nx,j,1:nz) 
cc     .              +by(1:nx,j+1,1:nz)*by(1:nx,j,1:nz)
cc     .                   *(tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
cc      enddo
cc
cc      do k = 0,nz
cc        call getMGmap(1,1,k,igx,igy,igz,ig,jg,kg)
cc
cc        Fluxz(1:nx,1:ny,k) = Fluxz(1:nx,1:ny,k)
cc     .              +bz(1:nx,1:ny,k+1)*bz(1:nx,1:ny,k)
cc     .                   *(tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
cc      enddo

      endif

c     End program

      contains

c     minmod
c     #################################################################
      function minmod(a,b)

      real(8) :: a,b,minmod

      if (a*b <= 0d0) then
        minmod = 0d0
        return
      endif

      if (a > 0d0) then
        minmod = min(a,b)
      else
        minmod = max(a,b)
      endif

      end function minmod

c     mc
c     #################################################################
      function mc(a,b)

      real(8) :: a,b,mc

      mc = minmod(2*minmod(a,b),0.5*(a+b))

      end function mc

      end subroutine par_flux_2nd

ccc     iso_flux_2nd
ccc     ################################################################
cc      subroutine iso_flux_2nd(nx,ny,nz,igx,igy,igz,tmp
cc     .                       ,Fluxx,Fluxy,Fluxz,symmetric)
ccc     ----------------------------------------------------------------
ccc     Calculates div.(grad) in 3D with 2nd order accuracy
ccc     ----------------------------------------------------------------
cc
cc      implicit none
cc
ccc     Call variables
cc
cc      integer  :: nx,ny,nz,igx,igy,igz
cc      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1)
cc      real(8)  :: Fluxx(0:nx,0:ny,0:nz)
cc     .           ,Fluxy(0:nx,0:ny,0:nz)
cc     .           ,Fluxz(0:nx,0:ny,0:nz)
cc
cc      logical,optional :: symmetric
cc
ccc     Local variables
cc
cc      real(8)  :: dtx,dty,dtz,a1,a2,sgn,bb,flx0,flxp
cc      integer  :: i,j,k,ig,jg,kg
cc      logical  :: limit,symm
cc
ccc     Begin program
cc
cc      Fluxx=0
cc      Fluxy=0
cc      Fluxz=0
cc
cc      if (PRESENT(symmetric)) then
cc        symm = symmetric
cc      else
cc        symm = .false.
cc      endif
cc
ccc     Symmetric discretization
cc
cc      if (symm) then
cc
cc        do k = 0,nz
cc          do j = 0,ny
cc            do i = 0,nx
cc
cc              call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc              dtx = .5*(tmp(i+1,j  ,k)-tmp(i,j  ,k)
cc     .                 +tmp(i+1,j+1,k)-tmp(i,j+1,k))/dx(ig)
cc              dty = .5*(tmp(i  ,j+1,k)-tmp(i  ,j,k)
cc     .                 +tmp(i+1,j+1,k)-tmp(i+1,j,k))/dy(jg)
cc              dtz = 0d0
cc
cc              Fluxx(i,j,k)= dtx
cc              Fluxy(i,j,k)= dty
cc              Fluxz(i,j,k)= dtz
cc            enddo
cc          enddo
cc        enddo
cc
cc        !Average fluxes to faces
cc        do k = nz,1,-1
cc          do j = ny,1,-1
cc            Fluxx(:,j,k) = 0.25*(Fluxx(:,j,k  )+Fluxx(:,j-1,k  )
cc     .                          +Fluxx(:,j,k-1)+Fluxx(:,j-1,k-1))
cc          enddo
cc        enddo
cc
cc        do k = nz,1,-1
cc          do i = nx,1,-1
cc            Fluxy(i,:,k) = 0.25*(Fluxy(i,:,k  )+Fluxy(i-1,:,k  )
cc     .                          +Fluxy(i,:,k-1)+Fluxy(i-1,:,k-1))
cc          enddo
cc        enddo
cc
cc        do j = ny,1,-1
cc          do i = nx,1,-1
cc            Fluxz(i,j,:) = 0.25*(Fluxz(i,j  ,:)+Fluxz(i-1,j  ,:)
cc     .                          +Fluxz(i,j-1,:)+Fluxz(i-1,j-1,:))
cc          enddo
cc        enddo
cc
ccc     Asymmetric discretization
cc
cc      else
cc
cc        !Add co-derivative terms: standard centered
cc        do i = 0,nx
cc          call getMGmap(i,1,1,igx,igy,igz,ig,jg,kg)
cc
cc          Fluxx(i,1:ny,1:nz) =
cc     .         (tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
cc        enddo
cc
cc        do j = 0,ny
cc          call getMGmap(1,j,1,igx,igy,igz,ig,jg,kg)
cc
cc          Fluxy(1:nx,j,1:nz) =
cc     .         (tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
cc        enddo
cc
cc        do k = 0,nz
cc          call getMGmap(1,1,k,igx,igy,igz,ig,jg,kg)
cc
cc          Fluxz(1:nx,1:ny,k) =
cc     .         (tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
cc        enddo
cc
cc      endif
cc
ccc     End program
cc
cc      end subroutine iso_flux_2nd

c     par_flux-4th
c     ################################################################
      subroutine par_flux_4th(i,j,k,nx,ny,nz,igx,igy,igz,F,Bx,By,Bz
     .                       ,Fluxx,Fluxy,Fluxz,limited,stencil)
c     ----------------------------------------------------------------
c     Calculates div.(BB.grad) in 3D with 4th order accuracy
c     Can deal with 1D, 2D, and 3D.
c     Checks if there are enough points for the stencil in all directions
c     CANNOT handle different dx, dy, dz, (except for reduced dimensions)
c     Need information about MG to do that, hence for now dx=dy=dz=h
c
c     Coded by Natalia Krashenninikova, 8/2008
c     Modified by L. Chacon, 9/2008
c     ----------------------------------------------------------------
      implicit none

c     Call variables

      integer  :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8)  :: F (0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
     .           ,By(0:nx+1,0:ny+1,0:nz+1)
     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)
      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
     .           ,Fluxz(i-1:i,j-1:j,k-1:k)

      integer,optional :: stencil
      logical,optional :: limited

c     Local variables

      real(8)  :: dx,dy,dz
      integer  :: tgx,tgy,tgz,xi,xf,yi,yf,zi,zf,ig,jg,kg,ii,jj,kk
     .           ,stncl
      logical  :: limit

      real(8),dimension(:,:,:),pointer :: jac

c     Begin program

      if (PRESENT(limited)) then
        limit = limited
      else
        limit = .true.
      endif

      if (PRESENT(stencil)) then
        stncl = stencil
      else
        stncl = 5
      endif

      Fluxx=0
      Fluxy=0
      Fluxz=0

c     Get (constant) mesh spacings

      call getMGmap(1,1,1,igx,igy,igz,ig,jg,kg)

      dx=grid_params%dxh(ig)
      dy=grid_params%dyh(jg)
      dz=grid_params%dzh(kg)

      jac => gmetric%grid(igx)%jac

c     Calculate fluxes

      !X flux
      if (nx > 1) then
        do ii=i-1,i
          tgx=IOR(min(ii-2,0),max(ii+2-nx,0))
          if (ii == 0) tgx = -2
          xi=max(ii-2,0)-max(ii+3-(nx+1),0)
          xf=min(ii+3,nx+1)+max(2-ii,0)
          Fluxx(ii,j,k)=Get_Flux(dx,tgx,F(xi:xf,j,k)
     .              ,Bx(xi:xf,j,k)*Bx(xi:xf,j,k)/jac(xi:xf,j,k))

          if (ny > 1) then
            tgx=IOR(min(ii-2,0),max(ii+2-(nx+1),0))
            tgy=IOR(min(j-2,0),max(j+2-(ny+1),0))
            if (ii == 0) tgx = -2
            xi=max(ii-2,0)-max(ii+2-(nx+1),0)
            xf=min(ii+2,nx+1)+max(2-ii,0)
            yi=max(j-2,0)-max(j+2-(ny+1),0)
            yf=min(j+2,ny+1)+max(2-j,0)
            Fluxx(ii,j,k)=Fluxx(ii,j,k)
     .                   +Get_Flux_Cross(dy,tgx,tgy,F(xi:xf,yi:yf,k)
     .                   ,Bx(xi:xf,j,k)*By(xi:xf,j,k)/jac(xi:xf,j,k))
          endif

          if (nz > 1) then
            tgx=IOR(min(ii-2,0),max(ii+2-(nx+1),0))
            if (ii == 0) tgx = -2
            tgz=IOR(min(k-2,0),max(k+2-(nz+1),0))
            xi=max(ii-2,0)-max(ii+2-(nx+1),0)
            xf=min(ii+2,nx+1)+max(2-ii,0)
            zi=max(k-2,0)-max(k+2-(nz+1),0)
            zf=min(k+2,nz+1)+max(2-k,0)

            Fluxx(ii,j,k)=Fluxx(ii,j,k)
     .                  +Get_Flux_Cross(dz,tgx,tgz,F(xi:xf,j,zi:zf)
     .                  ,Bx(xi:xf,j,k)*Bz(xi:xf,j,k)/jac(xi:xf,j,k))
         endif

        enddo
      endif

      !Y flux
      if (ny > 1) then
        do jj=j-1,j
          tgy=IOR(min(jj-2,0),max(jj+2-ny,0))
          if (jj == 0) tgy = -2
          yi=max(jj-2,0)-max(jj+3-(ny+1),0)
          yf=min(jj+3,ny+1)+max(2-jj,0)
          Fluxy(i,jj,k)=Get_Flux(dy,tgy,F(i,yi:yf,k)
     .                    ,By(i,yi:yf,k)*By(i,yi:yf,k)/jac(i,yi:yf,k))

          if (nx > 1) then
            tgx=IOR(min(i-2,0),max(i+2-(nx+1),0))
            tgy=IOR(min(jj-2,0),max(jj+2-(ny+1),0))
            if (jj == 0) tgy = -2
            xi=max(i-2,0)-max(i+2-(nx+1),0)
            xf=min(i+2,nx+1)+max(2-i,0)
            yi=max(jj-2,0)-max(jj+2-(ny+1),0)
            yf=min(jj+2,ny+1)+max(2-jj,0)
            Fluxy(i,jj,k)=Fluxy(i,jj,k)
     .                  +Get_Flux_Cross(dx,tgy,tgx
     .                     ,Transpose(F(xi:xf,yi:yf,k))
     .                     ,By(i,yi:yf,k)*Bx(i,yi:yf,k)/jac(i,yi:yf,k))
          endif

          if (nz > 1) then
            tgy=IOR(min(jj-2,0),max(jj+2-(ny+1),0))
            if (jj == 0) tgy = -2
            tgz=IOR(min(k-2,0),max(k+2-(nz+1),0))
            yi=max(jj-2,0)-max(jj+2-(ny+1),0)
            yf=min(jj+2,ny+1)+max(2-jj,0)
            zi=max(k-2,0)-max(k+2-(nz+1),0)
            zf=min(k+2,nz+1)+max(2-k,0)
            Fluxy(i,jj,k)=Fluxy(i,jj,k)
     .                   +Get_Flux_Cross(dz,tgy,tgz,F(i,yi:yf,zi:zf)
     .                     ,By(i,yi:yf,k)*Bz(i,yi:yf,k)/jac(i,yi:yf,k))
          endif
        enddo
      endif

      !Z flux
      if (nz > 1) then
        do kk=k-1,k
          tgz=IOR(min(kk-2,0),max(kk+2-nz,0))
          if (kk == 0) tgz = -2
          zi=max(kk-2,0)-max(kk+2-(nz+1),0)
          zf=min(kk+2,nz+1)+max(2-kk,0)
          Fluxz(i,j,kk)=Get_Flux(dz,tgz,F(i,j,zi:zf)
     .                    ,Bz(i,j,zi:zf)*Bz(i,j,zi:zf)/jac(i,j,zi:zf))

          if (ny > 1) then
            tgy=IOR(min(j-2,0),max(j+2-(ny+1),0))
            tgz=IOR(min(kk-2,0),max(kk+2-(nz+1),0))
            if (kk == 0) tgz = -2
            yi=max(j-2,0)-max(j+2-(ny+1),0)
            yf=min(j+2,ny+1)+max(2-j,0)
            zi=max(kk-2,0)-max(kk+2-(nz+1),0)
            zf=min(kk+2,nz+1)+max(2-kk,0)
            Fluxz(i,j,kk)=Fluxz(i,j,kk)
     .                  +Get_Flux_Cross(dy,tgz,tgy
     .                     ,Transpose(F(i,yi:yf,zi:zf))
     .                     ,Bz(i,j,zi:zf)*By(i,j,zi:zf)/jac(i,j,zi:zf))
          endif

          if (nx > 1) then
            tgx=IOR(min(i-2,0),max(i+2-(nx+1),0))
            tgz=IOR(min(kk-2,0),max(kk+2-(nz+1),0))
            if (kk == 0) tgz = -2
            xi=max(i-2,0)-max(i+2-(nx+1),0)
            xf=min(i+2,nx+1)+max(2-i,0)
            zi=max(kk-2,0)-max(kk+2-(nz+1),0)
            zf=min(kk+2,nz+1)+max(2-kk,0)
            Fluxz(i,j,kk)=Fluxz(i,j,kk)
     .                  +Get_Flux_Cross(dx,tgz,tgx
     .                     ,Transpose(F(xi:xf,j,zi:zf))
     .                     ,Bz(i,j,zi:zf)*Bx(i,j,zi:zf)/jac(i,j,zi:zf))
          endif
        enddo
      endif

      contains

c     ################################################################
      function Get_Flux_Cross(h,tagx,tagy,F,B2) result(val)
c     ----------------------------------------------------------------
c     Calculates the flux for d^2/dx dy part of div.(BB.grad) operator 
c     with 4th order accuracy at a point.
c     Can handle boundary points
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      integer  :: tagx,tagy
      real(8)  :: F(5,5),B2(5),val,h

c     Local variables
      integer  :: a(-1:1,5),cm(4),cl(5),cr(5)
      real(8)  :: twelve=12d0,bxby,flx0,flx1,flx
      logical  :: chk0,chk1,chk2

      val=0

      a(-1,:)=(/-3,-10,18,-6,1/)
      a(0,:)=(/1,-8,0,8,-1/)
      a(1,:)=(/-1,6,-18,10,3/)

      cm=(/-1,7,7,-1/)
      cl=(/2,17,-11,5,-1/)
      cr=(/-1,5,-11,17,2/)

      if((tagx+2)*(tagx+1)*tagx*(tagx-1)+(tagy+1)*tagy*(tagy-1)==0)then
         select case(tagx)
         case (-2)
            val=sum(cl*B2*matmul(F,a(tagy,:))/twelve/h)/twelve
         case (-1)
            val=sum(cm*B2(1:4)*matmul(F(1:4,:),a(tagy,:))
     .           /twelve/h)/twelve
         case (0)
            val=sum(cm*B2(2:5)*matmul(F(2:5,:),a(tagy,:))
     .           /twelve/h)/twelve
         case (1)
            val=sum(cr*B2*matmul(F,a(tagy,:))/twelve/h)/twelve
         end select
      else
         write(*,*)"Error in Get_Flux_Cross, wrong tags"
         write(*,*)"tagx=",tagx,"tagy=",tagy
         write(*,*)"Choices are -1, 0, or 1 for left boundary, "
         write(*,*)"interior or right boundary points"
         stop
      endif

      !Check for monotonicity (Sharma, Hammett, JCP 227 (2007))
      if (limit) then
        bxby = 0.5*(B2(4)+B2(3))

        flx0 = mc(F(3,3)-F(3,2),F(3,4)-F(3,3),chk0)
        flx1 = mc(F(4,3)-F(4,2),F(4,4)-F(4,3),chk1)
        flx  = mc(flx0         ,flx1         ,chk2)

cc        write (*,*) 'here',bxby,chk0,chk1,chk2
      
        if (chk0 .or. chk1 .or. chk2) then
          val = bxby/h*flx
cc          write (*,*) 'here'
cc          val = 0d0
        endif
      endif

      end function Get_Flux_Cross

c     ################################################################
      function Get_Flux(h,tag,F,B2) result (val)
c     ----------------------------------------------------------------
c     Calculates the flux for d^2/dx^2 part of div.(BB.grad) operator 
c     with 4th order accuracy at a point.
c     Can handle boundary points
c     ----------------------------------------------------------------

      implicit none

c     Call variables

      integer  :: tag
      real(8)  :: F(6),B2(6),val,h

c     Local variables
cc      integer  :: a(5,6)
      real(8)  :: sv2=72d0,shty=720d0
      integer  :: a(0:5,6),cm(4),cl(5),cr(5)
      real(8)  :: twelve=12d0,sixty=60d0

      select case(stncl)
      case(5)

        val=0
        a=0
        a(1,2:5)=(/0,1,4,-5/)
        a(2,2:5)=(/11,-54,33,10/)
        a(3,2:5)=(/-10,-33,54,-11/)
        a(4,2:5)=(/5,-4,-1,0/)

        select case(tag)
        case (-2)
           a(1,:)=(/-36, -185, 400, -230, 60, -9/)
           a(2,:)=(/-610, -290, 1770, -1340, 560, -90/)
           a(3,:)=(/116, 840, -1620, 970, -360, 54/)
           a(4,:)=(/-22, -430, 710, -360, 120, -18/)
           a(5,:)=(/12, 65, -120, 60, -20, 3/)
           val=sum(B2(1:5)*matmul(a(1:5,:),F)/h)/shty
           return
        case (-1)
           val=sum(B2(1:4)*matmul(a(1:4,2:5),F(1:4))/h)/sv2
           return
        case (0)
           val=sum(B2(2:5)*matmul(a(1:4,2:5),F(2:5))/h)/sv2
           return
        case (1)
           val=sum(B2(3:6)*matmul(a(1:4,2:5),F(3:6))/h)/sv2
           return
        case (2)
           a(1,:)=(/-3, 20, -60, 120, -65, -12/)
           a(2,:)=(/18, -120, 360, -710, 430, 22/)
           a(3,:)=(/-54, 360, -970, 1620, -840, -116/)
           a(4,:)=(/90, -560, 1340, -1770, 290, 610/)
           a(5,:)=(/9, -60, 230, -400, 185, 36/)
           val=sum(B2(2:6)*matmul(a(1:5,:),F)/h)/shty
           return
        case default
           write(*,*)"Error in Get_Flux, wrong tag=",tag
           write(*,*)"Choices are -2, -1, 0, 1 or 2 for left boundary, "
           write(*,*)"interior or right boundary points"
           return
        end select
      case(7)
        a(0,:)=(/-137,300,-300,200,-75,12/)
        a(1,:)=(/-12,-65,120,-60,20,-3/)
        a(2,:)=(/3,-30,-20,60,-15,2/)
        a(3,:)=(/-2,15,-60,20,30,-3/)
        a(4,:)=(/3,-20,60,-120,65,12/)
        a(5,:)=(/-12,75,-200,300,-300,137/)

        cm=(/-1,7,7,-1/)
        cl=(/2,17,-11,5,-1/)
        cr=(/-1,5,-11,17,2/)

        select case(tag)
        case (-2)
           val=sum(cl*B2(1:5)*matmul(a(0:4,:),F)/sixty/h)/twelve
           return
        case (-1)
           val=sum(cm*B2(1:4)*matmul(a(0:3,:),F)/sixty/h)/twelve
           return
        case (0)
           val=sum(cm*B2(2:5)*matmul(a(1:4,:),F)/sixty/h)/twelve
           return
        case (1)
           val=sum(cm*B2(3:6)*matmul(a(2:5,:),F)/sixty/h)/twelve
           return
        case (2)
           val=sum(cr*B2(2:6)*matmul(a(1:5,:),F)/sixty/h)/twelve
           return
        case default
           write(*,*)"Error in Get_Flux, wrong tag=",tag
           write(*,*)"Choices are -2, -1, 0, 1 or 2 for left boundary, "
           write(*,*)"interior or right boundary points"
           return
        end select
      case default
        call pstop('Get_Flux in par_flux_4th','Stencil not available')
      end select

      end function Get_Flux

c     minmod
c     #################################################################
      function minmod(a,b)

      real(8) :: a,b,minmod

      if (a*b <= 0d0) then
        minmod = 0d0
        return
      endif

      if (a > 0d0) then
        minmod = min(a,b)
      else
        minmod = max(a,b)
      endif

      end function minmod

c     mc
c     #################################################################
      function mc(a,b,chk_limited)

      real(8) :: a,b,mc,ho_flx

      logical :: chk_limited

      ho_flx = 0.5*(a+b)

      mc = minmod(2*minmod(a,b),ho_flx)

      chk_limited = (mc /= ho_flx)

      end function mc

      end subroutine par_flux_4th

      end module operators

c module equilibrium
c ######################################################################
      module equilibrium

        use problem_def

        real(8) :: dlambda,rshear,vparflow,vperflow,beta,E0(3),B0(3)

        real(8) :: eq_params(6)

        !Flux-surface averages
        real(8) :: bzz_avg,b2_avg,ssmax,x0,y0,z0,thmax

        logical :: closed_orb

        integer :: nh1,nh2,nh3,npw1=1,npw2=1,npw3=1

        logical :: odd(3),random

        real(8),dimension(10) :: pert

      end module equilibrium

c module nlfunction_setup
c ######################################################################
      module nlfunction_setup

        use grid

        use parameters

        use transport_params

        use equilibrium

        use auxiliaryVariables

        use operators

        use timeStepping

#if !defined(vec_pot)
        real(8),pointer,dimension(:,:,:):: rho,px,py,pz,bx,by,bz,tmp
#else
        real(8),pointer,dimension(:,:,:):: rho,px,py,pz,ax,ay,az,tmp
#endif

        integer :: advect               !Identifies advection scheme

        logical :: bootstrap =.false.   !Whether to calculate bootstrap current
cc     .            ,ion_hall  =.false.   !Whether to use ion EOM for Hall E
     .            ,fake_ve   =.true.    !Whether to use fake ve in div(Pe)
     .            ,solenoidal=.true.    !Whether we keep B solenoidal
     .            ,nc_eom_jxb=.false.   !Whether we use JxB instead of conserv. form
     .            ,nc_eom_gp =.false.   !Whether we use grad(p) instead of conserv. form
     .            ,nc_eom_v  =.false.   !Whether we use v.grad(v) instead of conserv. form
     .            ,solve_rho =.true.    !Whether we solve for rho or use rho=1
     .            ,alt_eom   =.false.   !Whether to use alternate EOM formulation in SP geometries
     .            ,sym_st    =.false.   !Whether to use symmetric stress tensor

        logical :: save_tn

        INTERFACE res
          module procedure res_mesh,res_ijk
        END INTERFACE

        INTERFACE vis
          module procedure vis_mesh,vis_ijk
        END INTERFACE

      contains

c     viscous_heat_src
c     ###################################################################
      function viscous_heat_src(i,j,k,nx,ny,nz,igx,igy,igz,visc
     .                         ,vcnv1,vcnv2) result(viscous)

c     -------------------------------------------------------------------
c     Finds viscous heat source  -Pi(v1):nabla(v2)
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,j,k,nx,ny,nz,igx,igy,igz
        real(8) :: vcnv1(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,vcnv2(0:nx+1,0:ny+1,0:nz+1,3)
     .            ,visc (0:nx+1,0:ny+1,0:nz+1)
     .            ,viscous

c     Local variables

        real(8) :: gsub(3,3),nabla_v(3,3),pi_tnsr(3,3)

c     Begin program

c     Viscous src: obtain cnv tensor pi, and dot with mixed tensor nabla_v (T_l^m)

        if (visc(i,j,k) > 0d0) then
          pi_tnsr=EOM_Pi_ijk(i,j,k,nx,ny,nz,igx,igy,igz,vcnv1,visc)

          nabla_v = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                      ,vcnv2(:,:,:,1)
     .                      ,vcnv2(:,:,:,2)
     .                      ,vcnv2(:,:,:,3),0)

          !Make nabla_v covariant on both sides, and do dot product
          gsub = gmetric%grid(igx)%gsub(i,j,k,:,:)

          viscous =-tensorScalarProduct(i,j,k,igx,matmul(nabla_v,gsub)
     .                                 ,pi_tnsr)
        else
          viscous = 0d0
        endif

c     End program

      end function viscous_heat_src

c     EOM_f
c     ###################################################################
      function EOM_f(i,j,k,nx,ny,nz,igx,igy,igz,rr,tt,bb,jj,vec_cnv
     $              ,nc_eom_jxb,nc_eom_gp,alt_eom,vol)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of the EOM force term
c     (jxb)-grad(ap*rho*tmp) at cell centers. Vectors are taken as
c     contravariant if vec_cnv=.true., and covariant otherwise
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nx,ny,nz,igx,igy,igz
        real(8)    :: cnv(3)
        real(8),target :: jj(0:nx+1,0:ny+1,0:nz+1,3)
     $                   ,bb(0:nx+1,0:ny+1,0:nz+1,3)
     $                   ,rr(0:nx+1,0:ny+1,0:nz+1)
     $                   ,tt(0:nx+1,0:ny+1,0:nz+1)
        logical    :: vec_cnv,nc_eom_jxb,nc_eom_gp,alt_eom
        logical,optional,intent(IN) :: vol

c     Local variables

        integer    :: igrid,ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: cov(3),jac,cnv1(3),dxx,dyy,dzz
        logical    :: vol_wgt

c     Begin program

        cnv = eom_jxb(i,j,k,nx,ny,nz,igx,igy,igz,bb,jj,vec_cnv
     $               ,nc_eom_jxb,alt_eom,vol)
     .       -eom_gp (i,j,k,nx,ny,nz,igx,igy,igz,rr,tt,nc_eom_gp
     .               ,alt_eom,vol)

c     End program

      end function EOM_f

c     EOM_jxb
c     ###################################################################
      function EOM_jxb(i,j,k,nx,ny,nz,igx,igy,igz,bb,jj,vec_cnv
     $                ,nc_eom_jxb,alt_eom,vol)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of the EOM force term
c     (jxb) at cell centers. Vectors are taken as
c     contravariant if vec_cnv=.true., and covariant otherwise
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nx,ny,nz,igx,igy,igz
        real(8)    :: cnv(3)
        real(8),target :: jj(0:nx+1,0:ny+1,0:nz+1,3)
     $                   ,bb(0:nx+1,0:ny+1,0:nz+1,3)
        logical    :: vec_cnv,nc_eom_jxb,nc_eom_gp,alt_eom
        logical,optional,intent(IN) :: vol

c     Local variables

        integer    :: igrid
        logical    :: vol_wgt

c     Begin program

        igrid = igx

        vol_wgt = .true.
        if (PRESENT(vol)) vol_wgt = vol

        !Lorentz force
        if (nc_eom_jxb) then

          cnv = crossProduct(i,j,k,igrid
     .                      ,jj(i,j,k,:),bb(i,j,k,:),vec_cnv)

          if (vec_cnv)
     .         cnv = matmul(gmetric%grid(igrid)%gsup(i,j,k,:,:),cnv)

          if (vol_wgt) cnv = cnv*gmetric%grid(igrid)%dvol(i,j,k)

        else
          if (.not.vec_cnv) call pstop('eom_jxb'
     .                ,'Magnetic field should be contravariant')

          vec1 => bb
          vec2 => jj
          cnv =-div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,eom_jxb_x,eom_jxb_y,eom_jxb_z,vol=vol_wgt)
          nullify(vec1,vec2)
        endif

c     End program

      end function EOM_jxb

ccc     EOM_jxb_x
ccc     #############################################################
cc      subroutine EOM_jxb_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t11,t12,t13,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t11-t13 for jxb force term in EOM.
ccc     In the call sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t11,t12,t13: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t11,t12,t13
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,ip
cc        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3),gsuper(3,3)
cc        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,pmag
cc
ccc     Begin program
cc
cc        ip = i+1
cc        if (flag == 0) ip = i
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
cc     .               +gmetric%grid(igx)%jac (i ,j,k))
cc        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
cc     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))
cc
cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv
cc
cc        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
cc     .      .and. bcSP()
cc     .      .and. flag /= 0           ) then
cc          jacp = gmetric%grid(igx)%jac(ip,j,k)
cc          jac0 = gmetric%grid(igx)%jac(i ,j,k)
cc        else
cc          jacp = jac
cc          jac0 = jac
cc        endif
cc
cc        ijac0 = 1d0/jac0
cc        ijacp = 1d0/jacp
cc
cc        cnv0 = vec1(i,j,k,:)
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                              ,cov0(1),cov0(2),cov0(3)
cc     .                              ,cnv0(1),cnv0(2),cnv0(3),.false.)
cc        cnvp = vec2(ip,j,k,:)
cc        call transformFromCurvToCurv(ip,j,k,igx,igy,igz
cc     .                              ,covp(1),covp(2),covp(3)
cc     .                              ,cnvp(1),cnvp(2),cnvp(3),.false.)
cc
cc        pmag = jac*(cnvp(1)*cov0(1)*ijacp
cc     .             +cnv0(1)*covp(1)*ijac0
cc     .             +cnvp(2)*cov0(2)*ijac0
cc     .             +cnv0(2)*covp(2)*ijacp
cc     .             +cnvp(3)*cov0(3)*ijacp
cc     .             +cnv0(3)*covp(3)*ijac0)*0.25
cc
cc        t11 = -0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,1)
cc     .              +vec1(ip,j,k,1)*vec2(i ,j,k,1))*ijacp*ijac0*jac**2
cc     .        +gsuper(1,1)*pmag
cc
cc        t12 = -0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijacp
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijac0 )*jac
cc     .        +gsuper(1,2)*pmag
cc
cc        t13 = -0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,3)
cc     .              +vec1(i ,j,k,1)*vec2(ip,j,k,3))*ijacp*ijac0*jac**2
cc     .        +gsuper(1,3)*pmag
cc
cc        if (flag /= 0) then
cc          ijac = 1d0/jac
cc
cc          t11 = t11*ijac
cc          if (.not.alt_eom) t12 = t12*ijac
cc          t13 = t13*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine EOM_jxb_x
cc
ccc     EOM_jxb_y
ccc     #############################################################
cc      subroutine EOM_jxb_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t21,t22,t23,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t21-t23 for jxb force term in EOM.
ccc     In the call sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t21,t22,t23: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t21,t22,t23
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,jp
cc        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3),gsuper(3,3)
cc        real(8)    :: jac,ijac,pmag
cc
ccc     Begin program
cc
cc
cccc        if (flag == 0) then
cccc          jp = j
cccc          jac    = gmetric%grid(igx)%jac (i,j ,k)
cccc          gsuper = gmetric%grid(igx)%gsup(i,j ,k,:,:)
cccc
cccc          cnv0 = vec1(i,j,k,:)
cccc          cov0 = matmul(gmetric%grid(igx)%gsup(i,j,k,:,:)
cccc     .                 ,vec1(i,j,k,:))
cccc        else
cccc          jp = j+1
cccc          jac    = 0.25*(gmetric%grid(igx)%jac_v(i-1,j,k-1)
cccc     .                  +gmetric%grid(igx)%jac_v(i-1,j,k  )
cccc     .                  +gmetric%grid(igx)%jac_v(i  ,j,k-1)
cccc     .                  +gmetric%grid(igx)%jac_v(i  ,j,k  ))
cccc          gsuper = 0.25*(gmetric%grid(igx)%gsup_v(i-1,j,k-1,:,:)
cccc     .                  +gmetric%grid(igx)%gsup_v(i-1,j,k  ,:,:)
cccc     .                  +gmetric%grid(igx)%gsup_v(i  ,j,k-1,:,:)
cccc     .                  +gmetric%grid(igx)%gsup_v(i  ,j,k  ,:,:))
cccc        endif
cc        jp = j+1
cc        if (flag == 0) jp = j
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
cc     .               +gmetric%grid(igx)%jac (i,j ,k))
cc        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
cc     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))
cc
cc        cnv0 = vec1(i,j,k,:)
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                              ,cov0(1),cov0(2),cov0(3)
cc     .                              ,cnv0(1),cnv0(2),cnv0(3),.false.)
cc        cnvp = vec2(i,jp,k,:)
cc        call transformFromCurvToCurv(i,jp,k,igx,igy,igz
cc     .                              ,covp(1),covp(2),covp(3)
cc     .                              ,cnvp(1),cnvp(2),cnvp(3),.false.)
cc
cc        pmag = (cnvp(1)*cov0(1)
cc     .         +cnv0(1)*covp(1)
cc     .         +cnvp(2)*cov0(2)
cc     .         +cnv0(2)*covp(2)
cc     .         +cnvp(3)*cov0(3)
cc     .         +cnv0(3)*covp(3))*0.25
cc
cc        t21 = -0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,1)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,1) )
cc     .        +gsuper(2,1)*pmag
cc
cc        t22 = -0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,2)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,2) )
cc     .        +gsuper(2,2)*pmag
cc
cc        t23 = -0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,3)
cc     .              +vec1(i,j ,k,2)*vec2(i,jp,k,3) )
cc     .        +gsuper(2,3)*pmag
cc
cc        if (flag /= 0) then
cc          ijac = 1d0/jac
cc
cc          t21 = t21*ijac
cc          if (.not.alt_eom) t22 = t22*ijac
cc          t23 = t23*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine EOM_jxb_y
cc
ccc     EOM_jxb_z
ccc     #############################################################
cc      subroutine EOM_jxb_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
cc     .                    ,t31,t32,t33,flag)
ccc     -------------------------------------------------------------
ccc     Calculates tensor components t31-t33 for jxb force term in EOM.
ccc     In the call sequence:
ccc       * i,j,k: grid position
ccc       * nx,ny,nz: grid size
ccc       * igx,igy,igz: grid level (for MG evaluations)
ccc       * alt_eom: whether to use alternate EOM in singular coord.
ccc                  systems or not.
ccc       * t31,t32,t33: tensor components.
ccc       * flag: whether evaluation is a cell center i,j,k (flag=0)
ccc               or at cell face i+1/2,j,k (flag /= 0)
ccc     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
ccc     head.
ccc     -------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
cc        real(8)    :: t31,t32,t33
cc        logical    :: alt_eom
cc
ccc     Local variables
cc
cc        integer    :: ig,jg,kg,kp
cc        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3),gsuper(3,3)
cc        real(8)    :: jac,ijac,pmag
cc
ccc     Begin program
cc
cc        kp = k+1
cc        if (flag == 0) kp = k
cc
cc        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
cc     .               +gmetric%grid(igx)%jac (i,j,k ))
cc        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
cc     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))
cc
cc        cnv0 = vec1(i,j,k,:)
cc        call transformFromCurvToCurv(i,j,k,igx,igy,igz
cc     .                              ,cov0(1),cov0(2),cov0(3)
cc     .                              ,cnv0(1),cnv0(2),cnv0(3),.false.)
cc        cnvp = vec2(i,j,kp,:)
cc        call transformFromCurvToCurv(i,j,kp,igx,igy,igz
cc     .                              ,covp(1),covp(2),covp(3)
cc     .                              ,cnvp(1),cnvp(2),cnvp(3),.false.)
cc
cc        pmag = (cnvp(1)*cov0(1)
cc     .         +cnv0(1)*covp(1)
cc     .         +cnvp(2)*cov0(2)
cc     .         +cnv0(2)*covp(2)
cc     .         +cnvp(3)*cov0(3)
cc     .         +cnv0(3)*covp(3))*0.25
cc
cc        t31 = -0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,1)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,1) )
cc     .        +gsuper(3,1)*pmag
cc
cc        t32 = -0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,2)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,2) )
cc     .        +gsuper(3,2)*pmag
cc
cc        t33 = -0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,3)
cc     .              +vec1(i,j,k ,3)*vec2(i,j,kp,3) )
cc     .        +gsuper(3,3)*pmag
cc
cc        if (flag /= 0) then
cc          ijac = 1d0/jac
cc
cc          t31 = t31*ijac
cc          if (.not.alt_eom) t32 = t32*ijac
cc          t33 = t33*ijac
cc        endif
cc
ccc     End program
cc
cc      end subroutine EOM_jxb_z

c     EOM_jxb_x
c     #############################################################
      subroutine EOM_jxb_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for jxb force term in EOM.
c     In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t11,t12,t13: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip
        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3),gsuper(3,3)
     .               ,car0(3),carp(3)
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,pmag,coeff

c     Begin program

cc        if (flag == 0) then
cc          ip = i
cc          jac    = gmetric%grid(igx)%jac (i,j ,k)
cc          gsuper = gmetric%grid(igx)%gsup(i,j ,k,:,:)
cc        else
cc          ip = i+1
cc          jac    = 0.25*(gmetric%grid(igx)%jac_v (i,j-1,k-1)
cc     .                  +gmetric%grid(igx)%jac_v (i,j-1,k  )
cc     .                  +gmetric%grid(igx)%jac_v (i,j  ,k-1)
cc     .                  +gmetric%grid(igx)%jac_v (i,j  ,k  ))
cc          gsuper = 0.25*(gmetric%grid(igx)%gsup_v(i,j-1,k-1,:,:)
cc     .                  +gmetric%grid(igx)%gsup_v(i,j-1,k  ,:,:)
cc     .                  +gmetric%grid(igx)%gsup_v(i,j  ,k-1,:,:)
cc     .                  +gmetric%grid(igx)%gsup_v(i,j  ,k  ,:,:))
cc        endif

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
     .               +gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))

cc        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

cc        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
cc     .      .and. bcSP()
cc     .      .and. flag /= 0           ) then
        if (flag /= 0) then
          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac0 = gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

        coeff = ijacp*ijac0*jac**2

cc        cnv0 = vec1(i,j,k,:)
cc        cov0 = matmul(gmetric%grid(igx)%gsub(i ,j,k,:,:),cnv0)
cc
cc        cnvp = vec2(ip,j,k,:)
cc        covp = matmul(gmetric%grid(igx)%gsub(ip,j,k,:,:),cnvp)
cc
cc        pmag = jac*(cnvp(1)*cov0(1)*ijacp
cc     .             +cnv0(1)*covp(1)*ijac0
cc     .             +cnvp(2)*cov0(2)*ijac0
cc     .             +cnv0(2)*covp(2)*ijacp
cc     .             +cnvp(3)*cov0(3)*ijacp
cc     .             +cnv0(3)*covp(3)*ijac0)*0.25
cc        pmag = (cnvp(1)*cov0(1)
cc     .         +cnv0(1)*covp(1)
cc     .         +cnvp(2)*cov0(2)
cc     .         +cnv0(2)*covp(2)
cc     .         +cnvp(3)*cov0(3)
cc     .         +cnv0(3)*covp(3))*0.25*sqrt(coeff)
cc        pmag = 0d0
cc        pmag = (cnvp(1)*covp(1)
cc     .         +cnv0(1)*cov0(1)
cc     .         +cnvp(2)*covp(2)
cc     .         +cnv0(2)*cov0(2)
cc     .         +cnvp(3)*covp(3)
cc     .         +cnv0(3)*cov0(3))*0.25

        car0 = XformToCar(i ,j,k,igx,vec1(i ,j,k,:),.false.)
        carp = XformToCar(ip,j,k,igx,vec2(ip,j,k,:),.false.)
        pmag = (carp(1)*car0(1)
     .         +car0(1)*carp(1)
     .         +carp(2)*car0(2)
     .         +car0(2)*carp(2)
     .         +carp(3)*car0(3)
     .         +car0(3)*carp(3))*0.25*jac

        t11 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,1)
     .               +vec1(i ,j,k,1)*vec2(ip,j,k,1) )*coeff
     .         +gsuper(1,1)*pmag)

        t12 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,2)
     .               +vec1(i ,j,k,1)*vec2(ip,j,k,2) )*coeff
     .         +gsuper(1,2)*pmag)

        t13 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,3)
     .               +vec1(i ,j,k,1)*vec2(ip,j,k,3) )*coeff
     .         +gsuper(1,3)*pmag)

        if (flag /= 0) then
          ijac = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine EOM_jxb_x

c     EOM_jxb_y
c     #############################################################
      subroutine EOM_jxb_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for jxb force term in EOM.
c     In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t21,t22,t23: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp
        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3),gsuper(3,3)
     .               ,car0(3),carp(3)
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,pmag,coeff

c     Begin program

cc        if (flag == 0) then
cc          jp = j
cc          jac    = gmetric%grid(igx)%jac (i,j ,k)
cc          gsuper = gmetric%grid(igx)%gsup(i,j ,k,:,:)
cc        else
cc          jp = j+1
cc          jac    = 0.25*(gmetric%grid(igx)%jac_v (i-1,j,k-1)
cc     .                  +gmetric%grid(igx)%jac_v (i-1,j,k  )
cc     .                  +gmetric%grid(igx)%jac_v (i  ,j,k-1)
cc     .                  +gmetric%grid(igx)%jac_v (i  ,j,k  ))
cc          gsuper = 0.25*(gmetric%grid(igx)%gsup_v(i-1,j,k-1,:,:)
cc     .                  +gmetric%grid(igx)%gsup_v(i-1,j,k  ,:,:)
cc     .                  +gmetric%grid(igx)%gsup_v(i  ,j,k-1,:,:)
cc     .                  +gmetric%grid(igx)%gsup_v(i  ,j,k  ,:,:))
cc        endif
        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))

        if (flag /= 0) then
          jacp = gmetric%grid(igx)%jac(i,jp,k)
          jac0 = gmetric%grid(igx)%jac(i,j ,k)
        else
          jacp = jac
          jac0 = jac
        endif

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

        coeff = ijacp*ijac0*jac**2

cc        cnv0 = vec1(i,j,k,:)
cc        cov0 = matmul(gmetric%grid(igx)%gsub(i,j,k,:,:),cnv0)
cc
cc        cnvp = vec2(i,jp,k,:)
cc        covp = matmul(gmetric%grid(igx)%gsub(i,jp,k,:,:),cnvp)
cc
cc        pmag = (cnvp(1)*cov0(1)
cc     .         +cnv0(1)*covp(1)
cc     .         +cnvp(2)*cov0(2)
cc     .         +cnv0(2)*covp(2)
cc     .         +cnvp(3)*cov0(3)
cc     .         +cnv0(3)*covp(3))*0.25

        car0 = XformToCar(i,j ,k,igx,vec1(i,j ,k,:),.false.)
        carp = XformToCar(i,jp,k,igx,vec2(i,jp,k,:),.false.)
        pmag = (carp(1)*car0(1)
     .         +car0(1)*carp(1)
     .         +carp(2)*car0(2)
     .         +car0(2)*carp(2)
     .         +carp(3)*car0(3)
     .         +car0(3)*carp(3))*0.25*jac

        t21 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,1)      
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,1) )*coeff
     .         +gsuper(2,1)*pmag)

        t22 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,2)      
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,2) )*coeff
     .         +gsuper(2,2)*pmag)

        t23 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,3)
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,3) )*coeff
     .         +gsuper(2,3)*pmag)

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine EOM_jxb_y

c     EOM_jxb_z
c     #############################################################
      subroutine EOM_jxb_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for jxb force term in EOM.
c     In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t31,t32,t33: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp
        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3),gsuper(3,3)
     .               ,car0(3),carp(3)
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,pmag,coeff

c     Begin program

cc        if (flag == 0) then
cc          kp = k
cc          jac    = gmetric%grid(igx)%jac (i,j ,k)
cc          gsuper = gmetric%grid(igx)%gsup(i,j ,k,:,:)
cc        else
cc          kp = k+1
cc          jac    = 0.25*(gmetric%grid(igx)%jac_v (i-1,j-1,k)
cc     .                  +gmetric%grid(igx)%jac_v (i-1,j  ,k)
cc     .                  +gmetric%grid(igx)%jac_v (i  ,j-1,k)
cc     .                  +gmetric%grid(igx)%jac_v (i  ,j  ,k))
cc          gsuper = 0.25*(gmetric%grid(igx)%gsup_v(i-1,j-1,k,:,:)
cc     .                  +gmetric%grid(igx)%gsup_v(i-1,j  ,k,:,:)
cc     .                  +gmetric%grid(igx)%gsup_v(i  ,j-1,k,:,:)
cc     .                  +gmetric%grid(igx)%gsup_v(i  ,j  ,k,:,:))
cc        endif
        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))

        if (flag /= 0) then
          jacp = gmetric%grid(igx)%jac(i,j,kp)
          jac0 = gmetric%grid(igx)%jac(i,j,k )
        else
          jacp = jac
          jac0 = jac
        endif

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

        coeff = ijacp*ijac0*jac**2

cc        cnv0 = vec1(i,j,k,:)
cc        cov0 = matmul(gmetric%grid(igx)%gsub(i,j,k,:,:),cnv0)
cc
cc        cnvp = vec2(i,j,kp,:)
cc        covp = matmul(gmetric%grid(igx)%gsub(i,j,kp,:,:),cnvp)
cc
cc        pmag = (cnvp(1)*cov0(1)
cc     .         +cnv0(1)*covp(1)
cc     .         +cnvp(2)*cov0(2)
cc     .         +cnv0(2)*covp(2)
cc     .         +cnvp(3)*cov0(3)
cc     .         +cnv0(3)*covp(3))*0.25

        car0 = XformToCar(i,j,k ,igx,vec1(i,j,k ,:),.false.)
        carp = XformToCar(i,j,kp,igx,vec2(i,j,kp,:),.false.)
        pmag = (carp(1)*car0(1)
     .         +car0(1)*carp(1)
     .         +carp(2)*car0(2)
     .         +car0(2)*carp(2)
     .         +carp(3)*car0(3)
     .         +car0(3)*carp(3))*0.25*jac

        t31 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,1)
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,1) )*coeff
     .         +gsuper(3,1)*pmag)

        t32 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,2)      
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,2) )*coeff
     .         +gsuper(3,2)*pmag)

        t33 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,3)
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,3) )*coeff
     .         +gsuper(3,3)*pmag)

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine EOM_jxb_z

c     EOM_gp
c     ###################################################################
      function EOM_gp(i,j,k,nx,ny,nz,igx,igy,igz,rr,tt,nc_eom_gp
     .               ,alt_eom,vol)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of the EOM force term
c     grad(ap*rho*tmp) at cell centers.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,nx,ny,nz,igx,igy,igz
        real(8)    :: cnv(3)
        real(8),target :: rr(0:nx+1,0:ny+1,0:nz+1)
     $                   ,tt(0:nx+1,0:ny+1,0:nz+1)
        logical    :: nc_eom_gp,alt_eom
        logical,optional,intent(IN) :: vol

c     Local variables

        integer    :: igrid,ig,jg,kg,ip,im,jp,jm,kp,km
        real(8)    :: cov(3),jac,cnv1(3),dxx,dyy,dzz
        logical    :: vol_wgt

c     Begin program

        igrid = igx

        vol_wgt = .true.
        if (PRESENT(vol)) vol_wgt = vol

        if (nc_eom_gp) then
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1

          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

          dxx = 2*grid_params%dxh(ig)
          dyy = 2*grid_params%dyh(jg)
          dzz = 2*grid_params%dzh(kg)

          if (adiabatic) then !tt is *electron* pressure
            cov(1) = (tt(ip,j,k)-tt(im,j,k))/dxx

            cov(2) = (tt(i,jp,k)-tt(i,jm,k))/dyy
          
            cov(3) = (tt(i,j,kp)-tt(i,j,km))/dzz
          else
            cov(1) = ( rr(i,j,k)*(tt(ip,j,k)-tt(im,j,k))/dxx
     .                +tt(i,j,k)*(rr(ip,j,k)-rr(im,j,k))/dxx)

            cov(2) = ( rr(i,j,k)*(tt(i,jp,k)-tt(i,jm,k))/dyy
     .                +tt(i,j,k)*(rr(i,jp,k)-rr(i,jm,k))/dyy)
          
            cov(3) = ( rr(i,j,k)*(tt(i,j,kp)-tt(i,j,km))/dzz
     .                +tt(i,j,k)*(rr(i,j,kp)-rr(i,j,km))/dzz)
          endif

          cov = cov*a_p         !Multiply by alpha_p=1 + Ti/Te to get total pressure

          call transformFromCurvToCurv(i,j,k,igx,igy,igz
     .                                ,cov(1),cov(2),cov(3)
     .                                ,cnv(1),cnv(2),cnv(3),.true.)

          if (vol_wgt) cnv = cnv*gmetric%grid(igrid)%dvol(i,j,k)

        else
          sc1 => rr
          sc2 => tt
          cnv = div_tensor(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,eom_grad_p_x,eom_grad_p_y,eom_grad_p_z
     .                    ,vol=vol_wgt)
          nullify(sc1,sc2)
        endif

c     End program

      end function EOM_gp

c     EOM_grad_p_x
c     #############################################################
      subroutine EOM_grad_p_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for force term in EOM.
c     In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t11,t12,t13: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,ptot,gsuper(3,3)

c     Begin program

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
     .               +gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(ip,j,k,:,:)
     .               +gmetric%grid(igx)%gsup(i ,j,k,:,:))

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        if (adiabatic) then    !p=T
          ptot = 0.5*jac*(sc2(ip,j,k)+sc2(i,j,k))
        else                   !p=nT
          ptot = 0.5*jac*(sc1(ip,j,k)*sc2(i ,j,k)
     .                   +sc1(i ,j,k)*sc2(ip,j,k))
        endif

        ptot = a_p*ptot         !Multiply by alpha_p=1 + Ti/Te to get total pressure

        t11 = gsuper(1,1)*ptot

        t12 = gsuper(1,2)*ptot

        t13 = gsuper(1,3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine EOM_grad_p_x

c     EOM_grad_p_y
c     #############################################################
      subroutine EOM_grad_p_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for force term in EOM.
c     In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t21,t22,t23: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp
        real(8)    :: x,y,z
        real(8)    :: jac,ijac,ptot,gsuper(3,3)

c     Begin program

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,jp,k,:,:)
     .               +gmetric%grid(igx)%gsup(i,j ,k,:,:))

        if (adiabatic) then    !p=T
          ptot = 0.5*jac*(sc2(i,jp,k)+sc2(i,j,k))
        else                   !p=nT
          ptot = 0.5*jac*(sc1(i,jp,k)*sc2(i,j ,k)
     .                   +sc1(i,j ,k)*sc2(i,jp,k))
        endif

        ptot = a_p*ptot         !Multiply by alpha_p=1 + Ti/Te to get total pressure

        t21 = gsuper(2,1)*ptot

        t22 = gsuper(2,2)*ptot

        t23 = gsuper(2,3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine EOM_grad_p_y

c     EOM_grad_p_z
c     #############################################################
      subroutine EOM_grad_p_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for force term in EOM.
c     In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t31,t32,t33: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp
        real(8)    :: x,y,z
        real(8)    :: jac,ijac,ptot,gsuper(3,3)

c     Begin program

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(gmetric%grid(igx)%gsup(i,j,kp,:,:)
     .               +gmetric%grid(igx)%gsup(i,j,k ,:,:))

        if (adiabatic) then    !p=T
          ptot = 0.5*jac*(sc2(i,j,kp)+sc2(i,j,k))
        else                   !p=nT
          ptot = 0.5*jac*(sc1(i,j,kp)*sc2(i,j,k )
     .                   +sc1(i,j,k )*sc2(i,j,kp))
        endif

        ptot = a_p*ptot         !Multiply by alpha_p=1 + Ti/Te to get total pressure

        t31 = gsuper(3,1)*ptot

        t32 = gsuper(3,2)*ptot

        t33 = gsuper(3,3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine EOM_grad_p_z

c     EOM_diff_x
c     #############################################################
      subroutine EOM_diff_x(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                     ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for EOM. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t11,t12,t13: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine assumes vec1->v (velocity), vec2->p (momentum).
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip,im,ipp,igpp,igp,igm
        real(8)    :: x,y,z,vxp,qip,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,vish,visp,vis0

c     Begin program

        if (     .not.associated(coef)
     .      .or. .not.associated(vec1)) then
          call pstop('EOM_diff_x','Pointers not associated')
        endif

        ip = i+1
        if (flag == 0) ip = i

        !Harmonic average for calculation of viscosity coeff. at faces
        if (coef(ip,j,k)*coef(i,j,k) > 0d0) then
          vish = 2.*coef(ip,j,k)*coef(i,j,k)/(coef(ip,j,k)+coef(i,j,k))
        else
          t11 = 0d0 ; t12 = 0d0; t13 = 0d0
          return
        endif

        jac    = 0.5*(gmetric%grid(igx)%jac (ip,j,k)
     .               +gmetric%grid(igx)%jac (i ,j,k))

        if ( i + grid_params%ilo(igx)-1 < grid_params%nxgl(igx)
     .      .and. bcSP()
     .      .and. flag /= 0           ) then
          jacp = gmetric%grid(igx)%jac(ip,j,k)
          jac0 = gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

        if (isSP(i+1,j,k,igx,igy,igz)) jac = SP_flsv

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

c     Diffusive part

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),1,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),0,cnv=.true.)
        endif

        if (sym_st) then
          nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
        endif

        t11 = -vish*nabla_v_cnv(1,1)
        t12 = -vish*nabla_v_cnv(1,2)
        t13 = -vish*nabla_v_cnv(1,3)

c       Postprocessing

        if (flag /= 0) then
          ijac  = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine EOM_diff_x

c     EOM_diff_y
c     #############################################################
      subroutine EOM_diff_y(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                     ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for EOM. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t21,t22,t23: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp,jm,jpp,jgpp,jgp,jgm
        real(8)    :: x,y,z,vyp,qjp,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,ijac,vish,visp,vis0

c     Begin program

        if (     .not.associated(coef)
     .      .or. .not.associated(vec1)) then
          call pstop('EOM_diff_y','Pointers not associated')
        endif

        jp = j+1
        if (flag == 0) jp = j

        !Harmonic average for calculation of viscosity coeff. at faces
        if (coef(i,jp,k)*coef(i,j,k) > 0d0) then
          vish = 2.*coef(i,jp,k)*coef(i,j,k)/(coef(i,jp,k)+coef(i,j,k))
        else
          t21 = 0d0 ; t22 = 0d0; t23 = 0d0
          return
        endif

        jac    = 0.5*(gmetric%grid(igx)%jac (i,jp,k)
     .               +gmetric%grid(igx)%jac (i,j ,k))

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c       Diffusive part

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),2,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),0,cnv=.true.)
        endif

        if (sym_st) then
          nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
        endif

        t21 = -vish*nabla_v_cnv(2,1)
        t22 = -vish*nabla_v_cnv(2,2)
        t23 = -vish*nabla_v_cnv(2,3)

c       Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine EOM_diff_y

c     EOM_diff_z
c     #############################################################
      subroutine EOM_diff_z(i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                     ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for EOM. In the call
c     sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t31,t32,t33: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     This routine has (rho,vx,vy,vz,bx,by,bz,tmp) passed via module
c     head.
c     -------------------------------------------------------------

        implicit none

c     Call variables

        integer    :: i,j,k,flag,igx,igy,igz,nx,ny,nz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp,km,kpp,kgpp,kgp,kgm
        real(8)    :: x,y,z,vzp,qkp,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,ijac,vish,visp,vis0

c     Begin program

        if (     .not.associated(coef)
     .      .or. .not.associated(vec1)) then
          call pstop('EOM_diff_z','Pointers not associated')
        endif

        kp = k+1
        if (flag == 0) kp = k

        !Harmonic average for calculation of viscosity coeff. at faces
        if (coef(i,j,kp)*coef(i,j,k) > 0d0) then
          vish = 2.*coef(i,j,kp)*coef(i,j,k)/(coef(i,j,kp)+coef(i,j,k))
        else
          t31 = 0d0 ; t32 = 0d0; t33 = 0d0
          return
        endif

        jac    = 0.5*(gmetric%grid(igx)%jac (i,j,kp)
     .               +gmetric%grid(igx)%jac (i,j,k ))

        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)

c       Diffusive part

        if (flag /= 0) then
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),3,cnv=.true.)
        else
          nabla_v_cnv = fnabla_v(i,j,k,nx,ny,nz,igx,igy,igz
     .                    ,vec1(:,:,:,1)
     .                    ,vec1(:,:,:,2)
     .                    ,vec1(:,:,:,3),0,cnv=.true.)
        endif

        if (sym_st) then
          nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
        endif

        t31 = -vish*nabla_v_cnv(3,1)
        t32 = -vish*nabla_v_cnv(3,2)
        t33 = -vish*nabla_v_cnv(3,3)

c       Postprocessing

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine EOM_diff_z

c     EOM_divPi_ijk
c     #####################################################################
      function EOM_divPi_ijk(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,diff)
     .         result(cnv)

c     ---------------------------------------------------------------------
c     Find divergence of ion pressure stress tensor. Return contravariant
c     components.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,nxx,nyy,nzz,igx,igy,igz
      real(8) :: cnv(3)
      real(8),target :: v0  (0:nxx+1,0:nyy+1,0:nzz+1,3)
      real(8),target :: diff(0:nxx+1,0:nyy+1,0:nzz+1)

c     Local variables

c     Begin program

      vec1 => v0
      coef => diff

      cnv = div_tensor(i,j,k,nxx,nyy,nzz,igx,igy,igz,alt_eom
     .                ,eom_diff_x,eom_diff_y,eom_diff_z
     .                ,vol=.false.)

      nullify(vec1,coef)

c     End program

      end function EOM_divPi_ijk

c     EOM_divPi
c     #####################################################################
      function EOM_divPi(nxx,nyy,nzz,igx,igy,igz,v0,diff) result(cnv)

c     ---------------------------------------------------------------------
c     Find divergence of electron pressure stress tensor. Returns
c     contravariant components.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nxx,nyy,nzz,igx,igy,igz
      real(8) :: cnv (0:nxx+1,0:nyy+1,0:nzz+1,3)
      real(8) :: v0  (0:nxx+1,0:nyy+1,0:nzz+1,3)
      real(8) :: diff(0:nxx+1,0:nyy+1,0:nzz+1)

c     Local variables

      integer :: i,j,k

c     Begin program

      cnv = 0d0

      do k=1,nzz
        do j=1,nyy
          do i=1,nxx
            cnv(i,j,k,:) = EOM_divPi_ijk(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                                  ,v0,diff)
          enddo
        enddo
      enddo

c     End program

      end function EOM_divPi

c     EOM_divPe_ijk
c     #####################################################################
      function EOM_divPe_ijk(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,diff)
     .         result(cnv)

c     ---------------------------------------------------------------------
c     Find divergence of electron pressure stress tensor. Returns
c     contravariant components.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,nxx,nyy,nzz,igx,igy,igz
      real(8) :: cnv(3)
      real(8) :: v0  (0:nxx+1,0:nyy+1,0:nzz+1,3)
      real(8) :: diff(0:nxx+1,0:nyy+1,0:nzz+1)

c     Local variables

c     Begin program

      cnv = EOM_divPi_ijk(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,diff)

cc      !Set div(Pi) x,y components to zero to match Slava's implementation
cc      cnv(1:2) = 0d0

c     End program

      end function EOM_divPe_ijk

c     EOM_divPe
c     #####################################################################
      function EOM_divPe(nxx,nyy,nzz,igx,igy,igz,v0,field) result(cnv)

c     ---------------------------------------------------------------------
c     Find divergence of electron pressure stress tensor. Returns
c     contravariant components.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nxx,nyy,nzz,igx,igy,igz
      real(8) :: cnv  (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .          ,v0   (0:nxx+1,0:nyy+1,0:nzz+1,3)
     .          ,field(0:nxx+1,0:nyy+1,0:nzz+1)

c     Local variables

      integer :: i,j,k

c     Begin program

      cnv = 0d0

      do k=1,nzz
        do j=1,nyy
          do i=1,nxx
            cnv(i,j,k,:) = EOM_divPe_ijk(i,j,k,nxx,nyy,nzz,igx,igy,igz
     .                                  ,v0,field)
          enddo
        enddo
      enddo

c     End program

      end function EOM_divPe

c     EOM_Pi_ijk
c     #####################################################################
      function EOM_Pi_ijk(i,j,k,nx,ny,nz,igx,igy,igz,v0,diff)
     .         result(pi_tnsr)

c     ---------------------------------------------------------------------
c     Find pressure stress tensor, in contravariant representation.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: pi_tnsr(3,3)
      real(8),target :: v0  (0:nx+1,0:ny+1,0:nz+1,3)
      real(8),target :: diff(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

c     Begin program

      vec1 => v0
      coef => diff
      call eom_diff_x(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .               ,pi_tnsr(1,1),pi_tnsr(1,2),pi_tnsr(1,3),0)
      call eom_diff_y(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .               ,pi_tnsr(2,1),pi_tnsr(2,2),pi_tnsr(2,3),0)
      call eom_diff_z(i,j,k,nx,ny,nz,igx,igy,igz,.false.
     .               ,pi_tnsr(3,1),pi_tnsr(3,2),pi_tnsr(3,3),0)
      nullify(vec1,coef)

c     End program

      end function EOM_Pi_ijk

c     E_fld_res
c     #####################################################################
      function E_fld_res(nx,ny,nz,igx,igy,igz,vcnv,bcnv,dumcov,eta2)
     .         result(E_res)

c     ---------------------------------------------------------------------
c     Computes ideal+resistive electric field, E=-(vxb)+eta*j.
c     On input, dumcov can be either jcov (for B-field version)
c     or acov (for vector potential version).
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,igx,igy,igz
      real(8) :: E_res (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,vcnv  (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,bcnv  (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,dumcov(0:nx+1,0:ny+1,0:nz+1,3)
     .          ,eta2  (0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer :: i,j,k,ig,jg,kg,ieq,bcnd(6,3)
      real(8) :: eps,betap,qa,jbs,cov(3),cnv(3),bnorm,x1,y1,z1
cc     .          ,vdummy(0:nx+1,0:ny+1,0:nz+1,3)

c     Begin program

      E_res(:,:,:,1) =-E0(1)
      E_res(:,:,:,2) =-E0(2)
      E_res(:,:,:,3) =-E0(3)

c     Resistive term

      if (eta > 0d0) then

cc      if (adiabatic) then
cc        pe = tmp
cc      else
cc        pe = rho*tmp
cc      endif
cc
cc        !Bootstrap current  (need to get out of here to treat BC separately)
cc        if (bootstrap) then
cc          !Get equilibrium configuration
cc          eps   = eq_params(1)
cc          betap = eq_params(2)
cc          qa    = eq_params(3)
cc
cc          do k=0,nz+1
cc            do j=0,ny+1
cc              do i=0,nx+1
cc                call getCurvilinearCoordinates(i,j,k,igx,igy,igz
cc     .                                        ,ig,jg,kg,x1,y1,z1)
cc
cc                !Find grad(p)
cc                cov = grad(i,j,k,nx,ny,nz,igx,igy,igz,pe)
cc
cc                bnorm = vectorNorm(i,j,k,igx,igy,igz
cc     .                       ,bcnv(i,j,k,1),bcnv(i,j,k,2),bcnv(i,j,k,3)
cc     .                       ,.false.)
cc
cc                jbs = -0.5*sqrt(eps*abs(x1)/bnorm)*betap/qa*(a_p*cov(1))
cc
cc                E_res(i,j,k,:) =-eeta(i,j,k)*jbs*bcov(i,j,k,:)    !Form J_parallel
cc              enddo
cc            enddo
cc          enddo
cc
cc#if !defined(vec_pot)
cc          !Impose topological BCs
cc          bcnd(:,1) = bcond
cc          bcnd(:,2) = bcond
cc          bcnd(:,3) = bcond
cc
cc          call setBC(IAX,3,nx,ny,nz,vdummy,E_res,vzeros,bcnd
cc     .              ,igx,igy,igz,iorder=2,is_cnv=.false.)
cc#endif
cc        endif

        !Ohmic current
#if !defined(vec_pot)
        E_res(:,:,:,1) = E_res(:,:,:,1)+eta2*dumcov(:,:,:,1)
        E_res(:,:,:,2) = E_res(:,:,:,2)+eta2*dumcov(:,:,:,2)
        E_res(:,:,:,3) = E_res(:,:,:,3)+eta2*dumcov(:,:,:,3)
#else
        do k=1,nz
          do j=1,ny
            do i=1,nx
              E_res(i,j,k,:) = E_res(i,j,k,:)
     .           - eta2(i,j,k)*veclap_cov(i,j,k,nx,ny,nz,igx,igy,igz
     .                                   ,dumcov,vol=.false.)
            enddo
          enddo
        enddo
#endif

      endif

c     VxB term

      E_res = E_res - crossProduct(igx,vcnv,bcnv,.true.)

c     End program

      end function E_fld_res

c     E_fld_Hall
c     #####################################################################
      function E_fld_Hall(nx,ny,nz,igx,igy,igz,rho,tmp) result(E_h)

c     ---------------------------------------------------------------------
c     Computes Hall electric field, in two forms: ion EOM, and e EOM
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,igx,igy,igz
      real(8) :: rho (0:nx+1,0:ny+1,0:nz+1)
     .          ,tmp (0:nx+1,0:ny+1,0:nz+1)
     .          ,E_h(0:nx+1,0:ny+1,0:nz+1,3)

c     Local variables

      integer :: i,j,k,ig,jg,kg,ieq,bcnd(6,3)
      real(8) :: pe(0:nx+1,0:ny+1,0:nz+1),cov(3),cnv(3)
cc     .          ,bcnv_star(0:nx+1,0:ny+1,0:nz+1,3)

c     Begin program

c     Find required auxiliary quantities

      if (adiabatic) then
        pe = tmp
      else
        pe = rho*tmp
      endif

cc      if (de > 0d0) then
cc        bcnv_star = bcnv -de**2/di*curl(nx,ny,nz,igx,igy,igz,vecov)
cc      else
cc        bcnv_star = bcnv
cc      endif

      !Standard Hall:  E_hall = (j x B* - grad(pe) - div(Pi_e))/rho - de^2/di^2*d(v_e)/dt
      if (.not.ion_hall) then

        E_h = crossProduct(igx,jcnv,bcnv,.true.)
     .       -grad(igx,pe)
     .       -XformVector(igx,div_pe,.false.)

        E_h(:,:,:,1) = E_h(:,:,:,1)/rho
        E_h(:,:,:,2) = E_h(:,:,:,2)/rho
        E_h(:,:,:,3) = E_h(:,:,:,3)/rho

        !Electron inertia term
        if (de > 0d0) then
          E_h = E_h - de**2/di**2/alpha*cnp/dt*(vecov-vecov_n)
        endif

      !Alternate Hall: E_hall = (dv/dt + v.grad(v) + div(Pi)/rho + grad(p_i)/rho)
      else

        !Temporal term: dv/dt
        if (save_tn.or.source_eval) then
          E_h = 0d0
        else
          if (.not.nc_eom_v) then
            if (cnm == 0d0) then
              E_h(:,:,:,1) = (vcnv(:,:,:,1)
     .                       -u_n %array_var(IVX )%array
     .                       /u_n %array_var(IRHO)%array)/alpha/dt
              E_h(:,:,:,2) = (vcnv(:,:,:,2)                          
     .                       -u_n %array_var(IVY )%array
     .                       /u_n %array_var(IRHO)%array)/alpha/dt
              E_h(:,:,:,3) = (vcnv(:,:,:,3)                          
     .                       -u_n %array_var(IVZ )%array
     .                       /u_n %array_var(IRHO)%array)/alpha/dt
            else
              E_h(:,:,:,1) = (cnp*vcnv(:,:,:,1)
     .                       +cn *u_n %array_var(IVX )%array
     .                           /u_n %array_var(IRHO)%array
     .                       +cnm*u_nm%array_var(IVX )%array
     .                           /u_nm%array_var(IRHO)%array)/alpha/dt
              E_h(:,:,:,2) = (cnp*vcnv(:,:,:,2)                          
     .                       +cn *u_n %array_var(IVY )%array
     .                           /u_n %array_var(IRHO)%array
     .                       +cnm*u_nm%array_var(IVY )%array
     .                           /u_nm%array_var(IRHO)%array)/alpha/dt
              E_h(:,:,:,3) = (cnp*vcnv(:,:,:,3)                          
     .                       +cn *u_n %array_var(IVZ )%array
     .                           /u_n %array_var(IRHO)%array
     .                       +cnm*u_nm%array_var(IVZ )%array
     .                           /u_nm%array_var(IRHO)%array)/alpha/dt
            endif
          else
            E_h(:,:,:,1) = (cnp*vcnv(:,:,:,1)
     .                     +cn *u_n %array_var(IVX)%array
     .                     +cnm*u_nm%array_var(IVX)%array)/alpha/dt
            E_h(:,:,:,2) = (cnp*vcnv(:,:,:,2)                      
     .                     +cn *u_n %array_var(IVY)%array          
     .                     +cnm*u_nm%array_var(IVY)%array)/alpha/dt
            E_h(:,:,:,3) = (cnp*vcnv(:,:,:,3)                      
     .                     +cn *u_n %array_var(IVZ)%array          
     .                     +cnm*u_nm%array_var(IVZ)%array)/alpha/dt
          endif
        endif

        !Stress
        E_h(:,:,:,1) = E_h(:,:,:,1) + div_pi(:,:,:,1)/rho
        E_h(:,:,:,2) = E_h(:,:,:,2) + div_pi(:,:,:,2)/rho
        E_h(:,:,:,3) = E_h(:,:,:,3) + div_pi(:,:,:,3)/rho

        !Advection
        do k=1,nz
          do j=1,ny
            do i=1,nx
              cnv = tensor_advc(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,vcnv
     .                         ,min(advect,2),cons=.false.)
     .             +E_h(i,j,k,:)

              E_h(i,j,k,:) = matmul(gmetric%grid(igx)%gsub(i,j,k,:,:)
     .                             ,cnv)
            enddo
          enddo
        enddo

        !Ion pressure term: grad(pi)/rho
        if (temp_ratio > 0d0) then
          do k=1,nz
            do j=1,ny
              do i=1,nx
                cov = grad(i,j,k,nx,ny,nz,igx,igy,igz,pe)

                E_h(i,j,k,:) = E_h(i,j,k,:) + temp_ratio*cov/rho(i,j,k)
              enddo
            enddo
          enddo
        endif

      endif

c     End program

      end function E_fld_Hall

c     res_ijk
c     #############################################################
      function res_ijk(i,j,k,nx,ny,nz,igx,igy,igz,field,rrr)
     .         result(res)
c     -------------------------------------------------------------
c     This function computes the resistivity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,nx,ny,nz,igx,igy,igz
      real(8) :: res,field(0:nx+1,0:ny+1,0:nz+1)

      real(8),optional :: rrr

c     Local variables

      integer    :: ig,jg,kg
      real(8)    :: rr

c     Begin program

      select case (equil)
      case ('ppnch','ppnsl','ppnst','ppn3d','p3nsl','2dhel')

        if (PRESENT(rrr)) then
          rr = rrr
        else
          call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
          rr = grid_params%xx(ig)
        endif

        !Resistivity profile eta*(1 + aa*x^bb)^cc (Daniele Bonfiglio)
        res = eta*(1. + aa_eta*rr**bb_eta)**cc_eta
cc        res = eta*(1. + aa_eta*grid_params%xx(ig)**nn_eta)

      case default

        res = eta

      end select

c     End program

      end function res_ijk

c     res_mesh
c     #############################################################
      function res_mesh(nx,ny,nz,igx,igy,igz,field,rrr) result(rsvty)
c     -------------------------------------------------------------
c     This function computes the resistivity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,igx,igy,igz
      real(8) :: field(0:nx+1,0:ny+1,0:nz+1)
     .          ,rsvty(0:nx+1,0:ny+1,0:nz+1)

      real(8),optional :: rrr

c     Local variables

      integer :: i,j,k

c     Begin program

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            rsvty(i,j,k)=res_ijk(i,j,k,nx,ny,nz,igx,igy,igz,psi
     .                          ,rrr=rrr)
          enddo
        enddo
      enddo

c     End program

      end function res_mesh

c     vis_ijk
c     #############################################################
      function vis_ijk(i,j,k,nx,ny,nz,igx,igy,igz) result(vis)
c     -------------------------------------------------------------
c     This function computes the viscosity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8)    :: vis
      integer    :: i,j,k,nx,ny,nz,igx,igy,igz

c     Local variables

      integer(4) :: ig,jg,kg

c     Begin program

      select case (equil)
      case ('ppnch','ppnsl','ppnst','ppn3d','p3nsl','2dhel')

        !Viscosity profile eta*(1 + aa*x^bb)^cc (Daniele Bonfiglio)
        call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
        vis = nu*(1. + aa_nu*grid_params%xx(ig)**bb_nu)**cc_nu

      case default

        vis = nu

      end select

c     End program

      end function vis_ijk

c     vis_mesh
c     #############################################################
      function vis_mesh(nx,ny,nz,igx,igy,igz) result(vsty)
c     -------------------------------------------------------------
c     This function computes the resistivity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      integer :: nx,ny,nz,igx,igy,igz
      real(8) :: vsty(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer :: i,j,k

c     Begin program

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            vsty(i,j,k) = vis_ijk(i,j,k,nx,ny,nz,igx,igy,igz)
          enddo
        enddo
      enddo

c     End program

      end function vis_mesh

c     hres
c     ###############################################################
      function hres(nx,ny,nz,igx,igy,igz,bcnv,rho) result(h_res)
c     ---------------------------------------------------------------
c     Finds hyper-resistivity coefficient (global on the grid)
c     ---------------------------------------------------------------

      use app_iosetup

      implicit none

c     Call variables

      integer :: nx,ny,nz,igx,igy,igz
      real(8) :: bcnv (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,rho  (0:nx+1,0:ny+1,0:nz+1)
     .          ,h_res(0:nx+1,0:ny+1,0:nz+1)

c     Local variables

      integer :: ig,jg,kg,i,j,k

      real(8) :: idx,idy,idz,k2,kk,bk_par,lheta,vol,lvol

c     Begin program

cc      hres = 0d0
cc      vol = 0d0
cc
cc      do k=1,nz
cc        do j=1,ny
cc          do i=1,nx
cc            call getMGmap(i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc            !kk
cc            idx  = pi/grid_params%dx(ig)
cc            if (nx == 1) idx = 0d0
cc            idy  = pi/grid_params%dy(jg)
cc            if (ny == 1) idy = 0d0
cc            idz  = pi/grid_params%dz(kg)
cc            if (nz == 1) idz = 0d0
cc
cc            k2=vectorNorm(i,j,k,igx,igy,igz,idx,idy,idz,.true.)
cc            kk=sqrt(k2)
cc
cc            !|B|*k_par
cc            bk_par=scalarProduct(i,j,k,igx,igy,igz,idx,idy,idz
cc     .                         ,bcnv(i,j,k,1)
cc     .                         ,bcnv(i,j,k,2)
cc     .                         ,bcnv(i,j,k,3))
cc
cccc            heta = heta + 0.1*di*sqrt(bnorm)/k2
cc            hres = hres +di/sqrt(rho(i,j,k))*abs(bk_par)
cc     .                     /kk**3
cccc     .                     /(kk**3 + dt*nu*kk**5)
cc     .                  *gmetric%grid(igx)%dvol(i,j,k)
cc            vol  = vol + gmetric%grid(igx)%dvol(i,j,k)
cc
cc          enddo
cc        enddo
cc      enddo
cc
cc#if defined(petsc) && !defined(samrai)
cc      lheta = hres
cc      call MPI_Allreduce(lheta,hres,1,MPI_DOUBLE_PRECISION
cc     .                        ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc      lvol = vol
cc      call MPI_Allreduce(lvol,vol,1,MPI_DOUBLE_PRECISION
cc     .                        ,MPI_SUM,MPI_COMM_WORLD,mpierr)
cc#else
cc#if defined(samrai)
cc      call pstop('applyAuxVarBC'
cc     .           ,'Need to fix hyperresistivity calculation in SAMRAI')
cc#endif
cc#endif
cc
cc      hres = 10*hres/vol/di**2  !Divide by di**2 to correct for multiplication later
cc
cc      if (my_rank == 0) write (*,*) 'Hyperresistivity=',hres

      if (di == 0d0 .or. test) then

        h_res = 0d0

      elseif (heta /= 0d0) then

cc        h_res = heta*rho
        h_res = heta

      else

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              call getMGmap(min(max(i,1),nx)
     .                     ,min(max(j,1),ny)
     .                     ,min(max(k,1),nz),igx,igy,igz,ig,jg,kg)

              !kk
              idx  = pi/grid_params%dx(ig)
              if (nx == 1) idx = 0d0
              idy  = pi/grid_params%dy(jg)
              if (ny == 1) idy = 0d0
              idz  = pi/grid_params%dz(kg)
              if (nz == 1) idz = 0d0

              k2=vectorNorm(i,j,k,igx,(/idx,idy,idz/),.true.)
              kk=sqrt(k2)

              !|B|*k_par
              bk_par=scalarProduct(i,j,k,igx,(/idx,idy,idz/)
     .                           ,bcnv(i,j,k,:))

              !h_res = C*rho*v_a*k||/k^3 (absent factor of di inherited later from v_e)
              h_res(i,j,k) = 1d1*sqrt(rho(i,j,k))*abs(bk_par)/kk**3
            enddo
          enddo
        enddo

      endif

c     End

      end function hres

      end module nlfunction_setup

c module precond_setup
c ######################################################################
      module precond_setup

        use grid

        use nlfunction_setup, ONLY: ion_hall

        integer :: pc_max_iter
     .            ,pc_debug_it
     .            ,pc_iter
     .            ,mg_ores
     .            ,mg_oprol
     .            ,mg_ores_si
     .            ,mg_oprol_si
     .            ,mg_coarse_size
     .            ,mg_vcyc
     .            ,sm_ncolors
     .            ,sm_iter
     .            ,sm_iter_si

        real(8) :: pc_tol
     .            ,sm_omega
     .            ,sm_omega_si

        character(2) :: pc_type
     .                 ,sm_type
     .                 ,sm_type_si

        logical :: mg_glrkin
     .            ,mg_gm_coarse
     .            ,mg_line_relax
     .            ,mg_zebra_relax
     .            ,mg_zebra_relax_si
     .            ,mg_debug
     .            ,mg_cvrg_test
     .            ,pc_asm
     .            ,pc_B_solve
     .            ,pc_v_solve
     .            ,pc_debug
     .            ,pc_ion_hall
     .            ,pc_divclean

cc        type :: pcsetup
cc          integer :: sm_iter,maxvcyc,pc_max_iter,debug_it
cc
cc          integer :: pc_iter,ores,oprol,mg_coarse_size,ncolors
cc
cc          real(8) :: mgtol,omeg_smth
cc
cc          character(2) :: smooth
cc          character(10):: pc_type
cc
cc          logical :: gm_smooth
cc     .              ,mg_glrkin
cc     .              ,gm_coarse
cc     .              ,line_relax
cc     .              ,zebra_relax
cc     .              ,inv_B
cc     .              ,v_solve
cc     .              ,mg_debug
cc     .              ,cvrg_test
cc     .              ,asm_PC
cc        end type pcsetup
cc
cc        type(pcsetup) :: pc__setup

      contains

c     PC_setup
c     ##################################################################
      subroutine PC_setup

      implicit none

c     ------------------------------------------------------------------
c     Hardwire setup variables for PC and perform consistency checks.
c     ------------------------------------------------------------------

c     Begin program

      pc_asm       =.false. !Whether to do additive Schwartz PC
      pc_B_solve   =.true.  !Whether to invert M_B or not
      pc_v_solve   =.false. !Whether to perform final v corrector step
      pc_ion_hall  = ion_hall  !Whether to use ion Ohm's law in PC

      mg_gm_coarse =.false. !Whether to use GMRES as the coarse solve
      mg_glrkin    =.false. !Whether to use Galerking MG
      mg_line_relax=.false. !Whether to perform line relaxation in MG
      mg_debug     =.false. !Whether to go into debugging mode in MG
      mg_cvrg_test =.false. !Whether to perform convergence test of smoother
      mg_coarse_size = 2    !Minimum coarse grid size (mg_ratio^mg_coarse_size)

      sm_type      = 'jb'   !Type of smoother ('jg','gs','gm')
      sm_omega     = 0.7    !Damping for MG smoother (only when 'jb' is chosen)
      sm_ncolors   = 8      !Number of colors for diagonal calc., smoothing

      !SCIDAC 08 conference data for Hall MHD
cc      pc_B_solve  = .true.
cc      mg_gm_coarse = .true.
cccc      ores   = 2
cccc      inv_B  = .false.
cccc      divCln = .false.

      !Specific setup
cc      inv_B  = .false.

      !Consistency checks
      if (pc_type == 's3') pc_B_solve = .true.    !Invert B for large-flow PC

      if (np > 1) mg_gm_coarse = .false.       !Do not perform GMRES coarse solve in parallel

      if (sm_type /= 'jb') sm_omega = 1d0

      !SI solver configuration
cc      sm_type_si     = 'gs'
cc      sm_omega_si    = 1d0
cc      sm_iter_si     = sm_iter/2
      sm_type_si     = sm_type
      sm_omega_si    = sm_omega
      sm_iter_si     = sm_iter
      mg_zebra_relax_si = mg_zebra_relax
      mg_ores_si     = mg_ores
      mg_oprol_si    = mg_oprol

      !Setup for singular point geometries
cc      if (bcSP()) then
cc        mg_zebra_relax_si = .true.
cc        mg_ores_si   = 1
cc        mg_oprol_si  = 1
cc
cc        sm_type_si   = 'gs'
cc        sm_omega_si  = 1d0
cc        sm_iter_si   = 2
cc      endif

      end subroutine PC_setup

      end module precond_setup
