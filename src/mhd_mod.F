c module app_iosetup
c ######################################################################
      module app_iosetup

        use var_io

        integer :: it_divcl=0

        logical :: dcon=.false.

        logical :: poincare=.false.,fourier=.false.,test=.false.

        integer :: u_equf=1000,u_prtf=2000
        character(100) :: equ_file,prt_file

        character(30) :: in_file  ='pixie3d.in' ! Default input file

      contains

c     write_bc_field_to_inputf
c     ###############################################################
      subroutine write_bc_field_to_inputf(var,desc)

c     ---------------------------------------------------------------
c     Writes BC fields to pixie3d.in for later reading.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: desc
      real(8) :: var(3)

c     Local variables

      integer :: ierr
      character(200):: command,line

c     Begin program

c     Transfer variable info to input file

      if (sum(var) /= 0d0 .and. my_rank == 0) then
        command = 'grep '//trim(desc)//' '//trim(in_file)
     .           // ' > /dev/null'
        ierr = shellcmd(trim(command))

        if (ierr == 0) ierr = delete_line(in_file,desc)

        write (line,*) '   '//trim(desc)//'=',var(1)
     .                                   ,',',var(2)
     .                                   ,',',var(3)

        command = 'grep datin ' // trim(in_file) // ' > /dev/null'
        ierr = shellcmd(trim(command))

        if (ierr == 0) then
          ierr = add_line_after_txt(in_file,'datin',line)
        else
          ierr = add_line_after_txt(in_file,'phys',line)
        endif

        if (ierr /= 0) then
          command = 'Could not write '//trim(desc)//' to input file'
          call pstop('setEq',command)
        endif
      endif

      end subroutine write_bc_field_to_inputf

c     write_real_field_to_inputf
c     ###############################################################
      subroutine write_real_field_to_inputf(var,desc)

c     ---------------------------------------------------------------
c     Writes real field to pixie3d.in for later reading.
c     ---------------------------------------------------------------

      implicit none

c     Call variables

      character(*) :: desc
      real(8) :: var

c     Local variables

      integer :: ierr
      character(200):: command,line

c     Begin program

c     Transfer variable info to input file

      if (my_rank == 0) then
        command = 'grep '//trim(desc)//' '//trim(in_file)
     .           // ' > /dev/null'
        ierr = shellcmd(trim(command))

        if (ierr == 0) ierr = delete_line(in_file,desc)

        write (line,*) '   '//trim(desc)//'=',var

        command = 'grep datin ' // trim(in_file) // ' > /dev/null'
        ierr = shellcmd(trim(command))

        if (ierr == 0) then
          ierr = add_line_after_txt(in_file,'datin',line)
        else
          ierr = add_line_after_txt(in_file,'phys',line)
        endif

        if (ierr /= 0) then
          command = 'Could not write '//trim(desc)//' to input file'
          call pstop('setEq',command)
        endif
      endif

      end subroutine write_real_field_to_inputf

      end module app_iosetup

c module problem_def
c ######################################################################
      module problem_def

        use variables

        implicit none

#if !defined(vec_pot)
        integer, parameter :: IRHO=1,IVX=2,IVY=3,IVZ=4,IBX=5,IBY=6,IBZ=7
     .                 ,IPRS=8,ITMPE=8,IJX=9,IJY=10,IJZ=11,IAX=12,IAY=13
     .                 ,IAZ=14
#else
        integer, parameter :: IRHO=1,IVX=2,IVY=3,IVZ=4,IAX=5,IAY=6,IAZ=7
     .                 ,IPRS=8,ITMPE=8,IBX=9,IBY=10,IBZ=11,IJX=12,IJY=13
     .                 ,IJZ=14
#endif

        character(5)  :: equil

        logical :: ion_hall     =.false.   !Whether to use ion EOM for Hall E
     .            ,electron_hall=.false.   !Whether to use e   EOM for Hall E
     .            ,pinch_flow   =.false.   !Whether equilibrium has pinch flow
     .            ,slava        =.false.   !Whether to use Slava Lukin's Hall MHD implementation
     .            ,noise        =.false.   !Whether to add noise to EOM (to break symmetry)
     .            ,adiabatic    =.false.   !Whether to use adiabatic EOS (no heat sources)
     .            ,adiabatic_rho=.false.   !Whether to use adiabatic formula for rho
     .            ,solve_rho    =.true.    !Whether we solve for rho or use rho=1
     .            ,solve_prs    =.true.    !Whether we solve for total prs or e tmp
     .            ,use_p_eom    =.false.   !Whether we use pressure or n*T to drive EOM
     .            ,use_p_diff   =.false.   !Whether we use pressure in diffusion term in pressure eq.
     .            ,use_p_classic_eq =.false.!Whether we use classic (p.div(v)) pressure evolution eq.
     .            ,nu_is_rhonu  =.false.   !Whether to use rho*nu in electron viscous term
     .            ,heta_is_rhoheta=.true.  !Whether to use rho*heta in ion viscoust term
     .            ,e_cons_vstar =.true.    !Whether we use energy conserving velocity in T eq.

        real(8) :: noise_lev

      end module problem_def

c module transport_params
c ######################################################################
      module transport_params

        use grid

        use problem_def

        integer :: xport_order=2
        
        real(8) :: nu,eta,dd,chi,chi_par,gamma,di,de,a_p,temp_ratio
     .            ,heta=0d0,c_hyper,beta=1d-3,te0_ref

        real(8) :: aa_eta,bb_eta,cc_eta
     .            ,aa_nu ,bb_nu ,cc_nu
     .            ,aa_chi,bb_chi ! Daniele, 02/07/2014
     .            ,cc_chi,dd_chi ! Daniele, 02/17/2014

        logical :: lagrangian=.false.     !Whether we perform Lagrangian step for Te
     .            ,spitzer   =.false.     !Whether we are using Spitzer resistivity
     .            ,fixed_Te  =.false.     !Whether we are fixing the T to T_0
     .            ,lagged_Te =.false.     !Whether we are lagging T to T_n
     .            ,chi_bgnskii =.false.   !Whether we are using Braginskii conductivities

      contains

c     res_rfx
c     #############################################################
      function res_rfx(rr) result(res)
c     -------------------------------------------------------------
c     This function computes the resistivity profile at a radius rr.
c     (Daniele's RFX profile).
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: res,rr

c     Local variables

      real(8) :: rf

c     Begin program

      rf = 1. + aa_eta*rr**bb_eta

      !Resistivity profile eta*(1 + aa*x^bb)^cc
      if (rf >= 0d0) then
         res = rf**cc_eta
      else
         res = 3e2
      endif

      res = min(res,3d2)

c     End program

      end function res_rfx

c     vis_rfx
c     #############################################################
      function vis_rfx(rr) result(vis)
c     -------------------------------------------------------------
c     This function computes the viscosity profile at a radius rr.
c     (Daniele's RFX profile).
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: vis,rr

c     Local variables

c     Begin program

      !Viscosity profile (1 + aa*x^bb)^cc
      vis = (1. + aa_nu*rr**bb_nu)**cc_nu

c     End program

      end function vis_rfx

 ! Daniele, 02/07/2014
c     chi_rfx
c     #############################################################
      function chi_rfx(rr) result(dfvty)
c     ------------------------------------------------------------- 
c     This function computes the heat diffusivity profile at a 
c     radius rr. (Daniele's RFX profile).
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: dfvty,rr

c     Begin program

      !Diffusivity profile 
cc      dfvty = (1.-(1.-aa_chi)*sin(pi*rr/bb_chi))
      !Daniele, 02/17/2014
      dfvty = (1.-(1.-aa_chi)
     .     *sin(pi*bb_chi*abs(rr)**cc_chi)
     .     *abs(rr)**dd_chi)

c     End program

      end function chi_rfx

c     res_spitzer
c     #############################################################
      function res_spitzer(te,te0) result(res)
c     -------------------------------------------------------------
c     This function computes the Spitzer resistivity for ELECTRON
c     temperature value te.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: eta0,res,te,te0

c     Local variables

cc      real(8) :: f0=0d0 ! Daniele, 02/07/2014

c     Begin program

cc      res = ((1+f0)/(abs(te)+f0*beta))**(1.5)
      res = (abs(te/te0))**(-1.5)

c     End program

      end function res_spitzer

c     diff_def
c     #############################################################
      function diff_def(te,rho,b2,te0) result(dfvty)
c     -------------------------------------------------------------
c     This function computes the particle diffusivity as a 
c     function of Te, rho, B^2.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: dfvty,te,rho,b2,te0

c     Local variables

c     Begin program

      dfvty = 1d0
cc      dfvty = sqrt(abs(te0/te))*rho*rho/b2

c     End program
      
      end function diff_def
      
c     chi_perp_brag
c     #############################################################
      function chi_perp_brag(te,te0,rho,b2) result(dfvty)
c     -------------------------------------------------------------
c     This function computes the perpendicular heat transport
c     diffusivity for ELECTRON temperature value te.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: dfvty,te,rho,b2,te0

c     Local variables

c     Begin program

      if (b2 > 0d0.and.te/=0d0) then
         dfvty = sqrt(abs(te0/te))*rho*rho/b2
      else
         dfvty = 0d0
      endif

c     End program
      
      end function chi_perp_brag

c     chi_par_brag
c     #############################################################
      function chi_par_brag(te,te0) result(dfvty)
c     -------------------------------------------------------------
c     This function computes the parallel heat transport diffusivity
c     for ELECTRON temperature value te.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      real(8) :: dfvty,te,te0

c     Local variables

c     Begin program

      dfvty = (abs(te/te0))**(2.5)

c     End program

      end function chi_par_brag


      end module transport_params

c module auxiliaryVariables
c ######################################################################
      module auxiliaryVariables

        use variables

        use transport_params

        integer :: advect,v_advect      !Identifies advection scheme for scalars, vectors
        logical :: nc_eom_v    =.false. !Whether we use v.grad(v) instead of conserv. form

#if !defined(vec_pot)
        integer, parameter :: IJCNV   =1
     .                       ,IJCNV_0 =2
     .                       ,IJCOV   =3
     .                       ,IJCOV_0 =4

        integer, parameter :: IBCNV   =5
     .                       ,IBCNV_0 =6
     .                       ,IBCOV   =7
     .                       ,IBCOV_0 =8
     .                       ,IBCOV_N =22

        integer, parameter :: IVCNV   =9
     .                       ,IVCNV_N =10
     .                       ,IVCNV_0 =11
     .                       ,IVCOV   =12

        integer, parameter :: IVECNV  =13
     .                       ,IVECOV  =14
     .                       ,IVECOV_N=15
     .                       ,IVECNV_0=21

        integer, parameter :: IENI    =16
     .                       ,IENI_0  =17
     .                       ,IEH     =18

        integer, parameter :: IDIVPI  =19
     .                       ,IDIVPE  =20
#else
        integer, parameter :: IJCNV   =1
     .                       ,IJCNV_0 =2
     .                       ,IJCOV   =3
     .                       ,IJCOV_0 =4

        integer, parameter :: IBCNV   =5
     .                       ,IBCNV_0 =6
     .                       ,IBCOV   =7
     .                       ,IBCOV_0 =8
     .                       ,IBCOV_N =25

        integer, parameter :: IACNV   =9
     .                       ,IACOV   =10
     .                       ,IACOV_0 =11
     .                       ,IACOV_N =26

        integer, parameter :: IVCNV   =12
     .                       ,IVCNV_0 =13
     .                       ,IVCNV_N =14
     .                       ,IVCOV   =15

        integer, parameter :: IVECNV  =16
     .                       ,IVECOV  =17
     .                       ,IVECOV_N=18
     .                       ,IVECNV_0=24

        integer, parameter :: IENI    =19
     .                       ,IEH     =20
     .                       ,IENI_0  =21

        integer, parameter :: IDIVPI  =22
     .                       ,IDIVPE  =23
#endif

        integer, parameter :: IETA     =1
     .                       ,IETA_0   =2
     .                       ,INU      =3
     .                       ,IHETA    =4
     .                       ,IRHO_0   =5
     .                       ,ITMP_0   =6
     .                       ,IRHO_AUX =7
     .                       ,ITMP_AUX =8
     .                       ,ITMPE_0  =9
     .                       ,ITMPE_AUX=10
     .                       ,ICHI     =11
     .                       ,ICHIPAR  =12
     .                       ,IDFF     =13

#if defined(vec_pot)
        real(8),pointer,dimension(:,:,:,:) :: acov
     .                                       ,acov_0
     .                                       ,acnv
     .                                       ,acov_n
#endif

        real(8),pointer,dimension(:,:,:)   :: nuu
     .                                       ,dff
     .                                       ,eeta
     .                                       ,eeta_0
     .                                       ,h_eta
     .                                       ,prs
     .                                       ,rho_0
     .                                       ,tmp_0
     .                                       ,tmpe_0
     .                                       ,rho
     .                                       ,tmp
     .                                       ,tmpe
     .                                       ,cchi
     .                                       ,cchi_par
     .                                       ,rho_advc
     .                                       ,prs_advc
     .                                       ,tmp_advc
     .                                       ,heat_xport
     .                                       ,heat_src

        real(8),pointer,dimension(:,:,:,:) :: bcnv
     .                                       ,bcnv_0
     .                                       ,bcov
     .                                       ,bcov_0
     .                                       ,bcov_n
     .                                       ,chi_tnsr

        real(8),pointer,dimension(:,:,:,:) :: vcnv
     .                                       ,vcnv_0
     .                                       ,vcnv_n
     .                                       ,vcov
     .                                       ,pcnv => null() !Used to check

        real(8),pointer,dimension(:,:,:,:) :: bhat
     .                                       ,v_advc
     .                                       ,E_ni
     .                                       ,E_ni_0
     .                                       ,E_h
     .                                       ,div_pi
     .                                       ,div_pe
     .                                       ,jcnv
     .                                       ,jcnv_0
     .                                       ,jcov
     .                                       ,jcov_0

        real(8),pointer,dimension(:,:,:,:) :: vecnv
     .                                       ,vecov
     .                                       ,vecov_n
     .                                       ,vecnv_0
     .                                       ,vscnv

        real(8),pointer,dimension(:,:,:,:) :: E_bc
        
      contains

c     createAuxVariables
c     ###################################################################
      subroutine createAuxVariables(vaux)

c     -------------------------------------------------------------------
c     Allocates auxiliary variable storage, sets up pointers.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(aux_array),pointer :: vaux

c     Local variables

        integer :: nx,ny,nz

c     Begin program

        if (.not.associated(vaux)) then
          call pstop('defineAuxVariables'
     .              ,'Auxiliary variable structured not allocated')
        endif

c     Setup pointers

        !Equilibrium variables
        vaux%var_list(IETA_0)%descr = 'eta_0'
        vaux%var_list(IETA_0)%time_0 = .true.
        eeta_0=>vaux%var_list(IETA_0)%array

        vaux%var_list(IRHO_0)%descr = 'rho_0'
        vaux%var_list(IRHO_0)%time_0 = .true.
        rho_0=> vaux%var_list(IRHO_0)%array

        vaux%var_list(ITMP_0)%descr = 'T_0'
        vaux%var_list(ITMP_0)%time_0 = .true.
        tmp_0=> vaux%var_list(ITMP_0)%array

        vaux%var_list(ITMPE_0)%descr = 'Te_0'
        vaux%var_list(ITMPE_0)%time_0 = .true.
        tmpe_0=> vaux%var_list(ITMPE_0)%array

#if defined(vec_pot)
        vaux%vec_list(IACOV_0)%cnv   = .false.
        vaux%vec_list(IACOV_0)%descr = 'A_0 cov'
        vaux%vec_list(IACOV_0)%time_0 = .true.
        acov_0 => vaux%vec_list(IACOV_0)%vec
#endif

        vaux%vec_list(IBCNV_0)%cnv   = .true.
        vaux%vec_list(IBCNV_0)%descr = 'B_0 cnv'
        vaux%vec_list(IBCNV_0)%time_0 = .true.
        bcnv_0 => vaux%vec_list(IBCNV_0)%vec

        vaux%vec_list(IBCOV_0)%cnv = .false.
        vaux%vec_list(IBCOV_0)%descr = 'B_0 cov'
        vaux%vec_list(IBCOV_0)%time_0 = .true.
        bcov_0 => vaux%vec_list(IBCOV_0)%vec

        vaux%vec_list(IJCNV_0)%cnv = .true.
        vaux%vec_list(IJCNV_0)%descr = 'J_0 cnv'
        vaux%vec_list(IJCNV_0)%time_0 = .true.
        jcnv_0 => vaux%vec_list(IJCNV_0)%vec

        vaux%vec_list(IJCOV_0)%cnv = .false.
        vaux%vec_list(IJCOV_0)%descr = 'J_0 cov'
        vaux%vec_list(IJCOV_0)%time_0 = .true.
        jcov_0 => vaux%vec_list(IJCOV_0)%vec

        vaux%vec_list(IVCNV_0)%cnv = .true.
        vaux%vec_list(IVCNV_0)%descr = 'V_0 cnv'
        vaux%vec_list(IVCNV_0)%time_0 = .true.
        vcnv_0 => vaux%vec_list(IVCNV_0)%vec

        vaux%vec_list(IENI_0)%cnv   = .false.
        vaux%vec_list(IENI_0)%descr = 'E_0 cov'
        vaux%vec_list(IENI_0)%time_0 = .true.
        E_ni_0 => vaux%vec_list(IENI_0)%vec

        vaux%vec_list(IVECNV_0)%cnv   = .true.
        vaux%vec_list(IVECNV_0)%descr = 'Ve_0 cnv'
        vaux%vec_list(IVECNV_0)%time_0 = .true.
        vecnv_0 => vaux%vec_list(IVECNV_0)%vec

        !Old-time variables
        vaux%vec_list(IVCNV_N)%cnv = .true.
        vaux%vec_list(IVCNV_N)%descr = 'V_n cnv'
        vaux%vec_list(IVCNV_N)%time_n= .true.
        vcnv_n => vaux%vec_list(IVCNV_N)%vec

        vaux%vec_list(IVECOV_N)%cnv   = .false.
        vaux%vec_list(IVECOV_N)%descr = 'Ve_n cov'
        vaux%vec_list(IVECOV_N)%time_n= .true.
        vecov_n => vaux%vec_list(IVECOV_N)%vec

        vaux%vec_list(IBCOV_N)%cnv = .false.
        vaux%vec_list(IBCOV_N)%descr = 'B_n cov'
        vaux%vec_list(IBCOV_N)%time_n= .true.
        bcov_n => vaux%vec_list(IBCOV_N)%vec

        !New-time variables
        vaux%var_list(IRHO_AUX)%descr = 'rho'
        rho => vaux%var_list(IRHO_AUX)%array

        vaux%var_list(ITMP_AUX)%descr = 'T=Te+Ti'
        tmp => vaux%var_list(ITMP_AUX)%array

        vaux%var_list(ITMPE_AUX)%descr = 'Te'
        tmpe => vaux%var_list(ITMPE_AUX)%array

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

        vaux%vec_list(IACOV_N)%cnv = .false.
        vaux%vec_list(IACOV_N)%descr = 'A_n cov'
        vaux%vec_list(IACOV_N)%time_n= .true.
        acov_n => vaux%vec_list(IACOV_N)%vec
#endif

        vaux%vec_list(IVECNV)%cnv   = .true.
        vaux%vec_list(IVECNV)%descr = 'Ve cnv'
        vecnv => vaux%vec_list(IVECNV)%vec

        vaux%vec_list(IVECOV)%cnv   = .false.
        vaux%vec_list(IVECOV)%descr = 'Ve cov'
        vecov => vaux%vec_list(IVECOV)%vec

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

        vaux%var_list(IDFF)%descr = 'diff'
        dff  => vaux%var_list(IDFF)%array

        vaux%var_list(ICHI)%descr = 'chi perp'
        cchi => vaux%var_list(ICHI)%array

        vaux%var_list(ICHIPAR)%descr = 'chi par'
        cchi_par => vaux%var_list(ICHIPAR)%array

cc        vaux%vec_list(IVSCNV)%cnv   = .true.
cc        vaux%vec_list(IVSCNV)%descr = 'Vstar cnv'
cc        vscnv => vaux%vec_list(IVSCNV)%vec

cc        vaux%vec_list(IPCNV)%cnv   = .true.
cc        vaux%vec_list(IPCNV)%descr = 'P=rho*v cnv'
cc        pcnv => vaux%vec_list(IPCNV)%vec

cc        vaux%vec_list(IBHAT)%cnv   = .true.
cc        vaux%vec_list(IBHAT)%descr = 'B normalized'
cc        bhat => vaux%vec_list(IBHAT)%vec
        
c     End program

      end subroutine createAuxVariables
      
c     deallocLocalAuxVars
c     ###################################################################
      subroutine deallocLocalAuxVars

c     -------------------------------------------------------------------
c     Deallocates LOCAL auxiliary variable storage.
c     -------------------------------------------------------------------

      integer :: istat

c     Begin program

      deallocate(pcnv,vscnv,prs,bhat,v_advc,rho_advc,prs_advc,tmp_advc
     .          ,heat_xport,heat_src,chi_tnsr,stat=istat)

      nullify(bcnv
     .       ,bcnv_0
     .       ,bcov
     .       ,bcov_0
     .       ,bcov_n
     .       ,jcnv
     .       ,jcnv_0
     .       ,jcov
     .       ,jcov_0
     .       ,vcnv
     .       ,vcnv_0
     .       ,vcnv_n
     .       ,vcov
#if defined(vec_pot)
     .       ,acnv
     .       ,acov
     .       ,acov_0
     .       ,acov_n
#endif
     .       ,vecnv
     .       ,vecov
     .       ,vecov_n
     .       ,vecnv_0
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
     .       ,tmpe_0
     .       ,rho
     .       ,tmp
     .       ,tmpe
     .       ,cchi
     .       ,cchi_par
cc     .       ,vecnv_0
cc     .       ,vscnv
cc     .       ,pcnv
cc     .       ,bhat
     .       )

c     End program

      end subroutine deallocLocalAuxVars

      end module auxiliaryVariables

c module auxPlotVariables
c ######################################################################
      module auxPlotVariables

        use variables

        use auxiliaryVariables, ONLY: heat_src

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

        type(var_array),pointer :: u_pert => null()

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
     .         ,p_tot (ilom:ihip,jlom:jhip,klom:khip)
     .         ,heat_src(ilom:ihip,jlom:jhip,klom:khip))

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
     .           ,jpar,jpsi,jperp,heat_src)

      deallocate (qfactor,lambda)

      if (associated(u_pert)) call deallocateDerivedType(u_pert)

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
      function curl_bxv(g_def,i,j,k,igr,vv,bb,half_elem) result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of curl(bb x vv) at the
c     grid node (i,j,k). One sided derivatives are employed when half_elem=1
c     (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,j,k,half_elem,igr
        real(8) :: cnv(3)
        real(8) :: vv(0:,0:,0:,:)
     $            ,bb(0:,0:,0:,:)
        type(grid_mg_def),pointer :: g_def

c     Local variables

        integer :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq,igx,igy,igz,nx,ny,nz
        integer :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

        real(8) :: idhx,idhy,idhz,a(3)
        real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

        real(8) :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .            ,jacp,jacm,jacph,jacmh,jach,jac0
     .            ,ijacip,ijacim,ijacjp,ijacjm,ijackp,ijackm
     .            ,ijacp,ijacm,ijacph,ijacmh,ijach,ijac0

        real(8) :: vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .            ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .            ,vzip,vzim,vzjp,vzjm,vzkp,vzkm

        real(8) :: bxip,bxim,bxjp,bxjm,bxkp,bxkm
     .            ,byip,byim,byjp,byjm,bykp,bykm
     .            ,bzip,bzim,bzjp,bzjm,bzkp,bzkm

        real(8) :: one=1d0,half=0.5d0

c     Begin program

c     Defaults

        nx = g_def%nxv(igr)
        ny = g_def%nyv(igr)
        nz = g_def%nzv(igr)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        igx = igr
        igy = igr
        igz = igr

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

        idhx = half/g_def%dxh(ig)
        idhy = half/g_def%dyh(jg)
        idhz = half/g_def%dzh(kg)

cc        jac  = g_def%gmetric%grid(igx)%jac(i,j,k)

cc        if (no_map) then
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
          idhx = one/g_def%dx(ig)
          im = i

          if (.not.no_map) then
            jacip  = g_def%gmetric%grid(igx)%jac(ip,j,k)
            jacim  = g_def%gmetric%grid(igx)%jac(i ,j,k)
            jacjp  = half*(g_def%gmetric%grid(igx)%jac(ip,jp,k)
     .                    +g_def%gmetric%grid(igx)%jac(i ,jp,k))
            jacjm  = half*(g_def%gmetric%grid(igx)%jac(ip,jm,k)
     .                    +g_def%gmetric%grid(igx)%jac(i ,jm,k))
            jackp  = half*(g_def%gmetric%grid(igx)%jac(ip,j,kp)
     .                    +g_def%gmetric%grid(igx)%jac(i ,j,kp))
            jackm  = half*(g_def%gmetric%grid(igx)%jac(ip,j,km)
     .                    +g_def%gmetric%grid(igx)%jac(i ,j,km))

            if (isSP2(g_def,i+1,igx)) then
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

          idhy = one/g_def%dy(jg)
          jm = j

          if (.not.no_map) then
            jacip  = half*(g_def%gmetric%grid(igx)%jac(ip,jp,k)
     .                    +g_def%gmetric%grid(igx)%jac(ip,j ,k))
            jacim  = half*(g_def%gmetric%grid(igx)%jac(im,jp,k)
     .                    +g_def%gmetric%grid(igx)%jac(im,j ,k))
            jacjp  = g_def%gmetric%grid(igx)%jac(i,jp,k)
            jacjm  = g_def%gmetric%grid(igx)%jac(i,j ,k)
            jackp  = half*(g_def%gmetric%grid(igx)%jac(i,jp,kp)
     .                    +g_def%gmetric%grid(igx)%jac(i,j ,kp))
            jackm  = half*(g_def%gmetric%grid(igx)%jac(i,jp,km)
     .                    +g_def%gmetric%grid(igx)%jac(i,j ,km))
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
          idhz = one/g_def%dz(kg)
          km = k

          if (.not.no_map) then
            jacip  = half*(g_def%gmetric%grid(igx)%jac(ip,j,kp)
     .                    +g_def%gmetric%grid(igx)%jac(ip,j,k ))
            jacim  = half*(g_def%gmetric%grid(igx)%jac(im,j,kp)
     .                    +g_def%gmetric%grid(igx)%jac(im,j,k ))
            jacjp  = half*(g_def%gmetric%grid(igx)%jac(i,jp,kp)
     .                    +g_def%gmetric%grid(igx)%jac(i,jp,k ))
            jacjm  = half*(g_def%gmetric%grid(igx)%jac(i,jm,kp)
     .                    +g_def%gmetric%grid(igx)%jac(i,jm,k ))
            jackp  = g_def%gmetric%grid(igx)%jac(i,j,kp)
            jackm  = g_def%gmetric%grid(igx)%jac(i,j,km)
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

          if (.not.no_map) then
            jacip  = g_def%gmetric%grid(igx)%jac(ip,j,k)
            jacim  = g_def%gmetric%grid(igx)%jac(im,j,k)
            jacjp  = g_def%gmetric%grid(igx)%jac(i,jp,k)
            jacjm  = g_def%gmetric%grid(igx)%jac(i,jm,k)
            jackp  = g_def%gmetric%grid(igx)%jac(i,j,kp)
            jackm  = g_def%gmetric%grid(igx)%jac(i,j,km)
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

        if (.not.no_map) then
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

c     curl_bxv_upwd
c     ###################################################################
      function curl_bxv_upwd(g_def,i,j,k,igr,vv,bb) result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of curl(bb x vv) at the
c     grid node (i,j,k). One sided derivatives are employed when half_elem=1
c     (i,i+1), half_elem=2 (j,j+1), and half_elem=3 (k,k+1).
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        integer :: i,j,k,igr
        real(8) :: cnv(3)
        real(8) :: vv(0:,0:,0:,:)
     $            ,bb(0:,0:,0:,:)
        type(grid_mg_def),pointer :: g_def

c     Local variables

        integer :: ig,jg,kg,ip,im,jp,jm,kp,km,ieq,nx,ny,nz,igx,igy,igz
        integer :: ijk,ijkg,ipjkg,imjkg,ijpkg,ijmkg,ijkpg,ijkmg

        real(8) :: idhx,idhy,idhz,a(3)
        real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm

        real(8) :: jacip,jacim,jacjp,jacjm,jackp,jackm
     .            ,jacp,jacm,jacph,jacmh,jach,jac0,jac
     .            ,ijacip,ijacim,ijacjp,ijacjm,ijackp,ijackm
     .            ,ijacp,ijacm,ijacph,ijacmh,ijach,ijac0,ijac

        real(8) :: vxip,vxim,vxjp,vxjm,vxkp,vxkm
     .            ,vyip,vyim,vyjp,vyjm,vykp,vykm
     .            ,vzip,vzim,vzjp,vzjm,vzkp,vzkm

        real(8) :: bxip,bxim,bxjp,bxjm,bxkp,bxkm
     .            ,byip,byim,byjp,byjm,bykp,bykm
     .            ,bzip,bzim,bzjp,bzjm,bzkp,bzkm

        real(8) :: one=1d0,half=0.5d0

c     Begin program

c     Defaults

        nx = g_def%nxv(igr)
        ny = g_def%nyv(igr)
        nz = g_def%nzv(igr)

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        igx = igr
        igy = igr
        igz = igr

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

        idhx = half/g_def%dxh(ig)
        idhy = half/g_def%dyh(jg)
        idhz = half/g_def%dzh(kg)

        jacip  = g_def%gmetric%grid(igr)%jac(ip,j,k)
        jacim  = g_def%gmetric%grid(igr)%jac(im,j,k)
        jacjp  = g_def%gmetric%grid(igr)%jac(i,jp,k)
        jacjm  = g_def%gmetric%grid(igr)%jac(i,jm,k)
        jackp  = g_def%gmetric%grid(igr)%jac(i,j,kp)
        jackm  = g_def%gmetric%grid(igr)%jac(i,j,km)
        jac    = g_def%gmetric%grid(igr)%jac(i,j,k)

        ijac    = one/jac
        ijacip  = one/jacip
        ijacim  = one/jacim
        ijacjp  = one/jacjp
        ijacjm  = one/jacjm
        ijackp  = one/jackp
        ijackm  = one/jackm

cc        !First component
cc        flxjp = 0.5/(jac+jacjp)*(
cc     .     (    (vv(i,j,k,2)+vv(i,jp,k,2))
cc     .      +abs(vv(i,j,k,2)+vv(i,jp,k,2)) )*bb(i,j ,k,1)
cc     .    +(    (vv(i,j,k,2)+vv(i,jp,k,2))         
cc     .      -abs(vv(i,j,k,2)+vv(i,jp,k,2)) )*bb(i,jp,k,1))
cc     .         -0.5*(vv(i,j ,k,1)*bb(i,j ,k,2)/jac
cc     .              +vv(i,jp,k,1)*bb(i,jp,k,2)/jacjp)
cc        flxjm = 0.5/(jac+jacjm)*(
cc     .     (    (vv(i,j,k,2)+vv(i,jm,k,2))
cc     .      +abs(vv(i,j,k,2)+vv(i,jm,k,2)) )*bb(i,jm,k,1)
cc     .    +(    (vv(i,j,k,2)+vv(i,jm,k,2))          
cc     .      -abs(vv(i,j,k,2)+vv(i,jm,k,2)) )*bb(i,j ,k,1))
cc     .         -0.5*(vv(i,j ,k,1)*bb(i,j ,k,2)/jac
cc     .              +vv(i,jm,k,1)*bb(i,jm,k,2)/jacjm)
cc
cc        flxkp = 0.5/(jac+jackp)*(
cc     .     (    (vv(i,j,k,3)+vv(i,j,kp,3))
cc     .      +abs(vv(i,j,k,3)+vv(i,j,kp,3)) )*bb(i,j,k ,1)
cc     .    +(    (vv(i,j,k,3)+vv(i,j,kp,3))          
cc     .      -abs(vv(i,j,k,3)+vv(i,j,kp,3)) )*bb(i,j,kp,1))
cc     .         -0.5*(vv(i,j,k ,1)*bb(i,j,k ,3)/jac
cc     .              +vv(i,j,kp,1)*bb(i,j,kp,3)/jackp)
cc        flxkm = 0.5/(jac+jackm)*(
cc     .     (    (vv(i,j,k,3)+vv(i,j,km,3))
cc     .      +abs(vv(i,j,k,3)+vv(i,j,km,3)) )*bb(i,j,km,1)
cc     .    +(    (vv(i,j,k,3)+vv(i,j,km,3))          
cc     .      -abs(vv(i,j,k,3)+vv(i,j,km,3)) )*bb(i,j,k ,1))
cc     .         -0.5*(vv(i,j,k ,1)*bb(i,j,k ,3)/jac
cc     .              +vv(i,j,km,1)*bb(i,j,km,3)/jackm)
cc
cc        cnv(1) =  (flxjp-flxjm)/g_def%dyh(jg)
cc     .           +(flxkp-flxkm)/g_def%dzh(kg)
cc
cc        !Second component
cc        flxip = 0.5/(jac+jacip)*(
cc     .     (    (vv(i,j,k,1)+vv(ip,j,k,1))
cc     .      +abs(vv(i,j,k,1)+vv(ip,j,k,1)) )*bb(i ,j,k,2)
cc     .    +(    (vv(i,j,k,1)+vv(ip,j,k,1))          
cc     .      -abs(vv(i,j,k,1)+vv(ip,j,k,1)) )*bb(ip,j,k,2))
cc     .         -0.5*(vv(i ,j,k,2)*bb(i ,j,k,1)/jac
cc     .              +vv(ip,j,k,2)*bb(ip,j,k,1)/jacip)
cc
cccc        if (isSP(i,j,k,igr,igr,igr)) then
cccc          flxim = 0d0
cccc        else
cc          flxim = 0.5/(jac+jacim+SP_flsv)*(
cc     .     (    (vv(i,j,k,1)+vv(im,j,k,1))
cc     .      +abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(im,j,k,2)
cc     .    +(    (vv(i,j,k,1)+vv(im,j,k,1))          
cc     .      -abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(i ,j,k,2))
cc     .         -0.5*(vv(i ,j,k,2)*bb(i ,j,k,1)/jac
cc     .              +vv(im,j,k,2)*bb(im,j,k,1)/jacim)
cccc        endif
cc
cc        flxkp = 0.5/(jac+jackp)*(
cc     .     (    (vv(i,j,k,3)+vv(i,j,kp,3))
cc     .      +abs(vv(i,j,k,3)+vv(i,j,kp,3)) )*bb(i,j,k ,2)
cc     .    +(    (vv(i,j,k,3)+vv(i,j,kp,3))          
cc     .      -abs(vv(i,j,k,3)+vv(i,j,kp,3)) )*bb(i,j,kp,2))
cc     .         -0.5*(vv(i,j,k ,2)*bb(i,j,k ,3)/jac
cc     .              +vv(i,j,kp,2)*bb(i,j,kp,3)/jackp)
cc        flxkm = 0.5/(jac+jackm)*(
cc     .     (    (vv(i,j,k,3)+vv(i,j,km,3))
cc     .      +abs(vv(i,j,k,3)+vv(i,j,km,3)) )*bb(i,j,km,2)
cc     .    +(    (vv(i,j,k,3)+vv(i,j,km,3))          
cc     .      -abs(vv(i,j,k,3)+vv(i,j,km,3)) )*bb(i,j,k ,2))
cc     .         -0.5*(vv(i,j,k ,2)*bb(i,j,k ,3)/jac
cc     .              +vv(i,j,km,2)*bb(i,j,km,3)/jackm)
cc
cc        cnv(2) =  (flxip-flxim)/g_def%dxh(ig)
cc     .           +(flxkp-flxkm)/g_def%dzh(kg)
cc
cc        !Third component
cc        flxip = 0.5/(jac+jacip)*(
cc     .     (    (vv(i,j,k,1)+vv(ip,j,k,1))
cc     .      +abs(vv(i,j,k,1)+vv(ip,j,k,1)) )*bb(i ,j,k,3)
cc     .    +(    (vv(i,j,k,1)+vv(ip,j,k,1))          
cc     .      -abs(vv(i,j,k,1)+vv(ip,j,k,1)) )*bb(ip,j,k,3))
cc     .         -0.5*(vv(i ,j,k,3)*bb(i ,j,k,1)/jac
cc     .              +vv(ip,j,k,3)*bb(ip,j,k,1)/jacip)
cc
cccc        if (isSP(i,j,k,igr,igr,igr)) then
cccc          flxim = 0d0
cccc        else
cc          flxim = 0.5/(jac+jacim+SP_flsv)*(
cc     .     (    (vv(i,j,k,1)+vv(im,j,k,1))
cc     .      +abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(im,j,k,3)
cc     .    +(    (vv(i,j,k,1)+vv(im,j,k,1))          
cc     .      -abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(i ,j,k,3))
cc     .         -0.5*(vv(i ,j,k,3)*bb(i ,j,k,1)/jac
cc     .              +vv(im,j,k,3)*bb(im,j,k,1)/jacim)
cccc        endif
cc
cc        flxjp = 0.5/(jac+jacjp)*(
cc     .     (    (vv(i,j,k,2)+vv(i,jp,k,2))
cc     .      +abs(vv(i,j,k,2)+vv(i,jp,k,2)) )*bb(i,j ,k,3)
cc     .    +(    (vv(i,j,k,2)+vv(i,jp,k,2))          
cc     .      -abs(vv(i,j,k,2)+vv(i,jp,k,2)) )*bb(i,jp,k,3))
cc     .         -0.5*(vv(i,j ,k,3)*bb(i,j ,k,2)/jac
cc     .              +vv(i,jp,k,3)*bb(i,jp,k,2)/jacjp)
cc        flxjm = 0.5/(jac+jacjm)*(
cc     .     (    (vv(i,j,k,2)+vv(i,jm,k,2))
cc     .      +abs(vv(i,j,k,2)+vv(i,jm,k,2)) )*bb(i,jm,k,3)
cc     .    +(    (vv(i,j,k,2)+vv(i,jm,k,2))          
cc     .      -abs(vv(i,j,k,2)+vv(i,jm,k,2)) )*bb(i,j ,k,3))
cc     .         -0.5*(vv(i,j ,k,3)*bb(i,j ,k,2)/jac
cc     .              +vv(i,jm,k,3)*bb(i,jm,k,2)/jacjm)
cc
cc        cnv(3) =  (flxip-flxim)/g_def%dxh(ig)
cc     .           +(flxjp-flxjm)/g_def%dyh(jg)
cc
cc ##################################################################
        !First component
        flxjp = 0.5/(jac+jacjp)*(
     .     (    (vv(i,j,k,2)+vv(i,jp,k,2))
     .      +abs(vv(i,j,k,2)+vv(i,jp,k,2)) )*bb(i,j ,k,1)
     .    +(    (vv(i,j,k,2)+vv(i,jp,k,2))         
     .      -abs(vv(i,j,k,2)+vv(i,jp,k,2)) )*bb(i,jp,k,1))
     .         -0.5*(vv(i,j ,k,1)*bb(i,j ,k,2)*ijac
     .              +vv(i,jp,k,1)*bb(i,jp,k,2)*ijacjp)
        flxjm = 0.5/(jac+jacjm)*(
     .     (    (vv(i,j,k,2)+vv(i,jm,k,2))
     .      +abs(vv(i,j,k,2)+vv(i,jm,k,2)) )*bb(i,jm,k,1)
     .    +(    (vv(i,j,k,2)+vv(i,jm,k,2))          
     .      -abs(vv(i,j,k,2)+vv(i,jm,k,2)) )*bb(i,j ,k,1))
     .         -0.5*(vv(i,j ,k,1)*bb(i,j ,k,2)*ijac
     .              +vv(i,jm,k,1)*bb(i,jm,k,2)*ijacjm)

        flxkp = 0.5/(jac+jackp)*(
     .     (    (vv(i,j,k,3)+vv(i,j,kp,3))
     .      +abs(vv(i,j,k,3)+vv(i,j,kp,3)) )*bb(i,j,k ,1)
     .    +(    (vv(i,j,k,3)+vv(i,j,kp,3))          
     .      -abs(vv(i,j,k,3)+vv(i,j,kp,3)) )*bb(i,j,kp,1))
     .         -0.5*(vv(i,j,k ,1)*bb(i,j,k ,3)*ijac
     .              +vv(i,j,kp,1)*bb(i,j,kp,3)*ijackp)
        flxkm = 0.5/(jac+jackm)*(
     .     (    (vv(i,j,k,3)+vv(i,j,km,3))
     .      +abs(vv(i,j,k,3)+vv(i,j,km,3)) )*bb(i,j,km,1)
     .    +(    (vv(i,j,k,3)+vv(i,j,km,3))          
     .      -abs(vv(i,j,k,3)+vv(i,j,km,3)) )*bb(i,j,k ,1))
     .         -0.5*(vv(i,j,k ,1)*bb(i,j,k ,3)*ijac
     .              +vv(i,j,km,1)*bb(i,j,km,3)*ijackm)

        cnv(1) =  (flxjp-flxjm)/g_def%dyh(jg)
     .           +(flxkp-flxkm)/g_def%dzh(kg)

        !Second component
        flxip = 0.5/(jac+jacip)*(
     .     (    (vv(i,j,k,1)+vv(ip,j,k,1))
     .      +abs(vv(i,j,k,1)+vv(ip,j,k,1)) )*bb(i ,j,k,2)
     .    +(    (vv(i,j,k,1)+vv(ip,j,k,1))          
     .      -abs(vv(i,j,k,1)+vv(ip,j,k,1)) )*bb(ip,j,k,2))
     .         -0.5*(vv(i ,j,k,2)*bb(i ,j,k,1)*ijac
     .              +vv(ip,j,k,2)*bb(ip,j,k,1)*ijacip)

        if (isSP(g_def,i,j,k,igr,igr,igr)) then
cc          flxim = 0d0
cc          flxim = 0.25*(
cc     .     (    (vv(i,j,k,1)*ijac+vv(im,j,k,1)*ijacim)
cc     .      +abs(vv(i,j,k,1)*ijac+vv(im,j,k,1)*ijacim) )*bb(im,j,k,2)
cc     .    +(    (vv(i,j,k,1)*ijac+vv(im,j,k,1)*ijacim)          
cc     .      -abs(vv(i,j,k,1)*ijac+vv(im,j,k,1)*ijacim) )*bb(i ,j,k,2))
cc     .         -0.5*(vv(i ,j,k,2)*bb(i ,j,k,1)*ijac
cc     .              +vv(im,j,k,2)*bb(im,j,k,1)*ijacim)
          flxim = 0.5*(
     .     (    (vv(i,j,k,1))
     .      +abs(vv(i,j,k,1)) )*bb(im,j,k,2)*ijacim
     .    +(    (vv(i,j,k,1))          
     .      -abs(vv(i,j,k,1)) )*bb(i ,j,k,2)*ijac)
     .         -0.5*(vv(i ,j,k,2)*bb(i ,j,k,1)*ijac
     .              +vv(im,j,k,2)*bb(im,j,k,1)*ijacim)
cc          flxim = -0.5*(vv(i ,j,k,2)*bb(i ,j,k,1)*ijac
cc     .                 +vv(im,j,k,2)*bb(im,j,k,1)*ijacim)
        else
          flxim = 0.5/(jac+jacim)*(
     .     (    (vv(i,j,k,1)+vv(im,j,k,1))
     .      +abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(im,j,k,2)
     .    +(    (vv(i,j,k,1)+vv(im,j,k,1))          
     .      -abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(i ,j,k,2))
     .         -0.5*(vv(i ,j,k,2)*bb(i ,j,k,1)*ijac
     .              +vv(im,j,k,2)*bb(im,j,k,1)*ijacim)
        endif

        flxkp = 0.5/(jac+jackp)*(
     .     (    (vv(i,j,k,3)+vv(i,j,kp,3))
     .      +abs(vv(i,j,k,3)+vv(i,j,kp,3)) )*bb(i,j,k ,2)
     .    +(    (vv(i,j,k,3)+vv(i,j,kp,3))          
     .      -abs(vv(i,j,k,3)+vv(i,j,kp,3)) )*bb(i,j,kp,2))
     .         -0.5*(vv(i,j,k ,2)*bb(i,j,k ,3)*ijac
     .              +vv(i,j,kp,2)*bb(i,j,kp,3)*ijackp)
        flxkm = 0.5/(jac+jackm)*(
     .     (    (vv(i,j,k,3)+vv(i,j,km,3))
     .      +abs(vv(i,j,k,3)+vv(i,j,km,3)) )*bb(i,j,km,2)
     .    +(    (vv(i,j,k,3)+vv(i,j,km,3))          
     .      -abs(vv(i,j,k,3)+vv(i,j,km,3)) )*bb(i,j,k ,2))
     .         -0.5*(vv(i,j,k ,2)*bb(i,j,k ,3)*ijac
     .              +vv(i,j,km,2)*bb(i,j,km,3)*ijackm)

        cnv(2) =  (flxip-flxim)/g_def%dxh(ig)
     .           +(flxkp-flxkm)/g_def%dzh(kg)

        !Third component
        flxip = 0.5/(jac+jacip)*(
     .     (    (vv(i,j,k,1)+vv(ip,j,k,1))
     .      +abs(vv(i,j,k,1)+vv(ip,j,k,1)) )*bb(i ,j,k,3)
     .    +(    (vv(i,j,k,1)+vv(ip,j,k,1))          
     .      -abs(vv(i,j,k,1)+vv(ip,j,k,1)) )*bb(ip,j,k,3))
     .         -0.5*(vv(i ,j,k,3)*bb(i ,j,k,1)*ijac
     .              +vv(ip,j,k,3)*bb(ip,j,k,1)*ijacip)

        if (isSP(g_def,i,j,k,igr,igr,igr)) then
cc          flxim = 0d0
cc          flxim = 0.25*(
cc     .     (    (vv(i,j,k,1)+vv(im,j,k,1))
cc     .      +abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(im,j,k,3)*ijacim
cc     .    +(    (vv(i,j,k,1)+vv(im,j,k,1))          
cc     .      -abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(i ,j,k,3)*ijac)
          flxim = 0.5*(
     .     (    (vv(i,j,k,1))
     .      +abs(vv(i,j,k,1)) )*bb(im,j,k,3)*ijacim
     .    +(    (vv(i,j,k,1))          
     .      -abs(vv(i,j,k,1)) )*bb(i ,j,k,3)*ijac)
     .         -0.5*(vv(i ,j,k,3)*bb(i ,j,k,1)*ijac
     .              +vv(im,j,k,3)*bb(im,j,k,1)*ijacim)
cc          flxim = -0.5*(vv(i ,j,k,3)*bb(i ,j,k,1)*ijac
cc     .                 +vv(im,j,k,3)*bb(im,j,k,1)*ijacim)
        else
          flxim = 0.5/(jac+jacim)*(
     .     (    (vv(i,j,k,1)+vv(im,j,k,1))
     .      +abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(im,j,k,3)
     .    +(    (vv(i,j,k,1)+vv(im,j,k,1))          
     .      -abs(vv(i,j,k,1)+vv(im,j,k,1)) )*bb(i ,j,k,3))
     .         -0.5*(vv(i ,j,k,3)*bb(i ,j,k,1)*ijac
     .              +vv(im,j,k,3)*bb(im,j,k,1)*ijacim)
        endif

        flxjp = 0.5/(jac+jacjp)*(
     .     (    (vv(i,j,k,2)+vv(i,jp,k,2))
     .      +abs(vv(i,j,k,2)+vv(i,jp,k,2)) )*bb(i,j ,k,3)
     .    +(    (vv(i,j,k,2)+vv(i,jp,k,2))          
     .      -abs(vv(i,j,k,2)+vv(i,jp,k,2)) )*bb(i,jp,k,3))
     .         -0.5*(vv(i,j ,k,3)*bb(i,j ,k,2)*ijac
     .              +vv(i,jp,k,3)*bb(i,jp,k,2)*ijacjp)
        flxjm = 0.5/(jac+jacjm)*(
     .     (    (vv(i,j,k,2)+vv(i,jm,k,2))
     .      +abs(vv(i,j,k,2)+vv(i,jm,k,2)) )*bb(i,jm,k,3)
     .    +(    (vv(i,j,k,2)+vv(i,jm,k,2))          
     .      -abs(vv(i,j,k,2)+vv(i,jm,k,2)) )*bb(i,j ,k,3))
     .         -0.5*(vv(i,j ,k,3)*bb(i,j ,k,2)*ijac
     .              +vv(i,jm,k,3)*bb(i,jm,k,2)*ijacjm)

        cnv(3) =  (flxip-flxim)/g_def%dxh(ig)
     .           +(flxjp-flxjm)/g_def%dyh(jg)

      end function curl_bxv_upwd

c     curla_x_v_upwd
c     ##############################################################
      function curla_x_v_upwd(g_def,i,j,k,igr,v0,aa) result(cov)

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def
      
        integer :: i,j,k,igr

        real(8),dimension(0:,0:,0:,:) :: aa
        real(8) :: v0(3)
        
c     Local variables

        integer :: ip,im,jp,jm,kp,km,ig,jg,kg
        real(8) :: cov(3),ijac
        real(8) :: flxip,flxim,flxjp,flxjm,flxkp,flxkm,idx,idy,idz

c     Begin program

        ip = i+1
        im = i-1
        jp = j+1
        jm = j-1
        kp = k+1
        km = k-1

        ijac = 1d0/g_def%gmetric%grid(igr)%jac(i,j,k)
        
        call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)

        idx = ijac/g_def%dxh(ig)
        idy = ijac/g_def%dyh(jg)
        idz = ijac/g_def%dzh(kg)

        !1st component
        flxip = -v0(2)*0.5*(aa(ip,j,k,2)+aa(i,j,k,2))
     .          -v0(3)*0.5*(aa(ip,j,k,3)+aa(i,j,k,3))

        flxim = -v0(2)*0.5*(aa(im,j,k,2)+aa(i,j,k,2))
     .          -v0(3)*0.5*(aa(im,j,k,3)+aa(i,j,k,3))
c$$$        flxip =+0.5*(    -v0(2)
c$$$     .              +abs(-v0(2)))*aa(i ,j,k,2)
c$$$     .         +0.5*(    -v0(2)
c$$$     .              -abs(-v0(2)))*aa(ip,j,k,2)
c$$$     .         +0.5*(    -v0(3)
c$$$     .              +abs(-v0(3)))*aa(i ,j,k,3)
c$$$     .         +0.5*(    -v0(3)
c$$$     .              -abs(-v0(3)))*aa(ip,j,k,3)
c$$$        flxim =+0.5*(    -v0(2)
c$$$     .              +abs(-v0(2)))*aa(im,j,k,2)
c$$$     .         +0.5*(    -v0(2)
c$$$     .              -abs(-v0(2)))*aa(i ,j,k,2)
c$$$     .         +0.5*(    -v0(3)
c$$$     .              +abs(-v0(3)))*aa(im,j,k,3)
c$$$     .         +0.5*(    -v0(3)
c$$$     .              -abs(-v0(3)))*aa(i ,j,k,3)

        flxjp = 0.5*(    v0(2)
     .              +abs(v0(2)))*aa(i,j ,k,1)
     .         +0.5*(    v0(2)
     .              -abs(v0(2)))*aa(i,jp,k,1)
        flxjm = 0.5*(    v0(2)
     .              +abs(v0(2)))*aa(i,jm,k,1)
     .         +0.5*(    v0(2)
     .              -abs(v0(2)))*aa(i,j ,k,1)

        flxkp = 0.5*(    v0(3)
     .              +abs(v0(3)))*aa(i,j,k ,1)
     .         +0.5*(    v0(3)            
     .              -abs(v0(3)))*aa(i,j,kp,1)
        flxkm = 0.5*(    v0(3)            
     .              +abs(v0(3)))*aa(i,j,km,1)
     .         +0.5*(    v0(3)            
     .              -abs(v0(3)))*aa(i,j,k ,1)

        cov(1) =  (flxip-flxim)*idx
     .           +(flxjp-flxjm)*idy
     .           +(flxkp-flxkm)*idz

        !2nd component
        flxip = 0.5*(    v0(1)
     .              +abs(v0(1)))*aa(i ,j,k,2)
     .         +0.5*(    v0(1)        
     .              -abs(v0(1)))*aa(ip,j,k,2)
        flxim = 0.5*(    v0(1)        
     .              +abs(v0(1)))*aa(im,j,k,2)
     .         +0.5*(    v0(1)        
     .              -abs(v0(1)))*aa(i ,j,k,2)

        flxjp = -v0(1)*0.5*(aa(i,jp,k,1)+aa(i,j,k,1))
     .          -v0(3)*0.5*(aa(i,jp,k,3)+aa(i,j,k,3))

        flxjm = -v0(1)*0.5*(aa(i,jm,k,1)+aa(i,j,k,1))
     .          -v0(3)*0.5*(aa(i,jm,k,3)+aa(i,j,k,3))
c$$$        flxjp =+0.5*(    -v0(1)
c$$$     .              +abs(-v0(1)))*aa(i,j ,k,1)
c$$$     .         +0.5*(    -v0(1)          
c$$$     .              -abs(-v0(1)))*aa(i,jp,k,1)
c$$$     .         +0.5*(    -v0(3)          
c$$$     .              +abs(-v0(3)))*aa(i,j ,k,3)
c$$$     .         +0.5*(    -v0(3)          
c$$$     .              -abs(-v0(3)))*aa(i,jp,k,3)
c$$$        flxjm =+0.5*(    -v0(1)
c$$$     .              +abs(-v0(1)))*aa(i,jm,k,1)
c$$$     .         +0.5*(    -v0(1)          
c$$$     .              -abs(-v0(1)))*aa(i,j ,k,1)
c$$$     .         +0.5*(    -v0(3)          
c$$$     .              +abs(-v0(3)))*aa(i,jm,k,3)
c$$$     .         +0.5*(    -v0(3)          
c$$$     .              -abs(-v0(3)))*aa(i,j ,k,3)

        flxkp = 0.5*(    v0(3)
     .              +abs(v0(3)))*aa(i,j,k ,2)
     .         +0.5*(    v0(3)            
     .              -abs(v0(3)))*aa(i,j,kp,2)
        flxkm = 0.5*(    v0(3)            
     .              +abs(v0(3)))*aa(i,j,km,2)
     .         +0.5*(    v0(3)            
     .              -abs(v0(3)))*aa(i,j,k ,2)

        cov(2) =  (flxip-flxim)*idx
     .           +(flxjp-flxjm)*idy
     .           +(flxkp-flxkm)*idz

        !3rd component
        flxip = 0.5*(    v0(1)
     .              +abs(v0(1)))*aa(i ,j,k,3)
     .         +0.5*(    v0(1)        
     .              -abs(v0(1)))*aa(ip,j,k,3)
        flxim = 0.5*(    v0(1)        
     .              +abs(v0(1)))*aa(im,j,k,3)
     .         +0.5*(    v0(1)        
     .              -abs(v0(1)))*aa(i ,j,k,3)

        flxjp = 0.5*(    v0(2)
     .              +abs(v0(2)))*aa(i,j ,k,3)
     .         +0.5*(    v0(2)
     .              -abs(v0(2)))*aa(i,jp,k,3)
        flxjm = 0.5*(    v0(2)
     .              +abs(v0(2)))*aa(i,jm,k,3)
     .         +0.5*(    v0(2)
     .              -abs(v0(2)))*aa(i,j ,k,3)

        flxkp = -v0(1)*0.5*(aa(i,j,kp,1)+aa(i,j,k,1))
     .          -v0(2)*0.5*(aa(i,j,kp,2)+aa(i,j,k,2))
                                               
        flxkm = -v0(1)*0.5*(aa(i,j,km,1)+aa(i,j,k,1))
     .          -v0(2)*0.5*(aa(i,j,km,2)+aa(i,j,k,2))
c$$$        flxkp =+0.5*(    -v0(1)
c$$$     .              +abs(-v0(1)))*aa(i,j,k ,1)
c$$$     .         +0.5*(    -v0(1)            
c$$$     .              -abs(-v0(1)))*aa(i,j,kp,1)
c$$$     .         +0.5*(    -v0(2)            
c$$$     .              +abs(-v0(2)))*aa(i,j,k ,2)
c$$$     .         +0.5*(    -v0(2)            
c$$$     .              -abs(-v0(2)))*aa(i,j,kp,2)
c$$$        flxkm =+0.5*(    -v0(1)
c$$$     .              +abs(-v0(1)))*aa(i,j,km,1)
c$$$     .         +0.5*(    -v0(1)            
c$$$     .              -abs(-v0(1)))*aa(i,j,k ,1)
c$$$     .         +0.5*(    -v0(2)            
c$$$     .              +abs(-v0(2)))*aa(i,j,km,2)
c$$$     .         +0.5*(    -v0(2)            
c$$$     .              -abs(-v0(2)))*aa(i,j,k ,2)

        cov(3) =  (flxip-flxim)*idx
     .           +(flxjp-flxjm)*idy
     .           +(flxkp-flxkm)*idz

      end function curla_x_v_upwd

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
cc        call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)
cc
cc        idhx = 0.5/gv%gparams%dxh(ig)
cc        idhy = 0.5/gv%gparams%dyh(jg)
cc        idhz = 0.5/gv%gparams%dzh(kg)
cc
cc        jac  = g_def%gmetric%grid(igx)%jac(i,j,k)
cc
cccc        sing_point = isSP(i,j,k,igx,igy,igz)
cc
ccc     Exceptions
cc
cc        !X
cccc        igp = ig+(ip-i)
cccc        igm = ig-(i-im)
cccc
cccc        dh(1)= (gv%gparams%xx(igp)-gv%gparams%xx(igm))
cc
cc
cc        if (hex == 1) then
cc          idhx = 1./dx(ig)
cc          im = i
cc        endif
cc
cc        jacip = face_add(ip,j,k,g_def%gmetric%grid(igx)%jac,1)
cc        jacim = face_add(im,j,k,g_def%gmetric%grid(igx)%jac,1)
cc
cc        jacjp = face_add(i,jp
cc          jacip  = g_def%gmetric%grid(igx)%jac(ip,j,k)
cc          jacim  = g_def%gmetric%grid(igx)%jac(i ,j,k)
cc          jacjp  = 0.5*(g_def%gmetric%grid(igx)%jac(ip,jp,k)
cc     .                 +g_def%gmetric%grid(igx)%jac(i ,jp,k))
cc          jacjm  = 0.5*(g_def%gmetric%grid(igx)%jac(ip,jm,k)
cc     .                 +g_def%gmetric%grid(igx)%jac(i ,jm,k))
cc          jackp  = 0.5*(g_def%gmetric%grid(igx)%jac(ip,j,kp)
cc     .                 +g_def%gmetric%grid(igx)%jac(i ,j,kp))
cc          jackm  = 0.5*(g_def%gmetric%grid(igx)%jac(ip,j,km)
cc     .                 +g_def%gmetric%grid(igx)%jac(i ,j,km))
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
cc          jacip  = 0.5*(g_def%gmetric%grid(igx)%jac(ip,jp,k)
cc     .                 +g_def%gmetric%grid(igx)%jac(ip,j ,k))
cc          jacim  = 0.5*(g_def%gmetric%grid(igx)%jac(im,jp,k)
cc     .                 +g_def%gmetric%grid(igx)%jac(im,j ,k))
cc          jacjp  = g_def%gmetric%grid(igx)%jac(i,jp,k)
cc          jacjm  = g_def%gmetric%grid(igx)%jac(i,j ,k)
cc          jackp  = 0.5*(g_def%gmetric%grid(igx)%jac(i,jp,kp)
cc     .                 +g_def%gmetric%grid(igx)%jac(i,j ,kp))
cc          jackm  = 0.5*(g_def%gmetric%grid(igx)%jac(i,jp,km)
cc     .                 +g_def%gmetric%grid(igx)%jac(i,j ,km))
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
cc          jacip  = 0.5*(g_def%gmetric%grid(igx)%jac(ip,j,kp)
cc     .                 +g_def%gmetric%grid(igx)%jac(ip,j,k ))
cc          jacim  = 0.5*(g_def%gmetric%grid(igx)%jac(im,j,kp)
cc     .                 +g_def%gmetric%grid(igx)%jac(im,j,k ))
cc          jacjp  = 0.5*(g_def%gmetric%grid(igx)%jac(i,jp,kp)
cc     .                 +g_def%gmetric%grid(igx)%jac(i,jp,k ))
cc          jacjm  = 0.5*(g_def%gmetric%grid(igx)%jac(i,jm,kp)
cc     .                 +g_def%gmetric%grid(igx)%jac(i,jm,k ))
cc          jackp  = g_def%gmetric%grid(igx)%jac(i,j,kp)
cc          jackm  = g_def%gmetric%grid(igx)%jac(i,j,km)
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
cc          jacip  = g_def%gmetric%grid(igx)%jac(ip,j,k)
cc          jacim  = g_def%gmetric%grid(igx)%jac(im,j,k)
cc          jacjp  = g_def%gmetric%grid(igx)%jac(i,jp,k)
cc          jacjm  = g_def%gmetric%grid(igx)%jac(i,jm,k)
cc          jackp  = g_def%gmetric%grid(igx)%jac(i,j,kp)
cc          jackm  = g_def%gmetric%grid(igx)%jac(i,j,km)
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
cc        jacp = g_def%gmetric%grid(igx)%jac(ip,j,k)
cc        jac0 = g_def%gmetric%grid(igx)%jac(i ,j,k)
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
cc        jacp = g_def%gmetric%grid(igx)%jac(i,jp,k)
cc        jac0 = g_def%gmetric%grid(igx)%jac(i,j ,k)
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
cc        jacp = g_def%gmetric%grid(igx)%jac(i,j,kp)
cc        jac0 = g_def%gmetric%grid(igx)%jac(i,j,k )
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

c$$$c     par_diff
c$$$c     ################################################################
c$$$      function par_diff(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
c$$$     .                 ,par_order,symmetric,limited,vol) result(par2)
c$$$c     ----------------------------------------------------------------
c$$$c     Calculates div.(BB.grad) in 2D with 2nd order accuracy. Vector
c$$$c     B *must* be unitary.
c$$$c     ----------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer  :: i,j,k,nx,ny,nz,igx,igy,igz,par_order
c$$$      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1),par2
c$$$      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
c$$$     .           ,By(0:nx+1,0:ny+1,0:nz+1)
c$$$     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)
c$$$
c$$$      logical,optional :: symmetric,limited,vol
c$$$
c$$$c     Local variables
c$$$
c$$$      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
c$$$     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
c$$$     .           ,Fluxz(i-1:i,j-1:j,k-1:k)
c$$$
c$$$      real(8)  :: bx2,by2,bxby,dtx,dty,a1,a2,sgn,bb,flx0,flxp
c$$$      integer  :: ii,jj,kk,ig,jg,kg
c$$$      logical  :: limit,symm,vol_wgt
c$$$
c$$$      real(8),dimension(:,:,:),pointer :: jac
c$$$
c$$$c     Begin program
c$$$
c$$$      Fluxx=0d0
c$$$      Fluxy=0d0
c$$$      Fluxz=0d0
c$$$
c$$$      if (PRESENT(limited)) then
c$$$        limit = limited
c$$$      else
c$$$        limit = .true.
c$$$      endif
c$$$
c$$$      vol_wgt = .false.
c$$$      if (PRESENT(vol)) vol_wgt = vol
c$$$
c$$$      if (PRESENT(symmetric)) then
c$$$        symm = symmetric
c$$$      else
c$$$        symm = .false.
c$$$      endif
c$$$
c$$$      jac => g_def%gmetric%grid(igx)%jac
c$$$
c$$$c     Find fluxes
c$$$
c$$$      select case(par_order)
c$$$      case(2) !second order
c$$$        call par_flux_2nd(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
c$$$     .                 ,Fluxx,Fluxy,Fluxz,symmetric=symm,limited=limit)
c$$$      case(4) !fourth order
c$$$        call par_flux_4th(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
c$$$     .                 ,Fluxx,Fluxy,Fluxz,limited=limit)
c$$$      case default
c$$$        write (*,*) 'Order not implemented'
c$$$        stop
c$$$      end select
c$$$
c$$$c     Assemble fluxes
c$$$
c$$$      call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$      par2 = (Fluxx(i,j,k)-Fluxx(i-1,j,k))/gv%gparams%dxh(ig)
c$$$     .      +(Fluxy(i,j,k)-Fluxy(i,j-1,k))/gv%gparams%dyh(jg)
c$$$     .      +(Fluxz(i,j,k)-Fluxz(i,j,k-1))/gv%gparams%dzh(kg)
c$$$
c$$$      par2 = par2/jac(i,j,k)
c$$$
c$$$      if (vol_wgt) par2=par2*gmetric%grid(igx)%dvol(i,j,k)
c$$$
c$$$      end function par_diff
c$$$
c$$$c     par_flux_2nd
c$$$c     ################################################################
c$$$      subroutine par_flux_2nd(i,j,k,nx,ny,nz,igx,igy,igz,tmp,Bx,By,Bz
c$$$     .                       ,Fluxx,Fluxy,Fluxz,symmetric,limited)
c$$$c     ----------------------------------------------------------------
c$$$c     Calculates div.(BB.grad) in 3D with 2nd order accuracy
c$$$c     ----------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer  :: i,j,k,nx,ny,nz,igx,igy,igz
c$$$      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1)
c$$$      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
c$$$     .           ,By(0:nx+1,0:ny+1,0:nz+1)
c$$$     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)
c$$$      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
c$$$     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
c$$$     .           ,Fluxz(i-1:i,j-1:j,k-1:k)
c$$$
c$$$      logical,optional :: limited,symmetric
c$$$
c$$$c     Local variables
c$$$
c$$$      real(8)  :: bx2,by2,bxby,dtx,dty,a1,a2,sgn,bb,flx0,flxp
c$$$      integer  :: ig,jg,kg,ii,jj,kk
c$$$      logical  :: limit,symm
c$$$
c$$$      real(8),dimension(:,:,:),pointer :: jac
c$$$
c$$$c     Begin program
c$$$
c$$$      Fluxx=0d0
c$$$      Fluxy=0d0
c$$$      Fluxz=0d0
c$$$
c$$$      if (PRESENT(limited)) then
c$$$        limit = limited
c$$$      else
c$$$        limit = .true.
c$$$      endif
c$$$
c$$$      if (PRESENT(symmetric)) then
c$$$        symm = symmetric
c$$$      else
c$$$        symm = .false.
c$$$      endif
c$$$
c$$$      jac => gmetric%grid(igx)%jac
c$$$
c$$$c     Symmetric discretization
c$$$
c$$$      if (symm) then
c$$$
c$$$        do kk = k-1,k
c$$$          do jj = j-1,j
c$$$            do ii = i-1,i
c$$$              call getMGmap(gv%gparams,ii,jj,kk,igx,igy,igz,ig,jg,kg)
c$$$
c$$$              bx2 = (.25*(bx(ii  ,jj  ,kk)/jac(ii  ,jj  ,kk)
c$$$     .                   +bx(ii+1,jj  ,kk)/jac(ii+1,jj  ,kk)
c$$$     .                   +bx(ii  ,jj+1,kk)/jac(ii  ,jj+1,kk)
c$$$     .                   +bx(ii+1,jj+1,kk)/jac(ii+1,jj+1,kk)))
c$$$     .             *(.25*(bx(ii  ,jj  ,kk)
c$$$     .                   +bx(ii+1,jj  ,kk)
c$$$     .                   +bx(ii  ,jj+1,kk)
c$$$     .                   +bx(ii+1,jj+1,kk)))
c$$$              bxby= (.25*(bx(ii  ,jj  ,kk)/jac(ii  ,jj  ,kk)
c$$$     .                   +bx(ii+1,jj  ,kk)/jac(ii+1,jj  ,kk)
c$$$     .                   +bx(ii  ,jj+1,kk)/jac(ii  ,jj+1,kk)
c$$$     .                   +bx(ii+1,jj+1,kk)/jac(ii+1,jj+1,kk)))
c$$$     .             *(.25*(by(ii  ,jj  ,kk)
c$$$     .                   +by(ii+1,jj  ,kk)
c$$$     .                   +by(ii  ,jj+1,kk)
c$$$     .                   +by(ii+1,jj+1,kk)))
c$$$              by2 = (.25*(by(ii  ,jj  ,kk)/jac(ii  ,jj  ,kk)
c$$$     .                   +by(ii+1,jj  ,kk)/jac(ii+1,jj  ,kk)
c$$$     .                   +by(ii  ,jj+1,kk)/jac(ii  ,jj+1,kk)
c$$$     .                   +by(ii+1,jj+1,kk)/jac(ii+1,jj+1,kk)))
c$$$     .             *(.25*(by(ii  ,jj  ,kk)
c$$$     .                   +by(ii+1,jj  ,kk)
c$$$     .                   +by(ii  ,jj+1,kk)
c$$$     .                   +by(ii+1,jj+1,kk)))
c$$$
c$$$              dtx = .5*(tmp(ii+1,jj  ,kk)-tmp(ii,jj  ,kk)
c$$$     .                 +tmp(ii+1,jj+1,kk)-tmp(ii,jj+1,kk))
c$$$     .                 /gv%gparams%dx(ig)
c$$$              dty = .5*(tmp(ii  ,jj+1,kk)-tmp(ii  ,jj,kk)
c$$$     .                 +tmp(ii+1,jj+1,kk)-tmp(ii+1,jj,kk))
c$$$     .                 /gv%gparams%dy(jg)
c$$$
c$$$              Fluxx(ii,jj,kk)= bx2*dtx + bxby*dty
c$$$              Fluxy(ii,jj,kk)= by2*dty + bxby*dtx
c$$$              Fluxz(ii,jj,kk) = 0d0
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$
c$$$        !Average fluxes from vertices to faces
c$$$        Fluxx(:,j,k) = 0.25*(Fluxx(:,j,k  )+Fluxx(:,j-1,k  )
c$$$     .                      +Fluxx(:,j,k-1)+Fluxx(:,j-1,k-1))
c$$$
c$$$        Fluxy(i,:,k) = 0.25*(Fluxy(i,:,k  )+Fluxy(i-1,:,k  )
c$$$     .                      +Fluxy(i,:,k-1)+Fluxy(i-1,:,k-1))
c$$$
c$$$        Fluxz(i,j,:) = 0.25*(Fluxz(i,j  ,:)+Fluxz(i-1,j  ,:)
c$$$     .                      +Fluxz(i,j-1,:)+Fluxz(i-1,j-1,:))
c$$$
c$$$c     Asymmetric discretization
c$$$
c$$$      else
c$$$
c$$$        if (limit) then !Limited differences (Sharma, Hammett, JCP 227 (2007))
c$$$
c$$$          !X-component cross flux
c$$$          do ii = i-1,i
c$$$            call getMGmap(gv%gparams,ii,j,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$            bxby = 0.5*(bx(ii  ,j,k)/jac(ii  ,j,k)
c$$$     .                 *by(ii  ,j,k)
c$$$     .                 +bx(ii+1,j,k)/jac(ii+1,j,k)
c$$$     .                 *by(ii+1,j,k))
c$$$
c$$$            flx0 = mc(tmp(ii  ,j  ,k)-tmp(ii  ,j-1,k)
c$$$     .               ,tmp(ii  ,j+1,k)-tmp(ii  ,j  ,k))
c$$$            flxp = mc(tmp(ii+1,j  ,k)-tmp(ii+1,j-1,k)
c$$$     .               ,tmp(ii+1,j+1,k)-tmp(ii+1,j  ,k))
c$$$
c$$$            Fluxx(ii,j,k)= bxby/gv%gparams%dyh(jg)*mc(flx0,flxp)
c$$$          enddo
c$$$
c$$$          !Y-component cross flux 
c$$$          do jj = j-1,j
c$$$            call getMGmap(gv%gparams,i,jj,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$            bxby = 0.5*(bx(i,jj  ,k)/jac(i,jj  ,k)
c$$$     .                 *by(i,jj  ,k)
c$$$     .                 +bx(i,jj+1,k)/jac(i,jj+1,k)
c$$$     .                 *by(i,jj+1,k))
c$$$
c$$$            flx0 = mc(tmp(i  ,jj  ,k)-tmp(i-1,jj  ,k)
c$$$     .               ,tmp(i+1,jj  ,k)-tmp(i  ,jj  ,k))
c$$$            flxp = mc(tmp(i  ,jj+1,k)-tmp(i-1,jj+1,k)
c$$$     .               ,tmp(i+1,jj+1,k)-tmp(i  ,jj+1,k))
c$$$
c$$$            Fluxy(i,jj,k)= bxby/gv%gparams%dxh(ig)*mc(flx0,flxp)
c$$$          enddo
c$$$
c$$$        else !Standard central average
c$$$
c$$$          !X-component cross flux
c$$$          do ii = i-1,i
c$$$            call getMGmap(gv%gparams,ii,j,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$            bxby = 0.5*(bx(ii  ,j,k)/jac(ii  ,j,k)*by(ii  ,j,k)
c$$$     .                 +bx(ii+1,j,k)/jac(ii+1,j,k)*by(ii+1,j,k))
c$$$
c$$$            Fluxx(ii,j,k)= 0.25/gv%gparams%dyh(jg)*bxby
c$$$     .             *(tmp(ii  ,j+1,k)-tmp(ii  ,j-1,k)
c$$$     .              +tmp(ii+1,j+1,k)-tmp(ii+1,j-1,k))
c$$$          enddo
c$$$
c$$$          !Y-component cross flux 
c$$$          do jj = j-1,j
c$$$            call getMGmap(gv%gparams,i,jj,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$            bxby = 0.5*(bx(i,jj  ,k)/jac(i,jj  ,k)*by(i,jj  ,k)
c$$$     .                 +bx(i,jj+1,k)/jac(i,jj+1,k)*by(i,jj+1,k))
c$$$
c$$$            Fluxy(i,jj,k)= 0.25/gv%gparams%dxh(ig)*bxby
c$$$     .                   *(tmp(i+1,jj  ,k)-tmp(i-1,jj  ,k)
c$$$     .                    +tmp(i+1,jj+1,k)-tmp(i-1,jj+1,k))
c$$$          enddo
c$$$
c$$$        endif
c$$$
c$$$        !Add co-derivative terms: standard centered
c$$$        do ii = i-1,i
c$$$          call getMGmap(gv%gparams,ii,j,k,igx,igy,igz,ig,jg,kg)
c$$$
c$$$          bx2 = (0.5*(bx(ii+1,j,k)/jac(ii+1,j,k)
c$$$     .               +bx(ii  ,j,k)/jac(ii  ,j,k)))**2
c$$$          Fluxx(ii,j,k) = Fluxx(ii,j,k)
c$$$     .                  + bx2*(tmp(ii+1,j,k)-tmp(ii,j,k))
c$$$     .                        /gv%gparams%dx(ig)
c$$$        enddo
c$$$
c$$$        do jj = j-1,j
c$$$          call getMGmap(gv%gparams,i,jj,k,igx,igy,igz,ig,jg,kg)
c$$$ 
c$$$          by2 = (0.5*(by(i,jj+1,k)/jac(i,jj+1,k)
c$$$     .               +by(i,jj  ,k)/jac(i,jj  ,k)))**2
c$$$
c$$$          Fluxy(i,jj,k) = Fluxy(i,jj,k)
c$$$     .                  + by2*(tmp(i,jj+1,k)-tmp(i,jj,k))
c$$$     .                        /gv%gparams%dy(jg)
c$$$        enddo
c$$$
c$$$        do kk = k-1,k
c$$$          call getMGmap(gv%gparams,i,j,kk,igx,igy,igz,ig,jg,kg)
c$$$
c$$$          Fluxz(i,j,kk) = Fluxz(i,j,kk)
c$$$     .                +(0.5*(bz(i,j,kk+1)/jac(i,j,kk+1)
c$$$     .                      +bz(i,j,kk  )/jac(i,j,kk  )))**2
c$$$     .                     *(tmp(i,j,kk+1)-tmp(i,j,kk))
c$$$     .                      /gv%gparams%dz(kg)
c$$$        enddo
c$$$
c$$$cc      !Add co-derivative terms: ZIP
c$$$cc      do i = 0,nx
c$$$cc        call getMGmap(gv%gparams,i,1,1,igx,igy,igz,ig,jg,kg)
c$$$cc
c$$$cc        Fluxx(i,1:ny,1:nz) = Fluxx(i,1:ny,1:nz) 
c$$$cc     .              +bx(i+1,1:ny,1:nz)*bx(i,1:ny,1:nz)
c$$$cc     .                   *(tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
c$$$cc      enddo
c$$$cc
c$$$cc      do j = 0,ny
c$$$cc        call getMGmap(gv%gparams,1,j,1,igx,igy,igz,ig,jg,kg)
c$$$cc
c$$$cc        Fluxy(1:nx,j,1:nz) = Fluxy(1:nx,j,1:nz) 
c$$$cc     .              +by(1:nx,j+1,1:nz)*by(1:nx,j,1:nz)
c$$$cc     .                   *(tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
c$$$cc      enddo
c$$$cc
c$$$cc      do k = 0,nz
c$$$cc        call getMGmap(gv%gparams,1,1,k,igx,igy,igz,ig,jg,kg)
c$$$cc
c$$$cc        Fluxz(1:nx,1:ny,k) = Fluxz(1:nx,1:ny,k)
c$$$cc     .              +bz(1:nx,1:ny,k+1)*bz(1:nx,1:ny,k)
c$$$cc     .                   *(tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
c$$$cc      enddo
c$$$
c$$$      endif
c$$$
c$$$c     End program
c$$$
c$$$      contains
c$$$
c$$$c     minmod
c$$$c     #################################################################
c$$$      function minmod(a,b)
c$$$
c$$$      real(8) :: a,b,minmod
c$$$
c$$$      if (a*b <= 0d0) then
c$$$        minmod = 0d0
c$$$        return
c$$$      endif
c$$$
c$$$      if (a > 0d0) then
c$$$        minmod = min(a,b)
c$$$      else
c$$$        minmod = max(a,b)
c$$$      endif
c$$$
c$$$      end function minmod
c$$$
c$$$c     mc
c$$$c     #################################################################
c$$$      function mc(a,b)
c$$$
c$$$      real(8) :: a,b,mc
c$$$
c$$$      mc = minmod(2*minmod(a,b),0.5*(a+b))
c$$$
c$$$      end function mc
c$$$
c$$$      end subroutine par_flux_2nd
c$$$
c$$$ccc     iso_flux_2nd
c$$$ccc     ################################################################
c$$$cc      subroutine iso_flux_2nd(nx,ny,nz,igx,igy,igz,tmp
c$$$cc     .                       ,Fluxx,Fluxy,Fluxz,symmetric)
c$$$ccc     ----------------------------------------------------------------
c$$$ccc     Calculates div.(grad) in 3D with 2nd order accuracy
c$$$ccc     ----------------------------------------------------------------
c$$$cc
c$$$cc      implicit none
c$$$cc
c$$$ccc     Call variables
c$$$cc
c$$$cc      integer  :: nx,ny,nz,igx,igy,igz
c$$$cc      real(8)  :: tmp(0:nx+1,0:ny+1,0:nz+1)
c$$$cc      real(8)  :: Fluxx(0:nx,0:ny,0:nz)
c$$$cc     .           ,Fluxy(0:nx,0:ny,0:nz)
c$$$cc     .           ,Fluxz(0:nx,0:ny,0:nz)
c$$$cc
c$$$cc      logical,optional :: symmetric
c$$$cc
c$$$ccc     Local variables
c$$$cc
c$$$cc      real(8)  :: dtx,dty,dtz,a1,a2,sgn,bb,flx0,flxp
c$$$cc      integer  :: i,j,k,ig,jg,kg
c$$$cc      logical  :: limit,symm
c$$$cc
c$$$ccc     Begin program
c$$$cc
c$$$cc      Fluxx=0
c$$$cc      Fluxy=0
c$$$cc      Fluxz=0
c$$$cc
c$$$cc      if (PRESENT(symmetric)) then
c$$$cc        symm = symmetric
c$$$cc      else
c$$$cc        symm = .false.
c$$$cc      endif
c$$$cc
c$$$ccc     Symmetric discretization
c$$$cc
c$$$cc      if (symm) then
c$$$cc
c$$$cc        do k = 0,nz
c$$$cc          do j = 0,ny
c$$$cc            do i = 0,nx
c$$$cc
c$$$cc              call getMGmap(gv%gparams,i,j,k,igx,igy,igz,ig,jg,kg)
c$$$cc
c$$$cc              dtx = .5*(tmp(i+1,j  ,k)-tmp(i,j  ,k)
c$$$cc     .                 +tmp(i+1,j+1,k)-tmp(i,j+1,k))/dx(ig)
c$$$cc              dty = .5*(tmp(i  ,j+1,k)-tmp(i  ,j,k)
c$$$cc     .                 +tmp(i+1,j+1,k)-tmp(i+1,j,k))/dy(jg)
c$$$cc              dtz = 0d0
c$$$cc
c$$$cc              Fluxx(i,j,k)= dtx
c$$$cc              Fluxy(i,j,k)= dty
c$$$cc              Fluxz(i,j,k)= dtz
c$$$cc            enddo
c$$$cc          enddo
c$$$cc        enddo
c$$$cc
c$$$cc        !Average fluxes to faces
c$$$cc        do k = nz,1,-1
c$$$cc          do j = ny,1,-1
c$$$cc            Fluxx(:,j,k) = 0.25*(Fluxx(:,j,k  )+Fluxx(:,j-1,k  )
c$$$cc     .                          +Fluxx(:,j,k-1)+Fluxx(:,j-1,k-1))
c$$$cc          enddo
c$$$cc        enddo
c$$$cc
c$$$cc        do k = nz,1,-1
c$$$cc          do i = nx,1,-1
c$$$cc            Fluxy(i,:,k) = 0.25*(Fluxy(i,:,k  )+Fluxy(i-1,:,k  )
c$$$cc     .                          +Fluxy(i,:,k-1)+Fluxy(i-1,:,k-1))
c$$$cc          enddo
c$$$cc        enddo
c$$$cc
c$$$cc        do j = ny,1,-1
c$$$cc          do i = nx,1,-1
c$$$cc            Fluxz(i,j,:) = 0.25*(Fluxz(i,j  ,:)+Fluxz(i-1,j  ,:)
c$$$cc     .                          +Fluxz(i,j-1,:)+Fluxz(i-1,j-1,:))
c$$$cc          enddo
c$$$cc        enddo
c$$$cc
c$$$ccc     Asymmetric discretization
c$$$cc
c$$$cc      else
c$$$cc
c$$$cc        !Add co-derivative terms: standard centered
c$$$cc        do i = 0,nx
c$$$cc          call getMGmap(gv%gparams,i,1,1,igx,igy,igz,ig,jg,kg)
c$$$cc
c$$$cc          Fluxx(i,1:ny,1:nz) =
c$$$cc     .         (tmp(i+1,1:ny,1:nz)-tmp(i,1:ny,1:nz))/dx(ig)
c$$$cc        enddo
c$$$cc
c$$$cc        do j = 0,ny
c$$$cc          call getMGmap(gv%gparams,1,j,1,igx,igy,igz,ig,jg,kg)
c$$$cc
c$$$cc          Fluxy(1:nx,j,1:nz) =
c$$$cc     .         (tmp(1:nx,j+1,1:nz)-tmp(1:nx,j,1:nz))/dy(jg)
c$$$cc        enddo
c$$$cc
c$$$cc        do k = 0,nz
c$$$cc          call getMGmap(gv%gparams,1,1,k,igx,igy,igz,ig,jg,kg)
c$$$cc
c$$$cc          Fluxz(1:nx,1:ny,k) =
c$$$cc     .         (tmp(1:nx,1:ny,k+1)-tmp(1:nx,1:ny,k))/dz(kg)
c$$$cc        enddo
c$$$cc
c$$$cc      endif
c$$$cc
c$$$ccc     End program
c$$$cc
c$$$cc      end subroutine iso_flux_2nd
c$$$
c$$$c     par_flux-4th
c$$$c     ################################################################
c$$$      subroutine par_flux_4th(i,j,k,nx,ny,nz,igx,igy,igz,F,Bx,By,Bz
c$$$     .                       ,Fluxx,Fluxy,Fluxz,limited,stencil)
c$$$c     ----------------------------------------------------------------
c$$$c     Calculates div.(BB.grad) in 3D with 4th order accuracy
c$$$c     Can deal with 1D, 2D, and 3D.
c$$$c     Checks if there are enough points for the stencil in all directions
c$$$c     CANNOT handle different dx, dy, dz, (except for reduced dimensions)
c$$$c     Need information about MG to do that, hence for now dx=dy=dz=h
c$$$c
c$$$c     Coded by Natalia Krashenninikova, 8/2008
c$$$c     Modified by L. Chacon, 9/2008
c$$$c     ----------------------------------------------------------------
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer  :: i,j,k,nx,ny,nz,igx,igy,igz
c$$$      real(8)  :: F (0:nx+1,0:ny+1,0:nz+1)
c$$$      real(8)  :: Bx(0:nx+1,0:ny+1,0:nz+1)
c$$$     .           ,By(0:nx+1,0:ny+1,0:nz+1)
c$$$     .           ,Bz(0:nx+1,0:ny+1,0:nz+1)
c$$$      real(8)  :: Fluxx(i-1:i,j-1:j,k-1:k)
c$$$     .           ,Fluxy(i-1:i,j-1:j,k-1:k)
c$$$     .           ,Fluxz(i-1:i,j-1:j,k-1:k)
c$$$
c$$$      integer,optional :: stencil
c$$$      logical,optional :: limited
c$$$
c$$$c     Local variables
c$$$
c$$$      real(8)  :: dx,dy,dz
c$$$      integer  :: tgx,tgy,tgz,xi,xf,yi,yf,zi,zf,ig,jg,kg,ii,jj,kk
c$$$     .           ,stncl
c$$$      logical  :: limit
c$$$
c$$$      real(8),dimension(:,:,:),pointer :: jac
c$$$
c$$$c     Begin program
c$$$
c$$$      if (PRESENT(limited)) then
c$$$        limit = limited
c$$$      else
c$$$        limit = .true.
c$$$      endif
c$$$
c$$$      if (PRESENT(stencil)) then
c$$$        stncl = stencil
c$$$      else
c$$$        stncl = 5
c$$$      endif
c$$$
c$$$      Fluxx=0d0
c$$$      Fluxy=0d0
c$$$      Fluxz=0d0
c$$$
c$$$c     Get (constant) mesh spacings
c$$$
c$$$      call getMGmap(gv%gparams,1,1,1,igx,igy,igz,ig,jg,kg)
c$$$
c$$$      dx=gv%gparams%dxh(ig)
c$$$      dy=gv%gparams%dyh(jg)
c$$$      dz=gv%gparams%dzh(kg)
c$$$
c$$$      jac => gmetric%grid(igx)%jac
c$$$
c$$$c     Calculate fluxes
c$$$
c$$$      !X flux
c$$$      if (nx > 1) then
c$$$        do ii=i-1,i
c$$$          tgx=IOR(min(ii-2,0),max(ii+2-nx,0))
c$$$          if (ii == 0) tgx = -2
c$$$          xi=max(ii-2,0)-max(ii+3-(nx+1),0)
c$$$          xf=min(ii+3,nx+1)+max(2-ii,0)
c$$$          Fluxx(ii,j,k)=Get_Flux(dx,tgx,F(xi:xf,j,k)
c$$$     .              ,Bx(xi:xf,j,k)*Bx(xi:xf,j,k)/jac(xi:xf,j,k))
c$$$
c$$$          if (ny > 1) then
c$$$            tgx=IOR(min(ii-2,0),max(ii+2-(nx+1),0))
c$$$            tgy=IOR(min(j-2,0),max(j+2-(ny+1),0))
c$$$            if (ii == 0) tgx = -2
c$$$            xi=max(ii-2,0)-max(ii+2-(nx+1),0)
c$$$            xf=min(ii+2,nx+1)+max(2-ii,0)
c$$$            yi=max(j-2,0)-max(j+2-(ny+1),0)
c$$$            yf=min(j+2,ny+1)+max(2-j,0)
c$$$            Fluxx(ii,j,k)=Fluxx(ii,j,k)
c$$$     .                   +Get_Flux_Cross(dy,tgx,tgy,F(xi:xf,yi:yf,k)
c$$$     .                   ,Bx(xi:xf,j,k)*By(xi:xf,j,k)/jac(xi:xf,j,k))
c$$$          endif
c$$$
c$$$          if (nz > 1) then
c$$$            tgx=IOR(min(ii-2,0),max(ii+2-(nx+1),0))
c$$$            if (ii == 0) tgx = -2
c$$$            tgz=IOR(min(k-2,0),max(k+2-(nz+1),0))
c$$$            xi=max(ii-2,0)-max(ii+2-(nx+1),0)
c$$$            xf=min(ii+2,nx+1)+max(2-ii,0)
c$$$            zi=max(k-2,0)-max(k+2-(nz+1),0)
c$$$            zf=min(k+2,nz+1)+max(2-k,0)
c$$$
c$$$            Fluxx(ii,j,k)=Fluxx(ii,j,k)
c$$$     .                  +Get_Flux_Cross(dz,tgx,tgz,F(xi:xf,j,zi:zf)
c$$$     .                  ,Bx(xi:xf,j,k)*Bz(xi:xf,j,k)/jac(xi:xf,j,k))
c$$$         endif
c$$$
c$$$        enddo
c$$$      endif
c$$$
c$$$      !Y flux
c$$$      if (ny > 1) then
c$$$        do jj=j-1,j
c$$$          tgy=IOR(min(jj-2,0),max(jj+2-ny,0))
c$$$          if (jj == 0) tgy = -2
c$$$          yi=max(jj-2,0)-max(jj+3-(ny+1),0)
c$$$          yf=min(jj+3,ny+1)+max(2-jj,0)
c$$$          Fluxy(i,jj,k)=Get_Flux(dy,tgy,F(i,yi:yf,k)
c$$$     .                    ,By(i,yi:yf,k)*By(i,yi:yf,k)/jac(i,yi:yf,k))
c$$$
c$$$          if (nx > 1) then
c$$$            tgx=IOR(min(i-2,0),max(i+2-(nx+1),0))
c$$$            tgy=IOR(min(jj-2,0),max(jj+2-(ny+1),0))
c$$$            if (jj == 0) tgy = -2
c$$$            xi=max(i-2,0)-max(i+2-(nx+1),0)
c$$$            xf=min(i+2,nx+1)+max(2-i,0)
c$$$            yi=max(jj-2,0)-max(jj+2-(ny+1),0)
c$$$            yf=min(jj+2,ny+1)+max(2-jj,0)
c$$$            Fluxy(i,jj,k)=Fluxy(i,jj,k)
c$$$     .                  +Get_Flux_Cross(dx,tgy,tgx
c$$$     .                     ,Transpose(F(xi:xf,yi:yf,k))
c$$$     .                     ,By(i,yi:yf,k)*Bx(i,yi:yf,k)/jac(i,yi:yf,k))
c$$$          endif
c$$$
c$$$          if (nz > 1) then
c$$$            tgy=IOR(min(jj-2,0),max(jj+2-(ny+1),0))
c$$$            if (jj == 0) tgy = -2
c$$$            tgz=IOR(min(k-2,0),max(k+2-(nz+1),0))
c$$$            yi=max(jj-2,0)-max(jj+2-(ny+1),0)
c$$$            yf=min(jj+2,ny+1)+max(2-jj,0)
c$$$            zi=max(k-2,0)-max(k+2-(nz+1),0)
c$$$            zf=min(k+2,nz+1)+max(2-k,0)
c$$$            Fluxy(i,jj,k)=Fluxy(i,jj,k)
c$$$     .                   +Get_Flux_Cross(dz,tgy,tgz,F(i,yi:yf,zi:zf)
c$$$     .                     ,By(i,yi:yf,k)*Bz(i,yi:yf,k)/jac(i,yi:yf,k))
c$$$          endif
c$$$        enddo
c$$$      endif
c$$$
c$$$      !Z flux
c$$$      if (nz > 1) then
c$$$        do kk=k-1,k
c$$$          tgz=IOR(min(kk-2,0),max(kk+2-nz,0))
c$$$          if (kk == 0) tgz = -2
c$$$          zi=max(kk-2,0)-max(kk+2-(nz+1),0)
c$$$          zf=min(kk+2,nz+1)+max(2-kk,0)
c$$$          Fluxz(i,j,kk)=Get_Flux(dz,tgz,F(i,j,zi:zf)
c$$$     .                    ,Bz(i,j,zi:zf)*Bz(i,j,zi:zf)/jac(i,j,zi:zf))
c$$$
c$$$          if (ny > 1) then
c$$$            tgy=IOR(min(j-2,0),max(j+2-(ny+1),0))
c$$$            tgz=IOR(min(kk-2,0),max(kk+2-(nz+1),0))
c$$$            if (kk == 0) tgz = -2
c$$$            yi=max(j-2,0)-max(j+2-(ny+1),0)
c$$$            yf=min(j+2,ny+1)+max(2-j,0)
c$$$            zi=max(kk-2,0)-max(kk+2-(nz+1),0)
c$$$            zf=min(kk+2,nz+1)+max(2-kk,0)
c$$$            Fluxz(i,j,kk)=Fluxz(i,j,kk)
c$$$     .                  +Get_Flux_Cross(dy,tgz,tgy
c$$$     .                     ,Transpose(F(i,yi:yf,zi:zf))
c$$$     .                     ,Bz(i,j,zi:zf)*By(i,j,zi:zf)/jac(i,j,zi:zf))
c$$$          endif
c$$$
c$$$          if (nx > 1) then
c$$$            tgx=IOR(min(i-2,0),max(i+2-(nx+1),0))
c$$$            tgz=IOR(min(kk-2,0),max(kk+2-(nz+1),0))
c$$$            if (kk == 0) tgz = -2
c$$$            xi=max(i-2,0)-max(i+2-(nx+1),0)
c$$$            xf=min(i+2,nx+1)+max(2-i,0)
c$$$            zi=max(kk-2,0)-max(kk+2-(nz+1),0)
c$$$            zf=min(kk+2,nz+1)+max(2-kk,0)
c$$$            Fluxz(i,j,kk)=Fluxz(i,j,kk)
c$$$     .                  +Get_Flux_Cross(dx,tgz,tgx
c$$$     .                     ,Transpose(F(xi:xf,j,zi:zf))
c$$$     .                     ,Bz(i,j,zi:zf)*Bx(i,j,zi:zf)/jac(i,j,zi:zf))
c$$$          endif
c$$$        enddo
c$$$      endif
c$$$
c$$$      contains
c$$$
c$$$c     ################################################################
c$$$      function Get_Flux_Cross(h,tagx,tagy,F,B2) result(val)
c$$$c     ----------------------------------------------------------------
c$$$c     Calculates the flux for d^2/dx dy part of div.(BB.grad) operator 
c$$$c     with 4th order accuracy at a point.
c$$$c     Can handle boundary points
c$$$c     ----------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer  :: tagx,tagy
c$$$      real(8)  :: F(5,5),B2(5),val,h
c$$$
c$$$c     Local variables
c$$$      integer  :: a(-1:1,5),cm(4),cl(5),cr(5)
c$$$      real(8)  :: twelve=12d0,bxby,flx0,flx1,flx
c$$$      logical  :: chk0,chk1,chk2
c$$$
c$$$      val=0d0
c$$$
c$$$      a(-1,:)=(/-3,-10,18,-6,1/)
c$$$      a(0,:)=(/1,-8,0,8,-1/)
c$$$      a(1,:)=(/-1,6,-18,10,3/)
c$$$
c$$$      cm=(/-1,7,7,-1/)
c$$$      cl=(/2,17,-11,5,-1/)
c$$$      cr=(/-1,5,-11,17,2/)
c$$$
c$$$      if((tagx+2)*(tagx+1)*tagx*(tagx-1)+(tagy+1)*tagy*(tagy-1)==0)then
c$$$         select case(tagx)
c$$$         case (-2)
c$$$            val=sum(cl*B2*matmul(F,a(tagy,:))/twelve/h)/twelve
c$$$         case (-1)
c$$$            val=sum(cm*B2(1:4)*matmul(F(1:4,:),a(tagy,:))
c$$$     .           /twelve/h)/twelve
c$$$         case (0)
c$$$            val=sum(cm*B2(2:5)*matmul(F(2:5,:),a(tagy,:))
c$$$     .           /twelve/h)/twelve
c$$$         case (1)
c$$$            val=sum(cr*B2*matmul(F,a(tagy,:))/twelve/h)/twelve
c$$$         end select
c$$$      else
c$$$         write(*,*)"Error in Get_Flux_Cross, wrong tags"
c$$$         write(*,*)"tagx=",tagx,"tagy=",tagy
c$$$         write(*,*)"Choices are -1, 0, or 1 for left boundary, "
c$$$         write(*,*)"interior or right boundary points"
c$$$         stop
c$$$      endif
c$$$
c$$$      !Check for monotonicity (Sharma, Hammett, JCP 227 (2007))
c$$$      if (limit) then
c$$$        bxby = 0.5*(B2(4)+B2(3))
c$$$
c$$$        flx0 = mc(F(3,3)-F(3,2),F(3,4)-F(3,3),chk0)
c$$$        flx1 = mc(F(4,3)-F(4,2),F(4,4)-F(4,3),chk1)
c$$$        flx  = mc(flx0         ,flx1         ,chk2)
c$$$
c$$$cc        write (*,*) 'here',bxby,chk0,chk1,chk2
c$$$      
c$$$        if (chk0 .or. chk1 .or. chk2) then
c$$$          val = bxby/h*flx
c$$$cc          write (*,*) 'here'
c$$$cc          val = 0d0
c$$$        endif
c$$$      endif
c$$$
c$$$      end function Get_Flux_Cross
c$$$
c$$$c     ################################################################
c$$$      function Get_Flux(h,tag,F,B2) result (val)
c$$$c     ----------------------------------------------------------------
c$$$c     Calculates the flux for d^2/dx^2 part of div.(BB.grad) operator 
c$$$c     with 4th order accuracy at a point.
c$$$c     Can handle boundary points
c$$$c     ----------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer  :: tag
c$$$      real(8)  :: F(6),B2(6),val,h
c$$$
c$$$c     Local variables
c$$$cc      integer  :: a(5,6)
c$$$      real(8)  :: sv2=72d0,shty=720d0
c$$$      integer  :: a(0:5,6),cm(4),cl(5),cr(5)
c$$$      real(8)  :: twelve=12d0,sixty=60d0
c$$$
c$$$      select case(stncl)
c$$$      case(5)
c$$$
c$$$        val=0
c$$$        a=0
c$$$        a(1,2:5)=(/0,1,4,-5/)
c$$$        a(2,2:5)=(/11,-54,33,10/)
c$$$        a(3,2:5)=(/-10,-33,54,-11/)
c$$$        a(4,2:5)=(/5,-4,-1,0/)
c$$$
c$$$        select case(tag)
c$$$        case (-2)
c$$$           a(1,:)=(/-36, -185, 400, -230, 60, -9/)
c$$$           a(2,:)=(/-610, -290, 1770, -1340, 560, -90/)
c$$$           a(3,:)=(/116, 840, -1620, 970, -360, 54/)
c$$$           a(4,:)=(/-22, -430, 710, -360, 120, -18/)
c$$$           a(5,:)=(/12, 65, -120, 60, -20, 3/)
c$$$           val=sum(B2(1:5)*matmul(a(1:5,:),F)/h)/shty
c$$$           return
c$$$        case (-1)
c$$$           val=sum(B2(1:4)*matmul(a(1:4,2:5),F(1:4))/h)/sv2
c$$$           return
c$$$        case (0)
c$$$           val=sum(B2(2:5)*matmul(a(1:4,2:5),F(2:5))/h)/sv2
c$$$           return
c$$$        case (1)
c$$$           val=sum(B2(3:6)*matmul(a(1:4,2:5),F(3:6))/h)/sv2
c$$$           return
c$$$        case (2)
c$$$           a(1,:)=(/-3, 20, -60, 120, -65, -12/)
c$$$           a(2,:)=(/18, -120, 360, -710, 430, 22/)
c$$$           a(3,:)=(/-54, 360, -970, 1620, -840, -116/)
c$$$           a(4,:)=(/90, -560, 1340, -1770, 290, 610/)
c$$$           a(5,:)=(/9, -60, 230, -400, 185, 36/)
c$$$           val=sum(B2(2:6)*matmul(a(1:5,:),F)/h)/shty
c$$$           return
c$$$        case default
c$$$           write(*,*)"Error in Get_Flux, wrong tag=",tag
c$$$           write(*,*)"Choices are -2, -1, 0, 1 or 2 for left boundary, "
c$$$           write(*,*)"interior or right boundary points"
c$$$           return
c$$$        end select
c$$$      case(7)
c$$$        a(0,:)=(/-137,300,-300,200,-75,12/)
c$$$        a(1,:)=(/-12,-65,120,-60,20,-3/)
c$$$        a(2,:)=(/3,-30,-20,60,-15,2/)
c$$$        a(3,:)=(/-2,15,-60,20,30,-3/)
c$$$        a(4,:)=(/3,-20,60,-120,65,12/)
c$$$        a(5,:)=(/-12,75,-200,300,-300,137/)
c$$$
c$$$        cm=(/-1,7,7,-1/)
c$$$        cl=(/2,17,-11,5,-1/)
c$$$        cr=(/-1,5,-11,17,2/)
c$$$
c$$$        select case(tag)
c$$$        case (-2)
c$$$           val=sum(cl*B2(1:5)*matmul(a(0:4,:),F)/sixty/h)/twelve
c$$$           return
c$$$        case (-1)
c$$$           val=sum(cm*B2(1:4)*matmul(a(0:3,:),F)/sixty/h)/twelve
c$$$           return
c$$$        case (0)
c$$$           val=sum(cm*B2(2:5)*matmul(a(1:4,:),F)/sixty/h)/twelve
c$$$           return
c$$$        case (1)
c$$$           val=sum(cm*B2(3:6)*matmul(a(2:5,:),F)/sixty/h)/twelve
c$$$           return
c$$$        case (2)
c$$$           val=sum(cr*B2(2:6)*matmul(a(1:5,:),F)/sixty/h)/twelve
c$$$           return
c$$$        case default
c$$$           write(*,*)"Error in Get_Flux, wrong tag=",tag
c$$$           write(*,*)"Choices are -2, -1, 0, 1 or 2 for left boundary, "
c$$$           write(*,*)"interior or right boundary points"
c$$$           return
c$$$        end select
c$$$      case default
c$$$        call pstop('Get_Flux in par_flux_4th','Stencil not available')
c$$$      end select
c$$$
c$$$      end function Get_Flux
c$$$
c$$$c     minmod
c$$$c     #################################################################
c$$$      function minmod(a,b)
c$$$
c$$$      real(8) :: a,b,minmod
c$$$
c$$$      if (a*b <= 0d0) then
c$$$        minmod = 0d0
c$$$        return
c$$$      endif
c$$$
c$$$      if (a > 0d0) then
c$$$        minmod = min(a,b)
c$$$      else
c$$$        minmod = max(a,b)
c$$$      endif
c$$$
c$$$      end function minmod
c$$$
c$$$c     mc
c$$$c     #################################################################
c$$$      function mc(a,b,chk_limited)
c$$$
c$$$      real(8) :: a,b,mc,ho_flx
c$$$
c$$$      logical :: chk_limited
c$$$
c$$$      ho_flx = 0.5*(a+b)
c$$$
c$$$      mc = minmod(2*minmod(a,b),ho_flx)
c$$$
c$$$      chk_limited = (mc /= ho_flx)
c$$$
c$$$      end function mc
c$$$
c$$$      end subroutine par_flux_4th

      end module operators

c module equilibrium
c ######################################################################
      module equilibrium

        use problem_def

        use mg_solver

        real(8) :: dlambda,rshear,vparflow,vperflow,E0(3),B0(3),M0(3)

        real(8) :: eq_params(6)

        character(3) :: efit_coords='tor'
        
        !Flux-surface averages
        real(8) :: bzz_avg,etab2_avg,ssmax,x0,y0,z0,thmax,iR2_avg
     $            ,bp2_avg,b2j_avg

        logical :: have_jparB=.false.,gstor_I_drive=.true.,closed_orb

        !Perturbations
        integer :: nh1,nh2,nh3,npw1=1,npw2=1,npw3=1
        logical :: odd(3),random

        real(8) :: br_pert_bc,br_pert_phase,br_pert_freq
     .            ,br_pert_bc_old,br_pert_phase_old

        logical :: irrot_br_pert=.false.
     .            ,eig_is_T=.true.

        real(8),dimension(20) :: pert

        type(mg_array),target :: gpsi0,gbz0

cc        real(8),pointer,dimension(:,:,:) :: jpar_B => null()

      contains

c     dump_1d_eq
c     ###############################################################
      subroutine dump_1d_eq(r,bth0,bz0,p0,n0,vr0,eta0)

c     ---------------------------------------------------------------
c     Dump 1d equilibrium for eigen analysis
c     ---------------------------------------------------------------

      implicit none

      real(8),dimension(:) :: r,bth0,bz0,p0
      real(8),optional,dimension(:) :: n0,vr0,eta0

      integer :: nr
      real(8) :: n_r=0d0,vr_r=0d0,eta_r=0d0

      if(PRESENT(n0)  ) n_r   = 1d0
      if(PRESENT(vr0) ) vr_r  = 1d0
      if(PRESENT(eta0)) eta_r = 1d0

      if (my_rank == 0) then
        nr = size(r)
        
        open(unit=123,file=trim(equil)//'-eq.txt',status='unknown')
        write (123,*) nr
        write (123,*) eta_r,vr_r,n_r
        write (123,*) r
        write (123,*) bth0
        write (123,*) bz0
        write (123,*) p0
        if(PRESENT(eta0)) write (123,*) eta0
        if(PRESENT(vr0 )) write (123,*) vr0
        if(PRESENT(n0  )) write (123,*) n0
        close(123)
      endif

      end subroutine dump_1d_eq
      
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

        use ts_setup

        logical :: bootstrap        =.false. !Whether to calculate bootstrap current
     .            ,solenoidal       =.true.  !Whether we keep B solenoidal
     .            ,nc_eom_jxb       =.false. !Whether we use JxB instead of conserv. form
     .            ,nc_eom_gp        =.false. !Whether we use grad(p) instead of conserv. form
     .            ,sym_st           =.false. !Whether to use symmetric stress tensor
     .            ,no_eom_divpe     =.false. !Whether we include div(Pe) in EOM
     .            ,vlap_etaj        =.false. !Whether we use eta*lap(A) instead of eta*j in E_res
     .            ,subtract_E0      =.false. !Whether we substract initial E in Ohm's law
     .            ,post_divclean    =.false. !Whether to perform divergence cleaning at postproc.
     .            ,post_smooth_B    =.false. !Whether to smooth at postprocessing stage
     .            ,limit_dt_flow_cfl=.false. !Whether to limit Dt according to flow CFL
     .            ,post_filter_hel  =.false. !Whether to filter helical simulations

        logical,private :: symm_pi

        INTERFACE res
          module procedure res_mesh
        END INTERFACE

        INTERFACE chi_perp
          module procedure chi_perp_mesh
        END INTERFACE

        INTERFACE chi_para
          module procedure chi_par_mesh
        END INTERFACE

        INTERFACE vis
          module procedure vis_mesh,vis_ijk
        END INTERFACE

      contains

c     form_ve
c     ###################################################################
      function form_ve(v,j,rho) result(ve)

c     -------------------------------------------------------------------
c     Electron velocity ve=v-di*j/rho
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: v(:,:,:,:)
     .            ,j(:,:,:,:)
     .            ,rho(:,:,:)
        real(8) :: ve(0:size(v,1)-1
     .               ,0:size(v,2)-1
     .               ,0:size(v,3)-1
     .               ,3)

c     Local variables

c     Begin program

        where (rho /= 0d0) 
          ve(:,:,:,1) = v(:,:,:,1)-di*j(:,:,:,1)/rho
          ve(:,:,:,2) = v(:,:,:,2)-di*j(:,:,:,2)/rho
          ve(:,:,:,3) = v(:,:,:,3)-di*j(:,:,:,3)/rho
cc          ve(:,:,:,1) = -di*j(:,:,:,1)/rho
cc          ve(:,:,:,2) = -di*j(:,:,:,2)/rho
cc          ve(:,:,:,3) = -di*j(:,:,:,3)/rho
        end where

c     End program

      end function form_ve

c$$$c     form_jfake
c$$$c     ###################################################################
c$$$      function form_jfake(g_def,igr,dv,div_pi,rho) result(jf)
c$$$
c$$$c     -------------------------------------------------------------------
c$$$c     Fake current: jf=-dt*curl(curl(E_fake))=dt*vlap(E_fake)
c$$$c     with
c$$$c          E_fake = di*(dv/dt + div_pi(v))
c$$$c     -------------------------------------------------------------------
c$$$
c$$$        implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$        type(grid_mg_def),pointer :: g_def
c$$$
c$$$        integer :: igr
c$$$        real(8) :: dv    (0:,0:,0:,:)
c$$$     .            ,div_pi(0:,0:,0:,:)
c$$$     .            ,rho   (0:,0:,0:)
c$$$     .            ,jf(0:size(dv,1)-1
c$$$     .               ,0:size(dv,2)-1
c$$$     .               ,0:size(dv,3)-1
c$$$     .               ,  size(dv,4)  )
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$        dv = di*dv
c$$$cc        dv(:,:,:,1) = di*(dv(:,:,:,1) + dt*div_pi(:,:,:,1)/rho)
c$$$cc        dv(:,:,:,2) = di*(dv(:,:,:,2) + dt*div_pi(:,:,:,2)/rho)
c$$$cc        dv(:,:,:,3) = di*(dv(:,:,:,3) + dt*div_pi(:,:,:,3)/rho)
c$$$
c$$$        jf = veclap(g_def,igr,dv)
c$$$
c$$$c     End program
c$$$
c$$$      end function form_jfake

c     viscous_heat_src
c     ###################################################################
      function viscous_heat_src(g_def,i,j,k,igr,visc,vcnv1,vcnv2,symm)
     .         result(vhs)

c     -------------------------------------------------------------------
c     Finds viscous heat source  -Pi(v1):nabla(v2)
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igr
        real(8) :: vcnv1(0:,0:,0:,:)
     .            ,vcnv2(0:,0:,0:,:)
     .            ,visc (0:,0:,0:)
     .            ,vhs
        
        logical :: symm

c     Local variables

        integer :: nx,ny,nz
        real(8) :: nabla_v(3,3),pi_tnsr(3,3)

c     Begin program

        nx = g_def%nxv(igr)
        ny = g_def%nyv(igr)
        nz = g_def%nzv(igr)

c     Viscous src: obtain cnv tensor pi, and dot with mixed tensor nabla_v (T_l^m)

        if (visc(i,j,k) > 0d0) then
          pi_tnsr=EOM_Pi_ijk(g_def,i,j,k,nx,ny,nz,igr,vcnv1,visc,symm)

          nabla_v = fnabla_v(g_def,i,j,k,igr,vcnv2,0,1,3)

          !Make nabla_v covariant on both sides
          nabla_v = matmul(nabla_v
     .                    ,g_def%gmetric%grid(igr)%gsub(i,j,k,:,:))

          vhs =-tensorScalarProduct(g_def,i,j,k,igr,nabla_v,pi_tnsr)
        else
          vhs = 0d0
        endif

c     End program

      end function viscous_heat_src

c     ext_heat_src
c     ###################################################################
      function ext_heat_src(g_def,i,j,k,igr) result(ehs)

c     -------------------------------------------------------------------
c     Provides external heat source
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igr
        real(8) :: ehs

c     Local variables

        integer :: nx,ny,nz,ig,jg,kg
        real(8) :: x1,x2,x3

c     Begin program

#if defined(RFX)
        call getCurvilinearCoordinates(g_def,i,j,k,igr,igr,igr,ig,jg,kg
     .                                ,x1,x2,x3)

        ehs = ext_heat_src_prof(x1)
#else
        ehs = 0d0
#endif

c     End program

      end function ext_heat_src

c     ext_heat_src_prof
c     ###################################################################
      function ext_heat_src_prof(rr) result(ehsp)

c     -------------------------------------------------------------------
c     Provides external heat source
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        real(8) :: ehsp,rr

c     Local variables

c     Begin program

! Daniele, 12 Nov. 2014
!        ehsp = 1d-6
!        ehsp = chi*exp(-(rr/0.5)**2)
! Luis, June 29, 2017
!        select case(equil)
!        case('ppnch','ppnsl','ppnst','ppn3d','p3nsl')
!          ehsp = 1d-6
!        case default
          ehsp = 0d0
!        end select

c     End program

      end function ext_heat_src_prof

c     ext_mom_src
c     ###################################################################
      function ext_mom_src(g_def,i,j,k,igr) result(ems)

c     -------------------------------------------------------------------
c     Provides external momentum source for selected equilibria, in
c     contravariant representation
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igr
        real(8) :: ems(3)

c     Local variables

        integer :: nx,ny,nz,ig,jg,kg
        real(8) :: rr,kk,mm

c     Begin program

        call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)
        rr = g_def%xx(ig)

        select case(equil)
        case('ppnch','ppnsl','ppnst')
          mm = g_def%params(1)
          kk = g_def%params(2)

          ems(1) = rr*M0(1)
          ems(2) = M0(2) + kk*rr/mm*M0(3)
          ems(3) = rr*M0(3)
        case ('ppn3d','p3nsl')
          ems(1) = rr*M0(1)
          ems(2) = M0(2)
          ems(3) = rr*M0(3)
        case default
          ems = 0d0
        end select

c     End program

      end function ext_mom_src

ccc     viscous_heat_src2
ccc     ###################################################################
cc      function viscous_heat_src2(i,j,k,igr,vcnv,tnsr) result(viscous)
cc
ccc     -------------------------------------------------------------------
ccc     Finds viscous heat source  -Pi(v1):nabla(v2)
ccc     -------------------------------------------------------------------
cc
cc        implicit none
cc
ccc     Call variables
cc
cc        integer :: i,j,k,igr
cc        real(8) :: vcnv(0:,0:,0:,:)
cc     .            ,tnsr(0:,0:,0:,:,:)
cc     .            ,viscous
cc
ccc     Local variables
cc
cc        integer :: nx,ny,nz
cc        real(8) :: nabla_v(3,3)
cc
ccc     Begin program
cc
cc        nx = gv%gparams%nxv(igr)
cc        ny = gv%gparams%nyv(igr)
cc        nz = gv%gparams%nzv(igr)
cc
ccc     Viscous src: obtain cnv tensor pi, and dot with mixed tensor nabla_v (T_l^m)
cc
cc        nabla_v = fnabla_v(i,j,k,nx,ny,nz,igr,igr,igr
cc     .                    ,vcnv(:,:,:,1)
cc     .                    ,vcnv(:,:,:,2)
cc     .                    ,vcnv(:,:,:,3),0)
cc
cc        !Make nabla_v covariant on both sides
cc        nabla_v = matmul(nabla_v,gmetric%grid(igr)%gsub(i,j,k,:,:))
cc
cc        viscous =-tensorScalarProduct(i,j,k,igr,nabla_v,tnsr)
cc
ccc     End program
cc
cc      end function viscous_heat_src2

c     int_heat_src
c     ###################################################################
      function int_heat_src(g_def,i,j,k,igr,eta,nu,heta,jcnv,jcov,vcnv
     .                     ,vecnv) result(ihs)

c     -------------------------------------------------------------------
c     Provides external momentum source for selected equilibria, in
c     contravariant representation
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer :: i,j,k,igr
        real(8) :: ihs
     .            ,jcnv(0:,0:,0:,:)
     .            ,jcov(0:,0:,0:,:)
     .            ,vcnv(0:,0:,0:,:)
     .            ,vecnv(0:,0:,0:,:)
     .            ,heta(0:,0:,0:)
     .            ,eta (0:,0:,0:)
     .            ,nu  (0:,0:,0:)

c     Local variables

c     Begin program

        !Joule heating
        ihs =eta(i,j,k)*scalarProduct_ijk(g_def,i,j,k,igr,jcnv(i,j,k,:)
     .                                                   ,jcov(i,j,k,:))
        
        !Ion viscosity
        ihs =ihs + viscous_heat_src(g_def,i,j,k,igr,nu,vcnv,vcnv,sym_st)

        !Electron viscosity
        !Comment out for old TM Hall MHD performance results
        if (di > 0d0) then
          ihs = ihs
     .       +viscous_heat_src(g_def,i,j,k,igr,heta,vecnv
     .                        ,vecnv,.false.)
        endif

c     End program

      end function int_heat_src

c     fillLocalAuxVars
c     ###################################################################
      subroutine fillLocalAuxVars(igr,varray,vaux)

      implicit none
c     -------------------------------------------------------------------
c     Fills local auxiliary variables (SAMRAI unaware)
c     -------------------------------------------------------------------

c     Call variables

      integer :: igr

      type(var_array),pointer :: varray

      type(aux_array),pointer :: vaux

c     Local variables

      integer :: i,j,k,ig,jg,kg,nx,ny,nz,bcs(6,3),bc(6)

      real(8) :: bnorm

      real(8),pointer,dimension(:,:,:) :: rho=>null(),tmp=>null()
     .                                   ,jac=>null(),tmpe=>null()
      real(8),pointer,dimension(:,:,:,:) :: dum1,dum2

c     Begin program

      nx = gv%gparams%nxv(igr)
      ny = gv%gparams%nyv(igr)
      nz = gv%gparams%nzv(igr)

c     Allocate LOCAL auxiliary variables (not SAMRAI-aware)

      if (.not.associated(pcnv)) then
        allocate(pcnv    (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,vscnv   (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,bhat    (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,v_advc  (0:nx+1,0:ny+1,0:nz+1,3)
     .          ,rho_advc(0:nx+1,0:ny+1,0:nz+1)
     .          ,prs_advc(0:nx+1,0:ny+1,0:nz+1)
     .          ,tmp_advc(0:nx+1,0:ny+1,0:nz+1)
     .          ,heat_xport(0:nx+1,0:ny+1,0:nz+1)
     .          ,heat_src(0:nx+1,0:ny+1,0:nz+1)
     .          ,prs     (0:nx+1,0:ny+1,0:nz+1)
     .          ,chi_tnsr(0:nx+1,0:ny+1,0:nz+1,6))
      endif

      pcnv = 0d0
      vscnv= 0d0
      bhat = 0d0
      prs  = 0d0
      v_advc = 0d0
      rho_advc = 0d0
      prs_advc = 0d0
      tmp_advc = 0d0
      heat_xport = 0d0
      heat_src   = 0d0
      chi_tnsr   = 0d0
      
      rho => vaux%var_list(IRHO_AUX )%array
      tmp => vaux%var_list(ITMP_AUX )%array
      tmpe=> vaux%var_list(ITMPE_AUX)%array
      jac => gv%gparams%gmetric%grid(igr)%jac

c$$$      !Store B/|B| unit vector (for parallel heat conduction and/or boostrap current)
c$$$      if (chi_par /= chi .and. (.not.lagrangian)) then
c$$$        do k=0,nz+1
c$$$          do j=0,ny+1
c$$$            do i=0,nx+1
c$$$              bnorm = vectorNorm(gv%gparams,i,j,k,igr,bcnv(i,j,k,:)
c$$$     .             ,.false.)
c$$$              if (bnorm == 0d0) then
c$$$                bhat(i,j,k,:) = 0d0
c$$$              else
c$$$                bhat(i,j,k,:) = bcnv(i,j,k,:)/sqrt(bnorm)
c$$$              endif
c$$$            enddo
c$$$          enddo
c$$$        enddo
c$$$      endif

      !Momentum
      pcnv(:,:,:,1) = rho*vcnv(:,:,:,1)
      pcnv(:,:,:,2) = rho*vcnv(:,:,:,2)
      pcnv(:,:,:,3) = rho*vcnv(:,:,:,3)

      !Total pressure
      if (solve_prs) then
        prs = varray%array_var(IPRS)%array/jac
      else
        prs = rho*tmp
      endif

      !Temperature transport velocity, v_star = v - di*j/rho/a_p
      if (.not.e_cons_vstar) then
        vscnv = vcnv            !Also to recover old TM Hall performance results
      else
        vscnv(:,:,:,1) = vcnv(:,:,:,1)-(di/a_p)*jcnv(:,:,:,1)/rho
        vscnv(:,:,:,2) = vcnv(:,:,:,2)-(di/a_p)*jcnv(:,:,:,2)/rho
        vscnv(:,:,:,3) = vcnv(:,:,:,3)-(di/a_p)*jcnv(:,:,:,3)/rho
      endif

      !Advection operators
#if !defined(flux_rhs)
      bcs = gv%aux%vec_list(IVCNV)%bconds
      if (nc_eom_v) then
        v_advc= tensor_nc_advc_mesh(gv%gparams,igr,vcnv,vcnv,bcs
     .                             ,v_advect,.false.)
      else
        v_advc= tensor_advc_mesh(gv%gparams,igr,vcnv,pcnv,bcs,v_advect)
      endif

      if (solve_rho) then
        bc = varray%array_var(IRHO)%bconds
        rho_advc =flx_advec_mesh(gv%gparams,igr,vcnv,rho,advect,bc) !,sp=bcSP())
      endif

      bc = varray%array_var(IPRS)%bconds
      if (solve_prs) then
        if (use_p_classic_eq) then !Classic p evolution eq
          prs_advc = flx_advec_mesh(gv%gparams,igr,vscnv,prs,advect,bc)
     .           +(gamma-1d0)*prs*div(gv%gparams,igr,vscnv)
        else
          prs_advc=
     .  gamma     *flx_advec_mesh   (gv%gparams,igr,vscnv,prs,advect,bc)
     .+(gamma-1d0)*flx_nc_advec_mesh(gv%gparams,igr,vscnv,prs,advect,bc
     .                              ,rvrs=.true.) !Reverse: v -> -v
        endif
      else
        if (use_p_classic_eq) then !Classic T evolution eq
          tmp_advc =
     .            flx_nc_advec_mesh(gv%gparams,igr,vscnv,tmpe,advect,bc)
     .           +(gamma-1d0)*tmpe*div(gv%gparams,igr,vscnv)
        else
          tmp_advc =
     . (gamma-1.)*flx_advec_mesh   (gv%gparams,igr,vscnv,tmpe,advect,bc)
     .+(2.-gamma)*flx_nc_advec_mesh(gv%gparams,igr,vscnv,tmpe,advect,bc)
        endif
      endif
#endif

      !Heat transport term
      if (gamma > 1d0) then
#if !defined(flux_rhs)
        if (chi_par == chi) then
          if (chi /= 0d0) then
            if (use_p_diff) then !Total pressure
              heat_xport = -lap(gv%gparams,igr,prs,dff=cchi)
            elseif (solve_prs) then
              heat_xport = -lap(gv%gparams,igr,tmp,dff=cchi)
            else
              heat_xport = -lap(gv%gparams,igr,tmpe,dff=cchi)
            endif
          endif
        else
          allocate(dum1(0:nx+1,0:ny+1,0:nz+1,1)
     .            ,dum2(0:nx+1,0:ny+1,0:nz+1,1))
          dum1(:,:,:,1) = cchi_par
          dum2(:,:,:,1) = cchi
          chi_tnsr = find_diff_tnsr(gv%gparams,igr,bcnv,dum1,dum2,6)
          heat_xport =-tnsr_diff(gv%gparams,igr,tmpe
     .                          ,varray%array_var(IPRS)%bconds,chi_tnsr
     .                          ,xport_order,.true.,limit_mthd=SMART)
cc     .                          ,2,.true.,limit_mthd=SMART)
          deallocate(dum1,dum2)
        endif
      endif
#endif
      
      !Heat source term
      if ((.not.adiabatic).and.(gamma > 1d0)) then
        do k=1,nz
          do j=1,ny
            do i=1,nx
              !Internal heat source
              heat_src(i,j,k) =
     .             int_heat_src(gv%gparams,i,j,k,igr,eeta,nuu,h_eta
     .                         ,jcnv,jcov,vcnv,vecnv)

              !External heat source
              heat_src(i,j,k) = heat_src(i,j,k)
     .                        + ext_heat_src(gv%gparams,i,j,k,igr)
            enddo
          enddo
        enddo
      endif

c     End program

      end subroutine fillLocalAuxVars

c     find_diff_tnsr
c     ########################################################################
      function find_diff_tnsr(g_def,igr,bb,chipar,chiper,neq)
     .         result(tnsr)

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr,neq
      real(8) :: bb(0:,0:,0:,:)
      real(8) :: chipar(0:,0:,0:,:)
      real(8) :: chiper(0:,0:,0:,:)
      real(8) :: tnsr(0:size(bb,1)-1
     .               ,0:size(bb,2)-1
     .               ,0:size(bb,3)-1,neq)

c     Local variables

      integer :: i,j,k,ig,jg,kg,nx,ny,nz,i1,j1,ieq
      real(8) :: x1,y1,z1,bnorm
      real(8) :: bhat  (0:size(bb,1)-1
     .                 ,0:size(bb,2)-1
     .                 ,0:size(bb,3)-1,3)

c     Begin program

      if (neq /= 6) then
        call pstop("find_diff_tnsr","Wrong dimension")
      endif

      nx = g_def%nxv(igr)
      ny = g_def%nyv(igr)
      nz = g_def%nzv(igr)

      bhat = find_bhat(g_def,igr,bb,3)

c     Diffusion tensors

      if (.not.lagrangian) then
        !Find diffusion tensor
        do i1 = 1,3
          do j1 = i1,3
            if (i1 == 1) then
              ieq = j1
            else
              ieq = i1+j1
            endif
            tnsr(:,:,:,ieq)=
     .        (chipar(:,:,:,1)-chiper(:,:,:,1))
     .                        *bhat(:,:,:,i1)
     .                        *bhat(:,:,:,j1)
     .                        /g_def%gmetric%grid(igr)%jac
     .        +chiper(:,:,:,1)*g_def%gmetric%grid(igr)%gsup(:,:,:,j1,i1)
          enddo
        enddo
      else
c$$$        !Find perpendicular diffusion tensor
c$$$        if (isotropic) then
c$$$          do i1 = 1,3
c$$$            do j1 = i1,3
c$$$              if (i1 == 1) then
c$$$                ieq = j1
c$$$              else
c$$$                ieq = i1+j1
c$$$              endif
c$$$              tnsr(:,:,:,ieq)= 
c$$$     .         chiper(:,:,:,1)*g_def%gmetric%grid(igr)%gsup(:,:,:,j1,i1)
c$$$            enddo
c$$$          enddo
c$$$        else
          do i1 = 1,3
            do j1 = i1,3
              if (i1 == 1) then
                ieq = j1
              else
                ieq = i1+j1
              endif
              tnsr(:,:,:,ieq)= 
     .            chiper(:,:,:,1)
     .                     *(g_def%gmetric%grid(igr)%gsup(:,:,:,j1,i1)
     .                     -bhat(:,:,:,i1)*bhat(:,:,:,j1)
     .                     /g_def%gmetric%grid(igr)%jac )
            enddo
          enddo
c$$$        endif
      endif

c     End program

      end function find_diff_tnsr

c     find_bhat
c     ########################################################################
      function find_bhat(g_def,igr,bb,neq) result(bhat)

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr,neq
      real(8) :: bb(0:,0:,0:,:)
      real(8) :: bhat(0:size(bb,1)-1
     .               ,0:size(bb,2)-1
     .               ,0:size(bb,3)-1,neq)

c     Local variables

      real(8) :: ibnorm(0:size(bb,1)-1
     .                 ,0:size(bb,2)-1
     .                 ,0:size(bb,3)-1)

c     Begin program

      if (neq /= 3) then
        call pstop("find_hat","Wrong dimension")
      endif

      ibnorm = sqrt(vectorNorm(g_def,igr,bb,.false.))
      where (ibnorm > 0d0)
        ibnorm = 1d0/ibnorm
        bhat(:,:,:,1) = bb(:,:,:,1)*ibnorm
        bhat(:,:,:,2) = bb(:,:,:,2)*ibnorm
        bhat(:,:,:,3) = bb(:,:,:,3)*ibnorm
      elsewhere
        bhat(:,:,:,1) = 0d0
        bhat(:,:,:,2) = 0d0
        bhat(:,:,:,3) = 0d0
      end where

c     End program

      end function find_bhat

c     EOM_jxb
c     ###################################################################
      function EOM_jxb(g_def,i,j,k,nx,ny,nz,igrid,b1,b2,vol,proj_parcmp)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components of the conservative EOM force term
c     (jxb) at cell centers.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,nx,ny,nz,igrid
        real(8)    :: cnv(3)
        real(8),pointer,dimension(:,:,:,:) :: b1,b2
        logical,optional,intent(IN) :: vol,proj_parcmp

c     Local variables

        logical    :: project
        real(8)    :: bbcov(3),mag

c     Begin program

        project = .false.
        if (PRESENT(proj_parcmp)) project = proj_parcmp

        vec1 => b1
        vec2 => b2
        cnv =-div_tensor(g_def,i,j,k,nx,ny,nz,igrid,igrid,igrid
     .                  ,alt__eom(),eom_jxb_x,eom_jxb_y,eom_jxb_z
     .                  ,vol=vol)
        nullify(vec1,vec2)

        !Project out parallel force (destroys conservation)
        if (project) then
          bbcov = XformToCov_ijk(g_def,i,j,k,igrid,b1(i,j,k,:))
          mag = dot_product(bbcov,b1(i,j,k,:))
          if (mag > 0d0) then
            cnv = cnv - b1(i,j,k,:)*dot_product(bbcov,cnv)/mag
          endif
        endif

c     End program

      end function EOM_jxb

c     EOM_jxb_x
c     #############################################################
      subroutine EOM_jxb_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip,igr
        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3),gsuper(3)
     .               ,car0(3),carp(3),vec1h(3),vec2h(3)
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,pmag,coeff

c     Begin program

        igr = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(g_def%gmetric%grid(igx)%jac (ip,j,k)
     .               +g_def%gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(g_def%gmetric%grid(igx)%gsup(ip,j,k,1,:)
     .               +g_def%gmetric%grid(igx)%gsup(i ,j,k,1,:))

        if (flag /= 0) then
          jacp = g_def%gmetric%grid(igx)%jac(ip,j,k)
          jac0 = g_def%gmetric%grid(igx)%jac(i ,j,k)
        else
          jacp = jac
          jac0 = jac
        endif

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

        !Magnetic pressure
c$$$  ===============================================        
c$$$        car0 = XformToCar(g_def,i ,j,k,igx,vec1(i ,j,k,:),.false.)
c$$$        carp = XformToCar(g_def,ip,j,k,igx,vec2(ip,j,k,:),.false.)
c$$$        pmag = (carp(1)*car0(1)
c$$$     .         +car0(1)*carp(1)
c$$$     .         +carp(2)*car0(2)
c$$$     .         +car0(2)*carp(2)
c$$$     .         +carp(3)*car0(3)
c$$$     .         +car0(3)*carp(3))*0.25
c$$$  ===============================================        
        cnv0 = vec1(i ,j,k,:)
        cnvp = vec2(ip,j,k,:)

        cov0 = XformToCov(g_def,i ,j,k,igx,cnv0)
        covp = XformToCov(g_def,ip,j,k,igx,cnvp)

        cnv0 = cnv0*ijac0
        cnvp = cnvp*ijacp

        pmag = (cnvp(1)*cov0(1)
     .         +cnv0(1)*covp(1)
     .         +cnvp(2)*cov0(2)
     .         +cnv0(2)*covp(2)
     .         +cnvp(3)*cov0(3)
     .         +cnv0(3)*covp(3))*0.25
        
        !Assemble tensor
c$$$  ===============================================                
c$$$        if (isSP(g_def,ip,j,k,igx,igy,igz).and.flag /= 0) jac = SP_flsv
c$$$
c$$$        coeff = ijacp*ijac0*jac**2
c$$$
c$$$        t11 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,1)
c$$$     .               +vec1(i ,j,k,1)*vec2(ip,j,k,1))*coeff
c$$$     .         +gsuper(1)*pmag*jac)
c$$$
c$$$        t12 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,2)
c$$$     .               +vec1(i ,j,k,1)*vec2(ip,j,k,2))*coeff
c$$$     .         +gsuper(2)*pmag*jac)
c$$$
c$$$        t13 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,3)
c$$$     .               +vec1(i ,j,k,1)*vec2(ip,j,k,3))*coeff
c$$$     .         +gsuper(3)*pmag*jac)
c$$$
c$$$        if (flag /= 0) then
c$$$          ijac = 1d0/jac
c$$$
c$$$          t11 = t11*ijac
c$$$          if (.not.alt_eom) t12 = t12*ijac
c$$$          t13 = t13*ijac
c$$$        endif
c$$$  ===============================================        

c$$$  ===============================================        
c$$$        coeff = ijacp*ijac0*jac
c$$$        t11 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,1)
c$$$     .               +vec1(i ,j,k,1)*vec2(ip,j,k,1))*coeff
c$$$     .         +gsuper(1)*pmag)
c$$$
c$$$        t12 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,2)
c$$$     .               +vec1(i ,j,k,1)*vec2(ip,j,k,2))*coeff
c$$$     .         +gsuper(2)*pmag)
c$$$
c$$$        t13 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,3)
c$$$     .               +vec1(i ,j,k,1)*vec2(ip,j,k,3))*coeff
c$$$     .         +gsuper(3)*pmag)
c$$$  ===============================================        
c$$$        vec1h = 0.5*(vec1(i,j,k,:)*ijac0+vec1(ip,j,k,:)*ijacp)
c$$$        vec2h = 0.5*(vec2(i,j,k,:)*ijac0+vec2(ip,j,k,:)*ijacp)
c$$$        
c$$$        t11 = (-vec1h(1)*vec2h(1)*jac+gsuper(1)*pmag)
c$$$        t12 = (-vec1h(1)*vec2h(2)*jac+gsuper(2)*pmag)
c$$$        t13 = (-vec1h(1)*vec2h(3)*jac+gsuper(3)*pmag)
c$$$  ===============================================        
        t11 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,1)*ijac0
     .               +vec1(i ,j,k,1)*vec2(ip,j,k,1)*ijacp)
     .         +gsuper(1)*pmag)

        t12 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,2)*ijac0
     .               +vec1(i ,j,k,1)*vec2(ip,j,k,2)*ijacp)
     .         +gsuper(2)*pmag)

        t13 = (-0.5*( vec1(ip,j,k,1)*vec2(i ,j,k,3)*ijac0
     .               +vec1(i ,j,k,1)*vec2(ip,j,k,3)*ijacp)
     .         +gsuper(3)*pmag)
c$$$  ===============================================        

        if (flag == 0) then
          t11 = t11*jac
          if (alt_eom) t12 = t12*jac
          t13 = t13*jac
        endif

c     End program

      end subroutine EOM_jxb_x

c     EOM_jxb_y
c     #############################################################
      subroutine EOM_jxb_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp,igr
        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3),gsuper(3)
     .               ,car0(3),carp(3),vec1h(3),vec2h(3)
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,pmag,coeff

c     Begin program

        igr = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(g_def%gmetric%grid(igx)%jac (i,jp,k)
     .               +g_def%gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(g_def%gmetric%grid(igx)%gsup(i,jp,k,2,:)
     .               +g_def%gmetric%grid(igx)%gsup(i,j ,k,2,:))

        if (flag /= 0) then
          jacp = g_def%gmetric%grid(igx)%jac(i,jp,k)
          jac0 = g_def%gmetric%grid(igx)%jac(i,j ,k)
        else
          jacp = jac
          jac0 = jac
        endif

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

        !Magnetic pressure
c$$$  ===============================================        
c$$$        car0 = XformToCar(g_def,i,j ,k,igx,vec1(i,j ,k,:),.false.)
c$$$        carp = XformToCar(g_def,i,jp,k,igx,vec2(i,jp,k,:),.false.)
c$$$        pmag = (carp(1)*car0(1)
c$$$     .         +car0(1)*carp(1)
c$$$     .         +carp(2)*car0(2)
c$$$     .         +car0(2)*carp(2)
c$$$     .         +carp(3)*car0(3)
c$$$     .         +car0(3)*carp(3))*0.25
c$$$  ===============================================        
        cnv0 = vec1(i,j ,k,:)
        cnvp = vec2(i,jp,k,:)

        cov0 = XformToCov(g_def,i,j ,k,igx,cnv0)
        covp = XformToCov(g_def,i,jp,k,igx,cnvp)

        cnv0 = cnv0*ijac0
        cnvp = cnvp*ijacp

        pmag = (cnvp(1)*cov0(1)
     .         +cnv0(1)*covp(1)
     .         +cnvp(2)*cov0(2)
     .         +cnv0(2)*covp(2)
     .         +cnvp(3)*cov0(3)
     .         +cnv0(3)*covp(3))*0.25

        !Assemble tensor
c$$$  ===============================================                
c$$$        coeff = ijacp*ijac0*jac**2
c$$$
c$$$        t21 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,1)      
c$$$     .               +vec1(i,j ,k,2)*vec2(i,jp,k,1) )*coeff
c$$$     .         +gsuper(1)*pmag*jac)
c$$$
c$$$        t22 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,2)      
c$$$     .               +vec1(i,j ,k,2)*vec2(i,jp,k,2) )*coeff
c$$$     .         +gsuper(2)*pmag*jac)
c$$$
c$$$        t23 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,3)
c$$$     .               +vec1(i,j ,k,2)*vec2(i,jp,k,3) )*coeff
c$$$     .         +gsuper(3)*pmag*jac)
c$$$
c$$$        if (flag /= 0) then
c$$$          ijac = 1d0/jac
c$$$
c$$$          t21 = t21*ijac
c$$$          if (.not.alt_eom) t22 = t22*ijac
c$$$          t23 = t23*ijac
c$$$        endif
c$$$  ===============================================        

c$$$  ===============================================        
c$$$        coeff = ijacp*ijac0*jac
c$$$        t21 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,1)      
c$$$     .               +vec1(i,j ,k,2)*vec2(i,jp,k,1) )*coeff
c$$$     .         +gsuper(1)*pmag)
c$$$
c$$$        t22 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,2)      
c$$$     .               +vec1(i,j ,k,2)*vec2(i,jp,k,2) )*coeff
c$$$     .         +gsuper(2)*pmag)
c$$$
c$$$        t23 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,3)
c$$$     .               +vec1(i,j ,k,2)*vec2(i,jp,k,3) )*coeff
c$$$     .         +gsuper(3)*pmag)
c$$$  ===============================================        
c$$$        vec1h = 0.5*(vec1(i,j,k,:)*ijac0+vec1(i,jp,k,:)*ijacp)
c$$$        vec2h = 0.5*(vec2(i,j,k,:)*ijac0+vec2(i,jp,k,:)*ijacp)
c$$$
c$$$        t21 = (-vec1h(2)*vec2h(1)*jac+gsuper(1)*pmag)
c$$$        t22 = (-vec1h(2)*vec2h(2)*jac+gsuper(2)*pmag)
c$$$        t23 = (-vec1h(2)*vec2h(3)*jac+gsuper(3)*pmag)
c$$$  ===============================================        
        t21 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,1)*ijac0      
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,1)*ijacp )
     .         +gsuper(1)*pmag)

        t22 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,2)*ijac0  
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,2)*ijacp )
     .         +gsuper(2)*pmag)

        t23 = (-0.5*( vec1(i,jp,k,2)*vec2(i,j ,k,3)*ijac0
     .               +vec1(i,j ,k,2)*vec2(i,jp,k,3)*ijacp )
     .         +gsuper(3)*pmag)

c$$$  ===============================================                

        if (flag == 0) then
          t21 = t21*jac
          if (alt_eom) t22 = t22*jac
          t23 = t23*jac
        endif

c     End program

      end subroutine EOM_jxb_y

c     EOM_jxb_z
c     #############################################################
      subroutine EOM_jxb_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp,igr
        real(8)    :: cnv0(3),cov0(3),cnvp(3),covp(3),gsuper(3)
     .               ,car0(3),carp(3)
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,pmag,coeff

c     Begin program

        igr = igx

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(g_def%gmetric%grid(igx)%jac (i,j,kp)
     .               +g_def%gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(g_def%gmetric%grid(igx)%gsup(i,j,kp,3,:)
     .               +g_def%gmetric%grid(igx)%gsup(i,j,k ,3,:))

        if (flag /= 0) then
          jacp = g_def%gmetric%grid(igx)%jac(i,j,kp)
          jac0 = g_def%gmetric%grid(igx)%jac(i,j,k )
        else
          jacp = jac
          jac0 = jac
        endif

        ijac0 = 1d0/jac0
        ijacp = 1d0/jacp

        !Magnetic pressure
c$$$  ===============================================        
c$$$        car0 = XformToCar(g_def,i,j,k ,igx,vec1(i,j,k ,:),.false.)
c$$$        carp = XformToCar(g_def,i,j,kp,igx,vec2(i,j,kp,:),.false.)
c$$$        pmag = (carp(1)*car0(1)
c$$$     .         +car0(1)*carp(1)
c$$$     .         +carp(2)*car0(2)
c$$$     .         +car0(2)*carp(2)
c$$$     .         +carp(3)*car0(3)
c$$$     .         +car0(3)*carp(3))*0.25
c$$$  ===============================================        
        cnv0 = vec1(i,j,k ,:)
        cnvp = vec2(i,j,kp,:)

        cov0 = XformToCov(g_def,i,j,k ,igx,cnv0)
        covp = XformToCov(g_def,i,j,kp,igx,cnvp)

        cnv0 = cnv0*ijac0
        cnvp = cnvp*ijacp

        pmag = (cnvp(1)*cov0(1)
     .         +cnv0(1)*covp(1)
     .         +cnvp(2)*cov0(2)
     .         +cnv0(2)*covp(2)
     .         +cnvp(3)*cov0(3)
     .         +cnv0(3)*covp(3))*0.25

        !Assemble tensor
c$$$  ===============================================                
c$$$        coeff = ijacp*ijac0*jac**2
c$$$
c$$$        t31 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,1)
c$$$     .               +vec1(i,j,k ,3)*vec2(i,j,kp,1) )*coeff
c$$$     .         +gsuper(1)*pmag*jac)
c$$$
c$$$        t32 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,2)      
c$$$     .               +vec1(i,j,k ,3)*vec2(i,j,kp,2) )*coeff
c$$$     .         +gsuper(2)*pmag*jac)
c$$$
c$$$        t33 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,3)
c$$$     .               +vec1(i,j,k ,3)*vec2(i,j,kp,3) )*coeff
c$$$     .         +gsuper(3)*pmag*jac)
c$$$
c$$$        if (flag /= 0) then
c$$$          ijac = 1d0/jac
c$$$
c$$$          t31 = t31*ijac
c$$$          if (.not.alt_eom) t32 = t32*ijac
c$$$          t33 = t33*ijac
c$$$        endif
c$$$  ===============================================

c$$$  ===============================================                
c$$$        coeff = ijacp*ijac0*jac
c$$$        t31 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,1)
c$$$     .               +vec1(i,j,k ,3)*vec2(i,j,kp,1) )*coeff
c$$$     .         +gsuper(1)*pmag)
c$$$
c$$$        t32 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,2)      
c$$$     .               +vec1(i,j,k ,3)*vec2(i,j,kp,2) )*coeff
c$$$     .         +gsuper(2)*pmag)
c$$$
c$$$        t33 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,3)
c$$$     .               +vec1(i,j,k ,3)*vec2(i,j,kp,3) )*coeff
c$$$     .         +gsuper(3)*pmag)
c$$$  ===============================================                
        t31 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,1)*ijac0
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,1)*ijacp )
     .         +gsuper(1)*pmag)                          
                                                         
        t32 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,2)*ijac0  
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,2)*ijacp )
     .         +gsuper(2)*pmag)                          
                                                         
        t33 = (-0.5*( vec1(i,j,kp,3)*vec2(i,j,k ,3)*ijac0
     .               +vec1(i,j,k ,3)*vec2(i,j,kp,3)*ijacp )
     .         +gsuper(3)*pmag)
c$$$  ===============================================                

        if (flag == 0) then
          t31 = t31*jac
          if (alt_eom) t32 = t32*jac
          t33 = t33*jac
        endif

c     End program

      end subroutine EOM_jxb_z

c     EOM_gp
c     ###################################################################
      function EOM_gp(g_def,i,j,k,nx,ny,nz,igrid,pp,nc_eom_gp,vol)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of the EOM force term
c     grad(prs) at cell centers.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: igrid,i,j,k,nx,ny,nz
        real(8)    :: cnv(3)
        real(8),target :: pp(0:nx+1,0:ny+1,0:nz+1)
        logical    :: nc_eom_gp
        logical,optional,intent(IN) :: vol

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igx,igy,igz
        real(8)    :: cov(3),jac,cnv1(3),dxx,dyy,dzz
        logical    :: vol_wgt

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        vol_wgt = .false.
        if (PRESENT(vol)) vol_wgt = vol

        if (nc_eom_gp) then
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1

          call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

          dxx = 2*g_def%dxh(ig)
          dyy = 2*g_def%dyh(jg)
          dzz = 2*g_def%dzh(kg)

          cov(1) = (pp(ip,j,k)-pp(im,j,k))/dxx
          cov(2) = (pp(i,jp,k)-pp(i,jm,k))/dyy
          cov(3) = (pp(i,j,kp)-pp(i,j,km))/dzz

          cnv = XformToCnv(g_def,i,j,k,igx,cov)

          if (vol_wgt) cnv = cnv*g_def%gmetric%grid(igrid)%dvol(i,j,k)
        else
          sc2 => pp
          cnv = div_tensor(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt__eom()
     .                    ,eom_grad_p_x,eom_grad_p_y,eom_grad_p_z
     .                    ,vol=vol_wgt)
          nullify(sc2)
        endif

c     End program

      end function EOM_gp
      
c     EOM_grad_p_x
c     #############################################################
      subroutine EOM_grad_p_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                       ,t11,t12,t13,flag)
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
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip,igr
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,ptot,gsuper(3)

c     Begin program

        igr = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(g_def%gmetric%grid(igx)%jac (ip,j,k)
     .               +g_def%gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(g_def%gmetric%grid(igx)%gsup(ip,j,k,1,:)
     .               +g_def%gmetric%grid(igx)%gsup(i ,j,k,1,:))

        if (isSP(g_def,ip,j,k,igx,igy,igz).and.flag /= 0) jac = SP_flsv

        ptot = 0.5*jac*(sc2(i,j,k)+sc2(ip,j,k))

        t11 = gsuper(1)*ptot
        t12 = gsuper(2)*ptot
        t13 = gsuper(3)*ptot

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
      subroutine EOM_grad_p_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp,igr
        real(8)    :: x,y,z
        real(8)    :: jac,ijac,ptot,gsuper(3)

c     Begin program

        igr = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(g_def%gmetric%grid(igx)%jac (i,jp,k)
     .               +g_def%gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(g_def%gmetric%grid(igx)%gsup(i,jp,k,2,:)
     .               +g_def%gmetric%grid(igx)%gsup(i,j ,k,2,:))

        ptot = 0.5*jac*(sc2(i,j,k)+sc2(i,jp,k))

        t21 = gsuper(1)*ptot
        t22 = gsuper(2)*ptot
        t23 = gsuper(3)*ptot

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
      subroutine EOM_grad_p_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                       ,t31,t32,t33,flag)
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
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp,igr
        real(8)    :: x,y,z
        real(8)    :: jac,ijac,ptot,gsuper(3)

c     Begin program

        igr = igx

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(g_def%gmetric%grid(igx)%jac (i,j,kp)
     .               +g_def%gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(g_def%gmetric%grid(igx)%gsup(i,j,kp,3,:)
     .               +g_def%gmetric%grid(igx)%gsup(i,j,k ,3,:))

        ptot = 0.5*jac*(sc2(i,j,k)+sc2(i,j,kp))

        t31 = gsuper(1)*ptot
        t32 = gsuper(2)*ptot
        t33 = gsuper(3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine EOM_grad_p_z

c     EOM_gnT
c     ###################################################################
      function EOM_gnT(g_def,i,j,k,nx,ny,nz,igrid,rr,tt,nc_eom_gp,vol)
     .         result(cnv)

c     -------------------------------------------------------------------
c     Finds contravariant components (a1,a2,a3) of the EOM force term
c     grad(prs) at cell centers.
c     -------------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: igrid,i,j,k,nx,ny,nz
        real(8)    :: cnv(3)
        real(8),target :: rr(0:nx+1,0:ny+1,0:nz+1)
     $                   ,tt(0:nx+1,0:ny+1,0:nz+1)
        logical    :: nc_eom_gp
        logical,optional,intent(IN) :: vol

c     Local variables

        integer    :: ig,jg,kg,ip,im,jp,jm,kp,km,igx,igy,igz
        real(8)    :: cov(3),jac,cnv1(3),dxx,dyy,dzz
        logical    :: vol_wgt

c     Begin program

        igx = igrid
        igy = igrid
        igz = igrid

        vol_wgt = .false.
        if (PRESENT(vol)) vol_wgt = vol

        if (nc_eom_gp) then
          ip = i+1
          im = i-1
          jp = j+1
          jm = j-1
          kp = k+1
          km = k-1

          call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

          dxx = 2*g_def%dxh(ig)
          dyy = 2*g_def%dyh(jg)
          dzz = 2*g_def%dzh(kg)

          !grad(nT) with ZIP
          cov(1) = (rr(i,j,k)*(tt(ip,j,k)-tt(im,j,k))
     .             +tt(i,j,k)*(rr(ip,j,k)-rr(im,j,k)))/dxx

          cov(2) = (rr(i,j,k)*(tt(i,jp,k)-tt(i,jm,k))
     .             +tt(i,j,k)*(rr(i,jp,k)-rr(i,jm,k)))/dyy
          
          cov(3) = (rr(i,j,k)*(tt(i,j,kp)-tt(i,j,km))
     .             +tt(i,j,k)*(rr(i,j,kp)-rr(i,j,km)))/dzz

          cnv = XformToCnv(g_def,i,j,k,igx,cov)

          if (vol_wgt) cnv = cnv*g_def%gmetric%grid(igrid)%dvol(i,j,k)

        else
          sc1 => rr
          sc2 => tt
          cnv = div_tensor(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt__eom()
     .                    ,eom_grad_nT_x,eom_grad_nT_y,eom_grad_nT_z
     .                    ,vol=vol_wgt)
          nullify(sc1,sc2)
        endif

c     End program

      end function EOM_gnT

c     EOM_grad_nT_x
c     #############################################################
      subroutine EOM_grad_nT_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                       ,t11,t12,t13,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t11-t13 for grad(nT) term in EOM.
c     In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t11,t12,t13: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip,igr
        real(8)    :: x,y,z
        real(8)    :: jac,jac0,jacp,ijac,ijac0,ijacp,ptot,gsuper(3)

c     Begin program

        igr = igx

        ip = i+1
        if (flag == 0) ip = i

        jac    = 0.5*(g_def%gmetric%grid(igx)%jac (ip,j,k)
     .               +g_def%gmetric%grid(igx)%jac (i ,j,k))
        gsuper = 0.5*(g_def%gmetric%grid(igx)%gsup(ip,j,k,1,:)
     .               +g_def%gmetric%grid(igx)%gsup(i ,j,k,1,:))

        if (isSP(g_def,ip,j,k,igx,igy,igz).and.flag /= 0) jac = SP_flsv

        ptot = 0.5*jac*(sc1(ip,j,k)*sc2(i ,j,k)
     .                 +sc1(i ,j,k)*sc2(ip,j,k))

        t11 = gsuper(1)*ptot
        t12 = gsuper(2)*ptot
        t13 = gsuper(3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine EOM_grad_nT_x

c     EOM_grad_nT_y
c     #############################################################
      subroutine EOM_grad_nT_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                    ,t21,t22,t23,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t21-t23 for grad(nT) term in EOM.
c     In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t21,t22,t23: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp,igr
        real(8)    :: x,y,z
        real(8)    :: jac,ijac,ptot,gsuper(3)

c     Begin program

        igr = igx

        jp = j+1
        if (flag == 0) jp = j

        jac    = 0.5*(g_def%gmetric%grid(igx)%jac (i,jp,k)
     .               +g_def%gmetric%grid(igx)%jac (i,j ,k))
        gsuper = 0.5*(g_def%gmetric%grid(igx)%gsup(i,jp,k,2,:)
     .               +g_def%gmetric%grid(igx)%gsup(i,j ,k,2,:))

        ptot = 0.5*jac*(sc1(i,jp,k)*sc2(i,j ,k)
     .                 +sc1(i,j ,k)*sc2(i,jp,k))

        t21 = gsuper(1)*ptot
        t22 = gsuper(2)*ptot
        t23 = gsuper(3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t21 = t21*ijac
          if (.not.alt_eom) t22 = t22*ijac
          t23 = t23*ijac
        endif

c     End program

      end subroutine EOM_grad_nT_y

c     EOM_grad_nT_z
c     #############################################################
      subroutine EOM_grad_nT_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
     .                       ,t31,t32,t33,flag)
c     -------------------------------------------------------------
c     Calculates tensor components t31-t33 for grad(nT) term in EOM.
c     In the call sequence:
c       * i,j,k: grid position
c       * nx,ny,nz: grid size
c       * igx,igy,igz: grid level (for MG evaluations)
c       * alt_eom: whether to use alternate EOM in singular coord.
c                  systems or not.
c       * t31,t32,t33: tensor components.
c       * flag: whether evaluation is a cell center i,j,k (flag=0)
c               or at cell face i+1/2,j,k (flag /= 0)
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp,igr
        real(8)    :: x,y,z
        real(8)    :: jac,ijac,ptot,gsuper(3)

c     Begin program

        igr = igx

        kp = k+1
        if (flag == 0) kp = k

        jac    = 0.5*(g_def%gmetric%grid(igx)%jac (i,j,kp)
     .               +g_def%gmetric%grid(igx)%jac (i,j,k ))
        gsuper = 0.5*(g_def%gmetric%grid(igx)%gsup(i,j,kp,3,:)
     .               +g_def%gmetric%grid(igx)%gsup(i,j,k ,3,:))

        ptot = 0.5*jac*(sc1(i,j,kp)*sc2(i,j,k )
     .                 +sc1(i,j,k )*sc2(i,j,kp))

        t31 = gsuper(1)*ptot
        t32 = gsuper(2)*ptot
        t33 = gsuper(3)*ptot

        if (flag /= 0) then
          ijac = 1d0/jac

          t31 = t31*ijac
          if (.not.alt_eom) t32 = t32*ijac
          t33 = t33*ijac
        endif

c     End program

      end subroutine EOM_grad_nT_z

c     EOM_diff_x
c     #############################################################
      subroutine EOM_diff_x(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t11,t12,t13
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,ip,im,ipp,igpp,igp,igm,igr
        real(8)    :: x,y,z,vxp,qip,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,ijac,vish,visp!,vis0,jac0,jacp,ijac0,ijacp

c     Begin program

        if (     .not.associated(coef)
     .      .or. .not.associated(vec1)) then
          call pstop('EOM_diff_x','Pointers not associated')
        endif

        igr = igx

        ip = i+1
        if (flag == 0) ip = i

        !Harmonic average for calculation of viscosity coeff. at faces
        if (coef(ip,j,k)*coef(i,j,k) > 0d0) then
          vish = 2.*coef(ip,j,k)*coef(i,j,k)/(coef(ip,j,k)+coef(i,j,k))
        else
          t11 = 0d0 ; t12 = 0d0; t13 = 0d0
          return
        endif

        jac = 0.5*(g_def%gmetric%grid(igx)%jac (ip,j,k)
     .            +g_def%gmetric%grid(igx)%jac (i ,j,k))

cc        if (isSP(ip,j,k,igx,igy,igz).and.flag /= 0) jac = SP_flsv
        if (isSP2(g_def,i+1,igx).and.flag /= 0) then
          ijac = 0d0
        else
          ijac = 1d0/jac
        endif

c     Diffusive part

        if (symm_pi) then
          if (flag /= 0) then
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igx,vec1,1,1,3,cnv=.true.)
          else
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igx,vec1,0,1,3,cnv=.true.)
          endif
          nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
        else
          if (flag /= 0) then
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igx,vec1,1,1,1,cnv=.true.)
          else
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igx,vec1,0,1,1,cnv=.true.)
          endif
        endif

        t11 = -vish*nabla_v_cnv(1,1)
        t12 = -vish*nabla_v_cnv(1,2)
        t13 = -vish*nabla_v_cnv(1,3)

c       Postprocessing

        if (flag /= 0) then
          t11 = t11*ijac
          if (.not.alt_eom) t12 = t12*ijac
          t13 = t13*ijac
        endif

c     End program

      end subroutine EOM_diff_x

c     EOM_diff_y
c     #############################################################
      subroutine EOM_diff_y(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t21,t22,t23
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,jp,jm,jpp,jgpp,jgp,jgm,igr
        real(8)    :: x,y,z,vyp,qjp,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,ijac,vish,visp

c     Begin program

        if (     .not.associated(coef)
     .      .or. .not.associated(vec1)) then
          call pstop('EOM_diff_y','Pointers not associated')
        endif

        igr = igx

        jp = j+1
        if (flag == 0) jp = j

        !Harmonic average for calculation of viscosity coeff. at faces
        if (coef(i,jp,k)*coef(i,j,k) > 0d0) then
          vish = 2.*coef(i,jp,k)*coef(i,j,k)/(coef(i,jp,k)+coef(i,j,k))
        else
          t21 = 0d0 ; t22 = 0d0; t23 = 0d0
          return
        endif

        jac = 0.5*(g_def%gmetric%grid(igx)%jac (i,jp,k)
     .            +g_def%gmetric%grid(igx)%jac (i,j ,k))

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

c       Diffusive part

        if (symm_pi) then
          if (flag /= 0) then
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igx,vec1,2,1,3,cnv=.true.)
          else
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igx,vec1,0,1,3,cnv=.true.)
          endif
          nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
        else
          if (flag /= 0) then
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igx,vec1,2,2,2,cnv=.true.)
          else
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igx,vec1,0,2,2,cnv=.true.)
          endif
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
      subroutine EOM_diff_z(g_def,i,j,k,nx,ny,nz,igx,igy,igz,alt_eom
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
c     -------------------------------------------------------------

        implicit none

c     Call variables

        type(grid_mg_def),pointer :: g_def

        integer    :: i,j,k,flag,nx,ny,nz,igx,igy,igz
        real(8)    :: t31,t32,t33
        logical    :: alt_eom

c     Local variables

        integer    :: ig,jg,kg,kp,km,kpp,kgpp,kgp,kgm,igr
        real(8)    :: x,y,z,vzp,qkp,ph,nabla_v_cnv(3,3)
        real(8)    :: jac,ijac,vish,visp

c     Begin program

        if (     .not.associated(coef)
     .      .or. .not.associated(vec1)) then
          call pstop('EOM_diff_z','Pointers not associated')
        endif

        igr = igx

        kp = k+1
        if (flag == 0) kp = k

        !Harmonic average for calculation of viscosity coeff. at faces
        if (coef(i,j,kp)*coef(i,j,k) > 0d0) then
          vish = 2.*coef(i,j,kp)*coef(i,j,k)/(coef(i,j,kp)+coef(i,j,k))
        else
          t31 = 0d0 ; t32 = 0d0; t33 = 0d0
          return
        endif

        jac = 0.5*(g_def%gmetric%grid(igr)%jac (i,j,kp)
     .            +g_def%gmetric%grid(igr)%jac (i,j,k ))

        call getMGmap(g_def,i,j,k,igx,igy,igz,ig,jg,kg)

c       Diffusive part

        if (symm_pi) then
          if (flag /= 0) then
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igr,vec1,3,1,3,cnv=.true.)
          else
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igr,vec1,0,1,3,cnv=.true.)
          endif
          nabla_v_cnv = nabla_v_cnv + transpose(nabla_v_cnv)
        else
          if (flag /= 0) then
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igr,vec1,3,3,3,cnv=.true.)
          else
            nabla_v_cnv =fnabla_v(g_def,i,j,k,igr,vec1,0,3,3,cnv=.true.)
          endif
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
      function EOM_divPi_ijk(g_def,i,j,k,nxx,nyy,nzz,igr,v0,dff,symm)
     .         result(cnv)

c     ---------------------------------------------------------------------
c     Find divergence of ion pressure stress tensor. Return contravariant
c     components.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nxx,nyy,nzz,igr
      real(8) :: cnv(3)
      real(8),target :: v0 (0:nxx+1,0:nyy+1,0:nzz+1,3)
      real(8),target :: dff(0:nxx+1,0:nyy+1,0:nzz+1)

      logical :: symm

c     Local variables

c     Begin program

      vec1 => v0
      coef => dff
      symm_pi = symm
      
      cnv = div_tensor(g_def,i,j,k,nxx,nyy,nzz,igr,igr,igr,alt__eom()
     .                ,eom_diff_x,eom_diff_y,eom_diff_z
     .                ,vol=.false.)

      nullify(vec1,coef)

c     End program

      end function EOM_divPi_ijk

c     EOM_divPi
c     #####################################################################
      function EOM_divPi(g_def,igr,v0,diff,symm) result(cnv)

c     ---------------------------------------------------------------------
c     Find divergence of electron pressure stress tensor. Returns
c     contravariant components.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr
      real(8) :: v0  (0:,0:,0:,:)
      real(8) :: diff(0:,0:,0:)

      real(8) :: cnv (0:size(v0,1)-1
     $               ,0:size(v0,2)-1
     $               ,0:size(v0,3)-1
     $               ,  size(v0,4)  )

      logical :: symm

c     Local variables

      integer :: i,j,k,nxx,nyy,nzz

c     Begin program

      cnv = 0d0

      nxx = g_def%nxv(igr)
      nyy = g_def%nyv(igr)
      nzz = g_def%nzv(igr)

      do k=1,nzz
        do j=1,nyy
          do i=1,nxx
            cnv(i,j,k,:)
     .        = EOM_divPi_ijk(g_def,i,j,k,nxx,nyy,nzz,igr,v0,diff,symm)
          enddo
        enddo
      enddo

c     End program

      end function EOM_divPi

c$$$c     EOM_divPe_ijk
c$$$c     #####################################################################
c$$$      function EOM_divPe_ijk(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,diff)
c$$$     .         result(cnv)
c$$$
c$$$c     ---------------------------------------------------------------------
c$$$c     Find divergence of electron pressure stress tensor. Returns
c$$$c     contravariant components.
c$$$c     ---------------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer :: i,j,k,nxx,nyy,nzz,igx,igy,igz
c$$$      real(8) :: cnv(3)
c$$$      real(8) :: v0  (0:nxx+1,0:nyy+1,0:nzz+1,3)
c$$$      real(8) :: diff(0:nxx+1,0:nyy+1,0:nzz+1)
c$$$
c$$$c     Local variables
c$$$
c$$$c     Begin program
c$$$
c$$$      cnv = EOM_divPi_ijk(i,j,k,nxx,nyy,nzz,igx,igy,igz,v0,diff)
c$$$
c$$$c     End program
c$$$
c$$$      end function EOM_divPe_ijk
c$$$
c$$$c     EOM_divPe
c$$$c     #####################################################################
c$$$      function EOM_divPe(nxx,nyy,nzz,igx,igy,igz,v0,field) result(cnv)
c$$$
c$$$c     ---------------------------------------------------------------------
c$$$c     Find divergence of electron pressure stress tensor. Returns
c$$$c     contravariant components.
c$$$c     ---------------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      integer :: nxx,nyy,nzz,igx,igy,igz
c$$$      real(8) :: cnv  (0:nxx+1,0:nyy+1,0:nzz+1,3)
c$$$     .          ,v0   (0:nxx+1,0:nyy+1,0:nzz+1,3)
c$$$     .          ,field(0:nxx+1,0:nyy+1,0:nzz+1)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: i,j,k
c$$$
c$$$c     Begin program
c$$$
c$$$      cnv = 0d0
c$$$
c$$$      do k=1,nzz
c$$$        do j=1,nyy
c$$$          do i=1,nxx
c$$$            cnv(i,j,k,:) = EOM_divPe_ijk(i,j,k,nxx,nyy,nzz,igx,igy,igz
c$$$     .                                  ,v0,field)
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$c     End program
c$$$
c$$$      end function EOM_divPe

c     EOM_Pi_ijk
c     #####################################################################
      function EOM_Pi_ijk(g_def,i,j,k,nx,ny,nz,igr,v0,diff,symm)
     .         result(pi_tnsr)

c     ---------------------------------------------------------------------
c     Find pressure stress tensor, in contravariant representation.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,nx,ny,nz,igr
      real(8) :: pi_tnsr(3,3)
      real(8),target :: v0  (0:nx+1,0:ny+1,0:nz+1,3)
      real(8),target :: diff(0:nx+1,0:ny+1,0:nz+1)

      logical :: symm

c     Local variables

c     Begin program

      vec1 => v0
      coef => diff
      symm_pi = symm

      call eom_diff_x(g_def,i,j,k,nx,ny,nz,igr,igr,igr,.false.
     .               ,pi_tnsr(1,1),pi_tnsr(1,2),pi_tnsr(1,3),0)
      call eom_diff_y(g_def,i,j,k,nx,ny,nz,igr,igr,igr,.false.
     .               ,pi_tnsr(2,1),pi_tnsr(2,2),pi_tnsr(2,3),0)
      call eom_diff_z(g_def,i,j,k,nx,ny,nz,igr,igr,igr,.false.
     .               ,pi_tnsr(3,1),pi_tnsr(3,2),pi_tnsr(3,3),0)
      nullify(vec1,coef)

c     End program

      end function EOM_Pi_ijk

c     E_fld_res_B
c     #####################################################################
      function E_fld_res_B(g_def,igr,vcnv,bcnv,jcov,eta2) result(E_res)

c     ---------------------------------------------------------------------
c     Computes ideal+resistive electric field, E=-(vxb)+eta*j.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr
      real(8) :: jcov  (0:,0:,0:,:)
     .          ,vcnv  (0:,0:,0:,:)
     .          ,bcnv  (0:,0:,0:,:)
     .          ,eta2  (0:,0:,0:)

      real(8) :: E_res(0:size(vcnv,1)-1
     .                ,0:size(vcnv,2)-1
     .                ,0:size(vcnv,3)-1
     .                ,  size(vcnv,4)  )

c     Local variables

      integer :: nx,ny,nz,i,j,k,ig,jg,kg,ieq,bcnd(6,3)
      real(8) :: eps,betap,qa,jbs,cov(3),cnv(3),bnorm,x1,y1,z1

c     Begin program

      nx = size(vcnv,1)-2
      ny = size(vcnv,2)-2
      nz = size(vcnv,3)-2

c     Resistive term

      if (eta > 0d0) then

cc        pe = rho*tmp/a_p
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
cc          !Impose topological BCs
cc          bcnd(:,1) = bcond
cc          bcnd(:,2) = bcond
cc          bcnd(:,3) = bcond
cc
cc          call setBC(IAX,3,nx,ny,nz,vdummy,E_res,vzeros,bcnd
cc     .              ,igx,igy,igz,iorder=2,is_cnv=.false.)
cc        endif

        !Ohmic current
        E_res(:,:,:,1) = eta2*jcov(:,:,:,1)
        E_res(:,:,:,2) = eta2*jcov(:,:,:,2)
        E_res(:,:,:,3) = eta2*jcov(:,:,:,3)

      else

        E_res = 0d0

      endif

c     VxB term

      E_res = E_res - crossProduct(g_def,igr,vcnv,bcnv,.true.)
      
c     External field

      E_res(:,:,:,1) = E_res(:,:,:,1) - E0(1)
      E_res(:,:,:,2) = E_res(:,:,:,2) - E0(2)
      E_res(:,:,:,3) = E_res(:,:,:,3) - E0(3)

c     End program

      end function E_fld_res_B

c     E_fld_res_A
c     #####################################################################
      function E_fld_res_A(g_def,igr,vcnv,acov,jcov,eta2) result(E_res)

c     ---------------------------------------------------------------------
c     Computes ideal+resistive electric field, E=-(v x curl(A))+eta*j.
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr
      real(8) :: jcov(0:,0:,0:,:)
     .          ,eta2(0:,0:,0:)
     .          ,vcnv(0:,0:,0:,:)
     .          ,acov(0:,0:,0:,:)

      real(8) :: E_res(0:size(vcnv,1)-1
     .                ,0:size(vcnv,2)-1
     .                ,0:size(vcnv,3)-1
     .                ,  size(vcnv,4)  )

c     Local variables

      integer :: nx,ny,nz,i,j,k,ig,jg,kg,ieq,bcnd(6,3)
      real(8) :: eps,betap,qa,jbs,cov(3),cnv(3),bnorm,x1,y1,z1

c     Begin program

      nx = size(vcnv,1)-2
      ny = size(vcnv,2)-2
      nz = size(vcnv,3)-2

c     Resistive term

      if (eta > 0d0) then

        !Ohmic current
        if (vlap_etaj) then
c$$$          E_res = -veclap(g_def,igr,acnv,alt_eom=.false.)
c$$$          E_res = XformToCov(g_def,igr,E_res)  !Worse than veclap_cov
c$$$          ============
c$$$          E_res = -veclap_cnv(g_def,igr,acnv)
c$$$          E_res = XformToCov(g_def,igr,E_res)  !Does not finish simulation
c$$$          ============
          E_res =-veclap_cov(g_def,igr,acov) !Finishes sim, but strong issues at origin in rho, p
c$$$          ============
c$$$          E_res = curlcurl_cov(g_def,igr,acov)  !Finishes sim, but strong issues at origin in rho, p
c$$$          ============
c$$$          E_res = curlcurl(g_def,igr,acnv,acov,.false.)
c$$$          ============
c$$$          E_res(:,:,:,1) = jcov(:,:,:,1)          !Fixes all issues at origin
        else
          E_res = jcov
        endif

        E_res(:,:,:,1) = eta2*E_res(:,:,:,1)
        E_res(:,:,:,2) = eta2*E_res(:,:,:,2)
        E_res(:,:,:,3) = eta2*E_res(:,:,:,3)

      else

        E_res = 0d0

      endif

c     VxB term

      !Topological BCs
      bcnd(:,1) = bcond
      bcnd(:,2) = bcond
      bcnd(:,3) = bcond

c$$$      if (bcSP()) then
        E_res = E_res 
     .  + tensor_nc_advc_mesh (g_def,igr,vcnv,acov,bcnd,ZIP,.true.)
     .  - tensor_nc_advc_mesh2(g_def,igr,vcnv,acov     ,ZIP,.true.
     .                        ,transpose=.true.)
c$$$      else
c$$$        E_res = E_res 
c$$$     .  + tensor_nc_advc_mesh (g_def,igr,vcnv,acov,bcnd,v_advect,.true.)
c$$$     .  - tensor_nc_advc_mesh2(g_def,igr,vcnv,acov     ,ZIP     ,.true.
c$$$     .                        ,transpose=.true.)
c$$$      endif

c     External field

      E_res(:,:,:,1) = E_res(:,:,:,1) - E0(1)
      E_res(:,:,:,2) = E_res(:,:,:,2) - E0(2)
      E_res(:,:,:,3) = E_res(:,:,:,3) - E0(3)

c     End program

      end function E_fld_res_A

c$$$c     E_fld_Hall
c$$$c     #####################################################################
c$$$      function E_fld_Hall(g_def,igr,nx,ny,nz) result(E_h)
c$$$
c$$$c     ---------------------------------------------------------------------
c$$$c     Computes Hall electric field, in two forms: ion EOM, and e EOM
c$$$c     ---------------------------------------------------------------------
c$$$
c$$$      implicit none
c$$$
c$$$c     Call variables
c$$$
c$$$      type(grid_mg_def),pointer :: g_def
c$$$
c$$$      integer :: nx,ny,nz,igr
c$$$      real(8) :: E_h(0:nx+1,0:ny+1,0:nz+1,3)
c$$$
c$$$c     Local variables
c$$$
c$$$      integer :: i,j,k,ig,jg,kg,ieq,bcnd(6,3),igx,igy,igz
c$$$      real(8) :: pe(0:nx+1,0:ny+1,0:nz+1),cov(3),cnv(3)
c$$$     .          ,bcnv_star(0:nx+1,0:ny+1,0:nz+1,3)
c$$$
c$$$      real(8),pointer,dimension(:,:,:) :: jac
c$$$
c$$$c     Begin program
c$$$
c$$$      igx = igr
c$$$      igy = igr
c$$$      igz = igr
c$$$
c$$$      pe = rho*tmpe            !Electron pressure
c$$$
c$$$c     Standard Hall:  E_hall = (j x B* - grad(pe) - div(Pi_e))/rho - de^2/di^2*d(v_e)/dt
c$$$
c$$$      if (electron_hall) then
c$$$
c$$$        if (de > 0d0) then
c$$$          bcnv_star = bcnv -de**2/di*curl(g_def,igx,vecov)
c$$$        else
c$$$          bcnv_star = bcnv
c$$$        endif
c$$$
c$$$        E_h = crossProduct(g_def,igx,jcnv,bcnv_star,.true.)
c$$$     .       -grad(g_def,igx,pe)
c$$$     .       -XformVector(g_def,igx,div_pe,.false.)
c$$$
c$$$        E_h(:,:,:,1) = E_h(:,:,:,1)/rho
c$$$        E_h(:,:,:,2) = E_h(:,:,:,2)/rho
c$$$        E_h(:,:,:,3) = E_h(:,:,:,3)/rho
c$$$
c$$$        !Electron inertia term
c$$$        if (de > 0d0) then
c$$$          E_h = E_h - de**2/di**2*cnp*(vecov-vecov_n)/alpha/dt  !Works for BDF2
c$$$        endif
c$$$
c$$$c     Alternate Hall: E_hall = (dv/dt + v.grad(v) + div(Pi)/rho + grad(p_i)/rho)
c$$$      else
c$$$
c$$$        !Temporal term: dv/dt
c$$$        if (save_tn.or.source_eval) then
c$$$          E_h = 0d0
c$$$        else
c$$$          E_h = cnp*(vcnv-vcnv_n)/alpha/dt  !Change constant in U_Av_ijk; works for BDF
c$$$        endif
c$$$
c$$$        !Stress
c$$$        E_h(:,:,:,1) = E_h(:,:,:,1) + div_pi(:,:,:,1)/rho
c$$$        E_h(:,:,:,2) = E_h(:,:,:,2) + div_pi(:,:,:,2)/rho
c$$$        E_h(:,:,:,3) = E_h(:,:,:,3) + div_pi(:,:,:,3)/rho
c$$$
c$$$        !Advection (covariant)
c$$$        E_h = XformVector(g_def,igx,E_h,.false.)
c$$$     .       +tensor_nc_advc_mesh(g_def,igx,vcnv,vcnv
c$$$     .                           ,gv%aux%vec_list(IVCNV)%bconds
c$$$     .                           ,v_advect,.true.)
c$$$c$$$        do k=1,nz
c$$$c$$$          do j=1,ny
c$$$c$$$            do i=1,nx
c$$$c$$$cc              cnv = tensor_nc_advc(i,j,k,nx,ny,nz,igx,igy,igz,vcnv,vcnv
c$$$c$$$cccc     .                            ,min(advect,2))
c$$$c$$$cc     .                            ,1)
c$$$c$$$              cnv = tensor_advc(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$c$$$     .                         ,vcnv,vcnv,v_advect)
c$$$c$$$     .             -vcnv(i,j,k,:)*div(g_def,i,j,k,nx,ny,nz,igx,igy,igz
c$$$c$$$     .                               ,vcnv)
c$$$c$$$     .             +E_h(i,j,k,:)
c$$$c$$$
c$$$c$$$              E_h(i,j,k,:)
c$$$c$$$     .             = matmul(g_def%gmetric%grid(igx)%gsub(i,j,k,:,:),cnv)
c$$$c$$$            enddo
c$$$c$$$          enddo
c$$$c$$$        enddo
c$$$
c$$$        !Ion pressure term: grad(pi)/rho
c$$$        if (temp_ratio > 0d0) then
c$$$cc          Eh = Eh + temp_ratio*grad(g_def,igx,pe)/rho
c$$$          do k=1,nz
c$$$            do j=1,ny
c$$$              do i=1,nx
c$$$                cov = temp_ratio
c$$$     .               *grad(g_def,i,j,k,nx,ny,nz,igx,igy,igz,pe)
c$$$
c$$$                E_h(i,j,k,:) = E_h(i,j,k,:) + cov/rho(i,j,k)
c$$$              enddo
c$$$            enddo
c$$$          enddo
c$$$
c$$$        endif
c$$$
c$$$      endif
c$$$
c$$$c     End program
c$$$
c$$$      end function E_fld_Hall

c     E_fld_eHall
c     #####################################################################
      function E_fld_eHall(g_def,igr,rho,tmpe,bcnv,div_pe,vecov,vecov_n)
     .         result(E_h)

c     ---------------------------------------------------------------------
c     Computes electron Hall electric field:
c       E_hall = (j x B* - grad(pe) - div(Pi_e))/rho - de^2/di^2*d(v_e)/dt
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr
      real(8) :: rho(0:,0:,0:),tmpe(0:,0:,0:),div_pe(0:,0:,0:,:)
     .          ,bcnv(0:,0:,0:,:),vecov(0:,0:,0:,:),vecov_n(0:,0:,0:,:)
      real(8) :: E_h(0:size(rho,1)-1
     .              ,0:size(rho,2)-1
     .              ,0:size(rho,3)-1,3)

c     Local variables

      integer :: i,j,k,ig,jg,kg,ieq,bcnd(6,3)
      real(8) :: cov(3),cnv(3)
      real(8) :: pe(0:size(rho,1)-1
     .             ,0:size(rho,2)-1
     .             ,0:size(rho,3)-1)
     .          ,bcnv_star(0:size(rho,1)-1
     .                    ,0:size(rho,2)-1
     .                    ,0:size(rho,3)-1,3)

      real(8),pointer,dimension(:,:,:) :: jac

c     Begin program

      pe = rho*tmpe            !Electron pressure 

      if (de > 0d0) then
        bcnv_star = bcnv -de**2/di*curl(g_def,igr,vecov)
      else
        bcnv_star = bcnv
      endif

      E_h = crossProduct(g_def,igr,jcnv,bcnv_star,.true.)
     .     -grad        (g_def,igr,pe)
     .     -XformVector (g_def,igr,div_pe,.false.)

      E_h(:,:,:,1) = E_h(:,:,:,1)/rho
      E_h(:,:,:,2) = E_h(:,:,:,2)/rho
      E_h(:,:,:,3) = E_h(:,:,:,3)/rho

      !Electron inertia term
      if (de > 0d0) then
        E_h = E_h - (de/di)**2*cnp*(vecov-vecov_n)/alpha/dt  !Does NOT work for BDF2
      endif

c     End program

      end function E_fld_eHall

c     E_fld_iHall
c     #####################################################################
      function E_fld_iHall(g_def,igr,rho,tmpe,div_pi,vcnv,vcnv_n) 
     .         result(E_h)

c     ---------------------------------------------------------------------
c     Computes ion Hall electric field:
c          E_hall = (dv/dt + v.grad(v) + div(Pi)/rho + grad(p_i)/rho)
c     ---------------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr
      real(8) :: rho(0:,0:,0:),tmpe(0:,0:,0:),div_pi(0:,0:,0:,:)
     .          ,vcnv(0:,0:,0:,:),vcnv_n(0:,0:,0:,:)
      real(8) :: E_h(0:size(rho,1)-1
     .              ,0:size(rho,2)-1
     .              ,0:size(rho,3)-1,3)

c     Local variables

      real(8) :: p_i(0:size(rho,1)-1
     .              ,0:size(rho,2)-1
     .              ,0:size(rho,3)-1),idt

c     Begin program

      !Pressure terms: (grad(pi)+div(Pi))/rho (contravariant)
      p_i = temp_ratio*rho*tmpe
      E_h = grad(g_def,igr,p_i)  !Covariant
      E_h = div_pi + XformVector(g_def,igr,E_h,.true.) !Contravariant

      E_h(:,:,:,1) = E_h(:,:,:,1)/rho
      E_h(:,:,:,2) = E_h(:,:,:,2)/rho
      E_h(:,:,:,3) = E_h(:,:,:,3)/rho

      !Temporal derivative: dv/dt (contravariant)
      idt = cnp/(alpha*dt)
      E_h = E_h + idt*(vcnv-vcnv_n)

      !Advection: v.grad(v) (contravariant)
      E_h = E_h
     .     +tensor_nc_advc_mesh(g_def,igr,vcnv,vcnv
     .                         ,gv%aux%vec_list(IVCNV)%bconds
     .                         ,v_advect,.false.)

      !Transform to covariant
      E_h = XformVector(g_def,igr,E_h,.false.)

c     End program

      end function E_fld_iHall

c     res_ijk
c     #################################################################
      function res_ijk(g_def,i,j,k,igrid,field) result(rsvty)
c     -----------------------------------------------------------------
c     This function computes the local value of the resistivity on the
c     node (i,j,k) of grid igrid according to field.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igrid
      real(8) :: field,rsvty

c     Local variables

      integer :: ig,jg,kg
      real(8) :: rr,th,zz

c     Begin program

      if (spitzer) then
        rsvty = min(eta*res_spitzer(field,te0_ref),1d-1)
      else
        select case (equil)
        case ('ppnch','ppnsl','ppnst','ppn3d','p3nsl','2dhel','ohtor')
          call find_polar_coords(g_def,igrid,i,j,k,rr,th,zz,coords)
          rsvty = eta*res_rfx(rr)
        case default
          rsvty = eta
        end select
      endif

c     End program

      end function res_ijk

c     res_mesh
c     #############################################################
      function res_mesh(g_def,igrid,field) result(rsvty)
c     -------------------------------------------------------------
c     This function computes the resistivity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: field(0:,0:,0:)
      real(8) :: rsvty(0:size(field,1)-1
     .                ,0:size(field,2)-1
     .                ,0:size(field,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k
      real(8) :: fmax,fmaxg

c     Begin program

      nx = size(field,1)-2
      ny = size(field,2)-2
      nz = size(field,3)-2

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            rsvty(i,j,k)=res_ijk(g_def,i,j,k,igrid,field(i,j,k))
          enddo
        enddo
      enddo

c     End program

      end function res_mesh

c     vis_ijk
c     #############################################################
      function vis_ijk(g_def,i,j,k,igrid,tmp,rho) result(vscty)
c     -------------------------------------------------------------
c     This function computes the viscosity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: i,j,k,igrid
      real(8) :: tmp,rho,vscty

c     Local variables

      integer :: ig,jg,kg
      real(8) :: rr,th,zz

c     Begin program

      select case (equil)
      case ('ppnch','ppnsl','ppnst','ppn3d','p3nsl','2dhel','ohtor')
        call find_polar_coords(g_def,igrid,i,j,k,rr,th,zz,coords)
        vscty = nu*vis_rfx(rr)
      case default
        vscty = nu
      end select

      if (nu_is_rhonu) vscty = rho*vscty

c     End program

      end function vis_ijk

c     vis_mesh
c     #############################################################
      function vis_mesh(g_def,igrid,tmp,rho) result(vscty)
c     -------------------------------------------------------------
c     This function computes the viscosity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: tmp(0:,0:,0:),rho(0:,0:,0:)
     .          ,vscty(0:size(tmp,1)-1
     .                ,0:size(tmp,2)-1
     .                ,0:size(tmp,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = size(tmp,1)-2
      ny = size(tmp,2)-2
      nz = size(tmp,3)-2

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
           vscty(i,j,k)=vis_ijk(g_def,i,j,k,igrid,tmp(i,j,k),rho(i,j,k))
          enddo
        enddo
      enddo

c     End program

      end function vis_mesh

c     hres
c     ###############################################################
      function hres(g_def,igr,bcnv,rho) result(h_res)
c     ---------------------------------------------------------------
c     Finds hyper-resistivity coefficient (global on the grid)
c     ---------------------------------------------------------------

      use app_iosetup

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igr
      real(8) :: bcnv (0:,0:,0:,:)
     .          ,rho  (0:,0:,0:)
      real(8) :: h_res(0:size(rho,1)-1
     .                ,0:size(rho,2)-1
     .                ,0:size(rho,3)-1)

c     Local variables

      integer :: nx,ny,nz,ig,jg,kg,i,j,k

      real(8) :: idx,idy,idz,k2,kk,bk_par,lheta,vol,lvol,maxv,maxvg

c     Begin program

      nx = size(rho,1)-2
      ny = size(rho,2)-2
      nz = size(rho,3)-2

cc      hres = 0d0
cc      vol = 0d0
cc
cc      do k=1,nz
cc        do j=1,ny
cc          do i=1,nx
cc            call getMGmap(g_def,i,j,k,igr,igr,igr,ig,jg,kg)
cc
cc            !kk
cc            idx  = pi/g_def%dx(ig)
cc            if (nx == 1) idx = 0d0
cc            idy  = pi/g_def%dy(jg)
cc            if (ny == 1) idy = 0d0
cc            idz  = pi/g_def%dz(kg)
cc            if (nz == 1) idz = 0d0
cc
cc            k2=vectorNorm(i,j,k,igr,igr,igr,idx,idy,idz,.true.)
cc            kk=sqrt(k2)
cc
cc            !|B|*k_par
cc            bk_par=scalarProduct(i,j,k,igr,igr,igr,idx,idy,idz
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

      if (di == 0d0.or.heta==0d0) then

        h_res = 0d0

      elseif (heta > 0d0) then

c$$$        if (slava) then
c$$$          h_res = heta
c$$$        else
c$$$          h_res = rho*heta
c$$$        endif
        if (heta_is_rhoheta) then
          h_res = rho*heta
        else
          h_res = heta
        endif

      else

        do k=0,nz+1
          do j=0,ny+1
            do i=0,nx+1
              call getMGmap(g_def,min(max(i,1),nx)
     .                     ,min(max(j,1),ny)
     .                     ,min(max(k,1),nz),igr,igr,igr,ig,jg,kg)

              !kk
              idx  = 2d0/g_def%dx(ig)
              if (nx == 1) idx = 0d0
              idy  = 2d0/g_def%dy(jg)
              if (ny == 1) idy = 0d0
              idz  = 2d0/g_def%dz(kg)
              if (nz == 1) idz = 0d0

              k2=vectorNorm(g_def,i,j,k,igr,(/idx,idy,idz/),.true.)
              kk=sqrt(k2)

              !|B|*k_par
              bk_par=scalarProduct(g_def,i,j,k,igr,(/idx,idy,idz/)
     .                           ,bcnv(i,j,k,:))

              !h_res = C*rho*v_a*k||/k^3/di (This will be multiplied by di^2 through ve and di*Ehall)
              h_res(i,j,k) = 5d0*sqrt(rho(i,j,k))*abs(bk_par)/kk**3/di
            enddo
          enddo
        enddo

        if (itime == 0) then
          maxv = maxval(h_res)
#if defined(petsc)
          call MPI_Reduce(maxv,maxvg,1,MPI_DOUBLE_PRECISION
     .                   ,MPI_MAX,0,MPI_COMM_WORLD,mpierr)
#else
          maxvg = maxv
#endif
          if (my_rank == 0) write (*,*) 'Hyperresistivity=',maxvg
        endif

      endif

c     End

      end function hres

c     chi_perp_ijk
c     #################################################################
      function chi_perp_ijk(g_def,i,j,k,igrid,te,rho,bcnv) result(dfvty)
c     -----------------------------------------------------------------
c     This function computes the local value of the heat diffusivity on
c     the node (i,j,k) of grid igrid according to field.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igrid
      real(8) :: te,rho,bcnv(3),dfvty
      type(grid_mg_def),pointer :: g_def

c     Local variables

      real(8) :: b2
! Daniele, 02/07/2014
      integer :: ig,jg,kg
      real(8) :: rr,th,zz

c     Begin program

      if (chi_bgnskii) then
        b2 = vectorNorm(g_def,i,j,k,igrid,bcnv,.false.)
        dfvty = min(chi*chi_perp_brag(te,te0_ref,rho,b2),1d-1)
      else
        select case (equil)
        case ('ppnch','ppnsl','ppnst','ppn3d','p3nsl','2dhel','ohtor')
          call find_polar_coords(g_def,igrid,i,j,k,rr,th,zz,coords)
          dfvty = chi*chi_rfx(rr)
        case default
          dfvty = chi
        end select
      endif

c     End program

      end function chi_perp_ijk

c     chi_perp_mesh
c     #############################################################
      function chi_perp_mesh(g_def,igrid,te,rho,bcnv) result(dfvty)
c     -------------------------------------------------------------
c     This function computes the diffusivity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: te(0:,0:,0:),rho(0:,0:,0:),bcnv(0:,0:,0:,:)
      real(8) :: dfvty(0:size(te,1)-1
     .                ,0:size(te,2)-1
     .                ,0:size(te,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = size(te,1)-2
      ny = size(te,2)-2
      nz = size(te,3)-2

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            dfvty(i,j,k)=chi_perp_ijk(g_def,i,j,k,igrid,te(i,j,k)
     .                               ,rho(i,j,k),bcnv(i,j,k,:))
          enddo
        enddo
      enddo

c     End program

      end function chi_perp_mesh

c     chi_par_ijk
c     #################################################################
      function chi_par_ijk(g_def,i,j,k,igrid,te,chi_per) result(dfvty)
c     -----------------------------------------------------------------
c     This function computes the local value of the heat diffusivity on
c     the node (i,j,k) of grid igrid according to field.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igrid
      real(8) :: te,dfvty,chi_per
      type(grid_mg_def),pointer :: g_def

c     Local variables

c     Begin program

      if (chi_bgnskii) then
        dfvty = max(chi_par*chi_par_brag(te,te0_ref),chi_per)
      else
c$$$        select case (equil)
c$$$        case ('ppnch','ppnsl','ppnst','ppn3d','p3nsl','2dhel','ohtor')
c$$$          call find_polar_coords(g_def,igrid,i,j,k,rr,th,zz,coords)
c$$$          dfvty = chi_par*chi_rfx(rr)
c$$$        case default
          dfvty = chi_par
c$$$        end select
      endif

c     End program

      end function chi_par_ijk

c     chi_par_mesh
c     #############################################################
      function chi_par_mesh(g_def,igrid,te,chi) result(dfvty)
c     -------------------------------------------------------------
c     This function computes the diffusivity on the grid level
c     igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: te(0:,0:,0:),chi(0:,0:,0:)
      real(8) :: dfvty(0:size(te,1)-1
     .                ,0:size(te,2)-1
     .                ,0:size(te,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = size(te,1)-2
      ny = size(te,2)-2
      nz = size(te,3)-2

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            dfvty(i,j,k)=chi_par_ijk(g_def,i,j,k,igrid,te(i,j,k)
     .                              ,chi(i,j,k))
          enddo
        enddo
      enddo

c     End program

      end function chi_par_mesh

c     diff_ijk
c     #################################################################
      function diff_ijk(g_def,i,j,k,igrid,te,rho,bcnv) result(dfvty)
c     -----------------------------------------------------------------
c     This function computes the local value of the particle diffusivity
c     on the node (i,j,k) of grid igrid according to fields te, rho.
c     -----------------------------------------------------------------

      implicit none

c     Call variables

      integer :: i,j,k,igrid
      real(8) :: te,rho,bcnv(3),dfvty
      type(grid_mg_def),pointer :: g_def

c     Local variables

      real(8) :: b2

c     Begin program

      if (chi_bgnskii) then
        b2 = vectorNorm(g_def,i,j,k,igrid,bcnv,.false.)
        dfvty = dd*diff_def(te,rho,b2,te0_ref)
      else
        dfvty = dd
      endif

c     End program

      end function diff_ijk

c     diff_mesh
c     #############################################################
      function diff_mesh(g_def,igrid,te,rho,bcnv) result(dfvty)
c     -------------------------------------------------------------
c     This function computes the particle diffusivity on the grid 
c     level igrid.
c     -------------------------------------------------------------

      implicit none

c     Call variables

      type(grid_mg_def),pointer :: g_def

      integer :: igrid
      real(8) :: te(0:,0:,0:),rho(0:,0:,0:),bcnv(0:,0:,0:,:)
      real(8) :: dfvty(0:size(te,1)-1
     .                ,0:size(te,2)-1
     .                ,0:size(te,3)-1)

c     Local variables

      integer :: nx,ny,nz,i,j,k

c     Begin program

      nx = size(te,1)-2
      ny = size(te,2)-2
      nz = size(te,3)-2

      do k=0,nz+1
        do j=0,ny+1
          do i=0,nx+1
            dfvty(i,j,k)=diff_ijk(g_def,i,j,k,igrid,te(i,j,k)
     .                           ,rho(i,j,k),bcnv(i,j,k,:))
          enddo
        enddo
      enddo

c     End program

      end function diff_mesh
      
      end module nlfunction_setup
