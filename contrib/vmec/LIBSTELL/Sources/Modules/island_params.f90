     MODULE island_params
     USE stel_kinds
!
!     CONTAINS VARIABLE DECLARATIONS FOR VMECPP (ISLAND) CODE
!     THAT WILL BE SHARED THROUGHOUT THE PROJECT
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
     INTEGER :: ns_i, nu_i, nv_i, nuv_i                                
!size of s, u, v dimensions of (island) metric tensors
     INTEGER :: mpol_i, ntor_i
     INTEGER :: nfp_i, mnmax_i
     INTEGER(iprec) :: nsh                                            
! Number of points in half mesh

     REAL(rprec) :: hs_i, ohs_i, dnorm_i, wb_i
     REAL(rprec), PARAMETER :: gamma =5._dp/3._dp                       ! Adiab. constant. Move to ISLAND_PARAMS?
     REAL(rprec), ALLOCATABLE, DIMENSION(:,:) :: cosmu, cosmum,        &
    &    cosmui, sinmu, sinmum, sinmui, cosnv, cosnvn, sinnv, sinnvn
     REAL(rprec), PARAMETER:: nfactor = 1.41421356237310d0            
! RS: Normalization of m=0/n=0 Fourier modes.
!      REAL(rprec), PARAMETER:: nfactor = 1.0d0
     REAL(rprec), ALLOCATABLE, DIMENSION(:)   :: phipf_i, iotaf_i,presf_i, vp_h
!-----------------------------------------------

     END MODULE island_params 
