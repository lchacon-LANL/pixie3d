program pixeq_xfer

!===========================================================================
! PERFORMS SOLUTION TRANSFER FROM ONE CONFIGURATION TO ANOTHER.
!===========================================================================

  use var_io

  use grid

  use app_iosetup
  
  implicit none

  character(Len=1024) :: ifile,ofile,irecfile,orecfile
  integer :: n_args,ntimelevels,itlevel,order,u_rec_in
  integer :: ierr,perr
  logical :: debug=.true.,extrude_dir(3)

  !namelist
  namelist /xfer/ ifile,ofile,irecfile,orecfile,order

  !State variables
  type(var_array),pointer :: vref=>null(),vout=>null(),vref_0=>null()

  integer :: ic,jc,kc,if,jf,kf,ig,jg,kg,ivar,igx,igy,igz,grf,nxf,nyf,nzf
  integer :: it,itm,send_buf(3),rec_buf(3),nx,ny,nz
  real(8) :: gammat,dt,tt,volt,vol

  !Interpolation
  real(8) :: xp,yp,zp,interp
  real(8),allocatable,dimension(:) :: xx,yy,zz

  integer :: kx,ky,kz,nxs,nys,nzs,dim,flg,bcmod(6)

  real(8), dimension(:),allocatable:: tx,ty,tz,work
  real(8), dimension(:,:,:),allocatable:: bcoef,volf

!!$  INTERFACE
!!$    subroutine xfer_varray(vref,vout)
!!$      use variables
!!$      type(var_array),pointer :: vref,vout
!!$    end subroutine xfer_varray
!!$  END INTERFACE
    
  !Begin program

#if defined(petsc)
  call PetscInitialize(PETSC_NULL_CHARACTER,ierr)
  call initMPI(MPI_COMM_WORLD,np,my_rank)
#endif

  debug = debug.and.(my_rank == 0)

  igx = 1
  igy = 1
  igz = 1

  !Defaults
  order = 2         !Order of interpolation
  ifile = ''        !Input PIXIE3D file
  ofile = ''        !Output PIXIE3D file
  irecfile = ''     !Reference PIXIE3D record file
  orecfile = ''     !Output PIXIE3D record file

  open(unit=uinput,file='pixeq_xfer.in',status='old')
  read(uinput,xfer,iostat=ierr)
  close(uinput)

  !Read application grid setup configuration
!!$  call readGridInput(trim(ifile))
  in_file = trim(ifile)
  call readInput

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Open ORIGIN file and read in configuration!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (debug) write (*,*) "Reading REFERENCE file ",trim(ifile)

  u_rec_in = openRestartFileForRead(file=irecfile)

!!$  if (NDEPS > 0 .and. NDEPV > 0 .and. NAUXS > 0 .and. NAUXV > 0) then
!!$    call setDepVarDims(NDEPS,NDEPV)
!!$    call setAuxVarDims(NAUXS,NAUXV)
!!$  else
!!$    call pstop("pixeq_xfer","No variable info available")
!!$  endif

  !Set global limits
!!$#if defined(petsc)
!!$  write (*,*) ihig,jhig,khig,nxl,nyl,nzl
!!$  stop
!!$  send_buf = (/ ihig,jhig,khig /)
!!$  call MPI_Allreduce(send_buf,rec_buf,3 &
!!$                    ,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpierr)
!!$
!!$  nxd = rec_buf(1)
!!$  nyd = rec_buf(2)
!!$  nzd = rec_buf(3)
!!$#else
!!$  nxd = nxl
!!$  nyd = nyl
!!$  nzd = nzl
!!$#endif

  if (debug) then
    write (*,*)
    write (*,*) 'REFERENCE problem dimensions'
    write (*,*) 'NEQ=',neqd
    write (*,*) 'NX=',nxd
    write (*,*) 'NY=',nyd
    write (*,*) 'NZ=',nzd
    write (*,*)
  endif

  !Set vector dimensions and allocate variables
  call allocateGlobalVar(gv)

  call allocateDerivedType(vref)
 
  call allocateDerivedType(vref_0)

  !Read first file until last time slice
  call readTimeStep(u_rec_in,itm,tt,dt,gammat,vref_0,ierr)

  if (ierr /= 0) call pstop("pixeq_xfer","Equilibrium file unreadable")

  do 
    call readTimeStep(u_rec_in,itm,tt,dt,gammat,vref,ierr)

    if (ierr /= 0) exit
  enddo

  call closeRestartFileForRead(u_rec_in)

  call finalize_IO
  
  !Init PIXIE3D variables
!  call setAppBCs(vref,gv%aux)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Spline setup of reference solution!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nxs = nxl+2
  nys = nyl+2
  nzs = nzl+2

  allocate(xx(nxs),yy(nys),zz(nzs))

  call getMGmap(gv%gparams,1,1,1,igx,igy,igz,ig,jg,kg)

  xx(1:nxs) = gv%gparams%xx(ig-1:ig+nxl)
  yy(1:nys) = gv%gparams%yy(jg-1:jg+nyl)
  zz(1:nzs) = gv%gparams%zz(kg-1:kg+nzl)

  flg = 0

  kx = min(order+1,nxs-1)
  ky = min(order+1,nys-1)
  kz = min(order+1,nzs-1)

  dim = nxs*nys*nzs + max(2*kx*(nxs+1),2*ky*(nys+1),2*kz*(nzs+1))

  allocate(tx(nxs+kx))
  allocate(ty(nys+ky))
  allocate(tz(nzs+kz))
  allocate(work(dim))
  allocate(bcoef(nxs,nys,nzs))

  !!!!!!!!!!!!!!
  !Perform Xfer!
  !!!!!!!!!!!!!!

  nullify(gv%gparams)
  call destroyGrid(gv%gparams)
  call deallocateGlobalVar(gv)

  if (debug) write (*,*) "Performing STATE transfer"

  !Read application grid setup configuration
  nxf = nxd
  nyf = nyd
  nzf = nzd
  
  in_file = trim(ofile)
  call readInput

  extrude_dir(1) = (nxf == 1) .and. (nxd /= 1)
  extrude_dir(2) = (nyf == 1) .and. (nyd /= 1)
  extrude_dir(3) = (nzf == 1) .and. (nzd /= 1)

  if (debug) then
    write (*,*)
    write (*,*) 'TARGET problem dimensions'
    write (*,*) 'NEQ=',neqd
    write (*,*) 'NX=',nxd
    write (*,*) 'NY=',nyd
    write (*,*) 'NZ=',nzd
    write (*,*)
    write (*,*) 'Extrude?',extrude_dir
  endif

  !Set vector dimensions and allocate variables
  call allocateGlobalVar(gv)

  call createGrid(nxd,nyd,nzd,xmin,xmax,ymin,ymax,zmin,zmax,gv%gparams)

  call setVectorDimensions

  call allocateDerivedType(vout)

  !Init PIXIE3D variables
!  call setAppBCs(vout,gv%aux)

  nx = nxl
  ny = nyl
  nz = nzl
  
  if (extrude_dir(1)) nx = 1
  if (extrude_dir(2)) ny = 1
  if (extrude_dir(3)) nz = 1

  !!!!!!!!!!!!!!!!!!!
  !Write output file!
  !!!!!!!!!!!!!!!!!!!

  if (debug) write (*,*) "Writing output ",trim(orecfile)

!!  if (debug) write (*,*) 'itime=',itm,'; time=',tt

  !Reset time counters
  itm = 0 ; tt = 0d0

  call xfer_varray(vref_0,vout)
!!  call applyBC(igx,vout,gv%aux)

#if defined(adios)
  ierr = init_ADIOS_IO()
#endif
  
  call writeRecordFile(orecfile,itm,tt,dt,gammat,vout,init=.true.)

  call xfer_varray(vref,vout)

  call writeRecordFile(orecfile,itm,tt,dt,gammat,vout,init=.false.)
  
!!  call writeDerivedType(vout,6,.true.)
  
  !!!!!!!!!!!!!!!!!!!!!!
  !Deallocate variables!
  !!!!!!!!!!!!!!!!!!!!!!

  call deallocateDerivedType(vout)
  call deallocateDerivedType(vref)
  call deallocateDerivedType(vref_0)

  deallocate(tx,ty,tz,work,bcoef,xx,yy,zz)

  call deallocateStructures
  call deallocateGlobalVar(gv)

#if defined(petsc)
  call PetscFinalize(mpierr)
#endif

contains

  ! xfer_varray
  !######################################################################
  subroutine xfer_varray(vref,vout)

    implicit none
    
    !Call variables
    type(var_array),pointer :: vref,vout

    !Local variables
    integer :: ivar

    real(8) :: db2val,db3val
    external   db2val,db3val

    real(8),pointer,dimension(:,:,:) :: arrayf,arrayc

    !Begin program
  
    do ivar=1,vout%nvar
      vout%array_var(ivar)%descr  &
             = vref%array_var(ivar)%descr

      vout%array_var(ivar)%bconds &
             = vref%array_var(ivar)%bconds

      vout%array_var(ivar)%bc_dep_list &
             = vref%array_var(ivar)%bc_dep_list
      vout%array_var(ivar)%dom_dep_list &
             = vref%array_var(ivar)%dom_dep_list

      arrayf => vref%array_var(ivar)%array
      arrayc => vout%array_var(ivar)%array

      !Spline reference component
      call db3ink(xx,nxs,yy,nys,zz,nzs,arrayf &
                 ,nxs,nys,kx,ky,kz,tx,ty,tz,bcoef,work,flg)

      !Interpolate to target grid
      do kc = 0,nz+1
        do jc = 0,ny+1
          do ic = 0,nx+1
            call getMGmap(gv%gparams,ic,jc,kc,igx,igy,igz,ig,jg,kg)

            xp = gv%gparams%xx(ig)
            yp = gv%gparams%yy(jg)
            zp = gv%gparams%zz(kg)

            interp = db3val(xp,yp,zp,0,0,0,tx,ty,tz,nxs,nys,nzs &
                           ,kx,ky,kz,bcoef,work)
!!            write (*,*) ic,jc,kc,xp,yp,zp,interp
           
            if     (extrude_dir(1).and.extrude_dir(2)) then
              arrayc(:,:,kc) = interp
            elseif (extrude_dir(1).and.extrude_dir(3)) then
              arrayc(:,jc,:) = interp
            elseif (extrude_dir(2).and.extrude_dir(3)) then
              arrayc(ic,:,:) = interp
            elseif (extrude_dir(1)) then
              arrayc(:,jc,kc) = interp
            elseif (extrude_dir(2)) then
              arrayc(ic,:,kc) = interp
            elseif (extrude_dir(3)) then
              arrayc(ic,jc,:) = interp
            else
              arrayc(ic,jc,kc) = interp
            endif

          enddo
        enddo
      enddo

    enddo

    nullify(arrayf,arrayc)

  end subroutine xfer_varray
  
end program pixeq_xfer

