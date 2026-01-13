! PFBplotter
! Converts parflow PFB formatted single variable outputs into a binary
! VTK format.
!
! The user guide for this describes how to compile and run so if
!   someone gave you this source code without the manual file
!   "PFBplotter_manual.pfb", you might call that person a jerk,
!   but I'd suggest asking them nicely to send the manual first.
!
! Version 0.1 - April 9, 2014
!   It's kind of a beta I guess
!
! Nick Engdahl (nengdahl@mines.edu)
! Colorado School of Mines
!
! ------------------------------------------------------------------------------
module global
! Sets up the globally available variables
! >>>>> A few options you might need to change at some point <<<<<
    integer,parameter :: Mode = 6
    logical,parameter :: ask_for_params = .false.

    logical,parameter :: use_lower_t = .true.
    logical,parameter :: use_upper_t = .true.
    logical,parameter :: use_extension = .false.

    integer,parameter :: clen=200
    logical,parameter :: fast_read = .true.
    integer,parameter :: interp_mode = 1
    integer,parameter :: num_len = 5


! The rest of this module should probably be left alone...
type gridinfo
    integer :: nx=1,ny=1,nz=1
    real*8 :: dx=1.0,dy=1.0,dz=1.0
    real*8 :: x0=0.0,y0=0.0,z0=0.0
end type gridinfo

    type(gridinfo) :: MstrGrid,OutGrid
    character(len=10) ext
    character(len=clen) rootname,MSKFile,Varname,DEMfile
    integer timezero,timelast
    real,allocatable :: Pnts(:,:)
    integer, allocatable :: Lnum(:,:,:),TopCell(:,:),BotCell(:,:)
    real*8,allocatable :: PFBdata(:,:,:),DEM(:,:,:),MASK(:,:,:),OutData(:,:,:)
    real*8 misval,mismsk
    logical use_a_mask,use_a_dem,grid_tfg,num_lyrs
    character(len=60) :: fmtstr

end module global
! ------------------------------------------------------------------------------
program pfhillout
use global
implicit none
integer i,j,k,n,nn,npnts,ncell
integer steps,IO,outs
character(len=10) extin
character(len=clen) testname,tmpnam,pfbfile,vtkfile
logical tmp_lgc

! ==============================================================================
! ==============================================================================
!           >>>>> User configurable options set at compile time <<<<<
!
! This section allows you to customize the inputs the program will require
! If an option is flagged .true. it will be asked for and/or assumed

! Descriptions:

! use_a_mask:       Look for and attempt to use a mask for your domain
! num_lyrs:         Add an attribute for the layer number in descending order
!                    (layer 1 is the top), requires a mask file
! grid_tfg:         Is grid (input or output) terrain following? If a mask is also used
!                    then the output is a structured grid, not structured points.
! use_a_dem:        Read in a DEM for the elevations. Only used if grid_tfg=.true.
! use_lower_t:      Specify integer number of the first time step to try (defaults to 0)
! use_upper_t:      Specify integer number of the last time step to try (defaults to 1e4)
! use_extension:    Specify an extension other than .pfb
! ask_for_params:   Prompt user for required parameters (dynamic) or use hard-wired parameters
!                    I put all the required inputs into a text file and pipe the contents in
! interp_mode:      Option for DEM interpolation mode. Doesn't do much right now.
! num_length:       How many digits make up the time series numbers? Defauts to 5, the parflow default
! fast_read:        Direct (.true.) or loop (.false.) reading of a PFB. Loop will always work but is SLOW


select case(Mode)
case(0)
! Standard grid, no mask
use_a_mask = .false.
num_lyrs = .false.
grid_tfg = .false.
use_a_dem = .false.
case(1)
! Standard grid, with mask, direct conversion
use_a_mask = .true.
num_lyrs = .false.
grid_tfg = .false.
use_a_dem = .false.
case(2)
! Standard grid, with mask, add layer numbers
use_a_mask = .true.
num_lyrs = .true.
grid_tfg = .false.
use_a_dem = .false.
case(3)
! Convert to terrain following, compute top elevations
use_a_mask = .true.
num_lyrs = .false.
grid_tfg = .true.
use_a_dem = .false.
case(4)
! Convert to terrain following, read elevations from DEM
use_a_mask = .true.
num_lyrs = .false.
grid_tfg = .true.
use_a_dem = .true.
case(5)
! Convert to terrain following, read elevations from DEM
use_a_mask = .false.
num_lyrs = .false.
grid_tfg = .true.
use_a_dem = .false.
case(6)
! Convert to terrain following, read elevations from DEM
use_a_mask = .false.
num_lyrs = .false.
grid_tfg = .true.
use_a_dem = .true.
case default
write(*,*) 'ERROR: Invalid operation mode'
stop

end select

!                   >>>>>> End of user options <<<<<<
! ==============================================================================
! ==============================================================================

if (ask_for_params .eqv. .false.) then
! -------------------------------------------------------------------------
! Internal parameters
!
! >>>>>> If ask_for_params=.false. the following values will be used <<<<<<
!
! But the parameters are only used if the corresponding flags above are set

rootname='neon_sa'      ! Root name of files, without extension
Varname='LULC'                  ! Name for variable
MSKFile='East_Inlet.out.mask.pfb'    ! Mask file, with .pfb extension
DEMfile='ER_dem.pfb'  ! DEM file, with .pfb extension
timezero=0                        ! First time to search for (integer > 0)
timelast=0                        ! Last time to search for
extin='.pfb'                        ! Extension that goes with "rootname"
! -------------------------------------------------------------------------
else
! Set a few defaults
timezero=0                        ! First time to search for (integer > 0)
timelast=1e4                      ! Last time to search for
extin='.pfb'
endif

if (ask_for_params) then
write(*,'(A)') 'PFB root name (e.g. runname.out.press):'
read(*,*) rootname

write(*,'(A)') 'Name for variable:'
read(*,*) Varname

if (use_a_mask) then
write(*,'(A)') 'Input Mask file (include the .pfb extension):'
read(*,*) MSKfile
endif

if (use_a_dem) then
write(*,'(A)') 'Input DEM file (include the .pfb extension):'
read(*,*) DEMfile
endif

if (use_lower_t) then
write(*,'(A)') 'First time step (t>=0):'
read(*,*) timezero
endif

if (use_upper_t) then
write(*,'(A)') 'Last time step:'
read(*,*) timelast
endif

if (use_extension) then
write(*,'(A)') 'Extension for files:'
read(*,*) extin
endif
! End of parameter specification
endif

!  -------- Some things that are always done --------

! These first two might need adjustment depending on different machines
mismsk = 0.0d+0         ! Inactive cell flag used in mask files
misval = -1.0e+37       ! Huge negative used to flag inactive cells in PFB

timelast=timelast + 1 ! Correct counting for time zero offset

! Set up the format string for the filenames
write(fmtstr,*) '(2A,I',num_len,'.',num_len,',A)'
ext=adjustl(trim(extin))

! check for the first file and pull the grid from it
write(tmpnam,fmtstr) adjustl(trim(rootname)),'.',timezero,ext
write(testname,'(A)') adjustl(trim(tmpnam))
write(*,*) testname

call get_grid(testname)

! Allocate the storage array
allocate(PFBdata(MstrGrid%nx,MstrGrid%ny,MstrGrid%nz))
PFBdata(:,:,:) = misval     ! Default to the missing value

! ---------- See how many files we've got to convert ----------
write(*,'(A,I6)') 'Time zero is: ',timezero
call nsteps(steps)
write(*,'(I6,A)') steps,' time steps found'

!  -------- Option specific tasks --------
! Check for the mask file, if it doesn't exist go on without it
tmp_lgc = use_a_mask
if (tmp_lgc) then
open(unit=50,file=trim(adjustl(MSKFile)),form='unformatted', &
access='stream',convert='BIG_ENDIAN',status='old',iostat=IO)
if (IO .ne. 0) then
write(*,*) 'WARNING: Mask not found, proceeding without it'
use_a_mask = .false.
endif
close(50)
endif

if (grid_tfg) then
! Check for the DEM file, if it doesn't exist go on without it
tmp_lgc = use_a_dem
if (tmp_lgc) then
open(unit=50,file=trim(adjustl(DEMfile)),form='unformatted', &
access='stream',convert='BIG_ENDIAN',status='old',iostat=IO)
if (IO .ne. 0) then
write(*,*) 'WARNING: DEM not found, proceeding without it'
use_a_dem = .false.
endif
close(50)
endif
else
! Uh, there's no need for a DEM if you're going all rectilinear on me
use_a_dem = .false.
endif

tmp_lgc = num_lyrs
if (grid_tfg .and. tmp_lgc) then
! And there's no need to compute layer numbers for a TFG
num_lyrs = .false.
endif

allocate(MASK(MstrGrid%nx,MstrGrid%ny,MstrGrid%nz))
if (use_a_mask) then
! OK, if we've got a mask, try to read it in
write(*,*) 'Reading mask...'
call readMASK()
else
! Define it anyway to a default
MASK(:,:,:)=1.0d+0
endif

! This only needs to be done if there are masked cells
if ((minval(MASK) .gt. 0.999d+0) .or. (use_a_mask .eqv. .false.)) then
! Huge time saver
if (num_lyrs) then
allocate(TopCell(MstrGrid%nx,MstrGrid%ny))
allocate(BotCell(MstrGrid%nx,MstrGrid%ny))
TopCell(:,:) = MstrGrid%nz
BotCell(:,:) = 1
endif

else
! Go get the top and bottom

! If the output is requested/needed compute the top of the model
if (grid_tfg .or. num_lyrs) then
! It should be impossible to get here if these are already allocated but just to be safe...
if (allocated(TopCell) .eqv. .false.) then
allocate(TopCell(MstrGrid%nx,MstrGrid%ny))
endif
if (allocated(BotCell) .eqv. .false.) then
allocate(BotCell(MstrGrid%nx,MstrGrid%ny))
endif

do j=1,MstrGrid%ny
  do i=1,MstrGrid%nx
   k=MstrGrid%nz
   do while ((MASK(i,j,k) .lt. 1.0) .and. (k .gt. 1))
   k=k-1
   enddo
    if ((k .eq. 1) .and. (MASK(i,j,k) .lt. 1.0))then
    ! Just some big negative integers to mark inactives
    TopCell(i,j) = -99999
    BotCell(i,j) = -88888
    else
    TopCell(i,j) = k

    n=1
    do while ((MASK(i,j,n) .lt. 1.0) .and. (n .lt. MstrGrid%nz))
    n=n+1
    enddo
    if (n .le. TopCell(i,j)) then
    BotCell(i,j) = n
    else
    BotCell(i,j) = -77777
    endif
    endif

  enddo ! loop on i
enddo ! loop on j

endif ! (grid_tfg .or. num_lyrs)
endif ! minval(MASK) test

! -------------------------------------------------------------------------------
!
! Setup is now done, moving on to the actual data read/conversion/write section
!
!   All the different cases are handled similarly with a few little changes
!   to give the desired output.
!
! -------------------------------------------------------------------------------

if (use_a_mask .and. grid_tfg) then
write(*,'(A)') ' -- Writing shifted, masked, TFG --'
! =====> Converting a masked rectilinear grid into a terrain following grid <=====
!  Uses DEM to compute top, flag should be 1
call compress_grid(1)

! Generate points for structured_grid output
if (use_a_dem) then
write(*,*) ' --- Reading DEM'
call readDEM()
else
write(*,*) '  --- Making DEM'
call makeDEM()
endif

call compute_points(npnts,ncell)

do outs=timezero+1,timezero+steps
! Loop over all the files
n=outs-1
write(pfbfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,ext
write(vtkfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,'.vtk'

! Initialize/blank this array every time
PFBdata(:,:,:)=0
call pfb_read(pfbfile)

! Shift it, flag should be 1
call shift_grid(1)

call struct_VTK(vtkfile,npnts,ncell)
enddo

elseif (use_a_mask) then
write(*,'(A)') ' -- Writing masked structured_points --'
! =====> Just re-write the standard grid as structured_points <=====
call compress_grid(0)

! --------------------------------------
! Compute layer numbers

! Allocate the layer numbers array
if (num_lyrs) then

if (OutGrid%nz .ne. MstrGrid%nz) then
write(*,*) 'ERROR: Grids not aligned.'
stop
endif

allocate(Lnum(OutGrid%nx,OutGrid%ny,OutGrid%nz))
Lnum(:,:,:)=0

do i=1,OutGrid%nx
 do j=1,OutGrid%ny
  n=1
  do k=OutGrid%nz,1,-1
   if (Mask(i,j,k) .gt. 0.5) then
   Lnum(i,j,k)=n
   n=n+1
   endif
  enddo

 enddo
enddo

endif
! --------------------------------------

! Time loop
do outs=timezero+1,timezero+steps
n=outs-1
write(pfbfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,ext
write(vtkfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,'.vtk'

! Initialize/blank this array every time
PFBdata(:,:,:)=0
call pfb_read(pfbfile)

call shift_grid(0)
call VTK_SPnts(vtkfile)

enddo ! End time loop

elseif (grid_tfg) then
! Turns out there are two cases in this one...

if (use_a_dem) then
write(*,'(A)') ' -- Writing TFG as structured_grid --'
! ====> Write terrain following grid as structured_grid <=====
!  Uses DEM to compute top
call compress_grid(1)
! NBE: this could be a slight problem, the shift won't work if the dimensions are different...
!       I added a check to the shifter for this. It would only happen if a terrain following
!       grid in parflow is masked such that ACTIVE EXTENTS are different from the grid extents
!   Solution: just run it without a mask

! Generate points for structured_grid output
call readDEM()
call compute_points(npnts,ncell)

! Time loop
do outs=timezero+1,timezero+steps
n=outs-1
write(pfbfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,ext
write(vtkfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,'.vtk'

! Initialize/blank this array every time
PFBdata(:,:,:)=0
call pfb_read(pfbfile)

call shift_grid(0)
call struct_VTK(vtkfile,npnts,ncell)

enddo ! End time loop

else
write(*,'(A)') ' -- Writing TFG as structured_point grid --'
! ====> Write terrain following grid as structured_points <=====
!  Uses DEM to compute top
call compress_grid(0)

! Time loop
do outs=timezero+1,timezero+steps
n=outs-1
write(pfbfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,ext
write(vtkfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,'.vtk'

! Initialize/blank this array every time
PFBdata(:,:,:)=0
call pfb_read(pfbfile)

call shift_grid(0)
call struct_VTK(vtkfile,npnts,ncell)

enddo ! End time loop
endif ! (use_a_dem)

elseif ((use_a_dem .eqv. .false.) .and. (use_a_dem .eqv. .false.) .and. &
    (use_a_dem .eqv. .false.) .and. (use_a_dem .eqv. .false.)) then
write(*,'(A)') ' -- Writing structured_points --'
! =====> Just re-write the standard grid as structured_points <=====
call compress_grid(0)

! Time loop
do outs=timezero+1,timezero+steps
n=outs-1
write(pfbfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,ext
write(vtkfile,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,'.vtk'

! Initialize/blank this array every time
PFBdata(:,:,:)=0
call pfb_read(pfbfile)

call shift_grid(0)
call VTK_SPnts(vtkfile)

enddo ! End time loop

else
write(*,*) 'ERROR: Invalid output configuration'
write(*,*) '       Sorry, you put one by me there!'
stop

endif

write(*,'(A)') ' --- Conversions complete ---'

end program
! -----------------------------------------------------------------------------------------
!                       ******** SUBROUTINE NATION ********
! -----------------------------------------------------------------------------------------
subroutine readMASK()
! Note: readMASK and readDEM are basically pfb_read, but I keep
!        them seperate anyway since the data goes to a module
use global
implicit none
integer nxa,nya,nza,ns
real*8 dx, dy, dz, x1, y1, z1
integer i,j,k, nni, nnj, nnk, ix, iy, iz
integer rx, ry, rz, nnx, nny, nnz,is,io,MSKF
MSKF=49

open(MSKF,file=trim(MSKfile),form='unformatted', access='stream', &
convert='BIG_ENDIAN',status='old',iostat=IO)
if (IO .ne. 0) then
! This check should be imposible to trigger, but I'll leave it in anyway
write(*,'(A)') 'ERROR: Mask file could not be opened. Check file name and try again'
stop
endif

! Main grid definition
read(MSKF) x1
read(MSKF) y1
read(MSKF) z1
read(MSKF) nxa
read(MSKF) nya
read(MSKF) nza
read(MSKF) dx
read(MSKF) dy
read(MSKF) dz
read(MSKF) ns

if ((nxa .ne. MstrGrid%nx) .or. (nya .ne. MstrGrid%ny) .or. (nza .ne. MstrGrid%nz)) then
write(*,'(A)') 'ERROR: Inconsistent data and mask dimensions'
write(*,'(A,3I8)') '  Grid:', MstrGrid%nx,MstrGrid%ny,MstrGrid%nz
write(*,'(A,3I8)') '  Mask:',nxa,nya,nza
stop
endif

do is = 0, (ns-1)
read(MSKF) ix
read(MSKF) iy
read(MSKF) iz
read(MSKF) nnx
read(MSKF) nny
read(MSKF) nnz
read(MSKF) rx
read(MSKF) ry
read(MSKF) rz

if (fast_read) then
read(MSKF) MASK(ix+1:ix+nnx,iy+1:iy+nny,iz+1:iz+nnz)
else
! This is a SLOW read operation but it will always work
do  k=iz +1 , iz + nnz
 do  j=iy +1 , iy + nny
  do  i=ix +1 , ix + nnx
    read(MSKF) MASK(i,j,k)
  end do
 end do
end do
endif

end do
close(MSKF)
end subroutine
! ---------------------------------------------------------------------------------------------------------------------------------------
subroutine readDEM()
use global
implicit none
integer nxa,nya,nza,ns,io
real*8 dx, dy, dz, x1, y1, z1
real*8,allocatable :: datemp(:,:,:)
integer i,j,k, nni, nnj, nnk, ix, iy, iz
integer rx, ry, rz, nnx, nny, nnz,is,DEMF
DEMF=51

open(DEMF,file=trim(DEMfile),form='unformatted', access='stream', &
convert='BIG_ENDIAN',status='old',iostat=IO)
if (IO .ne. 0) then
! This check should be imposible to trigger, but I'll leave it in anyway
write(*,'(A)') 'ERROR: DEM file could not be opened. Check file name and try again'
stop
endif

! Main grid definition
read(DEMF) x1
read(DEMF) y1
read(DEMF) z1
read(DEMF) nxa
read(DEMF) nya
read(DEMF) nza
read(DEMF) dx
read(DEMF) dy
read(DEMF) dz
read(DEMF) ns ! # of subgrids

if ((OutGrid%nx .ne. nxa) .or. (OutGrid%ny .ne. nya)) then
write(*,*) 'ERROR: DEM size differs from PFB data'
write(*,*) ''
write(*,*) "PFB size, DEM size:"
write(*,*) 'nx',OutGrid%nx,nxa
write(*,*) 'ny',OutGrid%ny,nya
stop
endif

!! Allocate the temporary data array, in case nz>1
allocate(datemp(nxa,nya,nza))

do is = 0, (ns-1)
! Sub-grid header
read(DEMF) ix
read(DEMF) iy
read(DEMF) iz
read(DEMF) nnx
read(DEMF) nny
read(DEMF) nnz
read(DEMF) rx
read(DEMF) ry
read(DEMF) rz
! Now grab the data
if (fast_read) then
read(DEMF) datemp(ix+1:ix+nnx,iy+1:iy+nny,iz+1:iz+nnz)
else
do  k=iz +1 , iz + nnz
 do  j=iy +1 , iy + nny
  do  i=ix +1 , ix + nnx
    read(DEMF) datemp(i,j,k)
  end do
 end do
end do
endif
end do

! Now store the DEM permanently in the global memory
!  Hang on to the TOP LAYER ONLY for this data
allocate(DEM(OutGrid%nx,OutGrid%ny,1))
DEM(:,:,1)=datemp(:,:,nza)
deallocate(datemp)

close(DEMF)
end subroutine
! ---------------------------------------------------------------------------------------------------------------------------------------
subroutine makeDEM()
use global
implicit none
integer io,i,j,MaxH,nx,ny,top
integer MSKF, elev(8),cnt

MSKF=49

allocate(DEM(OutGrid%nx,OutGrid%ny,1))

open(MSKF,file=trim(MSKFile),form='unformatted', access='stream', &
convert='BIG_ENDIAN',status='old',iostat=IO)
if (IO .ne. 0) then
write(*,'(A)') 'WARNING: Mask file not found. Using dz*nz for DEM'

DEM(:,:,1) = Mstrgrid%z0 + DBLE(MstrGrid%nz)*MstrGrid%dz
else
close(MSKF)
! Set a default value based on the UNCOMPRESSED grid
MaxH = nint(real(MstrGrid%nz)/real(2))

nx=OutGrid%nx
ny=OutGrid%ny

do j=1,ny
  do i=1,nx
if ((TopCell(i,j) .le. 0) .and. (BotCell(i,j) .le. 0)) then
elev(:) = 0
cnt = 0
top = 0
! Check for nearby active cells, being careful to avoid seg faults...
if (i+1 .le. nx) then
if (TopCell(i+1,j) .gt. 0) then
elev(1) = TopCell(i+1,j)
cnt = cnt + 1
endif
endif
if (i-1 .ge. 1) then
if (TopCell(i-1,j) .gt. 0) then
elev(2) = TopCell(i-1,j)
cnt = cnt + 1
endif
endif
if (j+1 .le. ny) then
if (TopCell(i,j+1) .gt. 0) then
elev(3) = TopCell(i,j+1)
cnt = cnt + 1
endif
endif
if (j-1 .ge. 1) then
if (TopCell(i,j-1) .gt. 0) then
elev(4) = TopCell(i,j-1)
cnt = cnt + 1
endif
endif
if ((i+1 .le. nx) .and. (j+1 .le. ny)) then
if (TopCell(i+1,j+1) .gt. 0) then
elev(5) = TopCell(i+1,j+1)
cnt = cnt + 1
endif
endif
if ((i-1 .ge. 1) .and. (j+1 .le. ny)) then
if (TopCell(i-1,j+1) .gt. 0) then
elev(6) = TopCell(i-1,j+1)
cnt = cnt + 1
endif
endif
if ((i+1 .le. nx) .and. (j-1 .ge. 1)) then
if (TopCell(i+1,j-1) .gt. 0) then
elev(7) = TopCell(i+1,j-1)
cnt = cnt + 1
endif
endif
if ((i-1 .ge. 1) .and. (j-1 .ge. 1)) then
if (TopCell(i-1,j-1) .gt. 0) then
elev(8) = TopCell(i-1,j-1)
cnt = cnt + 1
endif
endif

if (cnt .gt. 0) then
! We're an edge cell of the mask, add a value
top = NINT(real(sum(elev))/real(cnt))
DEM(i,j,1)=Mstrgrid%z0 + MstrGrid%dz*DBLE(top)
else
! Not near edge, inactive cell. Give 'er the default
DEM(i,j,1)=Mstrgrid%z0 + MstrGrid%dz*DBLE(MaxH)
endif
else
! Normal case, not an edge cell
DEM(i,j,1)=Mstrgrid%z0 + MstrGrid%dz*DBLE(TopCell(i,j))
endif

  enddo
enddo
endif

end subroutine
! -----------------------------------------------------------------------------------
subroutine compute_points(npnts,ncell)
use global
implicit none
integer n,steps,IO,m,outs,t
integer nx,ny,nz,i,j,k,ii,jj
integer nnx,nny,nnz,npnts,ncell
real x1,y1,z1,dx,dy,dz,maxz,zoff,tval
real,allocatable :: PntFlip(:,:)

nx=OutGrid%nx
ny=OutGrid%ny
nz=OutGrid%nz

dx=real(OutGrid%dx)
dy=real(OutGrid%dy)
dz=real(OutGrid%dz)

x1=real(OutGrid%x0)
y1=real(OutGrid%y0)
z1=real(OutGrid%z0)

nnx=nx+1
nny=ny+1
nnz=nz+1
npnts=nnx*nny*nnz
ncell=nx*ny*nz

allocate(Pnts(npnts,3))
Pnts(:,:)=0
m=1

! Duh, need the maximum height of the model
maxz=real(nz)*dz

! Should probably do bi-linear interpolation or something more
!  robust, but for visualization I'm OK with this for now

if (1 .eq. 0) then
! This is a simple way of handling the maximum edges
do k=1,nnz
do j=1,nny
do i=1,nnx
Pnts(m,1)=x1+real(i-1)*dx
Pnts(m,2)=y1+real(j-1)*dy
if (i .le. nx) then
ii=i
else
ii=nx
endif
if (j .le. ny) then
jj=j
else
jj=ny
endif
! This step translates the DEM
! The specified initial heights in the pfb (z1) are ignored and the
!  offset is computed based on the model thickness
zoff=DEM(ii,jj,1)-maxZ
Pnts(m,3)=zoff+real(k-1)*dz
m=m+1
enddo
enddo
enddo

else  ! "simple handling" switch
! Use the average of the surrounding nodes for the elevation
do k=1,nnz
do j=1,nny
do i=1,nnx
Pnts(m,1)=x1+real(i-1)*dx
Pnts(m,2)=y1+real(j-1)*dy

t=0
tval=0.0d+0

! Simple case, we're in the middle of the field
if ((i .gt. 1) .and. (i .lt. nnx) .and. (j .gt. 1) &
.and. (j .lt. nny)) then
tval=tval + DEM(i-1,j-1,1) + DEM(i,j,1) &
+ DEM(i-1,j,1) + DEM(i,j-1,1)
t=4
! The rest of the these logically handle the edges
elseif ((i .gt. 1) .and. (i .lt. nnx) .and. (j .eq. 1)) then
tval=tval + DEM(i-1,j,1) + DEM(i,j,1)
t=2
elseif ((i .gt. 1) .and. (i .lt. nnx) .and. (j .eq. nny)) then
tval=tval + DEM(i-1,j-1,1) + DEM(i,j-1,1)
t=2
elseif ((j .gt. 1) .and. (j .lt. nny) .and. (i .eq. 1)) then
tval=tval + DEM(i,j-1,1) + DEM(i,j,1)
t=2
elseif ((j .gt. 1) .and. (j .lt. nny) .and. (i .eq. nnx)) then
tval=tval + DEM(i-1,j-1,1) + DEM(i-1,j,1)
t=2
elseif ((i .eq. 1).and.(j.eq.1)) then
tval=DEM(i,j,1)
t=1
elseif ((i.eq.nnx).and.(j.eq.nny)) then
tval=DEM(i-1,j-1,1)
t=1
elseif ((i.eq.1).and.(j.eq.nny)) then
tval=DEM(i,j-1,1)
t=1
elseif ((i.eq.nnx).and.(j.eq.1)) then
tval=DEM(i-1,j,1)
t=1
else
write(*,*) 'Something went terribly wrong, oops.'
stop
endif

zoff = tval / real(t) - maxZ
Pnts(m,3)=zoff+real(k-1)*dz
m=m+1
enddo
enddo
enddo

endif

! Reorient the points for the VTK output order
allocate(PntFlip(3,npnts))
PntFlip(:,:)=0.0

do i=1,3
PntFlip(i,:)=real(Pnts(:,i))
enddo
deallocate(Pnts)

allocate(Pnts(3,npnts))
Pnts(:,:)=PntFlip(:,:)
deallocate(PntFlip)

end subroutine
! -----------------------------------------------------------------------------------
subroutine pfb_read(pfbname)
! Modified after pfb_read in /pftools/prepostproc
use global
implicit none
integer nxa,nya,nza
real*8 dx, dy, dz, x1, y1, z1                                                           
integer i,j,k, nni, nnj, nnk, ix, iy, iz,IO
integer ns,  rx, ry, rz, nnx, nny, nnz,is, PFID 
character(len=100),intent(IN) :: pfbname
PFID=50

open(PFID,file=trim(adjustl(pfbname)),form='unformatted', & 
access='stream',convert='BIG_ENDIAN',status='old',iostat=IO) 
if (IO .ne. 0) then
write(*,*) 'ERROR: problem opening PFB'
stop
endif
! Main grid definition
read(PFID) x1
read(PFID) y1
read(PFID) z1
read(PFID) nxa
read(PFID) nya
read(PFID) nza
read(PFID) dx
read(PFID) dy
read(PFID) dz
read(PFID) ns ! # of subgrids

do is = 0, (ns-1)
! Sub-grid header
read(PFID) ix
read(PFID) iy
read(PFID) iz
read(PFID) nnx
read(PFID) nny
read(PFID) nnz
read(PFID) rx
read(PFID) ry
read(PFID) rz
!! Now grab the data
if (fast_read) then
read(PFID) PFBdata(ix+1:ix+nnx,iy+1:iy+nny,iz+1:iz+nnz)
else
do  k=iz +1 , iz + nnz
 do  j=iy +1 , iy + nny
  do  i=ix +1 , ix + nnx
    read(PFID) PFBdata(i,j,k)
  end do
 end do
end do
endif

end do
close(PFID)
end subroutine
! -----------------------------------------------------------------------------------
subroutine nsteps(steps)
! This tests for and determines the number of time steps to be processed
! rootname is usedfor this and the filenames MUST be sequential
use global
implicit none
integer, intent(OUT) :: steps
integer i,n,IO,maxn
character(len=100) testname,tmpnam

maxn=timelast
maxn=maxn-1 ! This corrects the counting later on

IO=0
n=timezero-1
do while ((n.le.maxn).and.(IO .eq. 0))
n=n+1
write(tmpnam,'(2A,I5.5,A)') adjustl(trim(rootname)),'.',n,ext
write(testname,'(A)') adjustl(trim(tmpnam))
open(unit=1,file=trim(testname),status='old',action='read',iostat=IO)
if (IO .ne. 0) then
! File isn't there so bail out
exit
endif
if (IO .eq. 0) close(1)
enddo

if (n.eq.0) then
write(*,'(A)') 'ERROR: No time steps found. Check root name and try again.'
stop
else
endif

steps=n-timezero
! NOTE: The output is the NUMBER of steps!
end subroutine
! --------------------------------------------------------------------------------
subroutine get_grid(infile)
use global
implicit none
character(len=clen),intent(in) :: infile
integer IO,PFID
IO=0
PFID=50

open(PFID,file=trim(adjustl(infile)),form='unformatted', &
access='stream',convert='BIG_ENDIAN',status='old',iostat=IO)
if (IO .ne. 0) then
write(*,*) 'ERROR: First data file not found'
stop
endif

read(PFID) MstrGrid%x0
read(PFID) MstrGrid%y0
read(PFID) MstrGrid%z0
read(PFID) MstrGrid%nx
read(PFID) MstrGrid%ny
read(PFID) MstrGrid%nz
read(PFID) MstrGrid%dx
read(PFID) MstrGrid%dy
read(PFID) MstrGrid%dz

close(PFID)
end subroutine
! --------------------------------------------------------------------------------
subroutine shift_grid(flag)
! Shifts input data to the TOP of the output data array if flag == 1
! Otherwise it just does a direct copy of the grid
use global
implicit none
integer,intent(IN) :: flag
integer i,j,k,n

if (flag .eq. 1) then
k=OutGrid%nz
do j=1,OutGrid%ny
do i=1,OutGrid%nx

OutData(i,j,:) = misval

if ((TopCell(i,j) .gt. 0) .and. (BotCell(i,j) .gt. 0)) then
n=TopCell(i,j)-BotCell(i,j) + 1
OutData(i,j,(k-n+1):k) = PFBdata(i,j,BotCell(i,j):TopCell(i,j))
else
endif
enddo
enddo
else
if ((OutGrid%nx .eq. MstrGrid%nx) .and. &
    (OutGrid%ny .eq. MstrGrid%ny)  .and. &
    (OutGrid%nz .eq. MstrGrid%nz)) then
OutData(:,:,:)=PFBdata(:,:,:)
else
write(*,'(A)') 'ERROR: '
write(*,'(A)') ' You tried to write a compressed grid without compressing it'
write(*,'(A)') ' The attempted and actual grid dimensions follow:'
write(*,*) OutGrid%nx,MstrGrid%nx
write(*,*) OutGrid%ny,MstrGrid%ny
write(*,*) OutGrid%nz,MstrGrid%nz
write(*,'(A)') ' Check your compile options and try again, sorry!'
stop
endif
endif

end subroutine
! --------------------------------------------------------------------------------
subroutine compress_grid(flag)
use global
implicit none
integer,intent(IN) :: flag
integer thick(MstrGrid%nx,MstrGrid%ny)
integer tH,nlys
real zTop

! Version 1: Only use active "layers" but go ahead
!            and keep the areal extent for simplicity

if (flag .eq. 1) then
thick(:,:) = 0

if (allocated(TopCell) .and. allocated(BotCell)) then
! Only do this if Top and Bot are allocated, otherwise there is 
! no reduction in grid size so we don't need them.

! This forces the empty cells to be even more negative
thick(:,:) = TopCell(:,:) - abs(BotCell(:,:))
thick(:,:) = thick(:,:) + 1

tH = maxval(TopCell)
nlys = maxval(thick)
else
! Use normal extents
tH = MstrGrid%nz
nlys = MstrGrid%nz
endif

if (nlys .lt. 1) then
write(*,*) 'ERROR: Fewer than 1 layer'
stop
endif

zTop = MstrGrid%z0 + MstrGrid%dz*real(tH)
OutGrid%z0 = zTop - real(nlys)*MstrGrid%dz

if (nlys .le. MstrGrid%nz) then
OutGrid%nz = nlys
else
OutGrid%nz = MstrGrid%nz
endif
else
OutGrid%z0 = MstrGrid%z0
OutGrid%nz = MstrGrid%nz
endif

! These don't change for now, but they could...
OutGrid%x0 = MstrGrid%x0
OutGrid%y0 = MstrGrid%y0
OutGrid%nx = MstrGrid%nx
OutGrid%ny = MstrGrid%ny
OutGrid%dx = MstrGrid%dx
OutGrid%dy = MstrGrid%dy
OutGrid%dz = MstrGrid%dz

allocate(OutData(OutGrid%nx,OutGrid%ny,OutGrid%nz))

end subroutine
!-------------------------------------------------------------------------
subroutine VTK_SPnts(outfile)
use global
implicit none
character(len=clen),intent(IN) :: outfile
! Writes the data as a structured point dataset

!! This is a binary VTK output
open(Unit=1,file=outfile,action='write',status='replace')
write(1,'(A)')'# vtk DataFile Version 2.0'
write(1,'(A)') 'Requested output'
write(1,'(A)')'BINARY'
write(1,'(A)')'DATASET STRUCTURED_POINTS'
write(1,'(A,3I6)')'DIMENSIONS',OutGrid%nx+1,OutGrid%ny+1,OutGrid%nz+1
write(1,'(A,3F16.4)')'ORIGIN ',OutGrid%x0,OutGrid%y0,OutGrid%z0
write(1,'(A,3F16.4)')'SPACING ',OutGrid%dx,OutGrid%dy,OutGrid%dz

! ----- WRITE BINARY SCALAR DATA -----
open(Unit=1,file=Outfile,action='write',status='old',position='append')
write(1,'(A,I8)')'CELL_DATA ',OutGrid%nx*OutGrid%ny*OutGrid%nz
write(1,'(A)')'SCALARS ',adjustl(trim(Varname)),' float'
write(1,'(A)')'LOOKUP_TABLE  default'
close(1)
open(Unit=1,file=Outfile,action='write',status='old', &
form='unformatted',access='stream',position='append',convert='BIG_ENDIAN')
!write(1) (((real(OutData(i,j,k),kind=4), i=1,OutGrid%nx), j=1,OutGrid%ny), k=1,OutGrid%nz)
write(1) real(OutData(:,:,:))
close(1)

if (num_lyrs) then
if (allocated(Lnum)) then
! ----- WRITE LAYER NUMBERS -----
open(Unit=1,file=Outfile,action='write',status='old',position='append')
write(1,'(A)')'SCALARS LayerNum int'
write(1,'(A)')'LOOKUP_TABLE  default'
close(1)
open(Unit=1,file=Outfile,action='write',status='old', &
form='unformatted',access='stream',position='append',convert='BIG_ENDIAN')
!write(1) (((Lnum(i,j,k), i=1,OutGrid%nx), j=1,OutGrid%ny), k=1,OutGrid%nz)
write(1) Lnum(:,:,:)
endif ! allocated(Lnum)
endif ! (num_lyrs)
close(1)
end subroutine
! --------------------------------------------------------------------------------
subroutine struct_VTK(outfile,npnts,ncell)
use global
implicit none
integer,intent(IN) :: npnts,ncell
character(len=clen),intent(IN) :: outfile
! Writes the data as a structured grid dataset

!! This is a binary VTK output
open(Unit=1,file=outfile,action='write',status='replace')
write(1,'(A)')'# vtk DataFile Version 2.0'
write(1,'(A)') 'Requested output'
write(1,'(A)')'BINARY'
write(1,'(A)')'DATASET STRUCTURED_GRID'
write(1,'(A,3I5)')'DIMENSIONS ',OutGrid%nx+1,OutGrid%ny+1,OutGrid%nz+1
write(1,'(A,I10,A)')'POINTS ',npnts,'  float'
close(1)
open(Unit=1,file=Outfile,action='write',status='old', &
form='unformatted',access='stream',position='append',convert='BIG_ENDIAN')
!write(1) (real(pnts(i,1:3),kind=4), i=1,npnts)  ! This should work if the direct write doesn't
write(1) Pnts(:,:)
close(1)
! ----- SCALAR DATA -----
open(Unit=1,file=Outfile,action='write',status='old',position='append')
write(1,'(A,I10)')'CELL_DATA ',ncell
write(1,'(A)')'SCALARS ',adjustl(trim(Varname)),' float'
write(1,'(A)')'LOOKUP_TABLE  default'
close(1)
open(Unit=1,file=Outfile,action='write',status='old', &
form='unformatted',access='stream',position='append',convert='BIG_ENDIAN')
!write(1) (((real(Outdata(i,j,k),kind=4), i=1,OutGrid%nx), j=1,OutGrid%ny), k=1,OutGrid%nz)
write(1) real(OutData(:,:,:))
close(1)
end subroutine
