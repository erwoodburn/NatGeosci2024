Program Indicator_File
implicit none

INTEGER i,j,k,z,iaa,nx,ny,nz,nl,ndt,ii
REAL*8 dx,dy,x0,y0,z0
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: clm,satur
REAL*8, DIMENSION(:,:), ALLOCATABLE :: total_et_june,total_et_july,total_et_aug,total_et_annual,total_n_et_annual
REAL*8, DIMENSION(:,:), ALLOCATABLE :: total_infl_annual,peak_swe_value
INTEGER, DIMENSION(:,:), ALLOCATABLE :: peak_swe_hour
!REAL*8, DIMENSION(:,:), ALLOCATABLE :: total_satur_june,total_satur_july,total_satur_aug,total_satur_annual

Character (LEN=200) fid
  character*200 fname01,fname02
    nx=150
    ny=170
    !nz is number of CLM layers here
    nz=17
    nl=5
    ndt=17500
    !ndt=17496
    dx=100.0
    dy=100.0

    ALLOCATE (clm(nx,ny,nz),satur(nx,ny,5))
    ALLOCATE (total_et_june(nx,ny),total_et_july(nx,ny),total_et_aug(nx,ny),total_et_annual(nx,ny),total_n_et_annual(nx,ny))
    ALLOCATE (total_infl_annual(nx,ny),peak_swe_hour(nx,ny),peak_swe_value(nx,ny))
    !ALLOCATE (total_satur_june(nx,ny),total_satur_july(nx,ny),total_satur_aug(nx,ny),total_satur_annual(nx,ny))

    peak_swe_hour=0
    peak_swe_value=0.

      OPEN(29,FILE='et_june_map.bl_wy2015.txt')
      OPEN(30,FILE='et_july_map.bl_wy2015.txt')
      OPEN(31,FILE='et_aug_map.bl_wy2015.txt')
      OPEN(32,FILE='et_annual_map.bl_wy2015.txt')
      OPEN(932,FILE='n_et_annual_map.bl_wy2015.txt')
      OPEN(33,FILE='infl_annual_map.bl_wy2015.txt')
      OPEN(34,FILE='peak_swe_hour_map.bl_wy2015.txt')
      OPEN(35,FILE='peak_swe_value_map.bl_wy2015.txt')

      !OPEN(299,FILE='satur_june_map.bl_wy2015.txt')
      !OPEN(300,FILE='satur_july_map.bl_wy2015.txt')
      !OPEN(311,FILE='satur_aug_map.bl_wy2015.txt')      
      !OPEN(322,FILE='satur_annual_map.bl_wy2015.txt')

      !Main time loop here - 1
      DO i=5833,6552,1
        
      IF (i==6552) THEN
        write(*,*) 'End June'
      END IF

      !IF ( mod(i,500) == 0 ) THEN
	!        write(*,*) i
       ! END IF

        write(fid, '(I5.5)') i
        fname01='bl_wy2015.out.clm_output.' // trim(adjustl(fid)) // '.C.pfb'
        !fname02='bl_wy2015.out.satur.' // trim(adjustl(fid)) // '.pfb'
        !write(*,*) fname01
        !write(*,*) fname02

	call pfb_read(clm,fname01,nx,ny,nz)	
        !call pfb_read(satur,fname02,nx,ny,nl)

        !CLM array structure:
                !evap_tot,5]
                !evap_grnd,6]
                !evap_soi,7]
                !evap_veg,8]
                !tran_veg,9]
                !infl,10]
                !swe,11]

            DO j=1,ny
                DO k=1,nx
                        total_et_june(k,j)=total_et_june(k,j)+(clm(k,j,5)*3600.)
                        !5 above is total et of .C.pfb but 5 below is top layer of 5
                        !total_satur_june(k,j)=total_satur_june(k,j)+satur(k,j,5)
                END DO
            END DO
 END DO
 
             DO j=1,ny
                DO k=1,nx
                        WRITE(29,14) total_et_june(k,j)
                        !WRITE(299,14) total_satur_june(k,j)/719.                        
                END DO
             END DO


      !Main time loop here - 2
      DO i=6553,7296,1

      IF (i==7296) THEN
        write(*,*) 'End July'
      END IF

        write(fid, '(I5.5)') i
        fname01='bl_wy2015.out.clm_output.' // trim(adjustl(fid)) // '.C.pfb'
        !fname02='bl_wy2015.out.satur.' // trim(adjustl(fid)) // '.pfb'
        !write(*,*) fname01

        call pfb_read(clm,fname01,nx,ny,nz)
        !call pfb_read(satur,fname02,nx,ny,nl)

            DO j=1,ny
                DO k=1,nx
                        total_et_july(k,j)=total_et_july(k,j)+(clm(k,j,5)*3600.)
                        !total_satur_july(k,j)=total_satur_july(k,j)+satur(k,j,5)
                END DO
            END DO
 END DO

             DO j=1,ny
                DO k=1,nx
                        WRITE(30,14) total_et_july(k,j)
                        !WRITE(300,14) total_satur_july(k,j)/743.                        
                END DO
             END DO

             

      !Main time loop here - 3
      DO i=7297,8040,1

      IF (i==8040) THEN
        write(*,*) 'End Aug'
      END IF

        write(fid, '(I5.5)') i
        fname01='bl_wy2015.out.clm_output.' // trim(adjustl(fid)) // '.C.pfb'
        !fname02='bl_wy2015.out.satur.' // trim(adjustl(fid)) // '.pfb'
        !write(*,*) fname01

        call pfb_read(clm,fname01,nx,ny,nz)
        !call pfb_read(satur,fname02,nx,ny,nl)

            DO j=1,ny
                DO k=1,nx
                        total_et_aug(k,j)=total_et_aug(k,j)+(clm(k,j,5)*3600.)
                        !total_satur_aug(k,j)=total_satur_aug(k,j)+satur(k,j,5)
                END DO
            END DO
 END DO

             DO j=1,ny
                DO k=1,nx
                        WRITE(31,14) total_et_aug(k,j)
                        !WRITE(311,14) total_satur_aug(k,j)/743.
                END DO
             END DO

      !Main time loop here - 4
      !Also print a file with cumulative negative ET
      !Also do infl annual and calc peak swe hour and value
      DO i=1,8759,1

      IF (i==8759) THEN
              write(*,*) 'End Annual'
      END IF

        write(fid, '(I5.5)') i
        fname01='bl_wy2015.out.clm_output.' // trim(adjustl(fid)) // '.C.pfb'
        !fname02='bl_wy2015.out.satur.' // trim(adjustl(fid)) // '.pfb'
        !write(*,*) fname01

        call pfb_read(clm,fname01,nx,ny,nz)
        !call pfb_read(satur,fname02,nx,ny,nl)

            DO j=1,ny
                DO k=1,nx
                        total_et_annual(k,j)=total_et_annual(k,j)+(clm(k,j,5)*3600.)

                        IF (clm(k,j,5)<0.0) THEN
                                total_n_et_annual(k,j)=total_n_et_annual(k,j)+(clm(k,j,5)*3600.)
                        END IF

                        !total_satur_annual(k,j)=total_satur_annual(k,j)+satur(k,j,5)
                        total_infl_annual(k,j)=total_infl_annual(k,j)+(clm(k,j,10))

                        !Check for peak SWE hour for each cell        
                        IF ( (clm(k,j,11) > peak_swe_value(k,j)) ) THEN
                                peak_swe_value(k,j)=clm(k,j,11)
                                peak_swe_hour(k,j)=i
                        END IF
                END DO
            END DO
 END DO

             DO j=1,ny
                DO k=1,nx
                        WRITE(32,14) total_et_annual(k,j)
                        WRITE(932,14) total_n_et_annual(k,j)
                        !WRITE(322,14) total_satur_annual(k,j)/8759.
                        WRITE(33,14) total_infl_annual(k,j)
                        WRITE(34,15) peak_swe_hour(k,j)
                        WRITE(35,14) peak_swe_value(k,j)
                END DO
             END DO

 14 FORMAT(2000(1x,e13.7))
 15 FORMAT(2000(I5))

END
 
 subroutine pfb_read(value,fname,nx,ny,nz) 
  implicit none
  real*8 value(nx,ny,nz)
  real*8 dx, dy, dz, x1, y1, z1								
  integer*4 i,j,k, nni, nnj, nnk, ix, iy, iz,			&
            ns,  rx, ry, rz,nx,ny,nz, nnx, nny, nnz,is
  character*100 fname
  
  open(100,file=trim(adjustl(fname)),form='unformatted',   &
                 recordtype='stream',convert='BIG_ENDIAN',status='old') 
  ! Read in header info

! Start: reading of domain spatial information
  read(100) x1 !X
  read(100) y1 !Y 
  read(100) z1 !Z

  read(100) nx !NX
  read(100) ny !NY
  read(100) nz !NZ

  read(100) dx !DX
  read(100) dy !DY
  read(100) dz !DZ

  read(100) ns !num_subgrids
! End: reading of domain spatial information

! Start: loop over number of sub grids
  do is = 0, (ns-1)

! Start: reading of sub-grid spatial information
   read(100) ix
   read(100) iy
   read(100) iz

   read(100) nnx
   read(100) nny
   read(100) nnz
   read(100) rx
   read(100) ry
   read(100) rz

! End: reading of sub-grid spatial information

! Start: read in saturation data from each individual subgrid
  do  k=iz +1 , iz + nnz
   do  j=iy +1 , iy + nny
    do  i=ix +1 , ix + nnx
     read(100) value(i,j,k)
    end do
   end do
  end do
! End: read in saturation data from each individual subgrid

! End: read in saturation data from each individual subgrid

  end do
! End: loop over number of sub grids



  close(100)
  end subroutine


