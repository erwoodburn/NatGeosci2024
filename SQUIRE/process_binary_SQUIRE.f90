Program Indicator_File
implicit none
INTEGER i,j,k,nx,ny,nz,nx1,ndt,ii,ind1,nwrite,kk,d,a
REAL*8 dx, dy,x0,y0,z0,zzi,T_frz,snowfall_incr,rainfall_incr
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: T_air,precip
REAL*8, DIMENSION(:,:), ALLOCATABLE :: snowfall_sum_all,rainfall_sum_all

Character (LEN=200) fid,fid2
  character*200 fname01,fname02
    nx=150
    ny=170
    nz=1
    !ndt=8759
    nwrite=0 
    dx=100.0
    dy=100.0
 
    ALLOCATE (snowfall_sum_all(nx,ny))
    ALLOCATE (rainfall_sum_all(nx,ny))

    ALLOCATE (T_air(nx,ny,nz))
    ALLOCATE (precip(nx,ny,nz))


!Clear arrays for maps
snowfall_sum_all=0.
rainfall_sum_all=0.

 T_frz=273.66

! Main time loop here
! Sum only during SAIL SQUIRE period
! Dec 2021-March 2022
! Dec 1 = hour 1465
! March 31 = 4368

    !DO i=1465,4368
     DO i=1465,2568
        write(fid, '(I6.6)') i
        fname01='NLDAS.APCP.' // trim(adjustl(fid)) // '.pfb'
	call pfb_read(precip,fname01,nx,ny,nz)

        fname02='NLDAS.Temp.' // trim(adjustl(fid)) // '.pfb'
        call pfb_read(T_air,fname02,nx,ny,nz)

      DO j=1,ny
       DO k=1,nx

         rainfall_incr=0.
         snowfall_incr=0.
         
         IF (T_air(k,j,1)<T_frz) THEN
           snowfall_sum_all(k,j)=snowfall_sum_all(k,j)+precip(k,j,1) 
        ENDIF
        IF ((T_air(k,j,1)>=T_frz).AND.(T_air(k,j,1)<=(T_frz+2))) THEN
           rainfall_incr=(-54.632+(0.2*T_air(k,j,1)))*precip(k,j,1)
           rainfall_sum_all(k,j)=rainfall_sum_all(k,j)+rainfall_incr
           snowfall_incr=precip(k,j,1)-rainfall_incr
           snowfall_sum_all(k,j)=snowfall_sum_all(k,j)+snowfall_incr
        ENDIF
        IF (T_air(k,j,1)>T_frz+2) THEN
           rainfall_sum_all(k,j)=rainfall_sum_all(k,j)+precip(k,j,1)
        ENDIF


      END DO
     END DO

   END DO !Dec 1, 2021 - Jan 15, 2022 (midnight)

  ! 
  !4 day data gap in squire Jan 16 - Jan 19, 2022
  !

  !Continue onward with rest of winter 
  !Jan 20 - Mar 31, 2022
      DO i=2665,4368
        write(fid, '(I6.6)') i
        fname01='NLDAS.APCP.' // trim(adjustl(fid)) // '.pfb'
        call pfb_read(precip,fname01,nx,ny,nz)

        fname02='NLDAS.Temp.' // trim(adjustl(fid)) // '.pfb'
        call pfb_read(T_air,fname02,nx,ny,nz)

      DO j=1,ny
       DO k=1,nx

         rainfall_incr=0.
         snowfall_incr=0.

         IF (T_air(k,j,1)<T_frz) THEN
           snowfall_sum_all(k,j)=snowfall_sum_all(k,j)+precip(k,j,1)
        ENDIF
        IF ((T_air(k,j,1)>=T_frz).AND.(T_air(k,j,1)<=(T_frz+2))) THEN
           rainfall_incr=(-54.632+(0.2*T_air(k,j,1)))*precip(k,j,1)
           rainfall_sum_all(k,j)=rainfall_sum_all(k,j)+rainfall_incr
           snowfall_incr=precip(k,j,1)-rainfall_incr
           snowfall_sum_all(k,j)=snowfall_sum_all(k,j)+snowfall_incr
        ENDIF
        IF (T_air(k,j,1)>T_frz+2) THEN
           rainfall_sum_all(k,j)=rainfall_sum_all(k,j)+precip(k,j,1)
        ENDIF


      END DO
     END DO

   END DO !Dec 1, 2021 - Jan 15, 2022 (midnight)

        




      !!!!! CHANGE WY NAME HERE!

      !Print annual maps
      OPEN(95,FILE="SQUIRE_PERIOD.NLDAS_rainfall_sum_map.omitJan15-19.h.txt")
      OPEN(96,FILE="SQUIRE_PERIOD.NLDAS_snowfall_sum_map.omitJan15-19.h.txt")

      WRITE(95,*) "150 170 1 "
      WRITE(96,*) "150 170 1 "

      DO j=1,ny
                DO k=1,nx
                        WRITE(95,*) rainfall_sum_all(k,j)
                        WRITE(96,*) snowfall_sum_all(k,j)
                END DO
      END DO

 14 FORMAT(2000(1x,e13.7))
 15 FORMAT(2000(I4,',',A,',',F19.18))

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


