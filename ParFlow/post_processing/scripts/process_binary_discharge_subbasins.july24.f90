Program Indicator_File
implicit none
INTEGER i,j,k,iaa,nsub,nx,ny,nz,nx1,nz1,ndt,ii,ind1,ind2,ind3,ind4,ind,nwrite,kk
REAL*8 dx, dy,x0,y0,z0
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: press,slopex,slopey,slopet,manning
REAL*8, DIMENSION(:,:), ALLOCATABLE :: discharge_sub,press_sub
INTEGER,DIMENSION(:,:), ALLOCATABLE :: sub_ind
INTEGER,DIMENSION(:), ALLOCATABLE :: out_k,out_j
Character (LEN=200) Nom_Fich,fichier,Nom_Fich1,fichier1,Nom_Fich2,fichier2,Nom_Fich3,fichier3,Nom_Fich4,fichier4,fid
  character*200 fname01
    nx=150
    nx1=668
    ny=170
    nz=5
    !ndt=35100

    !TEMPORARY!
    ndt=8759
    !ndt=6369
    !ndt=8783
    nwrite=0 
    dx=100.0
    dy=100.0
    nsub=8
    nz1=1
 
    ALLOCATE (press(nx,ny,nz),slopex(nx,ny,nz1),slopey(nx,ny,nz1),slopet(nx,ny,nz1),manning(nx,ny,nz1))
    ALLOCATE (discharge_sub(nsub,ndt),out_k(nsub),out_j(nsub),sub_ind(nx,ny),press_sub(nsub,ndt))
   
    write(*,*) 'ici0' 

    discharge_sub=0.0

    Nom_Fich='bl_wy2015.out.press.00000' 
    ind = INDEX(Nom_Fich,' ')-1
        fichier(1:ind) = nom_fich(1:ind)
        ind = ind + 1

    Nom_Fich1='bl_wy2015.out.satur.00000'
    ind1 = INDEX(Nom_Fich1,' ')-1
    fichier1(1:ind1) = nom_fich1(1:ind1)
    ind1 = ind1 + 1

      fname01='ER_slope_x.pfb'
      call pfb_read(slopex,fname01,nx,ny,nz1)   

      fname01='ER_slope_y.pfb'
      call pfb_read(slopey,fname01,nx,ny,nz1)

      fname01='ER_mannings_100m.pfb'
      call pfb_read(manning,fname01,nx,ny,nz1)

      DO j=1,ny
        DO k=1,nx
                IF (slopex(k,j,nz1)==0.0) THEN
                        slopet(k,j,nz1)=slopey(k,j,nz1)
                ELSE IF (slopey(k,j,nz1)==0.0) THEN
                        slopet(k,j,nz1)=slopex(k,j,nz1)
                ELSE
                        WRITE(*,*) k,j
                END IF
        END DO
      END DO

      OPEN(20,FILE='subbasin_indicator.txt') 
      READ(20,*)
      DO j=1,ny
                READ(20,*) sub_ind(1:nx,j)
      END DO
      CLOSE (20)

      OPEN(20,FILE='subbasin_outlet_new.txt') 
      DO i=1,nsub
        READ(20,*) iaa,out_k(i),out_j(i)
      END DO
      CLOSE (20)

      OPEN(31,FILE='bl_wy2015.discharge_subbasin1-8.corr_july24.txt') 
      OPEN(32,FILE='bl_wy2015.discharge_check.txt')
      OPEN(33,FILE='bl_wy2015.pressure_check.txt')

      DO i=1,ndt

        IF ( mod(i,500) == 0 ) THEN
                write(*,*) i
        END IF

        write(fid, '(I5.5)') i
        fname01='bl_wy2015.out.press.' // trim(adjustl(fid)) // '.pfb'

        !write(*,*) i
        !IF (i < 10) THEN
        !    !WRITE(fichier(17:17),'(i1)') i
        !    !fname01=fichier(1:17)//'.pfb'
        !    WRITE(fichier(22:22),'(i1)') i
        !    fname01=fichier(1:22)//'.pfb'
        !ELSE IF (i>=10 .and. i<100) THEN
        !    !WRITE(fichier(16:17),'(i2)') i
	!    !fname01=fichier(1:17)//'.pfb'
        !    WRITE(fichier(21:22),'(i2)') i
        !    fname01=fichier(1:22)//'.pfb'
        !ELSE IF (i>=100 .and. i<1000) THEN
        !   ! WRITE(fichier(15:17),'(i3)') i
        !   ! fname01=fichier(1:17)//'.pfb'
        !    WRITE(fichier(20:22),'(i3)') i
        !    fname01=fichier(1:22)//'.pfb'           
        !ELSE IF (i>=1000 .and. i<10000) THEN
        !    !WRITE(fichier(14:17),'(i4)') i
        !    !fname01=fichier(1:17)//'.pfb'   
        !    WRITE(fichier(19:22),'(i4)') i
        !    fname01=fichier(1:22)//'.pfb'             
        !ELSE 
        !    !WRITE(fichier(13:17),'(i5)') i
        !    !fname01=fichier(1:17)//'.pfb'               
        !    WRITE(fichier(18:22),'(i5)') i
        !    fname01=fichier(1:22)//'.pfb'             
        !END IF
		
	call pfb_read(press,fname01,nx,ny,nz)	

	
        DO kk=1,nsub
                !discharge_sub(kk,i)=(press(out_k(kk),out_j(kk),nz)**(5./3.))*(sqrt(abs(slopet(out_k(kk),out_j(kk),nz1))))/(manning(out_k(kk),out_j(kk),nz1)*100/3600.)
                discharge_sub(kk,i)=(press(out_k(kk),out_j(kk),nz)**(5./3.))*(sqrt(abs(slopet(out_k(kk),out_j(kk),nz1))))/(manning(out_k(kk),out_j(kk),nz1))*100.
                press_sub(kk,i)=press(out_k(kk),out_j(kk),nz)
        END DO
                WRITE(31,14) i*1.0,discharge_sub(1:nsub,i)
                WRITE(32,14) i*1.0,press(out_k(1),out_j(1),nz),slopet(out_k(1),out_j(1),nz1),manning(out_k(1),out_j(1),nz1)
                WRITE(33,14) i*1.0,press_sub(1:nsub,i)

    END DO
 14 FORMAT(2000(1x,e13.7))

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


