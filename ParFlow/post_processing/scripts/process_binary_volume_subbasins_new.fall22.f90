Program Indicator_File
implicit none
INTEGER i,j,k,nsub,nx,ny,nz,nx1,ndt,ii,ind1,ind2,ind3,ind4,ind,nwrite,kk
REAL*8 dx, dy,x0,y0,z0,zzi,surf,surf_riv,surf_cv,surf_br,surf_riv_cv,surf_riv_br
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: satur,press,moist,stora,poros,ss
REAL*8, DIMENSION(:,:), ALLOCATABLE :: storage_gw,storage_vz,storage_riv
INTEGER,DIMENSION(:,:), ALLOCATABLE :: sub_ind
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: indi
DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: riv_stora,sat_stora,unsat_stora,dz
Character (LEN=100) Nom_Fich,fichier,Nom_Fich1,fichier1,Nom_Fich2,fichier2,Nom_Fich3,fichier3,Nom_Fich4,fichier4,fid
  character*300 fname01
    nx=150
    nx1=668
    ny=170
    nz=5
    ndt=8759
    nwrite=0 
    dx=100.0
    dy=100.0
    nsub=8
 
    ALLOCATE (indi(nx,ny,nz),satur(nx,ny,nz),press(nx,ny,nz),moist(nx,ny,nz),stora(nx,ny,nz),poros(nx,ny,nz),ss(nx,ny,nz))
    ALLOCATE (storage_gw(nsub,ndt),storage_vz(nsub,ndt),storage_riv(nsub,ndt),sub_ind(nx,ny))
    ALLOCATE (sat_stora(ndt),unsat_stora(ndt),dz(nz),riv_stora(ndt))    
   
    write(*,*) 'ici0' 

    sat_stora=0.0  ;	unsat_stora=0.0		;	riv_stora=0.0
    storage_gw=0.0  ;    storage_vz=0.0  ;    storage_riv=0.0  
     dz(1)=21.0   ;       dz(2)=8.0   ;       dz(3)=0.6   ;       dz(4)=0.3   ;       dz(5)=0.1

    Nom_Fich='bl_wy2015.out.press.000000' 
    ind = INDEX(Nom_Fich,' ')-1
        fichier(1:ind) = nom_fich(1:ind)
        ind = ind + 1

    Nom_Fich1='bl_wy2015.out.satur.000000'
    ind1 = INDEX(Nom_Fich1,' ')-1
    fichier1(1:ind1) = nom_fich1(1:ind1)
    ind1 = ind1 + 1

      fname01='bl_wy2015.out.porosity.pfb'
      call pfb_read(poros,fname01,nx,ny,nz)   

      fname01='bl_wy2015.out.specific_storage.pfb'
      call pfb_read(ss,fname01,nx,ny,nz)

      OPEN(20,FILE='subbasin_indicator.txt') 
      READ(20,*)
      DO j=1,ny
                READ(20,*) sub_ind(1:nx,j)
      END DO

      OPEN(31,FILE='bl_wy2015.groundwater_volume_subbasins.txt') 
      OPEN(32,FILE='bl_wy2015.river_volume_subbasins.txt')
      OPEN(33,FILE='bl_wy2015.vadosezone_volume_subbasins.txt')
         
      DO i=1,ndt

        IF ( mod(i,500) == 0 ) THEN
                write(*,*) i
        END IF

        write(fid, '(I5.5)') i
        fname01='bl_wy2015.out.press.' // trim(adjustl(fid)) // '.pfb'


        !write(*,*) i
        !IF (i < 10) THEN
        !    WRITE(fichier(17:17),'(i1)') i
	!    fname01=fichier(1:17)//'.pfb'
        !ELSE IF (i>=10 .and. i<100) THEN
        !    WRITE(fichier(16:17),'(i2)') i
	!    fname01=fichier(1:17)//'.pfb'
        !ELSE IF (i>=100 .and. i<1000) THEN
        !    WRITE(fichier(15:17),'(i3)') i
        !    fname01=fichier(1:17)//'.pfb'
        !ELSE IF (i>=1000 .and. i<10000) THEN
        !    WRITE(fichier(14:17),'(i4)') i
        !    fname01=fichier(1:17)//'.pfb'   
        !ELSE 
        !    WRITE(fichier(13:17),'(i5)') i
        !    fname01=fichier(1:17)//'.pfb'               
        !END IF
		
	call pfb_read(press,fname01,nx,ny,nz)	


        write(fid, '(I5.5)') i
        fname01='bl_wy2015.out.satur.' // trim(adjustl(fid)) // '.pfb'

        !IF (i < 10) THEN
        !    WRITE(fichier1(17:17),'(i1)') i
        !    fname01=fichier1(1:17)//'.pfb'
        !ELSE IF (i>=10 .and. i<100) THEN
        !    WRITE(fichier1(16:17),'(i2)') i
        !    fname01=fichier1(1:17)//'.pfb'
        !ELSE IF (i>=100 .and. i<1000) THEN
        !    WRITE(fichier1(15:17),'(i3)') i
        !    fname01=fichier1(1:17)//'.pfb'
        !ELSE IF (i>=1000 .and. i<10000) THEN
        !    WRITE(fichier1(14:17),'(i4)') i
        !    fname01=fichier1(1:17)//'.pfb'
        !ELSE
        !    WRITE(fichier1(13:17),'(i5)') i
        !    fname01=fichier1(1:17)//'.pfb'
        !END IF

	call pfb_read(satur,fname01,nx,ny,nz)

!write(*,*) i,'ici02'        
	DO ii=1,nz
            DO j=1,ny
                DO k=1,nx
                        DO kk=1,nsub
                                IF (sub_ind(k,j)==kk) THEN
                                        IF (press(k,j,ii)>0.01) THEN
                                                IF (ii<nz) THEN
                                                        stora(k,j,ii)=satur(k,j,ii)*ss(k,j,ii)*press(k,j,ii)*dx*dy*dz(ii)+(poros(k,j,ii)*satur(k,j,ii)*dx*dy*dz(ii))
                                                      !  stora(k,j,ii)=(poros(k,j,ii)*satur(k,j,ii)*dx*dy*dz(ii))
!                                                        stora(k,j,ii)=satur(k,j,ii)*ss(k,j,ii)*press(k,j,ii)*dx*dy*dz(ii)+(poros(k,j,ii)*satur(k,j,ii)*press(k,j,ii)*dx*dy*dz(ii))
                                                        storage_gw(kk,i)=storage_gw(kk,i)+stora(k,j,ii)
                                                ELSE
                                                        stora(k,j,ii)=press(k,j,ii)*dx*dy
                                                        storage_riv(kk,i)=storage_riv(kk,i)+stora(k,j,ii)
                                                END IF
                                        ELSE
                                                stora(k,j,ii)=poros(k,j,ii)*satur(k,j,ii)*dx*dy*dz(ii)
                                                storage_vz(kk,i)=storage_vz(kk,i)+stora(k,j,ii)
                                        END IF
                                END IF
                        END DO
               END DO
         END DO
    END DO
                WRITE(31,14) i*1.0,storage_gw(1:nsub,i)
                WRITE(32,14) i*1.0,storage_riv(1:nsub,i)
                WRITE(33,14) i*1.0,storage_vz(1:nsub,i)

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



