      subroutine o3read_amip

!  Read the AMIP2 lat-height ozone data (volume mixing ratio)
!  Based on the supplied example program to read the data.
!  This routine is fixed format f90

      use cc_mpi, only : myid ! MJT read
      use o3amip_m

      implicit none
      include 'mpif.h'       ! MJT read
      include 'filnames.h'   ! MJT read

      character(len=80)  :: finame
      character(len=120) :: label
      real, parameter :: amd = 28.9644, amo = 48.0000
      real, parameter :: massratio = amo / amd
      integer :: lato3d  ! = jg  number of data latitudinal grid
      integer :: layo3d  ! = kg  number of data vertical layers
      integer :: lvlo3d  ! = lg  number of data vertical layer interfaces
      integer :: j, k, month, ierr
      integer, parameter :: un = 16  ! Unit number for file
      real, dimension(kg) :: galt, gprs
      real, dimension(lg) :: gali
      real :: o3vubc  ! =     upper boundary o3 vmr (ppmv)
      real :: prsubc  ! =     upper boundary pressure (mb)
      real :: altubc  ! =     upper boundary altitude (km)


      !--------------------------------------------------------------
      ! MJT read
      if (myid==0) then

!       Opened in globpe.f
!        finame = "amip2o3.dat"
!        open (un,file=finame,status="old",action="read",iostat=ierr)
!        if ( ierr /= 0 ) then
!           print*, "Error opening file amip2o3.dat"
!           stop
!        end if
        open(un,file=o3file,form='formatted',status='old') ! MJT read

        read (un,"(i3)") lato3d
        read (un,"(i3)") layo3d
        lvlo3d  =    1 + layo3d
        if ( lato3d /= jg ) then
           print*, "Error in horizontal resolution of ozone data"
           stop
        end if
        if ( layo3d /= kg ) then
           print*, "Error in vertical resolution of ozone data"
           stop
        end if
!       print*, lato3d, layo3d

        read (un,"(3(1pe12.5))") o3vubc, prsubc, altubc
      
        read (un,"(a)") label   ! grid latitudes (deg) --------->
        read (un,"(10(1pe12.5))") (glat(j),j=1,lato3d)
      
        read (un,"(a)") label   ! layer pressure (mb) ---------->
        read (un,"(10(1pe12.5))") (gprs(k),k=1,layo3d)

        read (un,"(a)") label   ! layer altitude (km) ---------->
        read (un,"(10(1pe12.5))") (galt(k),k=1,layo3d)

        read (un,"(a)") label   ! interface pressure (mb) ------>
        read (un,"(10(1pe12.5))") (gpri(k),k=1,lvlo3d)

        read (un,"(a)") label   ! interface altitude (km) ------>
        read (un,"(10(1pe12.5))") (gali(k),k=1,lvlo3d)

        do month = 1,mo
           read (un,"(a)") label   ! monthly o3 vmr (ppmv) lat-alt
           read (un,"(10(1pe12.5))")
     $        ((gdat(j,k,month),j=1,lato3d),k=1,layo3d)
        enddo
        close(un)
      end if
      call MPI_Bcast(o3vubc,1,MPI_REAL,0,MPI_COMM_WORLD,ierr) ! not needed
      call MPI_Bcast(prsubc,1,MPI_REAL,0,MPI_COMM_WORLD,ierr) ! not needed
      call MPI_Bcast(altubc,1,MPI_REAL,0,MPI_COMM_WORLD,ierr) ! not needed
      call MPI_Bcast(glat,jg,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(gprs,kg,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(galt,kg,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(gpri,lg,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(gali,lg,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(gdat,jg*kg*mo,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !--------------------------------------------------------------

      dp = gpri(2:lvlo3d) - gpri(1:lvlo3d-1)

!     Convert from ppmv to mass mixing ratio
      gdat = gdat * 1e-6 * massratio

      end subroutine o3read_amip
