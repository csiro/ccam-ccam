      subroutine o3_read(sigma)
c
c     Reads the ozone data from the o3_datafile (filename set in namelist)
c
      use cc_mpi, only : myid ! MJT read

      include 'newmpar.h'
      include 'mpif.h'       ! MJT read
      parameter ( l=kl )

      parameter(sigtol=1.e-3)
      common /o3dat/ dduo3n(37,l),ddo3n2(37,l),ddo3n3(37,l),ddo3n4(37,l)
      real sigma(kl), sigin(kl)

      !--------------------------------------------------------------
      ! MJT read
      if (myid==0) then
        read(16,*) nlev
        if ( nlev.ne.kl ) then
          print*, ' ERROR - Number of levels wrong in o3_data file'
	    stop
        end if
c       Check that the sigma levels are the same
c       Note that the radiation data has the levels in the reverse order
        read(16,*) (sigin(i),i=kl,1,-1)
        do k=1,kl
	    if ( abs(sigma(k)-sigin(k)) .gt. sigtol ) then
	      print*, ' ERROR - sigma level wrong in o3_data file'
	      print*, k, sigma(k), sigin(k)
	      stop
          end if
        end do

c    Note that the data is written as MAM, SON, DJF, JJA. The arrays in
c    o3dat are in the order DJF, MAM, JJA, SON
        read(16,1000) ddo3n2
        read(16,1000) ddo3n4
        read(16,1000) dduo3n
        read(16,1000) ddo3n3
        close(16)
 1000   format(9f8.5)
      end if
      call MPI_Bcast(ddo3n2,37*l,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ddo3n4,37*l,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(dduo3n,37*l,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(ddo3n3,37*l,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      !--------------------------------------------------------------
      
      return
      end

