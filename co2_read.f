      block data co2_blk

      include 'newmpar.h'
      include 'rdparm.h'  ! needed before other radiation common blocks
      include 'co2dta.h'

c   The following coeffiecients don't depend on resolution or CO2 conc.
      data b0,b1,b2,b3/-.51926410e-4,-.18113332e-3,
     & -.10680132e-5,-.67303519e-7/
      end

c******************************************************************************

      subroutine co2_read(sigma)
c  This routine reads the CO2 transmission coefficients from the
c  co2_datafile (filename set in namelist)
c  was unit 15 for DARLAM, unit 17 for conformal-cubic

      use cc_mpi, only : myid
      include 'newmpar.h'
      include 'rdparm.h'  ! needed before other radiation common blocks
      include 'co2dta.h'
      include 'radisw.h' ! passes rrvco2 to radrive for use in swr89
      parameter(sigtol=1e-3)
      real sigma(kl), sigin(kl)
      data lu/15/
      read(lu,*) nlev
      if (myid==0) print *,'co2_read nlev=',nlev
c     Check that the number of levels is the same
      if ( nlev.ne.kl ) then
	  print*, ' ERROR - Number of levels wrong in co2_data file'
	  stop
      end if
c     Check that the sigma levels are the same
c     Note that the radiation data has the levels in the reverse order
      read(lu,*) (sigin(i),i=kl,1,-1)
      if (myid==0) print *,'co2_read sigin=',sigin
      do k=1,kl
	  if ( abs(sigma(k)-sigin(k)) .gt. sigtol ) then
	      print*, ' ERROR - sigma level wrong in co2_data file'
	      print*, k, sigma(k), sigin(k)
	      stop
          end if
      end do
      read(lu,*) rrvco2
      if (myid==0)
     &     write(*,*) ' CO2 mixing ratio is ', rrvco2*1e6,' ppmv'
      read(lu,*) stemp
      read(lu,*) gtemp
      read(lu,*) cdt51
      read(lu,*) co251
      read(lu,*) c2d51
      read(lu,*) cdt58
      read(lu,*) co258
      read(lu,*) c2d58
      read(lu,*) cdtm51
      read(lu,*) co2m51
      read(lu,*) c2dm51
      read(lu,*) cdtm58
      read(lu,*) co2m58
      read(lu,*) c2dm58
      read(lu,*) cdt31
      read(lu,*) co231
      read(lu,*) c2d31
      read(lu,*) cdt38
      read(lu,*) co238
      read(lu,*) c2d38
      read(lu,*) cdt71
      read(lu,*) co271
      read(lu,*) c2d71
      read(lu,*) cdt78
      read(lu,*) co278
      read(lu,*) c2d78
      read(lu,*) co211
      read(lu,*) co218
      close(lu)
      return
      end
