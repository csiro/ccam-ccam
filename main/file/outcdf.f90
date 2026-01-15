! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2026 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
! This file is part of the Conformal Cubic Atmospheric Model (CCAM)
!
! CCAM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! CCAM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with CCAM.  If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------

! CCAM NetCDF output routines

! itype=1,  iout=20  write outfile history file
! itype=-1, iout=19  write restart file (uncompressed)
! itype=-1, iout=21  write ensemble file (uncompressed)
! hp_output=0        compressed history file
! hp_output=1        uncompressed history file
! localhist=f        single file output 
! localhist=t        parallel output for a group of processors (e.g., for a node)

! Three output options exist for standard (all variables), cordex (all cordex
! variables) and freq (limited variables).  Standard output is typically
! six-hourly, cordex output is hourly and freq is intended for 10min output.
    
! Thanks to Paul Ryan for optimising netcdf routines.
    
module outcdf

use outcdf_common_m
use outcdf_cordex_m
use outcdf_freq_m
use outcdf_standard_m

private
public outfile, freqfile_cordex, freqfile_10, mslp

contains
    
subroutine outfile(iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
      
use cc_mpi         ! CC MPI routines
use dates_m        ! Date data
use parm_m         ! Model configuration
      
implicit none

integer, intent(in) :: iout
real, dimension(:), intent(in) :: psl_in
real, dimension(:,:), intent(in) :: u_in, v_in, t_in, q_in
character(len=*), intent(in) :: cdffile_in

call START_LOG(outfile_begin)
      
! Older text file for soil
if ( nrungcm==-2 .or. nrungcm==-3 .or. nrungcm==-5 ) then
  call soiltextfile  
endif      ! (nrungcm.eq.-2.or.nrungcm.eq.-3.or.nrungcm.eq.-5)

!---------------------------------------------------------------------------
! NetCDF files for history and restart
select case(iout)
  case(19)  
    select case(io_rest)  
      case(0)  ! No output
      case(1)  ! NetCDF 
        if ( myid==0 ) then
          write(6,*) "Restart write of data to netCDF"
        end if  
        call cdfout(-1,iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
      case default
        if ( myid==0 ) then
          write(6,*) "ERROR: unsupported file format io_rest ",io_rest
          write(6,*) "       valid options are 0=none, 1=NetCDF"
          call ccmpi_abort(-1)
        end if
    end select
  case(20)     
    select case(io_out)
      case(0)  ! No output
      case(1)  ! NetCDF
        call cdfout(1,iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
      case default
        if ( myid==0 ) then
          write(6,*) "ERROR: unsupported file format io_out ",io_out
          write(6,*) "       valid options are 0=none, 1=NetCDF"
          call ccmpi_abort(-1)
        end if
    end select
  case(21)    
     if ( myid==0 ) then
       write(6,*) "Ensemble write of data to netCDF"
     end if
    call cdfout(-1,iout,cdffile_in,psl_in,u_in,v_in,t_in,q_in)
  case default  
    if ( myid==0 ) then
      write(6,*) "ERROR: Unknown output file option iout=",iout
      call ccmpi_abort(-1)
    end if  
end select

call END_LOG(outfile_end)
      
return
end subroutine outfile    

subroutine soiltextfile

use arrays_m       ! Atmosphere dyamics prognostic arrays
use cc_mpi         ! CC MPI routines
use dates_m        ! Date data
use filnames_m     ! Filenames
use parm_m         ! Model configuration
use pbl_m          ! Boundary layer arrays
use soilsnow_m     ! Soil, snow and surface data

implicit none

character(len=1024) :: surfout
character(len=20) :: qgout

if ( ktau==nwrite/2 .or. ktau==nwrite ) then
!        usually after first 24 hours, save soil variables for next run
  if ( ktau==nwrite ) then  ! 24 hour write
    if ( ktime==1200 ) then
      surfout=surf_12   ! 'current.1200'
      qgout='qg_12'
    else
      surfout=surf_00   ! 'current.0000'
      qgout='qg_00'
    endif
  else                  ! 12 hour write
    if(ktime==1200)then
      surfout=surf_00   ! 'current.0000'
      qgout='qg_00'
    else
      surfout=surf_12   ! 'current.1200'
      qgout='qg_12'
    end if
  end if                ! (ktau.eq.nwrite)
  if ( myid==0 ) then
    write(6,*) "writing current soil & snow variables to ",surfout
    open(unit=77,file=surfout,form='formatted',status='unknown')
    write (77,*) kdate,ktime,' ktau = ',ktau
  end if
  call writeglobvar(77, wb, fmt='(14f6.3)')
  call writeglobvar(77, tgg, fmt='(12f7.2)')
  call writeglobvar(77, tss, fmt='(12f7.2)')
  call writeglobvar(77, snowd, fmt='(12f7.1)')
  call writeglobvar(77, sicedep, fmt='(12f7.1)')
  if ( myid==0 ) then
    close(77)
  end if
  if ( nrungcm==-2 .or. nrungcm==-5 ) then
    if ( myid==0 ) then
      write(6,*) "writing special qgout file: ",qgout
      open(unit=77,file=qgout,form='unformatted',status='unknown')
    end if
    call writeglobvar(77, qg)
    if ( myid==0 ) then
      close(77)
    end if
  endif  ! (nrungcm.eq.-2.or.nrungcm.eq.-5)
endif    ! (ktau.eq.nwrite/2.or.ktau.eq.nwrite)

return
end subroutine soiltextfile

end module outcdf
