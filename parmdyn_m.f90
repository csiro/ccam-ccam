! Conformal Cubic Atmospheric Model
    
! Copyright 2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module parmdyn_m

implicit none

private
public mex,mfix,mfix_qg,mfix_tr,mfix_aero,mfix_t
public nh,nritch_t,mspec,mup
public nstag,nstagu,ntbar,precon,helmmeth
public nstagoff
public epsp,epsu,epsf,epsh,restol

!     dynamics options

!     parameter (mfix_qg=1)   ! 1 "mass" fix for qg
!                               2 "mass" fix for qg and trace gases

!     mfix in namelist:        -1 on pslx in upglobal
!                               0 off
!                               1 cunning in adjust5
!                               2 more-cunning in adjust5

!            (ntbar=0)           ! 0 for standard
!            (ntbar=(kl+1)/2)    ! level# for tbar2d with T set in nonlin

integer, save :: mex=30,mfix=3,mfix_qg=1,mfix_tr=0,mfix_aero=0,mfix_t=0
integer, save :: nh=0,nritch_t=300,mspec,mup=1
integer, save :: nstag=-10,nstagu=-1,ntbar,precon=-2900,helmmeth=1
integer, save :: nstagoff=0
real, save :: epsp=-15.,epsu=0.,epsf=0.,epsh=1.,restol=4.e-7
    
end module parmdyn_m