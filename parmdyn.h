! Conformal Cubic Atmospheric Model
    
! Copyright 2015 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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


!     dynamics options (globpe, adjust5, nonlin, upglobal)

!     parameter (mfix_qg=1)   ! 1 "mass" fix for qg
!                               2 "mass" fix for qg and trace gases

!     mfix in namelist:        -1 on pslx in upglobal
!                               0 off
!                               1 cunning in adjust5
!                               2 more-cunning in adjust5

      integer         mex,mfix,mfix_qg,mspec,mup,mfix_tr
      integer         nh,nritch_t,mfix_aero
      integer         nstag,nstagu,ntbar,precon,helmmeth
      integer         nstagoff
      real            epsp,epsu,epsf,epsh,restol
      
      common/paramdyn/mex,mfix,mfix_qg,mspec,mup,                        &
     &                nh,nritch_t,                                       &
     &                nstag,nstagu,nstagoff,ntbar,precon,                &
     &                helmmeth,epsp,epsu,epsf,epsh,restol,mfix_tr,       &
     &                mfix_aero

!            (ntbar=0)           ! 0 for standard
!            (ntbar=(kl+1)/2)    ! level# for tbar2d with T set in nonlin

!            nvsplit    0  uses tendencies for vadv, radn & vertmix
!                       1  splits radn, vertmix, gwdrag, conjob (not vadv)
!                       2  splits radn, vertmix, gwdrag, conjob & vadv
!                       3  splits just vadv
!                      -1  splits just vertmix 
!                      N.B. qg always split for vadv
!                      N.B. always split for vadv called from adjust5
