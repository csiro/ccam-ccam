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


!     has possibility of imax=il*number
!     parameter settings for the longwave and shortwave radiation code: 
!          imax   =  no. points along the lat. circle used in calcs.
!          l      =  no. vertical levels (also layers) in model 
!***note: the user normally will modify only the imax and l parameters
!          nblw   =  no. freq. bands for approx computations. see 
!                      bandta for definition
!          nblx   =  no. freq bands for approx cts computations 
!          nbly   =  no. freq. bands for exact cts computations. see
!                      bdcomb for definition
!          inlte  =  no. levels used for nlte calcs.
!          nnlte  =  index no. of freq. band in nlte calcs. 
!          nb,ko2 are shortwave parameters; other quantities are derived
!                    from the above parameters. 
!   Source file must also include newmpar.h, before rdparm.h

      integer l, imax, nblw, nblx, nbly, nblm, lp1, lp2, lp3
      integer lm1, lm2, lm3, ll, llp1, llp2, llp3, llm1, llm2
      integer llm3, lp1m, lp1m1, lp1v, lp121, ll3p, nb, inlte
      integer inltep, nnlte, lp1i, llp1i, ll3pi, nb1, ko2, ko21, ko2m

      parameter (nblw=163,nblx=47,nbly=15)
      parameter (nblm=nbly-1)
      parameter (nb=12)
      parameter (inlte=3,inltep=inlte+1,nnlte=56) 
      parameter (nb1=nb-1)
      parameter (ko2=12)
      parameter (ko21=ko2+1,ko2m=ko2-1)
      
      common/rdparm/l,imax,lp1,lp2,lp3,lm1,lm2,lm3,ll,llp1,llp2,llp3,    &
     &              llm1,llm2,llm3,lp1m,lp1m1,lp1v,lp121,ll3p,lp1i,      &
     &              llp1i,ll3pi

