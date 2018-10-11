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
    
module soilv_m

use newmpar_m, only : mxst, mxvt, ms

implicit none

private
public i2bp3, ibp2, cnsd, hsbh
public swilt, ssat, sfc
public bch, css, hyds, rhos, sucs, clay, sand, silt
public rlaim44, rlais44, scveg44, rsunc44, slveg44
public froot, zse
public zshh, ww

integer, dimension(mxst), save :: i2bp3, ibp2
real, dimension(mxst), save :: cnsd, hsbh
real, dimension(0:mxst), save :: swilt = (/ 0., .072, .216, .286, .135, .219, .283, .175, .395, .216, .1142, .1547, .2864, &
                                                .2498/)
real, dimension(0:mxst), save :: ssat = (/ 2., .398, .479, .482, .443, .426, .482, .420, .451, .479, .435, .451, .482, .476/)
real, dimension(0:mxst), save :: sfc = (/ 1.,  .143, .301, .367, .218, .31 , .37 , .255, .45, .301, .22 , .25 , .367, .294/)
! bch for gravity term
real, dimension(mxst), save :: bch = (/ 4.2, 7.1, 11.4, 5.15, 10.4, 10.4, 7.12, 5.83, 7.1, 4.9, 5.39, 11.4, 8.52/)    
! heat capacity
real, dimension(mxst), save :: css = (/ 850., 850., 850., 850., 850., 850., 850., 1920., 2100., 850., 850., 850., 850./) 
real, dimension(mxst), save :: hyds = (/ 166.e-6, 4.e-6, 1.e-6, 21.e-6, 2.e-6, 1.e-6, 6.e-6,800.e-6, 1.e-6, 34.e-6, 7.e-6, &
                                              1.3e-6, 2.5e-6/)
real, dimension(mxst), save :: rhos = (/ 2600., 2600., 2600., 2600., 2600., 2600., 2600., 1300.,  910., 2600., 2600., 2600., &
                                              2600./)     ! soil density
real, dimension(mxst), save :: sucs = (/ -.106, -.591, -.405, -.348, -.153, -.49, -.299,-.356, -.153, -.218, -.478, -.405, &
                                              -.63/) ! phisat (m)
real, dimension(mxst), save :: clay = (/ .09, .3, .67, .2, .42, .48, .27, .17, .30, .2, .3, .3, .67/)    ! with mxst=13
real, dimension(mxst), save :: sand = (/ .83, .37, .16, .6, .52, .27, .58, .13, .37, .6, .37, .37, .17/) ! with mxst=13
real, dimension(mxst), save :: silt = (/ .08, .33, .17, .2, .06, .25, .15, .70, .33, .2, .33, .33, .17/) ! with mxst=13
real, dimension(44), parameter :: rlaim44 = (/ 4.8, 6.3, 5., 3.75, 2.78, 2.5, 3.9, 2.77, 2.04, 2.6,         & ! 1-10
                                               1.69, 1.9, 1.37, 1.5, 1.21, 1.58, 1.41, 2.3, 1.2, 1.71,      & ! 11-20
                                               1.21, 2.3, 2.3, 1.2, 1.2, 1.87, 1., 3., .01, .01, 1.2,       & ! 21-31
                                               6., 5.5, 5., 4.5, 5., 4., 3., 3.5, 1., 4., .5, 4., 0./)        ! 32-44
real, dimension(44), parameter :: rlais44 = (/ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,                      & ! 1-10
                                               1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,                      & ! 11-20
                                               1., 1., 1., 1., .6, .6, .5, 1., 0., 0., 1.,                  & ! 21-31
                                               2., 2., 2., 2., 2., 1.5, 1.5, 1.5, 1., .5, .5, .5, 0./)        ! 32-44
real, dimension(44), parameter :: scveg44 = (/ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                      & ! 1-10
                                               0., 0., 0., 0., 0., .1, .1, .1, .1, .1,                      & ! 11-20
                                               .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,                  & ! 21-31
                                               .05, 0., 0., 0., 0., .05, .05, .05, .1, 0., 0., .4, 0./)       ! 32-44
real, dimension(44), parameter :: rsunc44 = (/ 370., 330., 260., 200., 150., 130., 200., 150., 110., 160.,  & ! 1-10
                                               100., 120.,  90.,  90.,  80.,  90.,  90., 150.,  80., 100.,  & ! 11-20
                                               80.,  80.,  80.,  60.,  60., 120.,  80., 180., 995., 995.,   & ! 21-30
                                               80., 350., 300., 300., 300., 300., 230., 230., 230., 150.,   & ! 31-40
                                               230., 995., 150., 9900./)                                      ! 41-44
real, dimension(44), parameter :: slveg44 = (/ 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,                      & ! 1-10
                                               0., 0., 0., 0., 0., .1, .1, .1, .1, .1,                      & ! 11-20
                                               .1, .2, .4, .2, .1, .1, .1, 0., 0., 0., 0.,                  & ! 21-31
                                               1., 5.5, 3., 1., 3., 3., 3.5, 3., .5, 3.5, .1, 3.5, 0./)       ! 32-44
real, dimension(5), parameter :: froot = (/ .05, .10, .35, .40, .10 /) ! 10/02/99 veg. root distr.
real, dimension(ms), save :: zse = (/ .022, .058, .154, .409, 1.085, 2.872/) ! layer thickness
! so depths of centre of layers: .011, .051, .157, .4385, 1.1855, 3.164
! with base at 4.6     

real, dimension(ms+1), save :: zshh
real, dimension(ms), save :: ww

end module soilv_m
