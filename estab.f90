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
    
! Routines for saturated air used by cloud microphysics routines

module estab

implicit none

private
public establ,estabi,qsat,qsati,estab_init
public esdiffx,pow75

real, dimension(:), allocatable, save :: table, tablei
real, dimension(:), allocatable, save :: esdiff

interface qsat
  module procedure qsat_s, qsat_v
end interface qsat
interface qsati
  module procedure qsati_s, qsati_v
end interface qsati
interface establ
  module procedure establ_s, establ_v
end interface establ
interface estabi
  module procedure estabi_s, estabi_v
end interface estabi
interface tdiff
  module procedure tdiff_s, tdiff_v
end interface tdiff
interface tdiffx
  module procedure tdiffx_s, tdiffx_v
end interface tdiffx
interface esdiffx
  module procedure esdiffx_s, esdiffx_v
end interface esdiffx
interface pow75
  module procedure pow75_s, pow75_v
end interface

contains

subroutine estab_init

implicit none

allocate(table(0:220),tablei(0:220))
allocate(esdiff(-40:2))

table(0:99)=(/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9,                                    & !-146C
               6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9,                                 & !-141C
               36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9,                              & !-136C
               0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648,      & !-131C
               0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774,      & !-126C
               0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081,       & !-121C
               0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866,           & !-116C
               0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280,             & !-111C
               0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951,                & !-106C
               0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143,                 & !-101C
               .001403, .001719, .002101, .002561, .003117, .003784,                 & !-95C
               .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658,     & !-87C
               .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577,     & !-78C
               .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032,       & !-69C
               .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080,        & !-60C
               1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476/)         !-51C

table(100:220)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098,             & !-43C
                  10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88,     & !-34C
                  27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85,  & !-24C 
                  77.09, 85.02, 93.70, 103.20, 114.66, 127.20, 140.81, 155.67,       & !-16C
                  171.69, 189.03, 207.76, 227.96 , 249.67, 272.98, 298.00, 324.78,   & !-8C
                  353.41, 383.98, 416.48, 451.05, 487.69, 526.51, 567.52, 610.78,    & !0C
                  656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2,    & !8C
                  1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3,    & !16C
                  1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1,    & !24C
                  3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1,    & !32C
                  5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7,    & !40C
                  7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0,          & !47C
                  11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0,     & !54C
                  15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0,     & !61C
                  21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0,     & !68C
                  29845.0, 31169.0/)                                                   !70C

tablei(0:99)=(/ 1.e-9, 1.e-9, 2.e-9, 3.e-9, 4.e-9,                                   & !-146C
                6.e-9, 9.e-9, 13.e-9, 18.e-9, 26.e-9,                                & !-141C
                36.e-9, 51.e-9, 71.e-9, 99.e-9, 136.e-9,                             & !-136C
                0.000000188, 0.000000258, 0.000000352, 0.000000479, 0.000000648,     & !-131C
                0.000000874, 0.000001173, 0.000001569, 0.000002090, 0.000002774,     & !-126C
                0.000003667, 0.000004831, 0.000006340, 0.000008292, 0.00001081,      & !-121C
                0.00001404, 0.00001817, 0.00002345, 0.00003016, 0.00003866,          & !-116C
                0.00004942, 0.00006297, 0.00008001, 0.0001014, 0.0001280,            & !-111C
                0.0001613, 0.0002026, 0.0002538, 0.0003170, 0.0003951,               & !-106C
                0.0004910, 0.0006087, 0.0007528, 0.0009287, 0.001143,                & !-101C
                .001403, .001719, .002101, .002561, .003117, .003784,                & !-95C
                .004584, .005542, .006685, .008049, .009672,.01160,.01388,.01658,    & !-87C
                .01977, .02353, .02796,.03316,.03925,.04638,.05472,.06444,.07577,    & !-78C
                .08894, .1042, .1220, .1425, .1662, .1936, .2252, .2615, .3032,      & !-69C
                .3511, .4060, .4688, .5406, .6225, .7159, .8223, .9432, 1.080,       & !-60C
                1.236, 1.413, 1.612, 1.838, 2.092, 2.380, 2.703, 3.067, 3.476 /)       !-51C

tablei(100:220)=(/ 3.935,4.449, 5.026, 5.671, 6.393, 7.198, 8.097, 9.098,            & !-43C
                   10.21, 11.45, 12.83, 14.36, 16.06, 17.94, 20.02, 22.33, 24.88,    & !-34C
                   27.69, 30.79, 34.21, 37.98, 42.13, 46.69,51.70,57.20,63.23,69.85, & !-24C 
                   77.09, 85.02, 93.70, 103.06, 113.40, 124.68, 136.98, 150.39,      & !-16C
                   164.99, 180.88, 198.16, 216.94, 237.34, 259.47, 283.49, 309.51,   & !-8C
                   337.71, 368.23, 401.25, 436.96, 475.54, 517.21, 562.19, 610.70,   & !0C
                   656.62, 705.47, 757.53, 812.94, 871.92, 934.65, 1001.3, 1072.2,   & !8C
                   1147.4, 1227.2, 1311.9, 1401.7, 1496.9, 1597.7, 1704.4, 1817.3,   & !16C
                   1936.7, 2063.0, 2196.4, 2337.3, 2486.1, 2643.0, 2808.6, 2983.1,   & !24C
                   3167.1, 3360.8, 3564.9, 3779.6, 4005.5, 4243.0, 4492.7, 4755.1,   & !32C
                   5030.7, 5320.0, 5623.6, 5942.2, 6276.2, 6626.4, 6993.4, 7377.7,   & !40C
                   7780.2, 8201.5, 8642.3, 9103.4, 9585.5, 10089.0, 10616.0,         & !47C
                   11166.0, 11740.0, 12340.0, 12965.0, 13617.0, 14298.0, 15007.0,    & !54C
                   15746.0, 16516.0, 17318.0, 18153.0, 19022.0, 19926.0, 20867.0,    & !61C
                   21845.0, 22861.0, 23918.0, 25016.0, 26156.0, 27340.0, 28570.0,    & !68C
                   29845.0, 31169.0/)                                                  !70C

esdiff(-40:2)=(/ 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,  &
                 13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61, &
                 22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27, &
                 26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65, &
                 0.08, 0.0, 0.0 /)

return
end subroutine estab_init

function tdiff_s(t_) result(ans)
implicit none
real, intent(in) :: t_
real ans
! TDIFF is difference between T and 123.16, subject to 0 <= TDIFF <= 220
ans=min(max( t_-123.16, 0.), 219.)
end function tdiff_s

function tdiff_v(t_) result(ans)
implicit none
real, dimension(:), intent(in) :: t_
real, dimension(size(t_)) :: ans
! TDIFF is difference between T and 123.16, subject to 0 <= TDIFF <= 220
ans=min(max( t_-123.16, 0.), 219.)
end function tdiff_v

function tdiffx_s(tx_) result(ans)
implicit none
include 'const_phys.h'
real, intent(in) :: tx_
real ans
ans=min(max( tx_-tfrz, -40.), 1.)
end function tdiffx_s

function tdiffx_v(tx_) result(ans)
implicit none
include 'const_phys.h'
real, dimension(:), intent(in) :: tx_
real, dimension(size(tx_)) :: ans
ans=min(max( tx_-tfrz, -40.), 1.)
end function tdiffx_v

function establ_s(t_) result(ans)
implicit none
integer tpos
real, intent(in) :: t_
real ans
real tstore, tfrac
tstore = tdiff(t_)
tfrac = tstore-aint(tstore)
tpos = int(tstore)
! Arithmetic statement functions to replace call to establ.
! T is temp in Kelvin, which should lie between 123.16 and 343.16;
ans = (1.-tfrac)*table(tpos)+ tfrac*table(tpos+1)
end function establ_s

function establ_v(t_) result(ans)
implicit none
real, dimension(:), intent(in) :: t_
real, dimension(size(t_)) :: ans
real, dimension(size(t_)) :: tstore, tfrac
integer, dimension(size(t_)) :: tpos
tstore = tdiff(t_)
tfrac = tstore - aint(tstore)
tpos = int(tstore)
! Arithmetic statement functions to replace call to establ.
! T is temp in Kelvin, which should lie between 123.16 and 343.16;
ans = (1.-tfrac)*table(tpos)+ tfrac*table(tpos+1)
end function establ_v

function qsat_s(pp_,t_) result(ans)
implicit none
include 'const_phys.h'
real, intent(in) :: pp_, t_
real ans      
real estore
estore = establ(t_)
ans = epsil*estore/max(pp_-estore,.1) !jlm strato
! ans = epsil*establ(t_)/(pp_-establ(t_)) !Usual formula
! ans = epsil*establ(t_)/pp_ !Consistent with V4-5 to V4-7
end function qsat_s

function qsat_v(pp_,t_) result(ans)
implicit none
include 'const_phys.h'
real, dimension(:), intent(in) :: pp_, t_
real, dimension(size(pp_)) :: ans
real, dimension(size(pp_)) :: estore
estore = establ(t_)
ans = epsil*estore/max(pp_-estore,.1) !jlm strato
! ans = epsil*establ(t_)/(pp_-establ(t_)) !Usual formula
! ans = epsil*establ(t_)/pp_ !Consistent with V4-5 to V4-7
end function qsat_v

function estabi_s(t_) result(ans)
implicit none
integer tpos
real, intent(in) :: t_
real ans
real tstore, tfrac
tstore = tdiff(t_)
tfrac = tstore - aint(tstore)
tpos = int(tstore)
ans = (1.-tfrac)*tablei(tpos)+ tfrac*tablei(tpos+1)
end function estabi_s

function estabi_v(t_) result(ans)
implicit none
real, dimension(:), intent(in) :: t_
real, dimension(size(t_)) :: ans 
real, dimension(size(t_)) :: tstore, tfrac
integer, dimension(size(t_)) :: tpos
tstore = tdiff(t_)
tfrac = tstore - aint(tstore)
tpos = int(tstore)
ans = (1.-tfrac)*tablei(tpos)+ tfrac*tablei(tpos+1)
end function estabi_v

function qsati_s(pp_,t_) result(ans)
implicit none
include 'const_phys.h'
real, intent(in) :: pp_, t_
real ans
real estore
estore = estabi(t_)
ans = epsil*estore/max(pp_-estore,.1) !jlm strato
! ans = epsil*estabi(t_)/(pp_-estabi(t_)) !Usual formula
! ans = epsil*estabi(t_)/pp_ !Consistent with V4-5 to V4-7
end function qsati_s

function qsati_v(pp_,t_) result(ans)
implicit none
include 'const_phys.h'
real, dimension(:), intent(in) :: pp_, t_
real, dimension(size(pp_)) :: ans
real, dimension(size(pp_)) :: estore
estore = estabi(t_)
ans = epsil*estore/max(pp_-estore,.1) !jlm strato
! ans = epsil*estabi(t_)/(pp_-estabi(t_)) !Usual formula
! ans = epsil*estabi(t_)/pp_ !Consistent with V4-5 to V4-7
end function qsati_v

function esdiffx_s(tx_) result(ans)
implicit none
integer tpos
real, intent(in) :: tx_
real ans
real tstore, tfrac
tstore = tdiffx(tx_)
tfrac = tstore - aint(tstore)
tpos = int(tstore)
ans = (1.-tfrac)*esdiff(tpos)+tfrac*esdiff(tpos+1)
end function esdiffx_s

function esdiffx_v(tx_) result(ans)
implicit none
real, dimension(:), intent(in) :: tx_
real, dimension(size(tx_)) :: ans
real, dimension(size(tx_)) :: tstore, tfrac
integer, dimension(size(tx_)) :: tpos
tstore = tdiffx(tx_)
tfrac = tstore - aint(tstore)
tpos = int(tstore)
ans = (1.-tfrac)*esdiff(tpos)+tfrac*esdiff(tpos+1)
end function esdiffx_v

function pow75_s(x) result(ans)
implicit none
real, intent(in) :: x
real ans
ans=sqrt(x*sqrt(x))
end function pow75_s

function pow75_v(x) result(ans)
implicit none
real, dimension(:), intent(in) :: x
real, dimension(size(x)) :: ans
ans=sqrt(x*sqrt(x))
end function pow75_v

end module estab