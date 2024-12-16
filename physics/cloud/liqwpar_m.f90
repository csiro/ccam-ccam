! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2023 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module liqwpar_m

implicit none

private
public qlg,qfg ! cloud liquid water, cloud frozen water
public qrg,qsng,qgrg                  ! rain, snow, graupel
public nr, ni, ns                     ! 2nd moment terms
public stras_rliq, stras_rice, stras_rsno, stras_rrai ! droplet radius
public stras_cliq, stras_cice
public liqwpar_init,liqwpar_end

real, dimension(:,:), allocatable, save :: qlg,qfg
real, dimension(:,:), allocatable, save :: qrg,qsng,qgrg
real, dimension(:,:), allocatable, save :: nr, ni, ns
real, dimension(:,:), allocatable, save :: stras_rliq, stras_rice, stras_rsno, stras_rrai
real, dimension(:,:), allocatable, save :: stras_cliq, stras_cice

contains

subroutine liqwpar_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate( qlg(ifull+iextra,kl), qfg(ifull+iextra,kl) )
allocate( qrg(ifull+iextra,kl), qsng(ifull+iextra,kl), qgrg(ifull+iextra,kl) )
allocate( nr(ifull,kl), ni(ifull+iextra,kl), ns(ifull,kl) )
allocate(stras_rliq(ifull,kl),stras_rice(ifull,kl),stras_rsno(ifull,kl),stras_rrai(ifull,kl))
allocate(stras_cliq(ifull,kl),stras_cice(ifull,kl))
qlg = 0. ! liquid water for cloud
qfg = 0. ! frozen water for cloud
qrg = 0. ! precipitating rain
qsng = 0. ! precipitating snow
qgrg = 0. ! precipitation graupel
nr = 0. ! number concentration for rain
ni = 0. ! number concentration for ice/graupel
ns = 0. ! number concentration for snow
stras_rliq = 0.
stras_rice = 0.
stras_rsno = 0.
stras_rrai = 0.
stras_cliq = 0.
stras_cice = 0.

return
end subroutine liqwpar_init

subroutine liqwpar_end

implicit none

deallocate( qlg, qfg)
deallocate( qrg, qsng, qgrg )
deallocate( ni, nr, ns )
deallocate(stras_rliq,stras_rice,stras_rsno,stras_rrai)
deallocate(stras_cliq,stras_cice)

return
end subroutine liqwpar_end

end module liqwpar_m