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
    
module liqwpar_m

implicit none

private
public qlg,qfg ! cloud liquid water, cloud frozen water
public qrg,qsng,qgrg ! rain, snow, graupel
public liqwpar_init,liqwpar_end
public nc,nr,ni,ns
public cdrop_aerosol
public fluxr,fluxi,fluxs,fluxg
public stras_rliq,stras_rice,stras_rrai,stras_rsno
public psnow,psaut,psfw,psfi,praci, &
       piacr,psaci,psacw,psdep,     &
       pssub,pracs,psacr,psmlt,     &
       psmltevp,prain,praut,pracw,  &
       prevp,pgfr,pvapor,pclw,      &
       pladj,pcli,pimlt,pihom,      &
       pidw,piadj,qschg
public leo_pcaut,leo_psaut,leo_pgaut,leo_pgmlt,&
       leo_pgsub,leo_pgacw,leo_pgacr,leo_pgaci,&
       leo_pgacs,leo_psmlt,leo_pssub,leo_psacw,&
       leo_psacr,leo_psaci,leo_pimlt,leo_pisub,&
       leo_piacw,leo_piacr,leo_psure,leo_prevp,&
       leo_pracc,leo_pracs,leo_praci

real, dimension(:,:), allocatable, save :: qlg,qfg
real, dimension(:,:), allocatable, save :: qrg,qsng
real, dimension(:,:), allocatable, save :: qgrg
real, dimension(:,:), allocatable, save :: nc,nr,ni,ns
real, dimension(:,:), allocatable, save :: cdrop_aerosol
real, dimension(:,:), allocatable, save :: fluxr,fluxi,fluxs,fluxg
real, dimension(:,:), allocatable, save :: stras_rliq,stras_rice,stras_rrai,stras_rsno
real, dimension(:,:), allocatable, save :: psnow,psaut,psfw,psfi,praci, &
                                           piacr,psaci,psacw,psdep,     &
                                           pssub,pracs,psacr,psmlt,     &
                                           psmltevp,prain,praut,pracw,  &
                                           prevp,pgfr,pvapor,pclw,      &
                                           pladj,pcli,pimlt,pihom,      &
                                           pidw,piadj,qschg
real, dimension(:,:), allocatable, save :: leo_pcaut,leo_psaut,leo_pgaut,leo_pgmlt,&
                                           leo_pgsub,leo_pgacw,leo_pgacr,leo_pgaci,&
                                           leo_pgacs,leo_psmlt,leo_pssub,leo_psacw,&
                                           leo_psacr,leo_psaci,leo_pimlt,leo_pisub,&
                                           leo_piacw,leo_piacr,leo_psure,leo_prevp,&
                                           leo_pracc,leo_pracs,leo_praci

contains

subroutine liqwpar_init(ifull,iextra,kl,process_rate_mode)

implicit none

integer, intent(in) :: ifull,iextra,kl,process_rate_mode

allocate(qlg(ifull+iextra,kl),qfg(ifull+iextra,kl))
allocate(qrg(ifull,kl),qsng(ifull,kl))
allocate(qgrg(ifull,kl))
!rain, snow, graupel only fall vertically
allocate(nr(ifull,kl),ni(ifull+iextra,kl),ns(ifull,kl))
allocate(cdrop_aerosol(ifull,kl))
allocate(fluxr(ifull,kl),fluxi(ifull,kl),fluxs(ifull,kl),fluxg(ifull,kl))
allocate(stras_rliq(ifull,kl),stras_rice(ifull,kl),stras_rrai(ifull,kl),stras_rsno(ifull,kl))

qlg=0.
qfg=0.
qrg=0.
qsng=0.
qgrg=0.

nr=0.
ni=0.
ns=0.

cdrop_aerosol = 0.

fluxr = 0.
fluxi = 0.
fluxs = 0.
fluxg = 0.

stras_rliq = 0.
stras_rice = 0.
stras_rrai = 0.
stras_rsno = 0.

if (process_rate_mode == 1) then
  allocate(leo_pcaut(ifull,kl),leo_psaut(ifull,kl),leo_pgaut(ifull,kl),leo_pgmlt(ifull,kl),&
           leo_pgsub(ifull,kl),leo_pgacw(ifull,kl),leo_pgacr(ifull,kl),leo_pgaci(ifull,kl),&
           leo_pgacs(ifull,kl),leo_psmlt(ifull,kl),leo_pssub(ifull,kl),leo_psacw(ifull,kl),&
           leo_psacr(ifull,kl),leo_psaci(ifull,kl),leo_pimlt(ifull,kl),leo_pisub(ifull,kl),&
           leo_piacw(ifull,kl),leo_piacr(ifull,kl),leo_psure(ifull,kl),leo_prevp(ifull,kl),&
           leo_pracc(ifull,kl),leo_pracs(ifull,kl),leo_praci(ifull,kl))
  leo_pcaut = 0.
  leo_psaut = 0.
  leo_pgaut = 0.
  leo_pgmlt = 0.
  leo_pgsub = 0.
  leo_pgacw = 0.
  leo_pgacr = 0.
  leo_pgaci = 0.
  leo_pgacs = 0.
  leo_psmlt = 0.
  leo_pssub = 0.
  leo_psacw = 0.
  leo_psacr = 0.
  leo_psaci = 0.
  leo_pimlt = 0.
  leo_pisub = 0.
  leo_piacw = 0.
  leo_piacr = 0.
  leo_psure = 0.
  leo_prevp = 0.
  leo_pracc = 0.
  leo_pracs = 0.
  leo_praci = 0.
else if (process_rate_mode == 2) then
  allocate(psnow(ifull,kl),psaut(ifull,kl),psfw(ifull,kl),psfi(ifull,kl),praci(ifull,kl),&
           piacr(ifull,kl),psaci(ifull,kl),psacw(ifull,kl),psdep(ifull,kl)              ,&
           pssub(ifull,kl),pracs(ifull,kl),psacr(ifull,kl),psmlt(ifull,kl)              ,&
           psmltevp(ifull,kl),prain(ifull,kl),praut(ifull,kl),pracw(ifull,kl)           ,&
           prevp(ifull,kl),pgfr(ifull,kl),pvapor(ifull,kl),pclw(ifull,kl)               ,&
           pladj(ifull,kl),pcli(ifull,kl),pimlt(ifull,kl),pihom(ifull,kl)               ,&
           pidw(ifull,kl),piadj(ifull,kl),qschg(ifull,kl))
  psnow = 0.
  psaut = 0.
  psfw  = 0.
  psfi  = 0.
  praci = 0.
  piacr = 0.
  psaci = 0.
  psacw = 0.
  psdep = 0.
  pssub = 0.
  pracs = 0.
  psacr = 0.
  psmlt = 0.
  psmltevp = 0.
  prain = 0.
  praut = 0.
  pracw = 0.
  prevp = 0.
  pgfr  = 0.
  pvapor= 0.
  pclw  = 0.
  pladj = 0.
  pcli  = 0.
  pimlt = 0.
  pihom = 0.
  pidw  = 0.
  piadj = 0.
  qschg = 0.
end if

return
end subroutine liqwpar_init

subroutine liqwpar_end

implicit none

deallocate(qlg,qfg)
deallocate(qrg,qsng)
deallocate(qgrg)
deallocate(nr,ni,ns)
deallocate(cdrop_aerosol)
deallocate(fluxr,fluxi,fluxs,fluxg)
deallocate(stras_rliq,stras_rice,stras_rrai,stras_rsno)

if ( allocated(psnow) ) then
  deallocate(psnow,psaut,psfw,psfi,praci,&
           piacr,psaci,psacw,psdep              ,&
           pssub,pracs,psacr,psmlt              ,&
           psmltevp,prain,praut,pracw           ,&
           prevp,pgfr,pvapor,pclw               ,&
           pladj,pcli,pimlt,pihom               ,&
           pidw,piadj,qschg)
end if
if ( allocated(leo_pcaut) ) then
  deallocate(leo_pcaut,leo_psaut,leo_pgaut,leo_pgmlt,&
             leo_pgsub,leo_pgacw,leo_pgacr,leo_pgaci,&
             leo_pgacs,leo_psmlt,leo_pssub,leo_psacw,&
             leo_psacr,leo_psaci,leo_pimlt,leo_pisub,&
             leo_piacw,leo_piacr,leo_psure,leo_prevp,&
             leo_pracc,leo_pracs,leo_praci)
end if
return
end subroutine liqwpar_end

end module liqwpar_m
