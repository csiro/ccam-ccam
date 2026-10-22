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
    
module atmosphere_arrays_m

use arrays_m                               ! Atmosphere dyamics prognostic arrays
use dpsdt_m                                ! Vertical velocity
use epst_m                                 ! Off-centre terms
use nharrs_m                               ! Non-hydrostatic atmosphere arrays
use nlin_m                                 ! Atmosphere non-linear dynamics
use parmdyn_m                              ! Dynamics parameters
use parmhdff_m                             ! Horizontal diffusion parameters
use parmhor_m                              ! Horizontal advection parameters
use parmvert_m                             ! Vertical advection parameters
use savuvt_m                               ! Saved dynamic arrays
use savuv1_m                               ! Saved dynamic arrays
use sbar_m                                 ! Saved dynamic arrays
use tbar2d_m                               ! Atmosphere dynamics reference temperature
use unn_m                                  ! Saved dynamic arrays
use uvbar_m                                ! Saved dynamic arrays
use vecs_m                                 ! Eigenvectors for atmosphere dynamics
use vvel_m                                 ! Additional vertical velocity
use work3sav_m                             ! Water and tracer saved arrays
use xarrs_m                                ! Saved dynamic arrays

implicit none

private

public atmosphere_arrays_init, atmosphere_arrays_end

! arrays_m
public t,u,v,qg
public psl,ps,zs

! dpsdt
public dpsdt, dpsdtb, dpsdtbb

! epst
public epst

! nharrs
public phi,phi_nh,h_nh
public lrestart, lrestart_radiation, lrestart_tracer, always_mspeca
public lrestart_convection

! nlin
public tn,un,vn

! parmdyn
public mex,mfix,mfix_qg,mfix_tr,mfix_aero,mfix_t
public nh,nritch_t,mspec,mup
public nstag,nstagu,ntbar,precon,helmmeth
public nstagoff
public epsp,epsu,epsf,epsh,restol
public maxuv

! parmhdff
public nhor,nhorps,khor,khdif,nhorjlm
public hdiff

! parmhor
public mh_bs, nt_adv

! parmvert
public ntvd

! savuvt
public savt, savpsl
public savs, savu, savv

! savuv1
public savs1,savu1,savv1,savu2,savv2

! sbar
public sbar

! tbar
public tbar2d

! unn
public unn,vnn

! uvbar
public ubar,vbar

! vecs
public emat,einv,bam,bami,qvec,tmat

! vvel
public sdot, dpsldt, wvel
public updraft_helicity

! work3sav
public tsav, qgsav, qfgsav, qlgsav, trsav
public work3sav_init ! must occur after tracers_init

! xarrs
public ux,vx,tx,pslx

contains    
 
subroutine atmosphere_arrays_init(ifull,iextra,kl,nextout,nvmix,epsp)  

implicit none

integer, intent(in) :: ifull, iextra, kl
integer, intent(in) :: nextout, nvmix
real, intent(in) :: epsp

call arrays_init(ifull,iextra,kl)
call dpsdt_init(ifull,epsp)
call epst_init(ifull)
call nharrs_init(ifull,iextra,kl)
call nlin_init(ifull,kl)
call parmhdff_init(kl)
call savuvt_init(ifull,kl)
call savuv1_init(ifull,kl)
call sbar_init(ifull,kl)
call tbar2d_init(ifull)
call unn_init(ifull,kl)
call uvbar_init(ifull,kl)
call vecs_init(kl)
call vvel_init(ifull,iextra,kl)
call xarrs_init(ifull,iextra,kl)

return
end subroutine atmosphere_arrays_init

subroutine atmosphere_arrays_end

implicit none

call arrays_end
call dpsdt_end
call epst_end
call nharrs_end
call nlin_end
call parmhdff_end
call savuvt_end
call savuv1_end
call sbar_end
call tbar2d_end
call unn_end
call uvbar_end
call vecs_end
call vvel_end
call xarrs_end

return
end subroutine atmosphere_arrays_end

end module atmosphere_arrays_m