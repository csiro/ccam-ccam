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


      integer meso,ngwd,nrungcm,newtop,                                  &
     &        kountr,nrad,nvmix,nlocal,                                  &
     &        nhstest,namip,nspecial,newrough,newsoilm,nsib,nsoil,       &
     &        ntaft,ntsea,ntsur,ntsur2,lgwd,newztsea,nglacier,mbd,       &
     &        nbd,kbotdav,kbotu,nbox,nud_p,nud_q,nud_t,nud_uv,nud_hrs,   &
     &        nudu_hrs,ktau,ndi,ndi2,ntau,nperavg,nperday,nmaxpr,nlv,    &
     &        ia,ib,ja,jb,id,jd,idjd,                                    &
     &        io_clim,io_in,io_out,io_rest,io_spec,                      &
     &        nwt,nqg,nrun,nextout,nclim,m_fly,nsemble,tblock,tbave,     &
     &        nurban,nmr,nmlo,ktopdav,nud_sst,nud_sss,kbotmlo,ktopmlo,   &
     &        mloalpha,nud_ouv,nud_sfh,kblock,rescrn,knh,ccycle,iaero,   &
     &        nud_aero,mbd_maxscale,nriver
      real qgmin,                                                        &
     &     aleadfr,av_vmod,vmodmin,snmin,tss_sh,charnock,chn10,zobgin,   &
     &     rlongdn,rlongdx,rlatdn,rlatdx,ds,dt,dtin,timea,panfg,panzo,   &
     &     bpyear,helim,fc2,sigbot_gwd,alphaj,cgmap_offset,cgmap_scale,  &
     &     sigramplow,sigramphigh
      logical diag,localhist,amipo3
      common/parm1/meso,ngwd,nrungcm,newtop,bpyear,iaero,helim,fc2,      &
     &  sigbot_gwd,alphaj,qgmin     ! min value, esp. for stratosphere [1.e-6]

      common/parmradn/kountr,nrad,amipo3   

      common/parmvmix/nvmix,nlocal,cgmap_offset,cgmap_scale

      common/parmtest/nhstest,namip,nsemble,nspecial,rlongdn,rlongdx,    &
     &                rlatdn,rlatdx

      common/parmsfce/newrough,newsoilm,nsib,nsoil,ntsea,ntsur,ntsur2,   &
     &                lgwd,newztsea,aleadfr,av_vmod,vmodmin,snmin,       &
     &                tss_sh,nglacier,charnock,chn10,zobgin,ntaft,       &
     &                panfg,panzo,nurban,nmr,nmlo,rescrn,ccycle,nriver

      common/parmnudg/nbd,kbotdav,kbotu,nbox,nud_p,nud_q,nud_t,nud_uv,   &
     &                nud_hrs,nudu_hrs,mbd,ktopdav,nud_sst,nud_sss,      &
     &                kbotmlo,ktopmlo,mloalpha,nud_ouv,nud_sfh,kblock,   &
     &                nud_aero,mbd_maxscale,sigramplow,sigramphigh

      common/parmtime/ktau,ntau,nperavg,nperday,ds,dt,dtin,timea,nmaxpr, &
     &                diag,nlv,ia,ib,ja,jb,id,jd,idjd,ndi,ndi2,knh

      common/parmio/io_clim,io_in,io_out,io_rest,io_spec,                &  ! type of I/O
     &            nwt,nqg,nrun,nextout,nclim,m_fly,tblock,tbave,         &
     &            localhist

