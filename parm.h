      integer meso,nps,npsav,nsd,nem,ngwd,nrungcm,newtop,nhor,nhorps,
     .        khor,khdif,kountr,ndiur,nrad,nvmix,nlocal,
     .        nhstest,namip,nspecial,nsib,nsoil,newsoilm,
     .        ntsea,ntsur,ntsur2,lgwd,newztsea,nglacier,
     .        nbd,kbotdav,nbox,nud_p,nud_q,nud_t,nud_uv,nud_hrs,
     .        ktau,ndi,ndi2,ntau,nperday,nmaxpr,nlv,
     .        ia,ib,ja,jb,id,jd,idjd,
     .        io_clim,io_in,io_out,io_rest,io_spec,
     .        nwt,kwt,nqg,nrun,nrunx,nextout,nclim,nfly
      real qgmin,hdiff,hdifmax,rlong0,rlat0,schmidt,schm13,
     .     aleadfr,av_vmod,vmodmin,snmin,tss_sh,ds,dt,dtin,timea
      logical diag
      common/parm1/meso,nps,npsav,nsd,nem,ngwd,nrungcm,newtop 
     .  ,qgmin        ! min value, esp. for stratosphere           [1.e-6]
 
      common/parmhdff/nhor,nhorps,hdiff(kl),khor,khdif,hdifmax

      common/parmradn/kountr,ndiur,nrad   

      common/parmvmix/nvmix,nlocal

      common/parmtest/nhstest,namip,nspecial

      common/parmgeom/rlong0,rlat0,schmidt,schm13

      common/parmsfce/nsib,nsoil,newsoilm,ntsea,ntsur,ntsur2,
     .                lgwd,newztsea,aleadfr,av_vmod,vmodmin,snmin,
     .                tss_sh,nglacier

      common/parmnudg/nbd,kbotdav,nbox,nud_p,nud_q,nud_t,nud_uv,nud_hrs

      common/parmtime/ktau,ntau,nperday,ds,dt,dtin,timea,nmaxpr,    
     .                diag,nlv,ia,ib,ja,jb,id,jd,idjd,ndi,ndi2

      common/parmio/io_clim,io_in,io_out,io_rest,io_spec,    ! type of I/O
     .              nwt,kwt,nqg,nrun,nrunx,nextout,nclim,nfly  

