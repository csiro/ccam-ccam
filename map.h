      real, dimension(ifull+iextra) :: em, emu, emv, f, fu, fv, zs
      real, dimension(ifull) :: dmdx, dmdy, dmdxv, dmdyu
      logical land(ifull)

      common /map/ em,emu,emv,f,fu,fv,dmdx,dmdy,dmdxv,dmdyu,zs,land
