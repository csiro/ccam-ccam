      real, dimension(ifull+iextra) :: em, emu, emv, f, fu, fv, zs
      real, dimension(ifull) :: dmdx, dmdy
      logical land(ifull)

      common /map/ em,emu,emv,f,fu,fv,dmdx,dmdy,zs,land
