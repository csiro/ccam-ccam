module extraout_m

implicit none

private
public cloudlo,cloudmi,cloudhi,cloudtot
public rgsave,rtsave,sintsave,sgsave
public rtclsave,sgclsave,taux,tauy,ustar
public swrsave,fbeamvis,fbeamnir
public u10_3hr,v10_3hr,tscr_3hr,rh1_3hr
public extraout_init,extraout_end

real, dimension(:), allocatable, save :: cloudlo,cloudmi,cloudhi,cloudtot
real, dimension(:), allocatable, save :: rgsave,rtsave,sintsave,sgsave
real, dimension(:), allocatable, save :: rtclsave,sgclsave,taux,tauy,ustar
real, dimension(:), allocatable, save :: swrsave,fbeamvis,fbeamnir
real, dimension(:,:), allocatable, save :: u10_3hr,v10_3hr,tscr_3hr,rh1_3hr

contains

subroutine extraout_init(ifull,iextra,kl)

implicit none

integer, intent(in) :: ifull,iextra,kl

allocate(cloudlo(ifull),cloudmi(ifull),cloudhi(ifull),cloudtot(ifull))
allocate(rgsave(ifull),rtsave(ifull),sintsave(ifull),sgsave(ifull))
allocate(rtclsave(ifull),sgclsave(ifull),taux(ifull),tauy(ifull),ustar(ifull))
allocate(swrsave(ifull),fbeamvis(ifull),fbeamnir(ifull))
allocate(u10_3hr(ifull,8),v10_3hr(ifull,8),tscr_3hr(ifull,8),rh1_3hr(ifull,8))

return
end subroutine extraout_init

subroutine extraout_end

implicit none

deallocate(cloudlo,cloudmi,cloudhi,cloudtot)
deallocate(rgsave,rtsave,sintsave,sgsave)
deallocate(rtclsave,sgclsave,taux,tauy,ustar)
deallocate(swrsave,fbeamvis,fbeamnir)
deallocate(u10_3hr,v10_3hr,tscr_3hr,rh1_3hr)

return
end subroutine extraout_end

end module extraout_m