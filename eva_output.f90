SUBROUTINE eva_output(ktau, kstart, kend, ktauyear, dels, air, bgc, &
          canopy, met, bal, &
          rad, rough, soil, ssoil, sum_flux, veg)
  USE checks_module
  USE carbon_module
  IMPLICIT NONE
  INTEGER(i_d), INTENT(IN)		:: ktau ! integration step number
  INTEGER(i_d), INTENT(IN)	       	:: kstart ! starting value of ktau
  INTEGER(i_d), INTENT(IN)	       	:: kend 
  INTEGER(i_d), INTENT(IN)	       	:: ktauyear
  REAL(r_1), INTENT(IN)			:: dels ! integration time setp (s)
  TYPE (air_type), INTENT(INOUT)	:: air
  TYPE (bgc_pool_type), INTENT(INOUT)	:: bgc	
  TYPE (canopy_type), INTENT(INOUT)	:: canopy
  TYPE (met_type), INTENT(INOUT) 	:: met
  TYPE (balances_type), INTENT(INOUT)   :: bal
  TYPE (radiation_type), INTENT(INOUT) 	:: rad
  TYPE (roughness_type), INTENT(INOUT) 	:: rough
  TYPE (soil_parameter_type), INTENT(INOUT)	:: soil	
  TYPE (soil_snow_type), INTENT(INOUT)	:: ssoil
  TYPE (sum_flux_type), INTENT(INOUT)	:: sum_flux
  TYPE (veg_parameter_type), INTENT(INOUT)	:: veg
  REAL(r_1), dimension(mp)              :: xx1,xx2,xx3
  INTEGER(i_d) :: ijd  
  INTEGER(i_d)                 		:: u         ! output unit
  LOGICAL                      		:: is_open   ! Is file open?
  REAL(r_1), DIMENSION(mp)		:: owbtot
  REAL(r_1), DIMENSION(mp)		:: swnet
  REAL(r_1), DIMENSION(mp)		:: lwnet
  INTEGER :: k ! do loop counter
  INTEGER(i_d),DIMENSION(10) :: listp ! 


  include 'morepbl.h'
!  include 'parm.h'


  ! Call conservation routines:
!  CALL mass_balance(ktau,dels,ssoil,soil,canopy,met,air,bal)
!  CALL energy_balance(ktau,dels,met,rad,canopy,bal,ssoil)

  ijd = 1
  listp(1)=2726
  listp(2)=3002
  listp(3)=3009
  listp(4)=2468
  listp(5)=2560
  listp(6)=1547
  listp(7)=3688

  
!  u = 102
!  INQUIRE (u, opened=is_open)
!  IF (.not. is_open) THEN
!     OPEN (u, file='f102.txt', status='replace')
!     WRITE(u, '(a6,a7,a6,5a7,4a7,2a7,5a9,a12,3a7,2a7)') &
!          'ktau', 'fsd', 'fld', 'fnv', 'fev', 'fevc', 'fevw', 'fhv', & 
!          'fns', 'fes', 'fhs', 'ga', 'fe', 'fh', 'drybal', 'wetbal', &
!          'ebalt', 'ebal', 'ebal', 'ebal_tot', 'trad', 'tv', 'tss', 'vlaiw', 'transd'
!  END IF
!  WRITE (u, '(i6,f7.1,f6.1,5f7.1,4f7.1,2f7.1,5f9.1,f12.1,3f7.1,2f7.2)') &
!       ktau,met%fsd(ijd),met%fld(ijd), &
!       canopy%fnv(ijd),canopy%fev(ijd),canopy%fevc(ijd),canopy%fevw(ijd),canopy%fhv(ijd), &   !43-47
!       canopy%fns(ijd),canopy%fes(ijd),canopy%fhs(ijd),canopy%ga(ijd), &
!       canopy%fe(ijd),canopy%fh(ijd), &
!       bal%drybal(ijd),bal%wetbal(ijd), &
!       sum(rad%qcan(ijd,:,1))+sum(rad%qcan(ijd,:,2))+rad%qssabs(ijd) &
!       +met%fld(ijd) - sboltz*emsoil*rad%trad(ijd)**4    &
!       -canopy%fev(ijd) - canopy%fes(ijd)*ssoil%cls(ijd) -canopy%fh(ijd) -canopy%ga(ijd), &
!       met%fsd(ijd)*(1.-(rad%albedo(ijd,1)+rad%albedo(ijd,2))/2) &
!       +met%fld(ijd)-sboltz*emleaf*canopy%tv(ijd)**4*(1.-rad%transd(ijd))-rad%flws(ijd)*rad%transd(ijd) &
!       -canopy%fev(ijd) - canopy%fes(ijd)*ssoil%cls(ijd) -canopy%fh(ijd) -canopy%ga(ijd), &
!       bal%ebal(ijd),bal%ebal_tot(ijd), &
!       rad%trad(ijd),canopy%tv(ijd),ssoil%tss(ijd),veg%vlaiw(ijd),rad%transd(ijd)
!  ! New file-- - -- -
!  u = 101
!  INQUIRE (u, opened=is_open)
!  IF (.not. is_open) THEN
!     OPEN (u, file='f101.txt', status='replace')
!     WRITE(u, '(a6,23(1x,a12))') &
!          'ktau', 'precip', 'cansto', 'delwc', 'rnof1', 'rnof2', 'fe*dels/rlam', &
!          'fev*dels/rlam', 'fevc*dels/rlam', 'fevw*dels/rlam', 'fes/cls*dels/rlam', &
!          'owbtot', 'wbtot-owbtot', 'wbal', 'precip_tot', 'evap_tot', 'trnoff', &
!          'wbtot', 'wbal_tot', 'ebal_tot'
!  END IF
!  WRITE(u, '(i6,12f13.6,11f13.5)') &
!       ktau,met%precip(ijd),canopy%cansto(ijd),canopy%delwc(ijd), &
!       ssoil%rnof1(ijd),ssoil%rnof2(ijd),canopy%fe(ijd)*dels/air%rlam(ijd), &
!       canopy%fev(ijd)*dels/air%rlam(ijd),canopy%fevc(ijd)*dels/air%rlam(ijd), &
!       canopy%fevw(ijd)*dels/air%rlam(ijd),  &
!       canopy%fes(ijd)/ssoil%cls(ijd)*dels/air%rlam(ijd),owbtot(ijd), &
!       ssoil%wbtot(ijd)-owbtot(ijd),bal%wbal(ijd), &
!       bal%precip_tot(ijd), bal%evap_tot(ijd),bal%rnoff_tot(ijd),  &
!       ssoil%wbtot(ijd),bal%wbal_tot(ijd),bal%ebal_tot(ijd)
!
  u = 98
  INQUIRE (u, opened=is_open)
  IF (.not. is_open) THEN
     OPEN (u, file='f98.txt', status='replace')
     WRITE(u, '(a6,60a14)') &
          'ktau', 'hod', 'fsd', 'fld', 'precip', & !1-5
          'tk', 'pmb', 'ua', 'qv', 'ca', & !6-11
          'swnet', 'lwnet', &
          'fe', 'fh', & !12-14
          'fnv+fns', 'trad', 'tv', & !15-17
          'acond', 'wb(1,1)', 'wb(1,2)', & !18-20
          'wb(1,3)', 'wb(1,4)', 'wb(1,5)', 'wb(1,6)', & !21-24
          'xx1', 'xx2', 'xx3', & !25-27
          'tgg(1,1)', 'tgg(1,2)', 'tgg(1,3)', 'tgg(1,4)', & !28-31
          'tgg(1,5)', 'tgg(1,6)', & !32-33
          'fns', 'fes', 'fhs', 'ga', & !34-37
          'frp', 'frpw', 'frpr', 'fpn', 'frs', & !38-42
          'pfrs', &
          'fnv', 'fev', 'fevc', 'fevw', 'fhv', & !43-47
          'precis', 'vlai', 'ebal', 'ebalcan', 'ebalcan2', 'ebalsoil', & !48-53
          'runoff', 'rnof1', 'rnof2', 'delwc', 'snowd', & !54-58
          'wbal', 'albedo(:,1)', 'albedo(:,2)' !59-61
  END IF
  xx1=ssoil%tggsn(:,1)                
!  xx2=ssoil%tggsn(:,2)               
  xx3=ssoil%tggsn(:,3)
  xx2=ssoil%tggsn(:,2)               
  where (ssoil%isflag == 0 )
     xx1=ssoil%tgg(:,1)
     xx2=ssoil%tgg(:,2)
     xx3=ssoil%tgg(:,3)
  end where

  swnet=SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs
  lwnet=met%fld-sboltz*emleaf*canopy%tv**4 *(1-rad%transd)-rad%flws*rad%transd
!  do k=1,mp 
!    swnet(k)=SUM(rad%qcan(k,:,1))+SUM(rad%qcan(k,:,2))+rad%qssabs(k)
!    lwnet(k)=met%fld(k)-sboltz*emleaf*canopy%tv(k)**4 *(1-rad%transd(k))-mout%flws(k)*rad%transd(k)
!  enddo   
  do k=1,7
   ijd =listp(k)

     WRITE(u, '(i7,65e14.5)') &
          ktau,met%hod(ijd),met%fsd(ijd),met%fld(ijd),met%precip(ijd), & ! 1-5
          met%tk(ijd), met%pmb(ijd),met%ua(ijd),met%qv(ijd),met%ca(ijd), & ! 6-11
          swnet(ijd),lwnet(ijd), &
          canopy%fe(ijd),canopy%fh(ijd),&				  !12-14
          canopy%fnv(ijd)+canopy%fns(ijd), rad%trad(ijd),canopy%tv(ijd), & !15-17
          1./max(0.0001,rough%rt1(ijd)), ssoil%wb(ijd,1),ssoil%wb(ijd,2),& !18-20
          ssoil%wb(ijd,3),ssoil%wb(ijd,4), ssoil%wb(ijd,5),ssoil%wb(ijd,6),& !21-24
          xx1(ijd),xx2(ijd),xx3(ijd),      &   !25-27
          ! 0.0, 0.0, 0.0,	  &			  !25-27
          ssoil%tgg(ijd,1),ssoil%tgg(ijd,2), ssoil%tgg(ijd,3),ssoil%tgg(ijd,4),& !28-31
          ssoil%tgg(ijd,5),ssoil%tgg(ijd,6),	& !32-33
          canopy%fns(ijd),canopy%fes(ijd),canopy%fhs(ijd),canopy%ga(ijd), & !34-37
          canopy%frp(ijd),canopy%frpw(ijd),canopy%frpr(ijd),canopy%fpn(ijd),canopy%frs(ijd),& !38-42
          canopy%fnv(ijd),canopy%fev(ijd),canopy%fevc(ijd),canopy%fevw(ijd),canopy%fhv(ijd), & !43-47
          canopy%precis(ijd),veg%vlai(ijd),bal%ebal(ijd), &
          canopy%fnv(ijd)-canopy%fhv(ijd)-canopy%fev(ijd), &
          bal%drybal(ijd)+bal%wetbal(ijd), &
          canopy%fns(ijd)-canopy%fhs(ijd)-canopy%fes(ijd)*ssoil%cls(ijd)-canopy%ga(ijd), &
          ssoil%runoff(ijd),ssoil%rnof1(ijd),ssoil%rnof2(ijd),canopy%delwc(ijd),ssoil%snowd(ijd), &
          bal%wbal(ijd),rad%albedo(ijd,1),rad%albedo(ijd,2) , &
         1./max(0.001,ssoil%rtsoil(ijd)), canopy%us(ijd),canopy%cduv(ijd), pblh(ijd)
        enddo 

!  u = 100
!  INQUIRE (u, opened=is_open)
!  IF (.not. is_open) THEN
!     OPEN (u, file='f100.txt', status='replace')
!     WRITE(u, *) 'clitt ',  'frs ', 'cplant(:,1)', 'cplant(:,2)','cplant(:,3)', &
!          'csoil(:,1) ', 'csoil(:,2)' , 'coef_cd'
!  END IF
!  !	WRITE(u, '(11e13.5)') clitt, cfrts, cfwd, cfsf, cfrts, canopy%frs, &
!  !		bgc%csoil(:,1), bgc%csoil(:,2) ,coef_cd, coef_cdnew
!  WRITE(u, '(11e13.5)') clitt(ijd),canopy%frs(ijd),bgc%cplant(ijd,1),bgc%cplant(ijd,2),bgc%cplant(ijd,3), &
!       bgc%csoil(ijd,1), bgc%csoil(ijd,2) ,coef_cd(ijd), coef_cdnew(ijd)
! 
!  u = 103
!  INQUIRE (u, opened=is_open)
!  IF (.not. is_open) THEN
!     OPEN (u, file='f103.txt', status='replace')
!     WRITE(u, *) 'frs '
!  END IF
!  WRITE(u, '(9e13.5)') canopy%frs(ijd)
!
  !print*, canopy%fevc



!  IF(ktau==kend) THEN
!     print*, 'Total precip - evap - rnoff - change in wb - change in snow: ', &
!          bal%precip_tot(ijd)-bal%evap_tot(ijd)-bal%rnoff_tot(ijd)+  &
!          (bal%wbtot0(ijd)-ssoil%wbtot(ijd))+bal%osnowd0(ijd)-ssoil%snowd(ijd)
!     print*, 'Cumulative water balance: ', bal%wbal_tot
!     print*, 'Total components:'
!     print*, 'Total precip, evap, rnoff, soil water change, final wb, init wb : ',& 
!          bal%precip_tot(ijd), bal%evap_tot(ijd),bal%rnoff_tot(ijd), &
!          (bal%wbtot0(ijd)-ssoil%wbtot(ijd)),ssoil%wbtot(ijd),bal%wbtot0(ijd)
!     print*
!     print*, 'Total energy balance: ', bal%ebal_tot
!     print*
!     print*, 'Sum photosynthesis, plant respiration, soil respiration:'
!     print*, sum_flux%sumpn, sum_flux%sumrp, sum_flux%sumrs
!  END IF

END SUBROUTINE eva_output
