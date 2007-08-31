!!$ cable_output_text.f90
!!$
!!$ Sample text output module for CABLE land surface scheme offline text driver; 
!!$
!!$ Gab Abramowitz 2007 University of New South Wales/
!!$ CSIRO Marine and Atmospheric Research; gabsun@gmail.com
!!$
!!$ The subroutine in this file simply writes a collection of variables from a 
!!$ single time step to a text file. The choice of variables and matching 
!!$ gnuplot scripts are those used by Eva Kowalczyk for analysis.

SUBROUTINE text_output(ktau, kstart, kend, dels, air, bgc, &
     canopy, met, bal,rad, rough, soil, ssoil, sum_flux, veg, &
     filename_out)
  USE checks_module
  USE carbon_module
  USE physical_constants
  IMPLICIT NONE
  INTEGER(i_d), INTENT(IN)		:: ktau ! integration step number
  INTEGER(i_d), INTENT(IN)	       	:: kstart ! starting value of ktau
  INTEGER(i_d), INTENT(IN)	       	:: kend 
  REAL(r_1), INTENT(IN)			:: dels ! integration time setp (s)
  CHARACTER(LEN=*),INTENT(IN) :: filename_out ! name of file for CABLE output
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
  INTEGER(i_d)                          :: ijd  
  INTEGER(i_d)                 		:: u         ! output unit
  LOGICAL                      		:: is_open   ! Is file open?
  REAL(r_1), DIMENSION(mp)		:: owbtot
  REAL(r_1), DIMENSION(mp)		:: swnet
  REAL(r_1), DIMENSION(mp)		:: lwnet
  INTEGER :: k ! do loop counter

  ! Call conservation routines:
  CALL mass_balance(ktau,dels,ssoil,soil,canopy,met,air,bal)
  CALL energy_balance(ktau,dels,met,rad,canopy,bal,ssoil)

  ! which gridpoint to print results for?
  ijd = 1

  u = 98
  INQUIRE (u, opened=is_open)
  IF (.not. is_open) THEN
     OPEN (u, file=filename_out, status='replace')
     WRITE(u, '(a6,67a14)') &
          'ktau', 'hod', 'fsd', 'fld', 'precip', & !1-5
          'tk', 'pmb', 'ua', 'qv', 'ca', & !6-10
          'swnet', 'lwnet', 'fe', 'fh', & !11-14
          'fnv+fns', 'trad', 'tv', & !15-17
          'acond', 'wb(1,1)', 'wb(1,2)', & !18-20
          'wb(1,3)', 'wb(1,4)', 'wb(1,5)', 'wb(1,6)', & !21-24
          'xx1', 'xx2', 'xx3', & !25-27
          'tgg(1,1)', 'tgg(1,2)', 'tgg(1,3)', 'tgg(1,4)', & !28-31
          'tgg(1,5)', 'tgg(1,6)', & !32-33
          'fns', 'fes', 'fhs', 'ga', & !34-37
          'frp', 'frpw', 'frpr', 'fpn', 'frs', & !38-42
          'fnv','fev','fevc','fevw','fhv', & !43-47
          'precis', 'vlai', 'ebal', 'ebalcan', 'ebalcan2', 'ebalsoil', & !48-53
          'runoff', 'rnof1', 'rnof2', 'delwc', 'snowd', & !54-58
          'wbal', 'albedo(:,1)', 'albedo(:,2)', & !59-61
          'wbal_tot','ebal_tot','wbtot0','osnowd0', & ! 62-65
          'rnoff_tot','precip_tot','evap_tot' ! 66-68
  END IF
  xx1=ssoil%tggsn(:,1)                
  xx2=ssoil%tggsn(:,2)               
  xx3=ssoil%tggsn(:,3)
  where (ssoil%isflag == 0 ) ! no snow
     xx1=ssoil%tgg(:,1)
     xx2=ssoil%tgg(:,2)
     xx3=ssoil%tgg(:,3)
  ELSEWHERE
     xx1=ssoil%tggsn(:,1)                
     xx2=ssoil%tggsn(:,2)               
     xx3=ssoil%tggsn(:,3)
  end where

  swnet=SUM(rad%qcan(:,:,1),2)+SUM(rad%qcan(:,:,2),2)+rad%qssabs
  lwnet=met%fld-sboltz*emleaf*canopy%tv**4 *(1-rad%transd)-rad%flws*rad%transd
  WRITE(u, '(i7,67e14.5)') &
       ktau,met%hod(ijd),met%fsd(ijd),met%fld(ijd),met%precip(ijd), & ! 1-5
       met%tk(ijd), met%pmb(ijd),met%ua(ijd),met%qv(ijd),met%ca(ijd), & ! 6-10
       swnet(ijd),lwnet(ijd), canopy%fe(ijd),canopy%fh(ijd),&	    !11-14
       canopy%fnv(ijd)+canopy%fns(ijd), rad%trad(ijd),canopy%tv(ijd), & !15-17
       1./max(0.0001,rough%rt1(ijd)), ssoil%wb(ijd,1),ssoil%wb(ijd,2),& !18-20
       ssoil%wb(ijd,3),ssoil%wb(ijd,4), ssoil%wb(ijd,5),ssoil%wb(ijd,6),& !21-24
       xx1(ijd),xx2(ijd),xx3(ijd),      &   !25-27
       ssoil%tgg(ijd,1),ssoil%tgg(ijd,2), ssoil%tgg(ijd,3),ssoil%tgg(ijd,4),& !28-31
       ssoil%tgg(ijd,5),ssoil%tgg(ijd,6),	& !32-33
       canopy%fns(ijd),canopy%fes(ijd),canopy%fhs(ijd),canopy%ga(ijd), & !34-37
       canopy%frp(ijd),canopy%frpw(ijd),canopy%frpr(ijd),canopy%fpn(ijd),canopy%frs(ijd),& !38-42
       canopy%fnv(ijd),canopy%fev(ijd),canopy%fevc(ijd),canopy%fevw(ijd),canopy%fhv(ijd), & !43-47
       canopy%precis(ijd),veg%vlai(ijd),bal%ebal(ijd), & !48-50
       canopy%fnv(ijd)-canopy%fhv(ijd)-canopy%fev(ijd), & !51
       bal%drybal(ijd)+bal%wetbal(ijd), & !52
       canopy%fns(ijd)-canopy%fhs(ijd)-canopy%fes(ijd)*ssoil%cls(ijd)-canopy%ga(ijd), & !53
       ssoil%runoff(ijd),ssoil%rnof1(ijd),ssoil%rnof2(ijd),canopy%delwc(ijd),ssoil%snowd(ijd), &
       bal%wbal(ijd),rad%albedo(ijd,1),rad%albedo(ijd,2), & !59-61
       bal%wbal_tot(ijd),bal%ebal_tot(ijd),bal%wbtot0(ijd),bal%osnowd0(ijd), & ! 62-65
       bal%rnoff_tot(ijd), bal%precip_tot(ijd),bal%evap_tot(ijd) ! 66-68
END SUBROUTINE text_output
