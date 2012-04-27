
!========================================================================!
!=== PURPOSE: to enable accessibility of fundamental vars in global   ===!
!=== model thru out program (mainly CABLE).                           ===!
!=== USE: use module in subroutines (SRs) at top level of global model===!
!--- 2 define cable_timestep_ data with vars peculiar to global model ===!
!=== optionally call alias_ SR to give name smore familiar to cable   ===!
!=== people. then again use this module in  any subsequent SR in which===!
!=== you want to access this data.                                    ===!  
!========================================================================!

module cable_common_module
   implicit none 

   !---allows reference to "gl"obal timestep in run (from atm_step)
   !---total number of timesteps, and processing node 
   integer, save :: ktau_gl, kend_gl, knode_gl, kwidth_gl
   
   !---CABLE runtime switches def in this type
   type kbl_internal_switches
      logical :: um = .false., um_explicit = .false., um_implicit = .false., &
            um_radiation = .false., um_hydrology = .false.
      logical :: offline = .false., mk3l = .false.
   end type kbl_internal_switches 

   type (kbl_internal_switches), save :: cable_runtime

   !---CABLE runtime switches def in this type
   type kbl_user_switches
      character(len=200) :: VEG_PARS_FILE 
      character(len=20) :: DIAG_SOIL_RESP 
      character(len=200) :: LEAF_RESPIRATION
      character(len=200) :: FWSOIL_SWITCH 
      character(len=3) :: RUN_DIAG_LEVEL 
      logical :: INITIALIZE_MAPPING = .false. 
   end type kbl_user_switches

   type (kbl_user_switches), save :: cable_user

      type soilin_type
         real, dimension(:),allocatable :: silt
         real, dimension(:),allocatable :: clay
         real, dimension(:),allocatable :: sand
         real, dimension(:),allocatable :: swilt
         real, dimension(:),allocatable :: sfc
         real, dimension(:),allocatable :: ssat
         real, dimension(:),allocatable :: bch
         real, dimension(:),allocatable :: hyds
         real, dimension(:),allocatable :: sucs
         real, dimension(:),allocatable :: rhosoil
         real, dimension(:),allocatable :: css
         real, dimension(:),allocatable :: c3
      end type soilin_type

   type vegin_type
      real, dimension(:),allocatable :: canst1
      real, dimension(:),allocatable :: dleaf
      real, dimension(:),allocatable :: length 
      real, dimension(:),allocatable :: width
      real, dimension(:),allocatable :: vcmax
      real, dimension(:),allocatable :: ejmax
      real, dimension(:),allocatable :: hc
      real, dimension(:),allocatable :: xfang
      real, dimension(:),allocatable :: rp20
      real, dimension(:),allocatable :: rpcoef
      real, dimension(:),allocatable :: rs20
      real, dimension(:),allocatable :: wai 
      real, dimension(:),allocatable :: rootbeta 
      real, dimension(:),allocatable :: shelrb
      real, dimension(:),allocatable :: vegcf  !kdcorbin, 08/10
      real, dimension(:),allocatable :: frac4
      real, dimension(:,:),allocatable :: reflin
      real, dimension(:,:),allocatable :: taulin
      real, dimension(:),allocatable :: xalbnir
      real, dimension(:),allocatable :: extkn 
      real, dimension(:,:),allocatable :: froot
      real, dimension(:),allocatable :: tminvj
      real, dimension(:),allocatable :: tmaxvj
      real, dimension(:),allocatable :: vbeta
      real, dimension(:,:),allocatable :: cplant
      real, dimension(:,:),allocatable :: csoil
      real, dimension(:,:),allocatable :: ratecp
      real, dimension(:,:),allocatable :: ratecs
   end type vegin_type

   type(soilin_type) :: soilin
   type(vegin_type)  :: vegin

!   !---parameters, tolerances, etc. could be set in _directives.h
!   real, parameter :: RAD_TOLS = 1.0e-2

!jhan:temporary measure. improve hiding
!   real, dimension(:,:), pointer,save :: c1, rhoch
      
   contains


   SUBROUTINE get_type_parameters(logn,vegparmnew, classification)
      use define_dimensions
      use io_variables, only : filename, soil_desc, veg_desc
      implicit none
      ! Gets parameter values for each vegetation type and soil type.
      INTEGER(i_d),INTENT(IN) :: logn     ! log file unit number
      !opt. arg neccessary in offline CABLE IO
      character(len=4), intent(inout), optional :: classification
      LOGICAL,INTENT(IN)      :: vegparmnew ! new format input file (BP dec 2007)
      CHARACTER(LEN=80) :: comments 
      CHARACTER(LEN=10) :: vegtypetmp                   ! BP dec07
      CHARACTER(LEN=25) :: vegnametmp                   ! BP dec07
      !character(LEN=25), dimension(:),allocatable :: veg_desc 
      INTEGER(i_d) :: jveg                              ! BP dec07
      REAL(r_1)    :: notused                           ! BP dec07
      INTEGER(i_d) :: ioerror ! input error integer
      INTEGER(i_d) :: a ! do loop counter

         !================= Read in vegetation type specifications: ============
         OPEN(40,FILE=filename%veg,STATUS='old',ACTION='READ',IOSTAT=ioerror)
            IF(ioerror/=0) then 
               print *, 'CABLE_log: Cannot open veg type definitions.'
               stop
            endif
        
            IF (vegparmnew) THEN
               ! assume using IGBP/CSIRO vegetation types
               READ(40,*) comments
               READ(40,*) mvtype
               if( present(classification) ) &
                  WRITE(classification,'(a4)') comments(1:4)
            ELSE
               ! assume using CASA vegetation types
               !classification = 'CASA'
               READ(40,*)
               READ(40,*)
               READ(40,*) mvtype ! read # vegetation types
               READ(40,*)
               READ(40,*)
               comments = 'CASA'
            END IF
            
            WRITE(logn,'(A31,I3,1X,A10)') '  Number of vegetation types = ', mvtype, &
                 &  TRIM(comments)
        
            ! Allocate memory for type-specific vegetation parameters:
            ALLOCATE ( vegin%canst1(mvtype),vegin%dleaf(mvtype), &
                 vegin%vcmax(mvtype),vegin%ejmax(mvtype),vegin%hc(mvtype), &
                 vegin%xfang(mvtype),vegin%rp20(mvtype),vegin%rpcoef(mvtype), &
                 vegin%rs20(mvtype),vegin%shelrb(mvtype),vegin%frac4(mvtype), &
                 vegin%wai(mvtype),vegin%vegcf(mvtype),vegin%extkn(mvtype), &
                 vegin%tminvj(mvtype),vegin%tmaxvj(mvtype),vegin%vbeta(mvtype), &
                 vegin%rootbeta(mvtype),vegin%froot(ms,mvtype), &
                 vegin%cplant(ncp,mvtype),vegin%csoil(ncs,mvtype), &
                 vegin%ratecp(ncp,mvtype),vegin%ratecs(ncs,mvtype), &
                 vegin%xalbnir(mvtype),vegin%taulin(nrb,mvtype), &
                 vegin%reflin(nrb,mvtype), veg_desc(mvtype), & 
                 vegin%length(mvtype), vegin%width(mvtype) ) 
        
            IF (vegparmnew) THEN    ! added to read new format (BP dec 2007)
               
               ! Read in parameter values for each vegetation type:
               DO a = 1,mvtype 
                  READ(40,*) jveg, vegtypetmp, vegnametmp
                    
                     IF (jveg .GT. mvtype) then
                        write(6,*) 'jveg out of range in parameter file'
                        STOP 
                     endif   
                     veg_desc(jveg) = vegnametmp 
                  
                  READ(40,*) vegin%hc(jveg), vegin%xfang(jveg), vegin%width(jveg),  &
                           &   vegin%length(jveg), vegin%frac4(jveg)
                  READ(40,*) vegin%reflin(1:3,jveg) ! rhowood not used ! BP may2011
                  READ(40,*) vegin%taulin(1:3,jveg) ! tauwood not used ! BP may2011
                  READ(40,*) notused, notused, notused, vegin%xalbnir(jveg)
                  READ(40,*) notused, vegin%wai(jveg), vegin%canst1(jveg), &
                     vegin%shelrb(jveg), vegin%vegcf(jveg), vegin%extkn(jveg)
                  READ(40,*) vegin%vcmax(jveg), vegin%rp20(jveg), &
                             vegin%rpcoef(jveg), &
                             vegin%rs20(jveg)
                  READ(40,*) vegin%tminvj(jveg), vegin%tmaxvj(jveg), &
                             vegin%vbeta(jveg), vegin%rootbeta(jveg)
                  READ(40,*) vegin%cplant(1:3,jveg), vegin%csoil(1:2,jveg)
                  ! rates not currently set to vary with veg type
                  READ(40,*) vegin%ratecp(1:3,jveg), vegin%ratecs(1:2,jveg)
              END DO

            ELSE
   
               DO a = 1,mvtype 
                  READ(40,'(8X,A70)') veg_desc(a) ! Read description of each veg type
               END DO
   
               READ(40,*); READ(40,*) 
               READ(40,*) vegin%canst1
               READ(40,*) vegin%width
               READ(40,*) vegin%length
               READ(40,*) vegin%vcmax
               READ(40,*) vegin%hc
               READ(40,*) vegin%xfang
               READ(40,*) vegin%rp20
               READ(40,*) vegin%rpcoef
               READ(40,*) vegin%rs20
               READ(40,*) vegin%shelrb
               READ(40,*) vegin%frac4
               READ(40,*) vegin%wai
               READ(40,*) vegin%vegcf
               READ(40,*) vegin%extkn
               READ(40,*) vegin%tminvj
               READ(40,*) vegin%tmaxvj
               READ(40,*) vegin%vbeta
               READ(40,*) vegin%xalbnir
               READ(40,*) vegin%rootbeta
               READ(40,*) vegin%cplant(1,:)
               READ(40,*) vegin%cplant(2,:)
               READ(40,*) vegin%cplant(3,:)
               READ(40,*) vegin%csoil(1,:)
               READ(40,*) vegin%csoil(2,:)
               READ(40,*) 
               READ(40,*) vegin%ratecp(:,1)
               ! Set ratecp to be the same for all veg types:
               vegin%ratecp(1,:)=vegin%ratecp(1,1)
               vegin%ratecp(2,:)=vegin%ratecp(2,1)
               vegin%ratecp(3,:)=vegin%ratecp(3,1)
               READ(40,*) 
               READ(40,*) vegin%ratecs(:,1)
               vegin%ratecs(1,:)=vegin%ratecs(1,1)
               vegin%ratecs(2,:)=vegin%ratecs(2,1)
               ! old table does not have taul and refl ! BP may2011
               vegin%taulin(1,:) = 0.07
               vegin%taulin(2,:) = 0.425
               vegin%taulin(3,:) = 0.0
               vegin%reflin(1,:) = 0.07
               vegin%reflin(2,:) = 0.425
               vegin%reflin(3,:) = 0.0
   
            ENDIF
   
            write(6,*)'CABLE_log:Closing veg params file: ',trim(filename%veg)
         CLOSE(40)
         
         ! new calculation dleaf sin ce April 2012 
         vegin%dleaf = SQRT(vegin%width * vegin%length)
          
        !================= Read in soil type specifications: ============
         OPEN(40,FILE=filename%soil,STATUS='old',ACTION='READ',IOSTAT=ioerror)
            IF(ioerror/=0) then 
              print *, 'CABLE_log: Cannot open soil type definitions.'
              stop
            endif
   
            READ(40,*); READ(40,*)
            READ(40,*) mstype ! Number of soil types
            READ(40,*); READ(40,*)
        
            ALLOCATE ( soil_desc(mstype) )
            ALLOCATE ( soilin%silt(mstype), soilin%clay(mstype), soilin%sand(mstype) )
            ALLOCATE ( soilin%swilt(mstype), soilin%sfc(mstype), soilin%ssat(mstype) )
            ALLOCATE ( soilin%bch(mstype), soilin%hyds(mstype), soilin%sucs(mstype) )
            ALLOCATE ( soilin%rhosoil(mstype), soilin%css(mstype) )
        
            DO a = 1,mstype 
               READ(40,'(8X,A70)') soil_desc(a) ! Read description of each soil type
            END DO
   
            READ(40,*); READ(40,*) 
            READ(40,*) soilin%silt
            READ(40,*) soilin%clay
            READ(40,*) soilin%sand
            READ(40,*) soilin%swilt
            READ(40,*) soilin%sfc
            READ(40,*) soilin%ssat
            READ(40,*) soilin%bch
            READ(40,*) soilin%hyds
            READ(40,*) soilin%sucs
            READ(40,*) soilin%rhosoil
            READ(40,*) soilin%css
      CLOSE(40)

      return
   END SUBROUTINE get_type_parameters



end module cable_common_module

