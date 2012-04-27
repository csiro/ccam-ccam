!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++ USE this module in any subr. you wish to write vars from.             +++!
!+++ can write up to three vars of the form foo(x), per call               +++!
!+++ where x is typically the number of landpoints(tiles). binary file is  +++!
!+++ then appended every timestep with the new foo(x)                      +++!
!+++                                                                       +++! 
!+++ CALL syntax:                                                          +++!  
!+++                                                                       +++! 
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, var1 )     +++!
!+++                                                                       +++! 
!+++    OR                                                                 +++! 
!+++                                                                       +++! 
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, vname2, &  +++!
!+++              var1, var2 )                                             +++!      
!+++                                                                       +++! 
!+++    OR                                                                 +++! 
!+++                                                                       +++! 
!+++ cable_diag( Nvars, filename, dimx, dimy, timestep, vname1, vname2, &  +++!
!+++              vname3, var1, var2, var3 )                               +++!
!+++                                                                       +++! 
!+++ where:    Nvars       number of vars being written in call (1,2 or 3) +++!
!+++           filename    base of preferred filename.dat, filename.bin    +++!
!+++           dimx        length of x-dimension (# landpoints)            +++!      
!+++           dimy        length of y-dimension (# time steps )           +++!     
!+++           timestep    # of the timestep                               +++!     
!+++           vnameX      preferred/recognizable  name of var             +++!
!+++           varX        vsr to output                                   +++!   
!+++                                                                       +++! 
!+++ to plot the temp. of the first three soil layers in subr. bar()       +++!
!+++                                                                       +++! 
!+++ e.g.   in subr. bar()                                                 +++!
!+++     .                                                                 +++! 
!+++     .                                                                 +++! 
!+++      use cable_diag_mod                                               +++! 
!+++     .                                                                 +++! 
!+++     .                                                                 +++! 
!+++      call cable_diag( 3, 'bar_T', 2950, 48, ktau, 'tsoil1',           +++! 
!+++                'tsoil2','tsoil3',tsoil(:,1), tsoil(:,2), tsoil(:,3))  +++!
!+++     .                                                                 +++! 
!+++     .                                                                 +++!
!+++ following the run binaries can be interpreted from the command line   +++!
!+++ using the cable_diag.pl command (see cable_diag.pl for details)       +++!  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


module cable_diag_module
   use define_dimensions, only : i_d, r_1, r_2
   implicit none
   integer(i_d), parameter :: gok=0
   integer(i_d) :: gopenstatus,galloctest=1
  
   !--- subrs overloaded to respond to call cable_diag 
   interface cable_diag
      module procedure cable_diag1,cable_diag1i, cable_diag2, cable_diag3
   end interface cable_diag
  
   contains

   !=============================================================================!
   !--- cable generic print status
   !=============================================================================!

  subroutine cable_stat( routname)
      use cable_common_module, only : ktau_gl, knode_gl
      implicit none
      character(len=*), intent(in) :: routname
         if(knode_gl==1) & 
            write(6,*) 'CABLE@  ', routname, ktau_gl
      return 
   end subroutine cable_stat


   !=============================================================================!
   !--- cable_diag1/2/3 call subrs to write filename.dat which contains description
   !--- of data and format etc., and filename.bin containing the data   
   !=============================================================================!

   subroutine cable_diag1( Nvars, basename, dimx, dimy, timestep, node, vname1, var1 )
      implicit none
      integer(i_d), intent(in) :: Nvars,dimx, dimy, timestep,node
      real(r_1), intent(in), dimension(:) :: var1
      integer :: i=0
      character(len=*), intent(in) :: basename, vname1
      character(len=30) :: filename, chnode
         gopenstatus=1
         write(chnode,10) node
      10 format(i2.2)   
         filename=trim(trim(basename)//trim(chnode))
   
         call cable_diag_desc1( Nvars, trim(filename), dimx, dimy, vname1 )
         gopenstatus=1
         call cable_diag_data1( Nvars, trim(filename), dimx, timestep, dimy, var1 )
         i=i+1
      return 
   end subroutine cable_diag1

   subroutine cable_diag1i( Nvars, basename, dimx, dimy, timestep, node, vname1, var1 )
      implicit none
      integer(i_d), intent(in) :: Nvars,dimx, dimy, timestep,node
      integer, intent(in), dimension(:) :: var1
      integer :: i=0
      character(len=*), intent(in) :: basename, vname1
      character(len=30) :: filename, chnode
         gopenstatus=1
         write(chnode,10) node
      10 format(i2.2)   
         filename=trim(trim(basename)//trim(chnode))
         call cable_diag_desc1( Nvars, trim(filename), dimx, dimy, vname1 )
         gopenstatus=1
         call cable_diag_data1i( Nvars, trim(filename), dimx, timestep, dimy, var1 )
         i=i+1
      return 
   end subroutine cable_diag1i


   subroutine cable_diag2( Nvars, basename, dimx, dimy, timestep,node, vname1, vname2, var1, var2 )
      implicit none
      integer(i_d), intent(in) :: Nvars,dimx, dimy, timestep,node
      real(r_1), intent(in), dimension(:) :: var1, var2
      character(len=*), intent(in) :: basename, vname1, vname2
      character(len=30) :: filename, chnode
         gopenstatus=1
         write(chnode,10) node
      10 format(i2.2)   
         filename=trim(trim(basename)//trim(chnode))
         call cable_diag_desc2( Nvars, trim(filename), dimx, dimy, vname1, vname2 )
         gopenstatus=1
         call cable_diag_data2( Nvars, trim(filename), dimx, timestep, dimy, var1, var2 )
      return 
   end subroutine cable_diag2
   
   subroutine cable_diag3( Nvars, basename, dimx, dimy, timestep,node, vname1, vname2, vname3, &
                                             var1, var2, var3 )
      implicit none
      integer(i_d), intent(in) :: Nvars,dimx, dimy, timestep,node
      real(r_1), intent(in), dimension(:) :: var1, var2, var3
      character(len=*), intent(in) :: basename, vname1, vname2, vname3
      character(len=30) :: filename, chnode
         gopenstatus=1
         write(chnode,10) node
      10 format(i2.2)   
         filename=trim(trim(basename)//trim(chnode))
         call cable_diag_desc3( Nvars, trim(filename), dimx, dimy, vname1, vname2, vname3 )
         gopenstatus=1
         call cable_diag_data3( Nvars, trim(filename), dimx, timestep, dimy, var1, var2, var3 )
      return 
   end subroutine cable_diag3

!=============================================================================!
!=============================================================================!

   subroutine cable_diag_desc1( Nvars, filename, dimx, dimy, vname1 )
      implicit none
      integer(i_d), intent(in) :: Nvars,dimx,dimy 
      character(len=*), intent(in) :: filename, vname1
         open(unit=713941,file=filename//'.dat', status="replace",action="write", iostat=gopenstatus )
         if(gopenstatus==gok) then
               write (713941,*) 'Number of var(s): '
               write (713941,*) Nvars
               write (713941,*) 'Name of var(s): '
               write (713941,7139) vname1 
   7139        format(a)            
               write (713941,*) 'dimension of var(s) in x: '
               write (713941,*) dimx 
               write (713941,*) 'dimension of var(s) in y: '
               write (713941,*) dimy 
         else
            write (*,*), filename//'.dat',' Error: unable to write'
         endif
         close(713941)
      return 
   end subroutine cable_diag_desc1

   subroutine cable_diag_desc2( Nvars, filename, dimx, dimy, vname1, vname2 )
      implicit none
      integer(i_d), intent(in) :: Nvars,dimx,dimy 
      character(len=*), intent(in) :: filename, vname1, vname2
         open(unit=713941,file=filename//'.dat', status="replace",action="write", iostat=gopenstatus )
         if(gopenstatus==gok) then
               write (713941,*) 'Number of var(s): '
               write (713941,*) Nvars
               write (713941,*) 'Name of var(s): '
               write (713941,7139) vname1 
   7139        format(a)            
               write (713941,7139) vname2 
               write (713941,*) 'dimension of var(s) in x: '
               write (713941,*) dimx 
               write (713941,*) 'dimension of var(s) in y: '
               write (713941,*) dimy 
         else
            write (*,*) filename//'.dat',' Error: unable to write'
         endif
         close(713941)
      return 
   end subroutine cable_diag_desc2

   subroutine cable_diag_desc3( Nvars, filename, dimx, dimy, vname1, vname2, vname3 ) 
      implicit none
      integer(i_d), intent(in) :: Nvars,dimx,dimy 
      character(len=*), intent(in) :: filename, vname1, vname2, vname3 
         open(unit=713941,file=filename//'.dat', status="replace",action="write", iostat=gopenstatus )
         if(gopenstatus==gok) then
               write (713941,*) 'Number of var(s): '
               write (713941,*) Nvars
               write (713941,*) 'Name of var(s): '
               write (713941,7139) vname1 
   7139        format(a)            
               write (713941,7139) vname2 
               write (713941,7139) vname3 
               write (713941,*) 'dimension of var(s) in x: '
               write (713941,*) dimx 
               write (713941,*) 'dimension of var(s) in y: '
               write (713941,*) dimy 
         else
            write (*,*) filename//'.dat',' Error: unable to write'
         endif
         close(713941)
      return 
   end subroutine cable_diag_desc3

!=============================================================================!
!
   subroutine cable_diag_data1i( Nvars, filename, dimx, timestep, kend, var1  )
      implicit none
      integer(i_d), intent(in) :: Nvars, dimx, timestep, kend
      integer, intent(in), dimension(:) :: var1
      character(len=*), intent(in) :: filename
         open(unit=713942,file=filename//'.bin',status="unknown",action="write", &
            iostat=gopenstatus, form="unformatted", position='append' )
         if(gopenstatus==gok) then
               write (713942) var1
         else
            write (*,*) filename//'.bin',' NOT open for write. Error'
         endif
         close(713942)

      return 
   end subroutine cable_diag_data1i

!=============================================================================!

   subroutine cable_diag_data1( Nvars, filename, dimx, timestep, kend, var1  )
      implicit none
      integer(i_d), intent(in) :: Nvars, dimx, timestep, kend
      real(r_1), intent(in), dimension(:) :: var1
      character(len=*), intent(in) :: filename
         open(unit=713942,file=filename//'.bin',status="unknown",action="write", &
            iostat=gopenstatus, form="unformatted", position='append' )
         if(gopenstatus==gok) then
               write (713942) var1
         else
            write (*,*) filename//'.bin',' NOT open for write. Error'
         endif
         close(713942)

      return 
   end subroutine cable_diag_data1

   subroutine cable_diag_data2( Nvars, filename, dimx, timestep, kend, var1, var2  )
      implicit none
      integer(i_d), intent(in) :: Nvars, dimx, timestep, kend
      real(r_1), intent(in), dimension(:) :: var1, var2
      character(len=*), intent(in) :: filename
!      integer(i_d) :: frecl
         open(unit=713942,file=filename//'.bin',status="unknown",action="write", &
            iostat=gopenstatus, form="unformatted", position='append' )
         if(gopenstatus==gok) then
               write (713942) var1, var2
         else
            write (*,*) filename//'.bin',' NOT open for write. Error'
         endif
         close(713942)
#if 0
         frecl = Nvars * dimx* r_1 
         open(unit=713942,file=filename//'.bin',status="unknown",action="write", &
            iostat=gopenstatus, form="unformatted", access="direct", recl = frecl )
         if(gopenstatus==gok) then
               write (713942,rec=timestep) var1, var2
         else
            write (*,*) filename//'.bin',' NOT found'
         endif
         close(713942)
#endif
      return 
   end subroutine cable_diag_data2

   subroutine cable_diag_data3( Nvars, filename, dimx, timestep, kend, var1, var2, var3  )
      implicit none
      integer(i_d), intent(in) :: Nvars, dimx, timestep, kend
      real(r_1), intent(in), dimension(:) :: var1, var2, var3
      character(len=*), intent(in) :: filename
         open(unit=713942,file=filename//'.bin',status="unknown",action="write", &
            iostat=gopenstatus, form="unformatted", position='append' )
         if(gopenstatus==gok) then
               write (713942) var1, var2, var3
         else
            write (*,*) filename//'.bin',' NOT open for write. Error'
         endif
         close(713942)
      return 
   end subroutine cable_diag_data3

   subroutine cable_diag_data2r2( Nvars, filename, dimx, timestep, kend, var1, var2  )
      implicit none
      integer(i_d), intent(in) :: Nvars, dimx, timestep, kend
      real(r_2), intent(in), dimension(:) :: var1, var2
      character(len=*), intent(in) :: filename
         open(unit=713942,file=filename//'.bin',status="unknown",action="write", &
            iostat=gopenstatus, form="unformatted", position='append' )
         if(gopenstatus==gok) then
               write (713942) var1, var2
         else
            write (*,*) filename//'.bin',' NOT open for write. Error'
         endif
         close(713942)
      return 
   end subroutine cable_diag_data2r2

   subroutine cable_diag_data3r2( Nvars, filename, dimx, timestep, kend, var1, var2, var3  )
      implicit none
      integer(i_d), intent(in) :: Nvars, dimx, timestep, kend
      real(r_2), intent(in), dimension(:) :: var1, var2, var3
      character(len=*), intent(in) :: filename
         open(unit=713942,file=filename//'.bin',status="unknown",action="write", &
            iostat=gopenstatus, form="unformatted", position='append' )
         if(gopenstatus==gok) then
               write (713942) var1, var2, var3
         else
            write (*,*) filename//'.bin',' NOT open for write. Error'
         endif
         close(713942)
      return 
   end subroutine cable_diag_data3r2


end module cable_diag_module



