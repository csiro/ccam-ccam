! Conformal Cubic Atmospheric Model
    
! Copyright 2015-2024 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module getopt_m

   implicit none
   private


!!!  F90 version of getopt, modified from GNU glibc-2.2.2 getopt.c
!!!  Comments from original C code have just a leading !. New comments
!!!  have !!!

!!!  Translated by Martin.Dix@csiro.au.

!  This is also the maximum length of the command line and the optstring
   integer, parameter, public :: MAX_ARGLEN=256

!  For communication from `getopt' to the caller.
!  When `getopt' finds an option that takes an argument,
!  the argument value is returned here.
!   character(len=MAX_ARGLEN), public :: optarg

!  Index in ARGV of the next element to be scanned.
!  This is used for communication to and from the caller
!  and for communication between successive calls to `getopt'.
!
!  On entry to `getopt', zero means this is the first call; initialize.
!
!  When `getopt' returns -1, this is the index of the first of the
!  non-option elements that the caller should itself scan.
!

!  Formerly, initialization of getopt depended on optind==0, which
!  causes problems with re-calling getopt as programs generally don't
!  know that.



!  Callers store zero here to inhibit the error message
!  for unrecognized options.

   logical, public :: opterr = .true.

!  Set to an option character which was unrecognized.
!  This must be initialized on some systems to avoid linking in the
!  system's own getopt implementation.

   integer, public :: optopt = ichar("?")

   public :: getopt, getcline


!    Scan elements of ARGV (whose length is ARGC) for option characters
!    given in OPTSTRING.

!    If an element of ARGV starts with '-', and is not exactly "-" or "--",
!    then it is an option element.  The characters of this element
!    (aside from the initial '-') are option characters.  If `getopt'
!    is called repeatedly, it returns successively each of the option characters
!    from each of the option elements.

!    If `getopt' finds another option character, it returns that character,
!    updating `optind' and `nextchar' so that the next call to `getopt' can
!    resume the scan with the following option character or ARGV-element.

!    If there are no more option characters, `getopt' returns -1.
!    Then `optind' is the index in ARGV of the first ARGV-element
!    that is not an option.  (The ARGV-elements have been permuted
!    so that those that are not options now come last.)

!    OPTSTRING is a string containing the legitimate option characters.
!    If an option character is seen that is not listed in OPTSTRING,
!    return '?' after printing an error message.  If you set `opterr' to
!    zero, the error message is suppressed but we still return '?'.

!    If a char in OPTSTRING is followed by a colon, that means it wants an arg,
!    so the following text in the same ARGV-element, or the text of the following
!    ARGV-element, is returned in `optarg'.  Two colons mean an option that
!    wants an optional arg; if there is text in the current ARGV-element,
!    it is returned in `optarg', otherwise `optarg' is set to zero.

!    If OPTSTRING starts with `-' or `+', it requests different methods of
!    handling the non-option ARGV-elements.
!    See the comments about RETURN_IN_ORDER and REQUIRE_ORDER, above.

!    Long-named options begin with `--' instead of `-'.
!    Their names may be abbreviated as long as the abbreviation is unique
!    or is an exact match for some defined option.  If they have an
!    argument, it follows the option name in the same ARGV-element, separated
!    from the option name by a `=', or else the in next ARGV-element.
!    When `getopt' finds a long-named option, it returns 0 if that option's
!    `flag' field is nonzero, the value of the option's `val' field
!    if the `flag' field is zero.

!    The elements of ARGV aren't really const, because we permute them.
!    But we pretend they're const in the prototype to be compatible
!    with other systems.

contains

!   subroutine getopt (optstring, nopt, opt)
   subroutine getopt (optstring, optind, opt, optarg )
      character(len=*), intent(in) :: optstring
      integer, intent(out) :: optind
      integer, intent(out) :: opt
      character(len=*), intent(out) :: optarg

!  The next char to be scanned in the option-element
!  in which the last option character we returned was found.
!  This allows us to pick up the scan where we left off.
!
!  If this is zero, or a null string, it means resume the scan
!  by advancing to the next ARGV-element.
      integer, save :: nextchar=0

      logical, save :: getopt_initialized = .false.
      logical :: print_errors
      integer, save :: argc, i
      character(len=MAX_ARGLEN), dimension(:), allocatable, save :: argv
      character(len=MAX_ARGLEN) :: optname
      character(len=MAX_ARGLEN), save :: nextstr = ""
      character(len=1) :: c
      integer :: temp, temp_p1, temp_p2
      integer :: stderr = 6

      print_errors = opterr
      if (optstring(1:1) == ":") then
         print_errors = .false.
      end if

      if (command_argument_count() < 1) then
         opt = -1
         optind = 1
         return
      end if

      optarg = ""

      if ( .not. getopt_initialized ) then
         optind = 1
         !!! Need to use iargc()+1 to get the same result as with C
         argc = command_argument_count()+1
         allocate ( argv(0:argc-1) )
         do i=0,argc-1
            call get_command_argument(i,argv(i)) 
         end do
         getopt_initialized = .true.
      end if

!!!   nextchar is an index into nextstring.
!!!   Need to use max(1,nextchar) to avoid bounds error if nextchar=0
!!!   Short circuit of if isn't guaranteed.
      if (nextchar == 0 .or. len_trim(nextstr(max(1,nextchar):)) == 0 ) then

         ! Advance to the next ARGV-element.

         ! The special ARGV-element `--' means premature end of options.
	 ! Skip it like a null option,
	 ! then exchange with previous non-options as if it were an option,
	 ! then skip everything else like a non-option.

!!!      Split if tests to avoid bounds errors when optind=argc
         if (optind /= argc ) then
            if ( argv(optind) == "--") then
               optind = optind+1
            end if
         end if

      ! If we have done all the ARGV-elements, stop the scan
      ! and back over any non-options that we skipped and permuted.

         if (optind == argc) then
!!$	  /* Set the next-arg-index to point at the non-options
!!$	     that we previously skipped, so the caller will digest them.  */
            opt = -1
            return
         end if

      !  If we have come to a non-option and did not permute it,
      !  either stop the scan or describe it to the caller and pass it by.

         if ( argv(optind)(1:1) /= "-" ) then
            opt = -1
            return
         end if

      !  We have found another option-ARGV-element.
      !  Skip the initial punctuation.

         nextchar = 2
         nextstr = argv(optind)
      end if


  ! Decode the current option-ARGV-element.

  !  Check whether the ARGV-element is a long option.

  !  If long_only and the ARGV-element has the form "-f", where f is
  !  a valid short option, don't consider it an abbreviated form of
  !  a long option that starts with f.  Otherwise there would be no
  !  way to give the -f short option.

  !  On the other hand, if there's a long option "fubar" and
  !  the ARGV-element is "-fu", do consider that an abbreviation of
  !  the long option, just like "--fu", and not "-f" with arg "u".

  !  This distinction seems to be the most useful approach.

  ! Look at and handle the next short option-character.
      c = nextstr(nextchar:nextchar)
      nextchar = nextchar+1
      temp = index(optstring, c)

      ! Increment `optind' when we start to process its last character.
      if ( len_trim(nextstr(nextchar:)) == 0 ) then
         optind = optind+1
      end if

      if (temp == 0 .or. c == ":") then
         if (print_errors) then
            write(unit=stderr,fmt="(a,a,a)") trim(argv(0)), ": invalid option -- ", c
         end if
         optopt = ichar(c)
         opt = ichar("?")
         return
      end if

!   Does optstring have something appended to ensure this isn't off the end???
      temp_p1 = min( temp+1, len(optstring) )
      if (optstring(temp_p1:temp_p1) == ":" .and. temp_p1 == temp+1) then
         temp_p2 = min( temp+2, len(optstring) )
         if (optstring(temp_p2:temp_p2) == ":" .and. temp_p2 == temp+2) then
            ! This is an option that accepts an argument optionally.
            if (len_trim(nextstr(nextchar:)) /= 0 ) then
               optarg = trim(nextstr(nextchar:))
               optind = optind + 1
            else
               optarg = ""
               nextchar = 0
            end if
         else
            ! This is an option that requires an argument.
            if (len_trim(nextstr(nextchar:)) /= 0 ) then
               optarg = trim(nextstr(nextchar:))
               ! If we end this ARGV-element by taking the rest as an arg,
               ! we must advance to the next element now.
               optind = optind + 1
            else if (optind == argc) then
               if (print_errors) then
                  write(unit=stderr,fmt="(a,a,a)") trim(argv(0)), ": option requires an argument -- ", c
               end if
               optopt = ichar(c)
               if (optstring(1:1) == ":") then
                  c = ":"
               else
                  c = "?"
               end if
            else
               ! We already incremented `optind' once;
               ! increment it again when taking next ARGV-elt as argument.
               optarg = argv(optind)
               optind = optind + 1
            end if
            nextchar = 0  ! Where does this belong, perhaps on next line?
         end if
      end if
      opt = ichar(c)
      return
   end subroutine getopt

   subroutine getcline ( cline )

!     Get the complete program command line
      character(len=*), intent(out) :: cline
      integer :: iarg
      character(len=MAX_ARGLEN) :: arg

      cline = ''
      do iarg=0,command_argument_count()
         !call getarg(iarg,arg)
         call get_command_argument(iarg,arg)
!        Use >= here to allow for the extra space
         if ( len_trim(cline) + len_trim(arg) >= len(cline) ) then
            print*, "Error, increase length of command line variable"
            stop
         end if
         cline = cline(1:len_trim(cline)) // " " // trim(arg)
      end do

      !  The loop above adds a leading blank so adjustl
      cline = adjustl(cline)

   end subroutine getcline

end module getopt_m

