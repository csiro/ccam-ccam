! Conformal Cubic Atmospheric Model
    
! Copyright 2016 Commonwealth Scientific Industrial Research Organisation (CSIRO)
    
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
    
module filnames_m

implicit none

private
public albfile, icefile, sstfile, topofile, zofile, rsmfile
public soilfile, vegfile, restfile, surfile, surf_00, surf_12
public ifile, ofile, so4tfile, eigenv, radfile, o3file
public mesonest, laifile, albnirfile, urbanfile, bathfile
public vegprev, vegnext, vegnext2, cnsdir, salfile, oxidantfile
public phenfile, casafile

character(len=160), save :: albfile = ' '
character(len=160), save :: icefile = ' '
character(len=160), save :: sstfile = ' '
character(len=160), save :: topofile = ' '
character(len=160), save :: zofile = ' '
character(len=160), save :: rsmfile = ' '
character(len=160), save :: soilfile = ' '
character(len=160), save :: vegfile = ' '
character(len=160), save :: restfile = ' '
character(len=160), save :: surfile = ' '
character(len=160), save :: surf_00 = 's_00a'
character(len=160), save :: surf_12 = 's_12a'
character(len=160), save :: ifile = ' '
character(len=160), save :: ofile = ' '
character(len=160), save :: so4tfile = ' '
character(len=160), save :: eigenv = ' '
character(len=160), save :: radfile = ' '
character(len=160), save :: o3file = ' '
character(len=160), save :: mesonest = ' '
character(len=160), save :: laifile = ' '
character(len=160), save :: albnirfile = ' '
character(len=160), save :: urbanfile = ' '
character(len=160), save :: bathfile = ' '
character(len=160), save :: vegprev = ' '
character(len=160), save :: vegnext = ' '
character(len=160), save :: vegnext2 = ' '
character(len=160), save :: cnsdir = ' '
character(len=160), save :: salfile = ' '
character(len=160), save :: oxidantfile = ' '
character(len=160), save :: phenfile = ' '
character(len=160), save :: casafile = ' '

end module filnames_m