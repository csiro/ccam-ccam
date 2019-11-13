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
public phenfile, casafile, casapftfile
public ensembleoutfile
public solarfile, ch4file, n2ofile
public cfc11file, cfc12file, cfc113file, hcfc22file

character(len=1024), save :: albfile = ' '
character(len=1024), save :: icefile = ' '
character(len=1024), save :: sstfile = ' '
character(len=1024), save :: topofile = ' '
character(len=1024), save :: zofile = ' '
character(len=1024), save :: rsmfile = ' '
character(len=1024), save :: soilfile = ' '
character(len=1024), save :: vegfile = ' '
character(len=1024), save :: restfile = ' '
character(len=1024), save :: surfile = ' '
character(len=1024), save :: surf_00 = 's_00a'
character(len=1024), save :: surf_12 = 's_12a'
character(len=1024), save :: ifile = ' '
character(len=1024), save :: ofile = ' '
character(len=1024), save :: so4tfile = ' '
character(len=1024), save :: eigenv = ' '
character(len=1024), save :: radfile = ' '
character(len=1024), save :: o3file = ' '
character(len=1024), save :: mesonest = ' '
character(len=1024), save :: laifile = ' '
character(len=1024), save :: albnirfile = ' '
character(len=1024), save :: urbanfile = ' '
character(len=1024), save :: bathfile = ' '
character(len=1024), save :: vegprev = ' '
character(len=1024), save :: vegnext = ' '
character(len=1024), save :: vegnext2 = ' '
character(len=1024), save :: cnsdir = ' '
character(len=1024), save :: salfile = ' '
character(len=1024), save :: oxidantfile = ' '
character(len=1024), save :: phenfile = ' '
character(len=1024), save :: casafile = ' '
character(len=1024), save :: ensembleoutfile = ' '
character(len=1024), save :: casapftfile = ' '
character(len=1024), save :: solarfile = ' '
character(len=1024), save :: ch4file = ' '
character(len=1024), save :: n2ofile = ' '
character(len=1024), save :: cfc11file = ' '
character(len=1024), save :: cfc12file = ' '
character(len=1024), save :: cfc113file = ' '
character(len=1024), save :: hcfc22file = ' '

end module filnames_m
