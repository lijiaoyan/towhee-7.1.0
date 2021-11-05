c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004 Marcus G. Martin,                           *
c     * see the file license.gpl for the full license information      *
c     *                                                                *
c     * This program is free software; you can redistribute it and/or  *
c     * modify it under the terms of the GNU General Public License    *
c     * as published by the Free Software Foundation; either version 2 *
c     * of the License, or (at your option) any later version.         *
c     *                                                                *
c     * This program is distributed in the hope that it will be useful,*
c     * but WITHOUT ANY WARRANTY; without even the implied warranty of *
c     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
c     * GNU General Public License for more details.                   *
c     *                                                                *
c     * You should have received a copy of the GNU General Public      *
c     * License along with this program; if not, write to the Free     *
c     * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
c     * MA  02111-1307, USA.                                           *
c     ******************************************************************
      program xyz2towhee
c     ******************************************************************
c     * converts from xyz to towhee format                             *
c     * this file is rewritten to convert various files into Towhee    *
c     *                                                                *
c     * originally written 03-27-2003 by M.G. Martin                   *
c     * last modified 07-13-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none

      character*2 suffix
      character*30 filename
      integer maxtype
      parameter (maxtype=10)
      integer natoms,idim,iatom,idum,base,typecount,itype
      dimension typecount(maxtype)
      double precision boxdims,xxx,yyy,zzz,distmult
      dimension boxdims(3,2)

c     --- filename passed in when running the code
      write(6,*) 'Please enter a filename'
      read(5,*) filename
      open(unit=10,file=filename,form="formatted")
      write(6,*) 'Enter a distance multiplier (1.0 for no change)'
      read(5,*) distmult

      do itype = 1,maxtype
         typecount(itype) = 0
      enddo

c     --- get the number of atoms
      read(10,*)
      read(10,*) natoms
c     --- get the box dimensions
      read(10,*)
      do idim = 1,3
         read(10,*) boxdims(idim,1),boxdims(idim,2)
      enddo
      write(6,*) 'Total box lengths'
      do idim = 1,3
         write(6,*) 'dimension, length',idim,distmult*(boxdims(idim,2)
     &        -boxdims(idim,1))
      enddo

c     --- get the atom positions and types
      read(10,*)
      base = 40
      do iatom = 1,natoms
         read(10,*) idum,itype,xxx,yyy,zzz
         if ( itype .gt. maxtype ) then
            write(6,*) 'itype exceeds maxtype'
            write(6,*) 'itype,maxtype',itype,maxtype
            stop
         endif
c        --- see if we need to open a file
         if ( typecount(itype) .eq. 0 ) then
            filename = 'output.'
            write(suffix,'(i2.2)') itype
            filename(8:9) = suffix
            open(unit=base+itype,file=filename,form="formatted")
         endif
         typecount(itype) = typecount(itype) + 1
         write(base+itype,*) distmult*(xxx-boxdims(1,1))
     &        ,distmult*(yyy-boxdims(2,1))
     &        ,distmult*(zzz-boxdims(3,1))
      enddo
      close(10)
      do itype = 1,maxtype
         if ( typecount(itype) .ne. 0 ) then
            write(6,*)'type, typecount',itype,typecount(itype)
            close(base+itype)
         endif
      enddo

      end
