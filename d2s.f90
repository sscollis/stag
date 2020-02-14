program main

!     Purpose:  Convert a double precision grid to single precision
!
!     Author:   S. Scott Collis
!
!     Date:     12-21-98
!
!     Compile:  f90 -O2 d2s.f90 -o d2s
!

      use const
      implicit none

      integer :: nxi, neta, nz, i, j
      real*4, allocatable :: x4(:,:), y4(:,:)
      real*8, allocatable :: x8(:,:), y8(:,:)

!.... read single precision grid file

      open(10,file='grid.dbl',form='unformatted')
      read(10) nxi, neta, nz
      write(*,"('Read grid (',i3,',',i3')')") nxi, neta
      allocate( x4(nxi,neta), y4(nxi,neta), x8(nxi,neta), y8(nxi,neta))
      read(10) (( x8(i,j), i=1,nxi), j=1,neta), &
               (( y8(i,j), i=1,nxi), j=1,neta)
      close(10)
 
      x4 = x8; y4 = y8

!.... read single precision grid file

      open(10,file='grid.sgl',form='unformatted')
      write(10) nxi, neta, nz
      write(10) (( x4(i,j), i=1,nxi), j=1,neta), &
                (( y4(i,j), i=1,nxi), j=1,neta), &
                ((   1.0e0, i=1,nxi), j=1,neta)
      close(10)

      stop
end program main

