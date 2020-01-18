program mkgrid

! Simple two-dimensional grid generator
!
! S. Scott Collis
! Rice University
! MEMS, MS 321
! Houston, TX 77005
! 713-348-3617
! collis@rice.edu

implicit none

integer :: nx, ny, i, j
real, allocatable :: x(:), y(:), xi(:), eta(:)

!real :: xmin=-32, xmax=32, ymin=0, ymax=32
!real :: xmin=-20, xmax=20, ymin=0, ymax=30
real :: xmin=-20, xmax=20, ymin=0, ymax=40
real :: dx, dy, dxi, deta
real :: xs, ys, aa, bb

write(*,"('Enter number of cells (nx, ny) ==> ',$)") 
read(*,*) nx, ny
write(*,"('Enter Xs, Ys (Xs = Ys = 0 makes a uniform mesh) ==> ',$)")
read(*,*) xs, ys

allocate( x(nx+2), y(ny+2), xi(nx+2), eta(ny+2) )

dx = (xmax - xmin) / real(nx)
dy = (ymax - ymin) / real(ny)

if (xs.eq.0) then
  do i = 1, nx+1
    x(i) = xmin + dx * (i-1)
  end do
  x(nx+2) = x(nx+1)
else
  dxi = 1.0 / float(nx)
  write(*,"('Algebraic grid in x')")
  aa = (xmax-xmin) * xs / ( (xmax-xmin) - 2.0 * xs )
  bb = 1.0 + aa / (xmax-xmin)

  do i = 1, nx+1
    xi(i) = -0.5 + (i-1) * dxi
    x(i)  = 4.0 * ( aa * xi(i) / (bb - abs(xi(i))) )
  end do
  x(nx+2) = x(nx+1)
  do i = 1, nx+1
    write(11,"(2(1pe13.6,1x))") xi(i), x(i)
  end do
end if

if (ys.eq.0) then
  do j = 1, ny+1
    y(j) = ymin + dy * (j-1)
  end do
  y(ny+2) = y(ny+1)
else
  deta = 1.0 / float(ny)
  write(*,"('Algebraic grid in y')")
  aa = ymax * ys / ( ymax - 2.0 * ys )
  bb = 1.0 + aa / ymax
  do j = 1, ny+1
    eta(j) = (j-1) * deta
    y(j)   = aa * eta(j) / (bb - eta(j))
  end do
  y(ny+2) = y(ny+1)
  do j = 1, ny+1
    write(12,"(2(1pe13.6,1x))") eta(j), y(j)
  end do
end if

open(10,file='grid.in',form='formatted',status='unknown')
write(10,*) nx+2, ny+2
write(10,*) x
write(10,*) y
close(10)

stop
end program mkgrid
