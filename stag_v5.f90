!==============================================================================
!
! Program:  stag
!
! Purpose:  Solve the two-dimensional incompressible Navier-Stokes equations
!           on a staggered mesh system
!
! Author:   S. Scott Collis
!           Rice University
!           MEMS, MS 321
!           Houston, TX 77005
!           713-348-3617
!           collis@rice.edu
!
! Date:     12-23-99
!
! Version:  5.0 (Optimized for OpenMP)
!           Similar to version 4.0 but with REAL defined and without CG solver
!
!==============================================================================

#ifdef R8
#  define REAL real(kind=8)
#  define SGEFA dgefa
#  define SGESL sgesl
#endif

module global

!.... Module containing globally defined variables for the 
!     2-D Navier-Stokes solver

!.... spatial geometry
!
!     nx          := number of nodes in x-direction
!     ny          := number of nodes in y-direction
!     x,y         := nondimensional node locations
!     axis        := Axisymmetric flag

      integer :: nx, ny, nxm, nxmm, nym, nymm, i, j
      REAL, allocatable :: x(:), y(:), xc(:), yc(:), r(:), fx(:), fy(:)

!.... Code parameters
!
!     n       := current timestep index
!     nstep   := number of steps in this run
!     time    := current nondimensional time
!     dt      := current nondimensional timestep
!     dt_max  := maximum allowable timestep
!     cfl     := current CFL number
!     cfl_max := maximum allowable CFL number
!     tol     := convergence tolerance
!     method  := Discrete method
!     eps     := Fourth-order damping coefficient
!     ictype  := initial condition flag
!               0 = start using initial condition in ic() subroutine
!               1 = restart from file 'restart.dat'

      integer :: ictype, n, istep, nstep, method, k, kmax, pmethod, nout
      REAL    :: time, dt, cfl, dt_max, cfl_max, tol, eps, omega
      REAL    :: dtau, pnorm, pmax, psum
      logical :: delta_p

!.... field variables
!
!     u(nx,ny)      := x-velocity
!     v(nx,ny)      := y-velocity
!     p(nx,ny)      := presure
!     cx(nx,ny)     := mass flux out of the east face
!     cy(nx,ny)     := mass flux out of the north face

      REAL, allocatable :: u(:,:), v(:,:), p(:,:), cx(:,:), cy(:,:)
      REAL, allocatable :: u_old(:,:), v_old(:,:), p_old(:,:), ru(:,:), rv(:,:)
      REAL, allocatable :: uc(:,:),vc(:,:),pc(:,:),div(:,:),wz(:,:),ke(:,:)
      REAL :: tdiv, tke, enstrophy

!.... flow parameters
!
!     Re     := reference Reynolds number
!     rho    := fluid density
!     axis   := axisymmetry flag

      REAL :: Re, rho, vis
      logical :: axis

!.... units and filenames

      integer ::  ihist=9, io=10, ipres=11, istat=12
      character*80 :: base, filen
      integer :: iver, lfile=80
      integer, external :: igetver
      
end module global

module timemod

!.... Variables for monitoring CPU usage
!
!     Note that they must be explicitly declared as real*4

      real(kind=4), external :: second
      real(kind=4) :: cpu0, cpu, cpu2
end module timemod

module const

!.... Module defining some convenient constants to machine precision

      REAL, parameter :: zero    = 0.00000000000000000000
      REAL, parameter :: pt25    = 0.25000000000000000000
      REAL, parameter :: pt33    = 0.33333333333333333333
      REAL, parameter :: pt5     = 0.50000000000000000000
      REAL, parameter :: pt66    = 0.66666666666666666667
      REAL, parameter :: one     = 1.00000000000000000000
      REAL, parameter :: onept5  = 1.50000000000000000000
      REAL, parameter :: onept33 = 1.33333333333333333333
      REAL, parameter :: onept66 = 1.66666666666666666667
      REAL, parameter :: two     = 2.00000000000000000000
      REAL, parameter :: three   = 3.00000000000000000000
      REAL, parameter :: four    = 4.00000000000000000000
      REAL, parameter :: pi      = 3.14159265358979323846
end module const

!==============================================================================
program Navier_Stokes
!==============================================================================
!     Solve the incompressible, two-dimensional, non-dimensional 
!     Navier-Stokes equations on a staggered grid system
!
!     Author:  S. Collis
!
!     Date:    12-8-99
!
!     Compile: 
!               (SGI single precision) f90 -0 stag.f90 -o stag
!               (SGI double precision) f90 -r8 -DR8 -O stag.f90 -o stag
!
!     Variable Arrangement on Staggered Mesh:
!
!             u v
!            yc y                           Key
!                                        ---------
!             7 6  x+-+-+-+-+x            u   = -
!             6    -o-o-o-o-o-            v   = +
!               5  ++ + + + ++            u,v = x
!             5    -o-o-o-o-o-            p   = o
!               4  ++ + + + ++
!             4    -o-o-o-o-o-
!               3  ++ + + + ++
!             3    -o-o-o-o-o-
!               2  ++ + + + ++
!             2    -o-o-o-o-o-
!             1 1  x+-+-+-+-+x
!
!                  1 2 3 4 5 6   x , u
!                  12 3 4 5 67   xc, v
!==============================================================================
      use global
      implicit none

!$    integer, external :: omp_get_num_threads, omp_get_thread_num
!$    integer, external :: omp_get_num_procs

!     call omp_set_num_threads(nprocs)
!$omp parallel
!$    if (omp_get_thread_num() == 0) then
!$      write(*,*) 'Running on ',omp_get_num_threads(),' processor(s)'
!$      write(*,*) 'There are ',omp_get_num_procs(),' processor(s) available'
!$    end if
!$omp end parallel

      call setup
      call wgrid

      iver = igetver( base, 'h'//char(0) ) + 1
      call addver(base, 'h'//char(0), iver, filen, lfile)
      open(ihist,file=filen,status='unknown')

      iver = igetver( base, 'p'//char(0) ) + 1
      call addver(base, 'p'//char(0), iver, filen, lfile)
      open(ipres,file=filen,status='unknown')

      iver = igetver( base, 's'//char(0) ) + 1
      call addver(base, 's'//char(0), iver, filen, lfile)
      open(istat,file=filen,status='unknown')

      if (ictype.eq.0) then
        call ic
      else if (ictype.eq.1) then
        call restart
      else
        call error('Navier-Stokes$','Illegal value of ictype$')
      end if

      if (method .eq. 0) then
        call driver                   ! Runge-Kutta 2
      else
        call error('Navier-Stokes$','Illegal value of method$')
      end if

      stop
end program Navier_Stokes

subroutine driver

!.... Main driver for second order Runge-Kutta      

      use const
      use global
      use timemod
      implicit none

      REAL :: fact

      write(*,"(79('='))")

      call interp
      call setdt

      do istep = 1, nstep

!$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            u_old(i,j) = u(i,j)
            v_old(i,j) = v(i,j)
          end do
        end do

        if (delta_p) then 
!$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              p_old(i,j) = p(i,j)
              p(i,j) = zero
            end do
          end do
        endif
        
        n = n + 1

!.... first RK sub-step

        call calc_u( delta_p )
        call calc_v( delta_p )

        time = time + pt5 * dt
        fact = pt5 * dt
!$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            u(i,j) = u_old(i,j) + fact * ru(i,j)
            v(i,j) = v_old(i,j) + fact * rv(i,j)
          end do
        end do

        call calc_p( pt5 * dt, u, v, .true. )

!.... second RK sub-step

        call calc_u( delta_p )
        call calc_v( delta_p )

        time = time + pt5 * dt
        fact = dt
!$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            u(i,j) = u_old(i,j) + fact * ru(i,j)
            v(i,j) = v_old(i,j) + fact * rv(i,j)
          end do
        end do

        call calc_p ( dt, u, v, .true. )

        if (delta_p) then
!$omp parallel do private(i)
          do j = 1, ny
            do i = 1, nx
              p(i,j) = p_old(i,j) + p(i,j)
            end do
          end do
        endif

        cpu2 = second()
        write(*,"(i5,5(1x,1pe10.3))") n, time, dt, cfl, cpu2, cpu2-cpu
        write(ihist,"(i5,5(1x,1pe10.3))") n, time, dt, cfl, cpu2, cpu2-cpu
        cpu = cpu2

        call interp
        call setdt
        call stats

        call flush(ipres); call flush(istat); call flush(ihist)

!.... write restart file and Plot3d file

        if ( mod(n,nout) .eq. 0 .or. istep .eq. nstep) then
          iver = igetver( base, 'R'//char(0) ) + 1
          call addver(base, 'R'//char(0), iver, filen, lfile)
          write(*,"(79('='))")
          write (*,10) filen(1:index(filen,' ')-1), n
10        format('Writing restart file:  ',a,'  at time step ',i5)
          write(*,"(79('='))")
          open(io,file=filen,form='unformatted')
          write(io) n, nx, ny
          write(io) Re, time
          write(io) u, v, p, cx, cy
          close(io)
          iver = igetver( base, 'q'//char(0) ) + 1
          call addver(base, 'q'//char(0), iver, filen, lfile)
          call wdata( filen )
        end if

      end do

      close(ihist)

      return
end subroutine driver

subroutine calc_u( include_p )

!.... build the u-momentum equation over the u-control volumes

      use const
      use global
      implicit none

      logical :: include_p
      REAL :: dxwr, dxer, dysr, dynr, volr

      REAL :: ss, ds, cs, sn, dn, cn, se, de, ce, sw, dw, cw

!.... initialize quantities along the west boundary

      dxwr = one / (x(2)-x(1))
!$omp parallel do private(i,sn,ss,se,sw,dn,ds,de,dw,cn,cs,ce,cw) &
!$omp& private(dynr,dysr,dxer,volr)
      do j = 2, nym
        sw = pt5 * (r(j)+r(j-1)) * (y(j)-y(j-1))
        dw = vis * sw * dxwr
        cw = pt5 * ( cx(2,j) + cx(1,j) )

        se = pt5*(r(j)+r(j-1)) * (y(j)-y(j-1))
        dynr = one / ( yc(j+1) - yc(j) )
        dysr = one / ( yc(j) - yc(j-1) )

!.... main loop - north, south, and east sides

        do i = 2, nxmm

          ss = r(j-1) * (xc(i+1) - xc(i))
          sn = r(j) * (xc(i+1) - xc(i))

          dxer = one / ( x(i+1)-x(i) )

          cn = pt5 * ( cy(i,j) + cy(i+1,j) )
          cs = pt5 * ( cy(i,j-1) + cy(i+1,j-1) )
          ce = pt5 * ( cx(i,j) + cx(i+1,j) )

          dn = vis * sn * dynr
          ds = vis * ss * dysr
          de = vis * se * dxer

          volr = one / ((y(j)-y(j-1))*(xc(i+1)-xc(i)))

          ru(i,j) =  volr * ( &
                pt5 * ( (u(i,j)+u(i-1,j))*cw - (u(i+1,j)+u(i,j))*ce ) - &
                cn * ( u(i,j+1) * fy(j) + u(i,j) * (one-fy(j)) ) + &
                cs * ( u(i,j) * fy(j-1) + u(i,j-1) * (one-fy(j-1)) ) + &
                de * ( u(i+1,j)-u(i,j) ) - dw * ( u(i,j)-u(i-1,j) ) + &
                dn * ( u(i,j+1)-u(i,j) ) - ds * ( u(i,j)-u(i,j-1) ) )

          if (include_p) then
            ru(i,j) = ru(i,j) + volr * ( sw * p_old(i,j) - se * p_old(i+1,j) )
          end if

          sw = se
          dw = de
          cw = ce

        end do
      end do

end subroutine calc_u

subroutine calc_v( include_p )

!.... build the v-momentum equation over the v-control volumes

      use const
      use global
      implicit none

      logical :: include_p
      REAL :: dxwr, dxer, dysr, dynr, volr
      REAL :: ss, ds, cs, sn, dn, cn, se, de, ce, sw, dw, cw

!.... initialize quantities along the west boundary

      dxwr = one / (xc(2)-xc(1))
!$omp parallel do private(i,sn,ss,se,sw,dn,ds,de,dw,cn,cs,ce,cw) &
!$omp& private(dynr,dysr,dxer,volr)
      do j = 2, nymm
        sw = pt5 * (r(j+1)+r(j-1)) * (yc(j+1)-yc(j))
        dw = vis * sw * dxwr
        cw = pt5 * ( cx(1,j) + cx(1,j+1) )

        se = pt5 * (r(j+1)+r(j-1)) * (yc(j+1) - yc(j))
        dynr = one / ( y(j+1)-y(j) )
        dysr = one / ( y(j) - y(j-1) )

!.... main loop - north, south, and east sides

        do i = 2, nxm

          ss = pt5 * (r(j)+r(j-1)) * (x(i) - x(i-1))
          sn = pt5 * (r(j+1)+r(j)) * (x(i) - x(i-1))

          dxer = one / ( xc(i+1)-xc(i) )

          cn = pt5 * ( cy(i,j+1) + cy(i,j) )
          cs = pt5 * ( cy(i,j) + cy(i,j-1) )
          ce = pt5 * ( cx(i,j+1) + cx(i,j) )

          dn = vis * sn * dynr
          ds = vis * ss * dysr
          de = vis * se * dxer

          volr = one / ((x(i)-x(i-1))*(yc(j+1)-yc(j))) 

          rv(i,j) = volr * ( &
                pt5 * ( (v(i,j)+v(i,j-1))*cs - (v(i,j+1)+v(i,j))*cn ) - &
                ce * ( v(i+1,j) * fx(i) + v(i,j) * (one-fx(i)) ) + &
                cw * ( v(i,j) * fx(i-1) + v(i-1,j) * (one-fx(i-1)) ) + &
                de * ( v(i+1,j)-v(i,j) ) - dw * ( v(i,j)-v(i-1,j) ) + &
                dn * ( v(i,j+1)-v(i,j) ) - ds * ( v(i,j)-v(i,j-1) ) )

          if (include_p) then
            rv(i,j) = rv(i,j) + volr * (ss * p_old(i,j) - sn * p_old(i,j+1))
          end if

          sw = se
          dw = de
          cw = ce

        end do
      end do

end subroutine calc_v

subroutine calc_p( dtl, ul, vl, project )

!.... Find the pressure that makes the (ul, vl) field divergence free 
!.... by solving a Poisson equation

!.... If project is true, then project (ul, vl) to the divergence free
!.... space using the computed pressure

!.... Note that p here may be either (p - p_old) or the full
!.... pressure depending on the value of the delta_p flag

      use const
      use global
      implicit none

      REAL :: ul(nx,ny), vl(nx,ny)
      REAL :: dtl, dtlr
      logical :: project
      REAL :: sw, se, sn, ss, dw, de, dn, ds
      REAL :: dxwr, dxer, dysr, dynr
      REAL :: aw(nx,ny), ae(nx,ny), as(nx,ny), an(nx,ny), ap(nx,ny), rp(nx,ny)

!.... re-scale the pressure to be consistent

      if (dtl .ne. one) then
!$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            p(i,j) = p(i,j) * dtl
          end do
        end do
      end if

!$omp parallel do private(i)
      do j = 1, ny
        do i = 1, nx
          rp(i,j) = zero; aw(i,j) = zero; ae(i,j) = zero; 
          as(i,j) = zero; an(i,j) = zero; ap(i,j) = zero;
        end do
      end do

!.... initialize quantities along the south boundary

      dysr = one / ( yc(2) - yc(1) )
      do i = 2, nxm
        ss = r(1) * (x(i) - x(i-1))
        ds = ss * dysr
        cy(i,1) = rho * ss * vl(i,1)
      end do

!.... initialize quantities along the west boundary

      dxwr = one / (xc(2)-xc(1))
!$omp parallel do private(i,sn,ss,se,sw,dynr,dysr,dxer,dn,ds,de,dw)
      do j = 2, nym
        sw = pt5 * (r(j)+r(j-1)) * (y(j)-y(j-1))
        dw = sw * dxwr
        cx(1,j) = rho * sw * ul(1,j)

        se = pt5 * (r(j)+r(j-1)) * (y(j)-y(j-1))
        dynr = one / (yc(j+1) - yc(j))
        dysr = one / (yc(j) - yc(j-1))

!.... main loop - east and north sides

        do i = 2, nxm
        
          ss = r(j-1) * (x(i) - x(i-1))
          sn = r(j) * (x(i) - x(i-1))
          
          dxer = one / ( xc(i+1)-xc(i) )

          cx(i,j) = rho * se * ul(i,j)
          cy(i,j) = rho * sn * vl(i,j)

          dn = sn * dynr
          ds = ss * dysr
          de = se * dxer

          ae(i,j) = de
          aw(i,j) = dw
          an(i,j) = dn
          as(i,j) = ds

          dw = de

        end do
      end do

!.... form the RHS

      psum = zero
!$omp parallel do private(i) reduction(+:psum)
      do j = 2, nym
        do i = 2, nxm
          rp(i,j) = ( cx(i,j) + cy(i,j) - cx(i-1,j) - cy(i,j-1) )
          psum = psum + rp(i,j)**2
        end do
      end do

!.... satisfy wall boundary conditions (zero normal pressure gradient)
      
!$omp parallel do
      do j = 1, ny
        aw(2,j) = zero
        ae(nxm,j) = zero
      end do
!$omp parallel do
      do i = 1, ny
        as(i,2) = zero
        an(i,nym) = zero
      end do

!.... final assembly

!$omp parallel do private(i)
      do j = 2, nym
        do i = 2, nxm
          ap(i,j) = -(ae(i,j)+aw(i,j)+an(i,j)+as(i,j))
        end do
      end do

      if (pmethod .eq. 0) then
        call point_gs_sor( rp, aw, ae, as, an, ap )
      else if (pmethod .eq. 1) then
        call line_gs_sor( rp, aw, ae, as, an, ap )
      else if (pmethod .eq. 2) then
        call adi( rp, aw, ae, as, an, ap )
      else if (pmethod .eq. 3) then
        call direct( rp, aw, ae, as, an, ap )
      end if

      if (project) then

!.... correct u-velocity and Cx mass flux

!$omp parallel do private(i,se)
      do j = 2, nym
        do i = 2, nxmm
          se = pt5 * (r(j)+r(j-1)) * (y(j)-y(j-1))
          ul(i,j) = ul(i,j) - se * ( p(i+1,j) - p(i,j) ) / &
                   ( rho * (xc(i+1)-xc(i))*(y(j)-y(j-1)) )
          cx(i,j) = rho * se * ul(i,j)
        end do
      end do

!.... correct v-velocity and Cy mass flux

!$omp parallel do private(i,sn)
      do j = 2, nymm
        do i = 2, nxm
          sn = r(j) * ( x(i) - x(i-1) )
          vl(i,j) = vl(i,j) - sn * ( p(i,j+1) - p(i,j) ) / &
                   ( rho * (x(i)-x(i-1)) * (yc(j+1)-yc(j)) )
          cy(i,j) = rho * sn * vl(i,j)
        end do
      end do

      else

!$omp parallel do private(i,se)
      do j = 2, nym
        do i = 2, nxmm
          se = pt5 * (r(j)+r(j-1)) * (y(j)-y(j-1))
          cx(i,j) = rho * se * u(i,j)
        end do
      end do

!.... correct v-velocity and Cy mass flux

!$omp parallel do private(i,sn)
      do j = 2, nymm
        do i = 2, nxm
          sn = r(j) * ( x(i) - x(i-1) )
          cy(i,j) = rho * sn * v(i,j)
        end do
      end do

      end if

!.... re-scale the pressure to be consistent

      dtlr = one / dtl
      if (dtl .ne. one) then
!$omp parallel do private(i)
        do j = 1, ny
          do i = 1, nx
            p(i,j) = p(i,j) * dtlr
          end do
        end do
      end if

!.... compute pressure statistics

      psum = zero
      pmax = zero
      tdiv = zero
!$omp parallel do private(i) reduction(+: psum, tdiv) reduction(max: pmax)
      do j = 2, nym
        do i = 2, nxm
          psum = psum + p(i,j)**2
          pmax = max(abs(p(i,j)),pmax)
          tdiv = tdiv + (cx(i,j) + cy(i,j) - cx(i-1,j) - cy(i,j-1))**2
        end do
      end do
      psum = sqrt( psum / (nxmm*nymm) )
      tdiv = sqrt( tdiv / (nxmm*nymm) )

!.... write to history file

      write(*,"('----> ',4(1pe10.3,1x),'(',i5,')')") tdiv, pnorm, psum, pmax, k
      write(ipres,"(2(i5,1x),8(1pe10.3,1x))") n, k, tdiv, pnorm, psum, pmax

end subroutine calc_p

module direct_mod
      use const
      use global
      implicit none
      REAL, allocatable :: a(:,:)
contains
  subroutine init_direct
      allocate( a(nx*ny,nx*ny) )
  end subroutine init_direct
end module direct_mod

subroutine direct( rp, aw, ae, as, an, ap )

!.... Direct LU-factorization and solve of the Poisson equation for pressure

!.... Note that in general, one only has to factor the matrix operator once
!.... but this is currently done at every time step

      use const
      use global
      use direct_mod
      implicit none

      REAL :: rp(nx,ny), aw(nx,ny), ae(nx,ny), as(nx,ny), an(nx,ny), ap(nx,ny)
      integer :: row, ipvt(nx*ny), info

!.... hold a point fixed

      aw(2,2) = zero
      as(2,2) = zero
      ae(2,2) = zero
      an(2,2) = zero
      ap(2,2) = one
      rp(2,2) = zero

!$omp parallel do private(i)
      do j = 1, nx*ny
        do i = 1, nx*ny
          a(i,j)  = zero
        end do
      end do

!$omp parallel do private(i, row)
      do j = 2, nym
        do i = 2, nxm
          row = i + nx*(j-1)
          a(row,row-1)  = aw(i,j)
          a(row,row+1)  = ae(i,j)
          a(row,row)    = ap(i,j)
          a(row,row+nx) = an(i,j)
          a(row,row-nx) = as(i,j)
        end do
      end do

!$omp parallel do
      do i = 1, nx
        a(i,i) = one
        a(i+nx*nym,i+nx*nym) = one
      end do

!$omp parallel do
      do j = 1, ny
        a(1+nx*(j-1),1+nx*(j-1)) = one
        a(nx+nx*(j-1),nx+nx*(j-1)) = one
      end do

      call SGEFA( a, nx*ny, nx*ny, ipvt, info )
      if (info.ne.0) call error('direct$','Singular matrix$')
      call SGESL( A, nx*ny, nx*ny, ipvt, rp, 0)
   
!$omp parallel do private(i) 
      do j = 1, ny
        do i = 1, nx
          p(i,j) = rp(i,j)
        end do
      end do

!.... evaluate the residual norm

        pnorm = zero
!$omp parallel do private(i) reduction(+: pnorm)
        do j = 2, nym
          do i = 2, nxm
            pnorm = pnorm + ( ap(i,j) * p(i,j) + ae(i,j) * p(i+1,j) + &
                    an(i,j) * p(i,j+1) + aw(i,j) * p(i-1,j) + &
                    as(i,j) * p(i,j-1) - rp(i,j) )**2
          end do
        end do
        pnorm = sqrt(pnorm / (nxmm*nymm))

      return
end subroutine direct

subroutine point_gs_sor( rp, aw, ae, as, an, ap )

!.... Solve the pressure Poisson equation using point GS with SOR

      use const
      use global
      implicit none

      REAL :: rp(nx,ny), aw(nx,ny), ae(nx,ny), as(nx,ny), an(nx,ny), ap(nx,ny)
      integer :: pl

      do k = 1, kmax

!.... solve for pressure using point GS-SOR

!$omp parallel do private(i)
        do j = 2, nym
          do i = 2, nxm
            p(i,j) = (one-omega)*p(i,j) - omega * ( ae(i,j) * p(i+1,j) + &
                     an(i,j) * p(i,j+1) + aw(i,j) * p(i-1,j) + &
                     as(i,j) * p(i,j-1) - rp(i,j) ) / ap(i,j)
          end do
        end do

!.... ensure that boundary conditions are fixed

        p(1,1)         = p(2,2)
        p(nx,1)        = p(nxm,2)
        p(1,ny)        = p(2,nym)
        p(nx,ny)       = p(nxm,nym)

!$omp parallel do 
        do i = 2, nxm
          p(i,1)  = p(i,2)
          p(i,ny) = p(i,nym)
        end do

!$omp parallel do 
        do j = 2, nym
          p(1,j)  = p(2,j)
          p(nx,j) = p(nxm,j)
        end do

        pl = p(2,2)
!$omp parallel do private(i)
        do j = 2, nym
          do i = 2, nxm 
            p(i,j) = p(i,j) - pl  ! hold a point fixed
          end do
        end do

!.... evaluate the residual norm

        pnorm = zero
!$omp parallel do private(i) reduction(+: pnorm)
        do j = 2, nym
          do i = 2, nxm
            pnorm = pnorm + ( ap(i,j) * p(i,j) + ae(i,j) * p(i+1,j) + &
                    an(i,j) * p(i,j+1) + aw(i,j) * p(i-1,j) + &
                    as(i,j) * p(i,j-1) - rp(i,j) )**2
          end do
        end do
        pnorm = sqrt(pnorm / (nxmm*nymm))

        if (pnorm .lt. tol) exit

      end do

      return
end subroutine point_gs_sor

subroutine line_gs_sor( rp, aw, ae, as, an, ap )

!.... Solve the pressure Poisson equation using line GS with SOR

      use const
      use global
      implicit none

      REAL :: rp(nx,ny), aw(nx,ny), ae(nx,ny), as(nx,ny), an(nx,ny), ap(nx,ny)
      REAL :: dp(nx,ny), pl

      do k = 1, kmax

!.... solve for pressure using line GS-SOR in y

!!$     do i = 2, nxm
!!$       dp(i,2:nym) = -ae(i,2:nym) * p(i+1,2:nym) - &
!!$                      aw(i,2:nym) * p(i-1,2:nym) + rp(i,2:nym)
!!$       call TRDIAG( nymm, as(i,2:nym), ap(i,2:nym), &
!!$                       an(i,2:nym), dp(i,2:nym), dp(i,2:nym) )         
!!$       p(i,2:nym) = (one-omega)*p(i,2:nym) + omega*dp(i,2:nym)
!!$     end do
!!$
!!$     p(2:nxm,2:nym) = p(2:nxm,2:nym) - p(2,2)  ! hold a point fixed

!.... solve for pressure using line GS-SOR in x

!$omp parallel do shared(nxmm)
        do j = 2, nym
          dp(2:nxm,j) = -an(2:nxm,j) * p(2:nxm,j+1) - &
                         as(2:nxm,j) * p(2:nxm,j-1) + rp(2:nxm,j)
          call TRDIAG( nxmm, aw(2:nxm,j), ap(2:nxm,j), &
                       ae(2:nxm,j), dp(2:nxm,j), dp(2:nxm,j) )    
          p(2:nxm,j) = (one-omega)*p(2:nxm,j) + omega*dp(2:nxm,j)
        end do

        pl = p(2,2)
!$omp parallel do private(i)
        do j = 2, nym
          do i = 2, nxm 
            p(i,j) = p(i,j) - pl  ! hold a point fixed
          end do
        end do

!.... evaluate the residual norm

        pnorm = zero
!$omp parallel do private(i) reduction(+: pnorm)
        do j = 2, nym
          do i = 2, nxm
            pnorm = pnorm + ( ap(i,j) * p(i,j) + ae(i,j) * p(i+1,j) + &
                    an(i,j) * p(i,j+1) + aw(i,j) * p(i-1,j) + &
                    as(i,j) * p(i,j-1) - rp(i,j) )**2
          end do
        end do
        pnorm = sqrt(pnorm / (nxmm*nymm))

        if (pnorm .lt. tol) exit

      end do

      return
end subroutine line_gs_sor

subroutine adi( rp, aw, ae, as, an, ap )

!.... Solve the pressure Poisson equation using ADI

      use const
      use global
      implicit none
      
      REAL :: rp(nx,ny), aw(nx,ny), ae(nx,ny), as(nx,ny), an(nx,ny), ap(nx,ny)
      REAL :: dp(nx,ny), pl

      do k = 1, kmax

!.... solve for pressure using ADI

!$omp parallel do private(i)
        do j = 2, nym
          do i = 2, nxm
            dp(i,j) = dtau * ( ap(i,j) * p(i,j)   + ae(i,j) * p(i+1,j) + &
                      an(i,j) * p(i,j+1) + aw(i,j) * p(i-1,j) + &
                      as(i,j) * p(i,j-1) - rp(i,j) )
          end do
        end do

!$omp parallel do shared( nymm )
        do i = 2, nxm
          call TRDIAG( nymm, &
                       -dtau * as(i,2:nym), &
                        one + dtau * ( as(i,2:nym) + an(i,2:nym) ), &
                       -dtau * an(i,2:nym), dp(i,2:nym), dp(i,2:nym) ) 
        end do

!$omp parallel do shared( nxmm )
        do j = 2, nym
          call TRDIAG( nxmm, &
                       -dtau * aw(2:nxm,j), &
                       one + dtau * ( aw(2:nxm,j) + ae(2:nxm,j) ), &
                       -dtau * ae(2:nxm,j), dp(2:nxm,j), dp(2:nxm,j) )
        end do

        pl = dp(2,2)
!$omp parallel do private(i)
        do j = 2, nym
          do i = 2, nxm 
            p(i,j) = p(i,j) + dp(i,j) - pl
          end do
        end do

!.... evaluate the residual norm

        pnorm = zero
!$omp parallel do private(i) reduction(+: pnorm)
        do j = 2, nym
          do i = 2, nxm
            pnorm = pnorm + ( ap(i,j) * p(i,j) + ae(i,j) * p(i+1,j) + &
                    an(i,j) * p(i,j+1) + aw(i,j) * p(i-1,j) + &
                    as(i,j) * p(i,j-1) - rp(i,j) )**2
          end do
        end do
        pnorm = sqrt(pnorm / (nxmm*nymm))

        if (pnorm .lt. tol) exit
        
      end do

      return
end subroutine adi

subroutine setdt

!.... set the time-step based on a CFL condition

      use const
      use global
      implicit none

      REAL :: dtinv, dyinv, dxinv

!.... compute the convective CFL (diffusion effects are not included)

      dtinv = zero
!$omp parallel do private(i,dyinv,dxinv) reduction(max: dtinv)
      do j = 2, nym
        dyinv = one/(y(j) - y(j-1))
        do i = 2, nxm
          dxinv = one/(x(i) - x(i-1))
          dtinv = max( abs(uc(i,j)*dxinv)+ abs(vc(i,j)*dyinv), dtinv)
        end do
      end do

      dt = min( dt_max, cfl_max / dtinv )
      cfl = dt * dtinv

      return
end subroutine setdt

subroutine setup

!.... Input run parameters and form the grid

      use timemod
      use const
      use global
      use direct_mod
      implicit none

      integer :: idum
      REAL :: rdum, mu

      namelist /in/ Re, nstep, cfl_max, dt_max, tol, ictype, method, &
                    axis, kmax, omega, dtau, pmethod, base, nout, delta_p

!.... start the clock

      cpu  = second()
      cpu0 = cpu
      
!.... set default values for input variables

      Re=100.0; nstep=1000; cfl_max=0.8; dt_max=1.0e4; tol=1.0e-6; kmax=1000;
      ictype=0; method=0; axis=.false.; rho=one; omega = 1.0; dtau = 0.1;
      pmethod=0; base='output'; nout=100; delta_p = .true.
      
!.... get namelist input

      read(*,in)

      base = base(1:index(base,' ')-1)//char(0) ! null terminate
      filen = char(0)

      vis = one / Re

!.... read grid file

      open(io,file='grid.in',form='formatted',status='old',err=100)
      read(io,*,err=110,end=110) nx, ny
      write(*,"('Read grid with (',i3,',',i3,') CV''s ')") nx-2, ny-2
      allocate( x(nx), y(ny) )
      read(io,*,err=110,end=110) ( x(i), i=1,nx )
      read(io,*,err=110,end=110) ( y(j), j=1,ny )
      close(io)

      nxm  = nx - 1
      nym  = ny - 1
      nxmm = nxm - 1
      nymm = nym - 1

!.... allocate Matrix if using direct solver

      if (pmethod.eq.3) call init_direct

!.... dynamically set SOR parameter

      if (omega.eq.0) then
        mu = one - pt25 * pi**2 * ( one/nxmm**2 + one/nymm**2 )
        omega = two / ( one + sqrt( one - mu**2 ) )
        write(*,"('SOR parameter automatically set to ', 1pe13.6)") omega
      end if

!.... compute coordinates of CV-centers

      allocate ( xc(nx), yc(ny) )
      xc(1) = x(1)
      do i = 2, nxm
        xc(i) = pt5 * ( x(i) + x(i-1) )
      end do
      xc(nx) = x(nxm)

      yc(1) = y(1)
      do j = 2, nym
        yc(j) = pt5 * ( y(j) + y(j-1) )
      end do
      yc(ny) = y(nym)

!.... set the radius

      allocate ( r(ny), fx(nx), fy(ny) )
      if (axis) then
        r = y
      else
        r = one
      end if

!.... Interpolation factors (on scalar CVs)

      fx(1) = zero
      do i = 2, nxm
        fx(i) = (x(i)-x(i-1))/(x(i+1)-x(i-1))
      end do
      fx(nx) = zero

      fy(1) = zero
      do j = 2, nym
        fy(j) = (y(j)-y(j-1))/(y(j+1)-y(j-1))
      end do
      fy(ny) = zero

!.... allocate field variables

      allocate( u(nx,ny), v(nx,ny), p(nx,ny) )
      allocate( u_old(nx,ny), v_old(nx,ny), p_old(nx,ny) )
      allocate(uc(nx,ny),vc(nx,ny),pc(nx,ny),div(nx,ny),wz(nx,ny),ke(nx,ny))

!.... allocate mass flux variables

      allocate( cx(nx,ny), cy(nx,ny), ru(nx,ny), rv(nx,ny) )

!.... initialize all fields to zero (first touch policy for OpenMP)

!$omp parallel do private(i)
      do j = 1, ny
        do i = 1, ny
          u(i,j) = zero; v(i,j) = zero; p(i,j) = zero
          u_old(i,j) = zero; v_old(i,j) = zero; p_old(i,j) = zero
          uc(i,j) = zero; vc(i,j) = zero; pc(i,j) = zero; 
          div(i,j) = zero; wz(i,j) = zero; ke(i,j) = zero
          cx(i,j) = zero; cy(i,j) = zero; ru(i,j) = zero; rv(i,j) = zero
        end do
      end do

      return
100   call error('setup$','Unable to open grid.in$')
110   call error('setup$','Error reading grid.in$')
end subroutine setup

subroutine bc

!.... Set the boundary conditions

      use const
      use global
      implicit none

!.... Bottom

      u(1:nxm,1) = zero   ! no slip
      v(1:nx,1)  = zero

!.... Top

      u(1:nxm,ny) = zero  ! no slip
      v(1:nx,nym) = zero

!.... Left

      u(1,1:ny)  = zero   ! no slip
      v(1,1:nym) = zero

!.... Right

      u(nxm,1:ny) = zero  ! no slip
      v(nx,1:nym) = zero
      
      return
end subroutine bc

subroutine ic

!.... Set the initial condition
     
      use const
      use global
      implicit none
      REAL :: x0 = 2.0, y0 = 16.0

      n = 0         ! time-step zero
      time = zero

!.... Initialize to no-flow

!$omp parallel do private(i)
      do j = 1, ny
        do i = 1, nx
          u(i,j) = zero
          v(i,j) = zero
          p(i,j) = zero
        end do
      end do

!.... Add two Oseen vortices

      call Oseen( -one, -x0, y0, u, v )
      call Oseen(  one,  x0, y0, u, v )

!.... enforce the boundary conditions

      call bc

!.... compute mass fluxes and project to divergence free space

      call calc_p( one, u, v, .true. )

!.... compute the pressure for the initial condition from the divergence
!.... of the discrete momentum equations

      call calc_u( .false. )
      call calc_v( .false. )
      call calc_p( one,  ru, rv, .false. )

!.... interpolate to cell centered mesh and calculate statistics

      call interp
      call stats

!.... save the IC in a restart file and Plot3d file

      iver = 0
      call addver(base, 'R'//char(0), iver, filen, lfile)
      open(io,file=filen,form='unformatted')
      write(io) n, nx, ny
      write(io) Re, time
      write(io) u, v, p, cx, cy
      close(io)
      iver = 0
      call addver(base, 'q'//char(0), iver, filen, lfile)
      call wdata( filen )

      return
end subroutine ic

subroutine wdata( fname )

!.... Write a Plot3d solution file assuming that 
!.... pc, uc, vc, wz, and ke are valid

      use const
      use global
      implicit none

      character*80 :: fname

!.... output Plot3d solution-file at CV-centers

      open(io,file=fname,form='unformatted')
      write(io) nx, ny, 1
      write(io) zero, zero, Re, time
      write(io) (( pc(i,j), i=1,nx), j=1,ny), &
                (( uc(i,j), i=1,nx), j=1,ny), &
                (( vc(i,j), i=1,nx), j=1,ny), &
                (( wz(i,j), i=1,nx), j=1,ny), &
                (( ke(i,j), i=1,nx), j=1,ny)
      close(io)

      return
end subroutine wdata

subroutine wgrid

!.... write a Plot3d grid-file at CV-centers

      use const
      use global
      implicit none

      open(io,file='grid.xyz',form='unformatted')
      write(io) nx, ny, 1
      write(io) (( xc(i), i = 1, nx), j = 1, ny), &
                (( yc(j), i = 1, nx), j = 1, ny), &
                ((  zero, i = 1, nx), j = 1, ny)
      close(io)
      
      return
end subroutine wgrid

subroutine stats

!.... Computes divergence, vorticity, kinetic energy, and total enstrophy

      use const
      use global
      implicit none
      
      REAL :: area

!.... Vorticity must first be computed on the nodal mesh and then interpolated
!.... to the cell centered mesh

!$omp parallel do private(i)
      do j = 1, nym
        do i = 1, nxm
          div(i,j) = (v(i+1,j)-v(i,j))/(xc(i+1)-xc(i)) - &
                     (u(i,j+1)-u(i,j))/(yc(j+1)-yc(j))
        end do
      end do

!.... Now interpolate back to the cell centered mesh

      wz(1,1)      = div(1,1)
      wz(1,2:nym)  = pt5 * (div(1,2:nym) + div(1,1:nymm))
      wz(1,ny)     = div(1,nym)
      wz(2:nxm,ny) = pt5 * (div(2:nxm,nym) + div(1:nxmm,nym))
      wz(nx,ny)    = div(nxm,nym)
      wz(nx,2:nym) = pt5 * (div(nxm,2:nym) + div(nxm,1:nymm))
      wz(nx,1)     = div(nxm,1)
      wz(2:nxm,1)  = pt5 * (div(2:nxm,1) + div(1:nxmm,1))

!$omp parallel do private(i)
      do j = 2, nym
        do i = 2, nxm
          wz(i,j) = pt25 * ( div(i,j) + div(i-1,j) + &
                             div(i,j-1) + div(i-1,j-1) )
        end do
      end do

!.... Divergence and Kinetic energy are easily computed 
!.... on the cell centered mesh
      
!$omp parallel do private(i)
      do j = 2, nym
        do i = 2, nxm
          div(i,j) = rho * ( (u(i,j) - u(i-1,j)) * (y(j)-y(j-1)) + &
                             (v(i,j) - v(i,j-1)) * (x(i)-x(i-1)) )
          ke(i,j)  = pt5 * ( uc(i,j)**2 + vc(i,j)**2 )
        end do
      end do

!.... Integrate to get the total Energy and Enstrophy

      tke = zero; enstrophy = zero; tdiv = zero
!$omp parallel do private(i, area) reduction(+: tdiv, tke, enstrophy)
      do j = 2, nym
        do i = 2, nxm
          area = (x(i)-x(i-1)) * (y(j)-y(j-1))
          tdiv = tdiv + div(i,j)**2
          tke = tke + ke(i,j) * area
          enstrophy = enstrophy + wz(i,j)**2 * area
        end do
      end do
      enstrophy = pt5 * enstrophy
      tdiv = sqrt( tdiv / (nxmm*nymm) )

      write(istat,"(8(1pe13.6,1x))") time, tdiv, tke, enstrophy

      return
end subroutine stats

subroutine Oseen( gamma, x0, y0, ul, vl )

!.... Computes the velocity field for an Oseen vortex centered at
!.... (x0,y0) with circulation Gamma

      use global
      use const
      implicit none

      REAL :: ul(nx,ny), vl(nx,ny)
      REAL :: gamma, x0, y0, rr, theta, xi, alpha = 1.2564312086261696770
      REAL :: v_theta, v_r

!.... Oseen vortex at (x0, y0)

      xi = one - exp( -alpha )

!$omp parallel do private(i, rr, theta, v_theta, v_r)
      do j = 1, ny
        do i = 1, nxm
          rr = sqrt( (x(i)-x0)**2 + (yc(j)-y0)**2 )
          theta = atan2( yc(j)-y0, x(i)-x0 )

          if (rr.eq.zero) then
            v_theta = zero
            v_r = zero
          else
            v_theta = gamma * ( one - exp( -alpha * rr**2 ) ) / ( rr * xi )
            v_r = zero
          endif

          ul(i,j) = ul(i,j) + cos(theta) * v_r - sin(theta) * v_theta
        end do
      end do

!$omp parallel do private(i, rr, theta, v_theta, v_r)
      do j = 1, nym
        do i = 1, nx
          rr = sqrt( (xc(i)-x0)**2 + (y(j)-y0)**2 )
          theta = atan2( y(j)-y0, xc(i)-x0 )

          if (rr.eq.zero) then
            v_theta = zero
            v_r = zero
          else
            v_theta = gamma * ( one - exp( -alpha * rr**2 ) ) / ( rr * xi )
            v_r = zero
          endif

          vl(i,j) = vl(i,j) + sin(theta) * v_r + cos(theta) * v_theta
        end do
      end do

      return
end subroutine Oseen

subroutine interp

!.... Interpolate velocities and pressure to cell centered mesh
     
      use const
      use global
      implicit none

      pc(1,1)      = p(2,2)
      pc(nx,1)     = p(nxm,2)
      pc(1,ny)     = p(2,nym)
      pc(nx,ny)    = p(nxm,nym)
      pc(2:nxm,1)  = p(2:nxm,2)
      pc(2:nxm,ny) = p(2:nxm,nym)
      pc(1,2:nym)  = p(2,2:nym)
      pc(nx,2:nym) = p(nxm,2:nym)
!$omp parallel do private(i)
      do j = 2, nym
        do i = 2, nxm
          pc(i,j) = p(i,j)
        end do
      end do

      uc(1,1) = u(1,1)
      vc(1,1) = v(1,1)
      uc(nx,1) = u(nxm,1)
      vc(nx,1) = v(nx,1)

!$omp parallel do
      do j = 2, nym
        uc(1,j) = u(1,j)
        vc(1,j) = pt5*(v(1,j)+v(1,j-1))
        uc(nx,j) = u(nxm,j)
        vc(nx,j) = pt5*(v(nx,j)+v(nx,j-1))
      end do

      uc(1,ny) = u(1,ny)
      vc(1,ny) = v(1,nym)
      uc(nx,ny) = u(nxm,ny)
      vc(nx,ny) = v(nx,nym)

!$omp parallel do
      do i = 2, nxm
        uc(i,1) = pt5*(u(i,1)+u(i-1,1))
        vc(i,1) = v(i,1)
        uc(i,ny) = pt5*(u(i,ny)+u(i-1,ny))
        vc(i,ny) = v(i,nym)
      end do

!$omp parallel do private(i)
      do j = 2, nym
        do i = 2, nxm
          uc(i,j) = pt5*(u(i,j)+u(i-1,j))
          vc(i,j) = pt5*(v(i,j)+v(i,j-1))
        end do
      end do

      return
end subroutine interp

subroutine restart

!.... restart the computation using data from a previous run

!.... Assumes that the flow field is divergence free

      use const
      use global
      implicit none

      integer :: inx, iny, idum
      REAL :: rdum

!.... read a restart file

      iver = igetver( base, 'R'//char(0) )
      call addver(base, 'R'//char(0), iver, filen, lfile)
      write (*,10) filen(1:index(filen,' ')-1)
10    format(/,'Restarting from file:  ',a)

      open(io,file=filen,form='unformatted')
      read(io) n, inx, iny
      if (inx.ne.nx .or. iny.ne.ny) &
        call error('restart$','Incompatible restart file.$')
      read(io) rdum, time
      read(io) u, v, p, cx, cy
      close(io)

      write(*,"('Restarting at step ',i5,' at time = ',1pe13.6)") n, time

      return
end subroutine restart

subroutine error(name,msg)

!.... Error handler

      implicit none

      integer loc
      character*80 name, msg

      loc = index(name,'$')-1
      write(*,"(/,'*****************************************************')")
      write(*,"('Error in --> ',a)") name(1:loc)
      loc = index(msg,'$')-1
      write(*,"('-----------> ',a)") msg(1:loc)
      write(*,"('*****************************************************',/)")

      stop
end subroutine error

subroutine TRDIAG (N,A,B,C,X,G)    
      dimension A(N),B(N),C(N),X(N),G(N),BB(N)
!.....THIS SUBROUTINE SOLVES TRIDIAGONAL SYSTEMS OF EQUATIONS          
!.....BY GAUSS ELIMINATION     
!.....THE PROBLEM SOLVED IS MX=G where M=TRI(A,B,C)          
!.....THIS ROUTINE DOES NOT DESTROY THE ORIGINAL MATRIX      
!.....AND MAY BE CALLED A NUMBER OF TIMES WITHOUT REDEFINING 
!.....THE MATRIX     
!.....N = NUMBER OF EQUATIONS SOLVED
!.....FORWARD ELIMINATION
!.....BB IS A SCRATCH ARRAY NEEDED TO AVOID DESTROYING B ARRAY         
      do I=1,N     
        BB(I) = B(I)   
      end do
      do I=2,N     
        T = A(I)/BB(I-1)         
        BB(I) = BB(I) - C(I-1)*T 
        G(I) = G(I) - G(I-1)*T   
      end do
!.....BACK SUBSTITUTION        
      X(N) = G(N)/BB(N)        
      do I=1,N-1   
        J = N-I        
        X(J) = (G(J)-C(J)*X(J+1))/BB(J)    
      end do
      return         
end subroutine TRDIAG
