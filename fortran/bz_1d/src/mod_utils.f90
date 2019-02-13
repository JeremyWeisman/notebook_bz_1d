module mod_utils

use mod_precision
use mod_cartesian_grid

contains

subroutine read_data(neq, tini, tend, nt, xmin, xmax, nxib, method, tol, create_ref_sol)
  implicit none
  integer :: neq
  real(kind=dp) :: tini, tend
  integer :: nt
  real(kind=dp) :: xmin, xmax
  integer :: nxib
  character(len=20) :: method
  real(kind=dp) :: tol
  integer :: create_ref_sol
  logical :: existe

  open(unit=10, file="input.dat")
  read(10, *) neq
  read(10, *) tini
  read(10, *) tend
  read(10, *) nt
  read(10, *) xmin
  read(10, *) xmax
  read(10, *) nxib
  read(10, *) method
  read(10, *) tol
  read(10, *) create_ref_sol
  method = trim(method)
  close(10)

end subroutine read_data


subroutine save_sol(file_name, neq, grid, t, u)
  implicit none
  character(len=*) :: file_name
  integer :: neq
  type(cartesian_grid_type) :: grid
  real(kind=dp) :: t
  real(kind=dp), dimension(:) :: u
  ! local variables
  real(kind=dp) :: xi
  real(kind=dp) :: xmin
  real(kind=dp) :: dx
  integer :: nx
  integer :: inx, irow, i

  nx = grid%nx

  dx = grid%dx

  xmin = grid%xmin + dx

  open(unit=10, file=file_name)

  write(10,*) "#", neq, nx, t

  do inx = 1, nx
    irow = 1 + (inx-1)*neq
    !!print *, irow, irow+1
    xi = xmin  + (inx-1)*dx
    write(10,*) xi, (u(irow+i-1), i = 1, neq)
  end do

  close(10)

end subroutine save_sol

subroutine compute_error(neq, grid, u)
  implicit none
  integer :: neq
  type(cartesian_grid_type) :: grid
  real(kind=dp), dimension(:) :: u
  ! local variables
  real(kind=dp), allocatable, dimension(:) :: a, b, c
  real(kind=dp) :: norm_err_a, norm_err_b, norm_err_c
  integer :: nx
  real(kind=dp) :: xi
  integer :: ntot
  character(len=10) :: str
  integer :: read_nx, read_neq
  integer :: inx, irow
  integer :: ierr

  nx = grid%nx
  ntot = neq * nx

  ! read quasi exact solution
  open(unit=10, file="sol_ref.dat", action="read", iostat=ierr)
  if (ierr /= 0) then
    print *
    print *, "Cannot open sol_ref.dat to compute error"
    return
  end if

  read(10,*) str, read_neq, read_nx
  if (nx /= read_nx) then
    print *, "  Bad dimension of vector read "
    print *, "  Cannot compute norm of error "
    close(10)
    return
  end if

  if (neq==2) then
    allocate(b(nx), c(nx)) 
    do inx = 1, nx
      irow =  1 + (inx-1)*neq
      read(10,*) xi, b(inx), c(inx)
    end do
   
    norm_err_b = 0.d0
    norm_err_c = 0.d0
    do inx = 1, nx
      irow =  1 + (inx-1)*neq
      norm_err_b = norm_err_b + ((u(irow)   - b(inx)) * (u(irow)   - b(inx)))
      norm_err_c = norm_err_c + ((u(irow+1) - c(inx)) * (u(irow+1) - c(inx)))
    end do 
    norm_err_b = sqrt(norm_err_b/nx)
    norm_err_c = sqrt(norm_err_c/nx)
    print *
    print "(a, es20.12)", " Variable b : || unum - uexa || = ", norm_err_b
    print "(a, es20.12)", " Variable c : || unum - uexa || = ", norm_err_c
   
    deallocate(b, c)

  else if (neq==3) then 
    allocate(a(nx), b(nx), c(nx)) 
    do inx = 1, nx
      irow =  1 + (inx-1)*neq
      read(10,*) xi, a(inx), b(inx), c(inx)
    end do
   
    norm_err_a = 0.d0
    norm_err_b = 0.d0
    norm_err_c = 0.d0
    do inx = 1, nx
      irow =  1 + (inx-1)*neq
      norm_err_a = norm_err_a + ((u(irow)   - a(inx)) * (u(irow)   - a(inx)))
      norm_err_b = norm_err_b + ((u(irow+1) - b(inx)) * (u(irow+1) - b(inx)))
      norm_err_c = norm_err_c + ((u(irow+2) - c(inx)) * (u(irow+2) - c(inx)))
    end do
    norm_err_a = sqrt(norm_err_a/nx)
    norm_err_b = sqrt(norm_err_b/nx)
    norm_err_c = sqrt(norm_err_c/nx)
    print *
    print "(a, es20.12)", " Variable a : || unum - uexa || = ", norm_err_a
    print "(a, es20.12)", " Variable b : || unum - uexa || = ", norm_err_b
    print "(a, es20.12)", " Variable b : || unum - uexa || = ", norm_err_c
   
    deallocate(a, b, c)
  end if

end subroutine compute_error

end module mod_utils
