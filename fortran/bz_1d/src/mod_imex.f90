module mod_imex

use mod_precision
use mod_coo_matrix

contains

subroutine second_order_multistep_imex_integration(f_reac, build_diff_mat, n, u, nt, tini, tend)
  implicit none
  include 'dmumps_struc.h'
  ! arguments
  interface
    function f_reac(n, y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f_reac
    end function f_reac
    subroutine build_diff_mat(a) 
      use mod_coo_matrix 
      type(coo_matrix) :: a 
    end subroutine build_diff_mat
  end interface
  integer :: n
  real(kind=dp), dimension(:) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  ! local variables
  type(coo_matrix) :: diff_mat
  type(dmumps_struc) mumps
  real(kind=dp), dimension(n) :: un, unm1
  real(kind=dp) :: dt
  integer :: it, innz

  dt = (tend-tini)/(nt-1)

  unm1 = u

  call build_diff_mat(diff_mat) 

  mumps%comm = 0
  mumps%par  = 1
  mumps%sym  = 0

  ! initialize an instance of mumps
  mumps%job  = -1
  call dmumps(mumps)

  write(6, '(A13, ES14.7)') "dt used : ", dt
  print *, "  MUMPS version used : ", mumps%VERSION_NUMBER 

  ! first iteration : first order scheme

  ! allocate and initialize matrix of first order imex scheme
  mumps%n  = diff_mat%n
  mumps%nnz = diff_mat%nnz
           
  allocate(mumps%irn(mumps%nnz))
  allocate(mumps%jcn(mumps%nnz))
  allocate(mumps%a(mumps%nnz))
  allocate(mumps%rhs(mumps%n))

  mumps%irn = diff_mat%row
  mumps%jcn = diff_mat%col
  do innz = 1, mumps%nnz
    if (mumps%irn(innz) == mumps%jcn(innz)) then
      mumps%a(innz) = 1.d0 - dt*diff_mat%val(innz)
    else
      mumps%a(innz) =  -dt*diff_mat%val(innz)
    end if
    mumps%rhs = u + dt*f_reac(n,u)
  end do

  mumps%icntl(4) = 1

  ! compute solution
  mumps%job = 6
  call dmumps(mumps)

  un = mumps%rhs

  ! other iterations : second order scheme

  ! initialize second order imex matrix
  do innz = 1, mumps%nnz
    if (mumps%irn(innz) == mumps%jcn(innz)) then
      mumps%a(innz) = 1.5d0 - dt*diff_mat%val(innz)
    else
      mumps%a(innz) =  -dt*diff_mat%val(innz)
    end if
  end do

  ! factorize matrix
  mumps%job = 4
  call dmumps(mumps)

  do it = 2, nt-1

    !!print *, "  it : ", it

    mumps%rhs = 2.d0*un - 0.5d0*unm1 + 2.d0*dt*f_reac(n,un) - dt*f_reac(n,unm1)

    ! solve 
    mumps%job = 3
    call dmumps(mumps)

    unm1 = un
    un = mumps%rhs
    !!print *, it, un(1)

  end do

  u = mumps%rhs

  deallocate(mumps%irn)
  deallocate(mumps%jcn)
  deallocate(mumps%a)
  deallocate(mumps%rhs)

  ! destroy the instance of mumps
  mumps%job = -2
  call dmumps(mumps)

end subroutine second_order_multistep_imex_integration


subroutine second_order_2_stages_rk_imex_integration(f_reac, f_diff, build_diff_mat, n, u, nt, tini, tend)
  implicit none
  include 'dmumps_struc.h'
  ! arguments
  interface
    function f_reac(n, y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f_reac
    end function f_reac
    function f_diff(n, y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f_diff
    end function f_diff
    subroutine build_diff_mat(a)
      use mod_coo_matrix
      type(coo_matrix) :: a
    end subroutine build_diff_mat
  end interface
  integer :: n
  real(kind=dp), dimension(:) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  ! local variables
  type(coo_matrix) :: a_diff
  type(dmumps_struc) mumps
  real(kind=dp), dimension(n) :: un, u1, u2
  real(kind=dp) :: lambda
  real(kind=dp) :: dt
  integer :: it, innz

  lambda = 1.d0 - (sqrt(2.d0)/2.d0)
  dt = (tend-tini)/(nt-1)

  call build_diff_mat(a_diff)

  mumps%comm = 0
  mumps%par  = 1
  mumps%sym  = 0

  ! initialize an instance of mumps
  mumps%job  = -1
  call dmumps(mumps)

  write(6, '(A13, ES14.7)') "dt used : ", dt
  print *, "  MUMPS version used : ", mumps%VERSION_NUMBER

  ! allocate and initialize matrix of imex rk scheme
  mumps%n  = a_diff%n
  mumps%nnz = a_diff%nnz

  allocate(mumps%irn(mumps%nnz))
  allocate(mumps%jcn(mumps%nnz))
  allocate(mumps%a(mumps%nnz))
  allocate(mumps%rhs(mumps%n))

  mumps%irn = a_diff%row
  mumps%jcn = a_diff%col
  do innz = 1, mumps%nnz
    if (mumps%irn(innz) == mumps%jcn(innz)) then
      mumps%a(innz) = 1.d0 - lambda*dt*a_diff%val(innz)
    else
      mumps%a(innz) =  -lambda*dt*a_diff%val(innz)
    end if
  end do

  mumps%icntl(4) = 1

  ! factorize matrix
  mumps%job = 4
  call dmumps(mumps)

  un = u

  do it = 1, nt-1

    !!print *, "  it : ", it

    ! solve first linear system
    mumps%rhs = un

    mumps%job = 3
    call dmumps(mumps)

    u1 = mumps%rhs

    ! solve second linear system
    mumps%rhs = un + dt*f_reac(n,un) + dt*(1.d0-2.d0*lambda)*f_diff(n,u1)

    mumps%job = 3
    call dmumps(mumps)

    u2 = mumps%rhs

    un = un + (0.5*dt)*(f_reac(n,u1)+f_reac(n,u2)) + dt*f_diff(n,0.5*(u1+u2))

  end do

  u = un

  deallocate(mumps%irn)
  deallocate(mumps%jcn)
  deallocate(mumps%a)
  deallocate(mumps%rhs)

  ! destroy the instance of mumps
  mumps%job = -2
  call dmumps(mumps)

end subroutine second_order_2_stages_rk_imex_integration

subroutine second_order_3_stages_rk_imex_integration(f_reac, f_diff, build_diff_mat, n, u, nt, tini, tend)
  implicit none
  include 'dmumps_struc.h'
  ! arguments
  interface
    function f_reac(n, y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f_reac
    end function f_reac
    function f_diff(n, y)
       integer :: n
       real(kind=kind(1.0d0)), dimension(n) :: y, f_diff
    end function f_diff
    subroutine build_diff_mat(a)
      use mod_coo_matrix
      type(coo_matrix) :: a
    end subroutine build_diff_mat
  end interface
  integer :: n
  real(kind=dp), dimension(:) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  ! local variables
  type(coo_matrix) :: a_diff
  type(dmumps_struc) mumps
  real(kind=dp), dimension(n) :: un, u2, u3
  real(kind=dp) :: lambda
  real(kind=dp) :: dt
  integer :: it, innz

  lambda = 1.d0 - (sqrt(2.d0)/2.d0)
  dt = (tend-tini)/(nt-1)

  call build_diff_mat(a_diff)

  mumps%comm = 0
  mumps%par  = 1
  mumps%sym  = 0

  ! initialize an instance of mumps
  mumps%job  = -1
  call dmumps(mumps)

  write(6, '(A13, ES14.7)') "dt used : ", dt
  print *, "  MUMPS version used : ", mumps%VERSION_NUMBER

  ! allocate and initialize matrix of imex rk scheme
  mumps%n  = a_diff%n
  mumps%nnz = a_diff%nnz

  allocate(mumps%irn(mumps%nnz))
  allocate(mumps%jcn(mumps%nnz))
  allocate(mumps%a(mumps%nnz))
  allocate(mumps%rhs(mumps%n))

  mumps%irn = a_diff%row
  mumps%jcn = a_diff%col
  do innz = 1, mumps%nnz
    if (mumps%irn(innz) == mumps%jcn(innz)) then
      mumps%a(innz) = 1.d0 - lambda*dt*a_diff%val(innz)
    else
      mumps%a(innz) =  -lambda*dt*a_diff%val(innz)
    end if
  end do

  mumps%icntl(4) = 1

  ! factorize matrix
  mumps%job = 4
  call dmumps(mumps)

  un = u

  do it = 1, nt-1

    !!print *, "  it : ", it

    ! solve first linear system
    mumps%rhs = un + dt*lambda*f_reac(n,un)

    mumps%job = 3
    call dmumps(mumps)

    u2 = mumps%rhs

    ! solve second linear system
    mumps%rhs = un + dt*(lambda-1)*f_reac(n,un) + dt*2.d0*(1.d0-lambda)*f_reac(n,u2) + dt*(1.d0-2.d0*lambda)*f_diff(n,u2)

    mumps%job = 3
    call dmumps(mumps)

    u3 = mumps%rhs

    un = un + (0.5*dt)*(f_reac(n,u2)+f_reac(n,u3)) + dt*f_diff(n,0.5*(u2+u3))
    !!print *, it, un(1)

  end do

  u = un

  deallocate(mumps%irn)
  deallocate(mumps%jcn)
  deallocate(mumps%a)
  deallocate(mumps%rhs)

  ! destroy the instance of mumps
  mumps%job = -2
  call dmumps(mumps)

end subroutine second_order_3_stages_rk_imex_integration

end module mod_imex
