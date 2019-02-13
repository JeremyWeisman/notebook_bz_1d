module mod_bz_2eq_1d

use mod_precision
use mod_cartesian_grid
use mod_coo_matrix

integer, private, parameter :: neq = 2
real(kind=dp), private, parameter :: eps = 1.d-02
real(kind=dp), private, parameter :: f = 3.d0
real(kind=dp), private, parameter :: q = 2.d-03
real(kind=dp), private, parameter :: db = 1.d0/400.d0
real(kind=dp), private, parameter :: dc = 0.6d0/400.d0
real(kind=dp), private :: dboverdxdx
real(kind=dp), private :: dcoverdxdx
real(kind=dp), private :: oneovereps
integer, private :: nx
real(kind=dp), private :: xmin

contains

subroutine init_bz_2eq(grid)
  implicit none
  type(cartesian_grid_type) :: grid
  real(kind=dp) :: dx

  dx = grid%dx

  dboverdxdx = db / (dx*dx)
  dcoverdxdx = dc / (dx*dx)

  oneovereps = 1.d0 / eps

  nx = grid%nx

  xmin = grid%xmin

end subroutine init_bz_2eq

function f_bz_2eq(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_bz_2eq
  ! local 
  real(kind=dp) :: bim1, bi, bip1
  real(kind=dp) :: cim1, ci, cip1
  integer :: inx, irow

  ! left boundary
  bi   = u(1)
  ci   = u(2)
  bip1 = u(3)
  cip1 = u(4)

  f_bz_2eq(1) = (dboverdxdx * (bip1 - bi)) + oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
  f_bz_2eq(2) = (dcoverdxdx * (cip1 - ci)) + (bi - ci)

  do inx = 2, nx-1

    irow = 1 + (inx-1)*neq

    bim1 = u(irow-2)
    cim1 = u(irow-1)
    bi   = u(irow)
    ci   = u(irow+1)
    bip1 = u(irow+2)
    cip1 = u(irow+3)

    f_bz_2eq(irow)   = (dboverdxdx * (bim1 - 2.d0*bi + bip1)) + oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
    f_bz_2eq(irow+1) = (dcoverdxdx * (cim1 - 2.d0*ci + cip1)) + (bi - ci)

  end do 

  ! right boundary
  bim1 = u(n-3)
  cim1 = u(n-2)
  bi   = u(n-1)
  ci   = u(n)

  f_bz_2eq(n-1) = (dboverdxdx * (bim1 - bi)) +  oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
  f_bz_2eq(n)   = (dcoverdxdx * (cim1 - ci)) + (bi - ci)

end function f_bz_2eq

function f_loc_bz_2eq_reac(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_loc_bz_2eq_reac
  ! local
  real(kind=dp) :: bi
  real(kind=dp) :: ci

  bi   = u(1)
  ci   = u(2)

  f_loc_bz_2eq_reac(1) = oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
  f_loc_bz_2eq_reac(2) = (bi - ci)

end function f_loc_bz_2eq_reac

function f_bz_2eq_reac(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_bz_2eq_reac
  ! local 
  real(kind=dp) :: bi
  real(kind=dp) :: ci
  integer :: inx, irow

  do inx = 1, nx

   irow = 1 + (inx-1)*neq

    bi   = u(irow)
    ci   = u(irow+1)

    f_bz_2eq_reac(irow)   = oneovereps * (bi - bi*bi + ((f*(q-bi)*ci)/(q+bi)))
    f_bz_2eq_reac(irow+1) = (bi - ci)

  end do 

end function f_bz_2eq_reac

function f_bz_2eq_diff(n, u)
  implicit none
  integer :: n
  real(kind=dp), dimension(n) :: u
  real(kind=dp), dimension(n) :: f_bz_2eq_diff
  ! local 
  real(kind=dp) :: bim1, bi, bip1
  real(kind=dp) :: cim1, ci, cip1
  integer :: inx, irow

  ! left boundary
  bi   = u(1)
  ci   = u(2)
  bip1 = u(3)
  cip1 = u(4)

  f_bz_2eq_diff(1) = (dboverdxdx * (bip1 - bi))
  f_bz_2eq_diff(2) = (dcoverdxdx * (cip1 - ci))

  do inx = 2, nx-1

    irow = 1 + (inx-1)*neq

    bim1 = u(irow-2)
    cim1 = u(irow-1)
    bi   = u(irow)
    ci   = u(irow+1)
    bip1 = u(irow+2)
    cip1 = u(irow+3)

    f_bz_2eq_diff(irow)   = (dboverdxdx * (bim1 - 2.d0*bi + bip1))
    f_bz_2eq_diff(irow+1) = (dcoverdxdx * (cim1 - 2.d0*ci + cip1))

  end do 

  ! right boundary
  bim1 = u(n-3)
  cim1 = u(n-2)
  bi   = u(n-1)
  ci   = u(n)

  f_bz_2eq_diff(n-1) = (dboverdxdx * (bim1 - bi))
  f_bz_2eq_diff(n)   = (dcoverdxdx * (cim1 - ci))

end function f_bz_2eq_diff

subroutine build_bz_2eq_diff_matrix(a)
  implicit none
  type(coo_matrix) :: a
  ! local variables
  integer :: ntot, nnz
  integer :: innz, inx, irow

  ntot = neq*nx
  nnz = 2*2 + 2*3*(nx-2) + 2*2
  call create_coo_matrix(ntot, nnz, a)

  innz = 1

  ! left boundary
  a%row(innz) = 1
  a%col(innz) = 1
  a%val(innz) = -dboverdxdx
  innz = innz + 1
  a%row(innz) = 1
  a%col(innz) = 1+neq
  a%val(innz) = dboverdxdx
  innz = innz + 1

  a%row(innz) = 2
  a%col(innz) = 2
  a%val(innz) = -dcoverdxdx
  innz = innz + 1
  a%row(innz) = 2
  a%col(innz) = 2+neq
  a%val(innz) = dcoverdxdx
  innz = innz + 1

  do inx = 2, nx-1

    irow = 2*(inx-1) + 1

    a%row(innz) = irow
    a%col(innz) = irow - neq
    a%val(innz) = dboverdxdx
    innz = innz + 1
    a%row(innz) = irow
    a%col(innz) = irow
    a%val(innz) = -2.d0 * dboverdxdx
    innz = innz + 1
    a%row(innz) = irow
    a%col(innz) = irow + neq
    a%val(innz) = dboverdxdx
    innz = innz + 1

    a%row(innz) = irow + 1
    a%col(innz) = irow + 1 - neq
    a%val(innz) = dcoverdxdx
    innz = innz + 1
    a%row(innz) = irow + 1
    a%col(innz) = irow + 1
    a%val(innz) = -2.d0 * dcoverdxdx
    innz = innz + 1
    a%row(innz) = irow + 1
    a%col(innz) = irow + 1 + neq
    a%val(innz) = dcoverdxdx
    innz = innz + 1

  end do

  ! right boundary
  a%row(innz) = ntot - 1
  a%col(innz) = ntot - 1 - neq
  a%val(innz) = dboverdxdx
  innz = innz + 1
  a%row(innz) = ntot - 1
  a%col(innz) = ntot - 1
  a%val(innz) = -dboverdxdx
  innz = innz + 1

  a%row(innz) = ntot
  a%col(innz) = ntot - neq
  a%val(innz) = dcoverdxdx
  innz = innz + 1
  a%row(innz) = ntot
  a%col(innz) = ntot
  a%val(innz) = -dcoverdxdx

  !!do innz = 1, nnz
  !!  print *, a%val(innz)
  !!end do

end subroutine build_bz_2eq_diff_matrix

subroutine bz_2eq_init_sol(u)
  implicit none
  real(kind=dp), dimension(:) :: u
  ! local
  real(kind=dp), allocatable, dimension(:) :: b, c
  real(kind=dp) :: xcoor, ycoor, ylim
  real(kind=dp), parameter :: pi = 4.d0 * atan(1.d0)
  integer :: inx, irow

  allocate(b(nx))
  allocate(c(nx))

  ylim = 0.05d0

  do  inx = 1, int(nx/20)
    xcoor = 0.5d0
    ycoor = dfloat(inx)/dfloat(nx/20) - ylim

    if (ycoor >= 0.d0 .and. ycoor<= 0.3*xcoor) then
      b(inx) = 0.8d0
    else
      b(inx) = q*(f+1.d0)/(f-1.d0)
    endif

   if ( ycoor>=0.d0 ) then
     c(inx) = q*(f+1.d0)/(f-1.d0) + datan(ycoor/xcoor)/(8.d0*pi*f)
   else
     c(inx) = q*(f+1.d0)/(f-1.d0) + (datan(ycoor/xcoor) + 2.d0*pi)/(8.d0*pi*f)
   endif

  enddo

  do inx=int(nx/20)+1, nx
    b(inx)=b(int(nx/20))
    c(inx)=c(int(nx/20))
  end do

  do inx = 1, nx
    irow = 1 + (inx-1)*neq
    u(irow)   = b(inx)
    u(irow+1) = c(inx)
  end do

  deallocate(b, c)

end subroutine bz_2eq_init_sol

end module mod_bz_2eq_1d
