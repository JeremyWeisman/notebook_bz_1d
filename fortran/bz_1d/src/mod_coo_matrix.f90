module mod_coo_matrix

  implicit none

!*******************************************************************!
! define coo_matrix type                                            !
!*******************************************************************!
  type coo_matrix
     
  ! dimension of the matrix
  integer n
  ! total number of non zero coefficients
  integer nnz
  ! array of row indices
  integer, allocatable, dimension(:) :: row
  ! array of col indices
  integer, allocatable, dimension(:) :: col
  ! array of value 
  real(kind=kind(0.d0)), allocatable, dimension(:) :: val

  end type coo_matrix

contains

!*******************************************************************!
! create coo matrix                                                 !
!*******************************************************************!
  subroutine create_coo_matrix(n, nnz, a)
    integer :: n, nnz
    type(coo_matrix) :: a
    a%n = n
    a%nnz = nnz
    allocate(a%row(nnz))
    allocate(a%col(nnz))
    allocate(a%val(nnz))
  end subroutine create_coo_matrix

!*******************************************************************!
! destroy coo matrix                                                !
!*******************************************************************!
  subroutine destroy_coo_matrix(a)
    type(coo_matrix) :: a
    deallocate(a%val)
    deallocate(a%col)
    deallocate(a%row)
  end subroutine destroy_coo_matrix

end module mod_coo_matrix
