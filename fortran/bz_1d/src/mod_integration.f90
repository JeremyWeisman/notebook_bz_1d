module mod_integration

use mod_precision
use mod_radau
use mod_strang
use mod_imex
use mod_bz_2eq_1d
use mod_bz_3eq_1d

contains

subroutine integrate(method, tol, neq, nx, u, nt, tini, tend)
  implicit none
  ! arguments
  character(len=20) :: method
  real(kind=dp) :: tol
  integer :: neq
  integer :: nx
  real(kind=dp), dimension(:) :: u
  integer :: nt
  real(kind=dp) :: tini, tend
  integer :: ntot

  ntot = neq * nx

  select case (method) 
    case ("radau5")
      print *
      print *, "Radau5 integration"
      if (neq == 2) then
        call radau5_integration(f_bz_2eq, tol, neq, ntot, u, tini, tend, .true.)
      else
        call radau5_integration(f_bz_3eq, tol, neq, ntot, u, tini, tend, .true.)
      end if
    case ("strang")
      print *
      print *, "Strang splitting integration"
      print *, "  RADAU5 for reaction"
      print *, "  ROCK4  for diffusion"
      if (neq == 2) then
        call strang_integration(f_loc_bz_2eq_reac, f_bz_2eq_diff, tol, neq, nx, u, nt, tini, tend) 
      else
        call strang_integration(f_loc_bz_3eq_reac, f_bz_3eq_diff, tol, neq, nx, u, nt, tini, tend) 
      end if 
    case ("imex_ms_2")
      print *
      print *, "Second order multistep IMEX integration"
      if (neq == 2) then
        call second_order_multistep_imex_integration(f_bz_2eq_reac, build_bz_2eq_diff_matrix, ntot, u, nt, tini, tend) 
      else
        call second_order_multistep_imex_integration(f_bz_3eq_reac, build_bz_3eq_diff_matrix, ntot, u, nt, tini, tend) 
      end if 
    case ("imex_rk_2_2")
      print *
      print *, "Second order (with 2 stages) RK IMEX integration"
      if (neq == 2) then
        call second_order_2_stages_rk_imex_integration(f_bz_2eq_reac, f_bz_2eq_diff, build_bz_2eq_diff_matrix, &
                                                       ntot, u, nt, tini, tend) 
      else
        call second_order_2_stages_rk_imex_integration(f_bz_3eq_reac, f_bz_3eq_diff, build_bz_3eq_diff_matrix, &
                                                       ntot, u, nt, tini, tend) 
      end if 
    case ("imex_rk_2_3")
      print *
      print *, "Second order (with 3 stages) RK IMEX integration"
      if (neq == 2) then
        call second_order_3_stages_rk_imex_integration(f_bz_2eq_reac, f_bz_2eq_diff, build_bz_2eq_diff_matrix, &
                                                       ntot, u, nt, tini, tend) 
      else
        call second_order_3_stages_rk_imex_integration(f_bz_3eq_reac, f_bz_3eq_diff, build_bz_3eq_diff_matrix, &
                                                       ntot, u, nt, tini, tend) 
      end if 
    case default
      print *
      print *, method, "Unknown integration method"
      call exit(0)
  end select

end subroutine integrate

end module mod_integration
