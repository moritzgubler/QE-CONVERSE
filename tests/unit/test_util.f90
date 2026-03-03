! Unit tests for src/util.f90 using the test-drive framework.
!
! Tested routines:
!   trace              – sum of diagonal elements of an n×n matrix
!   principal_axis     – diagonalise 3×3 tensor, sort eigenvalues by |λ|
!   principal_axis_simpson – same but sort by |λ − Tr/3|
!
! Compilation: handled by tests/unit/Makefile.
! These tests require the QE libraries to be built (for cdiagh, hpsort, kinds).

module test_util
  use testdrive, only: new_unittest, unittest_type, error_type, check
  implicit none
  private
  public :: collect_util_tests

  ! Match QE's definition of double precision
  integer, parameter :: dp = selected_real_kind(14, 200)

contains

  subroutine collect_util_tests(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("trace_identity",          test_trace_identity),        &
      new_unittest("trace_diagonal",          test_trace_diagonal),        &
      new_unittest("trace_zero_matrix",       test_trace_zero_matrix),     &
      new_unittest("principal_axis_isotropic",test_principal_axis_iso),    &
      new_unittest("principal_axis_diagonal", test_principal_axis_diag),   &
      new_unittest("principal_axis_simpson",  test_principal_axis_simpson) &
    ]
  end subroutine collect_util_tests


  ! ----- trace tests -------------------------------------------------------

  subroutine test_trace_identity(error)
    ! Tr(I_3) = 3
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: a(3,3), tr
    integer  :: i, j

    a = 0.0_dp
    do i = 1, 3
      a(i,i) = 1.0_dp
    end do

    call trace(3, a, tr)
    call check(error, tr, 3.0_dp, thr=1.0e-14_dp)
  end subroutine test_trace_identity


  subroutine test_trace_diagonal(error)
    ! Tr(diag(1, 2, 3)) = 6
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: a(3,3), tr

    a      = 0.0_dp
    a(1,1) = 1.0_dp
    a(2,2) = 2.0_dp
    a(3,3) = 3.0_dp

    call trace(3, a, tr)
    call check(error, tr, 6.0_dp, thr=1.0e-14_dp)
  end subroutine test_trace_diagonal


  subroutine test_trace_zero_matrix(error)
    ! Tr(0) = 0
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: a(3,3), tr

    a = 0.0_dp
    call trace(3, a, tr)
    call check(error, tr, 0.0_dp, thr=1.0e-14_dp)
  end subroutine test_trace_zero_matrix


  ! ----- principal_axis tests ----------------------------------------------

  subroutine test_principal_axis_iso(error)
    ! Isotropic tensor c*I → all eigenvalues = c, sorted arbitrarily.
    type(error_type), allocatable, intent(out) :: error
    real(dp), parameter :: c = 5.0_dp
    real(dp) :: tens(3,3), eigs(3), eigv(3,3)
    integer  :: i

    tens      = 0.0_dp
    do i = 1, 3
      tens(i,i) = c
    end do

    call principal_axis(tens, eigs, eigv)

    ! All three eigenvalues must equal c
    do i = 1, 3
      call check(error, eigs(i), c, thr=1.0e-10_dp)
      if (allocated(error)) return
    end do
  end subroutine test_principal_axis_iso


  subroutine test_principal_axis_diag(error)
    ! diag(1, 2, 3): eigenvalues sorted ascending by |λ| → [1, 2, 3].
    !
    ! principal_axis calls cdiagh (which uses LAPACK's zheev and returns
    ! eigenvalues in ascending order), then hpsort ascending by |λ|.
    ! For all-positive eigenvalues the sort is a no-op.
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: tens(3,3), eigs(3), eigv(3,3)

    tens      = 0.0_dp
    tens(1,1) = 1.0_dp
    tens(2,2) = 2.0_dp
    tens(3,3) = 3.0_dp

    call principal_axis(tens, eigs, eigv)

    call check(error, eigs(1), 1.0_dp, thr=1.0e-10_dp)
    if (allocated(error)) return
    call check(error, eigs(2), 2.0_dp, thr=1.0e-10_dp)
    if (allocated(error)) return
    call check(error, eigs(3), 3.0_dp, thr=1.0e-10_dp)
    if (allocated(error)) return

    ! Eigenvectors must be orthonormal: e_i · e_j = δ_{ij}
    call check(error, &
      abs(dot_product(eigv(:,1), eigv(:,2))) < 1.0e-10_dp, &
      "eigenvectors 1 and 2 not orthogonal")
    if (allocated(error)) return
    call check(error, &
      abs(dot_product(eigv(:,1), eigv(:,3))) < 1.0e-10_dp, &
      "eigenvectors 1 and 3 not orthogonal")
    if (allocated(error)) return
    call check(error, &
      abs(dot_product(eigv(:,2), eigv(:,3))) < 1.0e-10_dp, &
      "eigenvectors 2 and 3 not orthogonal")
  end subroutine test_principal_axis_diag


  subroutine test_principal_axis_simpson(error)
    ! diag(1, 2, 4): the Simpson convention sorts by |λ − Tr/3|.
    !   Tr = 7,  Tr/3 ≈ 2.333
    !   |λ − Tr/3|: |1−2.333|=1.333, |2−2.333|=0.333, |4−2.333|=1.667
    !   Sort ascending → order [2, 1, 4]
    !   So eigs(1)=2, eigs(2)=1, eigs(3)=4
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: tens(3,3), eigs(3), eigv(3,3)

    tens      = 0.0_dp
    tens(1,1) = 1.0_dp
    tens(2,2) = 2.0_dp
    tens(3,3) = 4.0_dp

    call principal_axis_simpson(tens, eigs, eigv)

    call check(error, eigs(1), 2.0_dp, thr=1.0e-10_dp)
    if (allocated(error)) return
    call check(error, eigs(2), 1.0_dp, thr=1.0e-10_dp)
    if (allocated(error)) return
    call check(error, eigs(3), 4.0_dp, thr=1.0e-10_dp)
  end subroutine test_principal_axis_simpson

end module test_util
