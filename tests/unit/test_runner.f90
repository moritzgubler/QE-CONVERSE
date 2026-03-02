! Test runner for QE-CONVERSE unit tests.
!
! Uses the test-drive framework (https://github.com/fortran-lang/test-drive).
! Build with:   make -C tests/unit
! Run  with:    tests/unit/test_runner

program test_runner
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive,  only: run_testsuite, new_testsuite, testsuite_type
  use test_util,  only: collect_util_tests
  use mp_global,  only: mp_startup
  use mp_world,   only: mp_world_end
  implicit none

  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  call mp_startup()

  stat = 0

  testsuites = [ &
    new_testsuite("util", collect_util_tests) &
  ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  call mp_world_end()

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program test_runner
