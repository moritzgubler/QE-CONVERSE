# Tests

## Structure

```
tests/
├── Makefile              # top-level test driver
├── integration/          # physics integration tests
│   ├── CO+/              # EPR g-tensor of CO+ radical
│   ├── nacl_nmr/         # NMR chemical shifts of NaCl (PBEsol)
│   ├── quartz/           # NMR chemical shifts of α-SiO₂
│   ├── test_pseudos/     # pseudopotential files shared by all tests
│   ├── check_gtensor.py  # output parser and checker for g-tensor tests
│   ├── check_nmr.py      # output parser and checker for NMR tests
│   └── parse_output.py   # shared output parsing library
└── unit/                 # Fortran unit tests (testdrive framework)
```

## Prerequisites

- `qe-converse.x` built and installed to `bin/` (run `make` in the project root)
- `pw.x` from a Quantum ESPRESSO 7.2 build
- Python 3 (for integration test checkers)
- MPI launcher (e.g. `mpirun`) if running parallel tests

## Running the tests

From the project root:

```bash
make test
```

## Integration tests

Each test subdirectory contains:

| File | Purpose |
|---|---|
| `run.sh` | Runs pw.x SCF, then qe-converse, then the checker script |
| `*.in` | Input files for qe-converse |
| `pw_scf.in` / `pw.scf.in` | Input file for the pw.x ground-state calculation |
| `reference/` | Reference input and output files for comparison |

The checker scripts (`check_gtensor.py`, `check_nmr.py`) parse key physical
quantities from the output and compare them to the reference values. They exit
with a non-zero status if any value exceeds its tolerance, causing the test to
fail.

## Unit tests

The unit tests use a copy of the
[testdrive](https://github.com/fortran-lang/testdrive) framework (Apache-2.0 OR MIT).
The project must be configured (`./configure --with-qe-source=...`) before
building the unit tests, so that the required QE modules are available.

### Currently tested

**Suite `util`** (`test_util.f90`) — routines from `src/util.f90`:

| Test | What it checks |
|---|---|
| `trace_identity` | `trace` of the 3×3 identity matrix equals 3 |
| `trace_diagonal` | `trace` of diag(1,2,3) equals 6 |
| `trace_zero_matrix` | `trace` of the zero matrix equals 0 |
| `principal_axis_isotropic` | `principal_axis` on c·I returns all eigenvalues = c |
| `principal_axis_diagonal` | `principal_axis` on diag(1,2,3) returns eigenvalues [1,2,3] and orthonormal eigenvectors |
| `principal_axis_simpson` | `principal_axis_simpson` sorts eigenvalues by \|λ − Tr/3\| (Simpson convention) |

### Adding a new test

**New test in an existing module** — add a subroutine to the relevant `test_*.f90`
and register it in the suite's `collect_*_tests` array:

```fortran
subroutine test_my_routine(error)
  type(error_type), allocatable, intent(out) :: error
  ! ... set up inputs ...
  call check(error, actual_value, expected_value, thr=1.0e-10_dp)
end subroutine
```

Then in `collect_util_tests` (or whichever suite applies):

```fortran
testsuite = [ &
  ..., &
  new_unittest("my_routine", test_my_routine) &
]
```

**New module** — create `test_mymodule.f90` with a `collect_mymodule_tests`
subroutine, add a `new_testsuite(...)` entry in `test_runner.f90`, and add the
new source file to `SRCS` in `unit/Makefile`.

`check` accepts scalars, arrays, and logical values; see `testdrive.F90` for the
full interface.
