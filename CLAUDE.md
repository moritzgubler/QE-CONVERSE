# Agent Guide for QE-CONVERSE

This file provides guidance for AI agents working on the QE-CONVERSE codebase.
It is derived from the source code and the companion publication (Fioccola et al., 2025).

---

## Project Overview

QE-CONVERSE is a standalone Fortran package for Quantum ESPRESSO (QE 7.5) that computes
**orbital magnetization** non-perturbatively via the *converse* GIPAW method. It replaces the
original converse routines from the obsolete PWscf 3.2 code.

Derived properties:
- **NMR chemical shielding tensors** (¹⁷O, ²⁹Si, ²⁷Al, …)
- **EPR g-tensor** deviations Δg
- **Electric field gradient (EFG) tensors** and NMR/NQR quadrupolar parameters (`Cq`, `η`, `ν_Q`),
  via the standalone `qe-efg.x` driver (a ground-state property, no converse SCF needed)

Scope: isolated molecules and periodic solids. Functionals: LDA, GGA. Pseudopotentials: GIPAW
norm-conserving (Trouiller-Martins + GIPAW reconstruction required).

---

## Physics Background

### Orbital Magnetization (Modern Theory)

In periodic systems the position operator is ill-defined, so orbital magnetization is computed
via the Berry-phase formula in the Wannier representation:

```
M = -(α N_c / 2 N_k) Im Σ_{n,k} f_{n,k}
      <∂_k u_{n,k}| × (H_k + ε_{n,k} - 2ε_F) |∂_k u_{n,k}>
```

where `H_k` is the zero-field GIPAW Hamiltonian, `ε_{n,k}` and `u_{n,k}` its eigenvalues
and eigenvectors, `ε_F` the Fermi level, and `α = 1/c` the fine-structure constant.

### Total Orbital Magnetization (GIPAW decomposition)

The full orbital magnetization splits into five terms:

| Term | Meaning |
|------|---------|
| `M_bare` | Berry-phase integral with the GIPAW Hamiltonian (LC + IC parts) |
| `ΔM_NL` | Non-local pseudopotential correction |
| `ΔM_para` | Paramagnetic GIPAW correction (uses `F_R^NL` for EPR, `K_R^NL` for NMR) |
| `ΔM_dia` | Diamagnetic GIPAW correction (uses `E_R^NL` for EPR, `J_R^NL` for NMR) |
| `ΔM_Hub` | Hubbard U correction (active when `lda_plus_u=.true.` and `lhub_magnetization=.true.`) |

`M_total = M_LC + M_IC + ΔM_NL + ΔM_para + ΔM_dia + ΔM_Hub`. The output keyword `Delta_M`
groups `ΔM_NL + ΔM_para + ΔM_dia + ΔM_Hub`; `Delta_M_hubbard` is printed separately.

### Hubbard U Contribution to Orbital Magnetization (DFT+U)

The Hubbard potential on site `I` is separable:

```
V_I^Hub = Σ_{mm'} |φ_{Im}^k> V_{mm'}^I <φ_{Im'}^k|
```

where `V_{mm'}^I = U^I [δ_{mm'}/2 - n_{mm'}^{I,σ}]` is the Hubbard potential matrix
(stored as `v%ns(m,m',σ,na)` in QE) and `φ_{Im}` are the Hubbard projector orbitals
(atomic or orthogonalized atomic; `pseudo`-type is excluded).

Starting from the same GIPAW commutator expression as for the bare non-local term, one
arrives at an identical structure with KB projectors replaced by Hubbard projectors:

```
ΔM_Hub,γ = -(α/2) Σ_{n,k} f_{nk} Σ_I Σ_{mm'} Σ_{μν} ε_{γμν}
              Im[ <u_{nk}|∂_{k_μ} φ̃_{Im}^k> V_{mm'}^I <∂_{k_ν} φ̃_{Im'}^k|u_{nk}> ]
```

**No-phase convention for Hubbard projectors.** The tilde `φ̃` denotes the no-phase
projector, defined by stripping the global atomic-center phase `e^{-ik·R_I}` from the
full-phase projector produced by `orthoUwfc_k`. In reciprocal space:

```
φ̃^{k,I}_m(G) = e^{-iG·R_I} φ_m(k+G)        (no e^{-ik·R_I} factor)
```

This ensures that `∂_{k_μ} φ̃^{k,I}_m = -i(r_μ - R_{I,μ}) φ̃^{k,I}_m`, i.e., the
derivative generates the relative position operator from the atomic center, not the
absolute position. Without stripping the phase, the derivative would produce a spurious
term proportional to `R_{I,μ}`.

**Implementation in `calc_orbital_magnetization.f90`:**
- `compute_dhubbecp` — computes `dhubbecp(m,n,ipol) = d/dk_ipol <φ̃_{Im}^k | u_{nk}>`
  via central finite differences. Calls `orthoUwfc_k` at `k±δk`, then multiplies by the
  phase `e^{+i(k±δk)·R_I}` to strip the Bloch center phase, then calls `ZGEMM` to
  form the overlap with `evc`.
- `calc_delta_M_hub` — evaluates the double sum over `m,m'` and bands, contracting
  `conjg(dhubbecp(m,n,ii)) * dhubbecp(m',n,jj)` with `v%ns`, and accumulates
  `delta_M_hub(kk) -= 2·Im(tmp)` (same sign convention as `delta_M_bare`).

**Supported Hubbard projector types:** `atomic` and `ortho-atomic`. The `pseudo` type
is explicitly excluded (Kleinman-Bylander projectors already enter via `ΔM_NL`).

**DFT+U+V (inter-site interactions):** **not yet implemented**. The current code handles
on-site `U` only.

### EPR g-tensor

The GIPAW Hamiltonian at zero external field is:

- `H^(0,0)` = kinetic + local potential + KB non-local pseudopotential (Eq. 4)
- `H^(1,0)` = spin-orbit (SO) term + paramagnetic non-local `F_R^NL` (Eq. 5); activated by
  `lambda_so ≠ 0`

The g-tensor deviation from the free-electron value is:

```
Δg_{μν} = -(2 / α S) e_μ · M(e_ν)
```

where `S` is the total spin and `e_ν` is one of the three spin-axis directions.
Three independent SCF+magnetization runs are needed (one per Cartesian direction of the spin axis).

### NMR Chemical Shift

For NMR, an additional vector potential from a magnetic point dipole `m_s` at atom `s` is added
to the Hamiltonian (Eq. 17-26). The `H^(1,0)` term gains the vector-potential kinetic term plus
the non-local paramagnetic operator `K_R^NL` (activated by `m_0 ≠ 0`).

The chemical shielding tensor is:

```
σ_{s,αβ} = δ_{αβ} - Ω ∂M_β/∂m_s
```

where `Ω` is the cell volume. Again three runs are needed (one per dipole direction).

### Electric Field Gradient (EFG)

The EFG is the traceless rank-2 tensor `V_{αβ}` (the l=2 multipole of the electrostatic
potential) evaluated at each nucleus. Unlike the orbital-magnetization properties above, it is
a **pure ground-state property**: it is computed directly from the converged SCF charge density
and wavefunctions (`read_file` → `gipaw_setup` → `calc_efg`), with **no converse SCF**, no
`q_gipaw` k-derivative, and no `lambda_so`/`m_0` perturbation. A single SCF run gives the full
tensor at every site.

The total EFG is the sum of three contributions (all printed separately by `calc_efg`):

| Term | Meaning |
|------|---------|
| `efg_bare`  | Electronic valence contribution from the smooth pseudo charge density, G-space sum `V_{αβ}(R) = (4π/Ω) Σ_{G≠0} ρ(G)(G_α G_β/G² − δ_{αβ}/3) e^{iG·R}/G²` |
| `efg_ion`   | Ionic (point-nucleus) contribution from the surrounding lattice, via Ewald sum (`ewald_dipole`); self-term excluded |
| `efg_gipaw` | GIPAW/PAW reconstruction: on-site AE−PS `1/r³` radial integrals contracted with the l=2 Gaunt coefficients `ap(5..9, lm1, lm2)` — restores the aspherical near-nucleus density the pseudo density gets wrong |

`efg_tot = efg_bare + efg_ion + efg_gipaw`, then symmetrized with `symtensor` (which uses the
crystal symmetry and ties symmetry-equivalent atoms together).

**None of the three terms is generically small or negligible** — the ranking depends on
ionicity and site symmetry. For ionic, low-symmetry sites (e.g. O in α-quartz) `efg_ion` can
dominate; the GIPAW term is typically ~2× the bare term because of the `1/r³` weighting near the
nucleus. The result is **k-point insensitive** (it is a local charge-density property): a coarse
grid converges it; residual error vs. experiment comes from the pseudopotential/XC/geometry.

Spectroscopic parameters from the principal values (`|Vzz| ≥ |Vyy| ≥ |Vxx|`):

```
Cq   = e Q V_zz / h            (MHz; quadrupolar coupling constant)
η    = (V_xx − V_yy) / V_zz    (asymmetry)
ν_Q  = 3 Cq / (2I(2I−1))       (MHz; quadrupolar frequency, needs I ≥ 1)
```

`Q` (nuclear quadrupole moment) and `I` (nuclear spin) are inputs per atom type (`q_efg`,
`i_efg`). `Cq` is computed only when `Q ≠ 0`; `ν_Q` only when `I ≥ 1`.

**Symmetry note:** the EFG path does **not** require `nosym=.true.`/`noinv=.true.`/`nspin=2`.
A symmetry-reduced k-grid is fine and faster (`symtensor` handles symmetrization). Do not reuse
an EFG SCF for a converse EPR/NMR run, which still needs the no-symmetry, spin-polarized setup.

---

## Repository Layout

```
src/                    Core Fortran source files
doc/                    User manual
examples/               NMR and EPR example inputs (with Tutorial.wiki files)
applications/           Application-level EPR and NMR calculations
benchmarking/           Input files for EPR g-tensor benchmarks
tests/                  Integration and unit tests
```

The QE 7.5 source tree is at `/home/moritz/src/q-e/`. Do not modify files there without
explicit instruction.

---

## Build System

```bash
./configure --with-qe-source=<path-to-QE-make.inc>
make
```

Two executables are placed in `bin/`: `qe-converse.x` (orbital magnetization / NMR / EPR) and
`qe-efg.x` (electric field gradient). Both link against the pre-installed QE libraries.
Optional: scaLAPACK or ELPA for improved linear algebra performance.

---

## Code Architecture

### Workflow

```
Pre-requisite: pw.x SCF with nosym=.true., noinv=.true., nspin=2 (for EPR)

qe-converse.x
  ├─ read_file          read ground-state data from QE XML output
  ├─ gipaw_setup        compute GIPAW projectors and integrals
  ├─ init_nmr           build NMR vector potential A_s in G-space and real space
  ├─ newscf             re-run SCF with GIPAW perturbation (SO or dipole field)
  │    ├─ wfcinit_gipaw      initialise wavefunctions
  │    └─ electrons_gipaw    SCF loop (Davidson diagonalisation)
  │         └─ c_bands_gipaw diagonalise H(k) for each k-point
  └─ calc_orbital_magnetization
       ├─ compute_dudk_new    covariant du/dk via finite differences
       │    └─ compute_u_kq   diagonalise H(k+q) for each k and q=±q_gipaw·x̂,ŷ,ẑ
       ├─ calc_delta_M_bare   bare orbital magnetization (Berry curvature integral)
       │    + GIPAW correction terms (ΔM_NL, ΔM_para, ΔM_dia)
       ├─ compute_dhubbecp    k-derivative of no-phase Hubbard projectors (if DFT+U)
       └─ calc_delta_M_hub    Hubbard U correction to orbital magnetization (if DFT+U)
```

The EFG driver is much shorter — it is a single-shot ground-state post-processing step
(no converse SCF):

```
Pre-requisite: ordinary pw.x SCF (symmetry allowed; no nosym/noinv/nspin=2 needed)

qe-efg.x
  ├─ read_file            read ground-state data from QE XML output
  ├─ gipaw_setup          compute GIPAW projectors and radial integrals
  ├─ calc_efg             EFG tensor at each site = efg_bare + efg_ion + efg_gipaw
  │    ├─ efg_bare_el        valence electronic term (G-space sum over pseudo ρ)
  │    ├─ ewald_dipole       ionic term (Ewald sum over nuclear point charges)
  │    ├─ efg_correction     GIPAW PAW reconstruction (on-site 1/r³, l=2 Gaunt)
  │    └─ symtensor + principal_axis → Vxx,Vyy,Vzz, Cq, η, ν_Q
  └─ print_efg_summary    Cq/η/ν_Q table for all atoms
```

### Covariant Finite Difference for du/dk

The k-derivative of Bloch functions uses gauge-invariant "dual" states:

```
|∂̃_i u_{n,k}> = ½ (|ũ_{n,k+q}> - |ũ_{n,k-q}>)

|ũ_{n,k+q}> = Σ_{n'} (S_{k+q}^{-1})_{n'n} |u_{n',k+q}>

(S_{k+q})_{nn'} = <u_{n,k}|u_{n',k+q}>
```

Loop structure (implemented in `compute_dudk_new.f90` + `compute_u_kq.f90`):
1. For each k-point: read `|u_{n,k}>` from disk.
2. For each of the 6 displaced points ±q·x̂, ±q·ŷ, ±q·ẑ: diagonalise H(k+q), compute
   overlap S, invert, build dual state, compute finite difference.
3. Save du/dk to disk (or in-memory buffer via `dudk_storage`).

### Orbital Magnetization Parallelisation

`calc_orbital_magnetization` contains a new **band-group parallelisation** (vs the original
PWscf 3.2 code): the inner `ibnd` loop is split across MPI processors via the QE band-group
communicator. This improves scalability for large supercell calculations.

---

## Source File Summary

| File | Description |
|------|-------------|
| `qe-converse.f90` | Main program. Parses input namelist `&input_qeconverse`, broadcasts to MPI ranks, calls `read_file` → `gipaw_setup` → `init_nmr` → `newscf` → `calc_orbital_magnetization`, prints timing. |
| `qe-efg.f90` | Standalone EFG driver program. Parses namelist `&input_qeefg` (`prefix`, `outdir`, `q_efg`, `i_efg`), then `read_file` → `gipaw_setup` → `calc_efg` → `print_efg_summary`. No converse SCF. Built as `bin/qe-efg.x`. |
| `efg.f90` | EFG implementation. `calc_efg` assembles `efg_bare` + `efg_ion` + `efg_gipaw`, symmetrizes, diagonalises, prints `Cq`/`η`/`ν_Q`. Internal routines: `efg_bare_el` (G-space valence term), `efg_correction` (GIPAW PAW reconstruction, l=2 Gaunt), `print_efg_summary`. Uses QE `ewald_dipole` for the ionic term. Adapted from qe-gipaw (Ceresoli et al.). |
| `efg_module.f90` | Module `efg_mod`: stores the symmetrized `efg_tensor(3,3,nat)` for end-of-run summary printing. |
| `newscf.f90` | Re-runs SCF with `io_level=0`, `starting_wfc='atomic'`. Finds Fermi level for non-metallic systems. |
| `wfcinit_gipaw.f90` | Generates initial wavefunctions (atomic superposition or random); calls subspace diagonalisation to produce starting eigenstates for GIPAW SCF. |
| `electrons_gipaw.f90` | Full SCF convergence loop: charge mixing (Broyden), density updates, Hartree/XC, convergence check. |
| `c_bands_gipaw.f90` | Diagonalises H(k) for each k-point. **Must always call `save_buffer(evc, nwordwfc, iunwfc, ik)` after diagonalisation** — even when `nks=1` (k-pool parallel case). |
| `h_psi_gipaw.f90` | Applies the GIPAW Hamiltonian to wavefunctions: kinetic energy, local potential, non-local KB, spin-orbit (`H^(1,0)` SO term when `lambda_so≠0`), NMR vector potential (`H^(1,0)` dipole term when `m_0≠0`). Central routine analogous to QE's `h_psi`. |
| `rotate_wfc_gipaw.f90` | Dispatcher for wavefunction rotation (subspace diagonalisation): routes to serial or parallel, Γ-point or k-point routines. |
| `rotate_wfc_k_gipaw.f90` | Serial and parallel subspace Hamiltonian diagonalisation in the k-point basis; computes H and S matrices. |
| `compute_dudk_new.f90` | Computes covariant du/dk: calls `compute_u_kq` for each displaced k±q, builds overlap matrix, inverts, assembles dual states, saves result. Supports in-memory (`dudk_storage`) or disk I/O. After calling `compute_u_kq(ik,0)`, saves `evc` back to `iunwfc`. Forces in-memory storage when `nproc_pool > 1` (disk I/O is not safe with distributed G-vectors). |
| `compute_u_kq.f90` | Diagonalises H(k+q). Reads `evc` from `iunwfc` as initial guess (line 133); after diagonalisation, restores `evc` from `iunwfc` (lines 193-196) to not corrupt the k-point buffer. |
| `dudk_storage.f90` | Module providing in-memory buffer management for du/dk matrices; allocates/deallocates storage and handles retrieval for better performance than disk I/O. |
| `calc_orbital_magnetization.f90` | Computes all orbital magnetization terms: LC, IC (Berry curvature formula), and GIPAW corrections (`delta_M_bare`, `delta_M_para`, `delta_M_dia`, `delta_M_hub`) for both EPR and NMR modes. Contains internal subroutines `compute_dhubbecp` and `calc_delta_M_hub` for the DFT+U contribution. Contains band-group MPI parallelisation over `ibnd`. |
| `orbital_magnetization.f90` | Module storing orbital magnetization results (`M_LC`, `M_IC`, `Delta_M`, `delta_M_hub`, total) and physical constants (α, g_e = 2.002319, a2gp8). |
| `gipaw_module.f90` | Global GIPAW parameters: q-vector (`q_gipaw`), convergence thresholds, SO coupling strengths (`lambda_so`), radial integrals, L-operators, projector data. Also holds `lhub_magnetization` flag. |
| `gipaw_setup.f90` | Initialises GIPAW infrastructure: reads UPF pseudopotential data, sets up radial integrals for NMR/EPR (`f_Rnm`, `e_Rnm`, `k_Rnm`, `j_Rnm`), L-operator matrices, spline interpolations. |
| `init_gipaw_1.f90` | Constructs GIPAW projectors `ρ̃_{R,n}` from UPF data; computes projection coefficients and orthogonalisation matrices. |
| `init_gipaw_2.f90` | Computes GIPAW projector coefficients in reciprocal space with structure factors, phase factors, and spherical harmonics. |
| `init_us_2_no_phase.f90` | Computes Kleinman-Bylander projectors in G-space without phase factors (required for k-derivative calculations, where phase is handled separately). |
| `paw_gipaw.f90` | PAW reconstruction module: defines projector types, stores AE/PS partial waves (`φ_{R,n}`, `φ̃_{R,n}`), radial integrals, and UPF data readers. |
| `nmr_module.f90` | Defines NMR parameters: magnetic dipole location (`m_0_atom`) and direction (`m_0`), fine-structure constant, vector potential arrays `A_s(r)`. |
| `nmr_routine.f90` | `init_nmr`: builds NMR vector potential `A_s` in G-space and real space. Runs on all MPI ranks (no distribution issue). Computes induced current and chemical shift contributions from orbital magnetization. |
| `set_dvrs.f90` | Calculates the gradient of the local Kohn-Sham potential in real space via FFT (∇V_loc, needed for the SO coupling term in `h_psi_gipaw`). |
| `util.f90` | Utilities: tensor diagonalisation, spin component selection, matrix operations (trace, kinetic energy, radial integrals), spherical averaging, coordinate rotation. |
| `stop_code.f90` | Minimal cleanup: synchronises MPI processes and calls `stop`. |

---

## MPI / k-Pool Parallelism

QE-CONVERSE supports `npool > 1` (k-point pools). Each pool handles a subset of k-points.
When `npool > 1`, a single pool may hold only one k-point (`nks = 1`).

**Critical invariant:** `c_bands_gipaw` must **always** call
`save_buffer(evc, nwordwfc, iunwfc, ik)` after diagonalisation — even when `nks = 1`.
Skipping this causes `compute_u_kq` to restore stale ground-state wavefunctions from the `.wfc`
disk file, which zeros out `delta_M_bare`.

**du/dk with intra-pool parallelism:** When `nproc_pool > 1`, G-vectors are distributed across
MPI processes within a pool. `davcio` cannot safely write distributed data to a shared file
record. `compute_dudk_new` therefore forces `dudk_in_memory = .true.` when `nproc_pool > 1`.

**npw per k-point:** `wvfct%npw` must be set to `ngk(ik)` at the start of each k-point
iteration in `compute_dudk_new`. Neither `compute_u_kq` nor `diag_bands_gipaw` update the
global `wvfct%npw` (they use local variables). Failing to set it causes out-of-bounds reads
of `evc`/`evq` for k-points where `ngk(ik) ≠ ngk(nks)`.

### Buffer / I/O Conventions

- `io_level = 0`: in-memory buffers (`buiol` in `buffers.f90`).
- `get_buffer` on a virgin in-memory slot falls back to reading the `.wfc` disk file via `davcio`.
- `open_buffer(..., io_level=0, ...)` calls `diropn` with `recl=-1` (existence check only; does
  not connect the unit for direct-access I/O).
- `qe-converse.f90` sets `io_level=1` before `read_file`; `newscf.f90` resets it to `io_level=0`.

---

## Testing

### Directory layout

```
tests/
  Makefile                      top-level test runner (unit + integration)
  integration/
    check_gtensor.py            output parser and checker for EPR g-tensor tests
    check_nmr.py                output parser and checker for NMR shift tests
    parse_output.py             regex parser for QE-CONVERSE stdout
    CO+/                        serial EPR g-tensor test (CO+ molecule)
    CO+_kpools/                 parallel EPR g-tensor test (CO+, npool=2, 4 MPI ranks)
    nacl_nmr/                   serial NMR chemical shift test (NaCl)
    nacl_nmr_kpool/             parallel NMR chemical shift test (NaCl, npool=4, 8 MPI ranks)
    quartz/                     serial NMR chemical shift test (quartz Si, O)
    licoo2_U/                   parallel DFT+U NMR test (LiCoO₂, npool=9, 9 MPI ranks; checks Li and Co shieldings)
    <test>/
      run.sh                    runs pw.x SCF + qe-converse.x, then calls check script
      reference/                stored reference output files (generated with serial run.sh)
        run.sh                  regenerates reference values (run locally when physics changes)
  unit/
    Makefile                    downloads test-drive, builds and runs unit tests
    test_util.f90               test-drive test module (trace, principal_axis)
    test_runner.f90             test-drive runner program
    testdrive.F90               auto-downloaded from fortran-lang/test-drive
```

### Integration tests

**What is tested:** the full binary end-to-end against stored reference output values.

Each test directory has a `run.sh` that:
1. Runs `pw.x` SCF.
2. Runs `qe-converse.x` (serially or with `mpirun --oversubscribe`).
3. Calls `check_gtensor.py` or `check_nmr.py` to compare against `reference/`.

**Environment variables required:**
```bash
export QECONVERSE=/path/to/bin/qe-converse.x
export PW=pw.x   # or full path to pw.x
```

**Run all integration tests:**
```bash
make -C tests integration
# or individually:
cd tests/integration/CO+ && bash run.sh
```

**Parallel tests** use `mpirun --oversubscribe` so they run on CI machines with fewer cores
than MPI ranks requested.

**Regenerating reference values:** run `bash run.sh` inside the `reference/` subdirectory
of a test using the serial binary. Do this whenever a physics-changing bug fix alters results.

**Tolerances:** set per-test in the check scripts (typically `atol=1.0 ppm` for g-tensor,
`atol=2.0 ppm` for NMR shift).

### Unit tests (test-drive)

**Framework:** [test-drive](https://github.com/fortran-lang/test-drive) (single-file Fortran)
**What is tested:** pure utility routines in `src/util.f90` that do not require
a full QE initialisation at runtime.

Current test coverage: `trace`, `principal_axis`, `principal_axis_simpson`.

Build and run:
```bash
make -C tests/unit        # downloads testdrive.F90, compiles, links
make -C tests/unit run    # also executes
```

Note: the unit test binary **does** link against all QE libraries (required
because `util.f90` references QE modules at link time). It does NOT require
MPI and can be run as a plain serial executable.

### Adding new tests

**New integration test case:**
1. Create `tests/integration/<name>/` with `run.sh`, input files, and `reference/run.sh`.
2. Run `bash reference/run.sh` to generate reference output.
3. Add the test to `INTEGRATION_TESTS` in `tests/Makefile`.

**New unit test:**
1. Add a `subroutine test_<name>(error)` to `tests/unit/test_util.f90`.
2. Register it in `collect_util_tests`.

---

## Input Namelist

There are two input namelists: `&input_qeconverse` (orbital magnetization / NMR shift / EPR
g-tensor, driver `qe-converse.x`) and `&input_qeefg` (EFG, driver `qe-efg.x`, see below).

### `&input_qeconverse`

Key parameters:

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `prefix` | char | — | Must match the pw.x SCF `prefix`. |
| `outdir` | char | `./` | Must match the pw.x SCF `outdir`. |
| `diagonalization` | string | `david` | Only `david` (Davidson) supported. |
| `q_gipaw` | real | 0.01 | Finite-difference step in k-space (bohr⁻¹). Keep at default. |
| `dudk_method` | string | `covariant` | Only `covariant` supported. |
| `diago_thr_init` | real | 1e-7 | Initial diagonalisation convergence threshold (Ry²). |
| `conv_threshold` | real | 1e-8 | SCF convergence threshold (Ry²). |
| `mixing_beta` | real | 0.5 | Charge mixing factor. |
| `lambda_so(1..3)` | real | 0.0 | Spin-orbit coupling direction for EPR. Set one component ≠ 0 per run. |
| `m_0(1..3)` | real | 0.0 | Nuclear dipole moment direction for NMR. Set one component ≠ 0 per run. |
| `m_0_atom` | int | 0 | Index of the atom carrying the NMR dipole. |
| `lhub_magnetization` | logical | `.true.` | Compute DFT+U Hubbard contribution to orbital magnetization. Only active when `lda_plus_u=.true.` in the SCF and projector type is not `pseudo`. |

EPR example:
```fortran
&input_qeconverse
    prefix = 'example_EPR', outdir = './scratch/'
    q_gipaw = 0.01, dudk_method = 'covariant'
    conv_threshold = 1e-8, mixing_beta = 0.5
    lambda_so(1) = 1.0   ! run once per direction (1, 2, 3)
/
```

NMR example:
```fortran
&input_qeconverse
    prefix = 'example_NMR', outdir = './scratch/'
    q_gipaw = 0.01, dudk_method = 'covariant'
    conv_threshold = 1e-9, mixing_beta = 0.5
    m_0(1) = 1.0, m_0_atom = 3   ! run once per direction (1, 2, 3)
/
```

### `&input_qeefg`

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `prefix` | char | `'prefix'` | Must match the pw.x SCF `prefix`. |
| `outdir` | char | `$ESPRESSO_TMPDIR` or `./` | Must match the pw.x SCF `outdir`. |
| `q_efg(ntyp)` | real | 0.0 | Nuclear quadrupole moment per atom type, in units of `1e-30 m²` (= 10 mbarn). 1 barn = 100 of these. `Cq` is printed only for types with `q_efg ≠ 0`. |
| `i_efg(ntyp)` | real | 0.0 | Nuclear spin `I` per atom type. `ν_Q` is printed only for types with `I ≥ 1`. |

Atom types are indexed as in the pw.x `ATOMIC_SPECIES` block. A single ordinary SCF (symmetry
allowed) is the only prerequisite; one `qe-efg.x` run gives the EFG at every site.

EFG example (α-quartz: type 1 = Si, type 2 = O):
```fortran
&input_qeefg
    prefix = 'quartz', outdir = './scratch/'
    q_efg(1) = 0.0,    i_efg(1) = 0.5    ! 29Si: I=1/2, no quadrupole
    q_efg(2) = -2.558, i_efg(2) = 2.5    ! 17O:  Q = -0.02558 b, I=5/2
/
```

---

## Output Quantities

| Keyword | Description |
|---------|-------------|
| `M_LC` | Local Circulation orbital magnetization (a.u.) |
| `M_IC` | Itinerant Circulation orbital magnetization (a.u.) |
| `Delta_M` | Sum of non-local, paramagnetic, diamagnetic, and Hubbard GIPAW corrections (a.u.) |
| `Delta_M_hubbard` | Hubbard U contribution to `Delta_M` (a.u.), printed separately for diagnostics |
| `M_total` | Sum of `M_LC + M_IC + Delta_M` (includes Hubbard if active) |
| `delta_g RMC` | Relativistic mass correction to Δg (ppm) |
| `delta_g RMC (GIPAW)` | GIPAW correction to `delta_g RMC` (ppm) |
| `delta_g SO` | SO coupling contribution to Δg (ppm) |
| `delta_g tot` | Total Δg = RMC + RMC(GIPAW) + SO (ppm) |
| `Chemical shift (ppm)` | Full chemical shielding tensor |
| `Core shift (ppm)` | Core contribution to shielding |

EFG outputs (`qe-efg.x`):

| Keyword | Description |
|---------|-------------|
| `bare electronic EFG` | Valence electronic EFG tensor (Ha/bohr²) |
| `ionic EFG` | Ionic (Ewald) EFG tensor (Ha/bohr²) |
| `GIPAW correction` | PAW reconstruction EFG tensor (Ha/bohr²) |
| `total EFG (symmetrized)` | Sum of the three terms, symmetrized (Ha/bohr²) |
| `Vxx, Vyy, Vzz` | EFG principal values, ordered `|Vzz| ≥ |Vyy| ≥ |Vxx|`, with eigenvectors |
| `Cq` | Quadrupolar coupling constant `e Q V_zz / h` (MHz) |
| `eta` | Asymmetry parameter `(Vxx − Vyy)/Vzz` |
| `nu_Q` | Quadrupolar frequency `3 Cq / (2I(2I−1))` (MHz; printed for `I ≥ 1`) |

---

## Development Conventions

- **Language:** Fortran 90/95. Follow surrounding style: `implicit none`, `intent` declarations.
- **MPI:** Use QE's internal wrappers (`mp_sum`, `mp_bcast`, etc.). Never call raw MPI directly.
- **Testing:** Always run both serial (`npool=1`) and parallel (`npool=2`, `npool=4`) to catch
  pool-parallelism bugs. Use `make -C tests integration` to run the full test suite.
- **Do not commit** compiled artefacts (`.o`, `.mod`, `.x`).
- **External dependency:** QE 7.5 API. Avoid relying on QE internal routines outside the stable
  public interface.
- **Pseudopotentials:** GIPAW norm-conserving only. Library at
  https://sites.google.com/site/dceresoli/pseudopotentials

---

## Useful QE Source Locations

| Path | Contents |
|------|----------|
| `/home/moritz/src/q-e/PW/src/buffers.f90` | Buffer I/O: `get_buffer`, `save_buffer`, `open_buffer`. |
| `/home/moritz/src/q-e/PW/src/` | Main PW source: diagonalisation, SCF, k-point loops. |
| `/home/moritz/src/q-e/Modules/` | Shared QE modules: parallelism, I/O, units. |
