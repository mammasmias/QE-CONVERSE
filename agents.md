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

The full orbital magnetization splits into four terms:

| Term | Meaning |
|------|---------|
| `M_bare` | Berry-phase integral with the GIPAW Hamiltonian (LC + IC parts) |
| `ΔM_NL` | Non-local pseudopotential correction |
| `ΔM_para` | Paramagnetic GIPAW correction (uses `F_R^NL` for EPR, `K_R^NL` for NMR) |
| `ΔM_dia` | Diamagnetic GIPAW correction (uses `E_R^NL` for EPR, `J_R^NL` for NMR) |

The output keywords `M_LC`, `M_IC`, and `Delta_M` correspond to these terms.

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

The executable `qe-converse.x` is placed in `bin/`. Links against the pre-installed QE libraries.
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
       ├─ compute_dudk_new   covariant du/dk via finite differences
       │    └─ compute_u_kq  diagonalise H(k+q) for each k and q=±q_gipaw·x̂,ŷ,ẑ
       └─ calc_delta_M_bare  bare orbital magnetization (Berry curvature integral)
            + GIPAW correction terms (ΔM_NL, ΔM_para, ΔM_dia)
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
| `newscf.f90` | Re-runs SCF with `io_level=0`, `starting_wfc='atomic'`. Finds Fermi level for non-metallic systems. |
| `wfcinit_gipaw.f90` | Generates initial wavefunctions (atomic superposition or random); calls subspace diagonalisation to produce starting eigenstates for GIPAW SCF. |
| `electrons_gipaw.f90` | Full SCF convergence loop: charge mixing (Broyden), density updates, Hartree/XC, convergence check. |
| `c_bands_gipaw.f90` | Diagonalises H(k) for each k-point. **Must always call `save_buffer(evc, nwordwfc, iunwfc, ik)` after diagonalisation** — even when `nks=1` (k-pool parallel case). |
| `h_psi_gipaw.f90` | Applies the GIPAW Hamiltonian to wavefunctions: kinetic energy, local potential, non-local KB, spin-orbit (`H^(1,0)` SO term when `lambda_so≠0`), NMR vector potential (`H^(1,0)` dipole term when `m_0≠0`). Central routine analogous to QE's `h_psi`. |
| `rotate_wfc_gipaw.f90` | Dispatcher for wavefunction rotation (subspace diagonalisation): routes to serial or parallel, Γ-point or k-point routines. |
| `rotate_wfc_k_gipaw.f90` | Serial and parallel subspace Hamiltonian diagonalisation in the k-point basis; computes H and S matrices. |
| `compute_dudk_new.f90` | Computes covariant du/dk: calls `compute_u_kq` for each displaced k±q, builds overlap matrix, inverts, assembles dual states, saves result. Supports in-memory (`dudk_storage`) or disk I/O. After calling `compute_u_kq(ik,0)`, saves `evc` back to `iunwfc`. Forces in-memory storage when `nproc_pool > 1` (disk I/O is not safe with distributed G-vectors). |
| `converse_kdotp.f90` | Projected Sternheimer du/dk using QE `cgsolve_all` and the full GIPAW Hamiltonian; central differences of H include k-dependent Hubbard projectors. Restores k, projectors and kinetic energies. Uses the same occupied subspace as `covariant`. |
| `compute_u_kq.f90` | Diagonalises H(k+q). Reads `evc` from `iunwfc` as initial guess (line 133); after diagonalisation, restores `evc` from `iunwfc` (lines 193-196) to not corrupt the k-point buffer. |
| `dudk_storage.f90` | Module providing in-memory buffer management for du/dk matrices; allocates/deallocates storage and handles retrieval for better performance than disk I/O. |
| `calc_orbital_magnetization.f90` | Computes all orbital magnetization terms: LC, IC (Berry curvature formula), and GIPAW corrections (`delta_M_bare`, `delta_M_para`, `delta_M_dia`) for both EPR and NMR modes. Contains band-group MPI parallelisation over `ibnd`. |
| `orbital_magnetization.f90` | Module storing orbital magnetization results (`M_LC`, `M_IC`, `Delta_M`, total) and physical constants (α, g_e = 2.002319, a2gp8). |
| `gipaw_module.f90` | Global GIPAW parameters: q-vector (`q_gipaw`), convergence thresholds, SO coupling strengths (`lambda_so`), radial integrals, L-operators, projector data. |
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

The single input namelist is `&input_qeconverse`. Key parameters:

| Keyword | Type | Default | Description |
|---------|------|---------|-------------|
| `prefix` | char | — | Must match the pw.x SCF `prefix`. |
| `outdir` | char | `./` | Must match the pw.x SCF `outdir`. |
| `diagonalization` | string | `david` | Only `david` (Davidson) supported. |
| `q_gipaw` | real | 0.01 | Finite-difference step in k-space (bohr⁻¹). Keep at default. |
| `dudk_method` | string | `covariant` | `covariant`, `singlepoint`, or `kdotp` (projected Sternheimer). |
| `diago_thr_init` | real | 1e-7 | Initial diagonalisation convergence threshold (Ry²). |
| `conv_threshold` | real | 1e-8 | SCF convergence threshold (Ry²). |
| `mixing_beta` | real | 0.5 | Charge mixing factor. |
| `lambda_so(1..3)` | real | 0.0 | Spin-orbit coupling direction for EPR. Set one component ≠ 0 per run. |
| `m_0(1..3)` | real | 0.0 | Nuclear dipole moment direction for NMR. Set one component ≠ 0 per run. |
| `m_0_atom` | int | 0 | Index of the atom carrying the NMR dipole. |

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

---

## Output Quantities

| Keyword | Description |
|---------|-------------|
| `M_LC` | Local Circulation orbital magnetization (a.u.) |
| `M_IC` | Itinerant Circulation orbital magnetization (a.u.) |
| `Delta_M` | Non-local + paramagnetic + diamagnetic GIPAW correction (a.u.) |
| `M_total` | Sum of `M_LC + M_IC + Delta_M` |
| `delta_g RMC` | Relativistic mass correction to Δg (ppm) |
| `delta_g RMC (GIPAW)` | GIPAW correction to `delta_g RMC` (ppm) |
| `delta_g SO` | SO coupling contribution to Δg (ppm) |
| `delta_g tot` | Total Δg = RMC + RMC(GIPAW) + SO (ppm) |
| `Chemical shift (ppm)` | Full chemical shielding tensor |
| `Core shift (ppm)` | Core contribution to shielding |

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
