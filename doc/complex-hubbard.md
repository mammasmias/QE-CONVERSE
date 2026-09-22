# Complex Hubbard occupations in CONVERSE

This branch preserves complex collinear occupation matrices in the perturbed
CONVERSE SCF for Dudarev U and U+V, with one Hubbard manifold per atom.
Co-3d alone is supported, as are Co-3d and O-2p on different atoms.

## Implementation

`src/new_nsg.f90` is a GPL copy of `PW/src/new_nsg.f90` from QE tag `qe-7.5`.
The collinear complex-k accumulation no longer takes `DBLE` of the projection
product. Both upstream entry points (`new_nsg` and `new_nsg_nc`) are retained
in the same object: QE's `sum_band` references both, so retaining only one
would pull in the archive object and create a duplicate symbol. CONVERSE
links this object before `libpw.a`. The noncollinear implementation is
unchanged and is not enabled by this patch. The separate EFG executable
continues to link the upstream QE implementation.

`src/nsg_adj.f90` also overrides its QE 7.5 counterpart. When user-specified
`starting_ns` eigenvalues are applied, it reconstructs the full complex
matrix with the correct eigenvector conjugation and Hermitian partner.
Taking `DBLE` here, or retaining the real-case conjugation order, would lose
or conjugate the initial imaginary off-diagonal entries.

For U+V, QE already stores `nsg`, the potential `v_nsg`, and the Broyden
history as complex arrays. Its mixing norm intentionally takes the real
part of a Hermitian product; this does not discard imaginary residuals.

For U-only input, `prepare_complex_hubbard` runs after `potinit` and before
the first perturbed diagonalization. It sets `Hubbard_V(I,I,1)=Hubbard_U(type(I))`
and leaves all intersite interactions exactly zero. It calls QE's
`init_hubbard` again to initialize neighborhoods, dimensions, orbital labels,
and complex arrays, then transfers the initial on-site occupations and
rebuilds the potential. It checks that projector ordering and the initial
Hubbard energy are unchanged. This is an initialized representation change,
not just an assignment of `lda_plus_u_kind=2`.

Both inputs subsequently use QE's complex U+V SCF/mixing/potential path and
CONVERSE's existing U+V magnetization correction. In the U-only case,
self-neighbor displacements vanish, so the intersite displacement term is
zero. No uniform magnetic-field term is added to `h_psi_gipaw`.

For a Hermitian occupation matrix N the on-site functional is

```
E_U = U/2 [Tr(N) - Tr(N*N)]
Tr(N*N) = sum_ab |N_ab|^2
```

The imaginary off-diagonal matrix elements contribute to both its energy
and its derivative. Taking the real part of N first changes this functional.

## Build and input

Build against a normally compiled QE 7.5; no edits to the installed QE tree
or changes to its module ABI are required:

```sh
./configure --with-qe-source=/path/to/q-e-7.5
make clean
make -j4
```

Run `configure` again for an existing CONVERSE checkout so that the new
objects and module dependencies enter the generated Makefile.

Use `nosym=.true.` and `noinv=.true.` in the original pw.x input, with an
explicit complex k-point grid. `K_POINTS automatic` with `1 1 1 0 0 0` is
different from the optimized real `K_POINTS gamma`, which is rejected.
The original pw.x ground-state calculation is not modified by this branch;
complex occupations are updated during the subsequent CONVERSE SCF.
Restart density files from pw.x are not rewritten by this postprocessor.

Supported scope: collinear Dudarev U and U+V, one manifold per atom, and the
projectors already supported by CONVERSE. Background/additional manifolds,
Liechtenstein U, noncollinear spins, J0/beta, orbital-resolved U, and fixed
Hubbard potentials are rejected explicitly. These cases need separate
implementation and validation; they are not silently approximated.

The output includes the maximum absolute imaginary occupation and the
maximum on-site Hermiticity error after SCF. A nonzero imaginary part is
expected only when allowed by the state and perturbation; zero by itself
is not an error.

## Tests and validation status

```sh
python3 tests/test_complex_hubbard.py
PW=/path/to/pw.x QECONVERSE=/path/to/qe-converse.x \
  python3 tests/integration/check_complex_hubbard.py
make test
```

The NumPy tests check a complex pure-orbital counterexample, functional
derivatives, orbital-unitary covariance, the on-site limit of the complex
U+V operator, and mixing/reduction algebra. They do not execute QE.

The binary regression uses a CO+ fixture with one p manifold per atom,
first with U only and then with U+V. It checks retention of imaginary
occupations, Hermiticity, and agreement of energy and orbital magnetization
between serial, two-pool, and four-pool runs. Its purpose is implementation
validation, not a benchmark for a particular material. Outputs are kept
in the printed temporary directory. The dedicated GitHub Actions workflow
builds against QE 7.5 and runs these tests.

Full QE compilation and serial/MPI validation are performed by the Linux
GitHub Actions runner; the Windows editing environment only runs the
algebraic and syntax checks. Check the workflow result for the exact commit
being used. Existing Hubbard reference values must not be regenerated just
to hide discrepancies.
Passing these tests does not independently establish the completeness of
the Hubbard orbital-magnetization theory or convergence for Co compounds.
