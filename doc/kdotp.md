# k·p / Sternheimer derivatives

Select `dudk_method = 'kdotp'`. The default remains `covariant`.

The legacy velocity/Green-function algorithm is implemented in
`src/converse_kdotp.f90` using QE 7.5 interfaces. For every k point and
Cartesian direction it solves

```
(H_k - epsilon_n + alpha P) du_n = -Q (dH_k/dk) u_n,
P = sum_{j=1}^{occ} |u_j><u_j|, Q = 1-P.
```

`alpha` is at least 1 Ry and twice the energy span of the selected bands.
It shifts only the projected subspace. The solution is projected back onto Q.
The occupied-band count is the same as for `covariant`; with smearing this
means the existing energy cutoff in `find_nbnd_occ`, not a new metallic
occupation-weighted response prescription. A well-separated selected subspace
is needed. Failure to converge is an error, never a silent partial result.

The reference states are refreshed once at k. The right-hand side uses the
central difference of **the full GIPAW Hamiltonian**, applied to those states,
at k ± q. No displaced eigenvector diagonalizations are needed. This follows
the Hamiltonian-difference velocity approach already sketched in the legacy
source, and includes kinetic, nonlocal KB, SO, NMR and Hubbard terms. Hubbard
projectors and intersite phase factors are rebuilt at displaced k points.
The reference k, kinetic energy and projectors are restored before the solver
and on return. Existing SCF occupations and potentials are held fixed.

All energies and the Hamiltonian derivative are in Rydberg units; the old
Hartree conversion factors cancel and are not applied. The derivative is with
respect to Cartesian k in inverse bohr. `q_gipaw/tpiba` converts the physical
step to QE coordinates.

This is a Sternheimer solver with a finite-difference Hamiltonian derivative,
not an analytical velocity implementation. Reduce `q_gipaw` to verify step
convergence (the central difference has second-order truncation error).
`q_gipaw` and `conv_threshold` must be positive. The CG residual tolerance is
`min(sqrt(conv_threshold), 1e-8)`. QE's solver stops after at most 400 iterations.

Supported representation: complex collinear wavefunctions (`npol=1`) with
norm-conserving pseudopotentials, as used by this project's NMR/EPR paths.
Gamma-only real-wavefunction and ultrasoft representations are rejected.
A Gamma k point represented by complex wavefunctions is supported.

Storage uses the existing `dudk_in_memory` or disk mechanism; multiple ranks
per pool force memory storage as in the other methods. Scalar products reduce
over `intra_bgrp_comm`; QE's CG solver distributes/collects response bands.

## Tests

Build against QE 7.5, then run:

```sh
PW=/path/to/qe/bin/pw.x QECONVERSE=/path/to/bin/qe-converse.x \
    python3 tests/integration/kdotp/run.py
```

The runner compares kdotp to freshly computed covariant results using the
project's EPR/NMR physical tolerances, checks serial/MPI EPR agreement,
checks reported Sternheimer residuals, and requires successful completion. It uses existing CO+ and NaCl inputs
in temporary directories and keeps logs for diagnosis. A LiCoO2 DFT+U smoke
test compares the two methods using 30 Ry and one shifted k point; this reduced
case tests projector handling, not convergence of a material prediction. The MPI cases use
2 and 4 pools with 2 ranks per pool. `make test-integration` also runs the
existing regression suite. MPI transport settings can be supplied in the
environment; they are not hard-coded into the tests.
