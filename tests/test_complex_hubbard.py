"""Algebraic regressions for complex U and U+V (requires NumPy).

These complement, and do not replace, the binary/MPI integration test.
"""
import unittest

import numpy as np


def energy(n, u):
    return 0.5 * u * np.trace(n - n @ n).real


def potential(n, u):
    return u * (0.5 * np.eye(len(n)) - n)


class ComplexHubbardTests(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(473)

    def hermitian(self, size):
        a = self.rng.normal(size=(size, size)) + 1j * self.rng.normal(size=(size, size))
        return (a + a.conj().T) / 2

    def test_pure_complex_orbital_has_zero_u_energy(self):
        p = np.array([1, 1j, 0, 0, 0]) / np.sqrt(2)
        n = np.outer(p, p.conj())
        self.assertAlmostEqual(energy(n, 4), 0)
        self.assertAlmostEqual(energy(n.real, 4), 1)

    def test_u_potential_is_energy_derivative(self):
        n, dn = self.hermitian(5), self.hermitian(5)
        eps = 1e-6
        fd = (energy(n + eps * dn, 4) - energy(n - eps * dn, 4)) / (2 * eps)
        self.assertAlmostEqual(fd, np.trace(potential(n, 4) @ dn).real, places=7)

    def test_u_limit_of_directed_uv_operator(self):
        # QE stores nsg(m2,m1)=n(m1,m2), v_nsg=-U*conjg(nsg)+U/2.
        # Its two half-weighted Hermitian branches reconstruct the U operator.
        n = self.hermitian(5)
        v_nsg = -4 * n.T.conj() + 2 * np.eye(5)
        applied = 0.5 * (v_nsg.conj().T + v_nsg)
        np.testing.assert_allclose(applied, potential(n, 4), atol=1e-13)
        e_uv = 2 * (np.trace(n.T).real - np.sum(abs(n.T)**2))
        self.assertAlmostEqual(e_uv, energy(n, 4))

    def test_unitary_orbital_covariance(self):
        n = self.hermitian(5)
        q, _ = np.linalg.qr(self.hermitian(5))
        nr = q.conj().T @ n @ q
        self.assertAlmostEqual(energy(nr, 4), energy(n, 4))
        np.testing.assert_allclose(potential(nr, 4), q.conj().T @ potential(n, 4) @ q, atol=1e-13)

    def test_intersite_energy_derivative_with_phase(self):
        # A directed I-J pair and its conjugate J-I partner are counted once.
        p = self.rng.normal(size=3) + 1j * self.rng.normal(size=3)
        q = self.rng.normal(size=5) + 1j * self.rng.normal(size=5)
        n = np.outer(p, q.conj()) * np.exp(-0.37j)
        dn = self.rng.normal(size=n.shape) + 1j * self.rng.normal(size=n.shape)
        v = 2.0
        def ev(x):
            return -v * np.sum(abs(x)**2)
        eps = 1e-6
        fd = (ev(n + eps * dn) - ev(n - eps * dn)) / (2 * eps)
        analytic = -2 * v * np.vdot(n, dn).real
        self.assertAlmostEqual(fd, analytic, places=6)
        self.assertGreater(abs(ev(n) - ev(n.real)), 1e-3)

    def test_mixing_metric_retains_imaginary_residual(self):
        residual = np.array([[0, 0.2j], [-0.2j, 0]])
        metric = np.vdot(residual, residual).real
        self.assertAlmostEqual(metric, 0.08)
        mixed = 0.3 * residual
        self.assertGreater(np.max(abs(mixed.imag)), 0)
        np.testing.assert_allclose(mixed, mixed.conj().T)

    def test_partitioned_k_sum(self):
        p = self.rng.normal(size=(8, 5, 3)) + 1j * self.rng.normal(size=(8, 5, 3))
        nk = np.einsum('kmb,knb->kmn', p, p.conj()) / 8
        for pools in (1, 2, 4):
            reduced = sum(np.sum(nk[i::pools], axis=0) for i in range(pools))
            np.testing.assert_allclose(reduced, np.sum(nk, axis=0), atol=1e-13)


if __name__ == '__main__':
    unittest.main()
