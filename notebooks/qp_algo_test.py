import numpy as np
from qp_algo import qp_active_set, qp_eq, qp_eq_qr
from scipy import linalg as SLA

# equality constrained QP
G_mat = np.array([[2, 5, 0], [5, 2, 0], [0, 0, 4]])
A_mat = np.array([[1, 1], [1, -2], [1, -3]])
d_vec = np.array([0, -3, -7])
b_vec = np.array([1, -2])


# inequality constrained QP
G_mat_asm = 2 * np.eye(2, dtype=float)
A_mat_asm = np.array([[-1, 2], [1, 2], [1, -2], [-1, 0], [0, -1]], dtype=float).T
b_vec_asm = np.array([2, 6, 2, 0, 0], dtype=float)
d_vec_asm = np.array([-2, -5], dtype=float)


def test_qp_eq():
    sol, lam = qp_eq(G_mat, d_vec, A_mat, b_vec)
    assert np.allclose(sol, np.array([0.4, -0.6, 1.2]))


def test_qp_eq_qr():
    Q_mat, R_mat = SLA.qr(A_mat)
    sol, lam = qp_eq_qr(G_mat, d_vec, A_mat, Q_mat, R_mat, b_vec)
    assert np.allclose(sol, np.array([0.4, -0.6, 1.2]))


def test_qp_active_set():
    sol, trace_dict = qp_active_set(
        G_mat_asm,
        d_vec_asm,
        A_mat_asm,
        b_vec_asm,
        np.array([2, 0], dtype=float),
        0,
        verbose=True,
    )
    assert np.allclose(sol, np.array([1.4, 1.7]))
