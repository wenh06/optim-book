"""
"""

import numpy as np
from scipy import linalg as SLA


def qp_eq(G: np.ndarray, d: np.ndarray, A: np.ndarray, b: np.ndarray) -> np.ndarray:
    r"""
    Solve the equality-constrained quadratic program using QR decomposition.

        .. math::
            \min 1/2 x^T G x + d^T x,
            s.t. a_i^T x = b, i = 1, ..., m

    Parameters
    ----------
    G: ndarray, shape (n, n)
        Hessian matrix of the objective function.
    d: ndarray, shape (n,)
        Linear term of the objective function.
    A: ndarray, shape (n, m)
        Transpose of the constraint matrix.
    b: ndarray, shape (m,)
        Constraint vector.

    Returns
    -------
    x: ndarray, shape (n,)
        Optimal solution.
    lambda_: ndarray, shape (m,)
        Lagrange multipliers.

    """
    n = G.shape[0]
    m = A.shape[1]
    assert (G.T == G).all(), "G must be symmetric."
    assert A.shape[0] == n, "A must have shape (n, m)."
    assert b.shape == (m,), "b must have shape (m,)."
    # Compute the QR decomposition of A.
    # Q: ndarray, shape (n, n), orthonormal
    # R: ndarray, shape (n, m), upper triangular
    Q, R = SLA.qr(A, mode="full")
    Z = Q[:, m:]
    Q = Q[:, :m]
    R = R[:m, :]
    # Compute the solution.
    Yb = Q @ np.linalg.solve(R.T, b)
    x = Yb + Z @ np.linalg.solve(Z.T @ G @ Z, -Z.T @ (d + G @ Yb))
    lambda_ = np.linalg.solve(A[:m], -(G @ x + d)[:m])
    return x, lambda_


def qp_eq_qr(
    G: np.ndarray,
    d: np.ndarray,
    A: np.ndarray,
    Q: np.ndarray,
    R: np.ndarray,
    b: np.ndarray,
) -> np.ndarray:
    r"""
    Solve the equality-constrained quadratic program using QR decomposition,
    and the QR decomposition of the transpose of the constraint matrix is given.

        .. math::
            \min 1/2 x^T G x + d^T x,
            s.t. a_i^T x = b, i = 1, ..., m

    Parameters
    ----------
    G: ndarray, shape (n, n)
        Hessian matrix of the objective function.
    d: ndarray, shape (n,)
        Linear term of the objective function.
    A: ndarray, shape (n, m)
        Transpose of the constraint matrix.
    Q: ndarray, shape (n, n), orthonormal
        Orthogonal matrix of the QR decomposition of A.
    R: ndarray, shape (n, m), upper triangular
        Upper triangular matrix of the QR decomposition of A.
    b: ndarray, shape (m,)
        Constraint vector.

    Returns
    -------
    x: ndarray, shape (n,)
        Optimal solution.
    lambda_: ndarray, shape (m,)
        Lagrange multipliers.

    """
    n = G.shape[0]
    m = R.shape[1]
    assert (G.T == G).all(), "G must be symmetric."
    assert Q.shape == (n, n), "Q must have shape (n, n)."
    assert R.shape == (n, m), "R must have shape (n, m)."
    assert b.shape == (m,), "b must have shape (m,)."
    # Intermediate matrices.
    Z = Q[:, m:]
    Q = Q[:, :m]
    R = R[:m, :]
    # Compute the solution.
    Yb = Q @ np.linalg.solve(R.T, b)
    x = Yb + Z @ np.linalg.solve(Z.T @ G @ Z, -Z.T @ (d + G @ Yb))
    lambda_ = np.linalg.solve(A[:m], -(G @ x + d)[:m])
    return x, lambda_


def qp_active_set(
    G: np.ndarray,
    d: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    x0: np.ndarray,
    num_eq_constraints: int,
    max_iter: int = 100,
    precision: float = 1e-10,
    verbose: bool = False,
) -> np.ndarray:
    r"""
    Solve the quadratic program using the active-set method. (NOT finished yet!!!)

        .. math::
            \min 1/2 x^T G x + d^T x,
            s.t. a_i^T x = b, i = 1, ..., m_1
                 a_i^T x \le b, i = m_1, ..., m

    Parameters
    ----------
    G: ndarray, shape (n, n)
        Hessian matrix of the objective function.
    d: ndarray, shape (n,)
        Linear term of the objective function.
    A: ndarray, shape (n, m)
        Transpose of the constraint matrix.
    b: ndarray, shape (m,)
        Constraint vector.
    x0: ndarray, shape (n,)
        Initial guess.
    num_eq_constraints: int
        Number of equality constraints.
    max_iter: int, optional
        Maximum number of iterations.
    precision: float, default 1e-10
        Precision of the solution.
    verbose: bool, default False
        If True, the intermediate results are tracked and returned as a dict.

    Returns
    -------
    sol: ndarray, shape (n,)
        Optimal solution.

    """
    n = G.shape[0]
    m = A.shape[1]
    assert (G.T == G).all(), "G must be symmetric."
    assert A.shape[0] == n, "A must have shape (n, m)."
    assert b.shape == (m,), "b must have shape (m,)."
    assert x0.shape == (n,), "x0 must have shape (n,)."
    assert num_eq_constraints <= m, "num_eq_constraints must be <= m."
    # Compute the initial constraint violation.
    c = A.T @ x0 - b
    # Check the initial guess is feasible.
    assert (
        c[:num_eq_constraints] == 0
    ).all(), "x0 must satisfy the equality constraints."
    assert (
        c[num_eq_constraints:] <= 0
    ).all(), "x0 must satisfy the inequality constraints."
    # Initialize the solution.
    x = x0.copy()
    # Initialize the active set.
    active_set = np.where(A.T @ x0 == b)[0]

    # Compute the initial gradient vector.
    g = G @ x + d
    # Compute the initial QR decomposition of A_active.
    Q, R = SLA.qr(A[:, active_set], mode="full")

    trace_dict = {
        "x": [],
        "active_set": [],
        "direction": [],
        "lambda_": [],
        "q": [],
        "alpha": [],
        "j": [],
    }

    # Initialize the number of iterations.
    num_iter = 0
    while True:
        trace_dict["x"].append(x.copy())
        trace_dict["active_set"].append(active_set.copy())
        # Solve the reduced equality-constrained quadratic program,
        # get the search direction, and the Lagrange multipliers.
        direction, lambda_ = qp_eq_qr(G, g, Q, R, b[active_set])
        trace_dict["direction"].append(direction.copy())
        trace_dict["lambda_"].append(lambda_.copy())
        # if direction is zero
        if np.linalg.norm(direction) < precision:
            ind = np.argmin(lambda_)
            q = active_set[ind]
            trace_dict["q"].append(q)
            trace_dict["j"].append(np.nan)
            trace_dict["alpha"].append(np.nan)
            if lambda_[ind] < 0:
                active_set = np.delete(active_set, ind)
                Q, R = SLA.qr_delete(Q, R, ind, which="col")
            else:
                break
        else:
            # Compute the step size.
            candidates = np.array(
                [
                    (b[i] - A[:, i].T @ x) / (A[:, i].T @ direction)
                    for i in range(m)
                    if i not in active_set and A[:, i].T @ direction > 0
                ]
            )
            if candidates.size > 0:
                j = np.argmin(candidates)
                alpha = (b[j] - A[:, j].T @ x) / (A[:, j].T @ direction)
            else:
                j = np.nan
                alpha = 1
            trace_dict["q"].append(np.nan)
            trace_dict["j"].append(j)
            trace_dict["alpha"].append(alpha)
            # Update the solution.

            # TODO: finish the rest of the algorithm.

        num_iter += 1

        if num_iter >= max_iter:
            break

    if verbose:
        return x, trace_dict

    return x
