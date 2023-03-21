"""
"""

import numpy as np
from scipy import linalg as SLA


def qp_eq(
    G: np.ndarray,
    d: np.ndarray,
    A: np.ndarray,
    b: np.ndarray,
    return_lambda: bool = True,
) -> np.ndarray:
    """Solve the equality-constrained quadratic program
    using QR decomposition.

    .. math::

        \\begin{equation}
        \\begin{array}{cl}
        \\min & 1/2 x^T G x + d^T x,
        s.t.  & a_i^T x = b, i = 1, ..., m
        \\end{array}
        \\end{equation}

    Parameters
    ----------
    G : numpy.ndarray
        Hessian matrix of the objective function,
        of shape ``(n, n)``.
    d : numpy.ndarray
        Linear term of the objective function,
        of shape ``(n,)``.
    A : numpy.ndarray
        Transpose of the constraint matrix,
        of shape ``(n, m)``.
    b : numpy.ndarray
        Constraint vector, of shape ``(m,)``.
    return_lambda: bool, default True
        Whether to return the Lagrange multipliers.

    Returns
    -------
    x: numpy.ndarray
        Optimal solution, of shape ``(n,)``.
    lambda_: numpy.ndarray, optional
        Lagrange multipliers, of shape ``(m,)``.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        from scipy import linalg as SLA

        G_mat = np.array([[2, 5, 0], [5, 2, 0], [0, 0, 4]])
        A_mat = np.array([[1, 1], [1, -2], [1, -3]])
        d_vec = np.array([0, -3, -7])
        b_vec = np.array([1, -2])

        sol, lam = qp_eq(G_mat, d_vec, A_mat, b_vec)

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
    if return_lambda:
        lambda_ = np.linalg.solve(A[:m], -(G @ x + d)[:m])
        return x, lambda_
    return x


def qp_eq_qr(
    G: np.ndarray,
    d: np.ndarray,
    A: np.ndarray,
    Q: np.ndarray,
    R: np.ndarray,
    b: np.ndarray,
    return_lambda: bool = True,
) -> np.ndarray:
    """
    Solve the equality-constrained quadratic program using QR decomposition,
    with QR decomposition of the transpose of the constraint matrix being given.

    .. math::

        \\begin{equation}
        \\begin{array}{cl}
        \\min & 1/2 x^T G x + d^T x,
        s.t.  & a_i^T x = b, i = 1, ..., m
        \\end{array}
        \\end{equation}

    Parameters
    ----------
    G : numpy.ndarray
        Hessian matrix of the objective function,
        of shape ``(n, n)``.
    d : numpy.ndarray
        Linear term of the objective function,
        of shape ``(n,)``.
    A : numpy.ndarray
        Transpose of the constraint matrix,
        of shape ``(n, m)``.
    Q : numpy.ndarray
        Orthogonal matrix of the QR decomposition of A,
        of shape ``(n, n)``, orthonormal.
    R : numpy.ndarray
        Upper triangular matrix of the QR decomposition of A,
        of shape ``(n, m)``.
    b : numpy.ndarray
        Constraint vector, of shape ``(m,)``.
    return_lambda : bool, default True
        Whether to return the Lagrange multipliers.

    Returns
    -------
    x : numpy.ndarray
        Optimal solution, of shape ``(n,)``.
    lambda_ : numpy.ndarray, optional
        Lagrange multipliers, of shape ``(m,)``.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        from scipy import linalg as SLA

        G_mat = np.array([[2, 5, 0], [5, 2, 0], [0, 0, 4]])
        A_mat = np.array([[1, 1], [1, -2], [1, -3]])
        d_vec = np.array([0, -3, -7])
        b_vec = np.array([1, -2])

        Q_mat, R_mat = SLA.qr(A_mat)

        sol, lam = qp_eq_qr(G_mat, d_vec, A_mat, Q_mat, R_mat, b_vec)

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
    if return_lambda:
        lambda_ = np.linalg.solve(A[:m], -(G @ x + d)[:m])
        return x, lambda_
    return x


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
    """Solve the quadratic program using the active-set method.

    (NOT finished yet!!! mainly rank checking, feasibility checking, etc.)

    .. math::

        \\begin{equation}
        \\begin{array}{cl}
        \\min & 1/2 x^T G x + d^T x,
        s.t.  & a_i^T x = b, i = 1, ..., m_1
              & a_i^T x \\le b, i = m_1, ..., m
        \\end{array}
        \\end{equation}

    Parameters
    ----------
    G : numpy.ndarray
        Hessian matrix of the objective function,
        of shape ``(n, n)``.
    d : numpy.ndarray
        Linear term of the objective function,
        of shape ``(n,)``.
    A : numpy.ndarray
        Transpose of the constraint matrix,
        of shape ``(n, m)``.
    b : numpy.ndarray
        Constraint vector, of shape ``(m,)``.
    x0 : numpy.ndarray, shape (n,)
        Initial guess.
    num_eq_constraints : int
        Number of equality constraints.
    max_iter : int, optional
        Maximum number of iterations.
    precision : float, default 1e-10
        Precision of the solution.
    verbose : bool, default False
        If True, the intermediate results are tracked
        and returned as a dict.

    Returns
    -------
    sol : numpy.ndarray
        Optimal solution, of shape ``(n,)``.
    trace_dict : dict, optional,
        A dict containing the intermediate results.

    Examples
    --------
    .. code-block:: python

        import numpy as np
        from scipy import linalg as SLA

        G_mat_asm = 2 * np.eye(2, dtype=float)
        A_mat_asm = np.array([[-1, 2], [1, 2], [1, -2], [-1, 0], [0, -1]], dtype=float).T
        b_vec_asm = np.array([2, 6, 2, 0, 0], dtype=float)
        d_vec_asm = np.array([-2, -5], dtype=float)

        sol, trace_dict = qp_active_set(
            G_mat_asm,
            d_vec_asm,
            A_mat_asm,
            b_vec_asm,
            np.array([2, 0], dtype=float),
            0,
            verbose=True,
        )

    """
    n = G.shape[0]
    m = A.shape[1]
    assert (G.T == G).all(), "G must be symmetric."
    assert A.shape[0] == n, "A must have shape (n, m)."
    assert b.shape == (m,), "b must have shape (m,)."
    assert x0.shape == (n,), "x0 must have shape (n,)."
    assert num_eq_constraints <= m, "num_eq_constraints must be <= m."

    trace_dict = {
        "x": [],
        "active_set": [],
        "direction": [],
        "lambda_": [],
        "q": [],
        "alpha": [],
        "j": [],
    }

    if num_eq_constraints == m:
        # If all constraints are equality constraints, use `qp_eq`.
        # NOTE that in this case, we DO NOT check if x0 is feasible or not.
        x, lambda_ = qp_eq(G, d, A, b)
        trace_dict["x"].append(x)
        trace_dict["lambda_"].append(lambda_)
        if verbose:
            return x, trace_dict
        return x

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
    # Compute the initial QR decomposition of A_active.
    Q, R = SLA.qr(A[:, active_set], mode="full")

    # Initialize the number of iterations.
    num_iter = 0
    while True:
        trace_dict["x"].append(x.copy())
        trace_dict["active_set"].append(active_set.copy())

        # Compute the initial gradient vector.
        g = G @ x + d

        # Solve the reduced equality-constrained quadratic program,
        # get the search direction, and the Lagrange multipliers.
        direction = qp_eq_qr(
            G, g, A[:, active_set], Q, R, np.zeros(len(active_set)), return_lambda=False
        )
        trace_dict["direction"].append(direction.copy())

        # if direction is zero
        if np.linalg.norm(direction) < precision:
            lambda_ = np.linalg.lstsq(A[:, active_set], -g, rcond=None)[0]
            trace_dict["lambda_"].append(lambda_.copy())
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
            trace_dict["lambda_"].append(np.nan)
            # Compute the step size.
            candidates = np.array(
                [
                    (b[i] - A[:, i].T @ x) / (A[:, i].T @ direction)
                    for i in range(num_eq_constraints, m)
                    if i not in active_set and A[:, i].T @ direction > 0
                ]
            )
            if candidates.size > 0:
                j = np.argmin(candidates)
                alpha = min(1, candidates[j])
            else:
                j = np.nan
                alpha = 1
            trace_dict["q"].append(np.nan)
            trace_dict["j"].append(j)
            trace_dict["alpha"].append(alpha)
            # Update the solution.
            x += alpha * direction
            # Update the active set.
            if alpha < 1:
                active_set = np.append(active_set, j)
                Q, R = SLA.qr_insert(Q, R, A[:, j], j, which="col")
        num_iter += 1
        if num_iter >= max_iter:
            break

    if verbose:
        return x, trace_dict

    return x
