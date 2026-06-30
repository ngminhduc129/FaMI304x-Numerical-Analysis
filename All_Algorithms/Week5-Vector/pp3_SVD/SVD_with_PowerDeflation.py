# =============================================================================
# SVD_with_PowerDeflation.py - SVD dùng lũy thừa & xuống thang
#
# Chức năng: Phân rã SVD bằng cách dùng phương pháp lũy thừa kết
#            hợp xuống thang để tìm trị riêng (thay vì numpy.eig).
#
# Các hàm chính:
#   read_matrix(path)            - đọc ma trận
#   compute_svd(A_df, tol)       -> (U, Sigma, Vt)
#   compute_svd_2(A_df, tol)     -> (U, Sigma, Vt)
#
# Input: Đọc từ file SVD_input_A.txt
# Cách dùng: python SVD_with_PowerDeflation.py
# =============================================================================
from fractions import Fraction
from typing import List, Tuple, Union
import numpy as np
import pandas as pd

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column

def read_matrix(path: str) -> pd.DataFrame:
    """
    Reads a whitespace-delimited matrix from a text file (allowing decimals
    or integer fractions) and returns it as a pandas DataFrame of floats.
    """
    data = []
    with open(path, 'r') as f:
        for lineno, line in enumerate(f, start=1):
            tokens = line.strip().split()
            if not tokens:
                continue
            row = []
            for tok in tokens:
                try:
                    if tok.count('/') == 1:
                        num, den = tok.split('/')
                        val = int(num) / int(den)
                    else:
                        val = float(tok)
                except Exception as e:
                    raise ValueError(
                        f"Error parsing token '{tok}' on line {lineno} in {path}: {e}"
                    )
                row.append(val)
            data.append(row)

    if not data:
        # empty file → empty DataFrame
        return pd.DataFrame()

    # ensure all rows have same length
    ncol = len(data[0])
    if any(len(r) != ncol for r in data):
        lengths = [len(r) for r in data]
        raise ValueError(
            f"Row-length mismatch in {path}: expected {ncol} columns, got {lengths}"
        )

    return pd.DataFrame(data)


def print_matrix(df: pd.DataFrame, float_format: str = '{:,.9f}') -> None:
    """
    Nicely prints the DataFrame without row/col labels, formatting floats.
    """
    print(
        df.to_string(
            index=False,
            header=False,
            float_format=float_format.format
        )
    )
# input
A_df = read_matrix("SVD_input_A.txt")
print_matrix(A_df)
def deflate_once(A: np.ndarray, lam: float, v: np.ndarray, x0_left: np.ndarray, tol: float, max_iter: int, norm_ord: Union[int, float], precision: int) -> np.ndarray:
    """
    Deflate matrix A by removing the effect of eigenpair (lam, v) for non-symmetric A:
    1. Compute left eigenvector w of A: solve A^T w = lam w via power_method
    2. Normalize so that w^T v = 1
    3. A_new = A - lam * np.outer(v, w)
    """
    # Step 1: compute left eigenvector of A
    # Use power_method on A.T to get dominant eigenvector w
    _, w = power_method(
        A.T, x0_left, tol=tol, max_iter=max_iter,
        norm_ord=norm_ord, precision=precision
    )
    # Step 2: normalize w so that w^T v = 1
    scale = w.dot(v)
    if abs(scale) < 1e-12:
        raise ValueError("Left and right eigenvectors nearly orthogonal; cannot deflate.")
    w = w / scale
    # Step 3: rank-1 deflation
    print(A - lam * np.outer(v, w))
    return A - lam * np.outer(v, w)
def power_method(A,
                 x0,
                 tol=1e-6,
                 max_iter=1000,
                 norm_ord=2,
                 precision=12,
                 display=True):
    """
    Power Method for finding the dominant eigenvalue and eigenvector of A.
    Displays each iteration in a pandas DataFrame using the specified norm.

    Parameters:
      A (ndarray): Square matrix (n×n).
      x0 (ndarray): Initial guess vector (n,).
      tol (float): Convergence tolerance for eigenvalue change.
      max_iter (int): Maximum number of iterations.
      norm_ord (int, float, or str): Norm order for both eigenvalue estimation and vector normalization (e.g. 1, 2, np.inf).
      precision (int): Number of decimal places to display.

    Returns:
      eigenvalue (float): Approximate dominant eigenvalue.
      eigenvector (ndarray): Approximate eigenvector corresponding to the dominant eigenvalue.
    """
    # Normalize initial guess
    x = x0 / np.linalg.norm(x0, ord=norm_ord)

    # Step 0: record the very first y,x,λ
    y0 = x0           # or if you want y0=x0 then skip this dot
    lambda0 = np.linalg.norm(y0, ord=norm_ord)
    x0_norm = y0 / lambda0

    # history entries: [iteration, y1…yn, lambda, x1…xn]
    history = [
      [0] + y0.tolist() + [lambda0] + x0_norm.tolist()
    ]

    eigenvalue = lambda0
    n = A.shape[0]
    

    for k in range(1, max_iter+1):
        # Compute A x_{k-1}
        y = A.dot(x)
        # Estimate eigenvalue via the chosen norm
        lambda_new = np.linalg.norm(y, ord=norm_ord)
        if lambda_new == 0:
            raise ZeroDivisionError("Norm of A*x is zero; unable to continue.")
        # Normalize to get next eigenvector estimate
        x_new = y / lambda_new

        # Record iteration data: [iter, y1..yn, lambda, x1..xn]
        row = [k] + y.tolist() + [lambda_new] + x_new.tolist()
        history.append(row)

        # Check convergence on eigenvalue
        if abs(lambda_new - eigenvalue) < tol:
            break

        x = x_new
        eigenvalue = lambda_new

    # Build and display DataFrame
    cols = ['Iteration'] + [f'y{i+1}' for i in range(n)] + ['lambda'] + [f'x{i+1}' for i in range(n)]
    df = pd.DataFrame(history, columns=cols)
    if (display==True):
      print(df.to_string(
        index=False,
        float_format=lambda v: f"{v:.{precision}f}"
    ))

    return eigenvalue, x_new

def power_method_case2(
    A,
    x0,
    tol=1e-6,
    max_iter=1000,
    norm_ord=np.inf,
    precision=12,
    display=True
):
    """
    Power Method for when the two largest eigenvalues satisfy |λ1| = |λ2| > |λ3| and λ1 = -λ2.
    Uses even powers of A to separate eigencomponents, printing y and x in one table.

    Parameters:
      A (ndarray): Square matrix (n×n).
      x0 (ndarray): Initial guess vector (n,).
      tol (float): Convergence tolerance for λ^2 change.
      max_iter (int): Maximum number of iterations for k.
      norm_ord (int, float, or str): Norm order for vector normalization.
      precision (int): Decimal places for display.

    Returns:
      lambda1 (float): Approximate dominant eigenvalue (positive).
      x1 (ndarray): Eigenvector corresponding to λ1.
      x2 (ndarray): Eigenvector corresponding to λ2 = -λ1.
    """
    n = A.shape[0]

    # Step 0: input vector and its norm^2
    y0 = x0.astype(float)
    lambda0 = np.linalg.norm(y0, ord=norm_ord)
    lambda2_old = lambda0**2
    x0_norm = y0 / lambda0

    history = []
    # record iteration 0
    history.append(
      [0] + y0.tolist() + [lambda2_old] + x0_norm.tolist() + x0_norm.tolist()
    )

    # now do the usual even‐odd power iterations
    y_even = x0.astype(float)


    for k in range(1, max_iter + 1):
        if k > 1:
            y_even = A2.dot(y_even)
        # Compute next even power
        y_even2 = A2.dot(y_even)
        # Index of largest magnitude in y_even
        j = np.argmax(np.abs(y_even))
        # Estimate λ^2
        lambda2_new = y_even2[j] / y_even[j]
        # Compute λ
        lambda1 = np.sqrt(lambda2_new)
        # Compute y_{2k+1}
        y_odd = A.dot(y_even)
        # Form eigenvectors
        x1 = y_even + lambda1 * y_odd
        x2 = y_even - lambda1 * y_odd
        # Normalize
        x1 = x1 / np.linalg.norm(x1, ord=norm_ord)
        x2 = x2 / np.linalg.norm(x2, ord=norm_ord)

        # Record iteration: [k, y_even1..n, lambda2, x1_1..n, x2_1..n]
        row = [k] + y_even.tolist() + [lambda2_new] + x1.tolist() + x2.tolist()
        history.append(row)

        # Convergence check on λ^2
        if k > 1 and abs(lambda2_new - eigen2) < tol:
            break
        eigen2 = lambda2_new

    # Display results in pandas
    cols = ['Iteration'] + [f'y{i+1}' for i in range(n)] + ['lambda^2'] + [f'x1_{i+1}' for i in range(n)] + [f'x2_{i+1}' for i in range(n)]
    df = pd.DataFrame(history, columns=cols)
    if (display==True):
      print(df.to_string(
        index=False,
        float_format=lambda v: f"{v:.{precision}f}"
    ))

    # Final eigenvalue
    lambda1_final = np.sqrt(eigen2)
    return lambda1_final, x1, x2

def compute_all_eigenpairs(
    A: np.ndarray,
    x0: np.ndarray,
    tol: float = 1e-6,
    max_iter: int = 200,
    norm_ord: Union[int, float] = 2,
    precision: int = 8
) -> Tuple[List[float], List[np.ndarray]]:
    """
    Compute all real eigenvalues and eigenvectors of A via deflation.

    Returns:
      eigs: list of eigenvalues (sorted by absolute value descending)
      vecs: list of corresponding eigenvectors
    """
    n = A.shape[0]
    M = A.copy().astype(float)
    y = x0.astype(float).flatten()
    eigs = []
    vecs = []

    for k in range(n):
        # Determine which power-method to apply: always use case1 for real eigenvalues
        print(f"Finding eigenpair {k+1}: matrix size {M.shape}")
        # Call power_method (case1) for new matrix M
        lam, v = power_method(
            M, y, tol=tol, max_iter=max_iter,
            norm_ord=norm_ord, precision=precision, display=False
        )
        eigs.append(lam)
        vecs.append(v)
        # Deflate: pass y and solver parameters to compute left eigenvector
        M = deflate_once(M, lam, v, y, tol, max_iter, norm_ord, precision)
        # Reset y to a standard initial guess in the deflated space
        y = x0.astype(float).flatten()

    # Sort eigenpairs by descending magnitude of eigenvalue
    pairs = sorted(zip(eigs, vecs), key=lambda ev: -abs(ev[0]))
    eigs_sorted, vecs_sorted = zip(*pairs)

 # Display final results: eigenvalue and eigenvector on same line
    print("All eigenvalues and eigenvectors (sorted):")
    for i, (lam, vec) in enumerate(zip(eigs_sorted, vecs_sorted), start=1):
        vec_str = ", ".join(f"{x:.{precision}f}" for x in vec)
        print(f"  λ{i} = {lam:.{precision}f}, v = [{vec_str}]")

    return list(eigs_sorted), list(vecs_sorted)
def complete_orthonormal_basis(B: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    """
    Given B (n×k) with orthonormal columns, extend to a full n×n orthonormal basis.
    """
    n, k = B.shape
    Q = np.zeros((n, n))
    Q[:, :k] = B
    for j in range(k, n):
        # start with a random vector
        v = np.random.rand(n)
        # orthonormalize against existing columns
        for i in range(j):
            v -= np.dot(Q[:, i], v) * Q[:, i]
        norm = np.linalg.norm(v)
        # if vector is nearly zero, retry
        while norm < tol:
            v = np.random.rand(n)
            for i in range(j):
                v -= np.dot(Q[:, i], v) * Q[:, i]
            norm = np.linalg.norm(v)
        Q[:, j] = v / norm
    return Q
def compute_svd(A_df: pd.DataFrame, tol: float = 1e-10) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute the SVD of A (m×n) via the smaller of A·Aᵀ or Aᵀ·A.
    Prints intermediate steps using pandas formatting.

    Returns:
        U_df: m×m DataFrame of left singular vectors
        Sigma_df: m×n DataFrame of singular values on the diagonal
        Vt_df: n×n DataFrame of right singular vectors transposed
    """
    # Convert to NumPy
    A = A_df.values.astype(float)
    m, n = A.shape

    print("--- Input matrix A ---")
    print_matrix(A_df)

    # Choose smaller route
    if True:
        print("Route: m <= n, use M = A · Aᵀ ({}×{})".format(m, m))
        M = A.dot(A.T)
        M_df = pd.DataFrame(M)
        print("Matrix M:")
        print_matrix(M_df)

        # Eigen-decomposition of M
        x0_eig = np.ones(m)   # or any nonzero initial guess of length m
        eigvals_list, eigvecs_list = compute_all_eigenpairs(
            M,
            x0_eig,
            tol=tol,
            max_iter=500,
            norm_ord=2,
            precision=7      # set to 0 to suppress per-iteration printing here
        )
        eigvals = np.array(eigvals_list)
        eigvecs = np.column_stack(eigvecs_list)
        
        # Print eigenpairs
        for i, lam in enumerate(eigvals, start=1):
            print(f"Eigenvalue λ{i} = {lam}")
            print(f"Eigenvector u{i}:")
            print_matrix(pd.DataFrame(eigvecs[:, i-1]))

        # Singular values
        sigmas = np.sqrt(np.clip(eigvals, 0, None))
        print("Singular values:")
        for i, s in enumerate(sigmas, start=1):
            print(f"σ{i} = {s}")

        # Build Σ
        Sigma = np.zeros((m, n))
        for i in range(min(m, n)):
            Sigma[i, i] = sigmas[i]
        Sigma_df = pd.DataFrame(Sigma)
        print("Sigma matrix Σ:")
        print_matrix(Sigma_df)

        # Compute right singular vectors V
        V = np.zeros((n, n))
        for i in range(min(m, n)):
            if sigmas[i] > tol:
                v = A.T.dot(eigvecs[:, i]) / sigmas[i]
            else:
                v = np.zeros(n)
            V[:, i] = v
            print(f"Right singular vector v{i+1}:")
            print_matrix(pd.DataFrame(v))

        # Complete V to full orthonormal basis
        V_complete = complete_orthonormal_basis(V[:, :m], tol)
        Vt_df = pd.DataFrame(V_complete.T)

        # U is full eigenvector matrix of M
        U_complete = eigvecs
        U_df = pd.DataFrame(U_complete)

    return U_df, Sigma_df, Vt_df

def compute_svd_2(A_df: pd.DataFrame, tol: float = 1e-10) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute the SVD of A (m×n) via the smaller of Aᵀ·A.
    Prints intermediate steps using pandas formatting.

    Returns:
        U_df: m×m DataFrame of left singular vectors
        Sigma_df: m×n DataFrame of singular values on the diagonal
        Vt_df: n×n DataFrame of right singular vectors transposed
    """
    # Convert to NumPy
    A = A_df.values.astype(float)
    m, n = A.shape

    print("--- Input matrix A ---")
    print_matrix(A_df)

    # Choose smaller route
    if True:
        print("Route: m > n, use N = Aᵀ · A ({}×{})".format(n, n))
        N = A.T.dot(A)
        N_df = pd.DataFrame(N)
        print("Matrix N:")
        print_matrix(N_df)

        # Eigen-decomposition of M
        x0_eig = np.ones(n)   # or any nonzero initial guess of length m
        eigvals_list, eigvecs_list = compute_all_eigenpairs(
            N,
            x0_eig,
            tol=tol,
            max_iter=500,
            norm_ord=2,
            precision=7     # set to 0 to suppress per-iteration printing here
        )
        eigvals = np.array(eigvals_list)
        eigvecs = np.column_stack(eigvecs_list)

        for j, lam in enumerate(eigvals, start=1):
            print(f"Eigenvalue λ{j} = {lam}")
            print(f"Eigenvector v{j}:")
            print_matrix(pd.DataFrame(eigvecs[:, j-1]))

        sigmas = np.sqrt(np.clip(eigvals, 0, None))
        print("Singular values:")
        for j, s in enumerate(sigmas, start=1):
            print(f"σ{j} = {s}")

        Sigma = np.zeros((m, n))
        for j in range(min(m, n)):
            Sigma[j, j] = sigmas[j]
        Sigma_df = pd.DataFrame(Sigma)
        print("Sigma matrix Σ:")
        print_matrix(Sigma_df)

        # Compute left singular vectors U
        U = np.zeros((m, m))
        for j in range(min(m, n)):
            if sigmas[j] > tol:
                u = A.dot(eigvecs[:, j]) / sigmas[j]
            else:
                u = np.zeros(m)
            U[:, j] = u
            print(f"Left singular vector u{j+1}:")
            print_matrix(pd.DataFrame(u))

        # Complete U to full orthonormal basis
        U_complete = complete_orthonormal_basis(U[:, :n], tol)
        U_df = pd.DataFrame(U_complete)

        # V is full eigenvector matrix of N
        V_complete = eigvecs
        Vt_df = pd.DataFrame(V_complete.T)

    return U_df, Sigma_df, Vt_df

U, S, Vt = compute_svd(A_df, tol=1e-9)
print("--- Final U ---")
print_matrix(U)
print("--- Final Σ ---")
print_matrix(S)
print("--- Final V^T ---")
print_matrix(Vt)
U, S, Vt = compute_svd_2(A_df, tol=1e-9)
print("--- Final U ---")
print_matrix(U)
print("--- Final Σ ---")
print_matrix(S)
print("--- Final V^T ---")
print_matrix(Vt)




