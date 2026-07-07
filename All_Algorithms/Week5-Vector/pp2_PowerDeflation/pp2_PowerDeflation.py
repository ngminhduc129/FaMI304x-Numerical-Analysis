# =============================================================================
# pp2_PowerDeflation.py - Phương pháp lũy thừa & xuống thang
#
# Chức năng: Tìm tất cả giá trị riêng và vector riêng của ma trận
#            bằng kết hợp phương pháp lũy thừa (power method) và
#            xuống thang (deflation). Hỗ trợ trị riêng thực và phức.
#
# Các hàm chính:
#   power_method(A, x0, ...)                       - lũy thừa thường
#   power_method_case2(A, x0, ...)                 - 2 trị trội bằng nhau
#   power_method_case3(A, x0, ...)                 - trị phức liên hợp
#   compute_all_eigenpairs(A, x0)                  - tổng hợp tất cả
#   deflate_once(A, lam, v, x0_left, ...)          - xuống thang
#
# Input: Đọc từ PWDF_input_A.txt, PWDF_input_A3.txt
# Cách dùng: python pp2_PowerDeflation.py
# =============================================================================
from pathlib import Path
import contextlib
from fractions import Fraction
from typing import List, Tuple, Union
import numpy as np
import pandas as pd

__dir__ = Path(__file__).parent.resolve()

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def input_matrix(filename, convert_fractions=False):
    """
    Reads a matrix from a text file and returns it as a NumPy array.
    Supports fractional entries if present (though your files use decimals).
    """
    matrix = []

    with open(str(__dir__ / filename), 'r') as f:
        for line in f:
            tokens = line.strip().split()
            if not tokens:
                continue

            row = []
            for token in tokens:
                if '/' in token:
                    # convert fractions if any includes fraction sign
                    val = Fraction(token)
                    row.append(float(val) if convert_fractions else val)
                else:
                    # parse as float directly
                    row.append(float(token))

            matrix.append(row)
            
    dtype = float if convert_fractions else object
    return np.array(matrix, dtype=dtype)
def output_matrix(X: np.ndarray, precision: int = 7):
    """
    Prints a NumPy array (vector or matrix) in a clean tabular format using pandas.
    
    Parameters:
    - X: np.ndarray, 1D or 2D array.
    - precision: number of decimal places to round floats to.
    """
    # Wrap 1D arrays into a 2D DataFrame for consistent display
    if X.ndim == 1:
        df = pd.DataFrame(X, columns=["value"])
    elif X.ndim == 2:
        df = pd.DataFrame(X)
    else:
        raise ValueError("Only 1D or 2D arrays are supported.")
    
    # Round floats
    df = df.round(precision)
    # Print without index/header for cleaner look
    print(df.to_string(index=False, header=False))

def verify_eigenpair(A: np.ndarray, lam, v: np.ndarray) -> float:
    residual = np.linalg.norm(A @ v - lam * v) / np.linalg.norm(A)
    return float(abs(residual))

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
        norm_ord=norm_ord, precision=precision,
        display=True
    )
    # Step 2: normalize w so that w^T v = 1
    scale = w.dot(v)
    if abs(scale) < 1e-12:
        raise ValueError("Vector riêng trái và phải gần như trực giao; không thể xuống thang.")
    w = w / scale
    # Step 3: rank-1 deflation
    return A - lam * np.outer(v, w)
def power_method(A,
                 x0,
                 tol=1e-6,
                 max_iter=1000,
                 norm_ord=2,
                  precision=7,
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
            raise ZeroDivisionError("Chuẩn của A*x bằng 0; không thể tiếp tục.")
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
    cols = ['Lần lặp'] + [f'y{i+1}' for i in range(n)] + ['lambda'] + [f'x{i+1}' for i in range(n)]
    df = pd.DataFrame(history, columns=cols)
    if display:
        if len(df) > 6:
            df_show = pd.concat([df.head(3), df.tail(3)])
            print(df_show.to_string(index=False,
                                    float_format=lambda x: f"{x:.{precision}f}"))
            print("    ... (bỏ qua", len(df) - 6, "lần lặp) ...")
        else:
            print(df.to_string(index=False,
                               float_format=lambda x: f"{x:.{precision}f}"))

    return eigenvalue, x_new

def power_method_case2(
    A,
    x0,
    tol=1e-6,
    max_iter=1000,
    norm_ord=np.inf,
    precision=7,
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
    A2 = A @ A

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
    cols = ['Lần lặp'] + [f'y{i+1}' for i in range(n)] + ['lambda^2'] + [f'x1_{i+1}' for i in range(n)] + [f'x2_{i+1}' for i in range(n)]
    df = pd.DataFrame(history, columns=cols)
    if display:
        if len(df) > 6:
            df_show = pd.concat([df.head(3), df.tail(3)])
            print(df_show.to_string(index=False,
                                    float_format=lambda x: f"{x:.{precision}f}"))
            print("    ... (bỏ qua", len(df) - 6, "lần lặp) ...")
        else:
            print(df.to_string(index=False,
                               float_format=lambda x: f"{x:.{precision}f}"))

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
        print(f"Đang tìm cặp trị riêng {k+1}: kích thước ma trận {M.shape}")
        # Call power_method (case1) for new matrix M
        lam, v = power_method(
            M, y, tol=tol, max_iter=max_iter,
            norm_ord=norm_ord, precision=precision, display=True
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
    print("Tất cả trị riêng và vector riêng (sắp xếp):")
    for i, (lam, vec) in enumerate(zip(eigs_sorted, vecs_sorted), start=1):
        vec_str = ", ".join(f"{x:.{precision}f}" for x in vec)
        print(f"  λ{i} = {lam:.{precision}f}, v = [{vec_str}]")

    return list(eigs_sorted), list(vecs_sorted)
def power_method_case3(
    A: np.ndarray,
    x0: np.ndarray,
    tol: float = 1e-6,
    max_iter: int = 100,
    norm_ord: Union[int, float, str] = np.inf,
    precision: int = 7,
    display: bool = True
) -> (complex, complex, np.ndarray, np.ndarray):
    """
    Power‐Method variant for the case where the two largest eigenvalues of A
    are complex conjugates (|λ1| = |λ2| > |λ3| ≥ ...),
    i.e. λ2 = conj(λ1).  We estimate p and q in the quadratic Z^2 + pZ + q = 0
    and then recover λ1, λ2 and their eigenvectors.

    Parameters:
    - A (ndarray): square matrix (n×n), real or complex.
    - x0 (ndarray): initial guess vector (n,), complex or real.
    - tol (float): tolerance on p and q convergence.
    - max_iter (int): maximum number of successive powers k to attempt.
    - norm_ord: norm order for eigenvector normalization (1, 2, np.inf, etc.).
    - precision (int): number of decimal places to print in pandas output.
    - display (bool): whether to print iteration table.

    Returns:
    - lambda1, lambda2 (complex): the two conjugate eigenvalues of largest magnitude.
    - v1, v2 (ndarray): normalized eigenvectors corresponding to λ1 and λ2 (n,).
    """

    n = A.shape[0]
    # initial y vectors
    y_m = x0.astype(complex).reshape(n)
    y_m1 = A.dot(y_m)
    y_m2 = A.dot(y_m1)

    history = []
    # Step 0: solve for initial p0,q0 from rows 0,1
    M0 = np.array([[y_m1[0], y_m[0]], [y_m1[1], y_m[1]]], dtype=complex)
    b0 = -np.array([y_m2[0], y_m2[1]], dtype=complex)
    p0, q0 = np.linalg.solve(M0, b0)
    history.append([0] + y_m.tolist() + y_m1.tolist() + y_m2.tolist() + [p0, q0])

    p_old, q_old = p0, q0

    for m in range(max_iter):
        # Pick two distinct indices (e.g., 0 and 1) to form 2×2 system:
        i1, i2 = 0, 1
        # Solve for p, q in: y2 + p * y1 + q * y = 0  at rows i1 and i2
        M = np.array([[y_m1[i1], y_m[i1]],
                      [y_m1[i2], y_m[i2]]], dtype=complex)
        b = -np.array([y_m2[i1], y_m2[i2]], dtype=complex)

        # Solve 2×2 linear system
        try:
            p, q = np.linalg.solve(M, b)
        except np.linalg.LinAlgError:
            raise RuntimeError(
                f"Hệ 2×2 suy biến tại lần lặp m={m}; "
                "hãy chọn cặp chỉ số khác hoặc kiểm tra ma trận A."
            )

        # Record iteration data: [m, y_m, y_m1, y_m2, p, q]
        row = [m] + y_m.tolist() + y_m1.tolist() + y_m2.tolist() + [p, q]
        history.append(row)

        # Check convergence on p and q
        if m > 0 and (abs(p - p_old) < tol and abs(q - q_old) < tol):
            break

        p_old, q_old = p, q

        # Move to next power: shift (m→m+1)
        y_m = y_m1
        y_m1 = y_m2
        y_m2 = A.dot(y_m1)

        # Normalise all three vectors by the same factor to prevent
        # exponential overflow (|λ₁|^m grows fast).  The linear recurrence
        # y₂ + p·y₁ + q·y₀ = 0 is homogeneous → p,q invariant under scaling.
        nrm = float(np.max(np.abs(y_m1)))
        if nrm > 1e8:
            y_m  /= nrm
            y_m1 /= nrm
            y_m2 /= nrm

    # Build DataFrame columns
    cols = ['m'] + [f'y_m{i+1}' for i in range(n)] \
           + [f'y_m1_{i+1}' for i in range(n)] \
           + [f'y_m2_{i+1}' for i in range(n)] \
           + ['p', 'q']
    df = pd.DataFrame(history, columns=cols)

    # Print iteration history with specified precision
    if display:
        if len(df) > 6:
            df_show = pd.concat([df.head(3), df.tail(3)])
            print(df_show.to_string(index=False,
                                    float_format=lambda x: f"{x:.{precision}f}"))
            print("    ... (bỏ qua", len(df) - 6, "lần lặp) ...")
        else:
            print(df.to_string(index=False,
                               float_format=lambda x: f"{x:.{precision}f}"))

    # Now solve the quadratic for eigenvalues
    # z^2 + p z + q = 0  =>  z = [-p ± sqrt(p^2 - 4q)]/2
    discriminant = p * p - 4.0 * q
    sqrt_disc = np.sqrt(discriminant)
    lambda1 = (-p + sqrt_disc) / 2.0
    lambda2 = (-p - sqrt_disc) / 2.0

    # Build eigenvectors from formula:
    #   y_{m+1} - λ₁·y_m  ∝ v₂  (eigenvector of λ₂)
    #   y_{m+1} - λ₂·y_m  ∝ v₁  (eigenvector of λ₁)
    v1 = y_m1 - lambda2 * y_m
    v2 = y_m1 - lambda1 * y_m

    # Normalize each eigenvector
    v1 = v1 / np.linalg.norm(v1, ord=norm_ord)
    v2 = v2 / np.linalg.norm(v2, ord=norm_ord)

    # Display eigenvectors in a small DataFrame
    df_eigs = pd.DataFrame({
        'v1': v1,
        'v2': v2
    })

    return lambda1, lambda2, df_eigs

def run_all_cases(A: np.ndarray, x0: np.ndarray,
                  tol: float = 1e-6, max_iter: int = 200,
                  norm_ord: Union[int, float] = 2, precision: int = 8):
    n = A.shape[0]
    threshold = 1e-4

    print("=" * 70)
    print("CHẠY CẢ 3 CASE TRÊN CÙNG MA TRẬN")
    print("=" * 70)

    results = {}

    # ===== CASE 1 =====
    print("\n" + "=" * 70)
    print("CASE 1 - Lũy thừa thường (trị riêng trội thực, bội đơn)")
    print("=" * 70)
    try:
        M = A.copy().astype(float)
        y = x0.copy().astype(float).flatten()
        eigs1, vecs1 = [], []
        for k in range(n):
            lam, v = power_method(M, y, tol=tol, max_iter=max_iter,
                                  norm_ord=norm_ord, precision=precision,
                                  display=True)
            eigs1.append(lam)
            vecs1.append(v)
            M = deflate_once(M, lam, v, y, tol, max_iter, norm_ord, precision)
            y = x0.copy().astype(float).flatten()

        pairs = sorted(zip(eigs1, vecs1), key=lambda ev: -abs(ev[0]))
        eigs1, vecs1 = zip(*pairs)

        print(f"✅ Thành công - {n} trị riêng:")
        for i, (lam, v) in enumerate(zip(eigs1, vecs1), 1):
            r = verify_eigenpair(A, lam, v)
            stt = "✅" if r < threshold else "⚠️"
            print(f"  λ{i} = {lam:.{precision}f}  residual = {r:.2e}  {stt}")
        results['case1'] = {'status': '✅', 'eigs': list(eigs1), 'vecs': list(vecs1)}
    except Exception as e:
        print(f"❌ Thất bại: {e}")
        results['case1'] = {'status': '❌', 'error': str(e)}

    # ===== CASE 2 =====
    print("\n" + "=" * 70)
    print("CASE 2 - Hai trị trội thực, bằng trị tuyệt đối, trái dấu (λ₂ = -λ₁)")
    print("=" * 70)
    try:
        lam2a, v2a, v2b = power_method_case2(
            A, x0.copy().flatten(), tol=tol, max_iter=max_iter,
            norm_ord=norm_ord, precision=precision, display=False)
        r2a = verify_eigenpair(A, lam2a, v2a)
        r2b = verify_eigenpair(A, -lam2a, v2b)

        if r2a < threshold and r2b < threshold:
            print(f"✅ Phù hợp: λ₁ = {lam2a:.{precision}f}, λ₂ = {-lam2a:.{precision}f}")
            print(f"   residual₁ = {r2a:.2e}, residual₂ = {r2b:.2e}")

            M = A.copy().astype(float)
            y = x0.copy().astype(float).flatten()
            M = deflate_once(M, lam2a, v2a, y, tol, max_iter, norm_ord, precision)
            M = deflate_once(M, -lam2a, v2b, y, tol, max_iter, norm_ord, precision)

            eigs2 = [lam2a, -lam2a]
            vecs2 = [v2a, v2b]
            y = x0.copy().astype(float).flatten()

            for k in range(n - 2):
                try:
                    lam, v = power_method(M, y, tol=tol, max_iter=max_iter,
                                          norm_ord=norm_ord, precision=precision,
                                          display=True)
                    eigs2.append(lam)
                    vecs2.append(v)
                    M = deflate_once(M, lam, v, y, tol, max_iter, norm_ord, precision)
                    y = x0.copy().astype(float).flatten()
                except Exception:
                    break

            pairs = sorted(zip(eigs2, vecs2), key=lambda ev: -abs(ev[0]))
            eigs2, vecs2 = zip(*pairs)

            print(f"   → {len(eigs2)}/{n} trị riêng (deflation):")
            for i, (lam, v) in enumerate(zip(eigs2, vecs2), 1):
                r = verify_eigenpair(A, lam, v)
                stt = "✅" if r < threshold else "⚠️"
                print(f"     λ{i} = {lam:.{precision}f}  residual = {r:.2e}  {stt}")
            results['case2'] = {'status': '✅', 'eigs': list(eigs2), 'vecs': list(vecs2)}
        else:
            print(f"❌ Không phù hợp (residual: {r2a:.2e}, {r2b:.2e})")
            results['case2'] = {'status': '❌', 'error': f'residual lớn: {r2a:.2e}, {r2b:.2e}'}
    except Exception as e:
        print(f"❌ Thất bại: {e}")
        results['case2'] = {'status': '❌', 'error': str(e)}

    # ===== CASE 3 =====
    print("\n" + "=" * 70)
    print("CASE 3 - Hai trị trội là phức liên hợp (λ₂ = λ̅₁)")
    print("=" * 70)
    try:
        lam3a, lam3b, df3 = power_method_case3(
            A, x0.copy().astype(float).flatten(), tol=tol,
            max_iter=max_iter // 2, norm_ord=norm_ord,
            precision=precision, display=False)
        v3a = df3['v1'].values
        v3b = df3['v2'].values
        r3a = verify_eigenpair(A, lam3a, v3a)
        r3b = verify_eigenpair(A, lam3b, v3b)

        if r3a < threshold and r3b < threshold:
            print(f"✅ Phù hợp: λ₁ = {lam3a:.{precision}f}, λ₂ = {lam3b:.{precision}f}")
            print(f"   residual₁ = {r3a:.2e}, residual₂ = {r3b:.2e}")

            M = A.copy().astype(complex)
            y = x0.copy().astype(complex).flatten()
            M = deflate_once(M, lam3a, v3a, y, tol, max_iter, norm_ord, precision)
            M = deflate_once(M, lam3b, v3b, y, tol, max_iter, norm_ord, precision)

            eigs3 = [lam3a, lam3b]
            vecs3 = [v3a, v3b]
            y = x0.copy().astype(complex).flatten()

            for k in range(n - 2):
                try:
                    lam, v = power_method(M, y, tol=tol, max_iter=max_iter,
                                          norm_ord=norm_ord, precision=precision,
                                          display=True)
                    eigs3.append(lam)
                    vecs3.append(v)
                    M = deflate_once(M, lam, v, y, tol, max_iter, norm_ord, precision)
                    y = x0.copy().astype(complex).flatten()
                except Exception:
                    break

            pairs = sorted(zip(eigs3, vecs3), key=lambda ev: -abs(ev[0]))
            eigs3, vecs3 = zip(*pairs)

            print(f"   → {len(eigs3)}/{n} trị riêng (deflation):")
            for i, (lam, v) in enumerate(zip(eigs3, vecs3), 1):
                r = verify_eigenpair(A, lam, v)
                stt = "✅" if r < threshold else "⚠️"
                if np.isreal(lam):
                    fstr = f"{lam.real:.{precision}f}"
                else:
                    fstr = f"{lam.real:.{precision}f}{lam.imag:+.{precision}f}j"
                print(f"     λ{i} = {fstr}  residual = {r:.2e}  {stt}")
            results['case3'] = {'status': '✅', 'eigs': list(eigs3), 'vecs': list(vecs3)}
        else:
            print(f"❌ Không phù hợp (residual: {r3a:.2e}, {r3b:.2e})")
            results['case3'] = {'status': '❌', 'error': f'residual lớn: {r3a:.2e}, {r3b:.2e}'}
    except Exception as e:
        print(f"❌ Thất bại: {e}")
        results['case3'] = {'status': '❌', 'error': str(e)}

    print("\n" + "=" * 70)
    print("KẾT THÚC")
    print("=" * 70)
    return results


if __name__ == "__main__":
    output_path = str(__dir__ / "pp2_PowerDeflation_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        # ==== DEMO: Ma trận thực 5x5 (PWDF_input_A.txt) ====
        A = input_matrix('PWDF_input_A.txt', convert_fractions=False)
        x0 = np.array([1., 1., 1., 1., 1.])
        print("\nMa trận A (5x5 - trị riêng thực):")
        output_matrix(A, precision=4)
        run_all_cases(A, x0, precision=7)

        # ==== DEMO: Ma trận phức 4x4 (PWDF_input_A3.txt) ====
        A3 = input_matrix('PWDF_input_A3.txt', convert_fractions=False)
        x0_3 = np.array([-1., 1., 0., 0.])
        print("\n\nMa trận A3 (4x4 - trị riêng phức liên hợp):")
        output_matrix(A3, precision=4)
        run_all_cases(A3, x0_3, precision=7)
    print(f"Đã ghi kết quả vào {output_path}")




