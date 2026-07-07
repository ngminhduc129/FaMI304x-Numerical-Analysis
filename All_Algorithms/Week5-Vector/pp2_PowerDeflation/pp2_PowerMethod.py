# =============================================================================
# pp2_PowerMethod.py - Phương pháp lũy thừa (Power Method)
#
# Chức năng: Tìm giá trị riêng trội và vector riêng tương ứng của ma trận
#            bằng phương pháp lũy thừa. Hỗ trợ 3 trường hợp:
#               Case 1: trị riêng trội thực, bội đơn
#               Case 2: hai trị trội bằng trị tuyệt đối, trái dấu
#               Case 3: hai trị trội là phức liên hợp
#
# Input: Đọc từ input.txt
# Cách dùng: python pp2_PowerMethod.py
# =============================================================================
import contextlib
from pathlib import Path
from fractions import Fraction
from typing import List, Tuple, Union
import numpy as np
import pandas as pd

DIR = Path(__file__).parent.resolve()

pd.set_option('display.precision', 12)
pd.set_option('display.width', 300)
pd.set_option('display.max_columns', None)


def input_matrix(filename, convert_fractions=False):
    """
    Reads a matrix from a text file and returns it as a NumPy array.
    Supports fractional entries if present.
    """
    matrix = []

    with open(str(DIR / filename), 'r') as f:
        for line in f:
            tokens = line.strip().split()
            if not tokens:
                continue

            row = []
            for token in tokens:
                if '/' in token:
                    val = Fraction(token)
                    row.append(float(val) if convert_fractions else val)
                else:
                    row.append(float(token))

            matrix.append(row)

    return np.array(matrix, dtype=float)


def output_matrix(X: np.ndarray, precision: int = 7):
    """
    Prints a NumPy array (vector or matrix) in a clean tabular format using pandas.
    """
    if X.ndim == 1:
        df = pd.DataFrame(X, columns=["value"])
    elif X.ndim == 2:
        df = pd.DataFrame(X)
    else:
        raise ValueError("Only 1D or 2D arrays are supported.")

    df = df.round(precision)
    print(df.to_string(index=False, header=False))


def verify_eigenpair(A: np.ndarray, lam, v: np.ndarray) -> float:
    residual = np.linalg.norm(A @ v - lam * v) / np.linalg.norm(A)
    return float(abs(residual))


def power_method(A, x0, tol=1e-6, max_iter=1000, norm_ord=2, precision=7, display=True):
    """
    Power Method for finding the dominant eigenvalue and eigenvector of A.
    Displays each iteration in a pandas DataFrame using the specified norm.

    Parameters:
      A (ndarray): Square matrix (n×n).
      x0 (ndarray): Initial guess vector (n,).
      tol (float): Convergence tolerance for eigenvalue change.
      max_iter (int): Maximum number of iterations.
      norm_ord (int, float, or str): Norm order for both eigenvalue estimation and vector normalization.
      precision (int): Number of decimal places to display.

    Returns:
      eigenvalue (float): Approximate dominant eigenvalue.
      eigenvector (ndarray): Approximate eigenvector corresponding to the dominant eigenvalue.
    """
    x = x0 / np.linalg.norm(x0, ord=norm_ord)

    y0 = x0
    lambda0 = np.linalg.norm(y0, ord=norm_ord)
    x0_norm = y0 / lambda0

    history = [
      [0] + y0.tolist() + [lambda0] + x0_norm.tolist()
    ]

    eigenvalue = lambda0
    n = A.shape[0]

    for k in range(1, max_iter + 1):
        y = A.dot(x)
        lambda_new = np.linalg.norm(y, ord=norm_ord)
        if lambda_new == 0:
            raise ZeroDivisionError("Chuẩn của A*x bằng 0; không thể tiếp tục.")
        x_new = y / lambda_new

        row = [k] + y.tolist() + [lambda_new] + x_new.tolist()
        history.append(row)

        if abs(lambda_new - eigenvalue) < tol:
            break

        x = x_new
        eigenvalue = lambda_new

    cols = ['Lần lặp'] + [f'y{i+1}' for i in range(n)] + ['lambda'] + [f'x{i+1}' for i in range(n)]
    df = pd.DataFrame(history, columns=cols)
    if display:
        print(df.to_string(index=False,
                           float_format=lambda x: f"{x:.{precision}f}"))

    return eigenvalue, x_new


def power_method_case2(A, x0, tol=1e-6, max_iter=1000, norm_ord=np.inf, precision=7, display=True):
    """
    Power Method for when the two largest eigenvalues satisfy |λ1| = |λ2| > |λ3| and λ1 = -λ2.
    Uses even powers of A to separate eigencomponents.

    Returns:
      lambda1 (float): Approximate dominant eigenvalue (positive).
      x1 (ndarray): Eigenvector corresponding to λ1.
      x2 (ndarray): Eigenvector corresponding to λ2 = -λ1.
    """
    n = A.shape[0]

    y0 = x0.astype(float)
    lambda0 = np.linalg.norm(y0, ord=norm_ord)
    lambda2_old = lambda0 ** 2
    x0_norm = y0 / lambda0

    history = []
    history.append(
      [0] + y0.tolist() + [lambda2_old] + x0_norm.tolist() + x0_norm.tolist()
    )

    y_even = x0.astype(float)
    A2 = A @ A
    x1 = np.empty(n)
    x2 = np.empty(n)

    for k in range(1, max_iter + 1):
        if k > 1:
            y_even = A2.dot(y_even)
        y_even2 = A2.dot(y_even)
        j = np.argmax(np.abs(y_even))
        lambda2_new = y_even2[j] / y_even[j]
        lambda1 = np.sqrt(lambda2_new)
        y_odd = A.dot(y_even)
        x1 = y_even + lambda1 * y_odd
        x2 = y_even - lambda1 * y_odd
        x1 = x1 / np.linalg.norm(x1, ord=norm_ord)
        x2 = x2 / np.linalg.norm(x2, ord=norm_ord)

        row = [k] + y_even.tolist() + [lambda2_new] + x1.tolist() + x2.tolist()
        history.append(row)

        if k > 1 and abs(lambda2_new - eigen2) < tol:
            break
        eigen2 = lambda2_new

    cols = ['Lần lặp'] + [f'y{i+1}' for i in range(n)] + ['λ²'] + [f'x1_{i+1}' for i in range(n)] + [f'x2_{i+1}' for i in range(n)]
    df = pd.DataFrame(history, columns=cols)
    if display:
        print(df.to_string(index=False,
                           float_format=lambda x: f"{x:.{precision}f}"))

    lambda1_final = np.sqrt(eigen2)
    return lambda1_final, x1, x2


def power_method_case3(A, x0, tol=1e-6, max_iter=100, norm_ord=np.inf, precision=7, display=True):
    """
    Power‐Method variant for the case where the two largest eigenvalues of A
    are complex conjugates (|λ1| = |λ2| > |λ3| ≥ ...).
    We estimate p and q in the quadratic Z^2 + pZ + q = 0
    and then recover λ1, λ2 and their eigenvectors.

    Returns:
      lambda1, lambda2 (complex): the two conjugate eigenvalues of largest magnitude.
      v1, v2 (ndarray): normalized eigenvectors corresponding to λ1 and λ2.
    """
    n = A.shape[0]
    y_m = x0.astype(complex).reshape(n)
    y_m1 = A.dot(y_m)
    y_m2 = A.dot(y_m1)

    history = []
    M0 = np.array([[y_m1[0], y_m[0]], [y_m1[1], y_m[1]]], dtype=complex)
    b0 = -np.array([y_m2[0], y_m2[1]], dtype=complex)
    p0, q0 = np.linalg.solve(M0, b0)
    history.append([0] + y_m.tolist() + y_m1.tolist() + y_m2.tolist() + [p0, q0])

    p_old, q_old = p0, q0

    for m in range(max_iter):
        i1, i2 = 0, 1
        M = np.array([[y_m1[i1], y_m[i1]],
                      [y_m1[i2], y_m[i2]]], dtype=complex)
        b = -np.array([y_m2[i1], y_m2[i2]], dtype=complex)

        try:
            p, q = np.linalg.solve(M, b)
        except np.linalg.LinAlgError:
            raise RuntimeError(
                f"Hệ 2×2 suy biến tại lần lặp m={m}; "
                "hãy chọn cặp chỉ số khác hoặc kiểm tra A."
            )

        row = [m] + y_m.tolist() + y_m1.tolist() + y_m2.tolist() + [p, q]
        history.append(row)

        if m > 0 and (abs(p - p_old) < tol and abs(q - q_old) < tol):
            break

        p_old, q_old = p, q

        y_m = y_m1
        y_m1 = y_m2
        y_m2 = A.dot(y_m1)

        nrm = float(np.max(np.abs(y_m1)))
        if nrm > 1e8:
            y_m /= nrm
            y_m1 /= nrm
            y_m2 /= nrm

    cols = ['m'] + [f'y_m{i+1}' for i in range(n)] \
           + [f'y_m1_{i+1}' for i in range(n)] \
           + [f'y_m2_{i+1}' for i in range(n)] \
           + ['p', 'q']
    df = pd.DataFrame(history, columns=cols)

    if display:
        print(df.to_string(index=False,
                           float_format=lambda x: f"{x:.{precision}f}"))

    discriminant = p * p - 4.0 * q
    sqrt_disc = np.sqrt(discriminant)
    lambda1 = (-p + sqrt_disc) / 2.0
    lambda2 = (-p - sqrt_disc) / 2.0

    v1 = y_m1 - lambda2 * y_m
    v2 = y_m1 - lambda1 * y_m

    v1 = v1 / np.linalg.norm(v1, ord=norm_ord)
    v2 = v2 / np.linalg.norm(v2, ord=norm_ord)

    return lambda1, lambda2, v1, v2


# ==== DEMO ====
if __name__ == "__main__":
    output_path = str(DIR / "PowerMethod_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        A = input_matrix('input.txt', convert_fractions=False)
        x0 = np.ones(A.shape[0], dtype=float)

        print(f"Ma trận A ({A.shape[0]}x{A.shape[1]}):")
        output_matrix(A, precision=4)
        print()

        print("=" * 70)
        print("CASE 1 - Lũy thừa thường (trị riêng trội thực, bội đơn)")
        print("=" * 70)
        try:
            lam1, v1 = power_method(A, x0, precision=7)
            v1_str = ", ".join(f"{x:.7f}" for x in v1)
            print(f"\n  λ₁ ≈ {lam1:.7f}")
            print(f"  v₁ = [{v1_str}]")
        except Exception as e:
            print(f"  ❌ {e}")

        print()
        print("=" * 70)
        print("CASE 2 - Hai trị trội thực, bằng trị tuyệt đối, trái dấu")
        print("=" * 70)
        try:
            lam2a, v2a, v2b = power_method_case2(A, x0, precision=7)
            r2a = verify_eigenpair(A, lam2a, v2a)
            r2b = verify_eigenpair(A, -lam2a, v2b)
            v2a_str = ", ".join(f"{x:.7f}" for x in v2a)
            v2b_str = ", ".join(f"{x:.7f}" for x in v2b)
            print(f"\n  λ₁ ≈ {lam2a:.7f}, λ₂ ≈ {-lam2a:.7f}")
            print(f"  residual₁ = {r2a:.2e}, residual₂ = {r2b:.2e}")
            print(f"  v₁ = [{v2a_str}]")
            print(f"  v₂ = [{v2b_str}]")
        except Exception as e:
            print(f"  ❌ {e}")

        print()
        print("=" * 70)
        print("CASE 3 - Hai trị trội là phức liên hợp")
        print("=" * 70)
        try:
            lam3a, lam3b, v3a, v3b = power_method_case3(A, x0, precision=7)
            v3a_str = ", ".join(f"{x.real:.7f}{x.imag:+.7f}j" for x in v3a)
            v3b_str = ", ".join(f"{x.real:.7f}{x.imag:+.7f}j" for x in v3b)
            print(f"\n  λ₁ ≈ {lam3a:.7f}, λ₂ ≈ {lam3b:.7f}")
            print(f"  v₁ = [{v3a_str}]")
            print(f"  v₂ = [{v3b_str}]")
        except Exception as e:
            print(f"  ❌ {e}")

    print(f"Đã ghi kết quả vào {output_path}")
