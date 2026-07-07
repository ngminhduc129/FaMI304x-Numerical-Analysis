# =============================================================================
# pp6_Jacobi.py - Phương pháp lặp Jacobi
#
# Chức năng: Giải hệ phương trình tuyến tính Ax = B bằng phương
#            pháp lặp Jacobi (tách đường chéo).
#
# Các hàm chính:
#   check_dominance(A)                     - kiểm tra chéo trội
#   convert_to_iteration(A, B)             - đưa về dạng lặp
#   fixed_point_matrix_iteration_2(T,C,D,x0,domiType,eps,eta)
#
# Input: Đọc từ JCB_input_A1.txt, JCB_input_B1.txt
# Cách dùng: python pp6_Jacobi.py
# =============================================================================
from fractions import Fraction
from typing import Tuple, Union
from pathlib import Path
import numpy as np
import pandas as pd

__dir__ = Path(__file__).parent.resolve()
INPUT_A = str(__dir__ / 'JCB_input_A1.txt')
INPUT_B = str(__dir__ / 'JCB_input_B1.txt')

k = 4  # số lần lặp tối đa

pd.set_option('display.precision', 7)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def input_matrix(filename, convert_fractions=False):
    """
    Reads a matrix from a text file and returns it as a NumPy array.
    Supports fractional entries if present.
    """
    matrix = []

    with open(filename, 'r') as f:
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
def output_matrix(X: np.ndarray, precision: int = 12):
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
    # df = df.round(precision)
    # Print without index/header for cleaner look
    print(df.to_string(index=False, header=False))
def check_dominance(A: np.ndarray) -> int:
    """
    Returns:
      1 if A is row-dominant only,
      2 if A is column-dominant only,
      3 if both,
      0 if neither.
    """
    n = A.shape[0]
    row_dom = all(abs(A[i,i]) > np.sum(np.abs(A[i,:])) - abs(A[i,i]) for i in range(n))
    col_dom = all(abs(A[j,j]) > np.sum(np.abs(A[:,j])) - abs(A[j,j]) for j in range(n))
    if row_dom and col_dom:
        return 3
    if row_dom:
        return 1
    if col_dom:
        return 2
    return 0
def convert_to_iteration(A: np.ndarray, B: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Converts the system Ax = B into the iterative form x_new = C x + D.
    """

    # Step 1: Build diagonal matrix T where T[i,i] = 1 / A[i,i]
    T = np.diag(1.0 / np.diag(A))
    
    # Step 2: Compute C = I - T * A
    n = A.shape[0]
    C = np.eye(n) - T.dot(A)
    
    # Step 3: Compute D = T * B
    D = T.dot(B)
    
    return C, D
def fixed_point_matrix_iteration(
    A: np.ndarray,
    C: np.ndarray,
    D: np.ndarray,
    x0: Union[np.ndarray, list],
    domiType: int,
    eps: float,
    eta: float,
    max_iter: int = None
) -> pd.DataFrame:
    """
    Performs fixed-point iteration x_{k+1}=C x_k + D and returns history with errors.
    Convergence based on:
      - exact error if eps is provided (||x_new - x_old|| <= eps*(1-q)/(α*q))
      - relative error if eta is provided (||x_new-x_old||/||x_new|| <= eta*(1-q)/(α*q))

    Parameters:
    - C, D: iteration matrix and constant
    - x0: initial guess
    - domiType: 1 or 3 (row), 2 (column)
    - eps: absolute tolerance (use exact error)
    - eta: relative tolerance (use ||Δx||/||x_new||)
    """
    # Validate tolerance inputs
    if (eps is None) == (eta is None):
        raise ValueError("Specify exactly one of eps (exact) or eta (relative)")

    # Set norms and α, q based on dominance type
    if domiType in (1, 3):
        vec_norm = lambda x: np.max(np.abs(x))
        q = np.max(np.sum(np.abs(C), axis=1))
        alpha = 1.0
    else:
        raise ValueError("domiType must be 1 or 3 for row-dominant iteration")

    # Convergence threshold
    tol = (eps if eps is not None else eta) * (1 - q) / (alpha * q)
    print(f"Hệ số co (q): {q:.12f}, Hệ số phụ (α): {alpha:.12f}, Ngưỡng dừng: {tol:.12f}")

    # Initialize
    x_old = np.array(x0, dtype=float).flatten()
    history = [x_old.copy()]
    errors = [np.nan]
    iteration = 0

    # Iteration loop
    while True:
        x_new = C.dot(x_old) + D
        if eps is not None:
            err = vec_norm(x_new - x_old)
        else:
            err = vec_norm(x_new - x_old) / vec_norm(x_new)

        # Prepare DataFrame for display
        history.append(x_new.copy())
        errors.append(err)
        x_old = x_new
        iteration += 1
        if err <= tol:
            break
        if max_iter is not None and iteration >= max_iter:
            break
        
    # Compile DataFrame
    n = len(x_new)
    df = pd.DataFrame(history, columns=[f"x{i+1}" for i in range(n)])
    df["sai số"] = errors
    df.index.name = "Lần lặp"
    return df
#Original matrix Ax=B
A = input_matrix(INPUT_A, convert_fractions=False)
B = input_matrix(INPUT_B, convert_fractions=False).flatten()

print("\nMa trận A:"); output_matrix(A)
print("\nKiểm tra chéo trội của A:", check_dominance(A));
print("\nMa trận B:"); output_matrix(B)

domiType = check_dominance(A)

if domiType in (1, 3):
    print("\n=== Trường hợp chéo trội hàng ===")
    C, D = convert_to_iteration(A, B)
    print("\nMa trận C:"); output_matrix(C)
    print("\nMa trận D:"); output_matrix(D)
    x0 = np.zeros(A.shape[0])
    df_history = fixed_point_matrix_iteration(A, C, D, x0, domiType, eps=None, eta=1e-6, max_iter=k)
    print(df_history)
    solution_series = df_history.filter(regex=r'^x\d+$').iloc[-1]
    print("Nghiệm xấp xỉ:"),
    print(solution_series.to_string())
def convert_to_iteration_2(A: np.ndarray, 
                           B: np.ndarray, 
                          ) -> Tuple[np.ndarray, np.ndarray]:
    """
    Prepares iteration y_{k+1} = C y_k + D for column-dominant case.
    """

    # 1. T = diag of 1 / diagonal entries of A
    T = np.diag(1.0 / np.diag(A))

    # 2. C = I - A * T
    n = A.shape[0]
    C = np.eye(n) - A.dot(T)
    
    # 3. D = B
    D = B.copy()

    return T, C, D
def fixed_point_matrix_iteration_2(
    A: np.ndarray,
    T: np.ndarray,
    C: np.ndarray, #main input
    D: np.ndarray,
    x0: Union[np.ndarray, list],
    domiType: int,
    eps: float,
    eta: float,
    max_iter: int = None
) -> np.ndarray:
    """
    Iterates y_{k+1} = C y_k + D for column-dominant systems.
    
    Parameters:
    - C: (n×n) iteration matrix
    - D: (n,) constant vector
    - y0: initial guess (n,)
    - domiType: dominance type (must be 2 or 3 for column-dominant)
    - eps: relative error tolerance
    
    Uses:
    - vector norm: 1-norm
    - matrix norm: column-sum norm
    - alpha = max|diag(A)| / min|diag(A)| from original A
    
    Prints each iteration's y and error, stops when
      ||y_new - y_old||_1 <= eps * (1 - q) / (alpha * q)
    
    Returns:
    - y_new: approximate solution vector
    """
     # Validate tolerance inputs
    if (eps is None) == (eta is None):
        raise ValueError("Specify exactly one of eps (exact) or eta (relative)")

    # Set norms and α, q based on dominance type
    if domiType in (2,3):
        vec_norm = lambda x: np.sum(np.abs(x))
        q = np.max(np.sum(np.abs(C), axis=0))  # column-sum norm
        diag_vals = np.abs(np.diag(A))         # A should be defined globally
        alpha = np.max(diag_vals) / np.min(diag_vals)
    else:
        raise ValueError("domiType must be 2 for column-dominant iteration")

    # Convergence threshold
    tol = (eps if eps is not None else eta) * (1 - q) / (alpha * q)
    print(f"Hệ số co (q): {q:.12f}, Hệ số phụ (α): {alpha:.12f}, Ngưỡng dừng: {tol:.12f}")

    # Initialize
    x_old = np.array(x0, dtype=float).flatten()
    T_inv = np.diag(np.diag(A)); y_old = T_inv.dot(x_old)
    history = [np.concatenate([y_old, x_old])]
    errors = [np.nan]
    iteration = 0

    while True:
        y_new = C.dot(y_old) + D
        x_new = T.dot(y_new)
        if eps is not None:
            err = vec_norm(x_new - x_old)
        else:
            err = vec_norm(x_new - x_old) / vec_norm(x_new)

        # Prepare DataFrame for display
        history.append(np.concatenate([y_new, x_new]))
        errors.append(err)
        y_old = y_new; x_old = x_new;
        iteration += 1
        if err <= tol:
            break
        if max_iter is not None and iteration >= max_iter:
            break
        
    # Compile DataFrame
    n = len(x_new)
    df = pd.DataFrame(history, columns=[f"y{i+1}" for i in range(n)] + [f"x{i+1}" for i in range(n)])
    df["sai số"] = errors
    df.index.name = "Lần lặp"
    return df
#Original matrix Ax=B
A = input_matrix(INPUT_A, convert_fractions=False)
B = input_matrix(INPUT_B, convert_fractions=False).flatten()

print("\nMa trận A:"); output_matrix(A)
print("\nKiểm tra chéo trội của A:", check_dominance(A));
print("\nMa trận B:"); output_matrix(B)

domiType = check_dominance(A)

if domiType in (2, 3):
    print("\n=== Trường hợp chéo trội cột ===")
    T, C, D = convert_to_iteration_2(A, B)
    print("\nMa trận C:"); output_matrix(C)
    print("\nMa trận D:"); output_matrix(D)
    x0 = np.zeros(A.shape[0])
    df_history = fixed_point_matrix_iteration_2(A, T, C, D, x0, domiType, eps=1e-6, eta=None, max_iter=k)
    print(df_history)
    solution_series = df_history.filter(regex=r'^x\d+$').iloc[-1]
    print("Nghiệm xấp xỉ:"),
    print(solution_series.to_string())




