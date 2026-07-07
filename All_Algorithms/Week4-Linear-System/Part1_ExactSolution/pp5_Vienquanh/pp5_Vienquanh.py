# =============================================================================
# pp5_Vienquanh.py - Nghịch đảo ma trận khối (Vienquanh/Banachiewicz)
#
# Chức năng: Tính ma trận nghịch đảo bằng phương pháp đệ quy khối
#            hoặc thông qua (A^T A)^{-1} A^T.
#
# Các hàm chính:
#   block_matrix_recursion(A)  - nghịch đảo đệ quy
#   inverse_via_ata(A)         - nghịch đảo qua (A^T A)^{-1} A^T
#
# Input: Đọc từ file BLMT_input_A.txt
# Cách dùng: python pp5_Vienquanh.py
# =============================================================================
from pathlib import Path
from fractions import Fraction
from typing import Tuple, Union
import numpy as np
import pandas as pd
import contextlib

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
    df = df.round(precision)
    # Print without index/header for cleaner look
    print(df.to_string(index=False, header=False))
def compute_theta(A: np.ndarray, B_prev: np.ndarray, k: int) -> float:
    row = A[k-1, :k-1]
    col = A[:k-1, k-1]
    a_kk = A[k-1, k-1]
    theta = float(row @ B_prev @ col - a_kk)
    if np.isclose(theta, 0):
        raise ZeroDivisionError(f"θ_{k} is zero; block singular at size {k}.")
    return theta
def compute_b_nn(A: np.ndarray, B_prev: np.ndarray, k: int) -> float:
    # b_{k,k} = -1/θ_k
    theta = compute_theta(A, B_prev, k)
    return -1.0 / theta
def compute_beta_col(A: np.ndarray, B_prev: np.ndarray, k: int) -> np.ndarray:
    # β_{1,k-1}: last column block (excluding b_{k,k})
    alpha_col = A[:k-1, k-1]
    theta = compute_theta(A, B_prev, k)
    return B_prev.dot(alpha_col) / theta
def compute_beta_row(A: np.ndarray, B_prev: np.ndarray, k: int) -> np.ndarray:
    # β_{k-1,1}: last row block (excluding b_{k,k})
    alpha_row = A[k-1, :k-1]
    theta = compute_theta(A, B_prev, k)
    return alpha_row.dot(B_prev) / theta
def update_top_left_block(A: np.ndarray, B_prev: np.ndarray, k: int, beta_row: np.ndarray) -> np.ndarray:
    # B_new = B_prev @ (I - α_col ⊗ β_row)
    alpha_col = A[:k-1, k-1]
    I = np.eye(k-1)
    return B_prev.dot(I - np.outer(alpha_col, beta_row))
def block_matrix_recursion(A: np.ndarray) -> np.ndarray:
    """
    Computes A^{-1} via block recursive inversion on leading principal minors.
    """
    n = A.shape[0]
    # Base case k=1
    theta1 = A[0,0]
    if np.isclose(theta1, 0):
        raise ZeroDivisionError("θ_1 is zero; A[0,0] is singular.")
    B_prev = np.array([[1.0/theta1]])
    # Recursively build inverse for k=2..n
    for k in range(2, n+1):
        print(f"Theta {k}:", compute_theta(A, B_prev, k))
        # Compute new blocks
        beta_col = compute_beta_col(A, B_prev, k)
        beta_row = compute_beta_row(A, B_prev, k)
        b_kk = compute_b_nn(A, B_prev, k)
        # Update top-left (k-1)x(k-1)
        B_tl = update_top_left_block(A, B_prev, k, beta_row)
        # Assemble B_k
        Bk = np.zeros((k, k))
        Bk[:k-1, :k-1] = B_tl
        Bk[:k-1, k-1] = beta_col
        Bk[k-1, :k-1] = beta_row
        Bk[k-1, k-1] = b_kk
        B_prev = Bk
        print(f"Kích thước {k}:\n", B_prev, "\n")
    return B_prev
def inverse_via_ata(A: np.ndarray) -> np.ndarray:
    """
    Computes A^{-1} indirectly by inverting M = A^T A via block recursion,
    then A^{-1} = M^{-1} A^T.
    """
    M = A.T.dot(A)
    M_inv = block_matrix_recursion(M)
    return M_inv.dot(A.T)
if __name__ == "__main__":
    output_path = str(__dir__ / "pp5_Vienquanh_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        A = input_matrix(str(__dir__ / 'BLMT_input_A.txt'), convert_fractions=False)

        print("\nMa trận A:"); output_matrix(A)
        print("\nMa trận M = A^T * A:"); output_matrix(A.T.dot(A))
        A_inv = inverse_via_ata(A)
        print("Kết quả nghịch đảo của A:\n", A_inv)
    print(f"Đã ghi kết quả vào {output_path}")




