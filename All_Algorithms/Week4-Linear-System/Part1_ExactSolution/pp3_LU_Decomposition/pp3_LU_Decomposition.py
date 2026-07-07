# =============================================================================
# pp3_LU_Decomposition.py - Phân rã LU (Doolittle / Crout)
#
# Chức năng: Phân tích ma trận A thành tích L*U với L là ma trận
#            tam giác dưới, U là ma trận tam giác trên.
#
# Các hàm chính:
#   lu_decomposition(A)              - phân rã LU
#   verify_lu_decomposition(A, L, U) - kiểm tra A = L*U
#
# Input: Đọc từ file LU_input_A.txt
# Cách dùng: python pp3_LU_Decomposition.py
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
def lu_decomposition(A):
    """
    Perform LU decomposition on matrix A with U having unit diagonal.
    
    Parameters:
    A (numpy.ndarray): Input square matrix
    
    Returns:
    tuple: (L, U) matrices where U has diagonal elements = 1
    """
    n = len(A)
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    
    # Set diagonal of U to 1
    np.fill_diagonal(U, 1)
    
    # LU decomposition
    for i in range(n):
        # Lower triangular matrix
        L[i][i] = A[i][i] - sum(L[i][j] * U[j][i] for j in range(i))
        
        for k in range(i + 1, n):
            # Lower triangular elements
            L[k][i] = A[k][i] - sum(L[k][j] * U[j][i] for j in range(i))
            
            # Upper triangular elements (normalized)
            U[i][k] = (A[i][k] - sum(L[i][j] * U[j][k] for j in range(i))) / L[i][i]

        print(f"Bước {i + 1}:")
        print(pd.DataFrame(L).to_string(index=False, header=False))
        print(pd.DataFrame(U).to_string(index=False, header=False))
    
    return L, U
def verify_lu_decomposition(A, L, U):
    """
    Verify LU decomposition by multiplying L and U and comparing with original matrix.
    
    Parameters:
    A (numpy.ndarray): Original matrix
    L (numpy.ndarray): Lower triangular matrix
    U (numpy.ndarray): Upper triangular matrix
    """
    print("\nMa trận tam giác dưới L:")
    print(pd.DataFrame(L).to_string(index=False, header=False))
    print("\nMa trận tam giác trên U:")
    print(pd.DataFrame(U).to_string(index=False, header=False))

    A_np = np.array(A).astype(np.float64)
    L_np = np.array(L).astype(np.float64)
    U_np = np.array(U).astype(np.float64)
    
    # Multiply L and U
    LU = np.dot(L_np, U_np)
    print("\nL × U:")
    print(pd.DataFrame(LU).to_string(index=False, header=False))
    
    # Check if LU equals A
    if np.allclose(LU, A_np):
        print("\nKiểm tra thành công: LU = A")
    else:
        print("\nKiểm tra thất bại: LU ≠ A")
if __name__ == "__main__":
    output_path = str(__dir__ / "pp3_LU_Decomposition_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        A = input_matrix(str(__dir__ / 'LU_input_A.txt'), convert_fractions=False)

        print("\nMa trận A:"); output_matrix(A)
        if A is not None:
            L, U = lu_decomposition(A)
            verify_lu_decomposition(A, L, U)
    print(f"Đã ghi kết quả vào {output_path}")




