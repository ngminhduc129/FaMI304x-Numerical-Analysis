# =============================================================================
# pp5_FixedPointIteration.py - Lặp điểm cố định cho hệ tuyến tính
#
# Chức năng: Giải hệ Ax = B bằng cách đưa về dạng lặp
#            x = (I-A)x + B và lặp đến khi hội tụ.
#
# Các hàm chính:
#   check_norm(A)                          - tính chuẩn ma trận
#   fixed_point_matrix_iteration(A,B,x0,q,eps,eta,norm)
#
# Input: Đọc từ FXP_input_A.txt, FXP_input_B.txt, FXP_input_X0.txt
# Cách dùng: python pp5_FixedPointIteration.py
# =============================================================================
from pathlib import Path
from fractions import Fraction
from typing import Tuple, Union
import numpy as np
import pandas as pd
import contextlib
import warnings

__dir__ = Path(__file__).parent.resolve()

pd.set_option('display.precision', 12)  # Increase decimal precision
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
def check_norm (A: np.ndarray):
    tmp = np.abs(A);
    print("Chuẩn cột: ", res := np.max((np.sum(np.abs(A), axis=1) - np.diag(np.abs(A))) / np.diag(np.abs(A))))
    print("Chuẩn hàng: ", res := np.max((np.sum(np.abs(A), axis=0) - np.diag(np.abs(A))) / np.diag(np.abs(A))))
    print("Chuẩn 2: ", res := np.linalg.norm(tmp, ord=2))
    print("Chuẩn max: ", res := 3*np.max(tmp))
def fixed_point_matrix_iteration(A, B, initial_values, q, eps, eta, norm, max_iter=1000):
    """
    Perform fixed-point iteration for system of linear equations using x_new = Ax + B
    
    Parameters:
    A (numpy.ndarray): Coefficient matrix
    B (numpy.ndarray): Constant vector
    initial_values (numpy.ndarray): Initial guess for x
    q (float): Contraction coefficient
    eps (float): Tolerance for convergence
    norm (int): Type of norm to use (1 for L1 norm, -1 for max norm)
    
    Returns:
    numpy.ndarray: Solution vector
    """
    if (eps is None) == (eta is None):
        raise ValueError("Specify exactly one of eps (exact) or eta (relative)")
    
    if norm == 1:
        vec_norm = lambda x: np.sum(np.abs(x))
    elif norm == -1:
        vec_norm = lambda x: np.max(np.abs(x))
    
    n = len(initial_values)
    values = np.array(initial_values)
    results = [[0] + initial_values.tolist()]
    iteration = 0

    #Calculate the shrinking speed q:
    new_eps = (eps if eps is not None else eta) * (1-q) / q
    print(f"new_eps: {new_eps:.12f}")

    while True:
        # Calculate new values using matrix operation x_new = Ax + B
        new_values = np.dot(A, values) + B
        
        # Calculate the differences
        if eps is not None:
            total_diff = vec_norm(new_values - values)
        else:
            total_diff = vec_norm(new_values - values) / vec_norm(new_values)
        
        # Append results
        results.append(new_values.tolist() + [total_diff])
        values = new_values
        iteration += 1
        
        # Check for convergence
        if total_diff < new_eps:
            break
        if iteration >= max_iter:
            break
        
    columns = [f'x{j+1}' for j in range(n)] + ['total_diff']
    df = pd.DataFrame(results, columns=columns)
    df.index.name = "Lần lặp"
    return df
if __name__ == "__main__":
    output_path = str(__dir__ / "pp5_FixedPointIteration_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        #Original matrix Ax=B
        A = input_matrix(str(__dir__ / 'FXP_input_A.txt'), convert_fractions=True)
        B = input_matrix(str(__dir__ / 'FXP_input_B.txt'), convert_fractions=True).flatten() #remove flatten if B is multi-column matrix
        initial_guess = input_matrix(str(__dir__ / 'FXP_input_X0.txt'), convert_fractions=True).flatten() #remove flatten if initial guess is multi-column matrix

        print("\nMa trận A:"); output_matrix(A)
        check_norm(A)
        print("\nMa trận B:"); output_matrix(B)
        print("\nGiá trị khởi tạo:"); output_matrix(initial_guess)
        # Initial guess (zero vector or any other initial guess)
        q = np.linalg.norm(np.abs(A), ord=1)
        eps = 1e-6
        eta = None
        norm = 1

        # Call the function
        with np.errstate(all='ignore'):
            try:
                df_history = fixed_point_matrix_iteration (A, B, initial_guess, q, eps, eta, norm)
                print(df_history)
                solution_series = df_history.filter(regex=r'^x\d+$').iloc[-1]
                print("Nghiệm xấp xỉ:"),
                print(solution_series.to_string())
            except Exception as e:
                print(f"Lỗi: {e}")
    print(f"Đã ghi kết quả vào {output_path}")




