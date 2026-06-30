# =============================================================================
# pp7_Seidel.py - Phương pháp lặp Gauss-Seidel
#
# Chức năng: Giải hệ Ax = B bằng phương pháp Gauss-Seidel (cập nhật
#            giá trị mới ngay trong từng bước, hội tụ nhanh hơn Jacobi).
#
# Các hàm chính:
#   convert_to_iteration(A, B)
#   fixed_point_gauss_seidel_2(T, C, D, x0, domiType, eps, eta)
#
# Input: Đọc từ GS_input_A1.txt, GS_input_B1.txt
# Cách dùng: python pp7_Seidel.py
# =============================================================================
from fractions import Fraction
from typing import Tuple, Union
import numpy as np
import pandas as pd

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
def fixed_point_gauss_seidel(
    C: np.ndarray,
    D: np.ndarray,
    x0: Union[np.ndarray, list],
    domiType: int,
    eps: float,
    eta: float
) -> pd.DataFrame:
    """
    Performs Gauss–Seidel fixed-point iteration for row-dominant systems.
    Iteration form: x_new = L x_new + U x_old + D, where C = L + U and diag(C)=0.

    Parameters:
    - C: (n×n) iteration matrix
    - D: (n,) constant vector
    - x0: initial guess vector (length n)
    - domiType: 1 or 3 for row-dominance
    - eps: absolute error tolerance

    Returns:
    - DataFrame of iterates with columns x1...xn and 'error' (∞-norm of x_new−x_old).
    """

    if (eps is None) == (eta is None):
        raise ValueError("Specify exactly one of eps (exact) or eta (relative)")
    
    '''# Check dominance type
    if domiType not in (1, 3):
        raise ValueError("domiType must be 1 or 3 for row-dominant systems")
    else:'''
    vec_norm = lambda x: np.max(np.abs(x))

    # Split C into strict lower (L) and strict upper (U) parts
    n = C.shape[0]
    L = np.tril(C, k=-1)
    U = np.triu(C, k=1)

    # Compute convergence factor q = max_row_k sum(|L[k]|)/(1 - sum(|U[k]|))
    row_L = np.sum(np.abs(L), axis=1)
    row_U = np.sum(np.abs(U), axis=1)
    ratios = row_L / (1 - row_U)
    q = np.max(ratios)
    
    #Compute s
    s=0

    # Tolerance
    tol = (eps if eps is not None else eta) * (1 - q) * (1 - s) / q
    print(f"q: {q:.12f}, s: {s:.12f}, threshold: {tol:.12f}")

    # Initialize
    x_old = np.array(x0, dtype=float)
    history = [x_old.copy()]
    errors = [np.nan]

    # Iterative updates
    while True:
        x_new = np.zeros_like(x_old)
        for k in range(n):
            x_new[k] = (L[k].dot(x_new)) + (U[k].dot(x_old)) + D[k]

        if eps is not None:
            err = vec_norm(x_new - x_old)
        else:
            err = vec_norm(x_new - x_old) / vec_norm(x_new)

        history.append(x_new.copy())
        errors.append(err)
        x_old = x_new

        if err <= tol:
            break
        
    # Prepare DataFrame
    cols = [f"x{i+1}" for i in range(n)]
    df = pd.DataFrame(history, columns=cols)
    df["error"] = errors
    df.index.name = "Iteration"
    return df
#Original matrix Ax=B
A = input_matrix('GS_input_A1.txt', convert_fractions=False)
B = input_matrix('GS_input_B1.txt', convert_fractions=False) #remove flatten if B is multi-column matrix

print("\nMatrix A:"); output_matrix(A)
print("\nCheck dominace of A:", check_dominance(A));
print("\nMatrix B:"); output_matrix(B)
#Convert to recursion form x_new = Cx+D
C, D = convert_to_iteration(A, B)

print("\nMatrix C:"); output_matrix(C)
print("\nMatrix D:"); output_matrix(D)
C = input_matrix('GS_input_A1.txt', convert_fractions=False)
D = input_matrix('GS_input_B1.txt', convert_fractions=False) #remove flatten if B is multi-column matrix

print("\nMatrix C:"); output_matrix(C);
print("\nCheck dominace of C:", check_dominance(C));
print("\nMatrix D:"); output_matrix(D);
#Calculate the result
x0 = np.array([0,0,0,0,0]) #initial value
domiType = check_dominance(C)
eps = 1e-10
eta = None

print("\nMatrix C:"); output_matrix(C);
print("\nCheck dominace of C:", check_dominance(C));
print("\nMatrix D:"); output_matrix(D);
df_history = fixed_point_gauss_seidel(C, D, x0, domiType, eps, eta)
print(df_history.to_string(float_format="{: .4f}".format))

solution_series = df_history.filter(regex=r'^x\d+$').iloc[-1]
print("Approximate solution:"),
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
def fixed_point_gauss_seidel_2(
    T: np.ndarray,
    C: np.ndarray,
    D: np.ndarray,
    x0: Union[np.ndarray, list],
    domiType: int,
    eps: float,
    eta: float
) -> pd.DataFrame:
    """
    Performs Gauss–Seidel fixed-point iteration for column-dominant systems.
    Iteration form: x_new = L x_new + U x_old + D, where C = L + U and diag(C)=0.

    Parameters:
    - C: (n×n) iteration matrix
    - D: (n,) constant vector
    - x0: initial guess vector (length n)
    - domiType: 2 or 3 for column-dominance
    - eps: absolute error tolerance

    Returns:
    - DataFrame of iterates with columns x1...xn and 'error' (∞-norm of x_new−x_old).
    """

    if (eps is None) == (eta is None):
        raise ValueError("Specify exactly one of eps (exact) or eta (relative)")
    
    # Check dominance type
    if domiType not in (2, 3):
        raise ValueError("domiType must be 2 for column-dominant systems")
    else:
        vec_norm = lambda x: np.sum(np.abs(x))

    # Split C into strict lower (L) and strict upper (U) parts
    n = C.shape[0]
    L = np.tril(C, k=-1)
    U = np.triu(C, k=1)

    # Convergence factor q = max_col_k sum(|L[k]|)/(1 - sum(|U[k]|))
    col_L = np.sum(np.abs(L), axis=0)
    col_U = np.sum(np.abs(U), axis=0)
    ratios = col_U / (1 - col_L)
    q = np.max(ratios)

    # Compute s = max_col_k (sum|L[k])
    s = np.max(col_L)

    #Tolerance
    tol = (eps if eps is not None else eta) * (1 - q) * (1 - s) / q
    print(f"q: {q:.12f}, s: {s:.12f}, threshold: {tol:.12f}")

    # Initial y
    T_inv = np.diag(1.0 / np.diag(T))
    x_old = np.array(x0, dtype=float).flatten()
    y_old = T_inv.dot(np.array(x0, dtype=float).flatten())
    history = [np.concatenate([y_old, x_old])]
    errors = [np.nan]
    
    # Iterative updates
    while True:
        y_new = y_old.copy()
        for k in range(n):
            y_new[k] = L[k].dot(y_new) + U[k].dot(y_old) + D[k]

        x_new = T.dot(y_new)
        if eps is not None:
            err = vec_norm(x_new - x_old)
        else:
            err = vec_norm(x_new - x_old) / vec_norm(x_new)
            
        # Record both y and x
        history.append(np.concatenate([y_new, x_new]))
        errors.append(err)
        y_old = y_new; x_old = x_new;

        if err <= tol:
            break
            
    # Prepare DataFrame
    df = pd.DataFrame(history, columns=[f"y{i+1}" for i in range(n)] + [f"x{i+1}" for i in range(n)])
    df['error'] = errors
    df.index.name = 'Iteration'
    return df

#Original matrix Ax=B
A = input_matrix('GS_input_A1.txt', convert_fractions=False)
B = input_matrix('GS_input_B1.txt', convert_fractions=False) #remove flatten if B is multi-column matrix

print("\nMatrix A:"); output_matrix(A)
print("\nCheck dominace of A:", check_dominance(A));
print("\nMatrix B:"); output_matrix(B)
#Convert to recursion form x_new = Cx+D
T, C, D = convert_to_iteration_2(A, B)

print("\nMatrix C:"); output_matrix(C)
print("\nMatrix D:"); output_matrix(D)
print("\nMatrix T:"); output_matrix(T)
#Calculate the result
domiType = check_dominance(A)
x0 = [1,1,1,1,1,1,1] #initial value
eps = 1e-6
eta = None

df_history = fixed_point_gauss_seidel_2(T, C, D, x0, domiType, eps, eta)
print(df_history)

solution_series = df_history.filter(regex=r'^x\d+$').iloc[-1]
print("Approximate solution:"),
print(solution_series.to_string())




