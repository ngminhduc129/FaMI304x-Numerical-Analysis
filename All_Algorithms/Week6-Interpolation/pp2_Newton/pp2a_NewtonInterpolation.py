# =============================================================================
# pp2a_NewtonInterpolation.py - Nội suy Newton (chia sai phân)
#
# Chức năng: Xây dựng đa thức nội suy Newton dạng chia sai phân
#            (divided differences) qua các điểm (x_i, y_i).
#
# Các hàm chính:
#   divided_differences(points, condition) - bảng sai phân
#   newton_interpolation(points, condition)- đa thức Newton
#   all_derivatives(coeffs, c)            - đạo hàm tại c
#
# Cách dùng: python pp2a_NewtonInterpolation.py
# =============================================================================
import numpy as np
import pandas as pd
import math
from typing import Tuple

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def divided_differences(points, condition):
    """
    points: list of (x_i, y_i) with x_i strictly increasing
    condition: 1 -> forward (first elements of each column)
               0 -> backward (last elements of each column)
    returns: list of selected divided differences (length = len(points))
    """
    x = np.array([p[0] for p in points], dtype=float)
    y = np.array([p[1] for p in points], dtype=float)
    m = len(points)
    if m == 0:
        return []
    if not np.all(np.diff(x) > 0):
        raise ValueError("x values must be strictly increasing.")
    table = np.full((m, m), np.nan, dtype=float)
    table[:, 0] = y.copy()
    for j in range(1, m):
        for i in range(0, m - j):
            table[i, j] = (table[i+1, j-1] - table[i, j-1]) / (x[i+j] - x[i])

    # build DataFrame for display
    data = {'x_i': x, 'y_i': table[:, 0]}
    for j in range(1, m):
        col_vals = [table[i, j] if i < m - j else np.nan for i in range(m)]
        data[f'Order {j}'] = col_vals
    df = pd.DataFrame(data)

    # extract forward or backward selection
    result = []
    for j in range(m):
        col = table[:m - j, j]
        result.append(col[0] if condition == 1 else col[-1])
    return df, result
points = [(0.87, 2.4), (1.24, 2.2), (2.99, 2.0), (3.67, 1.8), (4.23, 1.6)]

df, result = divided_differences(points, condition = 1)

df.style
print(result)
def newton_interpolation(points, condition):
    """
    Build Newton interpolation polynomial coefficients (lowest -> highest).
    Returns numpy array of coefficients [a0, a1, ..., a_n] (constant first).
    """
    m = len(points)
    if m == 0:
        return np.array([])
    x_arr = np.array([p[0] for p in points], dtype=float)
    if not np.all(np.diff(x_arr) > 0):
        raise ValueError("x values must be strictly increasing for the input points.")
    
    # get divided differences (forward or backward)
    tmp, D_list = divided_differences(points, condition=condition)
    # for backward Newton, reverse x and D so loop is same shape
    if condition == 0:
        D_list = D_list[::-1]
        x_arr = x_arr[::-1]
    N_coeff = np.zeros(1, dtype=float)
    steps = []
    for i in range(m):
        D_i = float(D_list[i])
        # build B_{i-1}(x) using lowest-first coefficients
        if i == 0:
            B = np.array([1.0], dtype=float)
        else:
            B = np.array([1.0], dtype=float)
            for k in range(i):
                # (x - x_k) * B  -> x*B - x_k*B
                xB = np.concatenate(([0.0], B))               # x * B (length len(B)+1)
                aB = np.concatenate((x_arr[k] * B, [0.0]))    # x_k * B padded to same length
                B = xB - aB
        N_i = D_i * B
        # add to total polynomial (pad if needed)
        if len(N_coeff) < len(N_i):
            N_coeff = np.pad(N_coeff, (0, len(N_i) - len(N_coeff)), constant_values=0.0)
        N_coeff[:len(N_i)] += N_i
        steps.append({
            'i': i,
            'D_i': D_i,
            'B_(i-1) coeffs (low->high)': np.round(B, 8).tolist(),
            'N_i coeffs (low->high)': np.round(N_i, 8).tolist()
        })

    step_pd = pd.DataFrame(steps)   
    coeff_pd = pd.DataFrame({'Degree': list(range(len(N_coeff))), 'Coeff': N_coeff})

    return step_pd, coeff_pd

points = [(0.87, 2.4), (1.24, 2.2), (2.99, 2.0), (3.67, 1.8), (4.23, 1.6)]
step_pd, coeff_pd = newton_interpolation(points, condition=1)
step_pd.style
coeff_pd.style
#Horner Test
def synthetic_division(a, c):
    """
    Perform synthetic division for polynomial p(x) with coefficients a,
    evaluated at x = c.

    Parameters:
        a (list[float]): coefficients of p(x) from highest to lowest degree
        c (float): the value to evaluate p(c)

    Returns:
        df (pd.DataFrame): table with columns [i, a_i, b_i*c, b_i]
        p_c (float): value of p(c)
        q_coeff (list[float]): coefficients of q(x) = (p(x) - p(c)) / (x - c)
    """

    n = len(a) - 1
    b = [0.0] * (n + 1)
    bc_values = [""] * (n + 1)

    b[n] = a[n]
    for i in range(n - 1, -1, -1):
        b[i] = a[i] + c * b[i + 1]
        bc_values[i + 1] = b[i + 1] * c

    # Prepare table (i from n to 0)
    df = pd.DataFrame({
        "i": list(range(n, -1, -1)),
        "a_i": [a[i] for i in range(n, -1, -1)],
        "b_i*c": [bc_values[i] for i in range(n, -1, -1)],
        "b_i = a_i + b_(i+1)*c": [b[i] for i in range(n, -1, -1)]
    })

    p_c = b[0]
    q_coeff = b[1:]
    return df, p_c, q_coeff, b

def all_derivatives(a, c):
    """
    Compute all derivatives p^(i)(c) using repeated Horner division
    and display in transposed table format.
    """
    coeffs = a.copy()
    degree = len(a) - 1
    results = []
    b0_list = []
    derivative_list = []

    # Perform repeated synthetic division
    for i in range(degree + 1):
        df, b0, next_coeff, b_all = synthetic_division(coeffs, c)
        results.append(b_all)
        b0_list.append(b0)
        derivative_list.append(b0 * math.factorial(i))
        coeffs = next_coeff
        if len(coeffs) == 0:
            break

    # Pad b_i lists for equal column length
    max_len = max(len(b) for b in results)
    for b in results:
        b.extend([None] * (max_len - len(b)))

    # Create DataFrame horizontally
    df = pd.DataFrame(results).T
    df.columns = [f"i={i}" for i in range(len(results))]

    # Insert first column for original a coefficients
    a_col = a + [None] * (df.shape[0] - len(a))
    df.insert(0, "a_i", a_col)

    # Add b_0 and p^(i)(c) rows
    df.loc["b_0"] = [None] + b0_list
    df.loc["p^(i)(c)"] = [None] + derivative_list

    # Add a row on top showing the value of c
    df.loc["c"] = [c] + [None] * (df.shape[1] - 1)
    df = df.loc[["c"] + [idx for idx in df.index if idx != "c"]]  # Move row to top

    return df
coeff_list = coeff_pd['Coeff'].tolist()
df2 = all_derivatives(coeff_list, 3.5)
df2.style




