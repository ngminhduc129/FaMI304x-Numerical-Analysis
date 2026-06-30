# =============================================================================
# pp2b_NewtonFixedGap.py - Nội suy Newton bước đều
#
# Chức năng: Nội suy Newton cho các điểm cách đều (bước h = const)
#            dùng sai phân hữu hạn (finite differences).
#
# Các hàm chính:
#   EvenDifference(points, condition)  - bảng sai phân
#   NewtonInterpolation(points, condition=1) - đa thức Newton
#   all_derivatives(coeffs, t, h)     - đạo hàm
#
# Cách dùng: python pp2b_NewtonFixedGap.py
# =============================================================================
import numpy as np
import pandas as pd
import math
from typing import Tuple

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
# --- Helper Function to Get Difference Table and Coefficients (MODIFIED) ---
def EvenDifference(points, condition=1):
    """
    Calculates the difference table and extracts coefficients.

    Args:
        points (list of tuples): A list of (x, y) data points.
        condition (int): 1 for Forward (top diagonal), 0 for Backward (bottom diagonal).

    Returns:
        tuple: A tuple containing:
            - pandas.DataFrame: The full difference table.
            - list: The list of extracted coefficients.
    """
    x_values = [p[0] for p in points]
    y_values = [p[1] for p in points]
    n = len(y_values)

    # Internal calculation table
    diff_calc_table = np.full((n, n), np.nan)
    diff_calc_table[:, 0] = y_values
    for j in range(1, n):
        for i in range(n - j):
            diff_calc_table[i, j] = diff_calc_table[i+1, j-1] - diff_calc_table[i, j-1]

    # --- NEW: Format the output DataFrame ---
    data = {
        'x_i': x_values,
        'y_i': y_values
    }
    for j in range(1, n):
        data[f'Order {j}'] = diff_calc_table[:, j]
    
    df = pd.DataFrame(data)

    # Extract coefficients (logic is unchanged)
    if condition == 1:
        coefficients = diff_calc_table[0, :].tolist()
    else:
        coeffs = []
        for j in range(n):
            coeffs.append(diff_calc_table[n-1-j, j])
        coefficients = coeffs

    return df, coefficients
def NewtonInterpolation(points, condition=1):
    """
    Constructs the Newton interpolation polynomial and returns the steps
    and final coefficients as polished Pandas DataFrames.
    """
    n = len(points) - 1
    if n < 0:
        return pd.DataFrame(), pd.DataFrame()

    diff_table, diff_coeffs = EvenDifference(points, condition=condition)

    # --- NEW: Polished data structures for the output DataFrames ---
    steps_data = []
    N_coeffs_var = np.zeros(n + 1, dtype=float)
    B_coeffs_prev = np.array([1.0])

    for i in range(n + 1):
        D_i = diff_coeffs[i] / math.factorial(i)
        if i == 0:
            B_coeffs = np.array([1.0], dtype=float)
        else:
            k = i - 1
            varB = np.concatenate(([0.0], B_coeffs_prev))
            kB = np.concatenate((k * B_coeffs_prev, [0.0]))
            B_coeffs = varB - kB if condition == 1 else varB + kB
        
        B_coeffs_prev = B_coeffs.copy()
        Ni_coeffs = D_i * B_coeffs
        N_coeffs_var[:len(Ni_coeffs)] += Ni_coeffs

        # Append the polished row of data for the steps DataFrame
        steps_data.append({
            'i': i,
            'Diff Coeff': diff_coeffs[i],
            'D_i': D_i,
            'B_i Coeffs': B_coeffs.tolist(),
            'N_i Coeffs': Ni_coeffs.tolist()
        })

    step_pd = pd.DataFrame(steps_data)
    coeff_pd = pd.DataFrame({
        'Degree': np.arange(n + 1),
        'Coeff': N_coeffs_var
    })

    return step_pd, coeff_pd
points = [(1.4, 0.9523), (1.5, 0.9661), (1.6, 0.9753), (1.7, 0.9838), (1.8, 0.9891), (1.9, 0.9928), (2.0, 0.9)]

# Set the central node index
x0_index = 0
x0_val = points[x0_index][0] #-1 if use backward (condition = 0)
h = points[1][0] - points[0][0]
df, coeffs = EvenDifference(points, condition = 1)

df.style
step_df, final_coeff_df = NewtonInterpolation(points, condition=1)

step_df.style
final_coeff_df.style
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

        Actual differentation P'(k) (x) = P'(k) (c) / h^k
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

def all_derivatives(a, c, gap):
    """
    Compute all derivatives p^(i)(c) using repeated Horner division
    and display in transposed table format.
    """
    coeffs = a.copy()
    degree = len(a) - 1
    results = []
    b0_list = []
    derivative_list = []
    actual_dev_list = []

    # Perform repeated synthetic division
    for i in range(degree + 1):
        df, b0, next_coeff, b_all = synthetic_division(coeffs, c)
        results.append(b_all)
        b0_list.append(b0)
        derivative_list.append(b0 * math.factorial(i))
        actual_dev_list.append(b0 * math.factorial(i) / gap**i)
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
    df.loc["p^(i)(x)"] = [None] + actual_dev_list

    # Add a row on top showing the value of c
    df.loc["c"] = [c] + [None] * (df.shape[1] - 1)
    df = df.loc[["c"] + [idx for idx in df.index if idx != "c"]]  # Move row to top

    return df
coeff_list = final_coeff_df['Coeff'].tolist()

x_val = 1.43
t_val = (x_val - x0_val) / h

print(f"x0_val: {x0_val}; gap: {h}")

df2 = all_derivatives(coeff_list, t_val, h)
df2.style




