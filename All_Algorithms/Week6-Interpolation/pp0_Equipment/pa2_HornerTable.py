# =============================================================================
# pa2_HornerTable.py - Sơ đồ Horner (Synthetic Division)
#
# Chức năng: Tính giá trị đa thức và các đạo hàm tại một điểm
#            bằng sơ đồ Horner (chia đa thức tổng hợp).
#
# Các hàm chính:
#   synthetic_division(a, c)          - sơ đồ Horner
#   all_derivatives(a, c)             - tất cả đạo hàm tại c
#   reverse_horner(b, c)              - Horner ngược
#   w_function(x_points)              - hệ số w(x) = Π(x-xi)
#
# Cách dùng: python pa2_HornerTable.py
# =============================================================================
import numpy as np
import pandas as pd
import math
from typing import Tuple

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column

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
# Example polynomial: p(x) = 3x^3 - 2x - 1
a = [-46, -5, 4, -3, 2, 1]   # coefficients of p(x)
c = 2          # value at which to evaluate

df, p_c, q_coeff, b = synthetic_division(a, c)
df.style.hide(axis="index")
df2 = all_derivatives(a, c)
df2.style
def reverse_horner(b, c):
    """
    Reverse Horner method:
    Given b_i coefficients of q(x) and point c,
    reconstruct a_i coefficients of p(x) = (x - c) * q(x)
    """

    n = len(b)
    b = [0] + b;
    a = [0] * (n + 1)

    # highest term: a_n = b_n
    a[-1] = b[-1]

    # compute remaining coefficients
    for i in range(n-1, -1, -1):
        a[i] = b[i] - c * b[i+1]

    # Build DataFrame for clear display
    data = {
        "i": list(range(n, -1, -1)),
        "b_i": b[::-1],
        "b_i*c": [c * b[i] for i in range(n, -1, -1)],
        "a_i = b_i - c*b_(i+1)": a[::-1]  # exclude the last term a_n (it's beyond b_i range)
    }

    df = pd.DataFrame(data)

    return df, a
# Example usage
b = [2, 3, 5]  # coefficients of q(x)
c = 4

df, a = reverse_horner(b, c)
df.style.hide(axis="index")
def w_function(x_points):
    """
    Construct w_{n+1}(x) = Π(x - x_k)
    using Reverse Horner multiplications.
    """
    steps = {}
    
    # Start with first term (x - x_0)
    coeffs = [-x_points[0], 1]
    steps[0] = coeffs[:]
    
    # Multiply recursively for all remaining x_k
    for i in range(1, len(x_points)):
        df, coeffs = reverse_horner(coeffs, x_points[i])
        steps[i] = coeffs[:]
    
    # ---- Create the display table ----
    max_len = max(len(v) for v in steps.values())
    df_data = {}

    for i, coeff_list in steps.items():
        padded = coeff_list + [np.nan] * (max_len - len(coeff_list))
        df_data[f"i={i}"] = padded

    df = pd.DataFrame(df_data)
    df.insert(0, "a_i", [f"coeff_{j}" for j in range(max_len)])
    
    # Add top row for x_k
    df.loc["x_k"] = ["x_k"] + [x for x in x_points] + [np.nan]*(df.shape[1]-len(x_points)-1)
    df = df.loc[["x_k"] + [idx for idx in df.index if idx != "x_k"]]  # Move row to top

    return df, coeffs

x_points = [1.2, 1.5, 1.7, 1.8, 2.1, 2.3]

df, final_coeffs = w_function(x_points)
df.style




