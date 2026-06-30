# =============================================================================
# pp1_LagrangeInterpolation.py - Nội suy Lagrange
#
# Chức năng: Xây dựng đa thức nội suy Lagrange đi qua các điểm
#            (x_i, y_i) đã cho.
#
# Các hàm chính:
#   lagrange_interpolation(points) -> (w_df, sub_df, coeff_table)
#     - points: list các cặp [(x1,y1), (x2,y2), ...]
#
# Cách dùng: python pp1_LagrangeInterpolation.py
# =============================================================================
import numpy as np
import pandas as pd
import math
from typing import Tuple

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
# ------------------------------------------------------------
# Helper function 1: w_function
# Produce coefficient list of w_{n+1}(x) = Π (x - x_i)
# ------------------------------------------------------------
def w_function(x_points):
    coeff = np.array([1.0])  # start with polynomial 1
    for xi in x_points:
        coeff = np.convolve(coeff, np.array([1.0, -xi]))  # multiply by (x - xi)
    return coeff  # highest to lowest degree
# ------------------------------------------------------------
# Helper function 2: c_function
# Produce coefficient list of C_i(x) = w_{n+1}(x) / (x - x_i)
# ------------------------------------------------------------
def c_function(w_coeffs, xi):
    n = len(w_coeffs)
    c_coeffs = np.zeros(n - 1)
    remainder = 0.0
    # Synthetic division by (x - xi)
    for j in range(n - 1):
        if j == 0:
            c_coeffs[j] = w_coeffs[j]
        else:
            c_coeffs[j] = w_coeffs[j] + xi * c_coeffs[j - 1]
    remainder = w_coeffs[-1] + xi * c_coeffs[-1]
    # We expect remainder ≈ 0 if everything is correct
    return c_coeffs


# ------------------------------------------------------------
# Main function: lagrange_interpolation
# ------------------------------------------------------------
def lagrange_interpolation(points):
    x_points = np.array([p[0] for p in points], dtype=float)
    y_points = np.array([p[1] for p in points], dtype=float)
    n = len(points) - 1

    # Step 1: Compute w(x)
    w_coeffs = w_function(x_points)

    w_coeffs_df = pd.DataFrame({'Degree': list(range(len(w_coeffs)-1, -1, -1)), 'Coeff': w_coeffs})

    # Step 2: Compute each L_i(x)
    all_Li = []
    table_data = []

    for i in range(n + 1):
        xi, yi = x_points[i], y_points[i]
        
        # D_i calculation
        denom = np.prod([xi - xk for j, xk in enumerate(x_points) if j != i])
        D_i = yi / denom
        
        # C_i(x)
        C_i = c_function(w_coeffs, xi)
        
        # L_i(x) = D_i * C_i(x)
        L_i = D_i * C_i
        all_Li.append(L_i)
        
        # Store intermediate results
        table_data.append({
            'i': i,
            'x_i': xi,
            'y_i': yi,
            'denom': denom,
            'D_i': D_i,
            'C_i(x) coeffs': C_i.tolist(),
            'L_i(x) coeffs': L_i.tolist()
        })

    sub_lagrange_df = pd.DataFrame(table_data)

    # Step 3: Sum all L_i(x)
    max_len = max(len(Li) for Li in all_Li)
    L_coeffs = np.zeros(max_len)
    for Li in all_Li:
        L_coeffs[-len(Li):] += Li  # align degrees

    coeff_table = pd.DataFrame({
        'Degree': list(range(len(L_coeffs))),
        'Coeff': L_coeffs,
    })

    return w_coeffs_df, sub_lagrange_df, coeff_table
# ------------------------------------------------------------
# Example run
# ------------------------------------------------------------
#points = [(1, 17), (2, 17.5), (3, 76), (4, 210.5), (7, 1970)]
points = [(1.2, 0.892), (1.5, 1.179), (1.7, 1.358), (1.8, 1.445), (2.1, 1.688), (2.3, 1.839)]
w_coeffs_df, sub_lagrange_df, coeff_table = lagrange_interpolation(points)
# W-coeff (using Horner multiplication)

w_coeffs_df.style.hide(axis="index")
sub_lagrange_df.style.hide(axis="index")
#Step 3: L_coeff

coeff_table.style.hide(axis="index")




