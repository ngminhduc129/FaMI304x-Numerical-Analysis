# =============================================================================
# pp3b_GaussII.py - Nội suy Gauss lùi (Gauss Backward)
#
# Chức năng: Nội suy bằng công thức Gauss lùi dùng sai phân trung
#            tâm cho các điểm cách đều.
#
# Các hàm chính:
#   EvenDifference(points)              - bảng sai phân
#   gauss_2_interpolation(points, x0_index) - công thức Gauss lùi
#   all_derivatives(coeffs, t, h)      - đạo hàm
#
# Cách dùng: python pp3b_GaussII.py
# =============================================================================
import numpy as np
import pandas as pd
import math
from typing import Tuple

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def EvenDifference(points):
    """
    Calculates the finite difference table for evenly spaced nodes.
    
    Args:
        points (list of tuples): A list of (x, y) data points, sorted by x.

    Returns:
        pandas.DataFrame: The full difference table.
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

    # Format the output DataFrame
    data = {
        'x_i': x_values,
        'y_i': y_values
    }
    for j in range(1, n):
        # Pad with NaN to maintain table shape
        col_name = f'Order {j}'
        col_data = np.full(n, np.nan)
        col_data[:(n-j)] = diff_calc_table[:(n-j), j]
        data[col_name] = col_data
    
    df = pd.DataFrame(data)
    return df

def gauss_2_interpolation(points, x0_index):
    """
    Constructs the Gauss II interpolation polynomial in the variable t = (x - x0) / h.
    
    Args:
        points (list of tuples): List of (x, y) data points. Must be evenly spaced
                                 and have an odd number of points.
        x0_index (int): The index of the central node (x0).

    Returns:
        step_pd (pd.DataFrame): DataFrame showing the intermediate steps.
        coeff_pd (pd.DataFrame): DataFrame of the final polynomial coefficients for P(t).
    """
    
    n = len(points) - 1
    if (n + 1) % 2 == 0:
        raise ValueError("Gauss II formula requires an odd number of nodes (even degree n).")
    
    # --- 1. Generate Difference Table ---
    # We use the re-usable function from pp2b_NewtonFixedGap.ipynb
    diff_table_df = EvenDifference(points)

    # --- 2. Select Coefficients (Gauss II Path) ---
    C_coeffs = []
    for i in range(n + 1):
        order_col = f'Order {i}' if i > 0 else 'y_i'
        
        # This is the selection rule for Gauss II
        # y0, y-1, y-1, y-2, y-2, ...
        row_index = x0_index - (i + 1) // 2
        
        if row_index < 0 or row_index >= len(points):
             raise IndexError(f"Cannot access index {row_index} in difference table for order {i}.")
             
        C_coeffs.append(diff_table_df[order_col].iloc[row_index])

    # --- 3. Build Polynomial P(t) ---
    steps_data = []
    N_coeffs_total = np.zeros(n + 1, dtype=float)
    B_coeffs = np.array([1.0])  # B_0(t) = 1

    for i in range(n + 1):
        D_i = C_coeffs[i] / math.factorial(i)
        
        # Calculate B_i(t) from B_{i-1}(t)
        if i == 1:
            # B_1(t) = t
            B_coeffs = np.array([0.0, 1.0]) # [const, t]
        elif i > 1:
            if i % 2 == 0: # Even: B_2k = B_{2k-1} * (t + k)
                k = i // 2
                B_coeffs = np.convolve(B_coeffs, [k, 1])
            else: # Odd: B_{2k+1} = B_{2k} * (t - k)
                k = (i - 1) // 2
                B_coeffs = np.convolve(B_coeffs, [-k, 1])
        
        # N_i(t) = D_i * B_i(t)
        Ni_coeffs = D_i * B_coeffs
        
        # Add to total polynomial (pad with zeros)
        N_coeffs_total[:len(Ni_coeffs)] += Ni_coeffs

        # Store intermediate steps for printing
        steps_data.append({
            'i': i,
            'Difference Coeff': C_coeffs[i],
            'D_i = C_i / i!': D_i,
            'B_i(t) Coeffs (low->high)': B_coeffs.tolist(),
            'N_i(t) Coeffs (low->high)': Ni_coeffs.tolist()
        })

    step_pd = pd.DataFrame(steps_data)
    coeff_pd = pd.DataFrame({
        'Degree (t)': np.arange(n + 1),
        'Coeff': N_coeffs_total
    })

    return step_pd, coeff_pd
# 1. Define the data points from the image
points = [
    (9.2, 9.4341319),
    (9.3, 9.4307764),
    (9.4, 9.4261142),
    (9.5, 9.4211191),
    (9.6, 9.4170553),
    (9.7, 9.4147476),
    (9.8, 9.4152900),
    (9.9, 9.4196762),
    (10.0, 9.4288617)
]

# 2. Set the central node index
# 9 points (indices 0-8), the center is index 4
x0_index = 4
x0_val = points[x0_index][0]
h = points[1][0] - points[0][0]
print("--- Generated Finite Difference Table ---")
df = EvenDifference(points)

df.style
# 3. Calculate the Gauss II polynomial
step_df, final_coeff_df = gauss_2_interpolation(points, x0_index)

print("\n--- Polynomial Construction Steps (Gauss II) ---")
step_df.style
print("\n--- Final Polynomial P(t) Coefficients ---")
print(f"P(t) = a0 + a1*t + ... (where t = (x - {x0_val:.1f}) / {h:.1f})")
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

x_val = 9.68
t_val = (x_val - x0_val) / h

print(f"x0_val: {x0_val}; gap: {h}")

df2 = all_derivatives(coeff_list, t_val, h)
df2.style




