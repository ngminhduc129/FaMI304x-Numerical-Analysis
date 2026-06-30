# =============================================================================
# pp3d_Bessel.py - Nội suy Bessel
#
# Chức năng: Nội suy bằng công thức Bessel dùng sai phân trung tâm,
#            phù hợp khi điểm nội suy ở giữa hai nút.
#
# Các hàm chính:
#   parse_xy_data(filepath)             - đọc dữ liệu từ CSV
#   get_points_in_range(points, x0, n)  - chọn điểm
#   bessel_interpolation(points, x0_index) - công thức Bessel
#   all_derivatives(coeffs, t, h)      - đạo hàm
#
# Cách dùng: python pp3d_Bessel.py
# =============================================================================
import numpy as np
import pandas as pd
import math
import os
from typing import Tuple, List

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
#Reading file

def parse_xy_data(filepath, delimiter=None):
    """
    Reads a CSV-like file with x, y data and returns a list of (x, y) tuples.

    This function is designed to handle different delimiters (like ';' or ' ')
    and assumes that commas (',') are used as decimal separators, based on
    the provided image.

    Args:
        filepath (str): The path to the data file.
        delimiter (str, optional): The column delimiter (e.g., ';', ' '). 
                                   If None, the function will try to 
                                   auto-detect it.

    Returns:
        list: A list of (x, y) float tuples.
              Returns an empty list if the file cannot be read or is empty.
    """
    data_points = []
    detected_delimiter = delimiter
    
    # --- 1. Delimiter Sniffing (if not provided) ---
    if detected_delimiter is None:
        try:
            with open(filepath, 'r', encoding='utf-8') as f:
                # Read the first non-empty line to guess
                first_line = ""
                for line in f:
                    first_line = line.strip()
                    if first_line:
                        break
                
                if ';' in first_line:
                    detected_delimiter = ';'
                elif ' ' in first_line:
                    # Check if it's likely a space delimiter
                    parts = re.split(r'\s+', first_line)
                    if len(parts) == 2:
                        try:
                            # Try to parse to see if it makes sense
                            float(parts[0].replace(',', '.'))
                            float(parts[1].replace(',', '.'))
                            detected_delimiter = ' '
                        except (ValueError, IndexError):
                             # Not a valid 2-column space-delimited float line
                             pass
                
                if detected_delimiter is None and ',' in first_line:
                    # Comma is the last guess, as it's ambiguous with decimal
                    detected_delimiter = ','
                
                if detected_delimiter is None:
                    # Final fallback based on your image
                    print("Warning: Could not auto-detect delimiter. Falling back to ';'.")
                    detected_delimiter = ';'
        except Exception as e:
            print(f"Error opening/reading file for sniffing: {e}")
            return [] # Return empty list on error
    
    print(f"Using delimiter: '{detected_delimiter}'")

    # --- 2. File Parsing ---
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line_number, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue # Skip empty lines or comment lines

                # Split the line by the detected delimiter
                if detected_delimiter == ' ':
                    # Use regex split for spaces to handle multiple spaces
                    parts = re.split(r'\s+', line)
                else:
                    parts = line.split(detected_delimiter)

                # Ensure we have exactly two columns
                if len(parts) == 2:
                    x_str, y_str = parts
                    
                    try:
                        # KEY STEP: Replace comma with dot for float conversion
                        x_val = float(x_str.strip().replace(',', '.'))
                        y_val = float(y_str.strip().replace(',', '.'))
                        data_points.append((x_val, y_val))
                    except ValueError as e:
                        # Warn if conversion to float fails
                        print(f"Warning: Could not parse numbers on line {line_number}: '{line}'. Error: {e}")
                else:
                    # Warn if the line doesn't have exactly two parts
                    print(f"Warning: Skipping malformed line {line_number}: '{line}'. Expected 2 columns, found {len(parts)}")
    
    except FileNotFoundError:
        print(f"Error: File not found at '{filepath}'")
        return []
    except Exception as e:
        print(f"An error occurred while reading the file: {e}")
        return []

    return data_points
print("--- Example: Reading '19data.csv' ---")
    
# You would replace '19data.csv' with the path to your actual file
file_to_read = '20251GKtest2.csv' 
    
# Check if the file exists before trying to read it
if os.path.exists(file_to_read):
    # Call the function to parse the file.
    # It will try to auto-detect the delimiter.
    all_points = parse_xy_data(file_to_read)
        
    if all_points:
        print(f"Successfully parsed {len(all_points)} data points:")
        df_input = pd.DataFrame(all_points, columns=['x', 'y'])
        print(df_input.to_string(index=False))
    else:
        print(f"Could not parse any data points from '{file_to_read}'.")
        print("Please check the file format and warnings above.")
            
else:
    print(f"Error: The file '{file_to_read}' was not found.")
    print("Please create this file or change 'file_to_read' variable")
    print("to point to your existing data file.")
def get_points_in_range(points: List[Tuple[float, float]], 
                        x_start: float, 
                        x_end: float) -> List[Tuple[float, float]]:
    """
    Filters a list of (x, y) tuples to return only those within a 
    specified x-range [x_start, x_end].
    """
    subset_points = [
        (x, y) for (x, y) in points 
        if x_start <= x <= x_end
    ]
    return subset_points

# Select points for the second interval [9.7, 10.0]
# This interval is suitable for Newton's FORWARD method

x_start_demo = 2.150
x_end_demo = 2.955
    
monotonic_points_fwd = get_points_in_range(all_points, x_start_demo, x_end_demo)
    
print(f"Subset of points from x={x_start_demo} to x={x_end_demo}:")
df_subset = pd.DataFrame(monotonic_points_fwd, columns=['x (swap)', 'y (swap)'])
print(df_subset.to_string(index=False))
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

def bessel_interpolation(points, x0_index):
    """
    Constructs the Bessel interpolation polynomial in the variable t = (x - x0) / h.
    
    Args:
        points (list of tuples): List of (x, y) data points. Must be evenly spaced
                                 and have an even number of points.
        x0_index (int): The index of the first central node (x0).

    Returns:
        step_pd (pd.DataFrame): DataFrame showing the intermediate steps.
        coeff_pd (pd.DataFrame): DataFrame of the final polynomial coefficients for P(t).
    """
    
    n = len(points) - 1
    if (n + 1) % 2 != 0:
        raise ValueError("Bessel's formula requires an even number of nodes (odd degree n).")
    
    # --- 1. Generate Difference Table ---
    diff_table_df = EvenDifference(points)

    # --- 2. Select Coefficients (Bessel Path) ---
    C_coeffs = []
    for i in range(n + 1):
        order_col = f'Order {i}' if i > 0 else 'y_i'
        
        if i % 2 == 0: # Even order: C_{2k} = (d^{2k} y_{-k} + d^{2k} y_{-(k-1)}) / 2
            k = i // 2
            c1_idx = x0_index - k
            c2_idx = x0_index - k + 1
            if c1_idx < 0 or c2_idx >= len(points):
                raise IndexError(f"Cannot access indices {c1_idx}, {c2_idx} for order {i}.")
            
            c1 = diff_table_df[order_col].iloc[c1_idx]
            c2 = diff_table_df[order_col].iloc[c2_idx]
            C_coeffs.append((c1 + c2) / 2.0)
            
        else: # Odd order: C_{2k+1} = d^{2k+1} y_{-k}
            k = (i - 1) // 2
            c_idx = x0_index - k
            if c_idx < 0 or c_idx >= len(points):
                raise IndexError(f"Cannot access index {c_idx} for order {i}.")
            C_coeffs.append(diff_table_df[order_col].iloc[c_idx])

    # --- 3. Build Polynomial P(t) ---
    steps_data = []
    N_coeffs_total = np.zeros(n + 1, dtype=float)
    
    # B_i(t) basis polynomials
    B_coeffs_list = []
    
    # B_0(t) = 1
    B_0 = np.array([1.0])
    B_coeffs_list.append(B_0)
    
    # B_1(t) = (t - 0.5)
    B_1 = np.array([-0.5, 1.0]) # [const, t]
    B_coeffs_list.append(B_1)

    for i in range(2, n + 1):
        if i % 2 == 0: # Even: B_{2k} = B_{2k-2} * (t+k-1)(t-k)
            k = i // 2
            # (t - k)
            term1 = np.array([-k, 1.0])
            # (t + k - 1)
            term2 = np.array([k - 1, 1.0])
            # (t^2 - t - k(k-1))
            new_factor = np.convolve(term1, term2)
            # B_{2k-2} is at index i-2
            B_i = np.convolve(B_coeffs_list[i-2], new_factor)
        else: # Odd: B_{2k+1} = B_{2k} * (t - 0.5)
            # B_{2k} is at index i-1
            B_i = np.convolve(B_coeffs_list[i-1], [-0.5, 1.0])
        
        B_coeffs_list.append(B_i)

    # Calculate final polynomial
    for i in range(n + 1):
        D_i = C_coeffs[i] / math.factorial(i)
        B_i = B_coeffs_list[i]
        
        # N_i(t) = D_i * B_i(t)
        Ni_coeffs = D_i * B_i
        
        # Add to total polynomial (pad with zeros)
        N_coeffs_total[:len(Ni_coeffs)] += Ni_coeffs

        # Store intermediate steps for printing
        steps_data.append({
            'i': i,
            'Difference Coeff (C_i)': C_coeffs[i],
            'D_i = C_i / i!': D_i,
            'B_i(t) Coeffs (low->high)': B_i.tolist(),
            'N_i(t) Coeffs (low->high)': Ni_coeffs.tolist()
        })

    step_pd = pd.DataFrame(steps_data)
    coeff_pd = pd.DataFrame({
        'Degree (t)': np.arange(n + 1),
        'Coeff': N_coeffs_total
    })

    return step_pd, coeff_pd
# 1. Define the data points from the image
points = monotonic_points_fwd

# 2. Set the central node index
# 9 points (indices 0-8), the center is index 4
x0_index = 3
x0_val = (points[x0_index][0] + points[x0_index + 1][0])/2
h = points[1][0] - points[0][0]
print("--- Generated Finite Difference Table ---")
df = EvenDifference(points)

df.style
# 3. Calculate the Bessel polynomial
step_df, final_coeff_df = bessel_interpolation(points, x0_index)

print("\n--- Polynomial Construction Steps (Bessel) ---")
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

x_val = 2.55
t_val = (x_val - x0_val) / h

print(f"x0_val: {x0_val}; gap: {h}")

df2 = all_derivatives(coeff_list, t_val, h)
df2.style




