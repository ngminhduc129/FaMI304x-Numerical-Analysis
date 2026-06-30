# =============================================================================
# pp4a_InverseFunction.py - Nội suy hàm ngược
#
# Chức năng: Tìm x từ y cho trước bằng cách nội suy Newton trên
#            hàm ngược (đổi vai trò x và y).
#
# Các hàm chính:
#   parse_xy_data(filepath)               - đọc dữ liệu
#   find_monotonic_intervals(points)      - tìm khoảng đơn điệu
#   newton_interpolation(points, condition)- nội suy Newton
#   all_derivatives(coeffs, c)            - đạo hàm tại c
#
# Cách dùng: python pp4a_InverseFunction.py
# =============================================================================
import numpy as np
import pandas as pd
import math
from typing import Tuple, List
import re
import os

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
## Reading a file

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
file_to_read = '19data.csv' 
    
# Check if the file exists before trying to read it
if os.path.exists(file_to_read):
    # Call the function to parse the file.
    # It will try to auto-detect the delimiter.
    data_list = parse_xy_data(file_to_read)
        
    if data_list:
        print(f"Successfully parsed {len(data_list)} data points:")
        df_input = pd.DataFrame(data_list, columns=['x', 'y'])
        print(df_input.to_string(index=False))
    else:
        print(f"Could not parse any data points from '{file_to_read}'.")
        print("Please check the file format and warnings above.")
            
else:
    print(f"Error: The file '{file_to_read}' was not found.")
    print("Please create this file or change 'file_to_read' variable")
    print("to point to your existing data file.")
def find_monotonic_intervals(points: List[Tuple[float, float]]) -> pd.DataFrame:
    """
    Analyzes a list of (x, y) points and finds all continuous monotonic
    (non-decreasing or non-increasing) intervals.

    Args:
        points (List[Tuple[float, float]]): 
            A list of (x, y) data points. 
            It is ASSUMED that the x-values are already sorted in ascending order.

    Returns:
        pd.DataFrame: A DataFrame detailing each monotonic interval.
    """
    
    if len(points) < 2:
        return pd.DataFrame(columns=["Type", "Start (x, y)", "End (x, y)"])

    intervals = []
    
    # Start with the first point as the beginning of the first interval
    start_point = points[0]
    
    # Determine the initial direction (sign of the slope)
    # 1: increasing, -1: decreasing, 0: flat
    current_direction = np.sign(points[1][1] - points[0][1])
    
    # If the first segment is flat, we'll assign its direction on the first change
    if current_direction == 0:
        current_direction = 1 # Default to "non-decreasing"

    for i in range(1, len(points) - 1):
        x_i, y_i = points[i]
        x_ip1, y_ip1 = points[i+1]
        
        new_direction = np.sign(y_ip1 - y_i)
        
        # A change in direction (from + to -) or (- to +) marks a break.
        # We ignore flat segments (new_direction == 0) as they don't break monotonicity.
        if new_direction != 0 and new_direction != current_direction:
            # --- Trend has broken ---
            
            # 1. Save the previous interval
            interval_type = "Increasing" if current_direction > 0 else "Decreasing"
            intervals.append((interval_type, start_point, points[i]))
            
            # 2. Start a new interval
            start_point = points[i]
            current_direction = new_direction
            
    # After the loop, add the final interval
    interval_type = "Increasing" if current_direction > 0 else "Decreasing"
    intervals.append((interval_type, start_point, points[-1]))

    # Convert to a clean DataFrame for printing
    df = pd.DataFrame(intervals, columns=["Type", "Start Point (x, y)", "End Point (x, y)"])
    return df
intervals_df = find_monotonic_intervals(data_list)

print("Monotonic Intervals Found:")
print(intervals_df.to_string(index=False))
def get_points_in_range(points: List[Tuple[float, float]], 
                        x_start: float, 
                        x_end: float) -> List[Tuple[float, float]]:
    """
    Filters a list of (x, y) tuples to return only those within a 
    specified x-range [x_start, x_end].

    Args:
        points (List[Tuple[float, float]]): The complete list of (x, y) points.
        x_start (float): The starting x-value (inclusive).
        x_end (float): The ending x-value (inclusive).

    Returns:
        List[Tuple[float, float]]: A new list containing only the points
                                     where x_start <= x <= x_end.
    """
    
    # Use a list comprehension for an efficient, one-line filter
    subset_points = [
        (y, x) for (x, y) in points 
        if x_start <= x <= x_end
    ]
    
    return subset_points
x_start_demo = 1.3
x_end_demo = 1.7
    
subset_data = get_points_in_range(data_list, x_start_demo, x_end_demo)
    
print(f"Subset of points from x={x_start_demo} to x={x_end_demo}:")
df_subset = pd.DataFrame(subset_data, columns=['y (swap)', 'x (swap)'])
print(df_subset.to_string(index=False))

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

df, result = divided_differences(subset_data, condition = 1)

df.style
step_pd, coeff_pd = newton_interpolation(subset_data, condition=1)
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
df2 = all_derivatives(coeff_list, 2.5)
df2.style




