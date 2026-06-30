# =============================================================================
# pp4b_LoopInverse.py - Nội suy ngược bằng lặp
#
# Chức năng: Tìm x từ y cho trước bằng kết hợp nội suy Newton
#            và phương pháp lặp điểm cố định (fixed-point).
#
# Các hàm chính:
#   inverse_interpolation_newton_fixed_point(points, y_target, eps, condition)
#
# Cách dùng: python pp4b_LoopInverse.py
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
file_to_read = '20251GKtest11.csv' 
    
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
def find_monotonic_intervals(points: List[Tuple[float, float]]) -> pd.DataFrame:
    """
    Analyzes a list of (x, y) points and finds all continuous monotonic
    (non-decreasing or non-increasing) intervals.
    """
    if len(points) < 2:
        return pd.DataFrame(columns=["Type", "Start (x, y)", "End (x, y)"])
    intervals = []
    start_point = points[0]
    current_direction = np.sign(points[1][1] - points[0][1])
    if current_direction == 0:
        current_direction = 1 
    for i in range(1, len(points) - 1):
        x_i, y_i = points[i]
        x_ip1, y_ip1 = points[i+1]
        new_direction = np.sign(y_ip1 - y_i)
        if new_direction != 0 and new_direction != current_direction:
            interval_type = "Increasing" if current_direction > 0 else "Decreasing"
            intervals.append((interval_type, start_point, points[i]))
            start_point = points[i]
            current_direction = new_direction
    interval_type = "Increasing" if current_direction > 0 else "Decreasing"
    intervals.append((interval_type, start_point, points[-1]))
    df = pd.DataFrame(intervals, columns=["Type", "Start Point (x, y)", "End Point (x, y)"])
    return df
print("\n--- Finding Monotonic Intervals ---")
intervals_df = find_monotonic_intervals(all_points)
print(intervals_df.to_string(index=False))
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

x_start_demo = 2.265
x_end_demo = 2.610
    
monotonic_points_fwd = get_points_in_range(all_points, x_start_demo, x_end_demo)
    
print(f"Subset of points from x={x_start_demo} to x={x_end_demo}:")
df_subset = pd.DataFrame(monotonic_points_fwd, columns=['x (swap)', 'y (swap)'])
print(df_subset.to_string(index=False))
def EvenDifference(points, condition):
    """
    Calculates the finite difference table and extracts coefficients.
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
    data = {'x_i': x_values, 'y_i': y_values}
    for j in range(1, n):
        col_name = f'Order {j}'
        col_data = np.full(n, np.nan)
        col_data[:(n-j)] = diff_calc_table[:(n-j), j]
        data[col_name] = col_data
    
    df = pd.DataFrame(data)

    # Extract coefficients
    if condition == 1: # Newton-Gregory Forward (top diagonal)
        coefficients = diff_calc_table[0, :].tolist()
    else: # Newton-Gregory Backward (bottom diagonal)
        coeffs = []
        for j in range(n):
            coeffs.append(diff_calc_table[n-1-j, j])
        coefficients = coeffs

    return df, coefficients

def evaluate_basis_product(t, i, condition):
    if condition == 1:
        """ Helper function to calculate t(t-1)...(t-(i-1)) """
        if i == 0:
            return 1.0
        product = 1.0
        for j in range(i):
            product *= (t - j)
        return product
    
    if condition == 0:
        """ Helper function to calculate t(t+1)...(t+(i-1)) """
        if i == 0:
            return 1.0
        product = 1.0
        for j in range(i):
            product *= (t + j)
        return product
def phi_function(t, y_target, y_tmp, C_coeffs, condition):
    if condition == 1:
        """
        Calculates phi(t) based on the formula derived from Newton's Forward polynomial.
        phi(t) = (1/C_1) * [ (y_target - y_0) - sum_{i=2 to n}( C_i/i! * B_i(t) ) ]
        where C_i = Delta^i y_0
        """
        n = len(C_coeffs) - 1
        C_1 = C_coeffs[1]
    
        y_diff = y_target - y_tmp #y_tmp = y_0
    
        sum_higher_order = 0.0
        for i in range(2, n + 1):
            C_i = C_coeffs[i]
            D_i = C_i / math.factorial(i)
            B_i_t = evaluate_basis_product(t, i, condition)
            sum_higher_order += D_i * B_i_t

        return (y_diff - sum_higher_order) / C_1

    if condition == 0:
        """
        Calculates phi(t) based on the formula derived from Newton's Backward polynomial.
        phi(t) = (1/C_1) * [ (y_target - y_n) - sum_{i=2 to n}( C_i/i! * B_i(t) ) ]
        where C_i = Delta^i y_{n-i}
        """
        n = len(C_coeffs) - 1
        C_1 = C_coeffs[1] # This is Delta y_{n-1}
    
        y_diff = y_target - y_tmp #y_tmp = y_n
    
        sum_higher_order = 0.0
        for i in range(2, n + 1):
            C_i = C_coeffs[i] # This is Delta^i y_{n-i}
            D_i = C_i / math.factorial(i)
            B_i_t = evaluate_basis_product(t, i, condition)
            sum_higher_order += D_i * B_i_t
        
        return (y_diff - sum_higher_order) / C_1
def inverse_interpolation_newton_fixed_point(points, y_target, eps, condition):
    if condition == 1:
        """
        Finds the x value for a given y_target using the fixed-point iteration
        method derived from Newton's Forward formula.
    
        Args:
            points (list of tuples): List of (x, y) data points. 
                                 Must be evenly spaced and monotonic.
            y_target (float): The y-value to find the x for.
            eps (float): The tolerance for convergence |t_{k+1} - t_k| < eps.
        
        Returns:
            float: The interpolated x-value.
        """
    
        n = len(points) - 1
        if n < 1:
            raise ValueError("At least two points are required.")
    
        # --- 1. Setup and Difference Table ---
        x_0 = points[0][0]
        y_0 = points[0][1]
        h = points[1][0] - points[0][0]
    
        diff_table, C_coeffs = EvenDifference(points, condition=1)
    
        C_1 = C_coeffs[1]
        if np.isclose(C_1, 0):
            raise ValueError("Delta_y0 (C_1) is zero, this method cannot be used. "
                         "The function may not be monotonic.")
                         
        print("--- Inverse Interpolation Setup ---")
        print(f"Target y_target = {y_target}")
        print(f"Interval starts at x_0 = {x_0}, y_0 = {y_0}")
        print(f"h = {h:.2f}, C_1 (Delta_y0) = {C_1:.12f}")
    
        # --- 2. Initial Guess t_0 ---
        t_0 = (y_target - y_0) / C_1
        print(f"Initial guess: t_0 = (y_target - y_0) / C_1 = {t_0:.12f}")
        print("-" * 30)

        # --- 3. Perform Fixed-Point Iteration ---
        print("--- Iteration Steps ---")
    
        t_k = t_0
        results = []
    
        for k in range(100): # Max 100 iterations
            # t_{k+1} = phi(t_k)
            t_kplus1 = phi_function(t_k, y_target, y_0, C_coeffs, condition)
        
            error = abs(t_kplus1 - t_k)
        
            results.append({
                'k': k,
                't_k': t_k,
                't_k+1 = phi(t_k)': t_kplus1,
                'error = |t_k+1 - t_k|': error
            })
        
            t_k = t_kplus1
        
            if error < eps:
                break
            
        df_results = pd.DataFrame(results, columns=['k', 't_k', 't_k+1 = phi(t_k)', 'error = |t_k+1 - t_k|'])
    
        # --- 4. Final Result ---
        final_t = t_k
        final_x = x_0 + final_t * h
    
        print(f"\nConvergence reached in {k+1} iterations.")
        print(f"Final approximated x = x_0 + t*h = {x_0} + {final_t:.12f} * {h:.2f}")
    
        return diff_table, df_results, final_t, final_x
    
    if condition == 0:
        """
        Finds the x value for a given y_target using the fixed-point iteration
        method derived from Newton's Backward formula.
    
        Args:
            points (list of tuples): List of (x, y) data points. 
                                 Must be evenly spaced and monotonic.
            y_target (float): The y-value to find the x for.
            eps (float): The tolerance for convergence |t_{k+1} - t_k| < eps.
        
        Returns:
            float: The interpolated x-value.
        """
    
        n = len(points) - 1
        if n < 1:
            raise ValueError("At least two points are required.")
    
        # --- 1. Setup and Difference Table ---
        x_n = points[-1][0]
        y_n = points[-1][1]
        h = points[1][0] - points[0][0] # h is constant
    
        # Call EvenDifference with condition=0 to get BACKWARD coefficients
        diff_table, C_coeffs_bwd = EvenDifference(points, condition=0)
    
        C_1 = C_coeffs_bwd[1] # C_1 = Delta y_{n-1}
        if np.isclose(C_1, 0):
            raise ValueError("Delta_y_{n-1} (C_1) is zero, this method cannot be used. "
                         "The function may not be monotonic.")
                         
        print("--- Inverse Interpolation Setup (Newton-Backward) ---")
        print(f"Target y_target = {y_target}")
        print(f"Interval ends at x_n = {x_n}, y_n = {y_n}")
        print(f"h = {h:.2f}, C_1 (Delta y_n-1) = {C_1:.12f}")
    
        # --- 2. Initial Guess t_0 ---
        t_0 = (y_target - y_n) / C_1
        print(f"Initial guess: t_0 = (y_target - y_n) / C_1 = {t_0:.12f}")
        print("-" * 30)

        # --- 3. Perform Fixed-Point Iteration ---
        print("--- Iteration Steps ---")
    
        t_k = t_0
        results = []
    
        for k in range(100): # Max 100 iterations
            # t_{k+1} = phi(t_k)
            t_kplus1 = phi_function(t_k, y_target, y_n, C_coeffs_bwd, condition)
        
            error = abs(t_kplus1 - t_k)
        
            results.append({
                'k': k,
                't_k': t_k,
                't_k+1 = phi(t_k)': t_kplus1,
                'error = |t_k+1 - t_k|': error
            })
        
            t_k = t_kplus1
        
            if error < eps:
                break
            
        df_results = pd.DataFrame(results, columns=['k', 't_k', 't_k+1 = phi(t_k)', 'error = |t_k+1 - t_k|'])
    
        # --- 4. Final Result ---
        final_t = t_k
        final_x = x_n + final_t * h
    
        print(f"\nConvergence reached in {k+1} iterations.")
        print(f"Final approximated x = x_n + t*h = {x_n} + ({final_t:.12f}) * {h:.2f}")
    
        return diff_table, df_results, final_t, final_x
# 2. Define the problem
y_target_val = -1.43
tolerance = 1e-7

diff_table, df_results, t_solution, x_solution_fwd = inverse_interpolation_newton_fixed_point(
    points=monotonic_points_fwd,
    y_target=y_target_val,
    eps=tolerance,
    condition = 0
)
print("--- Generated Finite Difference Table ---")

diff_table.style
print("--- Iteration Steps ---")

df_results
print(f"Final approximated t = {t_solution:.12f}")
print(f"Final approximated x = {x_solution_fwd:.12f}")




