# =============================================================================
# pp6_LSQ.py - Phương pháp bình phương tối thiểu (Least Squares)
#
# Chức năng: Tìm hàm xấp xỉ tốt nhất cho tập dữ liệu theo nguyên lý
#            bình phương tối thiểu. Hỗ trợ 3 mô hình: tuyến tính,
#            exponential (mũ), power law (lũy thừa).
#
# Các hàm chính:
#   fit_linear_model(points)         - y = a0 + a1*x
#   fit_exponential_model(points)    - y = a*e^(b*x)
#   fit_power_law_model(points)      - y = a*x^b
#   least_squares_fit(x, y, basis)   - fit với hàm cơ sở tùy chỉnh
#
# Cách dùng: python pp6_LSQ.py
# =============================================================================
import numpy as np
import pandas as pd
import os
from typing import List, Tuple, Callable

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def calculate_least_squares_error(
    y_vals: np.ndarray, 
    b_prime: np.ndarray, 
    a_coeffs: np.ndarray
) -> Tuple[float, float, float]:
    """
    Calculates the SSE, MSE, and RMSE of the fit.
    
    Args:
        y_vals (np.ndarray): The n x 1 target vector y.
        b_prime (np.ndarray): The m x 1 vector (Phi^T * y).
        a_coeffs (np.ndarray): The m x 1 solution vector a.
        
    Returns:
        (SSE, MSE, RMSE)
    """
    n = len(y_vals)
    
    # 1. [y, y] = y^T * y
    y_t_y = y_vals.T.dot(y_vals).item()
    
    # 2. sum(a_j * [y, phi_j]) = a^T * b'
    a_t_b_prime = a_coeffs.T.dot(b_prime).item()
    
    # 3. SSE = [y, y] - sum(a_j * [y, phi_j])
    sse = y_t_y - a_t_b_prime
    
    # 4. MSE = SSE / n
    mse = sse / n
    
    # 5. RMSE = sqrt(MSE)
    rmse = np.sqrt(mse)
    
    return sse, mse, rmse
def least_squares_fit(
    points: List[Tuple[float, float]], 
    basis_funcs: List[Callable[[float], float]]
) -> np.ndarray:
    """
    Finds the coefficients for a set of basis functions that best fit
    the given (x, y) data points using the method of Normal Equations.

    Args:
        points: A list of (x, y) data tuples.
        basis_funcs: A list of lambda functions representing the basis [g1, g2, ..., gm].

    Returns:
        np.ndarray: The column vector of optimal coefficients [a1, a2, ..., am].
    """
    
    # --- 1. Construct Input Matrices and Vectors ---
    n = len(points)
    m = len(basis_funcs)
    
    if n < m:
        print("Warning: More basis functions (m) than data points (n). "
              "The system is underdetermined.")

    x_vals = np.array([p[0] for p in points])
    y_vals = np.array([p[1] for p in points]).reshape(-1, 1) # y must be n x 1

    df_input = pd.DataFrame(points, columns=['x_i', 'y_i'])

    # --- 2. Construct Design Matrix Phi (n x m) ---
    # Phi_ij = phi_j(x_i)
    phi_matrix = np.zeros((n, m))
    for j, func in enumerate(basis_funcs):
        phi_matrix[:, j] = func(x_vals)
    
    df_phi = pd.DataFrame(phi_matrix, columns=[f"\u03C6_{j+1}(x)" for j in range(m)])

    # --- 3. Form Gram Matrix M = Phi^T * Phi (m x m) ---
    phi_t = phi_matrix.T
    M = phi_t.dot(phi_matrix)
    
    df_M = pd.DataFrame(M)

    # --- 4. Form Target Vector b' = Phi^T * y (m x 1) ---
    b_prime = phi_t.dot(y_vals)
    
    df_b_prime = pd.DataFrame(b_prime)

    # --- 5. Solve the Normal Equations Ma = b' ---
    # Use np.linalg.solve, which is a robust and fast solver (like Gauss elimination)
    # as per your instruction to use core numpy/pandas.
    try:
        a_coeffs = np.linalg.solve(M, b_prime)
    except np.linalg.LinAlgError:
        print("Error: The matrix M is singular. "
              "The basis functions may be linearly dependent.")
        return None

    df_a = pd.DataFrame(a_coeffs, columns=['a_j'])
    df_a.index = [f"a_{j+1}" for j in range(m)]

    # --- 6. Calculate Error Metrics ---
    sse, mse, rmse = calculate_least_squares_error(y_vals, b_prime, a_coeffs)

    return df_phi, df_M, df_b_prime, df_a, a_coeffs, sse, mse, rmse
# --- Example Usage ---

points_data = [
    (np.deg2rad(30), 2.611),
    (np.deg2rad(60), 3.102),
    (np.deg2rad(90), 2.912),
    (np.deg2rad(120), 2.105),
    (np.deg2rad(150), 0.612),
    (np.deg2rad(180), -1.321),
    (np.deg2rad(210), -1.906),
    (np.deg2rad(240), -2.412),
    (np.deg2rad(270), -2.802),
    (np.deg2rad(300), -2.703),
    (np.deg2rad(330), -1.610),
    (np.deg2rad(360), -1.500),
]

basis_set = [
    lambda x: 1,      
    lambda x: np.cos(x),      
    lambda x: np.cos(2*x), 
    lambda x: np.sin(x), 
    lambda x: np.sin(2*x)
]

df_phi, df_M, df_b_prime, df_a, a_coeffs, sse, mse, rmse = least_squares_fit(points_data, basis_set)
# Construct Design Matrix Phi (n x m)
df_phi.style
## Form Gram Matrix M = Phi^T * Phi (m x m)
df_M.style
## Form Target Vector b' = Phi^T * y (m x 1)
df_b_prime.style
## Solve the Normal Equations Ma = b'
df_a.style
## Error
print(f"Sum of Squared Errors (S_min): {sse: .12f}")
print(f"Mean Squared Error (MSE):     {mse: .12f}")
print(f"Root Mean Square Error (RMSE):  {rmse: .12f}")
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
file_to_read = '0.csv' 
    
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
def modify_value(original_points: List[Tuple[float, float]],
                 transform: Callable[[float, float], Tuple[float, float]] = lambda x, y: (x, y)
):
    """Apply a user-provided transform to each (x, y) and return the transformed list.

    Args:
        original_points: list of (x, y) tuples.
        transform: callable that takes (x, y) and returns (X, Y).
                   Example: lambda x, y: (x, np.log(y))

    Returns:
        List of transformed (X, Y) tuples.
    """
    transformed_points: List[Tuple[float, float]] = []
    for x, y in original_points:
        try:
            X, Y = transform(x, y)
        except Exception as e:
            raise ValueError(f"Error applying transform to point (x={x}, y={y}): {e}")
        transformed_points.append((X, Y))

    return transformed_points
def fit_linear_model(
    transformed_points: List[Tuple[float, float]]
) -> Tuple[float, float]:

    # --- 2. Define the Linear Problem ---
    # We are fitting Y = A*X + B*SINX
    # Basis function for A is phi_1(X) = X
    # Basis function for B is phi_2(X) = SINX
    linear_basis = [
        lambda x: x,  # For A
        lambda x: np.sin(x)   # For B
    ]
    
    # --- 3. Solve for Linear Coefficients (A, B) ---
    # linear_coeffs will be [ [A], [B] ]
    df_phi, df_M, df_b_prime, df_a, linear_coeffs, sse, mse, rmse = least_squares_fit(
        transformed_points, 
        linear_basis
    )
    
    if linear_coeffs is None:
        return None, None
        
    # --- 4. Transform Coefficients Back ---
    A = linear_coeffs[0][0]  # A = a
    B = linear_coeffs[1][0]  # B = b
    
    a = A
    b = B
    
    return A, B, a, b
print("Original (x, y) converted to (X, Y) where X=x, Y=y")
transformed_points = modify_value(all_points, lambda x, y: (x,y))

df_transformed = pd.DataFrame(transformed_points, columns=['X_i', 'Y_i'])
# --- Run the Full Process ---
A, B, a,b = fit_linear_model(transformed_points)
print(f"Solved Linear A (a) = {A: .12f}")
print(f"Solved Linear B (b)     = {B: .12f}")
print(f"-> Non-Linear a = A = {a: .12f}")
print(f"-> Non-Linear b = B   = {b: .12f}")
print("-" * 30)
print(f"Final Equation: y = {a:.6f} * x  + {b:.6f} * sinx)")
def fit_exponential_model(
    transformed_points: List[Tuple[float, float]]
) -> Tuple[float, float]:
    """
    Fits data to y = a * exp(b*x) by linearizing to ln(y) = ln(a) + b*x.
    
    Args:
        original_points: List of (x, y) data tuples.
        
    Returns:
        (a, b): The calculated non-linear coefficients.
    """

    # --- 2. Define the Linear Problem ---
    # We are fitting Y = A + B*X
    # Basis function for A is phi_1(X) = 1
    # Basis function for B is phi_2(X) = X
    linear_basis = [
        lambda x: 1,  # For A
        lambda x: x   # For B
    ]
    
    # --- 3. Solve for Linear Coefficients (A, B) ---
    # linear_coeffs will be [ [A], [B] ]
    df_phi, df_M, df_b_prime, df_a, linear_coeffs, sse, mse, rmse = least_squares_fit(
        transformed_points, 
        linear_basis
    )
    
    if linear_coeffs is None:
        return None, None
        
    # --- 4. Transform Coefficients Back ---
    A = linear_coeffs[0][0]  # A = ln(a)
    B = linear_coeffs[1][0]  # B = b
    
    a = np.exp(A)
    b = B
    
    return A, B, a, b
#--- Example Data ---
example_points = [
    (1, 7.1),
    (2, 27.8),
    (3, 63.1),
    (4, 110),
    (5, 161)
]
print("Original (x, y) converted to (X, Y) where X=x, Y=ln(y)")
transformed_points = modify_value(example_points, lambda x, y: (x, np.log(y)))

df_transformed = pd.DataFrame(transformed_points, columns=['X_i', 'Y_i'])

df_transformed.style
# --- Run the Full Process ---
A, B, a,b = fit_exponential_model(transformed_points)
print(f"Solved Linear A (ln(a)) = {A: .12f}")
print(f"Solved Linear B (b)     = {B: .12f}")
print(f"-> Non-Linear a = e^A = {a: .12f}")
print(f"-> Non-Linear b = B   = {b: .12f}")
print("-" * 30)
print(f"Final Equation: y = {a:.6f} * e^({b:.6f} * x)")
def fit_power_law_model(
    transformed_points: List[Tuple[float, float]]
) -> Tuple[float, float]:
    """
    Fits data to y = a * x^b by linearizing to ln(y) = ln(a) + b*ln(x).
    
    Args:
        original_points: List of (x, y) data tuples.
        
    Returns:
        (a, b): The calculated non-linear coefficients.
    """
        
    # --- 2. Define the Linear Problem ---
    # We are fitting Y = A + B*X
    # Basis function for A is phi_1(X) = 1
    # Basis function for B is phi_2(X) = X
    linear_basis = [
        lambda x: 1,  # For A
        lambda x: x   # For B
    ]
    
    # --- 3. Solve for Linear Coefficients (A, B) ---
    # This function will print steps 1-6 for the *linear* problem
    df_phi, df_M, df_b_prime, df_a, linear_coeffs, sse, mse, rmse = least_squares_fit(
        transformed_points, 
        linear_basis,
    )
    
    if linear_coeffs is None:
        return None, None
        
    # --- 4. Transform Coefficients Back ---
    A = linear_coeffs[0][0]  # A = ln(a)
    B = linear_coeffs[1][0]  # B = b
    
    a = np.exp(A)
    b = B
    
    
    return A, B, a, b
#--- Example Data ---
example_points = [
    (1, 7.1),
    (2, 27.8),
    (3, 63.1),
    (4, 110),
    (5, 161)
]

print("Original (x, y) converted to (X, Y) where X=ln(x), Y=ln(y)")
transformed_points = modify_value(example_points, lambda x, y: (np.log(x), np.log(y)))

df_transformed = pd.DataFrame(transformed_points, columns=['X_i', 'Y_i'])

df_transformed.style
# --- Run the Full Process ---
A, B, a, b = fit_power_law_model(transformed_points)
print(f"Solved Linear A (ln(a)) = {A: .12f}")
print(f"Solved Linear B (b)     = {B: .12f}")
print(f"-> Non-Linear a = e^A = {a: .12f}")
print(f"-> Non-Linear b = B   = {b: .12f}")
print("-" * 30)
print(f"Final Equation: y = {a:.6f} * x^({b:.6f})")




