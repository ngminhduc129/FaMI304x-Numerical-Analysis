# =============================================================================
# pp5b_QuadraticSpline.py - Nội suy Spline bậc 2
#
# Chức năng: Nội suy bằng các đoạn parabol nối trơn giữa các điểm,
#            đảm bảo liên tục đến đạo hàm cấp 1 (C1).
#
# Các hàm chính:
#   quadratic_spline_interpolation(points, z0_condition) - tạo spline
#   evaluate_quadratic_spline(segments, x_val)            - tính giá trị
#
# Cách dùng: python pp5b_QuadraticSpline.py
# =============================================================================
import numpy as np
import pandas as pd
import math
from typing import Tuple

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def quadratic_spline_interpolation(points, z0_condition):
    """
    Calculates the set of quadratic spline functions S_i(x) for each interval
    in standard polynomial form: S_i(x) = a_2*x^2 + a_1*x + a_0.
    
    A boundary condition z_0 = S'(x_0) must be provided.
    
    Args:
        points (list of tuples): List of (x, y) data points, sorted by x.
        z0_condition (float): The boundary condition for the slope at x_0,
                              (e.g., z_0 = S'(x_0) = 0).

    Returns:
        pd.DataFrame: A DataFrame containing the coefficients for each
                      quadratic segment S_i(x).
    """
    n = len(points) - 1
    if n < 1:
        raise ValueError("At least two points are required for interpolation.")
        
    x = np.array([p[0] for p in points])
    y = np.array([p[1] for p in points])
    h = np.zeros(n)
    z = np.zeros(n + 1)
    
    z_steps_data = []

    # --- 1. Solve for all z_i (slopes at knots) ---
    z[0] = z0_condition
    z_steps_data.append({
        'i': 0,
        'x_i': x[0],
        'y_i': y[0],
        'h_i': 'N/A',
        'z_i = S\'(x_i)': z[0],
        'Calculation': 'Boundary Condition'
    })

    for i in range(n):
        h[i] = x[i+1] - x[i]
        if h[i] <= 0:
            raise ValueError(f"Points are not sorted or x values are not distinct at index {i}.")
            
        # z_{i+1} = 2*(y_{i+1} - y_i)/h_i - z_i
        z[i+1] = (2.0 * (y[i+1] - y[i]) / h[i]) - z[i]
        
        z_steps_data.append({
            'i': i + 1,
            'x_i': x[i+1],
            'y_i': y[i+1],
            'h_i': h[i],
            'z_i = S\'(x_i)': z[i+1],
            'Calculation': f"2*({y[i+1]:.4f}-{y[i]:.4f})/{h[i]:.2f} - {z[i]:.4f}"
        })

    z_steps_df = pd.DataFrame(z_steps_data)
    
    # --- 2. Calculate coefficients for each segment ---
    segments_data = []
    for i in range(n):
        # a_i = (z_{i+1} - z_i) / (2 * h_i)
        a_2 = (z[i+1] - z[i]) / (2.0 * h[i])
        
        # b_i = z_i - 2 * a_i * x_i
        a_1 = z[i] - 2.0 * a_2 * x[i]
        
        # c_i = y_i - a_i * x_i^2 - b_i * x_i
        a_0 = y[i] - a_2 * x[i]**2 - a_1 * x[i]
        
        segments_data.append({
            'i': i,
            'Interval': f"[{x[i]}, {x[i+1]}]",
            'a_0 (Degree 0)': a_0,
            'a_1 (Degree 1)': a_1,
            'a_2 (Degree 2)': a_2
        })
        
    segments_df = pd.DataFrame(segments_data, columns=['i', 'Interval', 'a_0 (Degree 0)', 'a_1 (Degree 1)', 'a_2 (Degree 2)'])
    
    return segments_df


# 1. Define the data points from the image
points = [
    (1.4, 2.4347052),
    (1.5, 2.5473627),
    (1.6, 2.6495552),
    (1.7, 2.7412610),
    (1.8, 2.8225627),
    (1.9, 2.8936474),
    (2.0, 2.9548000),
    (2.1, 3.0064203),
    (2.2, 3.0489801)
]

# 2. Define the boundary condition (z_0 = S'(x_0))
z0_boundary = 0.0
# 2. Generate the spline segments
segments_table = quadratic_spline_interpolation(points, z0_boundary)

print("--- Linear Spline Segments ---")
segments_table.style
def evaluate_quadratic_spline(segments_df, x_val):
    """
    Evaluates the piecewise quadratic spline at a specific x value.
    
    Args:
        segments_df (pd.DataFrame): The DataFrame of spline segments.
        x_val (float): The x-value to interpolate at.

    Returns:
        float: The interpolated y-value.
    """
    
    # 1. Find the correct segment
    segment_found = None
    last_x_ip1 = -np.inf
    
    for _, row in segments_df.iterrows():
        interval_str = row['Interval']
        # Parse interval string "[x_i, x_ip1]"
        x_i, x_ip1 = [float(x) for x in interval_str.strip('[]').split(', ')]
        last_x_ip1 = x_ip1
        
        # Use a small tolerance for the upper bound
        if x_i <= x_val <= (x_ip1 + 1e-14):
            segment_found = row
            break
            
    if segment_found is None:
        first_x_i = float(segments_df.iloc[0]['Interval'].strip('[]').split(', ')[0])
        if x_val < first_x_i or x_val > last_x_ip1:
            raise ValueError(f"x_val = {x_val} is outside the interpolation range [{first_x_i}, {last_x_ip1}]")
        raise ValueError(f"Could not find a segment for x_val = {x_val}")

    
    # 2. Get coefficients from the found segment
    i = segment_found['i']
    a_0 = segment_found['a_0 (Degree 0)']
    a_1 = segment_found['a_1 (Degree 1)']
    a_2 = segment_found['a_2 (Degree 2)']
    
    # 3. Apply the formula S_i(x) = a_2 * x^2 + a_1 * x + a_0
    y_val = a_2 * (x_val**2) + a_1 * x_val + a_0
    
    # --- Print Intermediate Steps ---
    print("\n" + "-"*30)
    print("--- Quadratic Spline Evaluation ---")
    print(f"Target value x = {x_val}")
    print(f"Value is in {segment_found['Interval']} (segment i={i})")
    print(f"Using formula: S_{i}(x) = a_2 * x^2 + a_1 * x + a_0")
    print(f"S_{i}({x_val}) = {a_2:.12f} * {x_val}^2 + {a_1:.12f} * {x_val} + {a_0:.12f}")
    print(f"S_{i}({x_val}) = {y_val:.12f}")
    
    return y_val
# 3. Evaluate the spline at a specific point (e.g., x = 1.75)
x_to_find = 1.75
y_found = evaluate_quadratic_spline(segments_table, x_to_find)




