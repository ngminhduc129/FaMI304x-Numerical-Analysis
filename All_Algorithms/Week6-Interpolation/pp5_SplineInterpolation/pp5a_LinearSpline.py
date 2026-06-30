# =============================================================================
# pp5a_LinearSpline.py - Nội suy Spline tuyến tính
#
# Chức năng: Nội suy bằng các đoạn thẳng nối giữa các điểm dữ liệu
#            (spline bậc 1), đảm bảo liên tục C0.
#
# Các hàm chính:
#   linear_spline_interpolation(points)     - tạo các đoạn spline
#   evaluate_linear_spline(segments, x_val) - tính giá trị tại x
#
# Cách dùng: python pp5a_LinearSpline.py
# =============================================================================
import numpy as np
import pandas as pd
import math
from typing import Tuple

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def linear_spline_interpolation(points):
    """
    Calculates the set of linear spline functions S_i(x) for each interval
    in standard polynomial form: S_i(x) = a_1*x + a_0.
    
    Args:
        points (list of tuples): List of (x, y) data points, sorted by x.

    Returns:
        pd.DataFrame: A DataFrame containing the coefficients for each
                      linear segment S_i(x).
    """
    n = len(points) - 1
    segments_data = []

    for i in range(n):
        x_i, y_i = points[i]
        x_ip1, y_ip1 = points[i+1]
        
        h = x_ip1 - x_i
        if h <= 0:
            raise ValueError(f"Points are not sorted or x values are not distinct at index {i}.")
            
        # S_i(x) = a_1 * x + a_0
        # a_1 = m_i = (y_{i+1} - y_i) / (x_{i+1} - x_i)
        a_1 = (y_ip1 - y_i) / h
        
        # a_0 = y_i - m_i * x_i
        a_0 = y_i - a_1 * x_i
        
        segments_data.append({
            'i': i,
            'Interval': f"[{x_i}, {x_ip1}]",
            'a_0 (Degree 0)': a_0,
            'a_1 (Degree 1)': a_1,
        })
        
    segments_df = pd.DataFrame(segments_data, columns=['i', 'Interval', 'a_0 (Degree 0)', 'a_1 (Degree 1)'])
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
# 2. Generate the spline segments
segments_table = linear_spline_interpolation(points)

print("--- Linear Spline Segments ---")
segments_table.style
def evaluate_linear_spline(segments_df, x_val):
    """
    Evaluates the piecewise linear spline at a specific x value.
    This function is now designed to work with the table of a_0 and a_1 coefficients.
    
    Args:
        segments_df (pd.DataFrame): The DataFrame of spline segments from
                                    linear_spline_interpolation.
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
        
        if x_i <= x_val <= x_ip1:
            segment_found = row
            break
            
    # Handle edge case: x_val is slightly larger than the last node
    # (e.g., due to floating point), or exactly on the last node.
    if segment_found is None and np.isclose(x_val, last_x_ip1):
         segment_found = segments_df.iloc[-1]
    
    if segment_found is None:
        # Check if x_val is outside the total range
        first_x_i = float(segments_df.iloc[0]['Interval'].strip('[]').split(', ')[0])
        if x_val < first_x_i or x_val > last_x_ip1:
            raise ValueError(f"x_val = {x_val} is outside the interpolation range [{first_x_i}, {last_x_ip1}]")
        # This case should not be hit if points are contiguous
        raise ValueError(f"Could not find a segment for x_val = {x_val}")

    
    # 2. Get coefficients from the found segment
    i = segment_found['i']
    a_0 = segment_found['a_0 (Degree 0)']
    a_1 = segment_found['a_1 (Degree 1)']
    
    # 3. Apply the formula S_i(x) = a_1 * x + a_0
    y_val = a_1 * x_val + a_0
    
    # --- Print Intermediate Steps ---
    print("\n" + "-"*30)
    print("--- Linear Spline Evaluation ---")
    print(f"Target value x = {x_val}")
    print(f"Value is in {segment_found['Interval']} (segment i={i})")
    print(f"Using formula: S_{i}(x) = a_1 * x + a_0")
    print(f"S_{i}({x_val}) = {a_1:.12f} * {x_val} + {a_0:.12f}")
    print(f"S_{i}({x_val}) = {y_val:.12f}")
    
    return y_val
# 3. Evaluate the spline at a specific point (e.g., x = 1.75)
x_to_find = 1.75
y_found = evaluate_linear_spline(segments_table, x_to_find)




