# =============================================================================
# pp1_LagrangeDerivative.py - Đạo hàm bằng Lagrange
#
# Chức năng: Tính đạo hàm của hàm số tại một điểm bằng phương pháp
#            trọng số Lagrange.
#
# Các hàm chính:
#   solve_derivative_fixed(points, target_c, num_points=3)
#     - points: list các cặp (x, y)
#     - target_c: điểm cần tính đạo hàm
#     - num_points: số điểm dùng để nội suy
#
# Cách dùng: python pp1_LagrangeDerivative.py
# =============================================================================
import numpy as np
import pandas as pd
import math

# Increase precision for display
pd.set_option('display.precision', 8)

def calculate_lagrange_weights_v2(x_vals, c):
    """
    Computes the derivative of Lagrange basis polynomials L'_j(c).
    FIXED: Corrected sign error in node-aligned case.
    """
    n = len(x_vals)
    weights = np.zeros(n)
    
    # Check if c is exactly one of the nodes
    node_index = -1
    for i, x in enumerate(x_vals):
        if math.isclose(x, c, rel_tol=1e-9):
            node_index = i
            break
            
    if node_index != -1:
        # Case 1: c is a node x_k
        k = node_index
        for i in range(n):
            if i == k:
                # Diagonal element: Sum of reciprocals 1/(xk - xj)
                for j in range(n):
                    if j != k:
                        weights[i] += 1.0 / (x_vals[k] - x_vals[j])
            else:
                # Off-diagonal elements: L'_i(xk)
                # Formula: (1 / (xi - xk)) * Product(...)
                product_term = 1.0
                for j in range(n):
                    if j != k and j != i:
                        product_term *= (x_vals[k] - x_vals[j]) / (x_vals[i] - x_vals[j])
                
                # FIXED: Pre-factor is 1/(xi - xk), NOT 1/(xk - xi)
                weights[i] = (1.0 / (x_vals[i] - x_vals[k])) * product_term
    else:
        # Case 2: c is not a node (General case)
        # This part was already correct, but we include it for completeness
        for i in range(n):
            Li_c = 1.0
            for j in range(n):
                if i != j:
                    Li_c *= (c - x_vals[j]) / (x_vals[i] - x_vals[j])
            
            sum_part = 0.0
            for j in range(n):
                if i != j:
                    sum_part += 1.0 / (c - x_vals[j])
            
            weights[i] = Li_c * sum_part
            
    return weights
def solve_derivative_fixed(points, target_c, num_points=3):
    # 1. Select Neighbors
    sorted_by_dist = sorted(points, key=lambda p: abs(p[0] - target_c))
    selected_points = sorted(sorted_by_dist[:num_points], key=lambda p: p[0])
    
    x_subset = np.array([p[0] for p in selected_points])
    y_subset = np.array([p[1] for p in selected_points])
    
    # 2. Calculate Weights (Using Fixed Function)
    weights = calculate_lagrange_weights_v2(x_subset, target_c)
    
    # 3. Compute Derivative
    terms = y_subset * weights
    derivative = np.sum(terms)
    
    # 4. Create DataFrame
    df = pd.DataFrame({
        'x_k': x_subset,
        'y_k': y_subset,
        "Weight L'_k(c)": weights,
        'Term (y_k * W)': terms
    })
    
    return df, derivative
# --- INPUT DATA FROM PDF PAGE 3 ---
points = [
    (1.0, -0.641), (1.1, -0.498), (1.2, -0.340), (1.3, -0.165), 
    (1.4, 0.028), (1.5, 0.241), (1.6, 0.477), (1.7, 0.737), 
    (1.8, 1.025), (1.9, 1.343), (2.0, 1.695)
]

# --- EXECUTION FOR PDF EXAMPLES ---
target = 1.5


df_result, dy_dx = solve_derivative_fixed(points, target, num_points=3)
df_result
print(f"Approximate Derivative f'({target}) = {dy_dx:.6f}")




