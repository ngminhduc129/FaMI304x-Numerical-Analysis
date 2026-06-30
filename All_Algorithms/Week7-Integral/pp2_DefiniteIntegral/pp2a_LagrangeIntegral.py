# =============================================================================
# pp2a_LagrangeIntegral.py - Tích phân bằng Lagrange
#
# Chức năng: Tính tích phân xác định bằng cách nội suy Lagrange
#            hàm dưới dấu tích phân rồi tích phân đa thức nội suy.
#
# Các hàm chính:
#   solve_lagrange_integral(points, a, b)
#     - points: list các cặp (x, f(x))
#     - a, b: cận tích phân
#
# Cách dùng: python pp2a_LagrangeIntegral.py
# =============================================================================
import numpy as np
import pandas as pd

# Set display options
pd.set_option('display.precision', 6)
pd.set_option('display.width', 1000)
def polynomial_integrate_definite(coeffs, a, b):
    """
    Integrates a polynomial P(x) defined by coeffs [c_n, ..., c_0] 
    from a to b.
    P(x) = c_n*x^n + ... + c_1*x + c_0
    """
    # Integration increases power by 1 and divides by new power
    # Coeffs: [c_n, c_{n-1}, ..., c_0]
    # Int Coeffs: [c_n/(n+1), ..., c_0/1, 0]
    
    n = len(coeffs) - 1
    int_coeffs = []
    
    for i, c in enumerate(coeffs):
        power = n - i
        int_coeffs.append(c / (power + 1))
    int_coeffs.append(0) # Constant C becomes 0 for definite integral
    
    # Define evaluation function
    def eval_poly(c_list, x_val):
        val = 0.0
        degree = len(c_list) - 1
        for i, c in enumerate(c_list):
            val += c * (x_val ** (degree - i))
        return val
    
    val_b = eval_poly(int_coeffs, b)
    val_a = eval_poly(int_coeffs, a)
    
    return val_b - val_a
def solve_lagrange_integral(points, a, b):
    """
    Calculates definite integral using General Lagrange Interpolation.
    
    Args:
        points (list): List of (x, y) tuples.
        a (float): Lower limit.
        b (float): Upper limit.
        
    Returns:
        pd.DataFrame: Table of weights.
        float: Integral result.
    """
    n = len(points)
    x_vals = [p[0] for p in points]
    y_vals = [p[1] for p in points]
    
    weights = []
    
    # 1. Calculate Weight w_i for each point
    for i in range(n):
        xi = x_vals[i]
        
        # Construct Numerator Polynomial: Product of (x - xj)
        # Represent polynomial as list of coefficients [highest_power, ..., constant]
        # (x - xj) is represented as [1, -xj]
        poly_numerator = [1.0] 
        denominator = 1.0
        
        for j in range(n):
            if i == j:
                continue
                
            xj = x_vals[j]
            
            # Update Denominator
            denominator *= (xi - xj)
            
            # Update Numerator Polynomial (Multiplication)
            # np.convolve computes the product of two polynomials
            poly_numerator = np.convolve(poly_numerator, [1.0, -xj])
            
        # 2. Integrate Numerator from a to b
        integral_val = polynomial_integrate_definite(poly_numerator, a, b)
        
        # 3. Final Weight
        w_i = integral_val / denominator
        weights.append(w_i)
        
    # 4. Compute Final Sum
    weighted_y = np.array(y_vals) * np.array(weights)
    I = np.sum(weighted_y)
    
    # Output Table
    df = pd.DataFrame({
        'x_i': x_vals,
        'y_i': y_vals,
        'Weight w_i': weights,
        'Term (y_i * w_i)': weighted_y
    })
    
    return df, I


# --- EXAMPLE USAGE ---
# Using the data from Page 3 of your PDF as a test case
# Let's integrate from x=1.0 to x=1.6 using the first 4 points (uneven check)
# or just the first 3 points to check against Simpson's rule logic.


# Case 1: 3 equidistant points (Should behave like Simpson's rule)
# Interval [1.0, 1.2], h=0.1. Simpson weights: h/3 * [1, 4, 1] = 0.1/3 * [1, 4, 1] 
# Expected Weights approx: [0.0333, 0.1333, 0.0333]
points_simpson = [
    (1.0, -0.641),
    (1.1, -0.498),
    (1.2, -0.340)
]
a1, b1 = 1.0, 1.2

print(f"\nTest 1: Integration on [1.0, 1.2] (Equidistant 3 points)")
df1, I1 = solve_lagrange_integral(points_simpson, a1, b1)

df1
print(f"Approximate Integral I = {I1:.8f}")

# Case 2: Unevenly spaced points (General Case)
# Points: 1.0, 1.1, 1.3 (Skip 1.2)
points_uneven = [
    (1.0, -0.641),
    (1.1, -0.498),
    (1.3, -0.165)
]
a2, b2 = 1.0, 1.3


print(f"\nTest 2: Integration on [1.0, 1.3] (Uneven points: 1.0, 1.1, 1.3)")
df2, I2 = solve_lagrange_integral(points_uneven, a2, b2)
df2

print(f"Approximate Integral I = {I2:.8f}")




