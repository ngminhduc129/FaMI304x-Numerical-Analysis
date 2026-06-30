# =============================================================================
# pp2c_Simpson.py - Công thức Simpson 1/3
#
# Chức năng: Tính tích phân xác định bằng công thức Simpson ghép,
#            độ chính xác cao hơn hình thang (bậc 3).
#
# Các hàm chính:
#   solve_simpson(f, d4f, a, b, eps)
#     - f: hàm cần tính tích phân
#     - d4f: đạo hàm cấp 4 (dùng ước lượng sai số)
#     - a, b: cận tích phân
#     - eps: sai số cho phép
#
# Cách dùng: python pp2c_Simpson.py
# =============================================================================
import numpy as np
import pandas as pd
import math

# Set display options
pd.set_option('display.precision', 6)
pd.set_option('display.width', 1000)
# Cell 2: Helper to estimate Max 4th Derivative (M4)
def get_m4_max_abs_fourth_derivative(func_d4, a, b, num_samples=1000):
    """
    Estimates M4 = max|f''''(x)| on [a, b] by sampling.
    Needed for Simpson's Rule error estimation.
    """
    x_test = np.linspace(a, b, num_samples)
    y_test = np.abs(func_d4(x_test))
    return np.max(y_test)                                      
# Cell 4: Calculation of n (loopTimes)
def loopTimes_simpson(f, d4f, a, b, eps):
    # Step 1: Calculate M4
    M4 = get_m4_max_abs_fourth_derivative(d4f, a, b)

    # Step 2: Determine number of segments n
    # Error <= (M4 * (b-a)^5) / (180 * n^4) < eps
    # n > [ (M4 * (b-a)^5) / (180 * eps) ] ^ (1/4)
    numerator = M4 * (b - a)**5
    denominator = 180 * eps
    n_float = (numerator / denominator)**0.25
    
    n = math.ceil(n_float)
    
    # Simpson's rule requires n to be EVEN
    if n % 2 != 0:
        n += 1

    return numerator, denominator, n_float, n
# Cell 5: Solve Function
def solve_simpson(f, d4f, a, b, eps):
    M4 = get_m4_max_abs_fourth_derivative(d4f, a, b)
    numerator, denominator, n_float, n = loopTimes_simpson(f, d4f, a, b, eps)
    h = (b - a) / n

    # Step 4: Generate points x_k and values y_k
    x_values = np.linspace(a, b, n + 1)
    y_values = f(x_values)

    # Display Table using Pandas
    data = {
        'k': range(n + 1),
        'x_k': x_values,
        'y_k': y_values
    }
    df = pd.DataFrame(data)

    # Step 5: Apply Simpson's Formula
    # I = h/3 * [y0 + yn + 4*Odd_Sum + 2*Even_Sum]
    
    y_start = y_values[0]
    y_end = y_values[-1]
    
    # Odd indices: 1, 3, 5, ..., n-1
    sum_odd = np.sum(y_values[1:n:2])
    
    # Even indices: 2, 4, 6, ..., n-2 (exclude 0 and n)
    sum_even = np.sum(y_values[2:n-1:2])
    
    I = (h / 3) * (y_start + y_end + 4 * sum_odd + 2 * sum_even)
    
    return df, I


# Cell 7: Define Problem (Using the same example: 1/(x^2+1))
f = lambda x: 1/(x**2+1)

# 4th derivative of 1/(x^2+1)
# Derived as: 24(5x^4 - 10x^2 + 1) / (x^2+1)^5
d4f = lambda x: (24 * (5*x**4 - 10*x**2 + 1)) / ((x**2 + 1)**5)

a = 0.0       # Start of interval
b = 2.0       # End of interval
eps = 1e-6    # Desired maximum error

print(f"--- Simpson's 1/3 Rule Calculation ---")
print(f"Goal: Integrate f(x) on [{a}, {b}] with error < {eps}\n")
# Cell 8: Step 1 Output
M4 = get_m4_max_abs_fourth_derivative(d4f, a, b)
print(f"Step 1: Estimate max|f''''(x)| (M4)")
print(f"M4 ≈ {M4:.6f}")
# Cell 9: Step 2 Output
numerator, denominator, n_float, n = loopTimes_simpson(f, d4f, a, b, eps)
print(f"\nStep 2: Calculate minimum n (must be even)")
print(f"n > ({numerator:.6f} / {denominator:.6f})^(1/4) = {n_float:.6f}")
print(f"Chosen n = {n}")

# Cell 10: Step 3 Output
h = (b - a) / n
print(f"\nStep 3: Calculate step size h")
print(f"h = ({b} - {a}) / {n} = {h:.6f}")
# Cell 11: Step 4 Output (Table)
df, I = solve_simpson(f, d4f, a, b, eps)
print(f"\nStep 4: Table of values (x_k, y_k)")
df
# Cell 12: Step 5 Output (Final)
print(f"\nStep 5: Final Calculation")
print(f"Approximate Integral I = {I:.8f}")




