# =============================================================================
# pp2b_Trapezoid.py - Công thức hình thang (Trapezoidal Rule)
#
# Chức năng: Tính tích phân xác định bằng công thức hình thang ghép,
#            tự động chọn số bước để đạt sai số eps cho trước.
#
# Các hàm chính:
#   solve_trapezoidal(f, d2f, a, b, eps)
#     - f: hàm cần tính tích phân
#     - d2f: đạo hàm cấp 2 (dùng ước lượng sai số)
#     - a, b: cận tích phân
#     - eps: sai số cho phép
#
# Cách dùng: python pp2b_Trapezoid.py
# =============================================================================
import numpy as np
import pandas as pd
import math
def get_m2_max_abs_second_derivative(func_d2, a, b, num_samples=1000):
    """
    Estimates M2 = max|f''(x)| on [a, b] by sampling.
    """
    x_test = np.linspace(a, b, num_samples)
    y_test = np.abs(func_d2(x_test))
    return np.max(y_test)
def loopTimes (f, d2f, a, b, eps):
    # Step 1: Calculate M2
    M2 = get_m2_max_abs_second_derivative(d2f, a, b)

    # Step 2: Determine number of segments n
    # Error <= (M2 * (b-a)^3) / (12 * n^2) < eps
    # n > sqrt( (M2 * (b-a)^3) / (12 * eps) )
    numerator = M2 * (b - a)**3
    denominator = 12 * eps
    n_float = np.sqrt(numerator / denominator)
    n = math.ceil(n_float)

    return numerator, denominator, n_float, n
def solve_trapezoidal(f, d2f, a, b, eps):
    M2 = get_m2_max_abs_second_derivative(d2f, a, b)
    numerator, denominator, n_float, n = loopTimes(f, d2f, a, b, eps)
    h = (b - a) / n

    # Step 4: Generate points x_k and values y_k
    # x points
    x_values = np.linspace(a, b, n + 1)
    # y points
    y_values = f(x_values)

    # Display Table using Pandas
    data = {
        'k': range(n + 1),
        'x_k': x_values,
        'y_k': y_values
    }
    df = pd.DataFrame(data)

    # Step 5: Apply Trapezoidal Formula
    # I = h/2 * [y0 + yn + 2 * sum(y_1 to y_{n-1})]
    sum_middle = np.sum(y_values[1:-1])
    y_start = y_values[0]
    y_end = y_values[-1]
    
    I = (h / 2) * (y_start + y_end + 2 * sum_middle)
    
    return df, I
f = lambda x: 1/(x**2+1)

d2f = lambda x: (2*(3*x**2-1))/(x**6+3*x**4+3*x**2+1)

a = 0.0       # Start of interval
b = 2.0       # End of interval
eps = 1e-6 # Desired maximum error

print(f"--- Trapezoidal Rule Calculation ---")
print(f"Goal: Integrate f(x) on [{a}, {b}] with error < {eps}\n")
M2 = get_m2_max_abs_second_derivative(d2f, a, b)
print(f"Step 1: Estimate max|f''(x)| (M2)")
print(f"M2 ≈ {M2:.6f}")
numerator, denominator, n_float, n = loopTimes(f, d2f, a, b, eps)
print(f"\nStep 2: Calculate minimum n")
print(f"n > sqrt({numerator:.6f} / {denominator:.6f}) = {n_float:.6f}")
print(f"Chosen n = {n}")
h = (b - a) / n
print(f"\nStep 3: Calculate step size h")
print(f"h = ({b} - {a}) / {n} = {h:.6f}")
df, I = solve_trapezoidal(f, d2f, a, b, eps)
print(f"\nStep 4: Table of values (x_k, y_k)")
df

print(f"\nStep 5: Final Calculation")
print(f"Approximate Integral I = {I:.8f}")




