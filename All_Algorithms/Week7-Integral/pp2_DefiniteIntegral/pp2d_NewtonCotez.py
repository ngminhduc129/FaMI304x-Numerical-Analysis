# =============================================================================
# pp2d_NewtonCotez.py - Công thức Newton-Cotes tổng quát
#
# Chức năng: Tính tích phân bằng công thức Newton-Cotes bậc n.
#            n=1: hình thang, n=2: Simpson, n=4: Boole.
#
# Các hàm chính:
#   solve_newton_cotes(f, a, b, n)
#     - f: hàm cần tính tích phân
#     - a, b: cận
#     - n: bậc Newton-Cotes
#
# Cách dùng: python pp2d_NewtonCotez.py
# =============================================================================
import numpy as np
import pandas as pd

# Set display options
pd.set_option('display.precision', 6)
pd.set_option('display.width', 1000)
# Cell 3: Integration Helper (Same as Lagrange file)
def polynomial_integrate_definite(coeffs, a, b):
    """
    Integrates a polynomial P(x) defined by coeffs [c_n, ..., c_0] 
    from a to b.
    """
    n = len(coeffs) - 1
    int_coeffs = []
    
    for i, c in enumerate(coeffs):
        power = n - i
        int_coeffs.append(c / (power + 1))
    int_coeffs.append(0) # Constant C becomes 0
    
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
# Cell 4: Newton-Cotes Logic
def solve_newton_cotes(f, a, b, n):
    """
    Calculates definite integral using Newton-Cotes method of order n.
    """
    # 1. Discretization
    h = (b - a) / n
    x_vals = np.linspace(a, b, n + 1)
    y_vals = f(x_vals)
    
    # 2. Calculate Coefficients H_i (Integration on [0, n])
    # We use the variable t ranging from 0 to n. Nodes are integers 0, 1, ..., n.
    H_coeffs = []
    
    nodes_t = list(range(n + 1)) # [0, 1, ..., n]
    
    for i in range(n + 1):
        # Construct Basis Polynomial for node i in terms of t
        # L_i(t) = Product( (t - j)/(i - j) )
        
        poly_numerator = [1.0]
        denominator = 1.0
        
        for j in range(n + 1):
            if i == j:
                continue
            
            # Update Denominator: (i - j)
            denominator *= (i - j)
            
            # Update Numerator: Multiply by (t - j) -> Poly [1, -j]
            poly_numerator = np.convolve(poly_numerator, [1.0, -j])
            
        # Integrate Numerator from 0 to n
        integral_val = polynomial_integrate_definite(poly_numerator, 0, n)
        
        # Final H_i
        H_i = integral_val / denominator
        H_coeffs.append(H_i)
        
    # 3. Compute Final Sum: I = h * Sum(y_i * H_i)
    weighted_sum = np.sum(y_vals * np.array(H_coeffs))
    I = h * weighted_sum
    
    # Output Table
    df = pd.DataFrame({
        'Index k': range(n + 1),
        'x_k': x_vals,
        'y_k': y_vals,
        'Coeff H_k': H_coeffs,
        'Weight (h * H_k)': h * np.array(H_coeffs),
        'Term': y_vals * h * np.array(H_coeffs)
    })
    
    return df, I, h
## Result
# Cell 6: Example Usage
# Define Function: f(x) = 1/(x^2 + 1)
f = lambda x: 1/(x**2+1)
a = 0.0
b = 2.0
n = 4  # Newton-Cotes of order 4 (Boole's Rule equivalent if n=4)

print(f"\nStep size h = {h:.6f}")
print(f"--- Newton-Cotes Integration (Order n={n}) ---")
print(f"Interval: [{a}, {b}]")


df, I, h = solve_newton_cotes(f, a, b, n)
print("Check Sum of H_i (should be n):", df['Coeff H_k'].sum())
print("\nTable of Values and Weights:")

df


print(f"\nApproximate Integral I = {I:.8f}")




