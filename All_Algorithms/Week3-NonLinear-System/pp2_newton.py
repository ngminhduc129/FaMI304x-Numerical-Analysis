# =============================================================================
# pp2_newton.py - Phương pháp Newton cho hệ phi tuyến
#
# Chức năng: Giải hệ phương trình phi tuyến f1=0, f2=0, ...
#            bằng phương pháp Newton với ma trận Jacobian.
#
# Các hàm chính:
#   newton_method(initial_values, tol, max_iter, *funcs)
#     - funcs bao gồm các hàm f1,f2,... và hàm Jacobian
#
# Cách dùng: python pp2_newton.py
# =============================================================================
import numpy as np
import pandas as pd
import contextlib
from pathlib import Path
__dir__ = Path(__file__).parent.resolve()

pd.set_option('display.precision', 7)  # Độ chính xác hiển thị
pd.set_option('display.width', 150)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def newton_method(initial_values, tol=1e-6, max_iter=100, *funcs):
    n = len(initial_values)  # Number of initial values (and functions)
    x = np.array(initial_values, dtype=float)  # Convert initial values to a NumPy array
    results = [[0] + list(x)]  # Initialize results with the initial values and a placeholder for the norm difference
    
    for i in range(max_iter):
        F = np.array([func(*x) for func in funcs[:n]])
        
        # Create the Jacobian matrix using pandas
        J = pd.DataFrame(funcs[n](*x))
        
        delta_x = np.linalg.solve(J.values, -F)
        x_new = x + delta_x
        
        # Calculate the L1 norm (sum of absolute differences)
        norm_diff = np.sum(np.abs(x_new - x))
        results.append([i+1] + list(x_new) + [norm_diff])
        x = x_new
        
        if norm_diff < tol:
            break
    
    columns = ['Lần lặp'] + [f'x{j+1}' for j in range(n)] + ['sai_số_chuẩn']
    df = pd.DataFrame(results, columns=columns)
    print(df.to_string(index=False))
    
    return x

if __name__ == "__main__":
    output_path = str(__dir__ / "pp2_newton_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        # Example functions
        def f1(x, y, z):
            return 3*x - np.cos(y*z) - 0.5

        def f2(x, y, z):
            return x**2 - 81*((y+0.1)**2) + np.sin(z) + 1.06

        def f3(x, y, z):
            return np.exp(-x*y) + 20*z + 9.1389 

        # Jacobian matrix function
        def jacobian(x, y, z):
            return [
                [3, z*np.sin(y*z), y*np.sin(y*z)],
                [2*x, -162*(y+0.1), np.cos(z)],
                [-y*np.exp(-x*y), -x*np.exp(-x*y), 20]
            ]

        # Store functions and Jacobian function
        funcs = [f1, f2, f3, jacobian]

        # Initial guess
        initial_values = [0, 0, 0]

        # Perform Newton-Raphson method
        try:
            solution = newton_method(initial_values, 1e-6, 100, *funcs)
            print("Nghiệm xấp xỉ:", solution)
        except np.linalg.LinAlgError:
            print("Ma trận suy biến – không thể giải hệ phương trình.")
        def f1(x, y):
            return x*x + 2*y*y -4

        def f2(x, y):
            return 5*y*y + x*y -4

        # Jacobian matrix function
        def jacobian(x, y):
            return [
                [2*x , 4*y],
                [y , (10*y + x)]
            ]

        # Store functions and Jacobian function
        funcs = [f1, f2, jacobian]

        # Initial guess
        initial_values = [0, 0]

        # Perform Newton-Raphson method
        try:
            solution = newton_method(initial_values, 1e-6, 100, *funcs)
            print("Nghiệm xấp xỉ:", solution)
        except np.linalg.LinAlgError:
            print("Ma trận suy biến – không thể giải hệ phương trình.")
    print(f"Đã ghi kết quả vào {output_path}")
