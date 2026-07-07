# =============================================================================
# pp5_fixed_point_iteration.py - Phương pháp lặp điểm cố định (1 biến)
#
# Chức năng: Tìm nghiệm của phương trình x = φ(x) bằng cách lặp
#            x(n+1) = φ(x(n)) cho đến khi hội tụ.
#
# Các hàm chính:
#   fixed_point_iteration_v1(φ, dφ, a, b, x0, q, n, rbl)
#   fixed_point_iteration_v2(φ, dφ, a, b, x0, q, n, rbl)
#
# Cách dùng: python pp5_fixed_point_iteration.py
# =============================================================================
import numpy as np
import pandas as pd

pd.set_option('display.precision', 7)  # Độ chính xác hiển thị
pd.set_option('display.width', 150)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
import contextlib
from pathlib import Path

__dir__ = Path(__file__).parent

def fixed_point_iteration_v1 (phi, dphi, a, b, x_0, q, n, rbl):
    #Initial value
    x = x_0; results = [];
    results.append({
        'n': 0,
        'x_n': x_0,
        'f(x_n)': f(x_0),
        'delta=q^n/(1-q) * |x_1-x_0|': None
        }) 
    
    # Fixed-point iteration
    delta = 0; diff = 0; temp_q = q;
    for i in range(n):
        # Compute next iteration
        x_new = phi(x)

        if (i==0): 
            diff = abs(x_new - x)
        delta = temp_q/(1-q) * diff
        temp_q *= q
        
        results.append({
            'n': i+1,
            'x_n': x_new,
            'f(x_n)': f(x_new),
            'delta=q^n/(1-q) * |x_1-x_0|': delta
        }) 

        # Prepare for next iteration
        x = x_new
        if (f(x_new) == 0): 
            break

	#Print the final result
    df_results = pd.DataFrame(results)
    print(df_results.to_string(index=False))
		
    if rbl == None:
        print(f"Giá trị nghiệm là: {x}")
    else:
        total_delta = delta + 0.5 * 10**(-rbl)
        print(f"Giá trị nghiệm với {rbl} chữ số thập phân là: {round(x, rbl)}")
        print(f"Sai số tương đối là: {total_delta}")

def fixed_point_iteration_v2 (phi, dphi, a, b, x_0, q, n, rbl):
    #Initial value
    x = x_0; results = [];
    results.append({
        'n': 0,
        'x_n': x_0,
        'f(x_n)': f(x_0),
        'delta_x = q/(1-q) * |x_n - x|': None
        }) 
    
    # Fixed-point iteration
    delta_x = 0;
    for i in range(n):
        # Compute next iteration
        x_new = phi(x)
        delta_x = q/(1-q) * abs(x_new - x)
        
        results.append({
            'n': i+1,
            'x_n': x_new,
            'f(x_n)': f(x_new),
            'delta_x = q/(1-q) * |x_n - x|': delta_x
        }) 

        # Prepare for next iteration
        x = x_new
        if (f(x_new) == 0): 
            break

	#Print the final result
    df_results = pd.DataFrame(results)
    print(df_results.to_string(index=False))
		
    if rbl == None:
        print(f"Giá trị nghiệm là: {x}")
    else:
        total_delta = delta_x + 0.5 * 10**(-rbl)
        print(f"Giá trị nghiệm với {rbl} chữ số thập phân là: {round(x, rbl)}")
        print(f"Sai số tương đối là: {total_delta}")

def fixed_point_recursion_absolute (f, phi, dphi, a, b, x_0, q, eps):
    #Initial value
    x = x_0; results = [];
    results.append({
        'n': 0,
        'x_n': x_0,
        'f(x_n)': f(x_0),
        'delta_x=|x_new - x|': None
        }) 
    
    # Fixed-point iteration
    i=0; new_eps = (1-q)/q * eps
    print(f"new_eps = {new_eps}")
    
    while True:
        # Compute next iteration
        x_new = phi(x)
        current_delta_x = abs(x_new - x)
        
        results.append({
            'n': i+1,
            'x_n': x_new,
            'f(x_n)': f(x_new),
            'delta_x=|x_new - x|': current_delta_x
        }) 

        x = x_new  
        #stop condition 
        if (f(x_new) == 0): 
            break
        elif current_delta_x < new_eps:
            break
        else:
            i += 1

	#Print the final result
    df_results = pd.DataFrame(results)
    print(df_results.to_string(index=False))
		
    print(f"Giá trị nghiệm với sai số tuyệt đối {eps} là: {x}")

def fixed_point_recursion_relative (f, phi, dphi, a, b, x_0, q, eta):
    #Initial value
    x = x_0; results = [];
    results.append({
        'n': 0,
        'x_n': x_0,
        'f(x_n)': f(x_0),
        'sigma_x=|x_n-x_(n-1)| / |x_n|': None
        }) 
    
    # Fixed-point iteration
    i=0; new_eta = (1-q)/q * eta
    print(f"new_eta = {new_eta}")
    
    while True:
        # Compute next iteration
        x_new = phi(x)
        current_sigma_x = abs(x_new - x) / abs(x_new)
        
        results.append({
            'n': i+1,
            'x_n': x_new,
            'f(x_n)': f(x_new),
            'sigma_x=|x_n-x_(n-1)| / |x_n|': current_sigma_x
        }) 

        #stop condition 
        if (f(x) == 0): 
            break
        elif current_sigma_x < new_eta:
            break
        else:
            i += 1
            x = x_new  

	#Print the final result
    df_results = pd.DataFrame(results)
    print(df_results.to_string(index=False))
		
    print(f"Giá trị nghiệm với sai số tương đối {eps} là: {x}")

if __name__ == "__main__":
    output_path = str(__dir__ / "pp5_fixed_point_iteration_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        f = lambda x: x**5 - 17*x + 2
        phi = lambda x: (x**5 + 2) / 17
        dphi = lambda x: 5*x**4 / 17

        a = 0
        b = 1
        x_0 = a
        q = max(np.abs(dphi(x)) for x in [a, b]) # Lipschitz constant (q<1)

        n = 15
        rbl = 9

        fixed_point_iteration_v1 (phi, dphi, a, b, x_0, q, n, rbl)

        f = lambda x: x**5 - 17*x + 2
        phi = lambda x: (x**5 + 2) / 17
        dphi = lambda x: 5*x**4 / 17

        a = 0
        b = 1
        x_0 = a
        q = max(np.abs(dphi(x)) for x in [a, b]) # Lipschitz constant (q<1)

        n = 15
        rbl = 9

        fixed_point_iteration_v2 (phi, dphi, a, b, x_0, q, n, rbl)

        f = lambda x: (1.4)**x - x
        phi = lambda x: (1.4)**x
        dphi = lambda x: (1.4)**x * np.log(1.4)

        a = 1
        b = 2
        x_0 = b
        q = dphi(2)

        eps = 0.5 * 10**(-7)

        fixed_point_recursion_absolute (f, phi, dphi, a, b, x_0, q, eps)
        f = lambda x: x**5 - 17*x + 2
        phi = lambda x: -pow(2 - 17*x , 0.2)
        dphi = lambda x: 0.2 * 17 * pow(2 - 17*x, -0.8)

        a = -3
        b = -2
        x_0 = b
        q = max(np.abs(dphi(x)) for x in [a, b]) # Lipschitz constant (q<1)

        eps = 10**(-8)

        fixed_point_recursion_absolute (f, phi, dphi, a, b, x_0, q, eps)

        f = lambda x: x**2 + 3*x - 1
        phi = lambda x: 1/(x+3)
        dphi = lambda x: -1/(x+3)**2

        a = 0
        b = 1
        x_0 = a
        q = max(np.abs(dphi(x)) for x in [a, b]) # Lipschitz constant (q<1)

        eta = 5 * 10**(-9)

        fixed_point_recursion_relative (f, phi, dphi, a, b, x_0, q, eta)
    print(f"Đã ghi kết quả vào {output_path}")
