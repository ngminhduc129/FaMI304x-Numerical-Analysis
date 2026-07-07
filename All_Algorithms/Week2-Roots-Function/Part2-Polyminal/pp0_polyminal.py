# =============================================================================
# pp0_polyminal.py - Tiện ích đa thức (Polynomial Utilities)
#
# Chức năng: Cung cấp các hàm xử lý đa thức: tính giá trị, đạo hàm,
#            chặn nghiệm, tìm nghiệm thực/phức.
#
# Các hàm chính:
#   polynomial(coeffs, x)                  - tính giá trị tại x
#   polynomial_derivative(coeffs, x)       - đạo hàm f'(x)
#   polynomial_second_derivative(coeffs, x)- đạo hàm f''(x)
#   complex_radius(coeffs)                 - bán kính chặn nghiệm phức
#   real_radius(coeffs)                    - bán kính chặn nghiệm thực
#   real_root_bound(coeffs, s)             - cận trên nghiệm thực
#   solve_complex_roots(degree, ...)       - giải nghiệm phức
#
# Cách dùng: python pp0_polyminal.py
# =============================================================================
import numpy as np
import pandas as pd
import contextlib
from pathlib import Path

__dir__ = Path(__file__).parent

def polynomial(coefficients, x):
    return sum(c * x**i for i, c in enumerate(coefficients))

def polynomial_derivative(coefficients, x):
    return sum(i * c * x**(i-1) for i, c in enumerate(coefficients) if i > 0)

def polynomial_second_derivative(coefficients, x):
    return sum(i * (i-1) * c * x**(i-2) for i, c in enumerate(coefficients) if i > 1)
def complex_radius(coefficients):
    R = 1 + max(abs(coeff) for coeff in coefficients) / abs(coefficients[-1])
    return R

def real_radius(coefficients):
    n = len(coefficients) - 1
    negative_coefficients = [(i, abs(x)) for i, x in enumerate(coefficients) if x<0]
    if not negative_coefficients:
        return None
    else:
        k = (max(negative_coefficients, key=lambda item: item[0]))[0]
        ratio = (max(negative_coefficients, key=lambda item: item[1]))[1] / abs(coefficients[-1])
    
        R = 1 + ratio**(1/(n-k)) if n!=k else 0
        return R
def find_extrema(coefficients, lower_bound, upper_bound, learning_rate=0.01, max_iterations=1000, tolerance=1e-6, num_starting_points=10):
    """
    Find local minima and maxima using gradient descent/ascent with adaptive learning rate
    Double-check the answer with your plotting tool. If answer is not correct, try changing the learning rate or number of starting points.
    """
    starting_points = np.linspace(lower_bound, upper_bound, num_starting_points)
    extrema = {'minima': [], 'maxima': []}
    
    def gradient_search(x0, direction=1):  # direction: 1 for maxima, -1 for minima
        x = x0
        current_lr = learning_rate
        prev_value = polynomial(coefficients, x)
        
        for _ in range(max_iterations):
            gradient = polynomial_derivative(coefficients, x)
            x_new = x + direction * current_lr * gradient
            
            # Stay within bounds
            x_new = np.clip(x_new, lower_bound, upper_bound)
            
            # Check if we're moving in the wrong direction
            current_value = polynomial(coefficients, x_new)
            if (direction == 1 and current_value < prev_value) or \
               (direction == -1 and current_value > prev_value):
                # We went too far, reduce learning rate and try again
                current_lr *= 0.5
                continue
                
            # Update position if improvement is made
            if abs(x_new - x) < tolerance:
                # Check if it's actually an extremum using second derivative
                second_deriv = polynomial_second_derivative(coefficients, x_new)
                if (direction == -1 and second_deriv > 0) or (direction == 1 and second_deriv < 0):
                    return x_new
                return None
                
            x = x_new
            prev_value = current_value
            
        return None

    # Rest of the function remains the same
    for x0 in starting_points:
        # Find minimum
        min_point = gradient_search(x0, direction=-1)
        if min_point is not None:
            if not any(abs(min_point - x) < tolerance for x in extrema['minima']):
                extrema['minima'].append(min_point)
        
        # Find maximum
        max_point = gradient_search(x0, direction=1)
        if max_point is not None:
            if not any(abs(max_point - x) < tolerance for x in extrema['maxima']):
                extrema['maxima'].append(max_point)
    
    # Sort the results
    extrema['minima'].sort()
    extrema['maxima'].sort()
    return extrema
if __name__ == "__main__":
    output_path = str(__dir__ / "pp0_polyminal_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        # Example usage
        #coefficients = [-15, 0, 3, 0, 3, 1] # Coefficients of the polynomial -15 + 3x^2 - ax^4 + x^5
        #hệ số của fx
        coefficients = [-23, 0, 3, 1]
        temp = [x if i%2==0 else -x for i, x in enumerate(coefficients)]
        if (temp[-1] < 0):
            revert_coefficients = [-i for i in temp]
        else:
            revert_coefficients = temp
        
        radius = complex_radius(coefficients)
        upper_radius = real_radius(coefficients)
        lower_radius = -real_radius(revert_coefficients)
        
        print(f"Bán kính phức cho nghiệm: {-radius} -> {radius}")
        print(f"Giá trị hàm tại mỗi đầu của bán kính phức: {polynomial(coefficients, -radius)} và {polynomial(coefficients, radius)}")
        print("\n")
        
        lower_result = polynomial(coefficients, lower_radius) if lower_radius else None
        print("Bán kính dưới cho nghiệm thực:", lower_radius)
        print("Giá trị hàm tại bán kính dưới:", lower_result)
        print("\n")
        
        upper_result = polynomial(coefficients, upper_radius) if upper_radius else None
        print("Bán kính trên cho nghiệm thực:", upper_radius)
        print("Giá trị hàm tại bán kính trên:", upper_result)
        
        #To optimize, first remove the free coefficent cause it doees't affect the extrema
        temp_coefficients = coefficients.copy()
        temp_coefficients[0] = 0
        extrema = find_extrema(temp_coefficients, lower_bound = lower_radius, upper_bound = upper_radius, learning_rate = 0.5)
        
        print("Cực tiểu địa phương:", [f"{x}" for x in extrema['minima']])
        print("Giá trị hàm tại cực tiểu:", [f"{polynomial(coefficients, x)}" for x in extrema['minima']])
        print("\n")
        
        print("Cực đại địa phương:", [f"{x}" for x in extrema['maxima']])
        print("Giá trị hàm tại cực đại:", [f"{polynomial(coefficients, x)}" for x in extrema['maxima']])
    print(f"Đã ghi kết quả vào {output_path}")
