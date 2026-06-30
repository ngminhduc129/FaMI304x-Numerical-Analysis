# =============================================================================
# pp1_bisection.py - Phương pháp chia đôi (Bisection Method)
#
# Chức năng: Tìm nghiệm của phương trình f(x) = 0 trên [a, b] với
#            f(a)*f(b) < 0 bằng cách chia đôi khoảng liên tục.
#
# Các hàm chính:
#   bisection_iteration_v1(f, a, b, n, rbl)    - n bước, làm tròn rbl chữ số
#   bisection_iteration_v2(f, a, b, n, rbl)    - dùng delta_x (sai số tuyệt đối)
#   bisection_iteration_v3(f, a, b, n, rbl)    - dùng sigma (sai số tương đối)
#   bisection_iteration_v4(f, a, b, eps, rbl)  - dùng epsilon làm điều kiện dừng
#   bisection_iteration_v5(f, a, b, sigma, rbl)- dùng sai số tương đối sigma
#
# Cách dùng:
#   1. Định nghĩa hàm f(x)
#   2. Gọi một trong các hàm trên
#   3. Chạy: python pp1_bisection.py
#
# Ví dụ (đã có sẵn cuối file):
#   f = lambda x: e**x - cos(2*x)
#   bisection_iteration_v1(f, -3, -2, 20, 5)
# =============================================================================
import numpy as np
import pandas as pd

pd.set_option('display.precision', 15)  # Increase decimal precision
pd.set_option('display.width', 150)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column

sign = lambda x: 1 if x > 0 else (-1 if x < 0 else 0)
def bisection_iteration_v1 (f,a,b,n,rbl):
    # Error function
    if (f(a) * f(b) >= 0):
        print("You have not assumed right a and b\n")
        return

    # Implementing Bisection Method
    x = 0; delta = 0; sign_f_a = sign(f(a));

    results = []; diff = b-a; temp_2 = 2;
    for i in range(n):
        # Find next value of x
        x_new = (a+b) / 2.0
        delta = diff / temp_2
        temp_2 *= 2

        results.append({
            'n': i,
            'a_n': a,
            'b_n': b,
            'x_(n+1)': x_new,
            'f(a_n)': f(a),
            'f(b_n)': f(b),
            'f(x_(n+1))': f(x_new),
            'delta=(b-a)/2^(n+1)': delta
        })

        # Prepare for next iteration
        x = x_new
        if f(x_new) == 0:
            break
        elif (f(x_new) * sign_f_a < 0):
            b = x_new
        elif (f(x_new) * sign_f_a > 0):
            a = x_new

    # Print the final result
    df_results = pd.DataFrame(results)
    print(df_results.to_string(index=False))

    if rbl == None:
        print(f"The value of root is: {x}")
    else:
        total_delta = delta + 0.5 * 10**(-rbl) #must calculate roundoff error
        print(f"The value of root with {rbl} decimal point is: {round(x, rbl)}")
        print(f"Relative error is: {total_delta}")
f = lambda x: np.e**x - np.cos(2*x)

a = -3
b = -2

n = 20
rbl = 5

bisection_iteration_v1 (f, a, b, n, rbl)
f = lambda x: np.log(x) - 1 #approximate e

a = 2
b = 3

eps = 0.5 * pow(10, -7) #epsilon should be small so that the root is accurate to 7 decimal places
n = int(np.floor(np.log2((b-a)/eps)) + 1)
rbl = 7

bisection_iteration_v1 (f, a, b, n, rbl)
f = lambda x: np.tan(x/4) - 1 #approximate pi

a = 3
b = 3.2

eps = 0.5 * pow(10, -7)
n = int(np.floor(np.log2((b-a)/eps)) + 1)
rbl = 7

bisection_iteration_v1 (f, a, b, n, rbl)
def bisection_iteration_v2 (f,a,b,n,rbl):
    # Error function
    if (f(a) * f(b) >= 0):
        print("You have not assumed right a and b\n")
        return

    # Implementing Bisection Method
    x = 0; delta_x = 0; sign_f_a = sign(f(a));

    results = []
    for i in range(n):
        # Find next value of x
        x_new = (a+b) / 2.0
        delta_x = abs(x_new - x)

        results.append({
            'n': i,
            'a_n': a,
            'b_n': b,
            'x_(n+1)': x_new,
            'f(a_n)': f(a),
            'f(b_n)': f(b),
            'f(x_(n+1))': f(x_new),
            'delta_x=|x_(n+1)-x_n|': delta_x
        })

        # Prepare for next iteration
        x =  x_new
        if f(x_new) == 0:
            break
        elif (f(x_new) * sign_f_a < 0):
            b = x_new
        elif (f(x_new) * sign_f_a > 0):
            a = x_new

    # Print the final result
    df_results = pd.DataFrame(results)
    print(df_results.to_string(index=False))

    if rbl == None:
        print(f"The value of root is: {x}")
    else:
        total_delta = delta_x + 0.5 * 10**(-rbl) #must calculate roundoff error
        print(f"The value of root with {rbl} decimal point is: {round(x, rbl)}")
        print(f"Relative error is: {total_delta}")
f = lambda x: np.e**x - np.cos(2*x)

a = -3
b = -2

n = 30
rbl = 5

bisection_iteration_v2 (f, a, b, n, rbl)
def bisection_recursion_absolute (f,a,b,eps):
    # Error function
    if (f(a) * f(b) >= 0):
        print("You have not assumed right a and b\n")
        return

    # Implementing Bisection Method
    x = 0; sign_f_a = sign(f(a));
    print(f"delta_x = {eps}")

    i=0; results = []
    while True:
        # Find next value of x
        x_new = (a+b) / 2.0
        current_delta_x = abs(x_new - x)

        results.append({
            'n': i,
            'a_n': a,
            'b_n': b,
            'x_(n+1)': x_new,
            'f(a_n)': f(a),
            'f(b_n)': f(b),
            'f(x_(n+1))': f(x_new),
            'delta_x=|x_(n+1)-x_n|': current_delta_x
        })

        # Prepare for next iteration
        x =  x_new
        if f(x_new) == 0:
            break
        elif (f(x_new) * sign_f_a < 0):
            b = x_new
        elif (f(x_new) * sign_f_a > 0):
            a = x_new

        # Stopping condition
        if current_delta_x < eps:
            break
        else:
            i += 1

    # Print the final result
    df_results = pd.DataFrame(results)
    print(df_results.to_string(index=False))

    print(f"The value of root with absolute error {eps} is: {x}")
f = lambda x: 3*np.sin(x) + x**3 - 8*x**2 + 8*x + 1

a = 6
b = 7

eps = 0.5 * pow(10, -7)

bisection_recursion_absolute (f, a, b, eps)
f = lambda x: x**5-17 #approximate 5th root of 17

a = 1
b = 2

eps = 0.5 * pow(10, -6)

bisection_recursion_absolute (f, a, b, eps)
def bisection_recursion_relative (f,a,b,eta):
    # Error function
    if (f(a) * f(b) >= 0):
        print("You have not assumed right a and b\n")
        return

    # Implementing Bisection Method
    x = 0; sign_f_a = sign(f(a));
    print(f"sigma_x = {eta}")

    i=0; results = []
    while True:
        # Find next value of x
        x_new = (a+b) / 2.0
        current_sigma_x = abs(x_new - x)/abs(x_new)

        results.append({
            'n': i,
            'a_n': a,
            'b_n': b,
            'x_(n+1)': x_new,
            'f(a_n)': f(a),
            'f(b_n)': f(b),
            'f(x_(n+1))': f(x_new),
            'sigma_x=|x_(n+1)-x_n|/|x_(n+1)': current_sigma_x
        })

        # Prepare for next iteration
        x = x_new
        if f(x_new) == 0:
            break
        elif (f(x_new) * sign_f_a < 0):
            b = x_new
        elif (f(x_new) * sign_f_a > 0):
            a = x_new

        # Stopping condition
        if current_sigma_x < eta:
            break
        else:
            i += 1

    # Print the final result
    df_results = pd.DataFrame(results)
    print(df_results.to_string(index=False))

    print(f"The value of root with relative error {eta} is: {x}")
f = lambda x: np.exp(x) - np.cos(2*x)

a = -3
b = -2

eta = 0.5 * pow(10, -6)

bisection_recursion_relative (f, a, b, eta)
f = lambda x: x**5 - 3*x**3 + 2*x**2 - x + 5

a = -3
b = -2

eta = 0.05 * pow(10, -2)

bisection_recursion_relative (f, a, b, eta)




