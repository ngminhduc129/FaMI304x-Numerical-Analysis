# =============================================================================
# pp2_secant.py - Phương pháp dây cung (Secant Method)
#
# Chức năng: Tìm nghiệm của f(x)=0 bằng cách thay đường cong f(x)
#            bằng dây cung qua hai điểm (a, f(a)) và (b, f(b)).
#
# Các hàm chính:
#   secant_iteration_v1(f, df, a, b, n, rbl)     - có công thức sai số
#   secant_iteration_v2(f, a, b, n, rbl)         - không cần df, dùng delta_x
#   secant_iteration_v3(f, a, b, eps, rbl)       - dùng epsilon
#
# Cách dùng: python pp2_secant.py
# =============================================================================
import numpy as np
import pandas as pd

pd.set_option('display.precision', 7)  # Độ chính xác hiển thị
pd.set_option('display.width', 150)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
import contextlib
from pathlib import Path

__dir__ = Path(__file__).parent

def secant_iteration_v1 (f, df, a, b, n, rbl):	
	M1 = max([np.abs(df(x)) for x in [a, b]]) #M1 is the maximum value of |f'(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print(f"m1 = {m1}, M1 = {M1}")

	# Error function
	if (f(a) * f(b) >= 0) or (M1 * m1 <= 0): 
		print("Không thể tìm nghiệm trong khoảng đã cho");
		return 0;

	# Implementing Secant Method
	x = a; mrk = b;

	results = [];
	for i in range(n):
		# Calculate the next value of x
		x_new = (mrk * f(x) - x * f(mrk)) / (f(x) - f(mrk))
		delta_x = abs(f(x_new))/ m1

		results.append({
            'n': i,
			'mrk': mrk,
            'x_(n+1)': x_new,
            'f(x_(n+1))': f(x_new),
            'delta_x=|f(x_n)| / m_1': delta_x
        })
		
		# Prepare for next iteration
		x = x_new
		if (f(x_new) == 0): 
			break
		elif (i == 0):
			if f(a) * f(x_new) < 0:
				mrk = a
			else:
				mrk = b

	# Print the final result
	df_results = pd.DataFrame(results)
	print(df_results.to_string(index=False))

	if rbl == None:
		print(f"Giá trị nghiệm là: {x}")
	else:
		total_delta = delta_x + 0.5 * 10**(-rbl) #must calculate roundoff error
		print(f"Giá trị nghiệm với {rbl} chữ số thập phân là: {round(x, rbl)}")
		print(f"Sai số tương đối là: {total_delta}")

def secant_iteration_v2 (f, df, a, b, n, rbl):	
	M1 = max([np.abs(df(x)) for x in [a, b]]) #M1 is the maximum value of |f'(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print(f"m1 = {m1}, M1 = {M1}")

	# Error function
	if (f(a) * f(b) >= 0) or (M1 * m1 <= 0): 
		print("Không thể tìm nghiệm trong khoảng đã cho");
		return 0;

	# Implementing Secant Method
	x = a; mrk = b;

	results = [];
	for i in range(n):
		# Calculate the next value of x
		x_new = (mrk * f(x) - x * f(mrk)) / (f(x) - f(mrk))
		delta_x = (M1 - m1) * abs(x - x_new) / m1

		results.append({
            'n': i,
			'mrk': mrk,
            'x_(n+1)': x_new,
            'f(x_(n+1))': f(x_new),
            'delta_x=(M1 - m1)*(x-x_new)/m1': delta_x
        })
		
		# Prepare for next iteration
		x = x_new
		if (f(x_new) == 0): 
			break
		elif (i == 0):
			if f(a) * f(x_new) < 0:
				mrk = a
			else:
				mrk = b

	# Print the final result
	df_results = pd.DataFrame(results)
	print(df_results.to_string(index=False))

	if rbl == None:
		print(f"Giá trị nghiệm là: {x}")
	else:
		total_delta = delta_x + 0.5 * 10**(-rbl) #must calculate roundoff error
		print(f"Giá trị nghiệm với {rbl} chữ số thập phân là: {round(x, rbl)}")
		print(f"Sai số tương đối là: {total_delta}")

def secant_recursion_absolute_v1 (f, df, a, b, eps):	
	M1 = max([np.abs(df(x)) for x in [a, b]]) #M1 is the maximum value of |f'(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print(f"m1 = {m1}, M1 = {M1}")

	# Error function
	if (f(a) * f(b) >= 0) or (M1 * m1 <= 0): 
		print("Không thể tìm nghiệm trong khoảng đã cho");
		return 0;

	# Implementing Secant Method
	x = a; mrk = b; new_eps = m1*eps;
	print(f"delta_x = {new_eps}")

	results = []; i = 0; 
	while True:
		# Calculate the next value of x
		x_new = (mrk * f(x) - x * f(mrk)) / (f(x) - f(mrk))
		current_delta_x = abs(f(x_new))

		results.append({
            'n': i,
			'mrk': mrk,
            'x_(n+1)': x_new,
            'f(x_(n+1))': f(x_new),
            'delta_x=|f(x_n)|': current_delta_x
        })
		
		# Prepare for next iteration
		x = x_new
		if (f(x_new) == 0): 
			break
		elif (i == 0):
			if f(a) * f(x_new) < 0:
				mrk = a
			else:
				mrk = b
				
        # Stopping condition
		if current_delta_x < new_eps:
			break
		else:
			i += 1

    # Print the final result
	df_results = pd.DataFrame(results) 
	print(df_results.to_string(index=False))

	print(f"Giá trị nghiệm với sai số tuyệt đối {eps} là: {x}")

def secant_recursion_absolute_v2 (f, df, a, b, eps):	
	M1 = max([np.abs(df(x)) for x in [a, b]]) #M1 is the maximum value of |f'(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print(f"m1 = {m1}, M1 = {M1}")

	# Error function
	if (f(a) * f(b) >= 0) or (M1 * m1 <= 0): 
		print("Không thể tìm nghiệm trong khoảng đã cho");
		return 0;

	# Implementing Secant Method
	x = a; mrk = b; new_eps = eps*m1/(M1 - m1);
	print(f"delta_x = {new_eps}")

	results = []; i = 0; 
	while True:
		# Calculate the next value of x
		x_new = (mrk * f(x) - x * f(mrk)) / (f(x) - f(mrk))
		current_delta_x = abs(x - x_new)

		results.append({
            'n': i,
			'mrk': mrk,
            'x_(n+1)': x_new,
            'f(x_(n+1))': f(x_new),
            'delta_x=|x_(n+1) - x_n|': current_delta_x
        })
		
		# Prepare for next iteration
		x = x_new
		if (f(x_new) == 0): 
			break
		elif (i == 0):
			if f(a) * f(x_new) < 0:
				mrk = a
			else:
				mrk = b
				
        # Stopping condition
		if current_delta_x < new_eps:
			break
		else:
			i += 1

    # Print the final result
	df_results = pd.DataFrame(results) 
	print(df_results.to_string(index=False))

	print(f"Giá trị nghiệm với sai số tuyệt đối {eps} là: {x}")

def secant_recursion_relative (f, df, a, b, eta):	
	M1 = max([np.abs(df(x)) for x in [a, b]]) #M1 is the maximum value of |f'(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print(f"m1 = {m1}, M1 = {M1}")

	# Error function
	if (f(a) * f(b) >= 0) or (M1 * m1 <= 0): 
		print("Không thể tìm nghiệm trong khoảng đã cho");
		return 0;

	# Implementing Secant Method
	x = a; mrk = b; new_eta = eta * m1/(M1 - m1);
	print(f"sigma_x = {new_eta}, m_1 = {m1}, M_1 = {M1}")

	results = []; i = 0;
	while True:
		# Calculate the next value of x
		x_new = (mrk * f(x) - x * f(mrk)) / (f(x) - f(mrk))
		current_sigma_x = abs(x - x_new)/abs(x_new)

		results.append({
            'n': i,
			'mrk': mrk,
            'x_(n+1)': x_new,
            'f(x_(n+1))': f(x_new),
            'sigma_x=|x_(n+1)-x_n|/|x_n|': current_sigma_x
        })
		
		# Prepare for next iteration
		x = x_new
		if (f(x_new) == 0): 
			break
		elif (i == 0):
			if f(a) * f(x_new) < 0:
				mrk = a
			else:
				mrk = b
				
        # Stopping condition
		if current_sigma_x < new_eta:
			break
		else:
			i += 1

    # Print the final result
	df_results = pd.DataFrame(results) 
	print(df_results.to_string(index=False))

	print(f"Giá trị nghiệm với sai số tương đối {eta} là: {x}")

if __name__ == "__main__":
    output_path = str(__dir__ / "pp2_secant_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        f = lambda x: x ** 5 - 12; 
        df = lambda x: 5 * (x ** 4)

        a = 1
        b = 2

        n = 20
        rbl = 9;
        secant_iteration_v1 (f, df, a, b, n, rbl)

        f = lambda x: x ** 5 - 12; 
        df = lambda x: 5 * (x ** 4)

        a = 1
        b = 2

        n = 20
        rbl = 9;
        secant_iteration_v2 (f, df, a, b, n, rbl)

        f = lambda x: x**5 - 81*x - 243 #approximate e
        df = lambda x: 20*(x**3)

        a = 3
        b = 6

        eps = 3 * 1e-6

        secant_recursion_absolute_v1 (f, df, a, b, eps)
        f = lambda x: np.tan(x/4) - 1 #approximate pi
        df = lambda x: (1/4) * 1/(np.cos(x/4)**2)

        a = 3
        b = 4

        eps = 0.5 * pow(10, -7)

        secant_recursion_absolute_v1 (f, df, a, b, eps)

        f = lambda x: x**5 - 17
        df = lambda x: 5 * x**4

        a = 1
        b = 2

        eps = 0.5 * pow(10, -6)

        secant_recursion_absolute_v2 (f, df, a, b, eps)
        f = lambda x: np.e**x - np.cos(2*x)
        df = lambda x: np.e**x + 2*np.sin(2*x)

        a = -0.5
        b = -0.1

        eps = 0.5 * pow(10, -8)

        secant_recursion_absolute_v2 (f, df, a, b, eps)

        f = lambda x: x**5 - 3*x**3 + 2*x**2 - x + 5
        df = lambda x: 5*x**4 - 9*x**2 + 4*x - 1

        a = -3
        b = -2

        eta = 0.05 * pow(10, -2)

        secant_recursion_relative (f, df, a, b, eta)
        f = lambda x: np.e**(-x) - x
        df = lambda x: -np.e**(-x) - 1

        a = 0
        b = 1

        eta = 0.05 * pow(10, -2)

        secant_recursion_relative (f, df, a, b, eta)
        f = lambda x: 2**x - 5*x + np.sin(x)
        df = lambda x: 2**x*np.log(2) - 5 + np.cos(x)

        a = 0
        b = 1

        eta = 0.05 * pow(10, -2)

        secant_recursion_relative (f, df, a, b, eta)
    print(f"Đã ghi kết quả vào {output_path}")
