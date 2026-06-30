# =============================================================================
# pp3_newton_raphson.py - Phương pháp tiếp tuyến (Newton-Raphson)
#
# Chức năng: Tìm nghiệm f(x)=0 bằng công thức lặp
#            x(n+1) = x(n) - f(x(n))/f'(x(n)).
#
# Các hàm chính:
#   newton_iteration_v1(f, df, d2f, a, b, n, rbl)
#   newton_iteration_v2(f, df, a, b, n, rbl)
#   newton_iteration_v3(f, df, a, b, eps, rbl)
#
# Cách dùng: python pp3_newton_raphson.py
# =============================================================================
import numpy as np
import pandas as pd

pd.set_option('display.precision', 15)  # Increase decimal precision
pd.set_option('display.width', 150)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def newton_iteration_v1 (f, df, d2f, a, b, n, rbl):	
	M2 = max([np.abs(d2f(x)) for x in [a, b]]) #M2 is the maximum value of |f''(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print (f"m1 = {m1}, M2 = {M2}")

	#Starting values
	if f(a) * d2f(a) > 0:
		x = a
	elif f(b) * d2f(b) > 0:
		x = b
	
	results = []
	results.append({
        'n': 0,
        'x_n': x,
        'f(x_n)': f(x),
        'delta_x=|f(x_n)| / m_1': abs(f(x)) / m1
    })

	#Newton's method
	x_new = 0; delta_x = 0;
	for i in range(n):
		# calculate the next value of x
		x_new = x - f(x) / df(x)
		delta_x = abs(f(x_new)) / m1

		results.append({
            'n': i+1,
            'x_n': x_new,
            'f(x_n)': f(x_new),
            'delta_x=|f(x_n)| / m_1': delta_x
        })

		# Prepare for next iteration
		x = x_new
		if (f(x_new) == 0): 
			break

	#Print the final result
	df_results = pd.DataFrame(results)
	print(df_results.to_string(index=False))

	if rbl == None:
		print(f"The value of root is: {x}")
	else:
		total_delta = delta_x + 0.5 * 10**(-rbl)
		print(f"The value of root with {rbl} decimal point is: {round(x, rbl)}")
		print(f"Relative error is: {total_delta}")
f = lambda x: x**4-27
df = lambda x: 4*x**3
d2f = lambda x: 12*x**2

a = 2
b = 3

n = 5
rbl = 9

newton_iteration_v1 (f, df, d2f, a, b, n, rbl);
def newton_iteration_v2 (f, df, d2f, a, b, n, rbl):	
	M2 = max([np.abs(d2f(x)) for x in [a, b]]) #M2 is the maximum value of |f''(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print (f"m1 = {m1}, M2 = {M2}")

	#Starting values
	if f(a) * d2f(a) > 0:
		x = a
	elif f(b) * d2f(b) > 0:
		x = b
	
	results = []
	results.append({
        'n': 0,
        'x_n': x,
        'f(x_n)': f(x),
        'delta_x=M_2/(2*m_1) * |x_n - x_(n-1)|^2': None
    })

	#Newton's method
	x_new = 0; delta_x = 0;
	for i in range(n):
		# calculate the next value of x
		x_new = x - f(x) / df(x)
		delta_x = M2 / (2 * m1) * (x_new - x)**2

		results.append({
            'n': i+1,
            'x_n': x_new,
            'f(x_n)': f(x_new),
            'delta=M_2/(2*m_1) * |x_n - x_(n-1)|^2': delta_x
        })

		# Prepare for next iteration
		x = x_new
		if (f(x_new) == 0): 
			break

	#Print the final result
	df_results = pd.DataFrame(results)
	print(df_results.to_string(index=False))

	if rbl == None:
		print(f"The value of root is: {x}")
	else:
		total_delta = delta_x + 0.5 * 10**(-rbl)
		print(f"The value of root with {rbl} decimal point is: {round(x, rbl)}")
		print(f"Relative error is: {total_delta}")
f = lambda x: x**4-27
df = lambda x: 4*x**3
d2f = lambda x: 12*x**2

a = 2
b = 3

n = 5
rbl = 9

newton_iteration_v2 (f, df, d2f, a, b, n, rbl);
def newton_recursion_absolute_v1 (f, df, d2f, a, b, eps):	
	M2 = max([np.abs(d2f(x)) for x in [a, b]]) #M2 is the maximum value of |f''(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print (f"m1 = {m1}, M2 = {M2}")

	#Starting values
	if f(a) * d2f(a) > 0:
		x = a
	elif f(b) * d2f(b) > 0:
		x = b
	
	results = []
	results.append({
        'n': 0,
        'x_n': x,
        'f(x_n)': f(x),
		"delta_x=|f(x_n)|": abs(f(x))
    })

	#Newton's method
	x_new = 0; new_eps = m1*eps;
	print(f"delta_x = {new_eps}")

	i=0;
	while True:
		# calculate the next value of x
		x_new = x - f(x) / df(x)
		current_delta_x = abs(f(x_new))

		results.append({
            'n': i+1,
            'x_n': x_new,
            'f(x_n)': f(x_new),
            'delta_x=|f(x_n)|': current_delta_x
        })
		
		# update the value of interval 
		x = x_new
		if (f(x_new) == 0): 
			break

		#stop condition
		if current_delta_x < new_eps:
			break
		else:
			i += 1

	#Print the final result
	df_results = pd.DataFrame(results)
	print(df_results.to_string(index=False))


	print(f"The value of root with absolute error {eps} is: {x}")
f = lambda x: x**5 - np.log(x) - 12
df = lambda x: 5*x**4 - 1/x
d2f = lambda x: 20*x**3 + 1/x**2

a = 1
b = 2

eps = 10**(-8)

newton_recursion_absolute_v1 (f, df, d2f, a, b, eps);
def newton_recursion_absolute_v2 (f, df, d2f, a, b, eps):	
	M2 = max([np.abs(d2f(x)) for x in [a, b]]) #M2 is the maximum value of |f''(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print (f"m1 = {m1}, M2 = {M2}")

	#Starting values
	if f(a) * d2f(a) > 0:
		x = a
	elif f(b) * d2f(b) > 0:
		x = b
	
	results = []
	results.append({
        'n': 0,
        'x_n': x,
        'f(x_n)': f(x),
        'delta_x=|x_n - x_(n-1)|': None
    })

	#Newton's method
	x_new = 0; new_eps = np.sqrt(eps * 2 * m1 / M2);
	print(f"delta_x = {new_eps}")

	i=0;
	while True:
		# calculate the next value of x
		x_new = x - f(x) / df(x)
		current_delta_x = abs(x-x_new)

		results.append({
            'n': i+1,
            'x_n': x_new,
            'f(x_n)': f(x_new),
            'delta_x=|x_n - x_(n-1)|': current_delta_x
        })
		
		# update the value of interval 
		x = x_new
		if (f(x_new) == 0): 
			break

		#stop condition
		if current_delta_x < new_eps:
			break
		else:
			i += 1

	#Print the final result
	df_results = pd.DataFrame(results)
	print(df_results.to_string(index=False))


	print(f"The value of root with absolute error {eps} is: {x}")
f = lambda x: x**4-27
df = lambda x: 4*x**3
d2f = lambda x: 12*x**2

a = 2
b = 3

eps = 0.5 * 10**(-7)

newton_recursion_absolute_v2 (f, df, d2f, a, b, eps);
def newton_recursion_relative (f, df, d2f, a, b, eta):	
	M2 = max([np.abs(d2f(x)) for x in [a, b]]) #M2 is the maximum value of |f''(x)| in the interval [a,b]
	m1 = min([np.abs(df(x)) for x in [a, b]]) #m1 is the minimum value of |f'(x)| in the interval [a,b]
	print (f"m1 = {m1}, M2 = {M2}")

	#Starting values
	if f(a) * d2f(a) > 0:
		x = a
	elif f(b) * d2f(b) > 0:
		x = b
	
	results = []
	results.append({
        'n': 0,
        'x_n': x,
        'f(x_n)': f(x),
        'sigma_x=|x_n-x_(n-1)|^2/|x_n|': None
    })

	#Newton's method
	x_new = 0; new_eta = 2 * eta * m1 / M2;
	print(f"sigma_x = {new_eta}")

	i=0;
	while True:
		# calculate the next value of x
		x_new = x - f(x) / df(x)
		current_sigma_x = (x-x_new)**2 / abs(x_new)

		results.append({
            'n': i+1,
            'x_n': x_new,
            'f(x_n)': f(x_new),
            'delta_x=|x_n - x_(n-1)|': current_sigma_x
        })
		
		# update the value of interval 
		x = x_new
		if (f(x_new) == 0): 
			break

		#stop condition
		if current_sigma_x < new_eta:
			break
		else:
			i += 1

	#Print the final result
	df_results = pd.DataFrame(results)
	print(df_results.to_string(index=False))


	print(f"The value of root with relative error {eta} is: {x}")
f = lambda x: x**4-27
df = lambda x: 4*x**3
d2f = lambda x: 12*x**2

a = 2
b = 3

eta = 0.5 * 10**(-7)

newton_recursion_relative (f, df, d2f, a, b, eta);




