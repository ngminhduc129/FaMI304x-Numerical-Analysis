# =============================================================================
# pp1_fixed_point_iteration.py - Lặp điểm cố định cho hệ phi tuyến
#
# Chức năng: Tìm nghiệm của hệ phương trình phi tuyến dạng
#            x = g1(x,y,z), y = g2(x,y,z), z = g3(x,y,z), ...
#            bằng phương pháp lặp điểm cố định.
#
# Các hàm chính:
#   fixed_point_recursion_absolute(init, q, eps, *funcs, norm)
#   fixed_point_recursion_relative(init, q, eta, *funcs, norm)
#
# Cách dùng:
#   1. Định nghĩa các hàm g1, g2, g3, ...
#   2. Gọi fixed_point_recursion_absolute([x0,y0,z0], q, eps, g1,g2,g3)
#   3. Chạy: python pp1_fixed_point_iteration.py
# =============================================================================
import numpy as np
import pandas as pd

pd.set_option('display.precision', 15)  # Increase decimal precision
pd.set_option('display.width', 150)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def fixed_point_recursion_absolute (initial_values, q, eps, *funcs, norm):
    n = len(initial_values)
    values = np.array(initial_values)
    results = [[0] + initial_values]; i=0;

    new_eps = eps * (1-q)/ q
    print(f"new_eps: {new_eps}")

    while True:
        new_values = np.array([funcs[j](*values) for j in range(n)])
        
        # Calculate the differences
        diffs = np.abs(new_values - values)
        if (norm == 1): total_diff = np.sum(diffs)
        elif (norm == -1): total_diff = np.max(diffs)
        
        # Append the iteration number, new values, differences, and total_diff to the results list
        results.append([i+1] + new_values.tolist() + diffs.tolist() + [total_diff])
        values = new_values
        
        
        # Check for convergence
        if total_diff < new_eps:
            break
        i += 1  # Increment the iteration counter
        
    columns = ['Iteration'] + [f'x{j+1}' for j in range(n)] + [f'diff_x{j+1}' for j in range(n)] + ['total_diff']
    df = pd.DataFrame(results, columns=columns)
    print(df.to_string(index=False))  # Print the DataFrame without the index
    
    return values  
def g1(x, y, z):
    return 1/3 * (np.cos(y*z) + 0.5)

def g2(x, y, z):
    return 1/9 * np.sqrt(x**2 + np.sin(z) + 1.06) - 0.1

def g3(x, y, z):
    return -1/20 * np.exp(-x*y) - 9.1389/20

# Store g functions
def store_g_funcs(*g_funcs):
    return g_funcs

g_funcs = store_g_funcs(g1, g2, g3)

# Initial guess
initial_values = [0, 0, 0]
eps = 0.5 * 1e-7
q = 0.56 #column-norm

# Perform fixed-point iteration
solution = fixed_point_recursion_absolute(initial_values, q, eps, *g_funcs, norm=1)
print("Approximate solution:", solution)
def g1(x, y, z):
    return 1/3 * (np.cos(y*z)) + 1/6

def g2(x, y, z):
    return -1/9 * np.sqrt(x**2 + np.sin(z) + 1.06) - 0.1

def g3(x, y, z):
    return -1/20 * np.exp(-x*y) -(10*np.pi-3)/60

# Store g functions
def store_g_funcs(*g_funcs):
    return g_funcs

g_funcs = store_g_funcs(g1, g2, g3)

# Initial guess
initial_values = [0, 0, 0]
eps = 1e-6
q = 0.416 #column-norm

# Perform fixed-point iteration
solution = fixed_point_recursion_absolute(initial_values, q, eps, *g_funcs, norm=1)
print("Approximate solution:", solution)
def g1(x, y, z):
    return (13 - y**2 + 1.5 * z**2) / 15

def g2(x, y, z):
    return (11 + z - x**2) / 10

def g3(x, y, z):
    return (20 + y**2) / 30

# Store g functions
def store_g_funcs(*g_funcs):
    return g_funcs

g_funcs = store_g_funcs(g1, g2, g3)

# Initial guess
initial_values = [0, 0, 0]
eps = 1 * 1e-6
q = 0.425 #column-norm

# Perform fixed-point iteration
solution = fixed_point_recursion_absolute(initial_values, q, eps, *g_funcs, norm=1)
print("Approximate solution:", solution)
def fixed_point_recursion_relative (initial_values, q, eta, *funcs, norm):
    n = len(initial_values)
    values = np.array(initial_values)
    results = [[0] + initial_values]; i=0;

    new_eta = eta * (1-q)/ q
    print(f"new_eta: {new_eta}")

    while True:
        new_values = np.array([funcs[j](*values) for j in range(n)])
        
        # Calculate the differences
        diffs = np.abs(new_values - values)
        itself = np.abs(new_values)
        if (norm == 1): 
            total_diff = np.sum(diffs)
            total_itself = np.sum(itself)
        elif (norm == -1): 
            total_diff = np.max(diffs)
            total_itself = np.max(itself)
       
        if (i!=0):
            relative_diff = total_diff / total_itself
        else:
            relative_diff = None
        
        # Append the iteration number, new values, differences, and total_diff to the results list
        results.append([i+1] + new_values.tolist() + diffs.tolist() + [total_diff] + [relative_diff])
        values = new_values
        
        # Check for convergence
        if i!=0 and relative_diff < new_eta:
            break    

        i += 1  # Increment the iteration counter 
    
    columns = ['Iteration'] + [f'x{j+1}' for j in range(n)] + [f'diff_x{j+1}' for j in range(n)] + ['total_diff'] + ['relative_diff']
    df = pd.DataFrame(results, columns=columns)
    print(df.to_string(index=False))  # Print the DataFrame without the index
    
    return values
def g1(x, y):
    return 1/3 * np.sin(x*y)

def g2(x, y):
    return 1/4 * np.cos(x**2+y**2)

# Store g functions
def store_g_funcs(*g_funcs):
    return g_funcs

g_funcs = store_g_funcs(g1, g2)

# Initial guess
initial_values = [1, 1]
eta = 1 * 1e-4
q = 1/3 * np.cos(1) + 1/2 * np.sin(2) #column-norm

# Perform fixed-point iteration
solution = fixed_point_recursion_relative (initial_values, q, eta, *g_funcs, norm=1)
print("Approximate solution:", solution)




