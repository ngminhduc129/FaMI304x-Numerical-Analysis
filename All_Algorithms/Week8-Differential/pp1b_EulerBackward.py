# =============================================================================
# pp1b_EulerBackward.py - Phương pháp Euler lùi (Backward Euler)
#
# Chức năng: Giải ODE bằng phương pháp Euler lùi (ẩn), dùng lặp
#            điểm cố định nội tại để giải phương trình phi tuyến.
#
# Các hàm chính:
#   backward_euler_system_solver(F, x0, Y0, X_end, N, tol, max_iter)
#   plot_solver_result(df, mode, col1, col2)
#
# Cách dùng: python pp1b_EulerBackward.py
# =============================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set display options
pd.set_option('display.precision', 6)
pd.set_option('display.width', 1000)
def backward_euler_system_solver(F, x0, Y0, X_end, N, tol=1e-6, max_iter=50):
    """
    Solves a system Y' = F(x, Y) using the Backward Euler method.
    Formula: Y_{n+1} = Y_n + h * F(x_{n+1}, Y_{n+1})
    
    Since Y_{n+1} is on the right side, we use Fixed-Point Iteration 
    to solve for it at every step.
    
    Parameters:
    F : function
        Vector function F(x, Y) -> returns numpy array of derivatives
    x0 : float
        Initial x
    Y0 : list/array
        Initial state vector Y
    X_end : float
        Target x
    N : int
        Number of steps
    tol : float
        Convergence tolerance for the inner loop
    max_iter : int
        Maximum iterations for the inner loop
    """
    
    # 1. Setup
    h = (X_end - x0) / N
    Y_curr = np.array(Y0, dtype=float)
    x_curr = x0
    
    # Storage
    results = []
    
    # Store Initial Condition
    row = {'t': x_curr}
    for i, val in enumerate(Y_curr):
        row[f'y_{i}'] = val
    results.append(row)
    
    # print(f"Backward Euler Initialized. Range: [{x0}, {X_end}], h: {h}")
    
    # 2. Iteration Loop (Time Stepping)
    for j in range(N):
        x_next = x_curr + h
        
        # --- Implicit Solver (Inner Loop) ---
        
        # A. Predictor: Use Forward Euler to get a good initial guess
        # This speeds up convergence significantly compared to guessing Y_curr
        K_pred = np.array(F(x_curr, Y_curr))
        Y_guess = Y_curr + h * K_pred
        
        # B. Corrector: Fixed-Point Iteration
        # We iterate Y = Y_prev + h * F(x_next, Y) until Y stops changing
        for k in range(max_iter):
            # Calculate slope at the guessed future point
            K_implicit = np.array(F(x_next, Y_guess))
            
            # Apply Backward Euler Formula
            Y_new = Y_curr + h * K_implicit
            
            # Check for convergence (Matrix Norm)
            # We check the maximum difference across all vector components
            error = np.max(np.abs(Y_new - Y_guess))
            
            # Update guess
            Y_guess = Y_new
            
            if error < tol:
                break
        
        # C. Update State
        Y_curr = Y_guess
        x_curr = x_next
        
        # Store results
        row = {'t': x_curr}
        for i, val in enumerate(Y_curr):
            row[f'y_{i}'] = val
        results.append(row)

    return pd.DataFrame(results)
# --- Usage Example: Stiff System ---
# Backward Euler is famous for handling "Stiff" equations where other methods fail.
# System: 
# y' = -20y  (Decays very fast)
# z' = y + z (Coupled)

def stiff_system(t, Y):
    y = Y[0]
    z = Y[1]
    
    dy = -20 * y
    dz = y + z
    return [dy, dz]

# Parameters
t_start = 0
Y_start = [1.0, 0.0] # y=1, z=0
t_stop = 1.0
steps = 20 # Step size h=0.05. Forward Euler would be unstable here!
# Run
df_backward = backward_euler_system_solver(stiff_system, t_start, Y_start, t_stop, steps)
df_backward
def plot_solver_result(df, mode, col1='x', col2='y_0', output_file='plot.png'):
    """
    Visualizes the solver output.
    
    Parameters:
    df : pd.DataFrame
        The solution table.
    mode : int
        1 for Time Series (x vs y).
        2 for Phase Plot (y_i vs y_j).
    col1 : str
        Name of the column for the X-axis.
    col2 : str
        Name of the column for the Y-axis.
    output_file : str
        Filename to save the image.
    """
    plt.figure(figsize=(10, 6))
    plt.grid(True, linestyle='--', alpha=0.6)
    
    # Plot Logic
    for item in col2:
        plt.plot(df[col1], df[item], label=f'{item} vs {col1}', linewidth=2, )
        plt.xlabel(col1, fontsize=12)
        plt.ylabel(item, fontsize=12)
    # Formatting
    plt.legend()
    
    if mode == 1:
        # Set the axis limits manually
        #plt.xlim(-10, 2000)       # Set x-axis range from 0 to 8
        #plt.ylim(-1.2, 1.2)  # Set y-axis range from -1.2 to 1.2

        plt.title(f"Mode 1: Time Series Plot ({col1} vs {col2})", fontsize=14)

    elif mode == 2:
        # Set the axis limits manually
        #plt.xlim(-10, 250)       # Set x-axis range from 0 to 8
        #plt.ylim(-1.2, 1.2)  # Set y-axis range from -1.2 to 1.2

        plt.title(f"Mode 2: Phase Plane Plot ({col1} vs {col2})", fontsize=14)
plot_solver_result(df_backward, mode=1, col1='t', col2=['y_0', 'y_1'])





