# =============================================================================
# pp1c_HeunModified.py - Phương pháp Heun (Euler cải tiến)
#
# Chức năng: Giải ODE bằng phương pháp Heun (Improved Euler):
#            k1 = F(x,Y), k2 = F(x+h, Y+h*k1), Y += h*(k1+k2)/2.
#
# Các hàm chính:
#   improved_euler_system_solver(F, x0, Y0, X_end, N) -> DataFrame
#   plot_solver_result(df, mode, col1, col2)
#
# Cách dùng: python pp1c_HeunModified.py
# =============================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
def improved_euler_system_solver(F, x0, Y0, X_end, N):
    """
    Solves a system of ODEs using the Improved Euler (Heun's) Method.
    Formula: Y_{n+1} = Y_n + (h/2) * [F(x_n, Y_n) + F(x_{n+1}, Y*_n+1)]
    where Y*_n+1 is the predictor (standard Euler step).
    
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
    """
    
    # 1. Setup
    h = (X_end - x0) / N
    Y_curr = np.array(Y0, dtype=float)
    x_curr = x0
    
    # Storage
    results = []
    
    # Store Initial Condition
    row = {'Step': 0, 'x': x_curr}
    for i, val in enumerate(Y_curr):
        row[f'y_{i}'] = val
    results.append(row)
    
    print(f"Improved Euler Initialized.\n Range: [{x0}, {X_end}], Steps: {N}, h: {h}")
    print("Iterating...")

    # 2. Iteration Loop
    for j in range(N):
        # --- Step A: Predictor (Standard Euler) ---
        # Calculate slope at current point (K1)
        K1 = np.array(F(x_curr, Y_curr))
        
        # Predict next value (Y_star)
        Y_star = Y_curr + h * K1
        x_next = x_curr + h
        
        # --- Step B: Corrector (Trapezoidal Rule) ---
        # Calculate slope at the predicted point (K2)
        K2 = np.array(F(x_next, Y_star))
        
        # Apply the Improved Euler formula (Eq 6.43)
        # Y_{next} = Y_{curr} + (h/2) * (K1 + K2)
        Y_next = Y_curr + (h / 2.0) * (K1 + K2)
        
        # Update state
        Y_curr = Y_next
        x_curr = x_next
        
        # Store results
        row = {'Step': j + 1, 'x': x_curr}
        for i, val in enumerate(Y_curr):
            row[f'y_{i}'] = val
        results.append(row)

    return pd.DataFrame(results)
# --- Usage Example: System of Equations ---

def system_derivatives(x, Y):
    y0 = Y[0] # This is y
    y1 = Y[1] # This is z
    y2 = Y[2] # This is w
    
    dy0 = 0.4*y0*(1-y0/20.0) + 0.4*y1 - 0.3*y0*y2
    dy1 = 0.7*y1*(1-y1/25.0) - 0.4*y1 - 0.4*y1*y2
    dy2 = -0.3*y2 + 0.35*(y0+y1)*y2
    
    return [dy0, dy1, dy2]

# Parameters
x_start = 0
Y_start = [12, 18, 8]  # y(0)=12, z(0)=18, w(0)=8
x_stop = 1500
steps = 15000
# Run Solver
df_system = improved_euler_system_solver(system_derivatives, x_start, Y_start, x_stop, steps)

# Display Result
df_system
## Visualize
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
        plt.xlim(-10, 1500)       # Set x-axis range from 0 to 8
        #plt.ylim(-1.2, 1.2)  # Set y-axis range from -1.2 to 1.2

        plt.title(f"Mode 1: Time Series Plot ({col1} vs {col2})", fontsize=14)

    elif mode == 2:
        # Set the axis limits manually
        #plt.xlim(-10, 250)       # Set x-axis range from 0 to 8
        #plt.ylim(-1.2, 1.2)  # Set y-axis range from -1.2 to 1.2

        plt.title(f"Mode 2: Phase Plane Plot ({col1} vs {col2})", fontsize=14)
# 2. Plot Mode 1: x vs y_0 (Time vs Position)
plot_solver_result(df_system, mode=1, col1='x', col2=['y_0', 'y_1', 'y_2'])
# 3D Phase Plot: y_0 vs y_1 vs y_2
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(9, 6))
ax = fig.add_subplot(111, projection='3d')
# Plot trajectory in 3D phase space
ax.plot(df_system['y_0'], df_system['y_1'], df_system['y_2'], color='tab:blue', linewidth=2.0)
ax.set_xlabel('y_0')
ax.set_ylabel('y_1')
ax.set_zlabel('y_2')
ax.set_title('3D Phase Plot: y_0 vs y_1 vs y_2')
plt.tight_layout()
plt.show()
# --- Usage Example: ---

def highord_equation(t, Y):
    # Unpack state vector
    y0 = Y[0] # y
    y1 = Y[1] # y'
    y2 = Y[2] #y''

    z_t = y1
    dz_dt = y2
    d2z_dt2 = (1+t*z_t)*np.sin(1+y0*z_t)/(float(1+y0**2+dz_dt**2))

    return [z_t, dz_dt, d2z_dt2]

# Parameters
t_start = 0
t_stop = 10
steps = 200
# Initial conditions: y(0) = 1, y'(0) = 0.5, y''(0) = -1
init_conds = [1.0, 0.5, -1.0]
# Run
df_result = improved_euler_system_solver(highord_equation, t_start, init_conds, t_stop, steps)

# Display
print("\nOutput Table:")
df_result





