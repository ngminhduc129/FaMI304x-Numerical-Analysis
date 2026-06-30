# =============================================================================
# pp1a_EulerForward.py - Phương pháp Euler tiến (Forward Euler)
#
# Chức năng: Giải hệ phương trình vi phân thường (ODE) bằng phương
#            pháp Euler tiến: Y(n+1) = Y(n) + h*F(x(n), Y(n)).
#
# Các hàm chính:
#   euler_system_solver(F, x0, Y0, X_end, N) -> DataFrame
#   plot_solver_result(df, mode, col1, col2) - vẽ đồ thị
#
# Cách dùng: python pp1a_EulerForward.py
# =============================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def euler_system_solver(F, x0, Y0, X_end, N):
    """
    Solves a system of Cauchy differential problems using Euler-Forward method.
    
    Parameters:
    F : function
        A vector function F(x, Y) that returns a numpy array of derivatives.
    x0 : float
        Initial value of x.
    Y0 : list or np.array
        Initial vector of values [y1_0, y2_0, ...].
    X_end : float
        The end value of x.
    N : int
        Number of steps.
        
    Returns:
    pd.DataFrame
        Table containing steps, x, and all solution components y_k.
    """
    
    # 1. Setup Step Size
    h = (X_end - x0) / N
    
    # 2. Initialize Vector State
    # Convert input list to numpy array for vector math
    Y_curr = np.array(Y0, dtype=float) 
    x_curr = x0
    
    # Prepare storage for results
    # We will store: Step, x, and then all components of Y
    results = []
    
    print(f"Input Parameters:\n x0 = {x0}, Y0 = {Y_curr}\n X_end = {X_end}, N = {N}")
    print(f"Calculated Step Size h = {h}\n")
    print("Iterating System...")
    
    # Store initial state (Step 0)
    # We create a dictionary for the row to easily handle dynamic columns
    row = {'Step': 0, 'x': x_curr}
    for i, val in enumerate(Y_curr):
        row[f'y_{i}'] = val
    results.append(row)
    
    # 3. Iteration Loop
    for j in range(N):
        # Calculate Derivative Vector: F(x, Y)
        # This returns a vector [y1', y2', ...]
        dY = np.array(F(x_curr, Y_curr))
        
        # Update Vector State: Y_{j+1} = Y_j + h * F(x_j, Y_j)
        Y_next = Y_curr + h * dY
        x_next = x_curr + h
        
        # Update current values for next loop
        Y_curr = Y_next
        x_curr = x_next
        
        # Store results
        row = {'Step': j + 1, 'x': x_curr}
        for i, val in enumerate(Y_curr):
            row[f'y_{i}'] = val
        results.append(row)

    # 4. Create Output Table
    df_result = pd.DataFrame(results)
    return df_result
# --- Usage Example: System of Equations ---

def system_derivatives(x, Y):
    y0 = Y[0] # This is y
    y1 = Y[1] # This is z
    
    dy0 = 0.5*y0*(1-y0) - 0.15*y0*y1       
    dy1 = -0.3*y1 + 0.2*y0*y1     
    
    return [dy0, dy1]

# Parameters
x_start = 0
Y_start = [0.7, 0.5]  # y(0)=0.7, z(0)=0.5
x_stop = 2000
steps = 20000
 # Run Solver
df_system = euler_system_solver(system_derivatives, x_start, Y_start, x_stop, steps)

# Display Result
df_system
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
        plt.xlim(-10, 2000)       # Set x-axis range from 0 to 8
        #plt.ylim(-1.2, 1.2)  # Set y-axis range from -1.2 to 1.2

        plt.title(f"Mode 1: Time Series Plot ({col1} vs {col2})", fontsize=14)

    elif mode == 2:
        # Set the axis limits manually
        #plt.xlim(-10, 250)       # Set x-axis range from 0 to 8
        #plt.ylim(-1.2, 1.2)  # Set y-axis range from -1.2 to 1.2

        plt.title(f"Mode 2: Phase Plane Plot ({col1} vs {col2})", fontsize=14)
# 2. Plot Mode 1: x vs y_0 (Time vs Position)
plot_solver_result(df_system, mode=1, col1='x', col2=['y_0', 'y_1'])

# 3. Plot Mode 2: y_0 vs y_1 (Position vs Velocity)
plot_solver_result(df_system, mode=2, col1='y_0', col2=['y_1'])
# --- Usage Example: ---

def highord_equation(t, Y):
    # Unpack state vector
    y0 = Y[0]      # y
    y1 = Y[1]     # y'

    z_t =  y1
    dz_dt = (t+y0)*np.cos(1+z_t)
    return [z_t, dz_dt]

# Parameters
t_start = 0
t_stop = 20
steps = 200
# Initial conditions: y(0) = 1, y'(0) = -1
init_conds = [1.0, -1.0]
# Run
df_result = euler_system_solver(highord_equation, t_start, init_conds, t_stop, steps)

# Display
print("\nOutput Table:")
df_result




