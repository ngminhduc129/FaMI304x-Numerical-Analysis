# =============================================================================
# pp2_RungeKetta.py - Phương pháp Runge-Kutta (RK2, RK3, RK4)
#
# Chức năng: Giải hệ ODE bằng phương pháp Runge-Kutta với bậc tùy
#            chọn (2, 3, 4). Hỗ trợ tạo bảng Butcher tùy chỉnh.
#
# Các hàm chính:
#   rk_system_solver(F, x0, Y0, X_end, N, rk_order=4)
#   generate_rk2_tableau(alpha2) / generate_rk3_tableau(...)
#   generate_robust_rk4(alpha2, alpha3, alpha4)
#   get_coefficients(order, nodes=None)
#   plot_solver_result(df, mode, col1, col2)
#
# Cách dùng: python pp2_RungeKetta.py
# =============================================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Set display options for high precision as requested
pd.set_option('display.precision', 15)
pd.set_option('display.width', 150)
pd.set_option('display.max_columns', None)
def generate_rk2_tableau(alpha2):
    """
    Derives the RK2 Butcher Tableau weights (r) and matrix (beta) 
    given the custom node alpha2.
    Based on Order Conditions from '26_Các phương pháp Runge – Kutta 20251.pdf' (Page 10).
    """
    print(f"--- Generating RK2 Tableau for node: alpha2={alpha2} ---")
    
    if alpha2 == 0:
        print("Error: alpha2 cannot be 0 for a 2nd-order method (division by zero in derivation).")
        return None

    # 1. Solve for Weights (r1, r2)
    # Equation 2: r2 * alpha2 = 1/2  => r2 = 1 / (2 * alpha2)
    r2 = 1.0 / (2 * alpha2)
    
    # Equation 1: r1 + r2 = 1
    r1 = 1.0 - r2

    # 2. Solve for Matrix Coefficients (beta)
    # Equation 3: r2 * beta11 = 1/2
    # Since r2 = 1/(2*alpha2), substituting gives: (1/(2*alpha2)) * beta11 = 1/2 => beta11 = alpha2
    beta11 = alpha2
    
    # 3. Construct the Butcher Tableau
    #   0  |
    #   a2 | b11
    #   ----------------
    #      | r1   r2
    
    tableau = {
        'nodes': np.array([0, alpha2]),
        'weights': np.array([r1, r2]),
        'matrix': np.array([
            [0, 0],       # Row 1 (k1)
            [beta11, 0]   # Row 2 (k2)
        ])
    }
    
    # Print Tableau
    print("\n[Butcher Tableau Result]")
    print(f"{0:<5} |")
    print(f"{alpha2:<5.4f} | {beta11:.4f}")
    print("-" * 20)
    print(f"      | {r1:.4f}  {r2:.4f}")
    
    return tableau

def print_rk2_algorithm(tableau):
    """
    Prints the raw mathematical algorithm based on the generated tableau.
    """
    if tableau is None: return

    r = tableau['weights']
    b = tableau['matrix']
    a = tableau['nodes']
    
    print("\n" + "="*30)
    print(" RAW MARKDOWN ALGORITHM ")
    print("="*30)
    
    md = f"""
### Generalized Runge-Kutta 2 Algorithm (Derived)

**Input:** $x_n, y_n, h, f(x,y)$

**Step 1: Calculate Slopes ($k_i$)**
* $k_1 = h \\cdot f(x_n, y_n)$
* $k_2 = h \\cdot f(x_n + {a[1]:.4g}h, y_n + {b[1,0]:.4g}k_1)$

**Step 2: Update Solution ($y_{{n+1}}$)**
$$y_{{n+1}} = y_n + {r[0]:.4g}k_1 + {r[1]:.4g}k_2$$

**Output:** $y_{{n+1}}$
"""
    print(md)
# --- EXECUTION ---

# Example 1: Heun's Method (Trapezoidal) -> alpha2 = 1 [cite: 47]
tableau_heun = generate_rk2_tableau(alpha2=1)
print_rk2_algorithm(tableau_heun)

# Example 2: Midpoint Method -> alpha2 = 1/2 [cite: 48]
#tableau_midpoint = generate_rk2_tableau(alpha2=0.5)
#print_rk2_algorithm(tableau_midpoint)

# Example 3: Ralston's Method (Minimum Error Bound) -> alpha2 = 2/3
#tableau_ralston = generate_rk2_tableau(alpha2=2/3)
#print_rk2_algorithm(tableau_ralston)
import numpy as np
import pandas as pd

pd.set_option('display.precision', 6)
pd.set_option('display.width', 150)

def generate_rk3_tableau(alpha2, alpha3):
    """
    Derives the RK3 Butcher Tableau weights (r) and matrix (beta) 
    given two custom nodes alpha2 and alpha3.
    Based on Order Conditions from '26_Các phương pháp Runge – Kutta 20251.pdf' (Page 13).
    """
    print(f"--- Generating RK3 Tableau for nodes: alpha2={alpha2}, alpha3={alpha3} ---")

    # 1. Solve for Weights (r1, r2, r3)
    # The system of linear equations for 3rd order accuracy:
    # Eq 1: r1 + r2 + r3 = 1
    # Eq 2: r2*a2 + r3*a3 = 1/2
    # Eq 3: r2*a2^2 + r3*a3^2 = 1/3  (derived from 1/2 * sum(r * a^2) = 1/6)
    
    A = np.array([
        [1, 1, 1],
        [0, alpha2, alpha3],
        [0, alpha2**2, alpha3**2]
    ], dtype=float)
    
    b = np.array([1, 0.5, 1/3], dtype=float)
    
    try:
        r = np.linalg.solve(A, b)
        r1, r2, r3 = r
    except np.linalg.LinAlgError:
        print("Error: The chosen nodes result in a singular matrix. Cannot derive a unique RK3 method.")
        return None

    # 2. Solve for Matrix Coefficients (beta)
    # Critical coupling equation (Page 13): r3 * beta22 * alpha2 = 1/6
    if r3 == 0 or alpha2 == 0:
        print("Error: r3 or alpha2 is zero, cannot solve coupling equation r3 * beta22 * alpha2 = 1/6.")
        return None
        
    beta22 = 1 / (6 * r3 * alpha2)
    
    # Use Row Sum assumption (alpha_i = sum(beta_ij))
    # Row 2 (index 1): beta11 = alpha2 (since explicit, only one prev term)
    beta11 = alpha2 
    
    # Row 3 (index 2): beta21 + beta22 = alpha3
    beta21 = alpha3 - beta22
    
    # 3. Construct the Butcher Tableau
    #   0  |
    #   a2 | b11
    #   a3 | b21  b22
    #   ----------------
    #      | r1   r2   r3
    
    tableau = {
        'nodes': np.array([0, alpha2, alpha3]),
        'weights': r,
        'matrix': np.array([
            [0, 0, 0],       # Row 1 (k1)
            [beta11, 0, 0],  # Row 2 (k2)
            [beta21, beta22, 0] # Row 3 (k3)
        ])
    }
    
    # Print Tableau
    print("\n[Butcher Tableau Result]")
    print(f"{0:<5} |")
    print(f"{alpha2:<5.4f} | {beta11:.4f}")
    print(f"{alpha3:<5.4f} | {beta21:.4f}  {beta22:.4f}")
    print("-" * 20)
    print(f"      | {r1:.4f}  {r2:.4f}  {r3:.4f}")
    
    return tableau

def print_rk3_algorithm(tableau):
    """
    Prints the raw mathematical algorithm based on the generated tableau.
    """
    if tableau is None: return

    r = tableau['weights']
    b = tableau['matrix']
    a = tableau['nodes']
    
    print("\n" + "="*30)
    print(" RAW MARKDOWN ALGORITHM ")
    print("="*30)
    
    md = f"""
### Generalized Runge-Kutta 3 Algorithm (Derived)

**Input:** $x_n, y_n, h, f(x,y)$

**Step 1: Calculate Slopes ($k_i$)**
* $k_1 = h \\cdot f(x_n, y_n)$
* $k_2 = h \\cdot f(x_n + {a[1]:.4g}h, y_n + {b[1,0]:.4g}k_1)$
* $k_3 = h \\cdot f(x_n + {a[2]:.4g}h, y_n + {b[2,0]:.4g}k_1 + {b[2,1]:.4g}k_2)$

**Step 2: Update Solution ($y_{{n+1}}$)**
$$y_{{n+1}} = y_n + {r[0]:.4g}k_1 + {r[1]:.4g}k_2 + {r[2]:.4g}k_3$$

**Output:** $y_{{n+1}}$
"""
    print(md)
# --- EXECUTION ---
# Example 1: Heun's RK3 (Nodes: 1/3, 2/3) - As seen on Page 15
tableau_heun = generate_rk3_tableau(alpha2=1/3, alpha3=2/3)
print_rk3_algorithm(tableau_heun)

# Example 2: Classic RK3/Simpson (Nodes: 1/2, 1) - As seen on Page 14
#tableau_classic = generate_rk3_tableau(alpha2=1/2, alpha3=1)
#print_rk3_algorithm(tableau_classic)
import numpy as np
import pandas as pd

pd.set_option('display.precision', 6)

def generate_robust_rk4(alpha2, alpha3, alpha4):
    """
    Derives RK4 Butcher Tableau for ANY nodes, including the singular 'Classic' case
    where alpha2 == alpha3.
    """
    print(f"--- Generating RK4 for nodes: [{alpha2}, {alpha3}, {alpha4}] ---")
    
    # --- STEP 1: Solve for Weights (r) ---
    # Standard System: sum(r*a^k) = 1/(k+1)
    # Check for singularity (Duplicate nodes)
    if abs(alpha2 - alpha3) < 1e-9:
        # SINGULAR CASE (e.g., Classic RK4: 0.5, 0.5, 1.0)
        # We assume symmetry: r2 = r3
        # Reduced system for (r1, r_sum, r4) where r_sum = r2 + r3
        
        # Eqs:
        # 1) r1 + r_sum + r4 = 1
        # 2) r_sum * a2 + r4 * a4 = 1/2
        # 3) r_sum * a2^2 + r4 * a4^2 = 1/3
        
        A_reduced = np.array([
            [1, 1, 1],
            [0, alpha2, alpha4],
            [0, alpha2**2, alpha4**2]
        ], dtype=float)
        b_reduced = np.array([1, 0.5, 1/3], dtype=float)
        
        try:
            res = np.linalg.solve(A_reduced, b_reduced)
            r1 = res[0]
            r4 = res[2]
            r2 = res[1] / 2.0  # Split r_sum equally
            r3 = res[1] / 2.0
            r = np.array([r1, r2, r3, r4])
        except np.linalg.LinAlgError:
             print("Error: Reduced system is still singular.")
             return None
             
    else:
        # STANDARD CASE (Distinct nodes, e.g., 3/8 Rule)
        A = np.array([
            [1, 1, 1, 1],
            [0, alpha2, alpha3, alpha4],
            [0, alpha2**2, alpha3**2, alpha4**2],
            [0, alpha2**3, alpha3**3, alpha4**3]
        ], dtype=float)
        b = np.array([1, 0.5, 1/3, 1/4], dtype=float)
    
        try:
            r = np.linalg.solve(A, b)
            r1, r2, r3, r4 = r
        except np.linalg.LinAlgError:
            print("Error: Nodes result in a singular matrix.")
            return None

    # --- STEP 2: Solve for Matrix Coefficients (beta) ---
    # We solve the coupling equations.
    # Analytic Denominator: r4 * a3 * (a3 - a2)
    
    numerator = (1/12) - (alpha2 * 1/6)
    denominator = r4 * alpha3 * (alpha3 - alpha2)

    if abs(denominator) < 1e-9:
        # Denominator is zero. Check numerator.
        if abs(numerator) < 1e-9:
            # 0/0 Case (Classic RK4)
            # Standard choice for Classic RK4 is beta43 = 1.0
            beta43 = 1.0  
            
            # Solve Eq 3: r4 * b43 * b32 * a2 = 1/24
            # b32 = 1 / (24 * r4 * b43 * a2)
            if abs(24 * r4 * beta43 * alpha2) < 1e-9:
                 print("Error: Cannot solve for beta32 (div by zero)")
                 return None
            beta32 = 1.0 / (24 * r4 * beta43 * alpha2)
            
            # Solve Eq 1 for beta42: r4*(b42*a2 + b43*a3) + ... = 1/6
            # For Classic RK4 (a2=a3=0.5), usually beta42 = 0
            # Let's solve specifically:
            # 1/6 = r3*b32*a2 + r4*b42*a2 + r4*b43*a3
            term1 = r3 * beta32 * alpha2
            term3 = r4 * beta43 * alpha3
            # r4 * b42 * a2 = 1/6 - term1 - term3
            rhs = (1/6) - term1 - term3
            beta42 = rhs / (r4 * alpha2)
            
        else:
            print("Error: Inconsistent nodes (Division by zero in coupling equations).")
            return None
    else:
        # General Case
        beta43 = numerator / denominator
        beta32 = 1.0 / (24 * r4 * beta43 * alpha2)
        rhs = (1/6) - (r3 * alpha2 * beta32) - (r4 * alpha3 * beta43)
        beta42 = rhs / (r4 * alpha2)

    # Fill remaining betas via Row Sums
    beta21 = alpha2
    beta31 = alpha3 - beta32
    beta41 = alpha4 - beta42 - beta43

    # --- STEP 3: Construct Tableau ---
    tableau = {
        'nodes': np.array([0, alpha2, alpha3, alpha4]),
        'weights': r,
        'matrix': np.array([
            [0, 0, 0, 0],
            [beta21, 0, 0, 0],
            [beta31, beta32, 0, 0],
            [beta41, beta42, beta43, 0]
        ])
    }
    
    # Print Tableau
    print("\n[Butcher Tableau Result]")
    print(f"{0:<6} |")
    print(f"{alpha2:<6.4f} | {beta21:.4f}")
    print(f"{alpha3:<6.4f} | {beta31:.4f}  {beta32:.4f}")
    print(f"{alpha4:<6.4f} | {beta41:.4f}  {beta42:.4f}  {beta43:.4f}")
    print("-" * 30)
    print(f"       | {r1:.4f}  {r2:.4f}  {r3:.4f}  {r4:.4f}")

    return tableau

def print_rk4_algorithm(tableau):
    if tableau is None: return

    r = tableau['weights']
    b = tableau['matrix']
    a = tableau['nodes']

    # Helper to clean up 0.0 and 1.0
    def fmt(val):
        if abs(val) < 1e-9: return "0"
        if abs(val - 1.0) < 1e-9: return "1"
        if abs(val - 0.5) < 1e-9: return "0.5" # Make 0.5 explicit
        return f"{val:.4g}"

    print("\n" + "="*30)
    print(" RAW MARKDOWN ALGORITHM ")
    print("="*30)
    
    md = f"""
### Generalized Runge-Kutta 4 Algorithm (Derived)

**Input:** $x_n, y_n, h, f(x,y)$

**Step 1: Calculate Slopes ($k_i$)**
* $k_1 = h \\cdot f(x_n, y_n)$
* $k_2 = h \\cdot f(x_n + {fmt(a[1])}h, y_n + {fmt(b[1,0])}k_1)$
* $k_3 = h \\cdot f(x_n + {fmt(a[2])}h, y_n + {fmt(b[2,0])}k_1 + {fmt(b[2,1])}k_2)$
* $k_4 = h \\cdot f(x_n + {fmt(a[3])}h, y_n + {fmt(b[3,0])}k_1 + {fmt(b[3,1])}k_2 + {fmt(b[3,2])}k_3)$

**Step 2: Update Solution ($y_{{n+1}}$)**
$$y_{{n+1}} = y_n + {fmt(r[0])}k_1 + {fmt(r[1])}k_2 + {fmt(r[2])}k_3 + {fmt(r[3])}k_4$$

**Output:** $y_{{n+1}}$
"""
    print(md)


# --- EXECUTION ---
# Trying your specific node selection:
tableau_classic = generate_robust_rk4(0.5, 0.5, 1.0)
print_rk4_algorithm(tableau_classic)
# Case 2: RK4 3/8 Rule (Page 17 of your file)
# Nodes: 0, 1/3, 2/3, 1
tableau_38 = generate_robust_rk4(alpha2=1/3, alpha3=2/3, alpha4=1.0)
print_rk4_algorithm(tableau_38)
def get_coefficients(order, nodes=None):
    """
    Generates weights (r) and matrix (beta) for explicit RK methods.
    Returns: nodes (alpha), matrix (beta), weights (r)
    """
    # Defaults if no nodes provided
    if nodes is None:
        if order == 2: nodes = [0.5] # Midpoint
        elif order == 3: nodes = [0.5, 1.0] # Classic RK3
        elif order == 4: nodes = [0.5, 0.5, 1.0] # Classic RK4

    # Always prepend 0 for the first node if not present (explicit methods start at x_n)
    full_nodes = np.array([0] + list(nodes)) if nodes[0] != 0 else np.array(nodes)
    
    s = len(full_nodes)
    r = np.zeros(s)
    beta = np.zeros((s, s))

    # --- RK2 ---
    if order == 2:
        alpha2 = full_nodes[1]
        r[1] = 1.0 / (2 * alpha2)
        r[0] = 1.0 - r[1]
        beta[1, 0] = alpha2

    # --- RK3 ---
    elif order == 3:
        a2, a3 = full_nodes[1], full_nodes[2]
        # Solve weights
        A = np.array([[1,1,1], [0, a2, a3], [0, a2**2, a3**2]])
        b = np.array([1, 0.5, 1/3])
        r = np.linalg.solve(A, b)
        # Solve Matrix
        beta[1, 0] = a2
        beta[2, 1] = 1 / (6 * r[2] * a2)
        beta[2, 0] = a3 - beta[2, 1]

    # --- RK4 ---
    elif order == 4:
        a2, a3, a4 = full_nodes[1], full_nodes[2], full_nodes[3]
        
        # Check for Singular Case (Classic RK4 where a2 == a3)
        if abs(a2 - a3) < 1e-9:
            # Assume Symmetry r2=r3
            A_red = np.array([[1,1,1], [0, a2, a4], [0, a2**2, a4**2]])
            b_red = np.array([1, 0.5, 1/3])
            res = np.linalg.solve(A_red, b_red)
            r[0], r[3] = res[0], res[2]
            r[1] = res[1] / 2.0
            r[2] = res[1] / 2.0
            
            # Classic Coupling
            beta[3, 2] = 1.0
            beta[2, 1] = 1.0 / (24 * r[3] * beta[3, 2] * a2) # beta32
            # beta42 derived from row sum or coupling equation
            beta[3, 1] = (1/6 - r[2]*a2*beta[2,1] - r[3]*a3*beta[3,2]) / (r[3]*a2)
        else:
            # Standard Solver
            A = np.array([
                [1, 1, 1, 1],
                [0, a2, a3, a4],
                [0, a2**2, a3**2, a4**2],
                [0, a2**3, a3**3, a4**3]
            ])
            b = np.array([1, 0.5, 1/3, 1/4])
            r = np.linalg.solve(A, b)
            
            # Coupling
            num = (1/12) - (a2 * 1/6)
            den = r[3] * a3 * (a3 - a2)
            beta[3, 2] = num / den # beta43
            beta[2, 1] = 1.0 / (24 * r[3] * beta[3, 2] * a2) # beta32
            beta[3, 1] = (1/6 - r[2]*a2*beta[2,1] - r[3]*a3*beta[3,2]) / (r[3]*a2) # beta42

        # Fill Row Sums
        beta[1, 0] = a2
        beta[2, 0] = a3 - beta[2, 1]
        beta[3, 0] = a4 - beta[3, 1] - beta[3, 2]

    return full_nodes, beta, r
def rk_system_solver(F, x0, Y0, X_end, N, rk_order=4, custom_nodes=None):
    """
    Solves a system of ODEs using a generalized Runge-Kutta method.
    
    Parameters:
    - F: Vector function F(x, Y) returning a numpy array.
    - x0, Y0: Initial conditions.
    - X_end: End of interval.
    - N: Number of steps.
    - rk_order: Order of the method (2, 3, or 4).
    - custom_nodes: List of nodes [alpha2, alpha3...] (optional).
    """
    
    # 1. Setup Step Size and Coefficients
    h = (X_end - x0) / N
    nodes, beta, r = get_coefficients(rk_order, custom_nodes)
    
    # Print Method Info
    print(f"--- Solving with RK{rk_order} ---")
    print(f"Nodes: {nodes}")
    print(f"Weights (r): {r}")
    print(f"Step size h: {h}\n")
    
    # 2. Initialize
    x = x0
    Y = np.array(Y0, dtype=float)
    num_vars = len(Y)
    
    # Storage for output
    results = []
    
    # Helper to store row
    def record_step(n, x_val, Y_val):
        row = {'n': n, 'x': x_val}
        for i in range(num_vars):
            row[f'y{i}'] = Y_val[i]
        results.append(row)

    # Record Initial State
    record_step(0, x, Y)
    
    # 3. Iteration Loop
    for n in range(N):
        # Calculate Slopes (k1...ks)
        k = np.zeros((len(nodes), num_vars))
        
        # k1 (Always standard)
        k[0] = h * F(x, Y)
        
        # k2 to ks
        for i in range(1, len(nodes)):
            # Calculate sum of previous k contributions
            k_sum = np.zeros(num_vars)
            for j in range(i):
                if beta[i, j] != 0:
                    k_sum += beta[i, j] * k[j]
            
            # Evaluate F at intermediate point
            k[i] = h * F(x + nodes[i]*h, Y + k_sum)
            
        # Update Y using weights (r)
        delta_Y = np.zeros(num_vars)
        for i in range(len(nodes)):
            delta_Y += r[i] * k[i]
            
        Y = Y + delta_Y
        x = x + h
        
        # Record Step
        record_step(n+1, x, Y)
        
    # 4. Final Output
    df_results = pd.DataFrame(results)
    return df_results
def circuit(x, Y):
    Q = Y[0]
    I = Y[1]
    
    # Arbitrary constants for demonstration
    R = 35
    L = 2.5
    C = 0.26
    E = 50*np.cos(5*x)

    dQ_dt = I
    dI_dt = ((E - R*I - (Q*1.0)/C) * 1.0) / L
    
    return np.array([dQ_dt, dI_dt])

# Define Parameters
x0 = 0
Y0 = [18, 15] # Initial population: 20 Prey, 5 Predators
X_end = 5
N_steps = 500
df_rk4 = rk_system_solver(circuit, x0, Y0, X_end, N_steps, rk_order=4, custom_nodes=[0.5, 0.5, 1.0])
df_rk4
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
plot_solver_result(df_rk4, mode=1, col1='x', col2=['y0','y1'])
plot_solver_result(df_rk4, mode=2, col1='y0', col2=['y1'])




