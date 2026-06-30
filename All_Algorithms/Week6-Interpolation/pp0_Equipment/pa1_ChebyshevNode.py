# =============================================================================
# pa1_ChebyshevNode.py - Nút Chebyshev (Chebyshev Nodes)
#
# Chức năng: Tạo các điểm nội suy Chebyshev trên [a, b] để tránh
#            hiện tượng Runge khi nội suy đa thức bậc cao.
#
# Các hàm chính:
#   choose_interpolation_points(a, b, n) -> (DataFrame, x_values)
#
# Cách dùng: python pa1_ChebyshevNode.py
# =============================================================================
import numpy as np
import pandas as pd
from typing import Tuple

pd.set_option('display.precision', 7)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column
def choose_interpolation_points(a, b, n):

    # Step 1: Check interval condition
    condition = 1 if (a == -1 and b == 1) else 0

    # Step 2: Compute t_i
    i_values = np.arange(0, n + 1)
    t_values = np.cos(np.pi * (2 * i_values + 1) / (2 * n + 2))

    # Step 3: Compute x_i
    if condition == 1:
        x_values = t_values
    else:
        x_values = 0.5 * ((b - a) * t_values + (b + a))

    # Step 4: Display results in pandas table
    df = pd.DataFrame({
        "i": i_values,
        "t_i (cos)": t_values,
        "x_i (chosen point)": x_values
    })

    return df, x_values
# Example run
df, points = choose_interpolation_points(a=-1, b=1, n=9)

df.style.hide(axis="index")
df, points = choose_interpolation_points(a=3, b=6, n=12)

df.style.hide(axis="index")




