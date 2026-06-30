# =============================================================================
# pp4_Choleski.py - Phân rã Cholesky (A = L*L^T)
#
# Chức năng: Phân tích ma trận đối xứng, xác định dương A thành
#            tích L*L^T, dùng để giải hệ Ax = b.
#
# Các hàm chính:
#   cholesky_decomposition(A) - phân rã Cholesky
#   forward_substitution(L,b) - thế xuôi
#   backward_substitution(U,y)- thế ngược
#   solve_cholesky(A, b)      - giải hệ hoàn chỉnh
#
# Cách dùng: python pp4_Choleski.py
# =============================================================================
import numpy as np

def read_matrix_from_file(filename):
    """
    Read a matrix from a text file and display it.
    
    Parameters:
    filename (str): Path to the text file containing the matrix
    
    Returns:
    numpy.ndarray: The matrix read from the file
    """
    try:
        # Read the matrix from the file
        matrix = np.loadtxt(filename)
        
        # Display the matrix
        print("Matrix read from file:")
        print(matrix)
        
        return matrix
    
    except FileNotFoundError:
        print(f"Error: File '{filename}' not found.")
        return None
    except Exception as e:
        print(f"Error reading file: {str(e)}")
        return None
def cholesky_decomposition(A):
    """
    Perform Cholesky decomposition on matrix A.
    
    Parameters:
    A (numpy.ndarray): Input symmetric positive-definite matrix
    
    Returns:
    numpy.ndarray: Lower triangular matrix L where A = L × L^T
    """
    n = len(A)
    L = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i + 1):
            if i == j:
                # Diagonal elements
                sum_sq = sum(L[i][k] ** 2 for k in range(j))
                L[i][j] = np.sqrt(A[i][i] - sum_sq)
            else:
                # Non-diagonal elements
                sum_mult = sum(L[i][k] * L[j][k] for k in range(j))
                L[i][j] = (A[i][j] - sum_mult) / L[j][j]
    
    return L

def verify_cholesky_decomposition(A, L):
    """
    Verify Cholesky decomposition by multiplying L and L^T and comparing with original matrix.
    
    Parameters:
    A (numpy.ndarray): Original matrix
    L (numpy.ndarray): Lower triangular matrix
    """
    print("\nLower triangular matrix L:")
    print(L)
    print("\nUpper triangular matrix L^T:")
    print(L.T)
    
    # Multiply L and L^T
    LLT = np.dot(L, L.T)
    print("\nL × L^T:")
    print(LLT)
    
    print("\nOriginal matrix A:")
    print(A)
    
    # Check if L×L^T equals A
    if np.allclose(LLT, A):
        print("\nVerification successful: L×L^T = A")
    else:
        print("\nVerification failed: L×L^T ≠ A")
# Example usage
A = read_matrix_from_file('CLSK_input_A.txt')
if A is not None:
    # Check if matrix is symmetric
    if not np.allclose(A, A.T):
        print("Error: Matrix must be symmetric for Cholesky decomposition")
    else:
        L = cholesky_decomposition(A)
        verify_cholesky_decomposition(A, L)
def forward_substitution(L, b):
    n = len(L)
    y = np.zeros(n)
    for i in range(n):
        s = sum(L[i][j] * y[j] for j in range(i))
        y[i] = (b[i] - s) / L[i][i]
    return y

def backward_substitution(U, y):
    n = len(U)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        s = sum(U[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (y[i] - s) / U[i][i]
    return x

def solve_cholesky(A, b):
    # Bước 1: Phân rã A = L*L^T
    L = cholesky_decomposition(A)
    
    # Bước 2: Giải Ly = b (Thế tiến)
    y = forward_substitution(L, b)
    
    # Bước 3: Giải L^T*x = y (Thế lùi)
    LT = L.T
    x = backward_substitution(LT, y)
    
    return x, L

# --- CHẠY THỬ VỚI MA TRẬN CỦA BẠN ---
b = np.array([8.0, 39.0])

x, L = solve_cholesky(A, b)

print("Ma trận L:\n", L)
print("\nNghiệm x của hệ phương trình:", x)

# Kiểm tra lại kết quả Ax có bằng b không
print("\nKiểm tra Ax:", np.dot(A, x))




