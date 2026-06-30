# =============================================================================
# pp1_Danielevsky.py - Phương pháp Danilevsky
#
# Chức năng: Tìm giá trị riêng và vector riêng của ma trận bằng cách
#            đưa ma trận về dạng Frobenius (dạng companion).
#
# Các hàm chính:
#   danilevsky_method(A, tol) -> (eigenvalues, eigenvectors)
#
# Input: Đọc từ file DVSK_input_A.txt
# Cách dùng: python pp1_Danielevsky.py
# =============================================================================
import numpy as np
from typing import Tuple, List
np.set_printoptions(linewidth=np.inf)
def danilevsky_method(A: np.ndarray, tol: float = 1e-10) -> Tuple[np.ndarray, List[np.ndarray]]:
    """
    Tìm trị riêng và vectơ riêng của ma trận A sử dụng phương pháp Danilevsky.
    
    Args:
        A: Ma trận vuông đầu vào
        tol: Ngưỡng sai số cho phép
        
    Returns:
        eigenvalues: Mảng các trị riêng
        eigenvectors: Danh sách các vectơ riêng tương ứng
    """
    n = len(A)
    if A.shape != (n, n):
        raise ValueError("Ma trận đầu vào phải là ma trận vuông")
    
    # Sao chép ma trận A để không thay đổi ma trận gốc
    P = A.copy()
    
    # Lưu trữ các ma trận biến đổi
    S_list = []
    
    # Thực hiện các phép biến đổi Danilevsky
    for k in range(n-1, 0, -1):
        print(f"\nBước k = {k}:")
        print("Ma trận hiện tại:")
        print(P)
        
        # Kiểm tra các trường hợp
        if abs(P[k,k-1]) >= tol:  # TH1: a_{n,n-1} ≠ 0
            print(f"\nTH1: a_{k+1},{k} ≠ 0")
            # Tạo ma trận biến đổi S
            S = np.eye(n)
            S[k-1] = -P[k] / P[k,k-1]
            S[k-1,k-1] = 1/P[k,k-1]
            
            print("\nMa trận biến đổi S:")
            print(S)
            
            # Lưu ma trận biến đổi
            S_list.append(S)
            
            # Cập nhật ma trận P
            P = np.linalg.inv(S) @ P @ S
            
            print("\nMa trận sau biến đổi:")
            print(P)
            
        else:  # TH2 hoặc TH3
            # Tìm phần tử khác 0 đầu tiên trong hàng k
            non_zero_idx = None
            for j in range(k-1, -1, -1):
                if abs(P[k,j]) >= tol:
                    non_zero_idx = j
                    break
            
            if non_zero_idx is not None:  # TH2: Tìm thấy phần tử khác 0
                print(f"\nTH2: a_{k+1},{k} = 0, a_{k+1},{non_zero_idx} ≠ 0")
                # Tạo ma trận hoán vị C
                C = np.eye(n)
                C[k-1,k-1] = 0
                C[non_zero_idx,non_zero_idx] = 0
                C[k-1,non_zero_idx] = 1
                C[non_zero_idx,k-1] = 1
                
                print(f"\nMa trận hoán vị C tại bước k={k}:")
                print(C)
                print(f"\nHoán vị hàng {k} và {non_zero_idx+1}")
                
                # Hoán vị hàng và cột
                P = C @ P @ C
                print("\nMa trận sau khi hoán vị:")
                print(P)
                
                S_list.append(C)
                i = k % (n-1) + 1
                # Tiếp tục với phép biến đổi thông thường
                S = np.eye(n)
                S[k-1] = -P[k] / P[k,k-1]
                S[k-1,k-1] = 1/P[k,k-1]
                
                print(f"\nMa trận nghịch đảo của M_{i}:")
                print(S)
                
                S_list.append(S)
                P = np.linalg.inv(S) @ P @ S
                
                print("\nMa trận sau biến đổi:")
                print(P)
                
            else:  # TH3: Tất cả phần tử đều bằng 0
                print(f"\nTH3: a_{k},j = 0 với mọi j ≤ {k-1}")
                print("Ma trận đã ở dạng Frobenius, không cần biến đổi")
                continue
    
    print("\nMa trận Frobenius cuối cùng:")
    print(P)
    
    # Lấy các hệ số của đa thức đặc trưng từ ma trận Frobenius
    coeffs = [1] + [-P[0,i] for i in range(n)]
    
    # In ra phương trình đặc trưng
    print("\nPhương trình đặc trưng P(λ) = 0:")
    poly_str = "λ^" + str(n)
    for i in range(1, n+1):
        if abs(coeffs[i]) >= tol:
            if coeffs[i] > 0:
                if abs(coeffs[i] - round(coeffs[i])) < tol:
                    poly_str += f" + {int(round(coeffs[i]))}λ^{n-i}"
                else:
                    poly_str += f" + {coeffs[i]:.15f}λ^{n-i}"
            else:
                if abs(coeffs[i] - round(coeffs[i])) < tol:
                    poly_str += f" - {abs(int(round(coeffs[i])))}λ^{n-i}"
                else:
                    poly_str += f" - {abs(coeffs[i]):.15f}λ^{n-i}"
    poly_str += " = 0"
    print(poly_str)
    
    # Tìm trị riêng từ ma trận Frobenius
    eigenvalues = np.roots(coeffs)
    
    # Tìm vectơ riêng
    eigenvectors = []
    for lambda_i in eigenvalues:
        # Tạo vectơ riêng cho ma trận Frobenius
        v = np.array([lambda_i**i for i in range(n)])
        
        # Biến đổi ngược để có vectơ riêng của ma trận ban đầu
        for S in reversed(S_list):
            v = S @ v
        
        # Chuẩn hóa vectơ riêng
        v = v / np.linalg.norm(v)
        eigenvectors.append(v)
    
    return eigenvalues, eigenvectors


def read_matrix_from_file(filename: str) -> np.ndarray:
    """
    Đọc ma trận từ file txt.
    
    Format file txt:
    - Mỗi hàng của ma trận trên một dòng
    - Các phần tử trong hàng cách nhau bởi dấu cách hoặc dấu phẩy
    
    Args:
        filename: Tên file txt chứa ma trận
        
    Returns:
        Ma trận numpy đọc được từ file
    """
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            
        # Loại bỏ các dòng trống và khoảng trắng thừa
        lines = [line.strip() for line in lines if line.strip()]
        
        # Đọc ma trận
        matrix = []
        for line in lines:
            # Hỗ trợ cả dấu cách và dấu phẩy làm dấu phân cách
            row = [float(x.strip()) for x in line.replace(',', ' ').split()]
            matrix.append(row)
            
        return np.array(matrix)
    except FileNotFoundError:
        print(f"Không tìm thấy file {filename}")
        return None
    except Exception as e:
        print(f"Lỗi khi đọc file: {str(e)}")
        return None

# Ví dụ đọc ma trận từ file

try:
    # Đọc ma trận từ file
    A_file = read_matrix_from_file("DVSK_input_A.txt")
    if A_file is not None:
        print("Ma trận đọc từ file:")
        print(A_file)
        
        # Tìm trị riêng và vectơ riêng
        eigenvalues_file, eigenvectors_file = danilevsky_method(A_file)
        
        print("\nTrị riêng:")
        for i, lambda_i in enumerate(eigenvalues_file):
            if abs(lambda_i.imag) < 1e-15:
                print(f"lambda_{i+1} = {lambda_i.real:.15f}")
            else:
                print(f"lambda_{i+1} = {lambda_i.real:.15f} + {lambda_i.imag:.15f}j")
                
        print("\nVectơ riêng:")
        for i, v in enumerate(eigenvectors_file):
            print(f"v_{i+1} = [")
            for j in range(len(v)):
                if abs(v[j].imag) < 1e-15:
                    print(f"    {v[j].real:.15f}")
                else:
                    print(f"    {v[j].real:.15f} + {v[j].imag:.15f}j")
            print("]")
except Exception as e:
    print(f"Lỗi khi xử lý ma trận: {str(e)}")




