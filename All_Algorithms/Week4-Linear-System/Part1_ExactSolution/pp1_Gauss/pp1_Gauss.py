# =============================================================================
# pp1_Gauss.py - Phương pháp khử Gauss (Gaussian Elimination)
#
# Chức năng: Giải hệ phương trình tuyến tính Ax = B bằng phương
#            pháp khử Gauss (forward elimination + back substitution).
#
# Các hàm chính:
#   gauss(A, B)       - đưa ma trận về dạng bậc thang
#   gauss_solve(A, B) - giải hệ hoàn chỉnh
#
# Input: Đọc từ file G_input_A.txt và G_input_B.txt
# Cách dùng: python pp1_Gauss.py
# =============================================================================
from pathlib import Path
from sympy import Matrix, zeros, N
import numpy as np
import contextlib

__dir__ = Path(__file__).parent.resolve()

# Hàm này được thêm vào để loại bỏ các giá trị nhỏ hơn 10⁻¹⁰, đảm bảo chúng được coi là 0.
def chop_small_values(mat):
    if isinstance(mat, Matrix):
        return mat.applyfunc(lambda x: 0 if abs(x) < 1e-10 else x)
    else:
        return 0 if abs(mat) < 1e-10 else mat
    
#Danh sách các vị trí khác 0 đầu tiên của từng hàng, dùng để sắp xếp AB theo cấu trúc bậc thang
def indexpivot(AB):
    rowAB, colAB = AB.shape
    index = []
    for row in range(rowAB):
        for col in range(colAB):
            if abs(AB[row, col]) > 1e-10: 
                index.append(col)
                break
            if col == colAB - 1:
                index.append(col + 1)
    return index
    
#Vị trí các cột có phần tử ở đường chéo chính bằng 0
def missingindex(A):
    colA = A.shape[1]
    index = indexpivot(A)
    return list(set(range(colA)) - set(index))
    
#Khử gauss
def gauss(A, B):
    A, B = Matrix(N(A)), Matrix(N(B))
    AB = chop_small_values(A.row_join(B))
    rowAB, colAB = AB.shape
    
    #Sắp xếp lại AB theo cấu trúc bậc thang trong trường hợp input không phải là bậc thang    
    index = indexpivot(AB)
    sorted_order = sorted(range(len(index)), key=lambda k: index[k])
    index = sorted(index)
    AB = chop_small_values(AB[list(sorted_order), :])
    display(AB)

    #Khử Gauss (Vừa khử, vừa sắp xếp lại AB theo cấu trúc bậc thang)
    for i in range(len(index)):
        if index[i] != colAB:
            for j in range(i + 1, len(index)):
                AB[j, :] = chop_small_values(Matrix(AB[j, :] - AB[i, :] / AB[i, index[i]] * AB[j, index[i]]))
            index = indexpivot(AB)
            sorted_order = sorted(range(len(index)), key=lambda k: index[k])
            index = sorted(index)
            AB = chop_small_values(AB[list(sorted_order), :])
            print("Bước khử")
            display(AB)
    print("Ma trận sau khi khử:")
    display(AB)
    return AB

#Giải AX=B
def gauss_solve(A, B):
    #Ma trận đã khử AB
    A, B = Matrix(N(A)), Matrix(N(B))
    AB = gauss(A, B)

    #Kiểm tra trường hợp vô nghiệm
    for i in indexpivot(AB):
        if i >= A.shape[1] and i < AB.shape[1]:
            print("Hệ vô nghiệm")
            return 

    #Ma trận A, B mới từ ma trận đã khử    
    A = AB[:, :A.shape[1]]
    B = AB[:, A.shape[1]:]
    colB = B.shape[1]

    #Chuyển những cột có phần tử trên đường chéo chính bằng 0 sang B
    missindex = missingindex(A)
    if missindex:
        print(f"index: {[x + 1 for x in indexpivot(A)]}")
        print(f"missing index: {[x + 1 for x in missindex]}")
        for i in missindex:
            move_col = -1 * A[:, i]
            B = B.row_join(move_col)
        A = A[:, [i for i in range(A.shape[1]) if i not in missindex]]
        print("Ma trận sau khi chuyển cột:")
        display(chop_small_values(A.row_join(B)))

    #Tính toán ma trận X    
    X = zeros(A.shape[1], B.shape[1])
    rowX, colX = X.shape
    index = indexpivot(A)
    
    for col in range(colX):
        for row in range(rowX - 1, -1, -1):
            sum_val = sum(A[row, k] * X[k, col] for k in range(index[row], rowX))
            X[row, col] = chop_small_values((B[row, col] - sum_val) / A[row, index[row]])
    
    print(f"Số phần tử tự do là {len(missindex)}: X{[x + 1 for x in missindex]}. Nghiệm:")
    return X
if __name__ == "__main__":
    output_path = str(__dir__ / "pp1_Gauss_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        A = Matrix(np.loadtxt(__dir__ / "G_input_A.txt"))
        print("A =")
        display(A)
        print("B =")
        B = Matrix(np.loadtxt(__dir__ / "G_input_B.txt"))
        display(B)
        print("A|B=")
        X = gauss_solve(A,B)
        display(X)
    print(f"Đã ghi kết quả vào {output_path}")




