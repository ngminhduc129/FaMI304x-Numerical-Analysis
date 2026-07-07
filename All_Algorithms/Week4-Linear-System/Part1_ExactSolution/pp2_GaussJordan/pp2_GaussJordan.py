# =============================================================================
# pp2_GaussJordan.py - Phương pháp Gauss-Jordan
#
# Chức năng: Giải hệ phương trình tuyến tính Ax = B bằng phương
#            pháp Gauss-Jordan (đưa về dạng ma trận đơn vị).
#
# Các hàm chính:
#   gauss_jordan(A, B)        - khử hoàn toàn về dạng rút gọn
#   gauss_jordan_solve(A, B)  - giải hệ kèm xử lý biến tự do
#
# Input: Đọc từ file GJ_input_A.txt và GJ_input_B.txt
# Cách dùng: python pp2_GaussJordan.py
# =============================================================================
from pathlib import Path
from sympy import Matrix, zeros, N
import numpy as np
import contextlib

__dir__ = Path(__file__).parent.resolve()

#Hàm này được thêm vào để loại bỏ các giá trị nhỏ hơn 10⁻¹⁰, đảm bảo chúng được coi là 0.
def chop_small_values(mat):
    if isinstance(mat, Matrix):
        return mat.applyfunc(lambda x: 0 if abs(x) < 1e-10 else x)
    else:
        return 0 if abs(mat) < 1e-10 else mat

#Hàm tìm kiếm phần tử khử, trả về giá trị và vị trí của phần tử khử
def findpivot(AB,rowA,colA,rowused,colused):
    #Hàng và cột chưa dùng
    row = list(set(range(rowA)) - set(rowused))
    col = list(set(range(colA)) - set(colused))
    
    #Ưu tiên 1 và ưu tiên 2
    firts_prioritize_element = [AB[i, j] for i in row for j in col
                                if np.abs(float(AB[i, j])) == float(int(1))]
    if firts_prioritize_element:
        pivot = min(firts_prioritize_element, key=abs)
    else:
        second_prioritize_element = [AB[i, j] for i in row for j in col]
        pivot = max(second_prioritize_element, key=abs)
    positions = [[i, j] for i in row for j in col if AB[i, j] == pivot]
    return pivot, positions
#Hàm khử gauss_jordan
def gauss_jordan(A,B):

    A, B = Matrix(N(A)), Matrix(N(B))
    rowused = []
    colused = []
    rowA, colA = A.shape
    AB = chop_small_values(A.row_join(B))
    display(AB)
    standardization = []

    #Tìm pivot và khử theo công thức (1)
    while ((len(rowused)!=rowA) and (len(colused)!=colA)):
        pivot = findpivot(AB,rowA,colA,rowused,colused)[0]
        if (abs(pivot) < 1e-10): break
        pivotrow, pivotcol = findpivot(AB,rowA,colA,rowused,colused)[1][0]
        for i in [x for x in range(rowA) if x!=pivotrow]:
            AB[i,:] = AB[i,:] - AB[pivotrow,:]*AB[i,pivotcol]/pivot
        rowused.append(pivotrow)
        colused.append(pivotcol)
        standardization.append([pivotrow,pivotcol])
        print("Khử")
        AB = chop_small_values(AB)
        display(AB)

    #Chuẩn hoá, chia các hàng cho pivot
    print("Chuẩn hoá")
    for i in range(len(standardization)):
        pivotrow, pivotcol = standardization[i]
        AB[pivotrow,:] = AB[pivotrow,:]/AB[pivotrow,pivotcol]
    return chop_small_values(AB), colused
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
#Giải AX=B
def gauss_jordan_solve(A, B):
    #Ma trận đã khử AB
    A, B = Matrix(N(A)), Matrix(N(B))
    AB,colused = gauss_jordan(A, B)
    display(AB)
    #Kiểm tra trường hợp vô nghiệm
    for i in indexpivot(AB):
        if i >= A.shape[1] and i < AB.shape[1]:
            print("Hệ vô nghiệm")
            return 
    
    print(f"Số phần tử tự do là {len(list(set(range(A.shape[1]))-set(colused)))}: X{sorted([x+1 for x in list(set(range(A.shape[1]))-set(colused))])}.")

    return
def print_detailed_solutions(AB_sorted, colused, n_vars):
    """
    Hàm trích xuất và in chi tiết các nghiệm x1, x2, ..., xn từ ma trận đã khử.
    
    Tham số:
    - AB_sorted: Ma trận mở rộng [A|B] đã được khử Gauss-Jordan và sắp xếp hàng.
    - colused: Danh sách các cột chứa phần tử cơ sở (pivot).
    - n_vars: Số lượng ẩn số của hệ phương trình (chính là số cột của ma trận A ban đầu).
    """
    # 1. Tạo tự động các biến x1, x2, ..., xn bằng Sympy
    X = symbols(f'x1:{n_vars + 1}') 
    
    # 2. Phân loại ẩn cơ sở và ẩn tự do
    free_cols = sorted(list(set(range(n_vars)) - set(colused)))
    
    if free_cols:
        print(f"=> Hệ có VÔ SỐ NGHIỆM. Số ẩn tự do là {len(free_cols)}: {[X[i] for i in free_cols]}")
    else:
        print("=> Hệ có NGHIỆM DUY NHẤT:")

    print("-" * 30)
    solutions = {}
    
    # 3. Gán các ẩn tự do bằng chính nó (ví dụ: x3 = x3)
    for c in free_cols:
        solutions[X[c]] = X[c]
        
    # 4. Trích xuất biểu thức cho các ẩn cơ sở
    # Duyệt qua từng hàng của ma trận
    for r in range(AB_sorted.shape[0]):
        # Tìm vị trí cột của phần tử khác 0 đầu tiên trong hàng r
        pivot_col = -1
        for c in range(n_vars):
            if abs(AB_sorted[r, c]) > 1e-10:
                pivot_col = c
                break
        
        # Nếu hàng r toàn số 0 (hoặc nằm ở cột B), bỏ qua
        if pivot_col == -1:
            continue
            
        # Biểu diễn nghiệm: x_pivot = Giá trị cột cuối (B) - Tổng(Hệ số * Ẩn tự do)
        expr = AB_sorted[r, -1] 
        for c in range(pivot_col + 1, n_vars):
            if abs(AB_sorted[r, c]) > 1e-10:
                expr = expr - AB_sorted[r, c] * X[c]
                
        solutions[X[pivot_col]] = expr

    # 5. In toàn bộ nghiệm ra màn hình theo thứ tự x1, x2, ...
    for i in range(n_vars):
        print(f"{X[i]} = {solutions[X[i]]}")
        
    return solutions
if __name__ == "__main__":
    output_path = str(__dir__ / "pp2_GaussJordan_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        A = Matrix(np.loadtxt(__dir__ / "GJ_input_A.txt"))
        print("A =")
        display(A)
        print("B =")
        B = Matrix(np.loadtxt(__dir__ / "GJ_input_B.txt"))
        display(B)
        print("A|B=")
        gauss_jordan_solve(A,B)
    print(f"Đã ghi kết quả vào {output_path}")





