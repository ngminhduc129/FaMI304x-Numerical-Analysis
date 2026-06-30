# Phân rã LU (Doolittle)

## Công thức toán học

### 1. Phân rã LU

Phân tích ma trận $A$ thành tích của ma trận tam giác dưới $L$ và ma trận tam giác trên $U$ (với đường chéo của $U$ bằng 1):

$$A = L \times U$$

Với $L$ là ma trận tam giác dưới và $U$ là ma trận tam giác trên có $U_{ii} = 1$.

### 2. Công thức tính

Với $i = 0, 1, \dots, n-1$:

**Phần tử đường chéo của $L$:**

$$L_{i,i} = A_{i,i} - \sum_{j=0}^{i-1} L_{i,j} \cdot U_{j,i}$$

**Phần tử tam giác dưới của $L$ ($k > i$):**

$$L_{k,i} = A_{k,i} - \sum_{j=0}^{i-1} L_{k,j} \cdot U_{j,i}$$

**Phần tử tam giác trên của $U$ ($k > i$):**

$$U_{i,k} = \frac{A_{i,k} - \sum_{j=0}^{i-1} L_{i,j} \cdot U_{j,k}}{L_{i,i}}$$

### 3. Kiểm tra

Sau khi phân rã, kiểm tra lại:

$$L \times U \approx A$$

## Thuật toán

**Đầu vào:** Ma trận vuông $A_{n \times n}$ đọc từ file `LU_input_A.txt`.  
**Đầu ra:** Ma trận $L$ (tam giác dưới) và $U$ (tam giác trên, $U_{ii}=1$) sao cho $A = L \times U$.

### Phần A: Khởi tạo

1. **Đọc dữ liệu đầu vào:**
   a. Mở file `LU_input_A.txt` và đọc ma trận $A$.
   b. Xác định kích thước $n$ của ma trận vuông $A$.
   c. Kiểm tra $n > 0$, nếu không hợp lệ thì thông báo lỗi và kết thúc.

2. **Khởi tạo ma trận $L$ và $U$:**
   a. Tạo ma trận $L$ kích thước $n \times n$, khởi tạo tất cả phần tử bằng $0$.
   b. Tạo ma trận $U$ kích thước $n \times n$, khởi tạo tất cả phần tử bằng $0$.
   c. Đặt đường chéo của $U$: với mọi $i$ từ $0$ đến $n-1$, gán $U[i][i] = 1$.

### Phần B: Phân rã LU

3. **Với mỗi $i$ từ $0$ đến $n-1$ (xử lý từng cột của $L$ và từng hàng của $U$):**

   a. **Tính phần tử đường chéo $L[i][i]$:**
      - Khởi tạo tổng $sum = 0$.
      - Với mỗi $j$ từ $0$ đến $i-1$:
        $$sum = sum + L[i][j] \times U[j][i]$$
      - Tính $L[i][i] = A[i][i] - sum$.

   b. **Kiểm tra phần tử đường chéo $L[i][i]$:**
      - Nếu $|L[i][i]| < 10^{-10}$ (xấp xỉ bằng 0):
        * Thông báo ma trận suy biến, không thể phân rã LU theo phương pháp Doolittle.
        * Kết thúc thuật toán.

   c. **Tính các phần tử tam giác dưới của $L$ tại cột $i$ (các hàng $k > i$):**
      - Với mỗi $k$ từ $i+1$ đến $n-1$:
        * Khởi tạo tổng $sum = 0$.
        * Với mỗi $j$ từ $0$ đến $i-1$:
          $$sum = sum + L[k][j] \times U[j][i]$$
        * Tính $L[k][i] = A[k][i] - sum$.

   d. **Tính các phần tử tam giác trên của $U$ tại hàng $i$ (các cột $k > i$):**
      - Với mỗi $k$ từ $i+1$ đến $n-1$:
        * Khởi tạo tổng $sum = 0$.
        * Với mỗi $j$ từ $0$ đến $i-1$:
          $$sum = sum + L[i][j] \times U[j][k]$$
        * Tính $U[i][k] = (A[i][k] - sum) / L[i][i]$.

   e. **Xuất trạng thái trung gian (tùy chọn):**
      - In ra ma trận $L$ và $U$ sau bước $i$ để theo dõi quá trình phân rã.

### Phần C: Kiểm tra kết quả

4. **Kiểm tra phân rã LU:**
   a. Tạo ma trận $A\_reconstructed$ kích thước $n \times n$, khởi tạo bằng $0$.
   b. Với mỗi $i$ từ $0$ đến $n-1$:
      - Với mỗi $j$ từ $0$ đến $n-1$:
        * Khởi tạo tổng $sum = 0$.
        * Với mỗi $k$ từ $0$ đến $n-1$:
          $$sum = sum + L[i][k] \times U[k][j]$$
        * Gán $A\_reconstructed[i][j] = sum$.
   c. So sánh $A\_reconstructed$ với $A$ gốc:
      - Với mỗi cặp $(i, j)$:
        * Nếu $|A\_reconstructed[i][j] - A[i][j]| > 10^{-10}$:
          - In cảnh báo sai số tại vị trí $(i, j)$.
   d. Nếu sai số tất cả đều nhỏ hơn $10^{-10}$, thông báo phân rã thành công.

### Phần D: Hoàn tất

5. **Trả về kết quả:**
   a. Trả về ma trận $L$ (tam giác dưới, đường chéo chính là các phần tử $L[i][i]$).
   b. Trả về ma trận $U$ (tam giác trên, đường chéo chính là $1$).
   c. Ghi kết quả ra file hoặc in ra màn hình theo yêu cầu.
