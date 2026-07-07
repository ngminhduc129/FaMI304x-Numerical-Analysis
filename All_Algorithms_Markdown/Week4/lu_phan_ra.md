# Phân rã LU (Doolittle)

## Công thức toán học

### 1. Phân rã LU

Phân tích ma trận $A$ thành tích của ma trận tam giác dưới $L$ và ma trận tam giác trên $U$ (với đường chéo của $U$ bằng 1):

$$A = L \times U$$

### 2. Công thức tính

Với $i = 0, 1, \dots, n-1$:

**Phần tử đường chéo của $L$:**

$$L_{i,i} = A_{i,i} - \sum_{j=0}^{i-1} L_{i,j} \cdot U_{j,i}$$

**Phần tử tam giác dưới của $L$ ($k > i$):**

$$L_{k,i} = A_{k,i} - \sum_{j=0}^{i-1} L_{k,j} \cdot U_{j,i}$$

**Phần tử tam giác trên của $U$ ($k > i$):**

$$U_{i,k} = \frac{A_{i,k} - \sum_{j=0}^{i-1} L_{i,j} \cdot U_{j,k}}{L_{i,i}}$$

### 3. Kiểm tra

Sau khi phân rã, kiểm tra lại: $L \times U \approx A$.

---

## Thuật toán

**Mục tiêu:** Phân rã ma trận $A$ thành $L \times U$ theo phương pháp Doolittle.
**Đầu vào:** Ma trận vuông $A_{n \times n}$ đọc từ file `LU_input_A.txt`.
**Đầu ra:** Ma trận $L$ (tam giác dưới) và $U$ (tam giác trên, $U_{ii}=1$) sao cho $A = L \times U$.

### Phần A: Khởi tạo

**Bước 1:** Đọc dữ liệu đầu vào
   - **Bước 1.1:** Mở file `LU_input_A.txt` và đọc ma trận $A$.
   - **Bước 1.2:** Xác định kích thước $n$ của ma trận vuông $A$.
   - **Bước 1.3:** Kiểm tra $n > 0$, nếu không hợp lệ thì thông báo lỗi và kết thúc.

**Bước 2:** Khởi tạo ma trận $L$ và $U$
   - **Bước 2.1:** Tạo ma trận $L$ kích thước $n \times n$, khởi tạo tất cả bằng $0$.
   - **Bước 2.2:** Tạo ma trận $U$ kích thước $n \times n$, khởi tạo tất cả bằng $0$.
   - **Bước 2.3:** Đặt đường chéo của $U$: với mọi $i$, gán $U[i][i] = 1$.

### Phần B: Phân rã LU

**Bước 3:** Với mỗi $i$ từ $0$ đến $n-1$ (xử lý từng cột của $L$ và từng hàng của $U$)

   **Bước 3.1:** Tính phần tử đường chéo $L[i][i]$
   - Khởi tạo $sum = 0$.
   - Với mỗi $j$ từ $0$ đến $i-1$: $sum = sum + L[i][j] \times U[j][i]$.
   - Tính $L[i][i] = A[i][i] - sum$.

   **Bước 3.2:** Kiểm tra $L[i][i]$
   - Nếu $|L[i][i]| < 10^{-10}$: thông báo ma trận suy biến, kết thúc.

   **Bước 3.3:** Tính các phần tử tam giác dưới của $L$ tại cột $i$ (hàng $k > i$)
   - Với mỗi $k$ từ $i+1$ đến $n-1$:
     * $sum = 0$.
     * Với mỗi $j$ từ $0$ đến $i-1$: $sum = sum + L[k][j] \times U[j][i]$.
     * $L[k][i] = A[k][i] - sum$.

   **Bước 3.4:** Tính các phần tử tam giác trên của $U$ tại hàng $i$ (cột $k > i$)
   - Với mỗi $k$ từ $i+1$ đến $n-1$:
     * $sum = 0$.
     * Với mỗi $j$ từ $0$ đến $i-1$: $sum = sum + L[i][j] \times U[j][k]$.
     * $U[i][k] = (A[i][k] - sum) / L[i][i]$.

   **Bước 3.5:** (Tùy chọn) In trạng thái trung gian của $L$ và $U$.

### Phần C: Kiểm tra kết quả

**Bước 4:** Kiểm tra phân rã LU
   - **Bước 4.1:** Tạo $A\_reconstructed$ kích thước $n \times n$, khởi tạo bằng $0$.
   - **Bước 4.2:** Với mỗi $i$, $j$ từ $0$ đến $n-1$:
     * $sum = 0$.
     * Với mỗi $k$ từ $0$ đến $n-1$: $sum = sum + L[i][k] \times U[k][j]$.
     * Gán $A\_reconstructed[i][j] = sum$.
   - **Bước 4.3:** So sánh $A\_reconstructed$ với $A$ gốc:
     * Với mỗi $(i, j)$: nếu $|A\_reconstructed[i][j] - A[i][j]| > 10^{-10}$, in cảnh báo.
   - **Bước 4.4:** Nếu tất cả sai số đều nhỏ hơn $10^{-10}$, thông báo phân rã thành công.

### Phần D: Hoàn tất

**Bước 5:** Trả về kết quả
   - **Bước 5.1:** Trả về ma trận $L$ (tam giác dưới).
   - **Bước 5.2:** Trả về ma trận $U$ (tam giác trên, đường chéo là $1$).
   - **Bước 5.3:** Ghi kết quả ra file hoặc in màn hình.
