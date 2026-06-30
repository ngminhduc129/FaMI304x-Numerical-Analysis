# Phân rã Cholesky

## Công thức toán học

### 1. Điều kiện áp dụng

Ma trận $A$ phải là ma trận **đối xứng** ($A = A^T$) và **xác định dương**.

### 2. Phân rã Cholesky

Phân tích ma trận $A$ thành tích của ma trận tam giác dưới $L$ và chuyển vị của nó:

$$A = L \times L^T$$

### 3. Công thức tính phần tử của $L$

Với $i = 0, 1, \dots, n-1$ và $j = 0, 1, \dots, i$:

**Phần tử đường chéo ($i = j$):**

$$L_{i,i} = \sqrt{A_{i,i} - \sum_{k=0}^{i-1} L_{i,k}^2}$$

**Phần tử ngoài đường chéo ($i > j$):**

$$L_{i,j} = \frac{A_{i,j} - \sum_{k=0}^{j-1} L_{i,k} \cdot L_{j,k}}{L_{j,j}}$$

### 4. Giải hệ phương trình

Sau khi có $L$, giải hệ $Ax = b$ qua hai bước:

**Thế tiến** (giải $Ly = b$):

$$y_i = \frac{b_i - \sum_{j=0}^{i-1} L_{i,j} \cdot y_j}{L_{i,i}}$$

**Thế lùi** (giải $L^T x = y$):

$$x_i = \frac{y_i - \sum_{j=i+1}^{n-1} L_{j,i} \cdot x_j}{L_{i,i}}$$

## Thuật toán

**Đầu vào:** Ma trận $A_{n \times n}$ đọc từ file `CLSK_input_A.txt`, vector $b_{n \times 1}$ đọc từ file đầu vào.  
**Đầu ra:** Ma trận $L$ (phân rã Cholesky) và vector nghiệm $x$ của hệ $Ax = b$.

### Phần A: Kiểm tra điều kiện đầu vào

1. **Đọc dữ liệu đầu vào:**
   a. Mở file `CLSK_input_A.txt`, đọc ma trận $A$.
   b. Xác định kích thước $n$ của $A$. Kiểm tra $n > 0$.
   c. Đọc vector $b$ có $n$ phần tử.

2. **Kiểm tra ma trận đối xứng:**
   a. Với mỗi $i$ từ $0$ đến $n-1$:
      - Với mỗi $j$ từ $0$ đến $n-1$:
        * Tính sai số $\delta = |A[i][j] - A[j][i]|$.
        * Nếu $\delta > 10^{-10}$:
          - Kết luận ma trận $A$ không đối xứng.
          - Thông báo lỗi và kết thúc thuật toán.
   b. Nếu tất cả sai số đều nhỏ hơn $10^{-10}$, kết luận $A$ là ma trận đối xứng.

### Phần B: Phân rã Cholesky

3. **Khởi tạo ma trận $L$:**
   a. Tạo ma trận $L$ kích thước $n \times n$, khởi tạo tất cả phần tử bằng $0$.

4. **Với mỗi $i$ từ $0$ đến $n-1$ (xử lý từng hàng của $L$):**
   a. **Với mỗi $j$ từ $0$ đến $i$ (xử lý từng cột từ trái sang phải trong hàng $i$):**

      b. **Trường hợp $i == j$ (phần tử đường chéo):**
         - Khởi tạo tổng $sum = 0$.
         - Với mỗi $k$ từ $0$ đến $i-1$:
           $$sum = sum + L[i][k] \times L[i][k]$$
         - Tính $L[i][i] = \sqrt{A[i][i] - sum}$.
         - **Kiểm tra:** nếu $A[i][i] - sum < 0$:
           * Kết luận ma trận $A$ không xác định dương.
           * Thông báo lỗi và kết thúc thuật toán.

      c. **Trường hợp $i > j$ (phần tử ngoài đường chéo):**
         - Khởi tạo tổng $sum = 0$.
         - Với mỗi $k$ từ $0$ đến $j-1$:
           $$sum = sum + L[i][k] \times L[j][k]$$
         - Kiểm tra $|L[j][j]| < 10^{-10}$:
           * Nếu đúng, thông báo lỗi (chia cho 0) và kết thúc.
         - Tính $L[i][j] = (A[i][j] - sum) / L[j][j]$.

### Phần C: Kiểm tra kết quả phân rã

5. **Kiểm tra phân rã Cholesky:**
   a. Tạo ma trận $A\_reconstructed$ kích thước $n \times n$, khởi tạo bằng $0$.
   b. Với mỗi $i$ từ $0$ đến $n-1$:
      - Với mỗi $j$ từ $0$ đến $n-1$:
        * Khởi tạo tổng $sum = 0$.
        * Với mỗi $k$ từ $0$ đến $n-1$:
          $$sum = sum + L[i][k] \times L[j][k]$$
        * Gán $A\_reconstructed[i][j] = sum$.
   c. So sánh $A\_reconstructed$ với $A$ gốc:
      - Nếu $|A\_reconstructed[i][j] - A[i][j]| > 10^{-10}$ với bất kỳ $(i, j)$ nào:
        * In cảnh báo sai số tại vị trí đó.
      - Nếu tất cả sai số đều nhỏ hơn $10^{-10}$, thông báo phân rã Cholesky thành công.

### Phần D: Giải hệ phương trình $Ax = b$

6. **Giải hệ $Ly = b$ bằng phương pháp thế tiến:**
   a. Tạo vector $y$ kích thước $n$, khởi tạo bằng $0$.
   b. Với mỗi $i$ từ $0$ đến $n-1$:
      - Khởi tạo tổng $sum = b[i]$.
      - Với mỗi $j$ từ $0$ đến $i-1$:
        $$sum = sum - L[i][j] \times y[j]$$
      - Kiểm tra $|L[i][i]| < 10^{-10}$:
        * Nếu đúng, thông báo lỗi (ma trận suy biến) và kết thúc.
      - Tính $y[i] = sum / L[i][i]$.

7. **Giải hệ $L^T x = y$ bằng phương pháp thế lùi:**
   a. Tạo vector $x$ kích thước $n$, khởi tạo bằng $0$.
   b. Với mỗi $i$ từ $n-1$ xuống $0$:
      - Khởi tạo tổng $sum = y[i]$.
      - Với mỗi $j$ từ $i+1$ đến $n-1$:
        $$sum = sum - L[j][i] \times x[j]$$
      - Kiểm tra $|L[i][i]| < 10^{-10}$:
        * Nếu đúng, thông báo lỗi và kết thúc.
      - Tính $x[i] = sum / L[i][i]$.

### Phần E: Kiểm tra nghiệm

8. **Kiểm tra nghiệm $x$:**
   a. Tính $Ax$:
      - Với mỗi $i$ từ $0$ đến $n-1$:
        * Khởi tạo $sum = 0$.
        * Với mỗi $j$ từ $0$ đến $n-1$:
          $$sum = sum + A[i][j] \times x[j]$$
        * So sánh $sum$ với $b[i]$:
          - Nếu $|sum - b[i]| > 10^{-10}$: in cảnh báo sai số.
   b. Nếu tất cả sai số đều nhỏ hơn $10^{-10}$, thông báo nghiệm chính xác.

### Phần F: Hoàn tất

9. **Trả về kết quả:**
   a. Trả về ma trận $L$ (tam giác dưới từ phân rã Cholesky).
   b. Trả về vector nghiệm $x$ của hệ $Ax = b$.
   c. Ghi kết quả ra file hoặc in ra màn hình theo yêu cầu.
