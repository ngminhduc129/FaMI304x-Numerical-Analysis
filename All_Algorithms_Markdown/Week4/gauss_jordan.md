# Phương pháp Gauss-Jordan

## Công thức toán học

### 1. Tìm phần tử trụ (Pivot Selection)

Hàm tìm pivot với hai mức ưu tiên:

- **Ưu tiên 1:** Chọn phần tử có giá trị bằng ±1 (đơn giản nhất).
- **Ưu tiên 2:** Nếu không có, chọn phần tử có trị tuyệt đối lớn nhất.

$$pivot = \begin{cases}
\min\{|a_{ij}| : a_{ij} \in \mathbb{Z}\} & \text{nếu tồn tại} \\
\max\{|a_{ij}|\} & \text{ngược lại}
\end{cases}$$

### 2. Công thức khử

Với pivot tại $(pivotrow, pivotcol)$ và hệ số $\lambda = A_{i, pivotcol} / pivot$, mỗi hàng $i \neq pivotrow$ được khử:

$$A_{i,:} = A_{i,:} - \frac{A_{i, pivotcol}}{pivot} \cdot A_{pivotrow,:}$$

### 3. Chuẩn hóa

Sau khi khử, các hàng chứa pivot được chia cho chính pivot:

$$A_{pivotrow,:} = \frac{A_{pivotrow,:}}{A_{pivotrow, pivotcol}}$$

---

## Thuật toán

**Mục tiêu:** Giải hệ phương trình tuyến tính $Ax = B$ bằng phương pháp Gauss-Jordan.
**Đầu vào:** Ma trận $A_{m \times n}$, vector $B_{m \times 1}$ đọc từ file đầu vào.
**Đầu ra:** Nghiệm của hệ phương trình (dạng tường minh nếu có nghiệm duy nhất, dạng tham số nếu vô số nghiệm, hoặc thông báo vô nghiệm).

### Phần A: Khởi tạo

**Bước 1:** Đọc dữ liệu đầu vào
   - **Bước 1.1:** Đọc ma trận $A$ từ file. Xác định số hàng $m$ và số cột $n$.
   - **Bước 1.2:** Đọc vector $B$ từ file. Kiểm tra $B$ có $m$ phần tử.
   - **Bước 1.3:** Nếu kích thước không khớp, thông báo lỗi và kết thúc.

**Bước 2:** Xây dựng ma trận mở rộng
   - **Bước 2.1:** Tạo ma trận $AB$ kích thước $m \times (n + 1)$.
   - **Bước 2.2:** Sao chép $A$ vào $AB[:, :n]$ và $B$ vào $AB[:, n]$.

**Bước 3:** Khởi tạo danh sách theo dõi
   - **Bước 3.1:** Tạo $row\_used$ rỗng để lưu các chỉ số hàng đã chọn làm pivot.
   - **Bước 3.2:** Tạo $col\_used$ rỗng để lưu các chỉ số cột đã chọn làm pivot.

### Phần B: Vòng lặp khử Gauss-Jordan

**Bước 4:** Lặp cho đến khi tất cả hàng hoặc cột đã được dùng

   **Bước 4.1:** Kiểm tra điều kiện dừng
   - Nếu $|row\_used| = m$: thoát vòng lặp.
   - Nếu $|col\_used| = n$: thoát vòng lặp.

   **Bước 4.2:** Tìm phần tử trụ (pivot)
   - Duyệt qua các hàng $i$ chưa trong $row\_used$ và cột $j$ chưa trong $col\_used$.
   - **Ưu tiên 1:** Tìm phần tử có $|AB[i][j]| = 1$. Nếu tìm thấy, chọn ngay.
   - **Ưu tiên 2:** Nếu không có phần tử $\pm 1$, chọn phần tử có trị tuyệt đối lớn nhất.
   - Nếu pivot tìm được có $|pivot| < 10^{-10}$: thoát vòng lặp.
   - Ghi nhận $pivotrow$ và $pivotcol$.

   **Bước 4.3:** Khử các hàng khác bằng pivot
   - Với mỗi hàng $i$ từ $0$ đến $m-1$ sao cho $i \neq pivotrow$:
     * Tính $factor = AB[i][pivotcol] / AB[pivotrow][pivotcol]$.
     * Nếu $|factor| < 10^{-10}$, bỏ qua.
     * Với mỗi cột $k$ từ $0$ đến $n$:
       $$AB[i][k] = AB[i][k] - factor \times AB[pivotrow][k]$$

   **Bước 4.4:** Cập nhật danh sách đã dùng
   - Thêm $pivotrow$ vào $row\_used$.
   - Thêm $pivotcol$ vào $col\_used$.

### Phần C: Chuẩn hóa

**Bước 5:** Chuẩn hóa các hàng chứa pivot
   - **Bước 5.1:** Với mỗi $pivotrow$ trong $row\_used$:
     * Tìm $pivotcol$ tương ứng trong $col\_used$.
     * Lấy $p = AB[pivotrow][pivotcol]$.
     * Nếu $|p| > 10^{-10}$: với mỗi cột $k$ từ $0$ đến $n$:
       $$AB[pivotrow][k] = AB[pivotrow][k] / p$$

### Phần D: Xử lý kết quả

**Bước 6:** Kiểm tra hệ vô nghiệm
   - **Bước 6.1:** Duyệt từng hàng $i$ từ $0$ đến $m-1$.
   - **Bước 6.2:** Nếu hàng $i$ có $|AB[i][j]| < 10^{-10}$ với mọi $j < n$ nhưng $|AB[i][n]| > 10^{-10}$:
     * Kết luận hệ vô nghiệm, kết thúc.

**Bước 7:** Xác định biến cơ sở và biến tự do
   - **Bước 7.1:** Biến cơ sở: các cột $j$ có trong $col\_used$.
   - **Bước 7.2:** Biến tự do: các cột $j$ từ $0$ đến $n-1$ không trong $col\_used$.
   - **Bước 7.3:** Ghi nhận $r$ (số biến cơ sở) và $f = n - r$ (số biến tự do).

**Bước 8:** Trích xuất nghiệm

   **Bước 8.1:** Nếu $f = 0$ (hệ có nghiệm duy nhất)
   - Với mỗi biến cơ sở $j$ tại hàng $i$ (hàng chứa pivot ở cột $j$):
     * $x_j = AB[i][n]$.
   - Trả về vector nghiệm $x = (x_0, x_1, \dots, x_{n-1})$.

   **Bước 8.2:** Nếu $f > 0$ (hệ vô số nghiệm)
   - Với mỗi biến tự do $j$: gán $x_j = t_j$ (tham số tự do).
   - Với mỗi biến cơ sở $j$ tại hàng $i$:
     * $x_j = AB[i][n] - \sum_{k \in free\_cols} AB[i][k] \times t_k$.
   - Trả về biểu diễn nghiệm dạng tham số.

**Bước 9:** Kết thúc thuật toán
   - In kết quả ra màn hình hoặc ghi file.
