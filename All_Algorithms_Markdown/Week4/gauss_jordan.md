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

## Thuật toán

**Đầu vào:** Ma trận $A_{m \times n}$, vector $B_{m \times 1}$ đọc từ file đầu vào.  
**Đầu ra:** Nghiệm của hệ phương trình (dạng tường minh nếu có nghiệm duy nhất, dạng tham số nếu vô số nghiệm, hoặc thông báo vô nghiệm).

### Phần A: Khởi tạo

1. **Đọc dữ liệu đầu vào:**
   a. Đọc ma trận $A$ từ file. Xác định số hàng $m$ và số cột $n$.
   b. Đọc vector $B$ từ file. Kiểm tra $B$ có $m$ phần tử.
   c. Nếu kích thước không khớp, thông báo lỗi và kết thúc.

2. **Xây dựng ma trận mở rộng:**
   a. Tạo ma trận $AB$ kích thước $m \times (n + 1)$.
   b. Sao chép $A$ vào $AB[:, :n]$ và $B$ vào $AB[:, n]$.

3. **Khởi tạo danh sách theo dõi:**
   a. Tạo danh sách $row\_used$ rỗng để lưu các chỉ số hàng đã được chọn làm pivot.
   b. Tạo danh sách $col\_used$ rỗng để lưu các chỉ số cột đã được chọn làm pivot.

### Phần B: Vòng lặp khử Gauss-Jordan

4. **Lặp cho đến khi tất cả các hàng hoặc tất cả các cột đã được dùng:**
   a. **Kiểm tra điều kiện dừng:**
      - Nếu số lượng hàng đã dùng $|row\_used| = m$: thoát vòng lặp (đã xử lý hết hàng).
      - Nếu số lượng cột đã dùng $|col\_used| = n$: thoát vòng lặp (đã xử lý hết cột).

   b. **Tìm phần tử trụ (pivot):**
      - Duyệt qua tất cả các hàng $i$ chưa có trong $row\_used$ và các cột $j$ chưa có trong $col\_used$.
      - Với mỗi cặp $(i, j)$, xét giá trị $|AB[i][j]|$.
      - **Ưu tiên 1:** Tìm phần tử có giá trị tuyệt đối bằng $1$ (tức $|AB[i][j]| = 1$). Nếu tìm thấy, chọn ngay pivot này.
      - **Ưu tiên 2:** Nếu không có phần tử nào bằng $\pm 1$, chọn phần tử có trị tuyệt đối lớn nhất trong số các phần tử còn lại.
      - Nếu pivot tìm được có trị tuyệt đối $< 10^{-10}$ (quá nhỏ), thoát vòng lặp vì không thể khử tiếp.
      - Ghi nhận $pivotrow$ và $pivotcol$ tìm được.

   c. **Khử các hàng khác bằng pivot:**
      - Với mỗi hàng $i$ từ $0$ đến $m-1$ sao cho $i \neq pivotrow$:
        * Tính hệ số khử: $factor = AB[i][pivotcol] / AB[pivotrow][pivotcol]$.
        * Nếu $|factor| < 10^{-10}$, bỏ qua hàng này (factor quá nhỏ, không cần khử).
        * Với mỗi cột $k$ từ $0$ đến $n$:
          $$AB[i][k] = AB[i][k] - factor \times AB[pivotrow][k]$$

   d. **Cập nhật danh sách đã dùng:**
      - Thêm $pivotrow$ vào $row\_used$.
      - Thêm $pivotcol$ vào $col\_used$.

### Phần C: Chuẩn hóa

5. **Chuẩn hóa các hàng chứa pivot:**
   a. Với mỗi $pivotrow$ trong $row\_used$ (theo thứ tự thời gian):
      - Tìm cột $pivotcol$ tương ứng trong $col\_used$.
      - Lấy giá trị pivot $p = AB[pivotrow][pivotcol]$.
      - Nếu $|p| > 10^{-10}$:
        * Với mỗi cột $k$ từ $0$ đến $n$:
          $$AB[pivotrow][k] = AB[pivotrow][k] / p$$
   b. Kết quả: mỗi hàng chứa pivot sẽ có pivot bằng $1$.

### Phần D: Xử lý kết quả

6. **Kiểm tra hệ vô nghiệm:**
   a. Duyệt từng hàng $i$ từ $0$ đến $m-1$.
   b. Nếu hàng $i$ có tất cả các phần tử từ cột $0$ đến $n-1$ đều bằng $0$ (hoặc rất nhỏ $< 10^{-10}$) nhưng $|AB[i][n]| > 10^{-10}$:
      - Kết luận hệ vô nghiệm.
      - Trả về thông báo vô nghiệm và kết thúc.

7. **Xác định biến cơ sở và biến tự do:**
   a. Biến cơ sở: các cột $j$ có trong $col\_used$.
   b. Biến tự do: các cột $j$ từ $0$ đến $n-1$ không có trong $col\_used$.
   c. Ghi nhận số biến cơ sở $r$ và số biến tự do $f = n - r$.

8. **Trích xuất nghiệm:**
   a. **Nếu $f = 0$ (hệ có nghiệm duy nhất):**
      - Với mỗi biến cơ sở $j$ tại hàng $i$ (hàng $i$ tương ứng chứa pivot ở cột $j$):
        * $x_j = AB[i][n]$
      - Trả về vector nghiệm $x = (x_0, x_1, \dots, x_{n-1})$.

   b. **Nếu $f > 0$ (hệ vô số nghiệm):**
      - Với mỗi biến tự do $j$, gán $x_j = t_j$ (tham số tự do).
      - Với mỗi biến cơ sở $j$ tại hàng $i$ (hàng chứa pivot ở cột $j$):
        * $x_j = AB[i][n] - \sum_{k \in free\_cols} AB[i][k] \times t_k$
      - Trả về biểu diễn nghiệm dạng tham số: mỗi biến cơ sở được biểu diễn qua các biến tự do.

9. **Kết thúc thuật toán:**
   a. In kết quả ra màn hình hoặc ghi ra file theo yêu cầu.
   b. Trả về ma trận nghiệm hoặc thông báo tương ứng.
