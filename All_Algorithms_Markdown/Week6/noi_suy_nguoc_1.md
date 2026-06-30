# Nội suy hàm ngược (Inverse Interpolation — Function Approach)

## Mô tả
Tìm $x$ từ $y$ cho trước bằng cách nội suy Newton trên hàm ngược (đổi vai trò $x$ và $y$). Phương pháp này hoạt động trên các khoảng đơn điệu của hàm số.

## Công thức toán học

**Đổi vai trò**: Nếu hàm $y = f(x)$ đơn điệu trên $[a, b]$, tồn tại hàm ngược $x = f^{-1}(y)$. Nội suy Newton trên các điểm $(y_i, x_i)$ thay vì $(x_i, y_i)$:

$$x(y) \approx P_n(y)$$

**Công thức Newton chia sai phân** trên các điểm $(y_i, x_i)$:

$$P_n(y) = f[y_0] + f[y_0, y_1](y - y_0) + \dots + f[y_0, \dots, y_n](y - y_0)\dots(y - y_{n-1})$$

## Thuật toán

**Đầu vào:** Tập dữ liệu $(x_i, y_i)$ từ file CSV, giá trị $y_{\text{target}}$ cần tìm $x$.

**Đầu ra:** Giá trị $x \approx f^{-1}(y_{\text{target}})$.

### Phần 1: Đọc và phân tích dữ liệu

1. **Đọc dữ liệu:** Đọc file CSV chứa các cột $x$ và $y$.
   a. Lưu danh sách các điểm $(x_i, y_i)$.
   b. Kiểm tra dữ liệu không thiếu và sắp xếp theo $x$ tăng dần.

2. **Tìm khoảng đơn điệu:** Duyệt toàn bộ dữ liệu, xác định các đoạn liên tục mà $y$ tăng hoặc giảm đơn điệu.
   a. Với mỗi đoạn $[x_k, x_{k+m}]$, kiểm tra $y$ tăng nghiêm ngặt hoặc giảm nghiêm ngặt.
   b. Chỉ các khoảng đơn điệu mới đảm bảo tồn tại hàm ngược.

### Phần 2: Đổi vai trò và nội suy

1. **Chọn khoảng đơn điệu chứa $y_{\text{target}}$:**
   a. Với mỗi khoảng đơn điệu, kiểm tra nếu $y_{\text{target}}$ nằm giữa $y_{\min}$ và $y_{\max}$ của khoảng đó.
   b. Nếu tìm thấy, lấy các điểm trong khoảng đó.

2. **Đổi vai trò $x$ và $y$:** Tạo danh sách mới $(y_i, x_i)$ từ khoảng đơn điệu đã chọn.
   a. $y_i$ đóng vai trò là biến độc lập (thay cho $x$).
   b. $x_i$ đóng vai trò là biến phụ thuộc (thay cho $y$).

3. **Xây dựng bảng chia sai phân trên $(y_i, x_i)$:**
   a. Cột 0: $f[y_i] = x_i$.
   b. Cột $j$: với $i$ từ 0 đến $n-j$:
      $$f[y_i, \dots, y_{i+j}] = \frac{f[y_{i+1}, \dots, y_{i+j}] - f[y_i, \dots, y_{i+j-1}]}{y_{i+j} - y_i}$$

4. **Xây dựng đa thức Newton $P_n(y)$:**
   a. Lấy các hệ số $D_i$ từ đầu mỗi cột: $D_i = f[y_0, \dots, y_i]$.
   b. Khởi tạo $P_n(y) = 0$.
   c. Với mỗi $i$ từ 0 đến $n$:
      - Xây dựng đa thức cơ sở $B_{i-1}(y) = \prod_{k=0}^{i-1} (y - y_k)$.
      - Cộng dồn: $P_n(y) = P_n(y) + D_i \cdot B_{i-1}(y)$.

### Phần 3: Đánh giá

1. **Tính giá trị:** Thay $y = y_{\text{target}}$ vào $P_n(y)$.
   a. Sử dụng sơ đồ Horner để tính hiệu quả.
   b. Kết quả là $x \approx P_n(y_{\text{target}})$.

2. **Trả về kết quả:** Giá trị $x$ tìm được.
