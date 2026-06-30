# Nội suy ngược bằng lặp (Inverse Interpolation — Fixed-Point Iteration)

## Mô tả
Tìm $x$ từ $y$ cho trước bằng kết hợp nội suy Newton bước đều và phương pháp lặp điểm cố định (fixed-point iteration).

## Công thức toán học

Cho $y_{\text{target}}$ cần tìm $x$. Đặt $t = \frac{x - x_0}{h}$, từ công thức Newton-Gregory tiến:

$$y_{\text{target}} = y_0 + \Delta y_0 \cdot t + \frac{\Delta^2 y_0}{2!} t(t-1) + \dots$$

Giải $t$ bằng lặp điểm cố định:

$$t_{k+1} = \phi(t_k) = \frac{1}{\Delta y_0}\left[y_{\text{target}} - y_0 - \sum_{i=2}^{n} \frac{\Delta^i y_0}{i!} \prod_{j=0}^{i-1} (t_k - j)\right]$$

Sau đó: $x = x_0 + t \cdot h$.

Tương tự cho Newton lùi:

$$t_{k+1} = \phi(t_k) = \frac{1}{\nabla y_n}\left[y_{\text{target}} - y_n - \sum_{i=2}^{n} \frac{\nabla^i y_n}{i!} \prod_{j=0}^{i-1} (t_k + j)\right]$$

Sau đó: $x = x_n + t \cdot h$.

## Thuật toán

**Đầu vào:** Tập dữ liệu $(x_i, y_i)$ với $x_i$ cách đều $h$, giá trị $y_{\text{target}}$, ngưỡng hội tụ $\varepsilon$, số lần lặp tối đa $K$.

**Đầu ra:** Giá trị $x$ thỏa $f(x) \approx y_{\text{target}}$.

### Phần 1: Chuẩn bị dữ liệu

1. **Đọc dữ liệu:** Đọc file CSV chứa $(x_i, y_i)$.
   a. Kiểm tra $x_i$ cách đều với bước $h = x_{i+1} - x_i$.

2. **Tìm khoảng đơn điệu:** Xác định đoạn $y$ đơn điệu chứa $y_{\text{target}}$.
   a. Nếu không tìm thấy, thông báo lỗi.

3. **Chọn điểm đầu:**
   a. **Newton tiến:** Chọn $x_0$ là điểm đầu của khoảng.
   b. **Newton lùi:** Chọn $x_n$ là điểm cuối của khoảng.
   c. Tính $h$ từ các điểm.

### Phần 2: Xây dựng bảng sai phân

1. **Tính bảng sai phân tiến hoặc lùi:**
   a. **Tiến:** Tính $\Delta^i y_0$ cho $i = 0, 1, \dots, n$.
      - $\Delta^0 y_0 = y_0$, $\Delta^1 y_0 = y_1 - y_0$, ...
   b. **Lùi:** Tính $\nabla^i y_n$ cho $i = 0, 1, \dots, n$.
      - $\nabla^0 y_n = y_n$, $\nabla^1 y_n = y_n - y_{n-1}$, ...

2. **Tính hệ số $C_i$:**
   a. Với Newton tiến: $C_i = \frac{\Delta^i y_0}{i!}$.
   b. Với Newton lùi: $C_i = \frac{\nabla^i y_n}{i!}$.

### Phần 3: Lặp điểm cố định

1. **Ước lượng ban đầu $t_0$:**
   a. Newton tiến: $t_0 = \frac{y_{\text{target}} - y_0}{\Delta y_0}$.
   b. Newton lùi: $t_0 = \frac{y_{\text{target}} - y_n}{\nabla y_n}$.

2. **Với mỗi bước lặp $k = 0, 1, 2, \dots$ đến $K$:**
   a. **Tính tổng các số hạng bậc cao** $S(t_k)$:
      - Khởi tạo $S = 0$.
      - Với $i = 2$ đến $n$:
        - Tính tích $P = \prod_{j=0}^{i-1} (t_k - j)$ (tiến) hoặc $\prod_{j=0}^{i-1} (t_k + j)$ (lùi).
        - $S = S + C_i \cdot P$.
   b. **Tính $t_{k+1}$:**
      - Newton tiến:
        $$t_{k+1} = \frac{y_{\text{target}} - y_0 - S}{\Delta y_0}$$
      - Newton lùi:
        $$t_{k+1} = \frac{y_{\text{target}} - y_n - S}{\nabla y_n}$$
   c. **Kiểm tra hội tụ:**
      - Nếu $|t_{k+1} - t_k| < \varepsilon$: dừng lặp.
      - Nếu $k = K$: thông báo chưa hội tụ.

3. **Cập nhật:** $t_k = t_{k+1}$ và tiếp tục vòng lặp.

### Phần 4: Kết quả

1. **Tính $x$:**
   a. Newton tiến: $x = x_0 + t \cdot h$.
   b. Newton lùi: $x = x_n + t \cdot h$.

2. **Trả về $x$** là giá trị cần tìm.
