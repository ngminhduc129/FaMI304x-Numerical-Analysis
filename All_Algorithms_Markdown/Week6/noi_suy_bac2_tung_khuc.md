# Nội suy Spline bậc 2 (Quadratic Spline Interpolation)

## Mô tả
Nội suy bằng các đoạn parabol nối trơn giữa các điểm, đảm bảo liên tục đến đạo hàm cấp 1 ($C^1$).

## Công thức toán học

Mỗi đoạn $[x_i, x_{i+1}]$ có hàm spline bậc 2:

$$S_i(x) = a_2^{(i)} x^2 + a_1^{(i)} x + a_0^{(i)}$$

**Điều kiện**:
- Nội suy: $S_i(x_i) = y_i$, $S_i(x_{i+1}) = y_{i+1}$
- Liên tục $C^1$: $S'_i(x_{i+1}) = S'_{i+1}(x_{i+1})$
- Điều kiện biên: $S'_0(x_0) = z_0$ cho trước

**Giải các $z_i = S'_i(x_i)$**:

$$z_{i+1} = \frac{2(y_{i+1} - y_i)}{h_i} - z_i, \quad i = 0, 1, \dots, n-1$$

**Hệ số spline**:

$$\begin{aligned}
a_2^{(i)} &= \frac{z_{i+1} - z_i}{2h_i} \\
a_1^{(i)} &= z_i - 2a_2^{(i)} x_i \\
a_0^{(i)} &= y_i - a_2^{(i)} x_i^2 - a_1^{(i)} x_i
\end{aligned}$$

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm $(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)$ với $x_0 < x_1 < \dots < x_n$, đạo hàm biên $z_0 = S'_0(x_0)$, điểm cần nội suy $x$.

**Đầu ra:** Giá trị nội suy $S(x)$.

### Phần 1: Tính đạo hàm tại các nút

1. **Khởi tạo:** Gán $z_0$ là điều kiện biên đã cho (đạo hàm tại $x_0$).

2. **Tính $h_i$ cho mỗi đoạn:**
   $$h_i = x_{i+1} - x_i, \quad i = 0, 1, \dots, n-1$$

3. **Tính $z_i$ cho các nút còn lại:** Với $i$ từ 0 đến $n-1$:
   a. **Tính hiệu $y$:** $\Delta y_i = y_{i+1} - y_i$.
   b. **Tính tỷ số:** $\frac{2\Delta y_i}{h_i}$.
   c. **Tính $z_{i+1}$:**
      $$z_{i+1} = \frac{2(y_{i+1} - y_i)}{h_i} - z_i$$
      - Nhân đôi $\Delta y_i$, chia cho $h_i$, trừ $z_i$.
   d. Lưu $z_{i+1}$ cho bước tiếp theo.

4. **Kết thúc:** Ta có danh sách $z_0, z_1, \dots, z_n$ là đạo hàm bậc 1 tại các nút.

### Phần 2: Tính hệ số spline cho từng đoạn

1. **Với mỗi đoạn $i$ từ 0 đến $n-1$:**
   a. **Tính $a_2^{(i)}$ (hệ số bậc 2):**
      $$a_2^{(i)} = \frac{z_{i+1} - z_i}{2h_i}$$
      - Hiệu các đạo hàm: $z_{i+1} - z_i$.
      - Chia cho $2h_i$.
   b. **Tính $a_1^{(i)}$ (hệ số bậc 1):**
      $$a_1^{(i)} = z_i - 2a_2^{(i)} x_i$$
      - Nhân $2a_2^{(i)}$ với $x_i$.
      - Lấy $z_i$ trừ kết quả trên.
   c. **Tính $a_0^{(i)}$ (hệ số hằng):**
      $$a_0^{(i)} = y_i - a_2^{(i)} x_i^2 - a_1^{(i)} x_i$$
      - Tính $a_2^{(i)} \cdot x_i^2$.
      - Tính $a_1^{(i)} \cdot x_i$.
      - Lấy $y_i$ trừ cả hai kết quả trên.
   d. **Lưu bộ ba $(a_0^{(i)}, a_1^{(i)}, a_2^{(i)})$** cho đoạn $i$.

### Phần 3: Đánh giá tại điểm $x$

1. **Xác định đoạn chứa $x$:**
   a. Tìm $i$ sao cho $x_i \le x \le x_{i+1}$.
   b. Nếu $x < x_0$: lấy đoạn $i = 0$.
   c. Nếu $x > x_n$: lấy đoạn $i = n-1$.

2. **Tính giá trị nội suy:**
   $$S_i(x) = a_2^{(i)} \cdot x^2 + a_1^{(i)} \cdot x + a_0^{(i)}$$
   a. Tính $x^2$, nhân với $a_2^{(i)}$.
   b. Nhân $a_1^{(i)}$ với $x$.
   c. Cộng với $a_0^{(i)}$.

3. **Trả về $S_i(x)$** là giá trị nội suy spline bậc 2 tại $x$.
