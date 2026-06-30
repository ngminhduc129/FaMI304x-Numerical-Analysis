# Sơ đồ Horner (Synthetic Division / Horner's Method)

## Mô tả
Tính giá trị đa thức $p(x)$ và các đạo hàm tại một điểm $c$ bằng sơ đồ Horner (chia đa thức tổng hợp).

## Công thức toán học

Cho đa thức $p(x) = a_n x^n + a_{n-1} x^{n-1} + \dots + a_1 x + a_0$, tính $p(c)$ bằng:

$$\begin{aligned}
b_n &= a_n \\
b_i &= a_i + c \cdot b_{i+1}, \quad i = n-1, \dots, 0 \\
p(c) &= b_0
\end{aligned}$$

Đa thức thương: $q(x) = b_n x^{n-1} + b_{n-1} x^{n-2} + \dots + b_1$

Đạo hàm cấp $i$ tại $c$: $p^{(i)}(c) = b_0^{(i)} \cdot i!$

**Horner ngược**: Cho $q(x)$ và $c$, khôi phục $p(x) = (x - c) q(x)$:

$$a_i = b_i - c \cdot b_{i+1}$$

**Hàm $w(x)$**: $w_{n+1}(x) = \prod_{k=0}^{n} (x - x_k)$ sử dụng Horner ngược lặp.

## Thuật toán

### Thuật toán 1: Synthetic Division (Tính $p(c)$)

**Đầu vào:** Mảng hệ số $a = [a_n, a_{n-1}, \dots, a_0]$ (bậc cao đến thấp), điểm $c$.

**Đầu ra:** Giá trị $p(c)$ và mảng hệ số $b$ của đa thức thương $q(x)$.

1. **Khởi tạo:** Gán $b_n = a_n$ (hệ số bậc cao nhất giữ nguyên).
   a. Tạo mảng $b$ có độ dài $n+1$.
   b. Đặt $b[n] = a[n]$.

2. **Vòng lặp xuôi:** Với $i$ chạy từ $n-1$ xuống 0:
   a. Lấy hệ số $a_i$ hiện tại.
   b. Lấy giá trị $b_{i+1}$ đã tính ở bước trước.
   c. Nhân $c$ với $b_{i+1}$.
   d. Cộng kết quả với $a_i$: $b_i = a_i + c \cdot b_{i+1}$.

3. **Kết quả:**
   a. $p(c) = b_0$ là giá trị đa thức tại $c$.
   b. $b_1, b_2, \dots, b_n$ là hệ số của đa thức thương $q(x)$ bậc $n-1$.

### Thuật toán 2: Tất cả đạo hàm tại $c$

**Đầu vào:** Mảng hệ số $a$, điểm $c$, số đạo hàm tối đa $k$.

**Đầu ra:** Danh sách $p(c), p'(c), p''(c), \dots, p^{(k)}(c)$.

1. **Khởi tạo:** Gán $b = a$ (bản sao mảng hệ số).
2. **Với mỗi bậc đạo hàm $i$ từ 0 đến $k$:**
   a. Áp dụng Synthetic Division lên mảng $b$ hiện tại để được $b_0$ mới.
   b. Lưu $b_0$ vào mảng kết quả tạm $r_i$.
   c. Nhân $r_i$ với $i!$ (giai thừa của $i$) để được $p^{(i)}(c)$.
   d. Cập nhật $b$ bằng mảng hệ số thương $[b_1, b_2, \dots]$ (bỏ $b_0$).
3. **Trả về** danh sách $[p(c), p'(c), p''(c), \dots]$.

### Thuật toán 3: Horner ngược (Khôi phục $p(x)$)

**Đầu vào:** Mảng hệ số $b$ của $q(x)$, điểm $c$, nhị thức $(x - c)$.

**Đầu ra:** Mảng hệ số $a$ của $p(x) = (x - c)q(x)$.

1. **Xác định bậc:** $m = \text{len}(b) - 1$ là bậc của $q(x)$.
2. **Tạo mảng $a$** có độ dài $m+2$ (bậc $m+1$).
3. **Gán:** $a_{m+1} = b_m$ (hệ số bậc cao nhất).
4. **Vòng lặp:** Với $i$ từ $m$ xuống 0:
   $$a_i = b_i - c \cdot a_{i+1}$$
   a. Lấy $b_i$ và $a_{i+1}$ đã biết.
   b. Nhân $c$ với $a_{i+1}$.
   c. Lấy $b_i$ trừ kết quả trên.
5. **Trả về** mảng $a$ là hệ số của $p(x)$.

### Thuật toán 4: Xây dựng hàm $w(x)$

**Đầu vào:** Danh sách các nút $x_0, x_1, \dots, x_n$.

**Đầu ra:** Mảng hệ số của $w_{n+1}(x) = \prod_{k=0}^{n} (x - x_k)$.

1. **Khởi tạo:** Đặt $w(x) = [1, -x_0]$ là đa thức $(x - x_0)$.
   a. Bậc 1: $w_1 x + w_0$, với $w_1 = 1$, $w_0 = -x_0$.
2. **Vòng lặp:** Với mỗi $k$ từ 1 đến $n$:
   a. Nhân $w(x)$ hiện tại với $(x - x_k)$.
   b. Sử dụng Horner ngược: $b = w$, $c = x_k$.
   c. Tính $a_i = b_i - x_k \cdot a_{i+1}$.
   d. Cập nhật $w$ thành mảng $a$ mới.
3. **Kết thúc:** $w(x)$ là đa thức bậc $n+1$ với hệ số từ bậc cao đến thấp.
