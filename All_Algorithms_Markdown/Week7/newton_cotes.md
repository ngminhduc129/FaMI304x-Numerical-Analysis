# Công thức Newton-Cotes tổng quát (Newton-Cotes)

## Công thức toán học

Chia đoạn $[a, b]$ thành $n$ đoạn con đều nhau, $h = \frac{b-a}{n}$. Các nút là $x_k = a + kh$ với $k = 0, 1, \ldots, n$.

Công thức Newton-Cotes bậc $n$:

$$\int_a^b f(x)\,dx \approx h \sum_{k=0}^{n} H_k \cdot f(x_k)$$

với $H_k$ là hệ số Newton-Cotes (còn gọi là hệ số Cotes), được tính bằng cách tích phân đa thức cơ sở Lagrange trên $[0, n]$:

$$H_k = \int_0^n L_k(t)\,dt = \int_0^n \prod_{j \neq k} \frac{t - j}{k - j}\,dt$$

### Các trường hợp đặc biệt

- $n = 1$: Công thức hình thang
- $n = 2$: Công thức Simpson 1/3
- $n = 4$: Công thức Boole

## Thuật toán

**Đầu vào:** Hàm $f(x)$; khoảng $[a, b]$; bậc $n$ của công thức Newton-Cotes.

**Đầu ra:** Giá trị xấp xỉ tích phân $I \approx \int_a^b f(x)\,dx$.

1. **Rời rạc hóa khoảng tích phân:**
   a. Tính $h = \dfrac{b-a}{n}$.
   b. Với $k = 0, 1, \ldots, n$:
      - $x_k = a + k \times h$.
      - Tính $y_k = f(x_k)$.

2. **Với mỗi $k$ từ $0$ đến $n$, xây dựng đa thức cơ sở Lagrange $L_k(t)$ trên $[0, n]$:**
   a. **Tử số:** $T_k(t) = \prod_{j \neq k} (t - j)$.
      - Khởi tạo $T_k(t) = 1$ (dạng đa thức).
      - Với mỗi $j \neq k$ từ $0$ đến $n$: nhân $T_k(t)$ với $(t - j)$.
   b. **Mẫu số:** $M_k = \prod_{j \neq k} (k - j)$.
      - Khởi tạo $M_k = 1$.
      - Với mỗi $j \neq k$ từ $0$ đến $n$: $M_k = M_k \times (k - j)$.

3. **Tính hệ số Newton-Cotes $H_k$:**
   a. Viết $T_k(t) = c_0 + c_1 t + c_2 t^2 + \cdots + c_n t^n$.
   b. Tính tích phân: $\int_0^n T_k(t)\,dt = \left[ \sum_{i=0}^{n} \dfrac{c_i}{i+1} t^{i+1} \right]_{t=0}^{t=n}$.
   c. $H_k = \dfrac{\int_0^n T_k(t)\,dt}{M_k}$.

4. **Tính tổng có trọng số:**
   a. Khởi tạo $I = 0$.
   b. Với mỗi $k$ từ $0$ đến $n$: $I = I + H_k \times y_k$.
   c. $I = h \times I$.

5. **Trả về** giá trị $I$.

### Nhận xét các trường hợp đặc biệt

- **$n = 1$ (hình thang):** $H_0 = H_1 = \dfrac{1}{2}$, $I \approx \dfrac{h}{2}(y_0 + y_1)$.
- **$n = 2$ (Simpson 1/3):** $H_0 = H_2 = \dfrac{1}{3}$, $H_1 = \dfrac{4}{3}$, $I \approx \dfrac{h}{3}(y_0 + 4y_1 + y_2)$.
