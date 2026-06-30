# Tích phân bằng Lagrange (Lagrange Integral)

## Công thức toán học

Cho các điểm nút $(x_0, y_0), (x_1, y_1), \ldots, (x_{n-1}, y_{n-1})$, tích phân xác định trên $[a, b]$ được xấp xỉ bằng:

$$I = \int_a^b f(x)\,dx \approx \sum_{i=0}^{n-1} y_i \cdot w_i$$

với trọng số $w_i$ được tính từ đa thức cơ sở Lagrange:

$$w_i = \int_a^b L_i(x)\,dx = \int_a^b \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}\,dx$$

## Thuật toán

**Đầu vào:** Bộ $n$ điểm dữ liệu $(x_0, y_0), (x_1, y_1), \ldots, (x_{n-1}, y_{n-1})$; khoảng tích phân $[a, b]$.

**Đầu ra:** Giá trị xấp xỉ tích phân $I \approx \int_a^b f(x)\,dx$.

1. **Với mỗi nút $i$ từ $0$ đến $n-1$, xây dựng đa thức tử số $P_i(x)$:**
   a. Khởi tạo $P_i(x) = 1$ (dưới dạng đa thức hệ số).
   b. Với mỗi $j \neq i$ từ $0$ đến $n-1$:
      - Nhân $P_i(x)$ với $(x - x_j)$ (khai triển và gộp hệ số).
   
2. **Tính mẫu số $D_i$ cho từng nút $i$:**
   a. Khởi tạo $D_i = 1$.
   b. Với mỗi $j \neq i$ từ $0$ đến $n-1$:
      - $D_i = D_i \times (x_i - x_j)$.

3. **Tính tích phân của đa thức tử số $P_i(x)$ trên $[a, b]$:**
   a. Viết $P_i(x) = c_0 + c_1 x + c_2 x^2 + \cdots + c_{n-1} x^{n-1}$.
   b. Tính nguyên hàm: $\int P_i(x)\,dx = c_0 x + \dfrac{c_1}{2} x^2 + \dfrac{c_2}{3} x^3 + \cdots + \dfrac{c_{n-1}}{n} x^{n}$.
   c. Tính $\int_a^b P_i(x)\,dx = \left[\sum_{k=0}^{n-1} \dfrac{c_k}{k+1} x^{k+1}\right]_{x=a}^{x=b}$.

4. **Tính trọng số $w_i$:**
   a. $w_i = \dfrac{\int_a^b P_i(x)\,dx}{D_i}$.

5. **Tính tích phân xấp xỉ:**
   a. Khởi tạo $I = 0$.
   b. Với mỗi $i$ từ $0$ đến $n-1$: $I = I + y_i \times w_i$.

6. **Trả về** giá trị $I$.
