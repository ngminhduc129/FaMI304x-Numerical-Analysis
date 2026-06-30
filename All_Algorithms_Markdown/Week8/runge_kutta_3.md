# Phương pháp Runge-Kutta bậc 3 (RK3)

## Công thức toán học

Bảng Butcher cho RK3 tổng quát:

$$
\begin{array}{c|ccc}
0 & 0 & 0 & 0 \\
\alpha_2 & \beta_{11} & 0 & 0 \\
\alpha_3 & \beta_{21} & \beta_{22} & 0 \\
\hline
& r_1 & r_2 & r_3
\end{array}
$$

Hệ phương trình điều kiện bậc 3:

$$
\begin{aligned}
r_1 + r_2 + r_3 &= 1 \\
r_2 \alpha_2 + r_3 \alpha_3 &= \frac{1}{2} \\
r_2 \alpha_2^2 + r_3 \alpha_3^2 &= \frac{1}{3} \\
r_3 \beta_{22} \alpha_2 &= \frac{1}{6}
\end{aligned}
$$

Công thức cập nhật:

$$
\begin{aligned}
k_1 &= h \cdot F(x_n, Y_n) \\
k_2 &= h \cdot F(x_n + \alpha_2 h, Y_n + \beta_{11} k_1) \\
k_3 &= h \cdot F(x_n + \alpha_3 h, Y_n + \beta_{21} k_1 + \beta_{22} k_2) \\
Y_{n+1} &= Y_n + r_1 k_1 + r_2 k_2 + r_3 k_3
\end{aligned}
$$

Với $\alpha_2 = 1/3$, $\alpha_3 = 2/3$ ta có phương pháp Heun RK3.

## Thuật toán

**Đầu vào:** Giá trị hiện tại $x_n$, $Y_n$; bước nhảy $h$; hàm vế phải $F(x, Y)$; bộ hệ số RK3 ($\alpha_2, \alpha_3, \beta_{11}, \beta_{21}, \beta_{22}, r_1, r_2, r_3$).

**Đầu ra:** Giá trị xấp xỉ $Y_{n+1}$ tại $x_{n+1} = x_n + h$.

1. **Tính độ dốc thứ nhất $k_1$:**
   a. $k_1 = h \times F(x_n, Y_n)$.

2. **Tính độ dốc thứ hai $k_2$:**
   a. $x^{(2)} = x_n + \alpha_2 \times h$.
   b. $Y^{(2)} = Y_n + \beta_{11} \times k_1$.
   c. $k_2 = h \times F(x^{(2)}, Y^{(2)})$.

3. **Tính độ dốc thứ ba $k_3$:**
   a. $x^{(3)} = x_n + \alpha_3 \times h$.
   b. $Y^{(3)} = Y_n + \beta_{21} \times k_1 + \beta_{22} \times k_2$.
   c. $k_3 = h \times F(x^{(3)}, Y^{(3)})$.

4. **Cập nhật nghiệm:**
   $$Y_{n+1} = Y_n + r_1 \times k_1 + r_2 \times k_2 + r_3 \times k_3$$

5. **Trả về** $Y_{n+1}$.

### Bộ tham số Heun RK3 ($\alpha_2 = 1/3$, $\alpha_3 = 2/3$)

- $\beta_{11} = \dfrac{1}{3}$, $\beta_{21} = 0$, $\beta_{22} = \dfrac{2}{3}$.
- $r_1 = \dfrac{1}{4}$, $r_2 = 0$, $r_3 = \dfrac{3}{4}$.
