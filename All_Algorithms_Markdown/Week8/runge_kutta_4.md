# Phương pháp Runge-Kutta bậc 4 (RK4)

## Công thức toán học

Bảng Butcher cho RK4 cổ điển (Classic RK4):

$$
\begin{array}{c|cccc}
0 & 0 & 0 & 0 & 0 \\
\frac{1}{2} & \frac{1}{2} & 0 & 0 & 0 \\
\frac{1}{2} & 0 & \frac{1}{2} & 0 & 0 \\
1 & 0 & 0 & 1 & 0 \\
\hline
& \frac{1}{6} & \frac{1}{3} & \frac{1}{3} & \frac{1}{6}
\end{array}
$$

Công thức RK4 cổ điển:

$$
\begin{aligned}
k_1 &= h \cdot F(x_n, Y_n) \\
k_2 &= h \cdot F\left(x_n + \frac{h}{2}, Y_n + \frac{k_1}{2}\right) \\
k_3 &= h \cdot F\left(x_n + \frac{h}{2}, Y_n + \frac{k_2}{2}\right) \\
k_4 &= h \cdot F(x_n + h, Y_n + k_3) \\
Y_{n+1} &= Y_n + \frac{1}{6}(k_1 + 2k_2 + 2k_3 + k_4)
\end{aligned}
$$

## Thuật toán

**Đầu vào:** Giá trị hiện tại $x_n$, $Y_n$; bước nhảy $h$; hàm vế phải $F(x, Y)$.

**Đầu ra:** Giá trị xấp xỉ $Y_{n+1}$ tại $x_{n+1} = x_n + h$.

1. **Tính độ dốc thứ nhất $k_1$:**
   a. $k_1 = h \times F(x_n, Y_n)$.

2. **Tính độ dốc thứ hai $k_2$:**
   a. $x_{mid1} = x_n + \dfrac{h}{2}$.
   b. $Y_{mid1} = Y_n + \dfrac{k_1}{2}$.
   c. $k_2 = h \times F(x_{mid1}, Y_{mid1})$.

3. **Tính độ dốc thứ ba $k_3$:**
   a. $x_{mid2} = x_n + \dfrac{h}{2}$.
   b. $Y_{mid2} = Y_n + \dfrac{k_2}{2}$.
   c. $k_3 = h \times F(x_{mid2}, Y_{mid2})$.

4. **Tính độ dốc thứ tư $k_4$:**
   a. $x_{end} = x_n + h$.
   b. $Y_{end} = Y_n + k_3$.
   c. $k_4 = h \times F(x_{end}, Y_{end})$.

5. **Tổ hợp các độ dốc để cập nhật nghiệm:**
   $$Y_{n+1} = Y_n + \frac{1}{6} \times k_1 + \frac{1}{3} \times k_2 + \frac{1}{3} \times k_3 + \frac{1}{6} \times k_4$$

6. **Trả về** $Y_{n+1}$.
