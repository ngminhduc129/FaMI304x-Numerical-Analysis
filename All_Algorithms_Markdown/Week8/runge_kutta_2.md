# Phương pháp Runge-Kutta bậc 2 (RK2)

## Công thức toán học

Phương pháp RK2 tổng quát được biểu diễn bằng bảng Butcher:

$$
\begin{array}{c|cc}
0 & 0 & 0 \\
\alpha_2 & \beta_{11} & 0 \\
\hline
& r_1 & r_2
\end{array}
$$

Các hệ số thỏa mãn điều kiện bậc:

- $r_1 + r_2 = 1$
- $r_2 \cdot \alpha_2 = \frac{1}{2}$
- $\beta_{11} = \alpha_2$

Công thức cập nhật:

$$
\begin{aligned}
k_1 &= h \cdot F(x_n, Y_n) \\
k_2 &= h \cdot F(x_n + \alpha_2 h, Y_n + \beta_{11} k_1) \\
Y_{n+1} &= Y_n + r_1 k_1 + r_2 k_2
\end{aligned}
$$

Với $\alpha_2 = 1$ ta có phương pháp Heun; $\alpha_2 = 0.5$ ta có phương pháp Midpoint; $\alpha_2 = 2/3$ ta có phương pháp Ralston.

## Thuật toán

**Đầu vào:** Giá trị hiện tại $x_n$, $Y_n$; bước nhảy $h$; hàm vế phải $F(x, Y)$; tham số $\alpha_2$ (chọn $1$, $0.5$ hoặc $2/3$ tùy phương pháp).

**Đầu ra:** Giá trị xấp xỉ $Y_{n+1}$ tại $x_{n+1} = x_n + h$.

1. **Tính độ dốc thứ nhất $k_1$:**
   a. $k_1 = h \times F(x_n, Y_n)$.

2. **Tính độ dốc thứ hai $k_2$:**
   a. Xác định điểm trung gian:
      - $x_{mid} = x_n + \alpha_2 \times h$.
      - $Y_{mid} = Y_n + \alpha_2 \times k_1$.
   b. $k_2 = h \times F(x_{mid}, Y_{mid})$.

3. **Tính các hệ số $r_1$, $r_2$:**
   a. Từ điều kiện bậc:
      - $r_2 = \dfrac{1}{2 \alpha_2}$.
      - $r_1 = 1 - r_2$.

4. **Cập nhật nghiệm:**
   $$Y_{n+1} = Y_n + r_1 \times k_1 + r_2 \times k_2$$

5. **Trả về** $Y_{n+1}$.

### Các bộ tham số cụ thể

- **Heun ($\alpha_2 = 1$):** $r_1 = \dfrac{1}{2}$, $r_2 = \dfrac{1}{2}$.
- **Midpoint ($\alpha_2 = 0.5$):** $r_1 = 0$, $r_2 = 1$.
- **Ralston ($\alpha_2 = 2/3$):** $r_1 = \dfrac{1}{4}$, $r_2 = \dfrac{3}{4}$.
