# Phương pháp Adams-Bashforth bậc 2 (AB2)

## Công thức toán học

Phương pháp Adams-Bashforth bậc 2 là phương pháp đa bước tường minh:

$$
Y_{n+1} = Y_n + h \left[ \beta_0 F_n + \beta_1 F_{n-1} \right]
$$

Hệ số cho bậc 2:

$$
\beta_0 = \frac{3}{2}, \quad \beta_1 = -\frac{1}{2}
$$

Công thức tường minh:

$$
Y_{n+1} = Y_n + h \left( \frac{3}{2} F(x_n, Y_n) - \frac{1}{2} F(x_{n-1}, Y_{n-1}) \right)
$$

## Thuật toán

**Đầu vào:** Hàm $F(x, Y)$; điều kiện đầu $x_0$, $Y_0$; bước nhảy $h$; số bước $steps$.

**Đầu ra:** Bảng giá trị xấp xỉ nghiệm tại các bước thời gian.

1. **Khởi tạo hệ số:**
   a. Đặt $s = 2$ (bậc của phương pháp).
   b. Hằng số Adams-Bashforth bậc 2: $\beta_0 = \dfrac{3}{2}$, $\beta_1 = -\dfrac{1}{2}$.

2. **Pha A - Khởi tạo (Bootstrapping) bằng RK4:**
   a. **Bước 0:**
      - Lưu $Y_0$.
      - Tính $F_0 = F(x_0, Y_0)$.
      - Lưu $(x_0, Y_0, F_0)$ vào lịch sử.
   b. **Bước 1 (dùng RK4 để có điểm lịch sử thứ hai):**
      - Dùng phương pháp RK4 để tính $Y_1$ từ $Y_0$ với bước nhảy $h$.
      - $x_1 = x_0 + h$.
      - Tính $F_1 = F(x_1, Y_1)$.
      - Lưu $(x_1, Y_1, F_1)$ vào lịch sử.

3. **Pha B - Vòng lặp Adams-Bashforth** (với $n = 2, 3, \ldots, steps$):
   a. Lấy hai giá trị $F$ gần nhất từ lịch sử: $F_{n-1}$ và $F_{n-2}$.
   b. Tính số gia:
      $$\Delta Y = h \times \left( \beta_0 \times F_{n-1} + \beta_1 \times F_{n-2} \right)$$
   c. Cập nhật nghiệm:
      $$Y_n = Y_{n-1} + \Delta Y$$
   d. Cập nhật thời gian:
      $$x_n = x_{n-1} + h$$
   e. Tính $F_n = F(x_n, Y_n)$ và thêm $(x_n, Y_n, F_n)$ vào lịch sử.

4. **Trả về** bảng kết quả chứa các bước, $x$ và $Y$.
