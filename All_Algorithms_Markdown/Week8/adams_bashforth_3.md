# Phương pháp Adams-Bashforth bậc 3 (AB3)

## Công thức toán học

Phương pháp Adams-Bashforth bậc 3:

$$
Y_{n+1} = Y_n + h \left[ \beta_0 F_n + \beta_1 F_{n-1} + \beta_2 F_{n-2} \right]
$$

Hệ số cho bậc 3:

$$
\beta_0 = \frac{23}{12}, \quad \beta_1 = -\frac{16}{12}, \quad \beta_2 = \frac{5}{12}
$$

Công thức tường minh:

$$
Y_{n+1} = Y_n + h \left( \frac{23}{12} F_n - \frac{16}{12} F_{n-1} + \frac{5}{12} F_{n-2} \right)
$$

## Thuật toán

**Đầu vào:** Hàm $F(x, Y)$; điều kiện đầu $x_0$, $Y_0$; bước nhảy $h$; số bước $steps$.

**Đầu ra:** Bảng giá trị xấp xỉ nghiệm tại các bước thời gian.

1. **Khởi tạo hệ số:**
   a. Đặt $s = 3$ (bậc của phương pháp).
   b. Hằng số Adams-Bashforth bậc 3:
      $\beta_0 = \dfrac{23}{12}$, $\beta_1 = -\dfrac{16}{12}$, $\beta_2 = \dfrac{5}{12}$.

2. **Pha A - Khởi tạo (Bootstrapping) bằng RK4:**
   a. **Bước 0:**
      - Lưu $Y_0$.
      - Tính $F_0 = F(x_0, Y_0)$.
      - Lưu $(x_0, Y_0, F_0)$ vào lịch sử.
   b. **Bước 1 (dùng RK4):**
      - Tính $Y_1$ từ $Y_0$ bằng RK4 với bước nhảy $h$.
      - $x_1 = x_0 + h$.
      - Tính $F_1 = F(x_1, Y_1)$.
      - Lưu $(x_1, Y_1, F_1)$ vào lịch sử.
   c. **Bước 2 (dùng RK4):**
      - Tính $Y_2$ từ $Y_1$ bằng RK4 với bước nhảy $h$.
      - $x_2 = x_1 + h$.
      - Tính $F_2 = F(x_2, Y_2)$.
      - Lưu $(x_2, Y_2, F_2)$ vào lịch sử.

3. **Pha B - Vòng lặp Adams-Bashforth** (với $n = 3, 4, \ldots, steps$):
   a. Lấy ba giá trị $F$ gần nhất từ lịch sử: $F_{n-1}$, $F_{n-2}$, $F_{n-3}$.
   b. Tính số gia:
      $$\Delta Y = h \times \left( \beta_0 \times F_{n-1} + \beta_1 \times F_{n-2} + \beta_2 \times F_{n-3} \right)$$
   c. Cập nhật nghiệm:
      $$Y_n = Y_{n-1} + \Delta Y$$
   d. Cập nhật thời gian:
      $$x_n = x_{n-1} + h$$
   e. Tính $F_n = F(x_n, Y_n)$ và thêm $(x_n, Y_n, F_n)$ vào lịch sử.

4. **Trả về** bảng kết quả chứa các bước, $x$ và $Y$.
