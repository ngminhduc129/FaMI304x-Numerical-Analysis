# Phương pháp Adams-Moulton (Predictor-Corrector)

## Công thức toán học

Phương pháp Adams-Moulton kết hợp **Adams-Bashforth (dự đoán)** và **Adams-Moulton (hiệu chỉnh)**.

**Dự đoán (Predictor)** - Adams-Bashforth bậc $s$:

$$
Y_{n+1}^P = Y_n + h \sum_{j=0}^{s-1} \beta_j^{AB} F_{n-j}
$$

**Hiệu chỉnh (Corrector)** - Adams-Moulton bậc $s$ (ẩn):

$$
Y_{n+1} = Y_n + h \left[ \beta_0^{AM} F(x_{n+1}, Y_{n+1}) + \sum_{j=1}^{s-1} \beta_j^{AM} F_{n+1-j} \right]
$$

Hệ số Adams-Moulton bậc 4:

$$
\beta_0^{AM} = \frac{9}{24}, \quad \beta_1^{AM} = \frac{19}{24}, \quad \beta_2^{AM} = -\frac{5}{24}, \quad \beta_3^{AM} = \frac{1}{24}
$$

## Thuật toán

**Đầu vào:** Hàm $F(x, Y)$; $x_0$, $Y_0$; bước nhảy $h$; số bước $steps$; bậc $s$; ngưỡng hội tụ $tol$; số lần lặp tối đa $max\_iter$.

**Đầu ra:** Bảng giá trị xấp xỉ nghiệm tại các bước thời gian.

1. **Tạo hệ số:**
   a. **Hệ số Adams-Bashforth (dự đoán):** $ab\_betas = [\beta_0^{AB}, \beta_1^{AB}, \ldots, \beta_{s-1}^{AB}]$ dựa trên bậc $s$.
   b. **Hệ số Adams-Moulton (hiệu chỉnh):** $am\_betas = [\beta_0^{AM}, \beta_1^{AM}, \ldots, \beta_{s-1}^{AM}]$ dựa trên bậc $s$.

2. **Pha khởi tạo (Bootstrapping) bằng RK4:**
   a. $Y_0$ đã biết, tính $F_0 = F(x_0, Y_0)$.
   b. Với $k = 1, 2, \ldots, s-1$:
      - Dùng RK4 để tính $Y_k$ từ $Y_{k-1}$ với bước nhảy $h$.
      - $x_k = x_{k-1} + h$.
      - Tính $F_k = F(x_k, Y_k)$.
      - Lưu $(x_k, Y_k, F_k)$ vào lịch sử.
   c. Đặt $x_{curr} = x_{s-1}$, $Y_{curr} = Y_{s-1}$, $n = s-1$.

3. **Vòng lặp chính** (với mỗi bước từ $n = s-1$ đến $steps$):
   a. **Xác định thời điểm tiếp theo:**
      - $x_{next} = x_{curr} + h$.
   b. **Bước dự đoán (Predictor) - Adams-Bashforth:**
      - Lấy $s$ giá trị $F$ gần nhất từ lịch sử: $F_{n}, F_{n-1}, \ldots, F_{n-s+1}$.
      - Tính tổng dự đoán: $sum\_pred = \sum_{j=0}^{s-1} \beta_j^{AB} \times F_{n-j}$.
      - Dự đoán: $Y_{pred} = Y_{curr} + h \times sum\_pred$.
   c. **Bước hiệu chỉnh (Corrector) - Adams-Moulton với lặp điểm cố định:**
      - Tính phần tường minh (không phụ thuộc vào $Y_{n+1}$):
        $$explicit\_sum = \sum_{j=1}^{s-1} \beta_j^{AM} \times F_{n+1-j}$$
      - Đặt $Y_{guess} = Y_{pred}$.
      - Lặp (với $iter = 0, 1, \ldots, max\_iter$):
        - Tính $F_{guess} = F(x_{next}, Y_{guess})$.
        - Tính $Y_{new} = Y_{curr} + h \times \left( \beta_0^{AM} \times F_{guess} + explicit\_sum \right)$.
        - Tính sai số: $error = \max(|Y_{new} - Y_{guess}|)$.
        - Gán $Y_{guess} = Y_{new}$.
        - Nếu $error < tol$: thoát khỏi vòng lặp.
   d. **Cập nhật nghiệm:**
      - $Y_{curr} = Y_{guess}$.
      - $x_{curr} = x_{next}$.
      - Tính $F_{curr} = F(x_{curr}, Y_{curr})$.
      - Thêm $(x_{curr}, Y_{curr}, F_{curr})$ vào lịch sử.

4. **Trả về** bảng kết quả chứa các bước, $x$ và $Y$.
