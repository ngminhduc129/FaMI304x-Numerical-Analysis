# Phương pháp Euler lùi (Backward Euler)

## Công thức toán học

Cho bài toán Cauchy:

$$
\frac{dY}{dx} = F(x, Y), \quad Y(x_0) = Y_0
$$

Công thức Euler lùi (ẩn):

$$
Y_{n+1} = Y_n + h \cdot F(x_{n+1}, Y_{n+1})
$$

Vì $Y_{n+1}$ xuất hiện ở cả hai vế, ta sử dụng **lặp điểm cố định** để giải:

1. **Dự đoán (Predictor):** Dùng Euler tiến để có $Y_{n+1}^{(0)}$
2. **Hiệu chỉnh (Corrector):** Lặp cho đến khi hội tụ

$$
Y_{next} = Y_{curr} + h \cdot F(x_{next}, Y_{guess})
$$

## Thuật toán

**Đầu vào:**
- $F(x, Y)$: hàm vế phải
- $x_0$: giá trị đầu
- $Y_0$: vector điều kiện đầu
- $X_{end}$: giá trị kết thúc
- $N$: số bước
- $tol$: ngưỡng hội tụ (mặc định $10^{-6}$)
- $max\_iter$: số lần lặp tối đa (mặc định 50)

**Đầu ra:** Bảng giá trị xấp xỉ nghiệm tại các bước thời gian.

1. **Tính bước nhảy:**
   $$h = \frac{X_{end} - x_0}{N}$$

2. **Khởi tạo:**
   a. Đặt $x_{curr} = x_0$.
   b. Đặt $Y_{curr} = Y_0$.
   c. Lưu dòng đầu tiên: bước $0$, $x = x_0$, $Y = Y_0$.

3. **Vòng lặp chính** (với $j = 0, 1, \ldots, N-1$):
   a. **Xác định thời điểm tiếp theo:**
      - $x_{next} = x_{curr} + h$.
   b. **Bước dự đoán (Predictor) - dùng Euler tiến:**
      - Tính $K_{pred} = F(x_{curr}, Y_{curr})$.
      - $Y_{guess} = Y_{curr} + h \times K_{pred}$.
   c. **Vòng lặp điểm cố định (Corrector)** (với $k = 0, 1, \ldots, max\_iter$):
      - Tính $K_{implicit} = F(x_{next}, Y_{guess})$.
      - Tính $Y_{new} = Y_{curr} + h \times K_{implicit}$.
      - Tính sai số: $error = \max(|Y_{new} - Y_{guess}|)$ (lấy max trên các thành phần).
      - Gán $Y_{guess} = Y_{new}$.
      - Nếu $error < tol$: thoát khỏi vòng lặp.
   d. **Cập nhật nghiệm:**
      - $Y_{curr} = Y_{guess}$.
      - $x_{curr} = x_{next}$.
   e. **Lưu kết quả:** bước $j+1$, $x = x_{curr}$, $Y = Y_{curr}$.

4. **Trả về** bảng kết quả.
