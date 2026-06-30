# Phương pháp Heun (Euler cải tiến / Improved Euler)

## Công thức toán học

Phương pháp Heun (còn gọi là Euler cải tiến) sử dụng công thức hình thang:

$$
Y_{n+1} = Y_n + \frac{h}{2}\left[F(x_n, Y_n) + F(x_{n+1}, Y^*_{n+1})\right]
$$

trong đó $Y^*_{n+1}$ là giá trị dự đoán từ Euler tiến:

$$
Y^*_{n+1} = Y_n + h \cdot F(x_n, Y_n)
$$

## Thuật toán

**Đầu vào:**
- $F(x, Y)$: hàm vế phải
- $x_0, Y_0$: điều kiện đầu
- $X_{end}$: giá trị kết thúc
- $N$: số bước

**Đầu ra:** Bảng giá trị xấp xỉ nghiệm tại các bước thời gian.

1. **Tính bước nhảy:**
   $$h = \frac{X_{end} - x_0}{N}$$

2. **Khởi tạo:**
   a. Đặt $x_{curr} = x_0$.
   b. Đặt $Y_{curr} = Y_0$.
   c. Lưu dòng đầu tiên: bước $0$, $x = x_0$, $Y = Y_0$.

3. **Vòng lặp chính** (với $j = 0, 1, \ldots, N-1$):
   a. **Bước dự đoán (Predictor):**
      - Tính độ dốc thứ nhất: $K_1 = F(x_{curr}, Y_{curr})$.
      - Dự đoán nghiệm tại $x_{next}$ bằng Euler tiến: $Y_{star} = Y_{curr} + h \times K_1$.
      - Xác định thời điểm tiếp theo: $x_{next} = x_{curr} + h$.
   b. **Bước hiệu chỉnh (Corrector):**
      - Tính độ dốc thứ hai tại điểm dự đoán: $K_2 = F(x_{next}, Y_{star})$.
      - Lấy trung bình hai độ dốc và cập nhật:
        $$Y_{next} = Y_{curr} + \frac{h}{2} \times (K_1 + K_2)$$
   c. **Cập nhật nghiệm:**
      - $Y_{curr} = Y_{next}$.
      - $x_{curr} = x_{next}$.
   d. **Lưu kết quả:** bước $j+1$, $x = x_{curr}$, $Y = Y_{curr}$.

4. **Trả về** bảng kết quả.
