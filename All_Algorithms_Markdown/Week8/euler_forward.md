# Phương pháp Euler tiến (Forward Euler)

## Công thức toán học

Cho bài toán Cauchy:

$$
\frac{dY}{dx} = F(x, Y), \quad Y(x_0) = Y_0
$$

Công thức Euler tiến:

$$
Y_{n+1} = Y_n + h \cdot F(x_n, Y_n)
$$

với $h = \frac{X_{end} - x_0}{N}$ là bước nhảy.

## Thuật toán

**Đầu vào:**
- Hàm $F(x, Y)$: hàm vế phải của hệ ODE
- $x_0$: giá trị đầu của $x$
- $Y_0$: vector điều kiện đầu $[y_{1,0}, y_{2,0}, \dots]$
- $X_{end}$: giá trị kết thúc của $x$
- $N$: số bước

**Đầu ra:** Bảng giá trị xấp xỉ nghiệm tại các bước thời gian.

1. **Tính bước nhảy:**
   $$h = \frac{X_{end} - x_0}{N}$$

2. **Khởi tạo:**
   a. Đặt $x_{curr} = x_0$.
   b. Đặt $Y_{curr} = Y_0$ (dạng vector).
   c. Lưu dòng đầu tiên: bước $0$, $x = x_0$, $Y = Y_0$.

3. **Vòng lặp chính** (với $j = 0, 1, 2, \ldots, N-1$):
   a. **Tính đạo hàm tại điểm hiện tại:**
      - $dY = F(x_{curr}, Y_{curr})$.
   b. **Cập nhật nghiệm theo công thức Euler tiến:**
      - $Y_{next} = Y_{curr} + h \times dY$.
   c. **Cập nhật thời gian:**
      - $x_{next} = x_{curr} + h$.
   d. **Gán giá trị mới thành giá trị hiện tại:**
      - $Y_{curr} = Y_{next}$.
      - $x_{curr} = x_{next}$.
   e. **Lưu kết quả:** bước $j+1$, $x = x_{curr}$, $Y = Y_{curr}$.

4. **Trả về** bảng kết quả chứa các bước, $x$ và các thành phần $y_k$.
