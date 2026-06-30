# Phương pháp bình phương tối thiểu (Least Squares Method)

## Mô tả
Tìm hàm xấp xỉ tốt nhất cho tập dữ liệu theo nguyên lý bình phương tối thiểu. Hỗ trợ mô hình tuyến tính tổng quát với các hàm cơ sở tùy chỉnh, cùng các mô hình exponential (mũ) và power law (lũy thừa).

## Công thức toán học

**Bài toán**: Cho dữ liệu $(x_i, y_i)$, $i=1,\dots,n$, tìm $y \approx \sum_{j=1}^m a_j \phi_j(x)$

**Phương pháp Normal Equations**:

$$\mathbf{M} \mathbf{a} = \mathbf{b}'$$

trong đó:
- Ma trận Gram $\mathbf{M} = \mathbf{\Phi}^T \mathbf{\Phi}$, với $\Phi_{ij} = \phi_j(x_i)$ ($n \times m$)
- Vector $\mathbf{b}' = \mathbf{\Phi}^T \mathbf{y}$
- Hệ số $\mathbf{a} = [a_1, a_2, \dots, a_m]^T$

**Sai số**:
- $SSE = \mathbf{y}^T \mathbf{y} - \mathbf{a}^T \mathbf{b}'$
- $MSE = SSE / n$
- $RMSE = \sqrt{MSE}$

**Mô hình con**:
- Tuyến tính: $y = a_0 + a_1 x$ (hoặc tổng quát với $\sin$, $\cos$, ...)
- Exponential: $y = a e^{bx}$ $\rightarrow$ tuyến tính hóa: $\ln y = \ln a + b x$
- Power law: $y = a x^b$ $\rightarrow$ tuyến tính hóa: $\ln y = \ln a + b \ln x$

## Thuật toán

**Đầu vào:** Tập dữ liệu $(x_i, y_i)$ với $i = 1, \dots, n$, các hàm cơ sở $\phi_1(x), \phi_2(x), \dots, \phi_m(x)$, loại mô hình (tuyến tính / exponential / power law).

**Đầu ra:** Hệ số $a_1, a_2, \dots, a_m$ của mô hình xấp xỉ và các chỉ số sai số.

### Phần 1: Xử lý tiền dữ liệu (cho mô hình phi tuyến)

1. **Nếu mô hình là exponential ($y = a e^{bx}$):**
   a. Biến đổi $y$: $Y_i = \ln y_i$.
   b. Giữ nguyên $x$: $X_i = x_i$.
   c. Mô hình tuyến tính hóa: $Y = A + bx$ với $A = \ln a$.
   d. Đặt các hàm cơ sở: $\phi_1(x) = 1$, $\phi_2(x) = x$.

2. **Nếu mô hình là power law ($y = a x^b$):**
   a. Biến đổi $y$: $Y_i = \ln y_i$.
   b. Biến đổi $x$: $X_i = \ln x_i$.
   c. Mô hình tuyến tính hóa: $Y = A + bX$ với $A = \ln a$.
   d. Đặt các hàm cơ sở: $\phi_1(x) = 1$, $\phi_2(x) = x$.

3. **Nếu mô hình tuyến tính tổng quát:** Giữ nguyên $(x_i, y_i)$.

### Phần 2: Xây dựng ma trận và giải hệ

1. **Xây dựng ma trận $\Phi$ ($n$ hàng $\times$ $m$ cột):**
   a. Với mỗi hàng $i$ từ 1 đến $n$:
      - Với mỗi cột $j$ từ 1 đến $m$:
        $\Phi_{ij} = \phi_j(x_i)$
        - $\phi_1(x_i) = 1$ (hệ số hằng).
        - $\phi_2(x_i) = x_i$ (hệ số bậc 1).
        - $\phi_j(x_i)$ tùy theo hàm cơ sở thứ $j$.

2. **Tính ma trận Gram $\mathbf{M}$:**
   $$\mathbf{M} = \mathbf{\Phi}^T \mathbf{\Phi}$$
   a. Chuyển vị $\Phi$ được $\Phi^T$ kích thước $m \times n$.
   b. Nhân $\Phi^T$ với $\Phi$: kết quả là ma trận vuông $m \times m$.
   c. $\mathbf{M}_{kl} = \sum_{i=1}^{n} \phi_k(x_i) \cdot \phi_l(x_i)$.

3. **Tính vector $\mathbf{b}'$:**
   $$\mathbf{b}' = \mathbf{\Phi}^T \mathbf{y}$$
   a. Nhân $\Phi^T$ (kích thước $m \times n$) với vector $\mathbf{y}$ (kích thước $n \times 1$).
   b. Kết quả là vector cột $m \times 1$.
   c. $\mathbf{b}'_k = \sum_{i=1}^{n} \phi_k(x_i) \cdot y_i$.

4. **Giải hệ phương trình Normal:**
   $$\mathbf{M} \mathbf{a} = \mathbf{b}'$$
   a. Giải hệ $m$ phương trình $m$ ẩn số để tìm $\mathbf{a} = [a_1, a_2, \dots, a_m]^T$.

### Phần 3: Tính sai số

1. **Tính $SSE$ (Sum of Squared Errors):**
   $$SSE = \mathbf{y}^T \mathbf{y} - \mathbf{a}^T \mathbf{b}'$$
   a. Tính $\mathbf{y}^T \mathbf{y} = \sum_{i=1}^{n} y_i^2$.
   b. Tính $\mathbf{a}^T \mathbf{b}' = \sum_{j=1}^{m} a_j \cdot \mathbf{b}'_j$.
   c. Lấy hiệu hai giá trị trên.

2. **Tính $MSE$ (Mean Squared Error):**
   $$MSE = \frac{SSE}{n}$$
   a. Chia $SSE$ cho số điểm dữ liệu $n$.

3. **Tính $RMSE$ (Root Mean Squared Error):**
   $$RMSE = \sqrt{MSE}$$

### Phần 4: Chuyển ngược tham số (cho mô hình phi tuyến)

1. **Nếu mô hình exponential:**
   a. $b = a_2$ (giữ nguyên).
   b. $a = e^{a_1}$ (lấy $e$ mũ hệ số hằng).
   c. Mô hình: $y = a e^{bx}$.

2. **Nếu mô hình power law:**
   a. $b = a_2$ (giữ nguyên).
   b. $a = e^{a_1}$ (lấy $e$ mũ hệ số hằng).
   c. Mô hình: $y = a x^b$.

3. **Kết quả:** Trả về các hệ số $a_j$ và các chỉ số sai số $SSE$, $MSE$, $RMSE$.
