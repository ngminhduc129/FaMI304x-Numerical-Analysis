# Nội suy Spline tuyến tính (Linear Spline Interpolation)

## Mô tả
Nội suy bằng các đoạn thẳng nối giữa các điểm dữ liệu (spline bậc 1), đảm bảo liên tục $C^0$.

## Công thức toán học

Mỗi đoạn $[x_i, x_{i+1}]$ có hàm spline tuyến tính:

$$S_i(x) = a_1^{(i)} x + a_0^{(i)}$$

với:

$$a_1^{(i)} = \frac{y_{i+1} - y_i}{x_{i+1} - x_i}$$

$$a_0^{(i)} = y_i - a_1^{(i)} x_i$$

Tính liên tục $C^0$: $S_i(x_{i+1}) = S_{i+1}(x_{i+1}) = y_{i+1}$.

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm $(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)$ với $x_0 < x_1 < \dots < x_n$, điểm cần nội suy $x$.

**Đầu ra:** Giá trị nội suy $S(x)$.

### Phần 1: Xây dựng các spline tuyến tính

1. **Với mỗi đoạn $i$ từ 0 đến $n-1$:**
   a. **Tính độ dài đoạn $h_i$:**
      $$h_i = x_{i+1} - x_i$$
   b. **Tính hệ số góc $a_1^{(i)}$:**
      $$a_1^{(i)} = \frac{y_{i+1} - y_i}{h_i}$$
      - Lấy hiệu các $y$: $y_{i+1} - y_i$.
      - Chia cho $h_i$.
   c. **Tính hệ số tự do $a_0^{(i)}$:**
      $$a_0^{(i)} = y_i - a_1^{(i)} \cdot x_i$$
      - Nhân $a_1^{(i)}$ với $x_i$.
      - Lấy $y_i$ trừ kết quả trên.
   d. **Lưu cặp hệ số $(a_0^{(i)}, a_1^{(i)})$** cho đoạn $i$.

### Phần 2: Đánh giá tại điểm $x$

1. **Xác định đoạn chứa $x$:**
   a. Tìm $i$ sao cho $x_i \le x \le x_{i+1}$.
   b. Nếu $x < x_0$: lấy đoạn đầu tiên ($i = 0$).
   c. Nếu $x > x_n$: lấy đoạn cuối cùng ($i = n-1$).

2. **Tính giá trị nội suy:**
   $$S_i(x) = a_1^{(i)} \cdot x + a_0^{(i)}$$
   a. Nhân $a_1^{(i)}$ với $x$.
   b. Cộng với $a_0^{(i)}$.

3. **Trả về $S_i(x)$** là giá trị nội suy tuyến tính tại $x$.
