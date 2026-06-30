# Đạo hàm bằng Lagrange (Lagrange Derivative)

## Công thức toán học

Cho các điểm nút $(x_0, y_0), (x_1, y_1), \ldots, (x_{n-1}, y_{n-1})$, đạo hàm tại điểm $c$ được xấp xỉ bằng:

$$f'(c) \approx \sum_{i=0}^{n-1} y_i L'_i(c)$$

với $L_i(x)$ là đa thức cơ sở Lagrange:

$$L_i(x) = \prod_{j \neq i} \frac{x - x_j}{x_i - x_j}$$

Đạo hàm của $L_i(x)$ tại $c$:

$$L'_i(c) = L_i(c) \sum_{j \neq i} \frac{1}{c - x_j}$$

## Trường hợp đặc biệt: $c = x_k$ (tại nút)

Khi $c$ trùng với nút $x_k$:

- $i = k$: $L'_k(x_k) = \sum_{j \neq k} \frac{1}{x_k - x_j}$
- $i \neq k$: $L'_i(x_k) = \frac{1}{x_i - x_k} \prod_{j \neq i,k} \frac{x_k - x_j}{x_i - x_j}$

## Thuật toán

**Đầu vào:** Bộ $(n+1)$ điểm dữ liệu $(x_0, y_0), (x_1, y_1), \ldots, (x_n, y_n)$ với $x_i$ phân biệt; điểm $c$ cần tính đạo hàm.

**Đầu ra:** Giá trị xấp xỉ $f'(c)$.

1. **Chọn tập con $m$ điểm gần nhất:**
   a. Tính khoảng cách $|c - x_i|$ cho mọi $i$.
   b. Chọn $m$ điểm có khoảng cách nhỏ nhất, sắp xếp chúng theo thứ tự $x$ tăng dần.
   c. Gọi các điểm đã chọn là $(\tilde{x}_0, \tilde{y}_0), (\tilde{x}_1, \tilde{y}_1), \ldots, (\tilde{x}_{m-1}, \tilde{y}_{m-1})$.

2. **Tính đa thức cơ sở $L_i(c)$ cho từng nút $i$:**
   a. Đặt $L_i(c) = 1$.
   b. Với mỗi $j \neq i$ trong khoảng $0$ đến $m-1$:
      - Nhân $L_i(c)$ với $\dfrac{c - \tilde{x}_j}{\tilde{x}_i - \tilde{x}_j}$.

3. **Tính đạo hàm $L'_i(c)$ cho từng nút $i$:**
   a. Nếu $c$ trùng với nút $\tilde{x}_k$ nào đó:
      - Với $i = k$: $L'_k(c) = \sum_{j \neq k} \dfrac{1}{\tilde{x}_k - \tilde{x}_j}$.
      - Với $i \neq k$: $L'_i(c) = \dfrac{1}{\tilde{x}_i - \tilde{x}_k} \prod_{j \neq i,k} \dfrac{\tilde{x}_k - \tilde{x}_j}{\tilde{x}_i - \tilde{x}_j}$.
   b. Nếu $c$ không trùng nút nào:
      - Tính tổng $S = \sum_{j \neq i} \dfrac{1}{c - \tilde{x}_j}$.
      - $L'_i(c) = L_i(c) \times S$.

4. **Tính đạo hàm xấp xỉ:**
   a. Khởi tạo $f'(c) = 0$.
   b. Với mỗi $i$ từ $0$ đến $m-1$: $f'(c) = f'(c) + \tilde{y}_i \times L'_i(c)$.

5. **Trả về** giá trị $f'(c)$.
