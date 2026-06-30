# Nội suy Lagrange (Lagrange Interpolation)

## Mô tả
Xây dựng đa thức nội suy Lagrange đi qua các điểm $(x_i, y_i)$ đã cho.

## Công thức toán học

Đa thức nội suy Lagrange bậc $n$:

$$P_n(x) = \sum_{i=0}^{n} y_i L_i(x)$$

với các đa thức cơ sở Lagrange:

$$L_i(x) = \frac{w_{n+1}(x)}{(x - x_i) w'_{n+1}(x_i)} = D_i \cdot C_i(x)$$

trong đó:

$$w_{n+1}(x) = \prod_{k=0}^{n} (x - x_k)$$

$$D_i = \frac{y_i}{w'_{n+1}(x_i)} = \frac{y_i}{\prod_{j \neq i} (x_i - x_j)}$$

$$C_i(x) = \frac{w_{n+1}(x)}{x - x_i} \quad \text{(chia tổng hợp)}$$

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm $(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)$.

**Đầu ra:** Đa thức nội suy Lagrange $P_n(x)$.

1. **Tính đa thức $w(x)$:**
   a. Khởi tạo $w(x) = 1$ (đa thức hằng).
   b. Với mỗi $k$ từ 0 đến $n$:
      - Nhân $w(x)$ hiện tại với $(x - x_k)$.
      - Kết quả là đa thức $w_{n+1}(x) = \prod_{k=0}^{n} (x - x_k)$ bậc $n+1$.
   c. Lưu các hệ số của $w(x)$ từ bậc $n+1$ đến bậc 0.

2. **Với mỗi $i$ từ 0 đến $n$:**
   a. **Tính mẫu số $w'(x_i)$:**
      - Với mỗi $j$ từ 0 đến $n$, nếu $j \neq i$:
        - Tính hiệu $(x_i - x_j)$.
        - Nhân dồn vào tích số $\text{denom} = \prod_{j \neq i} (x_i - x_j)$.
      - Tích số này chính là $w'_{n+1}(x_i)$.
   b. **Tính $D_i$:** Chia $y_i$ cho $\text{denom}$:
      $$D_i = \frac{y_i}{\prod_{j \neq i} (x_i - x_j)}$$
   c. **Tính $C_i(x)$ bằng chia tổng hợp (synthetic division):**
      - Lấy mảng hệ số của $w(x)$.
      - Chia $w(x)$ cho $(x - x_i)$:
        - Khởi tạo $c_0 = w_0$ (hệ số bậc cao nhất của $w$).
        - Với $j = 1$ đến $n-1$: $c_j = w_j + x_i \cdot c_{j-1}$.
      - Kết quả là đa thức $C_i(x)$ bậc $n-1$.
   d. **Tính $L_i(x)$:** Nhân $D_i$ với từng hệ số của $C_i(x)$:
      $$L_i(x) = D_i \cdot C_i(x)$$
   e. **Cộng dồn:** Cộng $L_i(x)$ vào đa thức tổng $P_n(x)$.

3. **Kết quả:** $P_n(x) = \sum_{i=0}^{n} L_i(x)$ là đa thức nội suy Lagrange đi qua tất cả điểm đã cho.
