# Nhân tử Lagrange tiến (Forward Lagrange Multipliers)

## Mô tả
Tính các nhân tử $w(x)$ và $C_i(x)$ trong nội suy Lagrange dạng tiến. Đây là các hàm phụ trợ để xây dựng đa thức nội suy Lagrange.

## Công thức toán học

**Hàm $w(x)$**: Đa thức cơ sở dạng tích:

$$w_{n+1}(x) = \prod_{i=0}^{n} (x - x_i) = (x - x_0)(x - x_1)\dots(x - x_n)$$

**Hàm $C_i(x)$**: Chia $w(x)$ cho $(x - x_i)$ bằng chia tổng hợp (synthetic division):

$$C_i(x) = \frac{w_{n+1}(x)}{x - x_i}$$

$C_i(x)$ là đa thức bậc $n-1$ với hệ số được tính bằng:

$$\begin{aligned}
c_0 &= w_0 \\
c_j &= w_j + x_i \cdot c_{j-1}, \quad j = 1, \dots, n-1
\end{aligned}$$

**Đa thức cơ sở Lagrange**:

$$L_i(x) = D_i \cdot C_i(x), \quad D_i = \frac{y_i}{\prod_{j \neq i} (x_i - x_j)}$$

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm $(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)$.

**Đầu ra:** Đa thức cơ sở $w(x)$, các nhân tử $C_i(x)$ và $L_i(x)$.

### Phần 1: Tính $w(x)$

1. **Khởi tạo:** Đặt $w(x) = [1]$ là đa thức hằng bậc 0.
2. **Vòng lặp nhân:** Với mỗi $i$ từ 0 đến $n$:
   a. Lấy đa thức $w(x)$ hiện tại (bậc $i$).
   b. Nhân $w(x)$ với $(x - x_i)$:
      - $w(x)$ mới có bậc $i+1$.
      - Hệ số bậc $i+1$ là hệ số bậc $i$ cũ.
      - Với $j$ từ $i$ xuống 0: $w_{j+1} = w_j$ (dịch chuyển), $w_j = w_j - x_i \cdot w_{j+1}$.
   c. Cập nhật $w(x)$ thành đa thức mới.
3. **Kết thúc:** $w(x)$ là đa thức bậc $n+1$ với hệ số $[w_0, w_1, \dots, w_{n+1}]$.

### Phần 2: Tính $C_i(x)$ và $L_i(x)$ cho từng $i$

1. **Với mỗi $i$ từ 0 đến $n$:**
   a. **Tính mẫu số $D_i$:**
      - Khởi tạo $\text{denom} = 1$.
      - Với mỗi $j$ từ 0 đến $n$, nếu $j \neq i$:
        $\text{denom} = \text{denom} \times (x_i - x_j)$.
      - $D_i = y_i / \text{denom}$.
   b. **Tính $C_i(x)$ bằng chia tổng hợp:**
      - Lấy mảng $w = [w_0, w_1, \dots, w_{n+1}]$ (bậc $n+1$).
      - Khởi tạo mảng $c = []$ độ dài $n$.
      - $c_0 = w_0$ (hệ số bậc $n$ của $C_i$).
      - Với $j = 1$ đến $n-1$:
        $c_j = w_j + x_i \cdot c_{j-1}$.
      - $C_i(x) = [c_0, c_1, \dots, c_{n-1}]$.
   c. **Tính $L_i(x)$:**
      - Nhân từng hệ số của $C_i(x)$ với $D_i$:
        $L_i(x) = D_i \cdot C_i(x)$.
      - $L_i(x)$ là đa thức bậc $n-1$.
2. **Kết quả:** Trả về $w(x)$, danh sách $C_i(x)$ và $L_i(x)$.
