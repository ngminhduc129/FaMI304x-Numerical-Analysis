# Nội suy Newton với chia sai phân (Newton Divided Differences)

## Mô tả
Xây dựng đa thức nội suy Newton dạng chia sai phân (divided differences) qua các điểm $(x_i, y_i)$. Hỗ trợ cả Newton tiến và Newton lùi.

## Công thức toán học

**Bảng chia sai phân**:

$$\begin{aligned}
f[x_i] &= y_i \\
f[x_i, x_{i+1}] &= \frac{f[x_{i+1}] - f[x_i]}{x_{i+1} - x_i} \\
f[x_i, x_{i+1}, \dots, x_{i+j}] &= \frac{f[x_{i+1}, \dots, x_{i+j}] - f[x_i, \dots, x_{i+j-1}]}{x_{i+j} - x_i}
\end{aligned}$$

**Đa thức Newton tiến** (lấy phần tử đầu mỗi cột):

$$P_n(x) = f[x_0] + f[x_0, x_1](x - x_0) + \dots + f[x_0, \dots, x_n](x - x_0)\dots(x - x_{n-1})$$

**Đa thức Newton lùi** (lấy phần tử cuối mỗi cột, đảo ngược $x$):

$$P_n(x) = f[x_n] + f[x_{n-1}, x_n](x - x_n) + \dots$$

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm $(x_0, y_0), (x_1, y_1), \dots, (x_n, y_n)$.

**Đầu ra:** Đa thức nội suy Newton $P_n(x)$ (dạng tiến hoặc lùi).

### Phần 1: Xây dựng bảng chia sai phân

1. **Khởi tạo bảng:** Tạo bảng hai chiều $F$ kích thước $(n+1) \times (n+1)$.
   a. Cột 0: $F[i][0] = y_i$ với $i = 0, 1, \dots, n$.
   b. Các ô còn lại khởi tạo bằng 0.

2. **Tính các cột sai phân:** Với $j$ từ 1 đến $n$:
   a. Với $i$ từ 0 đến $n - j$:
      - Lấy giá trị ô trên: $F[i+1][j-1]$.
      - Lấy giá trị ô dưới: $F[i][j-1]$.
      - Tính hiệu: $F[i+1][j-1] - F[i][j-1]$.
      - Tính mẫu số: $x_{i+j} - x_i$.
      - Chia hiệu cho mẫu số:
        $$F[i][j] = \frac{F[i+1][j-1] - F[i][j-1]}{x_{i+j} - x_i}$$

3. **Trích xuất hệ số:**
   a. **Newton tiến:** Lấy các phần tử đầu mỗi cột: $D_j = F[0][j]$ với $j = 0, 1, \dots, n$.
   b. **Newton lùi:** Lấy các phần tử cuối mỗi cột: $D_j = F[n-j][j]$ với $j = 0, 1, \dots, n$.

### Phần 2: Xây dựng đa thức Newton

1. **Khởi tạo:** Đặt $P_n(x) = 0$ (đa thức 0).
2. **Với mỗi $i$ từ 0 đến $n$:**
   a. **Lấy hệ số $D_i$** từ bảng chia sai phân (theo chế độ tiến hoặc lùi).
   b. **Xây dựng đa thức cơ sở $B_{i-1}(x)$:**
      - Nếu $i = 0$: $B_{-1}(x) = 1$.
      - Nếu $i \ge 1$: $B_{i-1}(x) = \prod_{k=0}^{i-1} (x - x_k)$.
        - Khởi tạo $B = [1, -x_0]$ cho $(x - x_0)$.
        - Với $k = 1$ đến $i-1$: nhân $B$ với $(x - x_k)$ bằng Horner ngược.
   c. **Tính số hạng:** $N_i(x) = D_i \cdot B_{i-1}(x)$ (nhân $D_i$ với từng hệ số của $B$).
   d. **Cộng dồn:** $P_n(x) = P_n(x) + N_i(x)$.

3. **Kết quả:** $P_n(x)$ là đa thức nội suy Newton.
