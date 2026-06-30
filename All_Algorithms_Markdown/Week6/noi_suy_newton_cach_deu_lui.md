# Nội suy Newton bước đều lùi (Newton Backward Fixed Gap)

## Mô tả
Nội suy Newton cho các điểm cách đều (bước $h = \text{const}$) dùng sai phân hữu hạn lùi.

## Công thức toán học

**Sai phân hữu hạn lùi**:

$$\begin{aligned}
\nabla y_i &= y_i - y_{i-1} \\
\nabla^k y_i &= \nabla^{k-1} y_i - \nabla^{k-1} y_{i-1}
\end{aligned}$$

**Công thức Newton-Gregory lùi**:

Đặt $x = x_n + th$, $t = \frac{x - x_n}{h}$:

$$P_n(x) = y_n + \frac{\nabla y_n}{1!} t + \frac{\nabla^2 y_n}{2!} t(t+1) + \dots + \frac{\nabla^n y_n}{n!} t(t+1)\dots(t+n-1)$$

Hay viết gọn:

$$P_n(x) = \sum_{k=0}^{n} \frac{\nabla^k y_n}{k!} \prod_{j=0}^{k-1} (t + j)$$

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm $(x_i, y_i)$ với $x_i$ cách đều $h = x_{i+1} - x_i$, điểm cần nội suy $x$.

**Đầu ra:** Giá trị nội suy $P_n(x)$.

### Phần 1: Xây dựng bảng sai phân lùi

1. **Khởi tạo bảng:** Tạo bảng hai chiều $\nabla$ kích thước $(n+1) \times (n+1)$.
   a. Cột 0: $\nabla[i][0] = y_i$ với $i = 0, 1, \dots, n$.
   b. Các ô còn lại khởi tạo bằng 0.

2. **Tính sai phân lùi các cấp:** Với $j$ từ 1 đến $n$:
   a. Với $i$ từ $j$ đến $n$:
      - Lấy giá trị ô hiện tại cùng cột trước: $\nabla[i][j-1]$.
      - Lấy giá trị ô trên cùng cột trước: $\nabla[i-1][j-1]$.
      - Tính hiệu: $\nabla[i][j] = \nabla[i][j-1] - \nabla[i-1][j-1]$.
   b. Cột $j$ chứa sai phân lùi cấp $j$: $\nabla^j y_i$.
   c. Lưu ý: Có thể tính từ bảng sai phân tiến rồi lấy $\nabla^k y_n = \Delta^k y_{n-k}$.

3. **Lấy đường chéo lùi:** Trích xuất $C_i = \nabla[n][i]$ với $i = 0, 1, \dots, n$.
   a. $C_0 = y_n = \nabla^0 y_n$.
   b. $C_1 = \nabla y_n$.
   c. $C_2 = \nabla^2 y_n$.
   d. ...
   e. $C_n = \nabla^n y_n$.

### Phần 2: Tính giá trị nội suy

1. **Tính $t$:** $t = \frac{x - x_n}{h}$.
   a. Tính hiệu $x - x_n$.
   b. Chia cho bước $h$ để được $t$.

2. **Khởi tạo:** Đặt $P_n(x) = 0$, $B = 1$ (đa thức cơ sở bậc 0).

3. **Với mỗi $i$ từ 0 đến $n$:**
   a. **Tính hệ số $D_i$:** $D_i = \frac{C_i}{i!}$.
      - Với $i = 0$: $D_0 = C_0 = y_n$.
      - Với $i \ge 1$: chia $C_i$ cho $i!$.
   b. **Xây dựng đa thức cơ sở $B_i(t)$:** (dùng dấu + thay vì -)
      - Nếu $i = 0$: $B_0(t) = 1$.
      - Nếu $i \ge 1$: $B_i(t) = B_{i-1}(t) \cdot (t + (i-1))$.
        - Nhân đa thức $B_{i-1}(t)$ với $(t + i - 1)$.
   c. **Tính số hạng:** $N_i = D_i \times B_i(t)$.
   d. **Cộng dồn:** $P_n(x) = P_n(x) + N_i$.

4. **Kết quả:** $P_n(x)$ là giá trị nội suy tại $x$.
