# Nội suy Newton bước đều tiến (Newton Forward Fixed Gap)

## Mô tả
Nội suy Newton cho các điểm cách đều (bước $h = \text{const}$) dùng sai phân hữu hạn tiến (finite differences).

## Công thức toán học

**Sai phân hữu hạn tiến**:

$$\begin{aligned}
\Delta y_i &= y_{i+1} - y_i \\
\Delta^2 y_i &= \Delta y_{i+1} - \Delta y_i \\
\Delta^k y_i &= \Delta^{k-1} y_{i+1} - \Delta^{k-1} y_i
\end{aligned}$$

**Công thức Newton-Gregory tiến**:

Đặt $x = x_0 + th$, $t = \frac{x - x_0}{h}$:

$$P_n(x) = y_0 + \frac{\Delta y_0}{1!} t + \frac{\Delta^2 y_0}{2!} t(t-1) + \dots + \frac{\Delta^n y_0}{n!} t(t-1)\dots(t-n+1)$$

Hay viết gọn:

$$P_n(x) = \sum_{k=0}^{n} \frac{\Delta^k y_0}{k!} \prod_{j=0}^{k-1} (t - j)$$

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm $(x_i, y_i)$ với $x_i$ cách đều $h = x_{i+1} - x_i$, điểm cần nội suy $x$.

**Đầu ra:** Giá trị nội suy $P_n(x)$.

### Phần 1: Xây dựng bảng sai phân tiến

1. **Khởi tạo bảng:** Tạo bảng hai chiều $\Delta$ kích thước $(n+1) \times (n+1)$.
   a. Cột 0: $\Delta[i][0] = y_i$ với $i = 0, 1, \dots, n$.
   b. Các ô còn lại khởi tạo bằng 0.

2. **Tính sai phân các cấp:** Với $j$ từ 1 đến $n$:
   a. Với $i$ từ 0 đến $n - j$:
      - Lấy giá trị ô dưới cùng cột trước: $\Delta[i+1][j-1]$.
      - Lấy giá trị ô hiện tại cùng cột trước: $\Delta[i][j-1]$.
      - Tính hiệu: $\Delta[i][j] = \Delta[i+1][j-1] - \Delta[i][j-1]$.
   b. Cột $j$ chứa sai phân cấp $j$: $\Delta^j y_i$.

3. **Lấy đường chéo tiến:** Trích xuất $C_i = \Delta[0][i]$ với $i = 0, 1, \dots, n$.
   a. $C_0 = y_0 = \Delta^0 y_0$.
   b. $C_1 = \Delta y_0$.
   c. $C_2 = \Delta^2 y_0$.
   d. ...
   e. $C_n = \Delta^n y_0$.

### Phần 2: Tính giá trị nội suy

1. **Tính $t$:** $t = \frac{x - x_0}{h}$.
   a. Tính hiệu $x - x_0$.
   b. Chia cho bước $h$ để được $t$.

2. **Khởi tạo:** Đặt $P_n(x) = 0$, $B = 1$ (đa thức cơ sở bậc 0).

3. **Với mỗi $i$ từ 0 đến $n$:**
   a. **Tính hệ số $D_i$:** $D_i = \frac{C_i}{i!}$.
      - Với $i = 0$: $D_0 = C_0 = y_0$.
      - Với $i \ge 1$: chia $C_i$ cho $i!$.
   b. **Xây dựng đa thức cơ sở $B_i(t)$:**
      - Nếu $i = 0$: $B_0(t) = 1$.
      - Nếu $i \ge 1$: $B_i(t) = B_{i-1}(t) \cdot (t - (i-1))$.
        - Nhân đa thức $B_{i-1}(t)$ với $(t - i + 1)$.
   c. **Tính số hạng:** $N_i = D_i \times B_i(t)$.
   d. **Cộng dồn:** $P_n(x) = P_n(x) + N_i$.

4. **Kết quả:** $P_n(x)$ là giá trị nội suy tại $x$.
