# Nội suy Stirling (Stirling Interpolation)

## Mô tả
Nội suy bằng công thức Stirling dùng sai phân trung tâm, tối ưu cho các điểm đối xứng quanh điểm nội suy. Yêu cầu số điểm lẻ.

## Công thức toán học

Đặt $t = \frac{x - x_0}{h}$. Công thức Stirling:

$$\begin{aligned}
P_n(x) &= y_0 + \frac{t}{1!}\left(\frac{\Delta y_{-1} + \Delta y_0}{2}\right) + \frac{t^2}{2!}\Delta^2 y_{-1} \\
&\quad + \frac{t(t^2 - 1^2)}{3!}\left(\frac{\Delta^3 y_{-2} + \Delta^3 y_{-1}}{2}\right) + \frac{t^2(t^2 - 1^2)}{4!}\Delta^4 y_{-2} + \dots
\end{aligned}$$

**Quy tắc chọn hệ số**:
- Bậc chẵn $2k$: $\Delta^{2k} y_{-k}$
- Bậc lẻ $2k-1$: $\left(\Delta^{2k-1} y_{-k} + \Delta^{2k-1} y_{-(k-1)}\right) / 2$

**Đa thức cơ sở**:
- $B_0(t) = 1$, $B_1(t) = t$
- $B_{2k}(t) = t \cdot B_{2k-1}(t)$
- $B_{2k+1}(t) = (t^2 - k^2) \cdot B_{2k-1}(t)$

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm cách đều $h$, với $n$ chẵn (số điểm lẻ), điểm nội suy $x$ gần tâm $x_0$.

**Đầu ra:** Giá trị nội suy $P_n(x)$.

### Phần 1: Xây dựng bảng sai phân

1. **Khởi tạo bảng:** Tạo bảng sai phân hai chiều.
   a. Cột 0: gán $y_i$.
   b. Chọn $x_0$ là điểm trung tâm.

2. **Tính sai phân tiến các cấp:** Với $j$ từ 1 đến $n$:
   $$\Delta^j y_i = \Delta^{j-1} y_{i+1} - \Delta^{j-1} y_i$$

### Phần 2: Chọn hệ số Stirling

1. **Tính $t$:** $t = \frac{x - x_0}{h}$.

2. **Với mỗi bậc $i$ từ 0 đến $n$:**
   a. **Nếu $i = 0$:** $C_0 = y_0$.
   b. **Nếu $i$ chẵn ($i = 2k$, $k \ge 1$):**
      - Chọn $\Delta^{2k} y_{-k}$ từ bảng sai phân.
      - Hàng: $x_0 - k$, cột $2k$.
      - $C_{2k} = \Delta^{2k} y_{-k}$.
   c. **Nếu $i$ lẻ ($i = 2k-1$, $k \ge 1$):**
      - Lấy $\Delta^{2k-1} y_{-k}$ (hàng $x_0 - k$, cột $2k-1$).
      - Lấy $\Delta^{2k-1} y_{-(k-1)}$ (hàng $x_0 - (k-1)$, cột $2k-1$).
      - Tính trung bình cộng:
        $$C_{2k-1} = \frac{\Delta^{2k-1} y_{-k} + \Delta^{2k-1} y_{-(k-1)}}{2}$$

### Phần 3: Xây dựng và tính đa thức

1. **Khởi tạo:** $P_n(x) = 0$.

2. **Với mỗi $i$ từ 0 đến $n$:**
   a. **Xây dựng đa thức cơ sở $B_i(t)$:**
      - $B_0(t) = 1$.
      - $B_1(t) = t$.
      - Với $i \ge 2$:
        - Nếu $i$ chẵn ($i = 2k$): $B_{2k}(t) = t \cdot B_{2k-1}(t)$.
        - Nếu $i$ lẻ ($i = 2k+1$): $B_{2k+1}(t) = (t^2 - k^2) \cdot B_{2k-1}(t)$.
   b. **Tính số hạng:** $N_i = \frac{C_i}{i!} \cdot B_i(t)$.
      - Chia $C_i$ cho $i!$.
      - Nhân với $B_i(t)$.
   c. **Cộng dồn:** $P_n(x) = P_n(x) + N_i$.

3. **Kết quả:** $P_n(x)$ là giá trị nội suy Stirling tại $x$.
