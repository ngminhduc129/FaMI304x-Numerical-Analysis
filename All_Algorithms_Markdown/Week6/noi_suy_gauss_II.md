# Nội suy Gauss lùi (Gauss Backward Interpolation)

## Mô tả
Nội suy bằng công thức Gauss lùi dùng sai phân trung tâm cho các điểm cách đều. Phù hợp khi điểm nội suy gần tâm của tập dữ liệu.

## Công thức toán học

**Công thức Gauss II (Gauss lùi)**:

Đặt $t = \frac{x - x_0}{h}$:

$$P_n(x) = y_0 + t\Delta y_{-1} + \frac{t(t+1)}{2!}\Delta^2 y_{-1} + \frac{t(t+1)(t-1)}{3!}\Delta^3 y_{-2} + \frac{t(t+1)(t-1)(t+2)}{4!}\Delta^4 y_{-2} + \dots$$

Tổng quát, các hệ số được chọn theo quy tắc: $y_0, \Delta y_{-1}, \Delta^2 y_{-1}, \Delta^3 y_{-2}, \Delta^4 y_{-2}, \dots$

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm $(x_i, y_i)$ với $x_i$ cách đều $h$, điểm cần nội suy $x$ (gần tâm).

**Đầu ra:** Giá trị nội suy $P_n(x)$.

### Phần 1: Xây dựng bảng sai phân hữu hạn

1. **Khởi tạo bảng:** Tạo bảng sai phân hai chiều với $n+1$ hàng và $n+1$ cột.
   a. Cột 0: gán $y_i$ cho hàng $i$.
   b. $x_0$ là điểm trung tâm.

2. **Tính sai phân tiến các cấp:** Với $j$ từ 1 đến $n$:
   a. Với $i$ từ 0 đến $n - j$:
      $$\Delta^j y_i = \Delta^{j-1} y_{i+1} - \Delta^{j-1} y_i$$
   b. Lưu giá trị vào ô $(i, j)$ của bảng.

### Phần 2: Chọn hệ số Gauss II

1. **Tính $t$:** $t = \frac{x - x_0}{h}$.
   a. Hiệu $x - x_0$.
   b. Chia cho $h$.

2. **Duyệt $i$ từ 0 đến $n$:**
   a. Tính chỉ số hàng: $\text{hàng} = x_0 - \lfloor (i+1)/2 \rfloor$.
   b. Với $i$ chẵn ($i = 2k$):
      - Hàng: $x_0 - k$, tức $x_0 - \lfloor (i+1)/2 \rfloor = x_0 - k$.
      - Cột: $i$.
      - Hệ số: $\Delta^i y_{-k}$.
   c. Với $i$ lẻ ($i = 2k+1$):
      - Hàng: $x_0 - (k+1)$, tức $x_0 - \lfloor (i+1)/2 \rfloor = x_0 - (k+1)$.
      - Cột: $i$.
      - Hệ số: $\Delta^i y_{-(k+1)}$.
   d. Cụ thể:
      - $i = 0$: $y_0$ (hàng 0, cột 0).
      - $i = 1$: $\Delta y_{-1}$ (hàng -1, cột 1).
      - $i = 2$: $\Delta^2 y_{-1}$ (hàng -1, cột 2).
      - $i = 3$: $\Delta^3 y_{-2}$ (hàng -2, cột 3).
      - $i = 4$: $\Delta^4 y_{-2}$ (hàng -2, cột 4).
      - ...

### Phần 3: Xây dựng và tính đa thức

1. **Khởi tạo:** $P_n(x) = 0$, $B_0(t) = 1$.

2. **Với mỗi $i$ từ 0 đến $n$:**
   a. **Lấy hệ số $C_i$** từ bảng theo quy tắc Gauss II.
   b. **Xây dựng đa thức cơ sở $B_i(t)$:**
      - $B_0(t) = 1$.
      - $B_1(t) = t$.
      - Với $i \ge 2$:
        - Nếu $i$ chẵn ($i = 2k$): $B_{2k}(t) = B_{2k-1}(t) \cdot (t + k)$.
        - Nếu $i$ lẻ ($i = 2k+1$): $B_{2k+1}(t) = B_{2k}(t) \cdot (t - k)$.
   c. **Tính số hạng:** $N_i(t) = \frac{C_i}{i!} \cdot B_i(t)$.
      - Chia $C_i$ cho $i!$ (giai thừa của $i$).
      - Nhân kết quả với $B_i(t)$.
   d. **Cộng dồn:** $P_n(x) = P_n(x) + N_i$.

3. **Kết quả:** $P_n(x)$ là giá trị nội suy Gauss lùi tại $x$.
