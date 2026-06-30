# Nội suy Bessel (Bessel Interpolation)

## Mô tả
Nội suy bằng công thức Bessel dùng sai phân trung tâm, phù hợp khi điểm nội suy ở giữa hai nút (giữa $x_0$ và $x_1$). Yêu cầu số điểm chẵn.

## Công thức toán học

Đặt $t = \frac{x - x_0}{h}$, với $x_0$ là điểm đầu của cặp trung tâm. Công thức Bessel:

$$\begin{aligned}
P_n(x) &= \frac{y_0 + y_1}{2} + \left(t - \frac{1}{2}\right)\Delta y_0 + \frac{t(t-1)}{2!}\left(\frac{\Delta^2 y_{-1} + \Delta^2 y_0}{2}\right) \\
&\quad + \frac{t(t-1)(t-1/2)}{3!}\Delta^3 y_{-1} + \dots
\end{aligned}$$

**Quy tắc chọn hệ số**:
- Bậc chẵn $2k$: $\left(\Delta^{2k} y_{-k} + \Delta^{2k} y_{-(k-1)}\right) / 2$
- Bậc lẻ $2k+1$: $\Delta^{2k+1} y_{-k}$

**Đa thức cơ sở**:
- $B_0(t) = 1$, $B_1(t) = t - 0.5$
- $B_{2k}(t) = B_{2k-2}(t) \cdot (t + k - 1)(t - k)$
- $B_{2k+1}(t) = B_{2k}(t) \cdot (t - 0.5)$

## Thuật toán

**Đầu vào:** Tọa độ $n+1$ điểm cách đều $h$, với $n$ lẻ (số điểm chẵn), điểm nội suy $x$ ở giữa $x_0$ và $x_1$.

**Đầu ra:** Giá trị nội suy $P_n(x)$.

### Phần 1: Xây dựng bảng sai phân

1. **Khởi tạo bảng:** Tạo bảng sai phân hai chiều.
   a. Cột 0: gán $y_i$.
   b. Chọn $x_0$ là điểm đầu của cặp trung tâm.

2. **Tính sai phân tiến các cấp:** Với $j$ từ 1 đến $n$:
   $$\Delta^j y_i = \Delta^{j-1} y_{i+1} - \Delta^{j-1} y_i$$

### Phần 2: Chọn hệ số Bessel

1. **Tính $t$:** $t = \frac{x - x_0}{h}$.

2. **Với mỗi bậc $i$ từ 0 đến $n$:**
   a. **Nếu $i = 0$:** $C_0 = \frac{y_0 + y_1}{2}$ (trung bình hai giá trị trung tâm).
   b. **Nếu $i$ lẻ ($i = 2k+1$, $k \ge 0$):**
      - Chọn $\Delta^{2k+1} y_{-k}$ từ bảng.
      - Hàng: $x_0 - k$, cột $2k+1$.
      - $C_{2k+1} = \Delta^{2k+1} y_{-k}$.
   c. **Nếu $i$ chẵn ($i = 2k$, $k \ge 1$):**
      - Lấy $\Delta^{2k} y_{-k}$ (hàng $x_0 - k$, cột $2k$).
      - Lấy $\Delta^{2k} y_{-(k-1)}$ (hàng $x_0 - (k-1)$, cột $2k$).
      - Tính trung bình cộng:
        $$C_{2k} = \frac{\Delta^{2k} y_{-k} + \Delta^{2k} y_{-(k-1)}}{2}$$

### Phần 3: Xây dựng và tính đa thức

1. **Khởi tạo:** $P_n(x) = 0$.

2. **Với mỗi $i$ từ 0 đến $n$:**
   a. **Xây dựng đa thức cơ sở $B_i(t)$:**
      - $B_0(t) = 1$.
      - $B_1(t) = t - 0.5$.
      - Với $i \ge 2$:
        - Nếu $i$ chẵn ($i = 2k$, $k \ge 1$):
          $B_{2k}(t) = B_{2k-2}(t) \cdot (t + k - 1)(t - k)$.
        - Nếu $i$ lẻ ($i = 2k+1$, $k \ge 1$):
          $B_{2k+1}(t) = B_{2k}(t) \cdot (t - 0.5)$.
   b. **Tính số hạng:** $N_i = \frac{C_i}{i!} \cdot B_i(t)$.
      - Chia $C_i$ cho $i!$.
      - Nhân với $B_i(t)$.
   c. **Cộng dồn:** $P_n(x) = P_n(x) + N_i$.

3. **Kết quả:** $P_n(x)$ là giá trị nội suy Bessel tại $x$.
