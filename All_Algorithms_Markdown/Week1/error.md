# Sai số (Error)

## Các khái niệm cơ bản

- **Sai số tuyệt đối giới hạn**: $\Delta_X$ là số thỏa mãn $|X - x| \leq \Delta_X$, với $X$ là giá trị đúng, $x$ là giá trị gần đúng.
- **Sai số tương đối giới hạn**: $\delta_X = \frac{\Delta_X}{|X|}$ hoặc $\delta_x = \frac{\Delta_x}{|x|}$.

## Công thức lan truyền sai số

Cho hàm $y = f(x_1, x_2, \dots, x_n)$ với các đối số $x_i$ có sai số $\Delta_{x_i}$:

Sai số tuyệt đối giới hạn:

$$
\Delta_y = \sum_{i=1}^{n} \left| \frac{\partial f}{\partial x_i} \right| \Delta_{x_i}
$$

Sai số tương đối giới hạn:

$$
\delta_y = \frac{\Delta_y}{|y|}
$$

## Ví dụ 1: Tính sai số thể tích hình cầu

Cho $V = \frac{1}{6} \pi d^3$, $d = 3,7 \pm 0,05$ cm, $\pi = 3,14 \pm 0,0016$.

**Các bước tính toán:**

1. Xác định hàm số và sai số từng đối số:
   - Hàm số: $V(\pi, d) = \dfrac{1}{6} \pi d^3$.
   - Đối số $\pi$: giá trị $3,14$, sai số $\Delta_\pi = 0,0016$.
   - Đối số $d$: giá trị $3,7$, sai số $\Delta_d = 0,05$.

2. Tính đạo hàm riêng tại giá trị đo được:
   - $\dfrac{\partial V}{\partial \pi} = \dfrac{1}{6} d^3$.
   - Thay $d = 3,7$: $\dfrac{\partial V}{\partial \pi} = \dfrac{1}{6} \times (3,7)^3 = 8,44$.
   - $\dfrac{\partial V}{\partial d} = \dfrac{1}{2} \pi d^2$.
   - Thay $\pi = 3,14$, $d = 3,7$: $\dfrac{\partial V}{\partial d} = \dfrac{1}{2} \times 3,14 \times (3,7)^2 = 21,5$.

3. Tính sai số tuyệt đối giới hạn:
   - Công thức: $\Delta_V = \left|\dfrac{\partial V}{\partial \pi}\right| \Delta_\pi + \left|\dfrac{\partial V}{\partial d}\right| \Delta_d$.
   - Thay số: $\Delta_V = 8,44 \cdot 0,0016 + 21,5 \cdot 0,05 = 1,088 \approx 1,1$ cm$^3$.

4. Tính giá trị $V$ tại các giá trị đo:
   - $V = \dfrac{1}{6} \times 3,14 \times (3,7)^3 = 27,4$ cm$^3$.

5. Kết luận:
   - Kết quả: $V = 27,4 \pm 1,1$ cm$^3$.
   - Sai số tương đối: $\delta_V = \dfrac{\Delta_V}{|V|} = \dfrac{1,088}{27,4} = 0,04 = 4\%$.

## Ví dụ 2: Sai số ngược — chia đều ảnh hưởng

Cho hình trụ $V = \pi R^2 H$, $R = 2$m, $H = 3$m. Yêu cầu $\Delta_V = 0,1$m$^3$.

**Các bước tính toán:**

1. Xác định hàm số và yêu cầu:
   - Hàm số: $V(\pi, R, H) = \pi R^2 H$.
   - Sai số tổng cho phép: $\Delta_V = 0,1$ m$^3$.
   - Số đối số có sai số: $n = 3$ ($\pi$, $R$, $H$).

2. Tính đạo hàm riêng tại giá trị đo được:
   - $\dfrac{\partial V}{\partial \pi} = R^2 H = 2^2 \times 3 = 12$.
   - $\dfrac{\partial V}{\partial R} = 2\pi R H = 2 \times 3,14 \times 2 \times 3 = 37,7$.
   - $\dfrac{\partial V}{\partial H} = \pi R^2 = 3,14 \times 2^2 = 12,6$.

3. Áp dụng nguyên tắc ảnh hưởng ngang nhau:
   - Mỗi đối số đóng góp $\dfrac{\Delta_V}{n}$ vào sai số tổng.
   - Sai số cho phép của từng đối số:
     - $\Delta_\pi \leq \dfrac{0,1}{3 \cdot 12} = \dfrac{0,1}{36} = 0,00278 < 0,003$.
     - $\Delta_R \leq \dfrac{0,1}{3 \cdot 37,7} = \dfrac{0,1}{113,1} = 0,000884 < 0,001$.
     - $\Delta_H \leq \dfrac{0,1}{3 \cdot 12,6} = \dfrac{0,1}{37,8} = 0,00265 < 0,003$.

4. Kết luận:
   - Cần đo $\pi$ với sai số $< 0,003$.
   - Cần đo $R$ với sai số $< 0,001$ m.
   - Cần đo $H$ với sai số $< 0,003$ m.

## Ví dụ 3: Sai số ngược điều chỉnh

Diện tích hình tròn $S = \pi R^2$, yêu cầu $\delta_S \leq 1\%$.

**Các bước tính toán:**

1. Xác định hàm số và sai số tương đối yêu cầu:
   - Hàm số: $S(\pi, R) = \pi R^2$.
   - Yêu cầu: $\delta_S = \dfrac{\Delta_S}{S} \leq 0,01$.

2. Thiết lập phương trình sai số tương đối:
   - $\dfrac{\Delta_S}{S} = \dfrac{\Delta_\pi}{\pi} + 2\dfrac{\Delta_R}{R} = 0,01$.

3. Áp dụng nguyên tắc ảnh hưởng ngang nhau (hai đối số):
   - Chia đều: $\dfrac{\Delta_\pi}{\pi} = \dfrac{0,01}{2} = 0,005$.
   - $2\dfrac{\Delta_R}{R} = \dfrac{0,01}{2} = 0,005$.
   - Suy ra $\dfrac{\Delta_R}{R} = \dfrac{0,005}{2} = 0,0025$.

4. Tính sai số tuyệt đối cho phép:
   - $\Delta_\pi \leq 0,005 \times \pi = 0,005 \times 3,14 = 0,0157 < 0,016$.
   - $\Delta_R \leq 0,0025 \times R$ (phụ thuộc vào $R$ cụ thể).

5. Kết luận:
   - Cần $\Delta_\pi \leq 0,016$.
   - Cần $\Delta_R \leq 0,07$ cm (với $R$ cụ thể).
