# Công thức Simpson 1/3 (Simpson's Rule)

## Công thức toán học

Chia đoạn $[a, b]$ thành $n$ đoạn con đều nhau ($n$ chẵn), $h = \frac{b-a}{n}$.

Công thức Simpson 1/3 ghép (Composite Simpson's Rule):

$$I = \int_a^b f(x)\,dx \approx \frac{h}{3} \left[ f(x_0) + f(x_n) + 4\sum_{k=1, \text{lẻ}}^{n-1} f(x_k) + 2\sum_{k=2, \text{chẵn}}^{n-2} f(x_k) \right]$$

Hay viết gọn:

$$I \approx \frac{h}{3} \left[ y_0 + y_n + 4(y_1 + y_3 + \ldots + y_{n-1}) + 2(y_2 + y_4 + \ldots + y_{n-2}) \right]$$

## Sai số

Sai số của công thức Simpson ghép:

$$|E| \leq \frac{M_4 (b-a)^5}{180 n^4}$$

trong đó $M_4 = \max_{x \in [a,b]} |f^{(4)}(x)|$.

Để đạt sai số $\varepsilon$, cần:

$$n > \left( \frac{M_4 (b-a)^5}{180 \varepsilon} \right)^{1/4}$$

và $n$ phải là số chẵn.

## Thuật toán

**Đầu vào:** Hàm $f(x)$; khoảng $[a, b]$; sai số cho phép $\varepsilon$ (hoặc số đoạn $n$ chẵn).

**Đầu ra:** Giá trị xấp xỉ tích phân $I \approx \int_a^b f(x)\,dx$.

1. **Xác định $M_4$ (nếu dùng sai số):**
   a. Tìm $f^{(4)}(x)$.
   b. Tính $M_4 = \max_{x \in [a,b]} |f^{(4)}(x)|$.

2. **Xác định số đoạn $n$:**
   a. Nếu đã biết $n$: đảm bảo $n$ là số chẵn, nếu lẻ thì tăng thêm $1$.
   b. Nếu dùng sai số $\varepsilon$:
      - Tính $n_{min} = \left( \dfrac{M_4 (b-a)^5}{180 \varepsilon} \right)^{1/4}$.
      - Làm tròn lên số nguyên gần nhất.
      - Nếu $n_{min}$ lẻ: $n = \lceil n_{min} \rceil + 1$.
      - Nếu $n_{min}$ chẵn: $n = \lceil n_{min} \rceil$.

3. **Tính độ dài mỗi đoạn:**
   $$h = \frac{b-a}{n}$$

4. **Tạo các điểm chia và tính giá trị hàm:**
   a. Với $k = 0, 1, 2, \ldots, n$:
      - $x_k = a + k \times h$.
      - Tính $y_k = f(x_k)$.

5. **Áp dụng công thức Simpson 1/3 ghép:**
   a. Khởi tạo $S = y_0 + y_n$.
   b. Với $k = 1$ đến $n-1$, bước nhảy $2$ (các chỉ số lẻ):
      - $S = S + 4 \times y_k$.
   c. Với $k = 2$ đến $n-2$, bước nhảy $2$ (các chỉ số chẵn):
      - $S = S + 2 \times y_k$.
   d. $I = \dfrac{h}{3} \times S$.

6. **Trả về** giá trị $I$.
