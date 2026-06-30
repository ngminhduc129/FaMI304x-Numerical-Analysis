# Công thức hình thang (Trapezoidal Rule)

## Công thức toán học

Xấp xỉ tích phân xác định bằng cách chia đoạn $[a, b]$ thành $n$ đoạn con đều nhau, mỗi đoạn có độ dài $h = \frac{b-a}{n}$.

Công thức hình thang ghép (Composite Trapezoidal Rule):

$$I = \int_a^b f(x)\,dx \approx \frac{h}{2} \left[ f(x_0) + f(x_n) + 2\sum_{k=1}^{n-1} f(x_k) \right]$$

## Sai số

Sai số của công thức hình thang ghép:

$$|E| \leq \frac{M_2 (b-a)^3}{12 n^2}$$

trong đó $M_2 = \max_{x \in [a,b]} |f''(x)|$.

Để đạt sai số $\varepsilon$, cần:

$$n > \sqrt{\frac{M_2 (b-a)^3}{12 \varepsilon}}$$

## Thuật toán

**Đầu vào:** Hàm $f(x)$; khoảng $[a, b]$; sai số cho phép $\varepsilon$ (hoặc số đoạn $n$).

**Đầu ra:** Giá trị xấp xỉ tích phân $I \approx \int_a^b f(x)\,dx$.

1. **Xác định $M_2$ (nếu dùng sai số):**
   a. Tìm $f''(x)$.
   b. Tính $M_2 = \max_{x \in [a,b]} |f''(x)|$ (khảo sát hàm hoặc tính tại các điểm nghi ngờ).

2. **Xác định số đoạn $n$:**
   a. Nếu đã biết $n$: chuyển sang Bước 3.
   b. Nếu dùng sai số $\varepsilon$:
      - Tính $n_{min} = \sqrt{\dfrac{M_2 (b-a)^3}{12 \varepsilon}}$.
      - Làm tròn lên: $n = \lceil n_{min} \rceil$.

3. **Tính độ dài mỗi đoạn:**
   $$h = \frac{b-a}{n}$$

4. **Tạo các điểm chia và tính giá trị hàm:**
   a. Với $k = 0, 1, 2, \ldots, n$:
      - $x_k = a + k \times h$.
      - Tính $y_k = f(x_k)$.

5. **Áp dụng công thức hình thang ghép:**
   a. Khởi tạo $S = y_0 + y_n$.
   b. Với $k = 1$ đến $n-1$:
      - $S = S + 2 \times y_k$.
   c. $I = \dfrac{h}{2} \times S$.

6. **Trả về** giá trị $I$.
