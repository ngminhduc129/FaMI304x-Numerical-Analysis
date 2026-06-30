# Phương pháp tiếp tuyến (Newton-Raphson)

## Mô tả

Phương pháp Newton-Raphson tìm nghiệm của $f(x) = 0$ bằng công thức lặp sử dụng tiếp tuyến của đồ thị tại điểm hiện tại.

## Công thức

$$
x_{n+1} = x_n - \frac{f(x_n)}{f'(x_n)}
$$

Sai số tiên nghiệm:

$$
\Delta_x = \frac{|f(x_n)|}{m_1}
$$

với $m_1 = \min|f'(x)|$, $M_2 = \max|f''(x)|$ trên $[a, b]$.

Sai số hậu nghiệm:

$$
\Delta_x = \frac{M_2}{2 m_1} (x_{n+1} - x_n)^2
$$

## Chọn giá trị khởi đầu

Nếu $f(a) \cdot f''(a) > 0$ thì $x_0 = a$; nếu $f(b) \cdot f''(b) > 0$ thì $x_0 = b$.

## Các biến thể

### v1 — Sai số tiên nghiệm $\Delta_x = |f(x_n)| / m_1$

**Thuật toán:**

**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$, khoảng $[a, b]$ chứa nghiệm, số bước $n$, số chữ số làm tròn $rbl$,
$m_1 = \min|f'(x)|$ và $M_2 = \max|f''(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ của $f(x) = 0$.

1. Chọn giá trị khởi đầu $x_0$:
   - Tính $f(a)$ và $f''(a)$.
   - Nếu $f(a) \cdot f''(a) > 0$: gán $x = a$.
   - Nếu $f(b) \cdot f''(b) > 0$: gán $x = b$.

2. Với mỗi bước $i$ từ $0$ đến $n-1$:
   a. Tính giá trị tiếp theo: $x_{\text{new}} = x - \dfrac{f(x)}{f'(x)}$.
   b. Tính sai số tiên nghiệm: $\Delta_x = \dfrac{|f(x_{\text{new}})|}{m_1}$.
   c. Lưu kết quả bước $i$ vào bảng $(x, f(x), f'(x), x_{\text{new}}, f(x_{\text{new}}), \Delta_x)$.
   d. Gán $x = x_{\text{new}}$.
   e. Nếu $f(x_{\text{new}}) = 0$: thoát vòng lặp.

3. In bảng kết quả.

4. Làm tròn:
   - Nếu $rbl = \text{None}$: in $x$.
   - Nếu $rbl$ có giá trị: in $\text{round}(x, rbl)$.

### v2 — Sai số hậu nghiệm $\Delta_x = \dfrac{M_2}{2 m_1} (x_{n+1} - x_n)^2$

**Thuật toán:**

**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$ và $f''(x)$, khoảng $[a, b]$ chứa nghiệm, số bước $n$, số chữ số làm tròn $rbl$,
$m_1 = \min|f'(x)|$ và $M_2 = \max|f''(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ của $f(x) = 0$.

1. Chọn giá trị khởi đầu $x_0$:
   - Nếu $f(a) \cdot f''(a) > 0$: gán $x = a$.
   - Nếu $f(b) \cdot f''(b) > 0$: gán $x = b$.

2. Với mỗi bước $i$ từ $0$ đến $n-1$:
   a. Tính giá trị tiếp theo: $x_{\text{new}} = x - \dfrac{f(x)}{f'(x)}$.
   b. Tính sai số hậu nghiệm: $\Delta_x = \dfrac{M_2}{2 \cdot m_1} \cdot (x_{\text{new}} - x)^2$.
   c. Lưu kết quả bước $i$ vào bảng $(x, f(x), f'(x), x_{\text{new}}, f(x_{\text{new}}), \Delta_x)$.
   d. Gán $x = x_{\text{new}}$.
   e. Nếu $f(x_{\text{new}}) = 0$: thoát vòng lặp.

3. In bảng kết quả.

4. Làm tròn:
   - Nếu $rbl = \text{None}$: in $x$.
   - Nếu $rbl$ có giá trị: in $\text{round}(x, rbl)$.

### v3 — Dừng theo sai số tuyệt đối

**Thuật toán (dùng công thức v1 — $\Delta_x = |f(x_n)| / m_1$):**

**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$, khoảng $[a, b]$, ngưỡng $\varepsilon > 0$,
$m_1 = \min|f'(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ của $f(x) = 0$.

1. Chọn $x_0$ dựa vào dấu $f \cdot f''$.
2. Lặp:
   a. $x_{\text{new}} = x - f(x) / f'(x)$.
   b. Nếu $|f(x_{\text{new}})| < m_1 \cdot \varepsilon$: thoát vòng lặp.
   c. Gán $x = x_{\text{new}}$.
3. In $\text{round}(x, rbl)$.

**Thuật toán (dùng công thức v2 — $\Delta_x = \dfrac{M_2}{2 m_1} (x_{n+1} - x_n)^2$):**

**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$ và $f''(x)$, khoảng $[a, b]$, ngưỡng $\varepsilon > 0$,
$m_1 = \min|f'(x)|$, $M_2 = \max|f''(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ của $f(x) = 0$.

1. Chọn $x_0$ dựa vào dấu $f \cdot f''$.
2. Lặp:
   a. $x_{\text{new}} = x - f(x) / f'(x)$.
   b. Nếu $|x_{\text{new}} - x| < \sqrt{\dfrac{2 \cdot m_1 \cdot \varepsilon}{M_2}}$: thoát vòng lặp.
   c. Gán $x = x_{\text{new}}$.
3. In $\text{round}(x, rbl)$.

### v4 — Dừng theo sai số tương đối

**Thuật toán:**

**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$ và $f''(x)$, khoảng $[a, b]$, ngưỡng $\eta > 0$,
$m_1 = \min|f'(x)|$, $M_2 = \max|f''(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ của $f(x) = 0$.

1. Chọn $x_0$ dựa vào dấu $f \cdot f''$.
2. Lặp:
   a. $x_{\text{new}} = x - f(x) / f'(x)$.
   b. Tính sai số tương đối: $\sigma_x = \dfrac{(x_{\text{new}} - x)^2}{|x_{\text{new}}|}$.
   c. Nếu $\sigma_x < \dfrac{2 \cdot \eta \cdot m_1}{M_2}$: thoát vòng lặp.
   d. Gán $x = x_{\text{new}}$.
3. In $\text{round}(x, rbl)$.
