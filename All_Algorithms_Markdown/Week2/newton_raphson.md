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

**Bước 0:** Chọn $x_0$
- Nếu $f(a) \cdot f''(a) > 0$ thì $x_0 = a$.
- Nếu $f(b) \cdot f''(b) > 0$ thì $x_0 = b$.

---

## v1 — Sai số tiên nghiệm $\Delta_x = |f(x_n)| / m_1$

**Mục tiêu:** Tìm nghiệm xấp xỉ của $f(x)=0$ với số bước cho trước, sai số tính theo công thức tiên nghiệm.
**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$, khoảng $[a, b]$ chứa nghiệm, số bước $n$, số chữ số làm tròn $rbl$, $m_1 = \min|f'(x)|$ và $M_2 = \max|f''(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Chọn giá trị khởi đầu $x_0$
- Tính $f(a)$ và $f''(a)$.
- Nếu $f(a) \cdot f''(a) > 0$: gán $x = a$.
- Nếu $f(b) \cdot f''(b) > 0$: gán $x = b$.

**Bước 2:** Lặp với $i = 0, 1, \dots, n-1$

   **Bước 2.1:** Tính giá trị tiếp theo:
   $$x_{\text{new}} = x - \dfrac{f(x)}{f'(x)}$$
   
   **Bước 2.2:** Tính sai số tiên nghiệm: $\Delta_x = \dfrac{|f(x_{\text{new}})|}{m_1}$.
   
   **Bước 2.3:** Lưu kết quả bước $i$ vào bảng $(x, f(x), f'(x), x_{\text{new}}, f(x_{\text{new}}), \Delta_x)$.
   
   **Bước 2.4:** Gán $x = x_{\text{new}}$.
   
   **Bước 2.5:** Nếu $f(x_{\text{new}}) = 0$: thoát vòng lặp.

**Bước 3:** In bảng kết quả.

**Bước 4:** Làm tròn và xuất kết quả
- Nếu $rbl = \text{None}$: in $x$.
- Nếu $rbl$ có giá trị: in $\text{round}(x, rbl)$.

---

## v2 — Sai số hậu nghiệm $\Delta_x = \dfrac{M_2}{2 m_1} (x_{n+1} - x_n)^2$

**Mục tiêu:** Tìm nghiệm xấp xỉ với sai số hậu nghiệm.
**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$ và $f''(x)$, khoảng $[a, b]$ chứa nghiệm, số bước $n$, số chữ số làm tròn $rbl$, $m_1 = \min|f'(x)|$ và $M_2 = \max|f''(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Chọn giá trị khởi đầu $x_0$
- Nếu $f(a) \cdot f''(a) > 0$: gán $x = a$.
- Nếu $f(b) \cdot f''(b) > 0$: gán $x = b$.

**Bước 2:** Lặp với $i = 0, 1, \dots, n-1$

   **Bước 2.1:** Tính giá trị tiếp theo:
   $$x_{\text{new}} = x - \dfrac{f(x)}{f'(x)}$$
   
   **Bước 2.2:** Tính sai số hậu nghiệm:
   $$\Delta_x = \dfrac{M_2}{2 \cdot m_1} \cdot (x_{\text{new}} - x)^2$$
   
   **Bước 2.3:** Lưu kết quả bước $i$ vào bảng.
   
   **Bước 2.4:** Gán $x = x_{\text{new}}$.
   
   **Bước 2.5:** Nếu $f(x_{\text{new}}) = 0$: thoát vòng lặp.

**Bước 3:** In bảng kết quả.

**Bước 4:** Làm tròn và xuất kết quả.

---

## v3 — Dừng theo sai số tuyệt đối

### Cách 1: Dùng công thức v1 ($\Delta_x = |f(x_n)| / m_1$)

**Mục tiêu:** Tìm nghiệm với độ chính xác tuyệt đối $\varepsilon$.
**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$, khoảng $[a, b]$, ngưỡng $\varepsilon > 0$, $m_1 = \min|f'(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Chọn $x_0$ dựa vào dấu $f \cdot f''$.
**Bước 2:** Lặp
   - **Bước 2.1:** $x_{\text{new}} = x - f(x) / f'(x)$.
   - **Bước 2.2:** Nếu $|f(x_{\text{new}})| < m_1 \cdot \varepsilon$: thoát vòng lặp.
   - **Bước 2.3:** Gán $x = x_{\text{new}}$.
**Bước 3:** In $\text{round}(x, rbl)$.

### Cách 2: Dùng công thức v2 ($\Delta_x = \dfrac{M_2}{2 m_1} (x_{n+1} - x_n)^2$)

**Mục tiêu:** Tìm nghiệm với độ chính xác tuyệt đối $\varepsilon$.
**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$ và $f''(x)$, khoảng $[a, b]$, ngưỡng $\varepsilon > 0$, $m_1 = \min|f'(x)|$, $M_2 = \max|f''(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Chọn $x_0$ dựa vào dấu $f \cdot f''$.
**Bước 2:** Lặp
   - **Bước 2.1:** $x_{\text{new}} = x - f(x) / f'(x)$.
   - **Bước 2.2:** Nếu $|x_{\text{new}} - x| < \sqrt{\dfrac{2 \cdot m_1 \cdot \varepsilon}{M_2}}$: thoát vòng lặp.
   - **Bước 2.3:** Gán $x = x_{\text{new}}$.
**Bước 3:** In $\text{round}(x, rbl)$.

---

## v4 — Dừng theo sai số tương đối

**Mục tiêu:** Tìm nghiệm với độ chính xác tương đối $\eta$.
**Đầu vào:** Hàm $f(x)$, đạo hàm $f'(x)$ và $f''(x)$, khoảng $[a, b]$, ngưỡng $\eta > 0$, $m_1 = \min|f'(x)|$, $M_2 = \max|f''(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Chọn $x_0$ dựa vào dấu $f \cdot f''$.
**Bước 2:** Lặp
   - **Bước 2.1:** $x_{\text{new}} = x - f(x) / f'(x)$.
   - **Bước 2.2:** Tính sai số tương đối: $\sigma_x = \dfrac{(x_{\text{new}} - x)^2}{|x_{\text{new}}|}$.
   - **Bước 2.3:** Nếu $\sigma_x < \dfrac{2 \cdot \eta \cdot m_1}{M_2}$: thoát vòng lặp.
   - **Bước 2.4:** Gán $x = x_{\text{new}}$.
**Bước 3:** In $\text{round}(x, rbl)$.
