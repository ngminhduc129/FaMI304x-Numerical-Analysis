# Phương pháp dây cung (Secant Method)

## Mô tả

Phương pháp dây cung tìm nghiệm của $f(x) = 0$ bằng cách thay đường cong $f(x)$ bằng dây cung qua hai điểm $(a, f(a))$ và $(b, f(b))$, sau đó tìm giao điểm của dây cung với trục $x$.

## Công thức

$$
x_{n+1} = \frac{b \cdot f(a) - a \cdot f(b)}{f(a) - f(b)}
$$

Sai số tiên nghiệm:

$$
\Delta_x = \frac{|f(x_n)|}{m_1}
$$

với $m_1 = \min|f'(x)|$ trên $[a, b]$, $M_1 = \max|f'(x)|$ trên $[a, b]$.

Sai số hậu nghiệm:

$$
\Delta_x = \frac{M_1 - m_1}{m_1} |x_{n+1} - x_n|
$$

---

## v1 — Sai số tiên nghiệm $\Delta_x = |f(x_n)| / m_1$

**Mục tiêu:** Tìm nghiệm xấp xỉ của $f(x)=0$ bằng phương pháp dây cung với số bước cho trước.
**Đầu vào:** Hàm $f(x)$, khoảng $[a, b]$ với $f(a) \cdot f(b) < 0$, số bước $n$, số chữ số làm tròn $rbl$, $m_1 = \min|f'(x)|$ và $M_1 = \max|f'(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Kiểm tra điều kiện đầu vào
- Tính $f(a)$ và $f(b)$.
- Nếu $f(a) \cdot f(b) \geq 0$: dừng, thông báo "Bạn chưa chọn đúng $a$ và $b$".
- Nếu $M_1 \cdot m_1 \leq 0$: dừng (đạo hàm đổi dấu hoặc bằng $0$).

**Bước 2:** Khởi tạo
- Gán $x = a$.
- Gán $\text{mrk} = b$.
- Gán $\text{sign}_a = 1$ nếu $f(a) > 0$, $-1$ nếu $f(a) < 0$.

**Bước 3:** Lặp với $i = 0, 1, \dots, n-1$

   **Bước 3.1:** Tính điểm cắt dây cung:
   $$x_{\text{new}} = \dfrac{\text{mrk} \cdot f(x) - x \cdot f(\text{mrk})}{f(x) - f(\text{mrk})}$$
   
   **Bước 3.2:** Tính sai số tiên nghiệm: $\Delta_x = \dfrac{|f(x_{\text{new}})|}{m_1}$.
   
   **Bước 3.3:** Nếu $i = 0$: xác định cặp điểm ban đầu
   - Nếu $f(a) \cdot f(x_{\text{new}}) < 0$: gán $\text{mrk} = x$ (nghiệm giữa $a$ và $x_{\text{new}}$).
   - Nếu $f(a) \cdot f(x_{\text{new}}) > 0$: giữ nguyên $\text{mrk} = b$ (nghiệm giữa $x_{\text{new}}$ và $b$).
   
   **Bước 3.4:** Lưu kết quả bước $i$ vào bảng $(x, \text{mrk}, x_{\text{new}}, f(x), f(\text{mrk}), f(x_{\text{new}}), \Delta_x)$.
   
   **Bước 3.5:** Gán $x = x_{\text{new}}$.
   
   **Bước 3.6:** Nếu $f(x_{\text{new}}) = 0$: thoát vòng lặp.

**Bước 4:** In bảng kết quả.

**Bước 5:** Làm tròn và xuất kết quả
- Nếu $rbl = \text{None}$: in $x$.
- Nếu $rbl$ có giá trị: in $\text{round}(x, rbl)$.

---

## v2 — Sai số hậu nghiệm $\Delta_x = \dfrac{M_1 - m_1}{m_1} |x_{n+1} - x_n|$

**Mục tiêu:** Tìm nghiệm xấp xỉ với sai số hậu nghiệm.
**Đầu vào:** Hàm $f(x)$, khoảng $[a, b]$ với $f(a) \cdot f(b) < 0$, số bước $n$, số chữ số làm tròn $rbl$, $m_1 = \min|f'(x)|$ và $M_1 = \max|f'(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Kiểm tra điều kiện đầu vào
- Tính $f(a)$ và $f(b)$.
- Nếu $f(a) \cdot f(b) \geq 0$: dừng.
- Nếu $M_1 \cdot m_1 \leq 0$: dừng.

**Bước 2:** Khởi tạo
- Gán $x = a$, $\text{mrk} = b$.
- Gán $\text{sign}_a = 1$ nếu $f(a) > 0$, $-1$ nếu $f(a) < 0$.

**Bước 3:** Lặp với $i = 0, 1, \dots, n-1$

   **Bước 3.1:** Tính điểm cắt dây cung:
   $$x_{\text{new}} = \dfrac{\text{mrk} \cdot f(x) - x \cdot f(\text{mrk})}{f(x) - f(\text{mrk})}$$
   
   **Bước 3.2:** Tính sai số hậu nghiệm: $\Delta_x = \dfrac{M_1 - m_1}{m_1} \cdot |x_{\text{new}} - x|$.
   
   **Bước 3.3:** Nếu $i = 0$: xác định cặp điểm ban đầu
   - Nếu $f(a) \cdot f(x_{\text{new}}) < 0$: gán $\text{mrk} = x$.
   - Nếu $f(a) \cdot f(x_{\text{new}}) > 0$: giữ nguyên $\text{mrk} = b$.
   
   **Bước 3.4:** Lưu kết quả bước $i$ vào bảng.
   
   **Bước 3.5:** Gán $x = x_{\text{new}}$.
   
   **Bước 3.6:** Nếu $f(x_{\text{new}}) = 0$: thoát vòng lặp.

**Bước 4:** In bảng kết quả.

**Bước 5:** Làm tròn và xuất kết quả.

---

## v3 — Dừng theo sai số tuyệt đối

### Cách 1: Dùng công thức v1 ($\Delta_x = |f(x_n)| / m_1$)

**Mục tiêu:** Tìm nghiệm với độ chính xác tuyệt đối $\varepsilon$.
**Đầu vào:** Hàm $f(x)$, khoảng $[a, b]$ với $f(a) \cdot f(b) < 0$, ngưỡng $\varepsilon > 0$, $m_1 = \min|f'(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Kiểm tra $f(a) \cdot f(b) < 0$ và $m_1 > 0$.
**Bước 2:** Gán $x = a$, $\text{mrk} = b$.
**Bước 3:** Lặp
   - **Bước 3.1:** Tính $x_{\text{new}} = \dfrac{\text{mrk} \cdot f(x) - x \cdot f(\text{mrk})}{f(x) - f(\text{mrk})}$.
   - **Bước 3.2:** Nếu $|f(x_{\text{new}})| < m_1 \cdot \varepsilon$: thoát vòng lặp.
   - **Bước 3.3:** Nếu $i = 0$: cập nhật $\text{mrk}$ theo dấu $f(a) \cdot f(x_{\text{new}})$.
   - **Bước 3.4:** Gán $x = x_{\text{new}}$.
**Bước 4:** In kết quả $\text{round}(x, rbl)$.

### Cách 2: Dùng công thức v2 ($\Delta_x = \dfrac{M_1 - m_1}{m_1} |x_{n+1} - x_n|$)

**Mục tiêu:** Tìm nghiệm với độ chính xác tuyệt đối $\varepsilon$.
**Đầu vào:** Hàm $f(x)$, khoảng $[a, b]$ với $f(a) \cdot f(b) < 0$, ngưỡng $\varepsilon > 0$, $m_1 = \min|f'(x)|$, $M_1 = \max|f'(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Kiểm tra $f(a) \cdot f(b) < 0$, $M_1 \cdot m_1 > 0$.
**Bước 2:** Gán $x = a$, $\text{mrk} = b$.
**Bước 3:** Lặp
   - **Bước 3.1:** Tính $x_{\text{new}} = \dfrac{\text{mrk} \cdot f(x) - x \cdot f(\text{mrk})}{f(x) - f(\text{mrk})}$.
   - **Bước 3.2:** Tính $\Delta_x = \dfrac{M_1 - m_1}{m_1} \cdot |x_{\text{new}} - x|$.
   - **Bước 3.3:** Nếu $\Delta_x < \varepsilon$: thoát vòng lặp.
   - **Bước 3.4:** Nếu $i = 0$: cập nhật $\text{mrk}$ theo dấu $f(a) \cdot f(x_{\text{new}})$.
   - **Bước 3.5:** Gán $x = x_{\text{new}}$.
**Bước 4:** In kết quả $\text{round}(x, rbl)$.

---

## v4 — Dừng theo sai số tương đối

**Mục tiêu:** Tìm nghiệm với độ chính xác tương đối $\eta$.
**Đầu vào:** Hàm $f(x)$, khoảng $[a, b]$ với $f(a) \cdot f(b) < 0$, ngưỡng $\eta > 0$, $m_1 = \min|f'(x)|$, $M_1 = \max|f'(x)|$ trên $[a, b]$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Kiểm tra $f(a) \cdot f(b) < 0$, $M_1 \cdot m_1 > 0$.
**Bước 2:** Gán $x = a$, $\text{mrk} = b$.
**Bước 3:** Lặp
   - **Bước 3.1:** Tính $x_{\text{new}} = \dfrac{\text{mrk} \cdot f(x) - x \cdot f(\text{mrk})}{f(x) - f(\text{mrk})}$.
   - **Bước 3.2:** Tính sai số tương đối: $\sigma_x = \dfrac{|x_{\text{new}} - x|}{|x_{\text{new}}|}$.
   - **Bước 3.3:** Nếu $\sigma_x < \dfrac{\eta \cdot m_1}{M_1 - m_1}$: thoát vòng lặp.
   - **Bước 3.4:** Nếu $i = 0$: cập nhật $\text{mrk}$ theo dấu $f(a) \cdot f(x_{\text{new}})$.
   - **Bước 3.5:** Gán $x = x_{\text{new}}$.
**Bước 4:** In kết quả $\text{round}(x, rbl)$.
