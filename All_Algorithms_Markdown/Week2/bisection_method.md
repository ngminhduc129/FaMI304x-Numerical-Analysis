# Phương pháp chia đôi (Bisection Method)

## Mô tả

Phương pháp chia đôi tìm nghiệm của phương trình $f(x) = 0$ trên $[a, b]$ với $f(a) \cdot f(b) < 0$ bằng cách liên tục chia đôi khoảng và chọn khoảng con có chứa nghiệm.

## Công thức

$$
x_{n+1} = \frac{a_n + b_n}{2}
$$

Sai số sau $n$ bước:

$$
\Delta = \frac{b - a}{2^{n+1}}
$$

---

## v1 — Sai số dùng công thức $\Delta = (b-a)/2^{n+1}$

**Mục tiêu:** Tìm nghiệm xấp xỉ của $f(x)=0$ trên $[a,b]$ với số bước cho trước.
**Đầu vào:** Hàm $f(x)$, khoảng $[a, b]$ với $f(a) \cdot f(b) < 0$, số bước $n$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Kiểm tra điều kiện đầu vào
- Tính $f(a)$ và $f(b)$.
- Nếu $f(a) \cdot f(b) \geq 0$: dừng, thông báo "Bạn chưa chọn đúng $a$ và $b$".

**Bước 2:** Khởi tạo
- Gán $x = 0$.
- Gán $\text{sign}_a = 1$ nếu $f(a) > 0$, $-1$ nếu $f(a) < 0$, $0$ nếu $f(a) = 0$.
- Gán $\text{diff} = b - a$.

**Bước 3:** Lặp với $i = 0, 1, \dots, n-1$

   **Bước 3.1:** Tính điểm chia đôi: $x_{\text{new}} = \dfrac{a + b}{2}$.
   
   **Bước 3.2:** Tính sai số: $\Delta = \dfrac{\text{diff}}{2^{\,i+1}}$.
   
   **Bước 3.3:** Lưu kết quả bước $i$ vào bảng $(a, b, x_{\text{new}}, f(a), f(b), f(x_{\text{new}}), \Delta)$.
   
   **Bước 3.4:** Gán $x = x_{\text{new}}$.
   
   **Bước 3.5:** Nếu $f(x_{\text{new}}) = 0$: thoát khỏi vòng lặp.
   
   **Bước 3.6:** Xác định khoảng con chứa nghiệm:
   - Nếu $f(x_{\text{new}}) \cdot \text{sign}_a < 0$: gán $b = x_{\text{new}}$ (nghiệm nằm nửa trái).
   - Nếu $f(x_{\text{new}}) \cdot \text{sign}_a > 0$: gán $a = x_{\text{new}}$ (nghiệm nằm nửa phải).

**Bước 4:** In bảng kết quả.

**Bước 5:** Làm tròn và xuất kết quả
- Nếu $rbl = \text{None}$: in $x$.
- Nếu $rbl$ có giá trị:
  - Tính $\text{total\_delta} = \Delta + 0.5 \times 10^{-rbl}$.
  - In $\text{round}(x, rbl)$ và $\text{total\_delta}$.

---

## v2 — Sai số dùng $\Delta_x = |x_{n+1} - x_n|$

**Mục tiêu:** Tìm nghiệm xấp xỉ của $f(x)=0$ trên $[a,b]$ với số bước cho trước, sai số tính từ chênh lệch hai lần lặp.
**Đầu vào:** Hàm $f(x)$, khoảng $[a, b]$ với $f(a) \cdot f(b) < 0$, số bước $n$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Kiểm tra điều kiện đầu vào
- Tính $f(a)$ và $f(b)$.
- Nếu $f(a) \cdot f(b) \geq 0$: dừng, thông báo "Bạn chưa chọn đúng $a$ và $b$".

**Bước 2:** Khởi tạo
- Gán $x = 0$.
- Gán $\text{sign}_a = 1$ nếu $f(a) > 0$, $-1$ nếu $f(a) < 0$, $0$ nếu $f(a) = 0$.

**Bước 3:** Lặp với $i = 0, 1, \dots, n-1$

   **Bước 3.1:** Tính điểm chia đôi: $x_{\text{new}} = \dfrac{a + b}{2}$.
   
   **Bước 3.2:** Tính sai số: $\Delta_x = |x_{\text{new}} - x|$.
   
   **Bước 3.3:** Lưu kết quả bước $i$ vào bảng $(a, b, x_{\text{new}}, f(a), f(b), f(x_{\text{new}}), \Delta_x)$.
   
   **Bước 3.4:** Gán $x = x_{\text{new}}$.
   
   **Bước 3.5:** Nếu $f(x_{\text{new}}) = 0$: thoát khỏi vòng lặp.
   
   **Bước 3.6:** Xác định khoảng con chứa nghiệm:
   - Nếu $f(x_{\text{new}}) \cdot \text{sign}_a < 0$: gán $b = x_{\text{new}}$.
   - Nếu $f(x_{\text{new}}) \cdot \text{sign}_a > 0$: gán $a = x_{\text{new}}$.

**Bước 4:** In bảng kết quả.

**Bước 5:** Làm tròn và xuất kết quả
- Nếu $rbl = \text{None}$: in $x$.
- Nếu $rbl$ có giá trị: in $\text{round}(x, rbl)$.

---

## v3 — Sai số tuyệt đối (dừng khi $\Delta_x < \varepsilon$)

**Mục tiêu:** Tìm nghiệm xấp xỉ của $f(x)=0$ với độ chính xác tuyệt đối $\varepsilon$.
**Đầu vào:** Hàm $f(x)$, khoảng $[a, b]$ với $f(a) \cdot f(b) < 0$, ngưỡng $\varepsilon > 0$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Kiểm tra điều kiện đầu vào
- Tính $f(a)$ và $f(b)$.
- Nếu $f(a) \cdot f(b) \geq 0$: dừng, thông báo "Bạn chưa chọn đúng $a$ và $b$".

**Bước 2:** Khởi tạo
- Gán $x = 0$.
- Gán $\text{sign}_a = 1$ nếu $f(a) > 0$, $-1$ nếu $f(a) < 0$, $0$ nếu $f(a) = 0$.
- Gán $\text{step} = 0$.

**Bước 3:** Lặp cho đến khi hội tụ

   **Bước 3.1:** Tính điểm chia đôi: $x_{\text{new}} = \dfrac{a + b}{2}$.
   
   **Bước 3.2:** Tính sai số: $\Delta_x = |x_{\text{new}} - x|$.
   
   **Bước 3.3:** Lưu kết quả bước $\text{step}$ vào bảng $(a, b, x_{\text{new}}, f(a), f(b), f(x_{\text{new}}), \Delta_x)$.
   
   **Bước 3.4:** Gán $x = x_{\text{new}}$.
   
   **Bước 3.5:** Kiểm tra điều kiện dừng:
   - Nếu $f(x_{\text{new}}) = 0$: thoát vòng lặp.
   - Nếu $\Delta_x < \varepsilon$: thoát vòng lặp.
   
   **Bước 3.6:** Xác định khoảng con chứa nghiệm:
   - Nếu $f(x_{\text{new}}) \cdot \text{sign}_a < 0$: gán $b = x_{\text{new}}$.
   - Nếu $f(x_{\text{new}}) \cdot \text{sign}_a > 0$: gán $a = x_{\text{new}}$.
   
   **Bước 3.7:** Tăng $\text{step}$ lên $1$.

**Bước 4:** In bảng kết quả.

**Bước 5:** Làm tròn: in $\text{round}(x, rbl)$ nếu $rbl$ có giá trị, ngược lại in $x$.

---

## v4 — Sai số tương đối (dừng khi $\sigma_x < \eta$)

**Mục tiêu:** Tìm nghiệm xấp xỉ của $f(x)=0$ với độ chính xác tương đối $\eta$.
**Đầu vào:** Hàm $f(x)$, khoảng $[a, b]$ với $f(a) \cdot f(b) < 0$, ngưỡng $\eta > 0$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Kiểm tra điều kiện đầu vào
- Tính $f(a)$ và $f(b)$.
- Nếu $f(a) \cdot f(b) \geq 0$: dừng, thông báo "Bạn chưa chọn đúng $a$ và $b$".

**Bước 2:** Khởi tạo
- Gán $x = 0$.
- Gán $\text{sign}_a = 1$ nếu $f(a) > 0$, $-1$ nếu $f(a) < 0$, $0$ nếu $f(a) = 0$.
- Gán $\text{step} = 0$.

**Bước 3:** Lặp cho đến khi hội tụ

   **Bước 3.1:** Tính điểm chia đôi: $x_{\text{new}} = \dfrac{a + b}{2}$.
   
   **Bước 3.2:** Tính sai số tương đối: $\sigma_x = \dfrac{|x_{\text{new}} - x|}{|x_{\text{new}}|}$.
   
   **Bước 3.3:** Lưu kết quả bước $\text{step}$ vào bảng $(a, b, x_{\text{new}}, f(a), f(b), f(x_{\text{new}}), \sigma_x)$.
   
   **Bước 3.4:** Gán $x = x_{\text{new}}$.
   
   **Bước 3.5:** Kiểm tra điều kiện dừng:
   - Nếu $f(x_{\text{new}}) = 0$: thoát vòng lặp.
   - Nếu $\sigma_x < \eta$: thoát vòng lặp.
   
   **Bước 3.6:** Xác định khoảng con chứa nghiệm:
   - Nếu $f(x_{\text{new}}) \cdot \text{sign}_a < 0$: gán $b = x_{\text{new}}$.
   - Nếu $f(x_{\text{new}}) \cdot \text{sign}_a > 0$: gán $a = x_{\text{new}}$.
   
   **Bước 3.7:** Tăng $\text{step}$ lên $1$.

**Bước 4:** In bảng kết quả.

**Bước 5:** Làm tròn: in $\text{round}(x, rbl)$ nếu $rbl$ có giá trị, ngược lại in $x$.
