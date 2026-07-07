# Phương pháp lặp điểm cố định (Fixed-Point Iteration)

## Mô tả

Phương pháp lặp điểm cố định tìm nghiệm của phương trình $x = \varphi(x)$ bằng cách xây dựng dãy lặp $x_{n+1} = \varphi(x_n)$ cho đến khi hội tụ.

Với $f(x) = 0$, ta biến đổi về dạng $x = \varphi(x)$.

## Công thức

$$
x_{n+1} = \varphi(x_n)
$$

Sai số tiên nghiệm ($q$ là hằng số Lipschitz, $q < 1$):

$$
\Delta_n = \frac{q^n}{1-q} |x_1 - x_0|
$$

Sai số hậu nghiệm:

$$
\Delta_n = \frac{q}{1-q} |x_{n+1} - x_n|
$$

---

## v1 — Sai số tiên nghiệm $\Delta_n = \dfrac{q^n}{1-q} \cdot |x_1 - x_0|$

**Mục tiêu:** Tìm nghiệm xấp xỉ của $x = \varphi(x)$ với số bước cho trước, sai số tiên nghiệm.
**Đầu vào:** Hàm $\varphi(x)$, khoảng $[a, b]$ chứa nghiệm, giá trị khởi đầu $x_0$, số bước $n$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Tính hằng số Lipschitz
- Tính $q = \max|\varphi'(x)|$ trên $[a, b]$.
- Nếu $q \geq 1$: dừng, thông báo "$q \geq 1$, phương pháp không hội tụ".

**Bước 2:** Khởi tạo
- Gán $x = x_0$.
- Gán $\text{dx} = 0$ (biến lưu $|x_1 - x_0|$).

**Bước 3:** Lặp với $i = 0, 1, \dots, n-1$

   **Bước 3.1:** Tính giá trị lặp: $x_{\text{new}} = \varphi(x)$.
   
   **Bước 3.2:** Nếu $i = 0$: gán $\text{dx} = |x_{\text{new}} - x|$ (lưu $|x_1 - x_0|$).
   
   **Bước 3.3:** Tính sai số tiên nghiệm: $\Delta = \dfrac{q^{\,i+1}}{1 - q} \cdot \text{dx}$.
   
   **Bước 3.4:** Lưu kết quả bước $i$ vào bảng $(x, \varphi(x), x_{\text{new}}, \Delta)$.
   
   **Bước 3.5:** Gán $x = x_{\text{new}}$.

**Bước 4:** In bảng kết quả.

**Bước 5:** Làm tròn
- Nếu $rbl = \text{None}$: in $x$.
- Nếu $rbl$ có giá trị: in $\text{round}(x, rbl)$.

---

## v2 — Sai số hậu nghiệm $\Delta_n = \dfrac{q}{1-q} \cdot |x_{n+1} - x_n|$

**Mục tiêu:** Tìm nghiệm xấp xỉ với sai số hậu nghiệm.
**Đầu vào:** Hàm $\varphi(x)$, khoảng $[a, b]$ chứa nghiệm, giá trị khởi đầu $x_0$, số bước $n$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Tính hằng số Lipschitz
- Tính $q = \max|\varphi'(x)|$ trên $[a, b]$.
- Nếu $q \geq 1$: dừng, thông báo "$q \geq 1$, phương pháp không hội tụ".

**Bước 2:** Khởi tạo: gán $x = x_0$.

**Bước 3:** Lặp với $i = 0, 1, \dots, n-1$

   **Bước 3.1:** Tính giá trị lặp: $x_{\text{new}} = \varphi(x)$.
   
   **Bước 3.2:** Tính sai số hậu nghiệm: $\Delta = \dfrac{q}{1 - q} \cdot |x_{\text{new}} - x|$.
   
   **Bước 3.3:** Lưu kết quả bước $i$ vào bảng $(x, \varphi(x), x_{\text{new}}, \Delta)$.
   
   **Bước 3.4:** Gán $x = x_{\text{new}}$.

**Bước 4:** In bảng kết quả.

**Bước 5:** Làm tròn.

---

## v3 — Dừng theo sai số tuyệt đối

**Mục tiêu:** Tìm nghiệm với độ chính xác tuyệt đối $\varepsilon$.
**Đầu vào:** Hàm $\varphi(x)$, khoảng $[a, b]$ chứa nghiệm, giá trị khởi đầu $x_0$, ngưỡng $\varepsilon > 0$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Tính $q = \max|\varphi'(x)|$ trên $[a, b]$. Yêu cầu $q < 1$.
**Bước 2:** Gán $x = x_0$.
**Bước 3:** Lặp
   - **Bước 3.1:** $x_{\text{new}} = \varphi(x)$.
   - **Bước 3.2:** Tính $\Delta_x = |x_{\text{new}} - x|$.
   - **Bước 3.3:** Nếu $\Delta_x < \dfrac{1 - q}{q} \cdot \varepsilon$: thoát vòng lặp.
   - **Bước 3.4:** Gán $x = x_{\text{new}}$.
**Bước 4:** In $\text{round}(x, rbl)$.

---

## v4 — Dừng theo sai số tương đối

**Mục tiêu:** Tìm nghiệm với độ chính xác tương đối $\eta$.
**Đầu vào:** Hàm $\varphi(x)$, khoảng $[a, b]$ chứa nghiệm, giá trị khởi đầu $x_0$, ngưỡng $\eta > 0$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ.

**Bước 1:** Tính $q = \max|\varphi'(x)|$ trên $[a, b]$. Yêu cầu $q < 1$.
**Bước 2:** Gán $x = x_0$.
**Bước 3:** Lặp
   - **Bước 3.1:** $x_{\text{new}} = \varphi(x)$.
   - **Bước 3.2:** Tính sai số tương đối: $\sigma_x = \dfrac{|x_{\text{new}} - x|}{|x_{\text{new}}|}$.
   - **Bước 3.3:** Nếu $\sigma_x < \dfrac{1 - q}{q} \cdot \eta$: thoát vòng lặp.
   - **Bước 3.4:** Gán $x = x_{\text{new}}$.
**Bước 4:** In $\text{round}(x, rbl)$.
