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

## Các biến thể

### v1 — Sai số tiên nghiệm $\Delta_n = \dfrac{q^n}{1-q} \cdot |x_1 - x_0|$

**Thuật toán:**

**Đầu vào:** Hàm $\varphi(x)$, khoảng $[a, b]$ chứa nghiệm, giá trị khởi đầu $x_0$, số bước $n$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ của $x = \varphi(x)$.

1. Tính hằng số Lipschitz:
   - Tính $q = \max|\varphi'(x)|$ trên $[a, b]$.
   - Nếu $q \geq 1$: dừng, thông báo "$q \geq 1$, phương pháp không hội tụ".

2. Khởi tạo:
   - Gán $x = x_0$.
   - Gán $\text{dx} = 0$ (biến lưu $|x_1 - x_0|$).

3. Với mỗi bước $i$ từ $0$ đến $n-1$:
   a. Tính giá trị lặp: $x_{\text{new}} = \varphi(x)$.
   b. Nếu $i = 0$: gán $\text{dx} = |x_{\text{new}} - x|$ (lưu $|x_1 - x_0|$).
   c. Tính sai số tiên nghiệm: $\Delta = \dfrac{q^{\,i+1}}{1 - q} \cdot \text{dx}$.
   d. Lưu kết quả bước $i$ vào bảng $(x, \varphi(x), x_{\text{new}}, \Delta)$.
   e. Gán $x = x_{\text{new}}$.

4. In bảng kết quả.

5. Làm tròn:
   - Nếu $rbl = \text{None}$: in $x$.
   - Nếu $rbl$ có giá trị: in $\text{round}(x, rbl)$.

### v2 — Sai số hậu nghiệm $\Delta_n = \dfrac{q}{1-q} \cdot |x_{n+1} - x_n|$

**Thuật toán:**

**Đầu vào:** Hàm $\varphi(x)$, khoảng $[a, b]$ chứa nghiệm, giá trị khởi đầu $x_0$, số bước $n$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ của $x = \varphi(x)$.

1. Tính hằng số Lipschitz:
   - Tính $q = \max|\varphi'(x)|$ trên $[a, b]$.
   - Nếu $q \geq 1$: dừng, thông báo "$q \geq 1$, phương pháp không hội tụ".

2. Khởi tạo:
   - Gán $x = x_0$.

3. Với mỗi bước $i$ từ $0$ đến $n-1$:
   a. Tính giá trị lặp: $x_{\text{new}} = \varphi(x)$.
   b. Tính sai số hậu nghiệm: $\Delta = \dfrac{q}{1 - q} \cdot |x_{\text{new}} - x|$.
   c. Lưu kết quả bước $i$ vào bảng $(x, \varphi(x), x_{\text{new}}, \Delta)$.
   d. Gán $x = x_{\text{new}}$.

4. In bảng kết quả.

5. Làm tròn:
   - Nếu $rbl = \text{None}$: in $x$.
   - Nếu $rbl$ có giá trị: in $\text{round}(x, rbl)$.

### v3 — Dừng theo sai số tuyệt đối

**Thuật toán:**

**Đầu vào:** Hàm $\varphi(x)$, khoảng $[a, b]$ chứa nghiệm, giá trị khởi đầu $x_0$, ngưỡng $\varepsilon > 0$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ của $x = \varphi(x)$.

1. Tính $q = \max|\varphi'(x)|$ trên $[a, b]$. Yêu cầu $q < 1$.
2. Gán $x = x_0$.
3. Lặp:
   a. $x_{\text{new}} = \varphi(x)$.
   b. Tính $\Delta_x = |x_{\text{new}} - x|$.
   c. Nếu $\Delta_x < \dfrac{1 - q}{q} \cdot \varepsilon$: thoát vòng lặp.
   d. Gán $x = x_{\text{new}}$.
4. In $\text{round}(x, rbl)$.

### v4 — Dừng theo sai số tương đối

**Thuật toán:**

**Đầu vào:** Hàm $\varphi(x)$, khoảng $[a, b]$ chứa nghiệm, giá trị khởi đầu $x_0$, ngưỡng $\eta > 0$, số chữ số làm tròn $rbl$.
**Đầu ra:** Nghiệm $x$ xấp xỉ của $x = \varphi(x)$.

1. Tính $q = \max|\varphi'(x)|$ trên $[a, b]$. Yêu cầu $q < 1$.
2. Gán $x = x_0$.
3. Lặp:
   a. $x_{\text{new}} = \varphi(x)$.
   b. Tính sai số tương đối: $\sigma_x = \dfrac{|x_{\text{new}} - x|}{|x_{\text{new}}|}$.
   c. Nếu $\sigma_x < \dfrac{1 - q}{q} \cdot \eta$: thoát vòng lặp.
   d. Gán $x = x_{\text{new}}$.
4. In $\text{round}(x, rbl)$.
