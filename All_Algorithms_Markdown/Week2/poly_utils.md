# Tiện ích đa thức (Polynomial Utilities)

## Mô tả

Cung cấp các hàm xử lý đa thức: tính giá trị, đạo hàm, chặn nghiệm thực/phức, tìm cực trị.

---

## Hàm `polynomial(coeffs, x)`

Tính giá trị đa thức tại $x$:

$$
P(x) = c_0 + c_1 x + c_2 x^2 + \dots + c_n x^n
$$

**Thuật toán (Horner):**

**Mục tiêu:** Tính $P(x)$ bằng sơ đồ Horner.
**Đầu vào:** Mảng hệ số $\text{coeffs} = [c_0, c_1, \dots, c_n]$, giá trị $x$.
**Đầu ra:** $P(x)$.

**Bước 1:** Khởi tạo $\text{result} = c_n$ (hệ số bậc cao nhất).
**Bước 2:** Với mỗi $i$ từ $n-1$ xuống $0$:
   - $\text{result} = \text{result} \cdot x + c_i$.
**Bước 3:** Trả về $\text{result}$.

---

## Hàm `polynomial_derivative(coeffs, x)`

Tính đạo hàm $P'(x)$:

$$
P'(x) = c_1 + 2 c_2 x + 3 c_3 x^2 + \dots + n c_n x^{n-1}
$$

**Thuật toán:**

**Mục tiêu:** Tính $P'(x)$.
**Đầu vào:** Mảng hệ số $\text{coeffs} = [c_0, c_1, \dots, c_n]$, giá trị $x$.
**Đầu ra:** $P'(x)$.

**Bước 1:** Khởi tạo $\text{result} = 0$.
**Bước 2:** Với mỗi $i$ từ $1$ đến $n$:
   - $\text{result} = \text{result} + i \cdot c_i \cdot x^{\,i-1}$.
**Bước 3:** Trả về $\text{result}$.

---

## Hàm `polynomial_second_derivative(coeffs, x)`

Tính đạo hàm cấp hai $P''(x)$.

**Thuật toán:**

**Mục tiêu:** Tính $P''(x)$.
**Đầu vào:** Mảng hệ số $\text{coeffs} = [c_0, c_1, \dots, c_n]$, giá trị $x$.
**Đầu ra:** $P''(x)$.

**Bước 1:** Khởi tạo $\text{result} = 0$.
**Bước 2:** Với mỗi $i$ từ $2$ đến $n$:
   - $\text{result} = \text{result} + i \cdot (i-1) \cdot c_i \cdot x^{\,i-2}$.
**Bước 3:** Trả về $\text{result}$.

---

## Hàm `complex_radius(coeffs)`

Bán kính chặn nghiệm phức:

$$
R = 1 + \frac{\max(|c_0|, |c_1|, \dots, |c_{n-1}|)}{|c_n|}
$$

Mọi nghiệm phức đều thỏa $|z| \leq R$.

**Thuật toán:**

**Mục tiêu:** Tính bán kính chặn nghiệm phức.
**Đầu vào:** Mảng hệ số $\text{coeffs} = [c_0, c_1, \dots, c_n]$ với $c_n \neq 0$.
**Đầu ra:** Bán kính $R$.

**Bước 1:** Xác định bậc đa thức: $n = \text{len(coeffs)} - 1$.
**Bước 2:** Lấy hệ số bậc cao nhất: $c_n = \text{coeffs}[n]$.
**Bước 3:** Nếu $c_n = 0$: báo lỗi "Hệ số bậc cao nhất không được bằng $0$".
**Bước 4:** Tìm giá trị tuyệt đối lớn nhất của các hệ số từ $c_0$ đến $c_{n-1}$:
   - $\text{max\_abs} = \max(|c_0|, |c_1|, \dots, |c_{n-1}|)$.
**Bước 5:** Tính $R = 1 + \dfrac{\text{max\_abs}}{|c_n|}$.
**Bước 6:** Trả về $R$.

---

## Hàm `real_radius(coeffs)`

Bán kính chặn nghiệm thực dương.

Nếu $c_k$ là hệ số âm có chỉ số lớn nhất, $B$ là giá trị tuyệt đối lớn nhất trong các hệ số âm:

$$
R = 1 + \left(\frac{B}{|c_n|}\right)^{1/(n-k)}
$$

Nếu không có hệ số âm, trả về $\text{None}$.

**Thuật toán:**

**Mục tiêu:** Tính bán kính chặn nghiệm thực dương.
**Đầu vào:** Mảng hệ số $\text{coeffs} = [c_0, c_1, \dots, c_n]$ với $c_n \neq 0$.
**Đầu ra:** Bán kính $R$ hoặc $\text{None}$.

**Bước 1:** Xác định bậc đa thức: $n = \text{len(coeffs)} - 1$.
**Bước 2:** Lấy $c_n = \text{coeffs}[n]$.
**Bước 3:** Duyệt qua các hệ số $c_0, c_1, \dots, c_{n-1}$:
   - Lọc ra các hệ số âm.
   - Nếu không có hệ số âm nào: trả về $\text{None}$.
**Bước 4:** Tìm $B$ là giá trị tuyệt đối lớn nhất trong các hệ số âm.
**Bước 5:** Tìm $k$ là chỉ số lớn nhất của hệ số âm có giá trị tuyệt đối bằng $B$.
**Bước 6:** Tính $R = 1 + \left(\dfrac{B}{|c_n|}\right)^{1/(n-k)}$.
**Bước 7:** Trả về $R$.

---

## Hàm `find_extrema(coeffs, lower_bound, upper_bound, lr, tolerance)`

Tìm cực trị địa phương bằng gradient descent/ascent với learning rate thích ứng.

**Thuật toán:**

**Mục tiêu:** Tìm các điểm cực trị của đa thức trên $[\text{lower\_bound}, \text{upper\_bound}]$.
**Đầu vào:** Hệ số đa thức $\text{coeffs} = [c_0, c_1, \dots, c_n]$, cận dưới $\text{lower\_bound}$, cận trên $\text{upper\_bound}$, tốc độ học $\text{lr}$, ngưỡng hội tụ $\text{tolerance}$.
**Đầu ra:** Danh sách các điểm cực trị $(x, P(x), \text{loại})$.

**Bước 1:** Khởi tạo tập điểm xuất phát
- Chia đoạn $[\text{lower\_bound}, \text{upper\_bound}]$ thành $N$ điểm cách đều.
- Gọi $\text{points}$ là danh sách các điểm xuất phát này.

**Bước 2:** Với mỗi điểm $p$ trong $\text{points}$:

   **Bước 2.1:** Gán $x = p$.
   
   **Bước 2.2:** Lặp cho đến khi hội tụ
   - Tính gradient $g = P'(x)$ bằng $\text{polynomial\_derivative}$.
   - Nếu tìm cực đại: $\text{direction} = +1$.
   - Nếu tìm cực tiểu: $\text{direction} = -1$.
   - Tính $x_{\text{new}} = x + \text{direction} \cdot \text{lr} \cdot g$.
   - Nếu $P(x_{\text{new}})$ không cải thiện so với $P(x)$:
     * Giảm tốc độ học: $\text{lr} = \text{lr} / 2$.
   - Nếu $|x_{\text{new}} - x| < \text{tolerance}$:
     * Thêm $x$ vào danh sách điểm dừng.
     * Thoát vòng lặp.
   - Gán $x = x_{\text{new}}$.

**Bước 3:** Xác định loại cực trị
- Với mỗi điểm $x$ trong danh sách điểm dừng:
  - Tính $P''(x)$ bằng $\text{polynomial\_second\_derivative}$.
  - Nếu $P''(x) > 0$: phân loại là "cực tiểu".
  - Nếu $P''(x) < 0$: phân loại là "cực đại".
  - Nếu $P''(x) = 0$: phân loại là "không xác định".

**Bước 4:** Hậu xử lý
- Loại bỏ các điểm trùng lặp: nếu $|x_i - x_j| < \text{tolerance}$, giữ lại một điểm.
- Sắp xếp danh sách theo $x$ tăng dần.

**Bước 5:** Trả về danh sách các điểm cực trị $(x, P(x), \text{loại})$.
