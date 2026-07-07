# Phương pháp Lũy thừa và Xuống thang (Power Method & Deflation)

## Công thức toán học

### Phương pháp lũy thừa — Case 1: Trị riêng trội thực, bội đơn

$$y^{(k)} = A \cdot x^{(k-1)}$$

$$\lambda^{(k)} = \|y^{(k)}\|$$

$$x^{(k)} = \frac{y^{(k)}}{\lambda^{(k)}}$$

Khi $k \to \infty$: $\lambda^{(k)} \to \lambda_1$, $x^{(k)} \to v_1$.

### Case 2: Hai trị riêng trội bằng trị tuyệt đối, trái dấu ($\lambda_2 = -\lambda_1$)

$$A_2 = A \cdot A$$

$$\lambda^2 = \frac{y_{even2}[j]}{y_{even}[j]}$$

$$\lambda_1 = \sqrt{\lambda^2}, \quad \lambda_2 = -\lambda_1$$

### Case 3: Hai trị riêng trội là phức liên hợp

Giải $y_{m+2} + p \cdot y_{m+1} + q \cdot y_m = 0$:

$$z^2 + p z + q = 0 \Rightarrow \lambda_{1,2} = \frac{-p \pm \sqrt{p^2 - 4q}}{2}$$

### Xuống thang (Deflation)

Cho $(\lambda, v)$, tìm vector riêng trái $w$:

$$A^T w = \lambda w, \quad w^T v = 1$$

$$A_{new} = A - \lambda \cdot v \cdot w^T$$

---

## Thuật toán `power_method_case1`

**Mục tiêu:** Tìm trị riêng trội nhất $\lambda_1$ và vector riêng $v_1$.
**Đầu vào:** Ma trận $A$ ($n \times n$), vector khởi tạo $x_0$, ngưỡng hội tụ $tol$.
**Đầu ra:** $(\lambda_1, v_1)$.

**Bước 1:** Chuẩn hóa vector khởi tạo: $x = \dfrac{x_0}{\|x_0\|}$.
**Bước 2:** Khởi tạo $\lambda = 0$.
**Bước 3:** Lặp cho đến khi hội tụ
   - **Bước 3.1:** Tính $y = A \cdot x$.
   - **Bước 3.2:** $\lambda_{new} = \|y\|$.
   - **Bước 3.3:** $x_{new} = \dfrac{y}{\lambda_{new}}$.
   - **Bước 3.4:** Nếu $|\lambda_{new} - \lambda| < tol$: dừng vòng lặp.
   - **Bước 3.5:** Gán $x = x_{new}$, $\lambda = \lambda_{new}$ và tiếp tục lặp.
**Bước 4:** Trả về $(\lambda, x)$.

---

## Thuật toán `power_method_case2`

**Mục tiêu:** Tìm 2 trị riêng $\lambda_1, \lambda_2$ với $\lambda_2 = -\lambda_1$.
**Đầu vào:** Ma trận $A$ ($n \times n$), vector khởi tạo $x_0$, ngưỡng hội tụ $tol$.
**Đầu ra:** $(\lambda_1, x_1, \lambda_2, x_2)$.

**Bước 1:** Tính $A_2 = A \cdot A$.
**Bước 2:** Khởi tạo $y_{even} = \dfrac{x_0}{\|x_0\|}$.
**Bước 3:** Lặp cho đến khi hội tụ
   - **Bước 3.1:** $y_{even\_old} = y_{even}$.
   - **Bước 3.2:** $y_{even} = A_2 \cdot y_{even}$.
   - **Bước 3.3:** Chuẩn hóa $y_{even} = \dfrac{y_{even}}{\|y_{even}\|}$.
   - **Bước 3.4:** Nếu $\|y_{even} - y_{even\_old}\| < tol$: dừng vòng lặp.
**Bước 4:** Tính $\lambda^2$
   - $y_{even2} = A_2 \cdot y_{even}$.
   - Chọn chỉ số $j$ là vị trí phần tử lớn nhất trong $y_{even}$.
   - $\lambda^2 = \dfrac{y_{even2}[j]}{y_{even}[j]}$.
**Bước 5:** Suy ra 2 trị riêng: $\lambda_1 = \sqrt{\lambda^2}$, $\lambda_2 = -\lambda_1$.
**Bước 6:** Tính 2 vector riêng
   - $y_{odd} = A \cdot y_{even}$.
   - $x_1 = y_{even} + \lambda_1 \cdot y_{odd}$, chuẩn hóa.
   - $x_2 = y_{even} - \lambda_1 \cdot y_{odd}$, chuẩn hóa.
**Bước 7:** Trả về $(\lambda_1, x_1, \lambda_2, x_2)$.

---

## Thuật toán `power_method_case3`

**Mục tiêu:** Tìm 2 trị riêng phức liên hợp $\lambda_1, \lambda_2 = \overline{\lambda_1}$.
**Đầu vào:** Ma trận $A$ ($n \times n$), vector khởi tạo $x_0$, ngưỡng hội tụ $tol$.
**Đầu ra:** $(\lambda_1, x_1, \lambda_2, x_2)$.

**Bước 1:** Khởi tạo dãy lũy thừa
   - $y_m = \dfrac{x_0}{\|x_0\|}$.
   - $y_{m+1} = A \cdot y_m$, chuẩn hóa.
   - $y_{m+2} = A \cdot y_{m+1}$, chuẩn hóa.

**Bước 2:** Tìm $p, q$ sơ bộ từ 2 chỉ số $i, j$ bất kỳ
   - Lập hệ:
     $$y_{m+2}[i] + p \cdot y_{m+1}[i] + q \cdot y_m[i] = 0$$
     $$y_{m+2}[j] + p \cdot y_{m+1}[j] + q \cdot y_m[j] = 0$$
   - Giải tìm $p_0, q_0$.

**Bước 3:** Lặp cho đến khi hội tụ
   - **Bước 3.1:** Chọn 2 chỉ số $i, j$ (thường là phần tử lớn nhất).
   - **Bước 3.2:** Lập hệ $2 \times 2$ và giải tìm $p, q$.
   - **Bước 3.3:** Nếu $|p - p_{old}| < tol$ và $|q - q_{old}| < tol$: dừng.
   - **Bước 3.4:** Cập nhật $y_m = y_{m+1}$, $y_{m+1} = y_{m+2}$, $y_{m+2} = A \cdot y_{m+1}$ (chuẩn hóa).

**Bước 4:** Giải phương trình bậc 2 tìm trị riêng
   $$z^2 + p z + q = 0$$
   $$\lambda_1 = \frac{-p + \sqrt{p^2 - 4q}}{2}, \quad \lambda_2 = \frac{-p - \sqrt{p^2 - 4q}}{2}$$

**Bước 5:** Tính vector riêng
   - $x_1 = y_{m+1} - \lambda_2 \cdot y_m$, chuẩn hóa.
   - $x_2 = y_{m+1} - \lambda_1 \cdot y_m$, chuẩn hóa.

**Bước 6:** Trả về $(\lambda_1, x_1, \lambda_2, x_2)$.

---

## Thuật toán `deflate_once`

**Mục tiêu:** Loại bỏ cặp trị riêng $(\lambda, v)$ khỏi ma trận $A$.
**Đầu vào:** Ma trận $A$, trị riêng $\lambda$, vector riêng $v$.
**Đầu ra:** Ma trận $A_{new}$ không còn chứa $\lambda$.

**Bước 1:** Tìm vector riêng trái $w$
   - Giải $(A^T - \lambda I)w = 0$ tìm $w$.
   - Chuẩn hóa: $w = \dfrac{w}{w^T v}$ để thỏa $w^T v = 1$.

**Bước 2:** Tính $A_{new} = A - \lambda \cdot v \cdot w^T$.
**Bước 3:** Trả về $A_{new}$.

---

## Thuật toán `compute_all_eigenpairs`

**Mục tiêu:** Tìm tất cả các cặp trị riêng - vector riêng của ma trận $A$.
**Đầu vào:** Ma trận $A$ ($n \times n$), vector khởi tạo $x_0$, ngưỡng hội tụ $tol$.
**Đầu ra:** Danh sách các cặp $(\lambda_i, v_i)$ sắp xếp theo $|\lambda_i|$ giảm dần.

**Bước 1:** $M = A$ (bản sao của $A$).
**Bước 2:** $y = x_0$.
**Bước 3:** Lặp với $k = 1$ đến $n$
   - **Bước 3.1:** Gọi $(\lambda, v) = \text{power\_method}(M, y, tol)$ (chọn case phù hợp).
   - **Bước 3.2:** Thêm $(\lambda, v)$ vào danh sách kết quả.
   - **Bước 3.3:** Nếu $k < n$: $M = \text{deflate\_once}(M, \lambda, v)$; đặt lại $y = x_0$.
**Bước 4:** Sắp xếp danh sách theo $|\lambda|$ giảm dần.
**Bước 5:** Trả về danh sách các cặp $(\lambda_i, v_i)$.
