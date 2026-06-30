# Phương pháp Lũy thừa và Xuống thang (Power Method & Deflation)

## Công thức toán học

### Phương pháp lũy thừa (Power Method) — Case 1: Trị riêng trội thực, bội đơn

$$y^{(k)} = A \cdot x^{(k-1)}$$

$$\lambda^{(k)} = \|y^{(k)}\|$$

$$x^{(k)} = \frac{y^{(k)}}{\lambda^{(k)}}$$

Khi $k \to \infty$: $\lambda^{(k)} \to \lambda_1$ (trị riêng trội) và $x^{(k)} \to v_1$ (vector riêng tương ứng).

### Case 2: Hai trị riêng trội bằng trị tuyệt đối, trái dấu ($\lambda_2 = -\lambda_1$)

Sử dụng lũy thừa bậc chẵn của $A$:

$$A^2 = A \cdot A$$

$$y_{even} = A^2 \cdot y_{even}$$

$$\lambda^2 = \frac{y_{even2}[j]}{y_{even}[j]}$$

$$\lambda_1 = \sqrt{\lambda^2}, \quad \lambda_2 = -\lambda_1$$

$$x_1 = y_{even} + \lambda_1 \cdot y_{odd}, \quad x_2 = y_{even} - \lambda_1 \cdot y_{odd}$$

### Case 3: Hai trị riêng trội là phức liên hợp ($\lambda_2 = \overline{\lambda_1}$)

Giải hệ phương trình từ dãy lũy thừa:

$$y_{m+2} + p \cdot y_{m+1} + q \cdot y_m = 0$$

Giải phương trình bậc hai:

$$z^2 + p z + q = 0 \Rightarrow \lambda_{1,2} = \frac{-p \pm \sqrt{p^2 - 4q}}{2}$$

### Xuống thang (Deflation)

Cho cặp trị riêng-vector riêng $(\lambda, v)$, tìm vector riêng trái $w$:

$$A^T w = \lambda w, \quad w^T v = 1$$

Ma trận mới sau khi xuống thang:

$$A_{new} = A - \lambda \cdot v \cdot w^T$$

## Điều kiện áp dụng từng Case

| Case | Điều kiện | Mô tả |
|------|-----------|-------|
| 1 | $|\lambda_1| > |\lambda_2|$ | Một trị riêng trội nhất, thực, đơn |
| 2 | $\lambda_2 = -\lambda_1$, $|\lambda_1| = |\lambda_2|$ | Hai trị riêng đối nhau, cùng modulus lớn nhất |
| 3 | $\lambda_2 = \overline{\lambda_1}$ | Hai trị riêng phức liên hợp, cùng modulus lớn nhất |

## Thuật toán chi tiết

### power_method_case1($A$, $x_0$, $tol$)

**Mục tiêu:** Tìm trị riêng trội nhất $\lambda_1$ và vector riêng $v_1$ tương ứng.

**Input:** Ma trận $A$ ($n \times n$), vector khởi tạo $x_0$, ngưỡng hội tụ $tol$.

**Output:** $(\lambda_1, v_1)$.

**Các bước:**

1. Chuẩn hóa vector khởi tạo:
   $$x = \frac{x_0}{\|x_0\|}$$

2. Khởi tạo $\lambda = 0$.

3. **Lặp** cho đến khi hội tụ:
   - Tính $y = A \cdot x$
   - $\lambda_{new} = \|y\|$
   - $x_{new} = \dfrac{y}{\lambda_{new}}$
   - **Nếu** $|\lambda_{new} - \lambda| < tol$: **dừng** vòng lặp.
   - **Ngược lại**: gán $x = x_{new}$, $\lambda = \lambda_{new}$ và tiếp tục lặp.

4. **Trả về** $(\lambda, x)$

---

### power_method_case2($A$, $x_0$, $tol$)

**Mục tiêu:** Tìm 2 trị riêng $\lambda_1, \lambda_2$ với $\lambda_2 = -\lambda_1$ và vector riêng tương ứng.

**Input:** Ma trận $A$ ($n \times n$), vector khởi tạo $x_0$, ngưỡng hội tụ $tol$.

**Output:** $(\lambda_1, x_1, \lambda_2, x_2)$.

**Các bước:**

1. Tính $A_2 = A \cdot A$ (bình phương ma trận $A$).

2. Khởi tạo $y_{even} = \dfrac{x_0}{\|x_0\|}$.

3. **Lặp** cho đến khi hội tụ:
   - $y_{even\_old} = y_{even}$
   - $y_{even} = A_2 \cdot y_{even}$
   - Chuẩn hóa $y_{even} = \dfrac{y_{even}}{\|y_{even}\|}$
   - **Nếu** $\|y_{even} - y_{even\_old}\| < tol$: **dừng** vòng lặp.
   - **Ngược lại**: tiếp tục lặp.

4. Tính $\lambda^2$:
   - $y_{even2} = A_2 \cdot y_{even}$
   - Chọn chỉ số $j$ là vị trí phần tử có trị tuyệt đối lớn nhất trong $y_{even}$.
   - $\lambda^2 = \dfrac{y_{even2}[j]}{y_{even}[j]}$

5. Suy ra 2 trị riêng:
   - $\lambda_1 = \sqrt{\lambda^2}$
   - $\lambda_2 = -\lambda_1$

6. Tính 2 vector riêng tương ứng:
   - $y_{odd} = A \cdot y_{even}$
   - $x_1 = y_{even} + \lambda_1 \cdot y_{odd}$, chuẩn hóa $x_1 = \dfrac{x_1}{\|x_1\|}$
   - $x_2 = y_{even} - \lambda_1 \cdot y_{odd}$, chuẩn hóa $x_2 = \dfrac{x_2}{\|x_2\|}$

7. **Trả về** $(\lambda_1, x_1, \lambda_2, x_2)$

---

### power_method_case3($A$, $x_0$, $tol$)

**Mục tiêu:** Tìm 2 trị riêng phức liên hợp $\lambda_1, \lambda_2 = \overline{\lambda_1}$ và vector riêng tương ứng.

**Input:** Ma trận $A$ ($n \times n$), vector khởi tạo $x_0$, ngưỡng hội tụ $tol$.

**Output:** $(\lambda_1, x_1, \lambda_2, x_2)$.

**Các bước:**

1. Khởi tạo dãy lũy thừa:
   - $y_m = \dfrac{x_0}{\|x_0\|}$
   - $y_{m+1} = A \cdot y_m$, chuẩn hóa
   - $y_{m+2} = A \cdot y_{m+1}$, chuẩn hóa

2. Tìm $p, q$ sơ bộ từ 2 chỉ số $i, j$ bất kỳ:
   - Lập hệ:
     $$y_{m+2}[i] + p \cdot y_{m+1}[i] + q \cdot y_m[i] = 0$$
     $$y_{m+2}[j] + p \cdot y_{m+1}[j] + q \cdot y_m[j] = 0$$
   - Giải tìm $p_0, q_0$.

3. **Lặp** cho đến khi hội tụ:
   - Chọn 2 chỉ số $i, j$ (thường là phần tử lớn nhất).
   - Lập hệ $2 \times 2$ từ 2 chỉ số đó:
     $$\begin{cases}
     p \cdot y_{m+1}[i] + q \cdot y_m[i] = -y_{m+2}[i] \\
     p \cdot y_{m+1}[j] + q \cdot y_m[j] = -y_{m+2}[j]
     \end{cases}$$
   - Giải hệ tìm $p, q$.
   - **Nếu** $|p - p_{old}| < tol$ và $|q - q_{old}| < tol$: **dừng** vòng lặp.
   - **Ngược lại**:
     - $y_m = y_{m+1}$
     - $y_{m+1} = y_{m+2}$
     - $y_{m+2} = A \cdot y_{m+1}$, chuẩn hóa
     - Tiếp tục lặp.

4. Giải phương trình bậc 2 tìm trị riêng:
   $$z^2 + p z + q = 0$$
   $$\lambda_1 = \frac{-p + \sqrt{p^2 - 4q}}{2}, \quad \lambda_2 = \frac{-p - \sqrt{p^2 - 4q}}{2}$$

5. Tính vector riêng:
   - $x_1 = y_{m+1} - \lambda_2 \cdot y_m$, chuẩn hóa
   - $x_2 = y_{m+1} - \lambda_1 \cdot y_m$, chuẩn hóa

6. **Trả về** $(\lambda_1, x_1, \lambda_2, x_2)$

---

### deflate_once($A$, $\lambda$, $v$)

**Mục tiêu:** Loại bỏ cặp trị riêng $(\lambda, v)$ khỏi ma trận $A$ để tìm trị riêng tiếp theo.

**Input:** Ma trận $A$, trị riêng $\lambda$, vector riêng $v$.

**Output:** Ma trận $A_{new}$ không còn chứa $\lambda$.

**Các bước:**

1. Tìm vector riêng trái $w$ bằng cách giải hệ:
   $$A^T w = \lambda w, \quad w^T v = 1$$
   - Giải $(A^T - \lambda I)w = 0$ tìm $w$ (sai khác hằng số).
   - Chuẩn hóa: $w = \dfrac{w}{w^T v}$ để thỏa $w^T v = 1$.

2. Tính ma trận mới:
   $$A_{new} = A - \lambda \cdot v \cdot w^T$$

3. **Trả về** $A_{new}$

---

### compute_all_eigenpairs($A$, $x_0$, $tol$)

**Mục tiêu:** Tìm tất cả các cặp trị riêng - vector riêng của ma trận $A$.

**Input:** Ma trận $A$ ($n \times n$), vector khởi tạo $x_0$, ngưỡng hội tụ $tol$.

**Output:** Danh sách các cặp $(\lambda_i, v_i)$ sắp xếp theo $|\lambda_i|$ giảm dần.

**Các bước:**

1. $M = A$ (bản sao của ma trận $A$)
2. $y = x_0$
3. **Lặp** với $k = 1$ đến $n$:
   - Gọi $(\lambda, v) = \text{power\_method}(M, y, tol)$ (chọn case phù hợp)
   - Thêm $(\lambda, v)$ vào danh sách kết quả
   - **Nếu** $k < n$:
     - $M = \text{deflate\_once}(M, \lambda, v)$
     - Đặt lại $y = x_0$
4. Sắp xếp danh sách theo $|\lambda|$ giảm dần.
5. **Trả về** danh sách các cặp $(\lambda_i, v_i)$
