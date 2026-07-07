# Phương pháp lặp Gauss-Seidel

## Công thức toán học

Xét hệ phương trình tuyến tính $A x = B$.

Đưa về dạng lặp điểm cố định:

$$
T = \text{diag}\left(\frac{1}{a_{11}}, \frac{1}{a_{22}}, \dots, \frac{1}{a_{nn}}\right)
$$

$$
C = I - T A, \quad D = T B
$$

Phân rã ma trận $C$ thành $C = L + U$, trong đó:
- $L$ là ma trận tam giác dưới chặt (strict lower triangular)
- $U$ là ma trận tam giác trên chặt (strict upper triangular)

Công thức lặp Gauss-Seidel (cập nhật từng thành phần):

$$
x_i^{(k+1)} = \sum_{j=1}^{i-1} c_{ij} x_j^{(k+1)} + \sum_{j=i+1}^{n} c_{ij} x_j^{(k)} + d_i, \quad i = 1, 2, \dots, n
$$

Hệ số hội tụ cho trường hợp chéo trội hàng:

$$
q = \max_i \frac{\sum_j |l_{ij}|}{1 - \sum_j |u_{ij}|}
$$

Ngưỡng hội tụ:

$$
\text{tol} = \varepsilon \cdot \frac{(1-q)(1-s)}{q}
$$

---

## Thuật toán — Trường hợp ma trận chéo trội hàng

**Mục tiêu:** Giải hệ $Ax = B$ bằng Gauss-Seidel với ma trận chéo trội hàng.
**Đầu vào:**
- Ma trận $A$ kích thước $n \times n$ (đọc từ file `GS_input_A1.txt`)
- Vector $B$ kích thước $n$ (đọc từ file `GS_input_B1.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:** Vector nghiệm $x^{(k)}$ và lịch sử lặp.

**Bước 1:** Đọc dữ liệu đầu vào
   - **Bước 1.1:** Đọc ma trận vuông $A$ từ file `GS_input_A1.txt`.
   - **Bước 1.2:** Đọc vector vế phải $B$ từ file `GS_input_B1.txt`.
   - **Bước 1.3:** Kiểm tra kích thước.

**Bước 2:** Kiểm tra tính chéo trội của $A$
   - **Bước 2.1:** Chéo trội hàng: với mỗi $i$, kiểm tra $|a_{ii}| > \sum_{j \neq i} |a_{ij}|$.
   - **Bước 2.2:** Chéo trội cột: với mỗi $j$, kiểm tra $|a_{jj}| > \sum_{i \neq j} |a_{ij}|$.

**Bước 3:** Xây dựng ma trận đường chéo nghịch đảo $T$
   - Với mỗi $i = 1, 2, \dots, n$: $t_{ii} = 1 / a_{ii}$.

**Bước 4:** Xây dựng ma trận lặp $C$ và vector hằng $D$
   - **Bước 4.1:** $C = I - T \cdot A$: $c_{ij} = \delta_{ij} - t_{ii} \cdot a_{ij}$.
   - **Bước 4.2:** $D = T \cdot B$: $d_i = t_{ii} \cdot b_i = b_i / a_{ii}$.

**Bước 5:** Tính hệ số hội tụ $q$
   - Với mỗi hàng $i = 1, 2, \dots, n$:
      * $L_i = \sum_{j < i} |c_{ij}|$, $U_i = \sum_{j > i} |c_{ij}|$.
     * $\text{ratio}_i = \dfrac{L_i}{1 - U_i}$.
   - Gán $q = \max_i \text{ratio}_i$, $s = 0$.

**Bước 6:** Tính ngưỡng hội tụ
   $$\text{tol} = \varepsilon \cdot \frac{(1 - q)(1 - s)}{q}$$

**Bước 7:** Khởi tạo vòng lặp
   - $x^{(0)} = \mathbf{0}$, $k = 0$, $\text{history} = [x^{(0)}]$.

**Bước 8:** Thực hiện vòng lặp Gauss-Seidel

   **Bước 8.1:** Sao chép $x^{(k+1)} = x^{(k)}$.

   **Bước 8.2:** Với mỗi $i = 1, 2, \dots, n$:
   - $\displaystyle x_i^{(k+1)} = \sum_{j < i} c_{ij} x_j^{(k+1)} + \sum_{j > i} c_{ij} x_j^{(k)} + d_i$.

   **Bước 8.3:** Tính sai số: $\Delta = \| x^{(k+1)} - x^{(k)} \|_\infty = \max_i |x_i^{(k+1)} - x_i^{(k)}|$.

   **Bước 8.4:** Lưu $x^{(k+1)}$ vào $\text{history}$.

   **Bước 8.5:** Nếu $\Delta \le \text{tol}$: chuyển sang Bước 9.

   **Bước 8.6:** $k = k + 1$. Nếu $k \ge k_{\max}$: chuyển sang Bước 9. Ngược lại, quay lại Bước 8.1.

**Bước 9:** Trả về $x^{(k)}$ và $\text{history}$.

---

## Thuật toán — Trường hợp ma trận chéo trội cột

**Mục tiêu:** Giải hệ $Ax = B$ bằng Gauss-Seidel với ma trận chéo trội cột.
**Đầu vào:**
- Ma trận $A$ kích thước $n \times n$ (đọc từ file `GS_input_A1.txt`)
- Vector $B$ kích thước $n$ (đọc từ file `GS_input_B1.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:** Vector nghiệm $x^{(k)}$ và lịch sử lặp.

**Bước 1:** Đọc dữ liệu đầu vào (giống bước 1 của trường hợp chéo trội hàng).

**Bước 2:** Xây dựng ma trận $T$, $C$, $D$ theo biến đổi cột
   - **Bước 2.1:** $T = \text{diag}(1/a_{11}, 1/a_{22}, \dots, 1/a_{nn})$.
   - **Bước 2.2:** $C = I - A \cdot T$: $c_{ij} = \delta_{ij} - a_{ij} \cdot t_{jj}$.
   - **Bước 2.3:** $D = B$ (không nhân với $T$).

**Bước 3:** Tính hệ số hội tụ $q$ và $s$ (dùng chuẩn cột)
   - Với mỗi cột $j = 1, 2, \dots, n$:
     * $L_j = \sum_{i < j} |c_{ij}|$, $U_j = \sum_{i > j} |c_{ij}|$.
     * $\text{ratio}_j = \dfrac{U_j}{1 - L_j}$.
   - Gán $q = \max_j \text{ratio}_j$.
   - Gán $s = \max_j L_j$.

**Bước 4:** Tính ngưỡng hội tụ
   $$\text{tol} = \varepsilon \cdot \frac{(1 - q)(1 - s)}{q}$$

**Bước 5:** Khởi tạo vòng lặp
   - $y^{(0)} = \mathbf{0}$, $x^{(0)} = \mathbf{0}$, $k = 0$, $\text{history} = [].

**Bước 6:** Thực hiện vòng lặp Gauss-Seidel (dạng $y$)

   **Bước 6.1:** Sao chép $y^{(k+1)} = y^{(k)}$.

   **Bước 6.2:** Với mỗi $i = 1, 2, \dots, n$:
   - $\displaystyle y_i^{(k+1)} = \sum_{j < i} c_{ij} y_j^{(k+1)} + \sum_{j > i} c_{ij} y_j^{(k)} + d_i$.

   **Bước 6.3:** Tính $x^{(k+1)} = T \cdot y^{(k+1)}$: $x_i^{(k+1)} = y_i^{(k+1)} / a_{ii}$.

   **Bước 6.4:** Tính sai số: $\Delta = \|x^{(k+1)} - x^{(k)}\|_\infty$.

   **Bước 6.5:** Lưu $x^{(k+1)}$ vào $\text{history}$.

   **Bước 6.6:** Nếu $\Delta \le \text{tol}$: chuyển sang Bước 7.

   **Bước 6.7:** $k = k + 1$.

**Bước 7:** Trả về $x^{(k)}$ và $\text{history}$.
