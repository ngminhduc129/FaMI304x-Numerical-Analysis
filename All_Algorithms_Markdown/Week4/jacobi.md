# Phương pháp lặp Jacobi

## Công thức toán học

Xét hệ phương trình tuyến tính $A x = B$.

Phân rã ma trận $A$ thành $A = D + L + U$, trong đó:
- $D$ là ma trận đường chéo của $A$
- $L$ là ma trận tam giác dưới (không gồm đường chéo)
- $U$ là ma trận tam giác trên (không gồm đường chéo)

Công thức lặp Jacobi:

$$
x^{(k+1)} = D^{-1}\left(B - (L+U)x^{(k)}\right)
$$

Hay dưới dạng ma trận lặp:

$$
T = \text{diag}\left(\frac{1}{a_{11}}, \frac{1}{a_{22}}, \dots, \frac{1}{a_{nn}}\right)
$$

$$
C = I - T A, \quad D = T B
$$

$$
x^{(k+1)} = C x^{(k)} + D
$$

---

## Thuật toán — Trường hợp ma trận chéo trội hàng

**Mục tiêu:** Giải hệ $Ax = B$ bằng Jacobi với ma trận chéo trội hàng.
**Đầu vào:**
- Ma trận $A$ kích thước $n \times n$ (đọc từ file `JCB_input_A1.txt`)
- Vector $B$ kích thước $n$ (đọc từ file `JCB_input_B1.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:** Vector nghiệm $x^{(k)}$ và lịch sử lặp.

**Bước 1:** Đọc dữ liệu đầu vào
   - **Bước 1.1:** Đọc ma trận vuông $A$ từ file `JCB_input_A1.txt`.
   - **Bước 1.2:** Đọc vector vế phải $B$ từ file `JCB_input_B1.txt`.
   - **Bước 1.3:** Kiểm tra kích thước: $A$ là $n \times n$, $B$ là $n$.

**Bước 2:** Kiểm tra tính chéo trội hàng
   - Với mỗi hàng $i = 1, 2, \dots, n$:
     * Tính $S_i = \sum_{j \neq i} |a_{ij}|$.
     * Nếu $|a_{ii}| \le S_i$: ma trận không chéo trội hàng, có thể không hội tụ.

**Bước 3:** Xây dựng ma trận đường chéo nghịch đảo $T$
   - Với mỗi $i = 1, 2, \dots, n$: $t_{ii} = 1 / a_{ii}$.

**Bước 4:** Xây dựng ma trận lặp $C$ và vector hằng $D$
   - **Bước 4.1:** $C = I - T \cdot A$: $c_{ij} = \delta_{ij} - t_{ii} \cdot a_{ij}$.
   - **Bước 4.2:** $D = T \cdot B$: $d_i = t_{ii} \cdot b_i = b_i / a_{ii}$.

**Bước 5:** Xác định hệ số co $q$ (dùng chuẩn max)
   - Với mỗi hàng $i$: $s_i = \sum_{j=1}^{n} |c_{ij}|$.
   - Gán $q = \max_i s_i$, $\alpha = 1$.

**Bước 6:** Tính ngưỡng hội tụ
   $$\text{tol} = \varepsilon \cdot \frac{1 - q}{\alpha \cdot q}$$

**Bước 7:** Khởi tạo vòng lặp
   - $x^{(0)} = \mathbf{0}$, $k = 0$, $\text{history} = [x^{(0)}]$.

**Bước 8:** Thực hiện vòng lặp Jacobi

   **Bước 8.1:** Tính $x^{(k+1)} = C \cdot x^{(k)} + D$.
   - Với mỗi $i$: $\displaystyle x_i^{(k+1)} = \sum_{j=1}^{n} c_{ij} x_j^{(k)} + d_i$.

   **Bước 8.2:** Tính sai số: $\Delta = \| x^{(k+1)} - x^{(k)} \|_\infty = \max_i |x_i^{(k+1)} - x_i^{(k)}|$.

   **Bước 8.3:** Lưu $x^{(k+1)}$ vào $\text{history}$.

   **Bước 8.4:** Nếu $\Delta \le \text{tol}$: chuyển sang Bước 9.

   **Bước 8.5:** $k = k + 1$. Nếu $k \ge k_{\max}$: chuyển sang Bước 9. Ngược lại, quay lại Bước 8.1.

**Bước 9:** Trả về $x^{(k)}$ và $\text{history}$.

---

## Thuật toán — Trường hợp ma trận chéo trội cột

**Mục tiêu:** Giải hệ $Ax = B$ bằng Jacobi với ma trận chéo trội cột.
**Đầu vào:**
- Ma trận $A$ kích thước $n \times n$ (đọc từ file `JCB_input_A1.txt`)
- Vector $B$ kích thước $n$ (đọc từ file `JCB_input_B1.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:** Vector nghiệm $x^{(k)}$ và lịch sử lặp.

**Bước 1:** Đọc dữ liệu đầu vào (giống bước 1 của trường hợp chéo trội hàng).

**Bước 2:** Kiểm tra tính chéo trội cột
   - Với mỗi cột $j = 1, 2, \dots, n$:
     * Tính $S_j = \sum_{i \neq j} |a_{ij}|$.
     * Nếu $|a_{jj}| \le S_j$: ma trận không chéo trội cột.

**Bước 3:** Xây dựng ma trận $T$, $C$, $D$ theo biến đổi cột
   - **Bước 3.1:** $T = \text{diag}(1/a_{11}, 1/a_{22}, \dots, 1/a_{nn})$.
   - **Bước 3.2:** $C = I - A \cdot T$: $c_{ij} = \delta_{ij} - a_{ij} \cdot t_{jj}$.
   - **Bước 3.3:** $D = B$ (không nhân với $T$).

**Bước 4:** Xác định hệ số co $q$ (dùng chuẩn L1)
   - Với mỗi cột $j$: $s_j = \sum_{i=1}^{n} |c_{ij}|$.
   - Gán $q = \max_j s_j$.
   - Tính $\displaystyle \alpha = \frac{\max_i |a_{ii}|}{\min_i |a_{ii}|}$.

**Bước 5:** Tính ngưỡng hội tụ
   $$\text{tol} = \varepsilon \cdot \frac{1 - q}{\alpha \cdot q}$$

**Bước 6:** Khởi tạo vòng lặp
   - $y^{(0)} = \mathbf{0}$, $k = 0$, $\text{history} = []$.

**Bước 7:** Thực hiện vòng lặp Jacobi (dạng $y$)

   **Bước 7.1:** Tính $y^{(k+1)} = C \cdot y^{(k)} + D$.

   **Bước 7.2:** Tính $x^{(k+1)} = T \cdot y^{(k+1)}$: $x_i^{(k+1)} = y_i^{(k+1)} / a_{ii}$.

   **Bước 7.3:** Tính sai số: $\Delta = \|x^{(k+1)} - x^{(k)}\|$ (chuẩn L1).

   **Bước 7.4:** Lưu $x^{(k+1)}$ vào $\text{history}$.

   **Bước 7.5:** Nếu $\Delta \le \text{tol}$: chuyển sang Bước 8.

   **Bước 7.6:** $k = k + 1$. Nếu $k \ge k_{\max}$: chuyển sang Bước 8.

**Bước 8:** Trả về $x^{(k)}$ và $\text{history}$.
