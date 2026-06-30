# Phương pháp lặp Jacobi

## Công thức toán học

Xét hệ phương trình tuyến tính:

$$
A x = B
$$

Phân rã ma trận $A$ thành:

$$
A = D + L + U
$$

trong đó:
- $D$ là ma trận đường chéo của $A$
- $L$ là ma trận tam giác dưới (không gồm đường chéo)
- $U$ là ma trận tam giác trên (không gồm đường chéo)

Công thức lặp Jacobi:

$$
x^{(k+1)} = D^{-1}\left(B - (L+U)x^{(k)}\right)
$$

Hay viết dưới dạng ma trận lặp:

$$
T = \text{diag}\left(\frac{1}{a_{11}}, \frac{1}{a_{22}}, \dots, \frac{1}{a_{nn}}\right)
$$

$$
C = I - T A, \quad D = T B
$$

$$
x^{(k+1)} = C x^{(k)} + D
$$

Với ma trận $A$ có **chéo trội hàng**, chuẩn sử dụng là chuẩn max ($\infty$-norm):

$$
q = \max_i \sum_{j} |c_{ij}|
$$

Ngưỡng hội tụ:

$$
\text{tol} = \varepsilon \cdot \frac{1-q}{\alpha \cdot q}
$$

trong đó $\alpha = 1$ cho trường hợp chéo trội hàng.

Với ma trận $A$ có **chéo trội cột**, biến đổi khác:

$$
T = \text{diag}\left(\frac{1}{a_{11}}, \dots, \frac{1}{a_{nn}}\right), \quad
C = I - A T, \quad D = B
$$

$$
y^{(k+1)} = C y^{(k)} + D, \quad x^{(k+1)} = T y^{(k+1)}
$$

Chuẩn sử dụng là chuẩn L1 (tổng trị tuyệt đối):

$$
q = \max_j \sum_{i} |c_{ij}|, \quad
\alpha = \frac{\max_i |a_{ii}|}{\min_i |a_{ii}|}
$$

---

## Thuật toán — Trường hợp ma trận chéo trội hàng

**Đầu vào:**
- Ma trận $A$ kích thước $n \times n$ (đọc từ file `JCB_input_A1.txt`)
- Vector $B$ kích thước $n$ (đọc từ file `JCB_input_B1.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:**
- Vector nghiệm xấp xỉ $x^{(k)}$
- DataFrame lịch sử lặp (danh sách các vector $x^{(0)}, x^{(1)}, \dots, x^{(k)}$)

**Các bước thực hiện:**

1. **Đọc dữ liệu đầu vào:**
   a. Đọc ma trận vuông $A$ từ file `JCB_input_A1.txt`.
   b. Đọc vector vế phải $B$ từ file `JCB_input_B1.txt`.
   c. Kiểm tra kích thước: $A$ có kích thước $n \times n$, $B$ có kích thước $n$.

2. **Kiểm tra tính chéo trội hàng của ma trận $A$:**
   a. Với mỗi hàng $i = 1, 2, \dots, n$:
      - Tính tổng trị tuyệt đối các phần tử ngoài đường chéo: $S_i = \sum_{j \neq i} |a_{ij}|$.
      - So sánh: nếu $|a_{ii}| \le S_i$ với bất kỳ $i$ nào, ma trận không chéo trội hàng.
   b. Nếu không chéo trội hàng, phương pháp có thể không hội tụ.

3. **Xây dựng ma trận đường chéo nghịch đảo $T$:**
   a. Với mỗi $i = 1, 2, \dots, n$: $t_{ii} = 1 / a_{ii}$.
   b. Các phần tử ngoài đường chéo của $T$ bằng $0$.
   c. Dạng ma trận: $T = \text{diag}(1/a_{11}, 1/a_{22}, \dots, 1/a_{nn})$.

4. **Xây dựng ma trận lặp $C$ và vector hằng $D$:**
   a. Tính ma trận lặp: $C = I - T \cdot A$.
      - Với mỗi $i, j$: $c_{ij} = \delta_{ij} - t_{ii} \cdot a_{ij}$, trong đó $\delta_{ij}$ là ký hiệu Kronecker ($\delta_{ij} = 1$ nếu $i=j$, $0$ nếu $i \neq j$).
   b. Tính vector hằng: $D = T \cdot B$.
      - Với mỗi $i$: $d_i = t_{ii} \cdot b_i = b_i / a_{ii}$.

5. **Xác định hệ số co $q$ (dùng chuẩn max):**
   a. Với mỗi hàng $i = 1, 2, \dots, n$:
      - Tính tổng trị tuyệt đối: $s_i = \sum_{j=1}^{n} |c_{ij}|$.
   b. Gán $q = \max_{i} s_i$.
   c. Đặt $\alpha = 1$ (cho trường hợp chéo trội hàng).

6. **Tính ngưỡng hội tụ:**
   $$
   \text{tol} = \varepsilon \cdot \frac{1 - q}{\alpha \cdot q}
   $$

7. **Khởi tạo vòng lặp:**
   a. Gán vector xấp xỉ ban đầu $x^{(0)} = \mathbf{0}$ (vector không), hoặc dùng giá trị đọc từ file nếu có.
   b. Gán $k = 0$.
   c. Khởi tạo lịch sử lặp: `history = [x^{(0)}]`.

8. **Thực hiện vòng lặp Jacobi:** Với $k = 0, 1, 2, \dots$ cho đến khi hội tụ hoặc đạt $k_{\max}$.
   a. **Tính xấp xỉ mới:** $x^{(k+1)} = C \cdot x^{(k)} + D$.
      - Với mỗi thành phần $i$: $\displaystyle x_i^{(k+1)} = \sum_{j=1}^{n} c_{ij} x_j^{(k)} + d_i$.
   b. **Tính sai số tuyệt đối giữa hai bước lặp:**
      $$
      \Delta = \left\| x^{(k+1)} - x^{(k)} \right\|_\infty = \max_{i} \left| x_i^{(k+1)} - x_i^{(k)} \right|
      $$
   c. **Lưu lịch sử:** Thêm $x^{(k+1)}$ vào `history`.
   d. **Kiểm tra điều kiện dừng:**
      - Nếu $\Delta \le \text{tol}$: chuyển sang bước 9.
   e. **Cập nhật:** Gán $k = k + 1$.
      - Nếu $k \ge k_{\max}$: chuyển sang bước 9.

9. **Trả về kết quả:**
   a. Vector nghiệm xấp xỉ cuối cùng $x^{(k)}$.
   b. DataFrame lịch sử lặp (mỗi hàng là một vector $x^{(k)}$).

---

## Thuật toán — Trường hợp ma trận chéo trội cột

**Đầu vào:**
- Ma trận $A$ kích thước $n \times n$ (đọc từ file `JCB_input_A1.txt`)
- Vector $B$ kích thước $n$ (đọc từ file `JCB_input_B1.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:**
- Vector nghiệm xấp xỉ $x^{(k)}$
- DataFrame lịch sử lặp

**Các bước thực hiện:**

1. **Đọc dữ liệu đầu vào:** (giống bước 1 của trường hợp chéo trội hàng)

2. **Kiểm tra tính chéo trội cột của ma trận $A$:**
   a. Với mỗi cột $j = 1, 2, \dots, n$:
      - Tính tổng trị tuyệt đối các phần tử ngoài đường chéo: $S_j = \sum_{i \neq j} |a_{ij}|$.
      - So sánh: nếu $|a_{jj}| \le S_j$ với bất kỳ $j$ nào, ma trận không chéo trội cột.

3. **Xây dựng ma trận $T$, $C$, $D$ theo biến đổi cột:**
   a. $T = \text{diag}(1/a_{11}, 1/a_{22}, \dots, 1/a_{nn})$.
   b. Tính ma trận lặp: $C = I - A \cdot T$.
      - Với mỗi $i, j$: $c_{ij} = \delta_{ij} - a_{ij} \cdot t_{jj}$.
   c. Gán vector hằng: $D = B$ (không nhân với $T$).

4. **Xác định hệ số co $q$ (dùng chuẩn L1 - chuẩn cột):**
   a. Với mỗi cột $j = 1, 2, \dots, n$:
      - Tính tổng trị tuyệt đối theo cột: $s_j = \sum_{i=1}^{n} |c_{ij}|$.
   b. Gán $q = \max_{j} s_j$.
   c. Tính $\displaystyle \alpha = \frac{\max_i |a_{ii}|}{\min_i |a_{ii}|}$.

5. **Tính ngưỡng hội tụ:**
   $$
   \text{tol} = \varepsilon \cdot \frac{1 - q}{\alpha \cdot q}
   $$

6. **Khởi tạo vòng lặp:**
   a. Gán $y^{(0)} = \mathbf{0}$.
   b. Gán $k = 0$.
   c. Khởi tạo lịch sử lặp: `history = []`.

7. **Thực hiện vòng lặp Jacobi (dạng $y$):** Với $k = 0, 1, 2, \dots$
   a. **Tính $y$ mới:** $y^{(k+1)} = C \cdot y^{(k)} + D$.
      - Với mỗi thành phần $i$: $\displaystyle y_i^{(k+1)} = \sum_{j=1}^{n} c_{ij} y_j^{(k)} + d_i$.
   b. **Tính $x$ tương ứng:** $x^{(k+1)} = T \cdot y^{(k+1)}$.
      - Với mỗi $i$: $x_i^{(k+1)} = t_{ii} \cdot y_i^{(k+1)} = y_i^{(k+1)} / a_{ii}$.
   c. **Tính sai số:** $\Delta = \|x^{(k+1)} - x^{(k)}\|$ (dùng chuẩn L1).
   d. **Lưu lịch sử:** Thêm $x^{(k+1)}$ vào `history`.
   e. **Kiểm tra điều kiện dừng:**
      - Nếu $\Delta \le \text{tol}$: chuyển sang bước 8.
   f. **Cập nhật:** $k = k + 1$.

8. **Trả về kết quả:**
   a. Vector nghiệm xấp xỉ cuối cùng $x^{(k)}$.
   b. DataFrame lịch sử lặp.
