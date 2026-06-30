# Phương pháp lặp Gauss-Seidel

## Công thức toán học

Xét hệ phương trình tuyến tính:

$$
A x = B
$$

Đưa về dạng lặp điểm cố định:

$$
T = \text{diag}\left(\frac{1}{a_{11}}, \frac{1}{a_{22}}, \dots, \frac{1}{a_{nn}}\right)
$$

$$
C = I - T A, \quad D = T B
$$

Phân rã ma trận $C$ thành:

$$
C = L + U
$$

trong đó:
- $L$ là ma trận tam giác dưới chặt (strict lower triangular)
- $U$ là ma trận tam giác trên chặt (strict upper triangular)

Công thức lặp Gauss-Seidel (cập nhật từng thành phần):

$$
x^{(k+1)} = L x^{(k+1)} + U x^{(k)} + D
$$

Hay viết theo từng thành phần:

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

với $s = 0$ cho trường hợp chéo trội hàng.

Cho trường hợp chéo trội cột, biến đổi:

$$
T = \text{diag}\left(\frac{1}{a_{11}}, \dots, \frac{1}{a_{nn}}\right), \quad
C = I - A T, \quad D = B
$$

$$
q = \max_j \frac{\sum_i |u_{ij}|}{1 - \sum_i |l_{ij}|}, \quad
s = \max_j \sum_i |l_{ij}|
$$

---

## Thuật toán — Trường hợp ma trận chéo trội hàng

**Đầu vào:**
- Ma trận $A$ kích thước $n \times n$ (đọc từ file `GS_input_A1.txt`)
- Vector $B$ kích thước $n$ (đọc từ file `GS_input_B1.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:**
- Vector nghiệm xấp xỉ $x^{(k)}$
- DataFrame lịch sử lặp (danh sách các vector $x^{(0)}, x^{(1)}, \dots, x^{(k)}$)

**Các bước thực hiện:**

1. **Đọc dữ liệu đầu vào:**
   a. Đọc ma trận vuông $A$ từ file `GS_input_A1.txt`.
   b. Đọc vector vế phải $B$ từ file `GS_input_B1.txt`.
   c. Kiểm tra kích thước: $A$ có kích thước $n \times n$, $B$ có kích thước $n$.

2. **Kiểm tra tính chéo trội của ma trận $A$:**
   a. **Chéo trội hàng:** Với mỗi $i$, kiểm tra $|a_{ii}| > \sum_{j \neq i} |a_{ij}|$.
   b. **Chéo trội cột:** Với mỗi $j$, kiểm tra $|a_{jj}| > \sum_{i \neq j} |a_{ij}|$.
   c. Ma trận có thể chéo trội hàng, chéo trội cột, hoặc cả hai.

3. **Xây dựng ma trận đường chéo nghịch đảo $T$:**
   a. Với mỗi $i = 1, 2, \dots, n$: $t_{ii} = 1 / a_{ii}$.
   b. $T$ là ma trận đường chéo với các phần tử ngoài đường chéo bằng $0$.

4. **Xây dựng ma trận lặp $C$ và vector hằng $D$:**
   a. Tính $C = I - T \cdot A$.
      - Với mỗi $i, j$: $c_{ij} = \delta_{ij} - t_{ii} \cdot a_{ij}$.
   b. Tính $D = T \cdot B$.
      - Với mỗi $i$: $d_i = t_{ii} \cdot b_i = b_i / a_{ii}$.

5. **Phân rã ma trận $C$ thành $L$ và $U$:**
   a. **$L$ (tam giác dưới chặt):** $l_{ij} = c_{ij}$ nếu $i > j$, ngược lại $l_{ij} = 0$.
      - Ký hiệu: $L = \text{tril}(C, -1)$ (lấy phần dưới đường chéo chính, không lấy đường chéo).
   b. **$U$ (tam giác trên chặt):** $u_{ij} = c_{ij}$ nếu $i < j$, ngược lại $u_{ij} = 0$.
      - Ký hiệu: $U = \text{triu}(C, 1)$ (lấy phần trên đường chéo chính, không lấy đường chéo).

6. **Tính hệ số hội tụ $q$:**
   a. Với mỗi hàng $i = 1, 2, \dots, n$:
      - Tính tổng trị tuyệt đối các phần tử của $L$ trên hàng $i$: $L_i = \sum_{j=1}^{n} |l_{ij}|$.
      - Tính tổng trị tuyệt đối các phần tử của $U$ trên hàng $i$: $U_i = \sum_{j=1}^{n} |u_{ij}|$.
      - Tính tỷ số: $\displaystyle \text{ratio}_i = \frac{L_i}{1 - U_i}$.
   b. Gán $q = \max_i \text{ratio}_i$.
   c. Đặt $s = 0$ (cho trường hợp chéo trội hàng).

7. **Tính ngưỡng hội tụ:**
   $$
   \text{tol} = \varepsilon \cdot \frac{(1 - q)(1 - s)}{q}
   $$

8. **Khởi tạo vòng lặp:**
   a. Gán vector xấp xỉ ban đầu $x^{(0)} = \mathbf{0}$ (vector không), hoặc giá trị từ file nếu có.
   b. Gán $k = 0$.
   c. Khởi tạo lịch sử lặp: `history = [x^{(0)}]`.

9. **Thực hiện vòng lặp Gauss-Seidel:** Với $k = 0, 1, 2, \dots$ cho đến khi hội tụ hoặc đạt $k_{\max}$.
   a. **Sao chép vector hiện tại:** $x^{(k+1)} = x^{(k)}$ (khởi tạo vector mới).
   b. **Với mỗi thành phần $i = 1, 2, \dots, n$:**
      - Tính tổng đóng góp từ các thành phần đã cập nhật (phần $L$):
        $$
        \text{sumL} = \sum_{j=1}^{i-1} c_{ij} \cdot x_j^{(k+1)}
        $$
      - Tính tổng đóng góp từ các thành phần chưa cập nhật (phần $U$):
        $$
        \text{sumU} = \sum_{j=i+1}^{n} c_{ij} \cdot x_j^{(k)}
        $$
      - Cập nhật thành phần thứ $i$:
        $$
        x_i^{(k+1)} = \text{sumL} + \text{sumU} + d_i
        $$
   c. **Tính sai số tuyệt đối giữa hai bước lặp:**
      $$
      \Delta = \left\| x^{(k+1)} - x^{(k)} \right\|_\infty = \max_{i} \left| x_i^{(k+1)} - x_i^{(k)} \right|
      $$
   d. **Lưu lịch sử:** Thêm $x^{(k+1)}$ vào `history`.
   e. **Kiểm tra điều kiện dừng:**
      - Nếu $\Delta \le \text{tol}$: chuyển sang bước 10.
   f. **Cập nhật:** Gán $k = k + 1$.
      - Nếu $k \ge k_{\max}$: chuyển sang bước 10.

10. **Trả về kết quả:**
    a. Vector nghiệm xấp xỉ cuối cùng $x^{(k)}$.
    b. DataFrame lịch sử lặp (mỗi hàng là một vector $x^{(k)}$).

---

## Thuật toán — Trường hợp ma trận chéo trội cột

**Đầu vào:**
- Ma trận $A$ kích thước $n \times n$ (đọc từ file `GS_input_A1.txt`)
- Vector $B$ kích thước $n$ (đọc từ file `GS_input_B1.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:**
- Vector nghiệm xấp xỉ $x^{(k)}$
- DataFrame lịch sử lặp

**Các bước thực hiện:**

1. **Đọc dữ liệu đầu vào:** (giống bước 1 của trường hợp chéo trội hàng)

2. **Xây dựng ma trận $T$, $C$, $D$ theo biến đổi cột:**
   a. $T = \text{diag}(1/a_{11}, 1/a_{22}, \dots, 1/a_{nn})$.
   b. Tính $C = I - A \cdot T$.
      - Với mỗi $i, j$: $c_{ij} = \delta_{ij} - a_{ij} \cdot t_{jj}$.
   c. Gán $D = B$ (không nhân với $T$).

3. **Phân rã ma trận $C$ thành $L$ và $U$:**
   a. $L = \text{tril}(C, -1)$ (tam giác dưới chặt).
   b. $U = \text{triu}(C, 1)$ (tam giác trên chặt).

4. **Tính hệ số hội tụ $q$ và $s$ (dùng chuẩn cột):**
   a. Với mỗi cột $j = 1, 2, \dots, n$:
      - Tính tổng trị tuyệt đối các phần tử của $L$ trên cột $j$: $L_j = \sum_{i=1}^{n} |l_{ij}|$.
      - Tính tổng trị tuyệt đối các phần tử của $U$ trên cột $j$: $U_j = \sum_{i=1}^{n} |u_{ij}|$.
      - Tính tỷ số: $\displaystyle \text{ratio}_j = \frac{U_j}{1 - L_j}$.
   b. Gán $q = \max_j \text{ratio}_j$.
   c. Gán $s = \max_j L_j$ (giá trị lớn nhất của tổng trị tuyệt đối các cột của $L$).

5. **Tính ngưỡng hội tụ:**
   $$
   \text{tol} = \varepsilon \cdot \frac{(1 - q)(1 - s)}{q}
   $$

6. **Khởi tạo vòng lặp:**
   a. Gán $y^{(0)} = \mathbf{0}$, $x^{(0)} = \mathbf{0}$.
   b. Gán $k = 0$.
   c. Khởi tạo lịch sử lặp: `history = []`.

7. **Thực hiện vòng lặp Gauss-Seidel (dạng $y$):** Với $k = 0, 1, 2, \dots$
   a. **Sao chép vector:** $y^{(k+1)} = y^{(k)}$.
   b. **Với mỗi thành phần $i = 1, 2, \dots, n$:**
      - Tính tổng từ phần $L$ (cập nhật): $\displaystyle \text{sumL} = \sum_{j=1}^{i-1} c_{ij} \cdot y_j^{(k+1)}$.
      - Tính tổng từ phần $U$ (chưa cập nhật): $\displaystyle \text{sumU} = \sum_{j=i+1}^{n} c_{ij} \cdot y_j^{(k)}$.
      - Cập nhật: $y_i^{(k+1)} = \text{sumL} + \text{sumU} + d_i$.
   c. **Tính $x$ tương ứng:** $x^{(k+1)} = T \cdot y^{(k+1)}$.
      - Với mỗi $i$: $x_i^{(k+1)} = y_i^{(k+1)} / a_{ii}$.
   d. **Tính sai số:** $\Delta = \|x^{(k+1)} - x^{(k)}\|_\infty$.
   e. **Lưu lịch sử:** Thêm $x^{(k+1)}$ vào `history`.
   f. **Kiểm tra điều kiện dừng:**
      - Nếu $\Delta \le \text{tol}$: chuyển sang bước 8.
   g. **Cập nhật:** $k = k + 1$.

8. **Trả về kết quả:**
   a. Vector nghiệm xấp xỉ cuối cùng $x^{(k)}$.
   b. DataFrame lịch sử lặp.
