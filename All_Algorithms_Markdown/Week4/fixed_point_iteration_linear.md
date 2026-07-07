# Phương pháp lặp điểm cố định cho hệ phương trình tuyến tính

## Công thức toán học

Xét hệ phương trình tuyến tính:

$$
A x = B
$$

Đưa về dạng lặp điểm cố định:

$$
x = (I - A)x + B
$$

Công thức lặp:

$$
x^{(k+1)} = A x^{(k)} + B
$$

Điều kiện hội tụ: $\|A\| < 1$

Sai số tiên nghiệm:

$$
\|x^{(k)} - x^*\| \le \frac{q}{1-q} \|x^{(k)} - x^{(k-1)}\|
$$

với $q = \|A\|$ là hệ số co.

Ngưỡng hội tụ:

$$
\varepsilon_{\text{new}} = \varepsilon \cdot \frac{1-q}{q}
$$

---

## Thuật toán

**Mục tiêu:** Giải hệ phương trình tuyến tính $Ax = B$ bằng phương pháp lặp điểm cố định.
**Đầu vào:**
- Ma trận hệ số $A$ kích thước $n \times n$ (đọc từ file `FXP_input_A.txt`)
- Vector vế phải $B$ kích thước $n$ (đọc từ file `FXP_input_B.txt`)
- Vector xấp xỉ ban đầu $x^{(0)}$ kích thước $n$ (đọc từ file `FXP_input_X0.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:**
- Vector nghiệm xấp xỉ $x^{(k)}$
- Lịch sử các bước lặp

**Bước 1:** Đọc dữ liệu đầu vào từ file
   - **Bước 1.1:** Đọc ma trận vuông $A$ từ file `FXP_input_A.txt`.
   - **Bước 1.2:** Đọc vector cột $B$ từ file `FXP_input_B.txt`.
   - **Bước 1.3:** Đọc vector xấp xỉ ban đầu $x^{(0)}$ từ file `FXP_input_X0.txt`.
   - **Bước 1.4:** Kiểm tra kích thước: số hàng của $A$ phải bằng số phần tử của $B$ và của $x^{(0)}$.

**Bước 2:** Tính các chuẩn của ma trận $A$ để đánh giá điều kiện hội tụ
   - **Bước 2.1:** Chuẩn cột (chuẩn 1): $\|A\|_1 = \max_j \sum_{i=1}^{n} |a_{ij}|$.
   - **Bước 2.2:** Chuẩn hàng (chuẩn $\infty$): $\|A\|_\infty = \max_i \sum_{j=1}^{n} |a_{ij}|$.
   - **Bước 2.3:** Chuẩn 2: $\|A\|_2 = \sqrt{\lambda_{\max}(A^T A)}$.
   - **Bước 2.4:** Chuẩn max: $\|A\|_{\max} = 3 \cdot \max_{i,j} |a_{ij}|$.

**Bước 3:** Xác định hệ số co $q$
   - **Bước 3.1:** Gán $q = \|A\|_1$ (chuẩn cột).
   - **Bước 3.2:** Nếu $q \ge 1$, phương pháp có thể không hội tụ; vẫn tiếp tục nhưng cảnh báo.

**Bước 4:** Tính ngưỡng hội tụ hiệu chỉnh
   $$\varepsilon_{\text{new}} = \varepsilon \cdot \frac{1 - q}{q}$$

**Bước 5:** Khởi tạo vòng lặp
   - **Bước 5.1:** Gán $k = 0$.
   - **Bước 5.2:** Khởi tạo lịch sử lặp: $\text{history} = [x^{(0)}]$.

**Bước 6:** Thực hiện vòng lặp lặp điểm cố định

   **Bước 6.1:** Tính xấp xỉ mới: $x^{(k+1)} = A \cdot x^{(k)} + B$.
   - Với mỗi $i$: $\displaystyle x_i^{(k+1)} = \sum_{j=1}^{n} a_{ij} x_j^{(k)} + b_i$.

   **Bước 6.2:** Tính sai số giữa hai bước lặp:
   $$\Delta = \left\| x^{(k+1)} - x^{(k)} \right\|$$
   (Sử dụng chuẩn L1 hoặc chuẩn max.)

   **Bước 6.3:** Lưu lịch sử: Thêm $x^{(k+1)}$ vào $\text{history}$.

   **Bước 6.4:** Kiểm tra điều kiện dừng
   - Nếu $\Delta < \varepsilon_{\text{new}}$: chuyển sang Bước 7.

   **Bước 6.5:** Cập nhật $k = k + 1$.
   - Nếu $k \ge k_{\max}$: chuyển sang Bước 7.
   - Ngược lại, quay lại Bước 6.1.

**Bước 7:** Trả về kết quả
   - **Bước 7.1:** Vector nghiệm xấp xỉ cuối cùng $x^{(k)}$.
   - **Bước 7.2:** Lịch sử lặp $\text{history}$.
