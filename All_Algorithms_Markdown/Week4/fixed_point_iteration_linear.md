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

## Thuật toán

**Đầu vào:**
- Ma trận hệ số $A$ kích thước $n \times n$ (đọc từ file `FXP_input_A.txt`)
- Vector vế phải $B$ kích thước $n$ (đọc từ file `FXP_input_B.txt`)
- Vector xấp xỉ ban đầu $x^{(0)}$ kích thước $n$ (đọc từ file `FXP_input_X0.txt`)
- Ngưỡng sai số $\varepsilon > 0$
- Số bước lặp tối đa $k_{\max}$

**Đầu ra:**
- Vector nghiệm xấp xỉ $x^{(k)}$
- Lịch sử các bước lặp (danh sách các vector $x^{(0)}, x^{(1)}, \dots, x^{(k)}$)

**Các bước thực hiện:**

1. **Đọc dữ liệu đầu vào từ file:**
   a. Đọc ma trận vuông $A$ từ file `FXP_input_A.txt`.
   b. Đọc vector cột $B$ từ file `FXP_input_B.txt`.
   c. Đọc vector xấp xỉ ban đầu $x^{(0)}$ từ file `FXP_input_X0.txt`.
   d. Kiểm tra kích thước: số hàng của $A$ phải bằng số phần tử của $B$ và của $x^{(0)}$.

2. **Tính các chuẩn của ma trận $A$ để đánh giá điều kiện hội tụ:**
   a. **Chuẩn cột (chuẩn 1):** $\displaystyle \|A\|_1 = \max_j \sum_{i=1}^{n} |a_{ij}|$
   b. **Chuẩn hàng (chuẩn $\infty$):** $\displaystyle \|A\|_\infty = \max_i \sum_{j=1}^{n} |a_{ij}|$
   c. **Chuẩn 2:** $\|A\|_2 = \sqrt{\lambda_{\max}(A^T A)}$, với $\lambda_{\max}$ là trị riêng lớn nhất của $A^T A$.
   d. **Chuẩn max:** $\|A\|_{\max} = 3 \cdot \max_{i,j} |a_{ij}|$

3. **Xác định hệ số co $q$:**
   a. Gán $q = \|A\|_1$ (chuẩn cột).
   b. Kiểm tra điều kiện hội tụ: nếu $q \ge 1$, phương pháp có thể không hội tụ; vẫn tiếp tục nhưng cần cảnh báo.

4. **Tính ngưỡng hội tụ hiệu chỉnh:**
   $$
   \varepsilon_{\text{new}} = \varepsilon \cdot \frac{1 - q}{q}
   $$
   (Ngưỡng này dùng để dừng vòng lặp sớm hơn dựa trên sai số tiên nghiệm.)

5. **Khởi tạo vòng lặp:**
   a. Gán $k = 0$.
   b. Khởi tạo danh sách lịch sử lặp: `history = [x^{(0)}]`.

6. **Thực hiện vòng lặp lặp điểm cố định:** Lặp với $k = 0, 1, 2, \dots$ cho đến khi hội tụ hoặc đạt $k_{\max}$.
   a. **Tính xấp xỉ mới:** $x^{(k+1)} = A \cdot x^{(k)} + B$.
      - Với mỗi thành phần $i$: $\displaystyle x_i^{(k+1)} = \sum_{j=1}^{n} a_{ij} x_j^{(k)} + b_i$.
   b. **Tính sai số giữa hai bước lặp liên tiếp:**
      $$
      \Delta = \left\| x^{(k+1)} - x^{(k)} \right\|
      $$
      (Sử dụng chuẩn L1 hoặc chuẩn max.)
   c. **Lưu lịch sử:** Thêm $x^{(k+1)}$ vào danh sách `history`.
   d. **Kiểm tra điều kiện dừng:**
      - Nếu $\Delta < \varepsilon_{\text{new}}$: dừng vòng lặp, chuyển sang bước 7.
   e. **Cập nhật:** Gán $k = k + 1$.
      - Nếu $k \ge k_{\max}$: dừng vòng lặp, chuyển sang bước 7.

7. **Trả về kết quả:**
   a. Vector nghiệm xấp xỉ cuối cùng $x^{(k)}$.
   b. Lịch sử lặp `history` (danh sách tất cả các vector $x^{(0)}, x^{(1)}, \dots, x^{(k)}$).
