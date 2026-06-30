# Phương Pháp Newton Cho Hệ Phi Tuyến

## Cơ Sở Lý Thuyết

Cho hệ phương trình phi tuyến:

$$
\begin{cases}
f_1(x_1, x_2, \dots, x_n) = 0 \\
f_2(x_1, x_2, \dots, x_n) = 0 \\
\vdots \\
f_n(x_1, x_2, \dots, x_n) = 0
\end{cases}
$$

Khai triển Taylor bậc nhất cho $\mathbf{F}(\mathbf{x})$ tại $\mathbf{x}^{(k)}$:

$$
\mathbf{F}(\mathbf{x}) \approx \mathbf{F}(\mathbf{x}^{(k)}) + \mathbf{J}(\mathbf{x}^{(k)})(\mathbf{x} - \mathbf{x}^{(k)})
$$

với $\mathbf{J}$ là ma trận Jacobian:

$$
\mathbf{J}(\mathbf{x}) =
\begin{bmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \dots & \frac{\partial f_1}{\partial x_n} \\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \dots & \frac{\partial f_2}{\partial x_n} \\
\vdots & \vdots & \ddots & \vdots \\
\frac{\partial f_n}{\partial x_1} & \frac{\partial f_n}{\partial x_2} & \dots & \frac{\partial f_n}{\partial x_n}
\end{bmatrix}
$$

## Công Thức Lặp Newton

Đặt $\mathbf{F}(\mathbf{x}^{(k+1)}) = 0$, ta có:

$$
\mathbf{J}(\mathbf{x}^{(k)}) \, \Delta\mathbf{x}^{(k)} = -\mathbf{F}(\mathbf{x}^{(k)})
$$

với $\Delta\mathbf{x}^{(k)} = \mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}$.

Giải hệ phương trình tuyến tính trên để tìm $\Delta\mathbf{x}^{(k)}$, sau đó cập nhật:

$$
\mathbf{x}^{(k+1)} = \mathbf{x}^{(k)} + \Delta\mathbf{x}^{(k)}
$$

## Điều Kiện Hội Tụ

Nếu $\mathbf{F}$ khả vi liên tục, $\mathbf{J}(\mathbf{x}^*)$ khả nghịch, và $\mathbf{x}^{(0)}$ đủ gần $\mathbf{x}^*$, thì phương pháp Newton hội tụ với tốc độ bậc hai:

$$
\|\mathbf{x}^{(k+1)} - \mathbf{x}^*\| \le C \|\mathbf{x}^{(k)} - \mathbf{x}^*\|^2
$$

## Tiêu Chuẩn Dừng

Dừng khi:

$$
\|\mathbf{x}^{(k)} - \mathbf{x}^{(k-1)}\|_1 = \sum_{i=1}^{n} |x_i^{(k)} - x_i^{(k-1)}| < \text{tol}
$$

với $\text{tol}$ là ngưỡng sai số cho phép.

## Thuật Toán

### newton_method

**Đầu vào:**
- `initial_values`: vector khởi tạo $\mathbf{x}^{(0)} = (x_1^{(0)}, x_2^{(0)}, \dots, x_n^{(0)})^T$
- `tol`: sai số cho phép (mặc định $10^{-6}$), dùng để kiểm tra điều kiện dừng
- `max_iter`: số lần lặp tối đa (mặc định 100), dùng để tránh lặp vô hạn
- `*funcs`: danh sách các hàm $f_1, f_2, \dots, f_n$ và hàm Jacobian

**Đầu ra:** Nghiệm xấp xỉ $\mathbf{x}^*$ của hệ phương trình.

**Các bước:**

1. **Khởi tạo giá trị lặp ban đầu và bảng kết quả:**
   - Gán $\mathbf{x} = \text{initial\_values}$, đây là $\mathbf{x}^{(0)}$
   - Tạo một bảng (danh sách) để lưu kết quả từng bước lặp, bao gồm: số thứ tự lần lặp, vector $\mathbf{x}$ hiện tại, vector $\Delta\mathbf{x}$, chuẩn của $\Delta\mathbf{x}$ và vector giá trị hàm $\mathbf{F}(\mathbf{x})$

2. **Với mỗi lần lặp $i = 0, 1, 2, \dots$ (thực hiện tối đa `max_iter` lần):**
   - **Bước 2.1:** Tính giá trị của các hàm thành phần tại $\mathbf{x}$:
     - Với mỗi $j = 1, 2, \dots, n$, tính $f_j(\mathbf{x})$
     - Tạo vector cột $\mathbf{F}(\mathbf{x}) = \big(f_1(\mathbf{x}), f_2(\mathbf{x}), \dots, f_n(\mathbf{x})\big)^T$
   - **Bước 2.2:** Tính ma trận Jacobian $\mathbf{J}(\mathbf{x})$:
     - Gọi hàm Jacobian (đã được cung cấp trong `*funcs`) để tính ma trận đạo hàm riêng cấp 1 tại $\mathbf{x}$
     - Ma trận $\mathbf{J}(\mathbf{x})$ có kích thước $n \times n$, phần tử tại hàng $i$, cột $j$ là $\dfrac{\partial f_i}{\partial x_j}(\mathbf{x})$
   - **Bước 2.3:** Xây dựng hệ phương trình tuyến tính:
     - Lập hệ $\mathbf{J}(\mathbf{x}) \, \Delta\mathbf{x} = -\mathbf{F}(\mathbf{x})$
     - Trong đó $\Delta\mathbf{x} = (\Delta x_1, \Delta x_2, \dots, \Delta x_n)^T$ là ẩn số cần tìm
   - **Bước 2.4:** Giải hệ phương trình tuyến tính:
     - Sử dụng phương pháp khử Gauss, phân rã LU, hoặc các phương pháp giải hệ tuyến tính khác để tìm $\Delta\mathbf{x}$
     - Kết quả là vector $\Delta\mathbf{x}$ thỏa mãn hệ phương trình trên
   - **Bước 2.5:** Cập nhật nghiệm xấp xỉ:
     - Tính $\mathbf{x}_\text{new} = \mathbf{x} + \Delta\mathbf{x}$
     - Cụ thể, với mỗi thành phần $j = 1, 2, \dots, n$: $x_{\text{new}, j} = x_j + \Delta x_j$
   - **Bước 2.6:** Tính chuẩn 1 của vector chênh lệch:
     - $\text{norm\_diff} = \|\mathbf{x}_\text{new} - \mathbf{x}\|_1 = \sum_{j=1}^{n} |x_{\text{new}, j} - x_j| = \sum_{j=1}^{n} |\Delta x_j|$
   - **Bước 2.7:** Lưu kết quả của lần lặp hiện tại vào bảng (số thứ tự $i$, $\mathbf{x}$, $\Delta\mathbf{x}$, $\text{norm\_diff}$, $\mathbf{F}(\mathbf{x})$)
   - **Bước 2.8:** Cập nhật giá trị lặp: $\mathbf{x} = \mathbf{x}_\text{new}$
   - **Bước 2.9:** Kiểm tra điều kiện dừng:
     - Nếu $\text{norm\_diff} < \text{tol}$ thì **kết thúc vòng lặp** (nghiệm đã hội tụ với độ chính xác yêu cầu)
     - Nếu $i = \text{max\_iter} - 1$ (đã đạt số lần lặp tối đa) thì **kết thúc vòng lặp** (tránh lặp vô hạn)
     - Ngược lại, tăng $i$ lên 1 và quay lại Bước 2.1 để tiếp tục lặp

3. **In bảng kết quả:** Hiển thị bảng kết quả lưu được trong quá trình lặp (gồm tất cả các lần lặp đã thực hiện)

4. **Trả về nghiệm xấp xỉ:** Giá trị $\mathbf{x}$ hiện tại chính là nghiệm xấp xỉ $\mathbf{x}^*$ của hệ phương trình phi tuyến
