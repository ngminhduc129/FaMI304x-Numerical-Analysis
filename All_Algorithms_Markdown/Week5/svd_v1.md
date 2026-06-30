# Phân rã Giá trị Kỳ dị (Singular Value Decomposition) — v1

## Công thức toán học

Phân rã SVD của ma trận $A$ kích thước $m \times n$:

$$A = U \cdot \Sigma \cdot V^T$$

Trong đó:
- $U$: ma trận trực giao $m \times m$ (các vector suy biến trái)
- $\Sigma$: ma trận đường chéo $m \times n$ chứa các giá trị kỳ dị $\sigma_i$
- $V$: ma trận trực giao $n \times n$ (các vector suy biến phải)

Tính toán qua ma trận $M = A \cdot A^T$ (khi $m \le n$):

$$M = A \cdot A^T \quad (m \times m)$$

Giải bài toán trị riêng của $M$:

$$M \cdot u_i = \lambda_i \cdot u_i$$

Giá trị kỳ dị: $\sigma_i = \sqrt{\lambda_i}$

Vector suy biến phải:

$$v_i = \frac{A^T \cdot u_i}{\sigma_i}$$

## Thuật toán

**Đầu vào:** Ma trận $A$ kích thước $m \times n$, ngưỡng sai số $tol$

**Đầu ra:** $U$, $\Sigma$, $V^T$

1. Chọn hướng tính: nếu $m \le n$, dùng $M = A \cdot A^T$
2. Tính $M = A \cdot A^T$
3. Tìm trị riêng và vector riêng của $M$ bằng `np.linalg.eigh`
4. Sắp xếp trị riêng giảm dần
5. Tính giá trị kỳ dị: $\sigma_i = \sqrt{\lambda_i}$
6. Xây dựng ma trận $\Sigma$ với $\sigma_i$ trên đường chéo
7. Tính vector suy biến phải $v_i = A^T \cdot u_i / \sigma_i$
8. Hoàn thiện $V$ thành cơ sở trực chuẩn đầy đủ
9. $U$ là ma trận các vector riêng của $M$, $V^T$ là chuyển vị của $V$

## Hàm bổ sung: complete_orthonormal_basis

Mở rộng $B$ ($n \times k$) có các cột trực chuẩn thành cơ sở trực chuẩn $n \times n$ đầy đủ:

1. Khởi tạo $Q$ kích thước $n \times n$, gán $Q[:, :k] = B$
2. Với $j = k$ đến $n-1$:
   - Tạo vector ngẫu nhiên $v$
   - Trực giao hóa $v$ với các cột $Q[:, :j]$
   - Chuẩn hóa $v$ và gán vào $Q[:, j]$
