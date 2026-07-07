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

---

## Thuật toán

**Mục tiêu:** Tính phân rã SVD của ma trận $A$.
**Đầu vào:** Ma trận $A$ kích thước $m \times n$, ngưỡng sai số $tol$.
**Đầu ra:** $U$, $\Sigma$, $V^T$.

**Bước 1:** Chọn hướng tính
   - Nếu $m \le n$: dùng $M = A \cdot A^T$.
   - Nếu $m > n$: dùng $M = A^T \cdot A$.

**Bước 2:** Tính $M = A \cdot A^T$ (hoặc $A^T \cdot A$).

**Bước 3:** Tìm trị riêng và vector riêng của $M$ bằng `np.linalg.eigh`.

**Bước 4:** Sắp xếp trị riêng giảm dần.

**Bước 5:** Tính giá trị kỳ dị: $\sigma_i = \sqrt{\lambda_i}$.

**Bước 6:** Xây dựng ma trận $\Sigma$ với $\sigma_i$ trên đường chéo.

**Bước 7:** Tính vector suy biến phải $v_i = A^T \cdot u_i / \sigma_i$.

**Bước 8:** Hoàn thiện $V$ thành cơ sở trực chuẩn đầy đủ (dùng `complete_orthonormal_basis`).

**Bước 9:** $U$ là ma trận các vector riêng của $M$, $V^T$ là chuyển vị của $V$.

**Bước 10:** Trả về $U$, $\Sigma$, $V^T$.

---

## Hàm bổ sung: `complete_orthonormal_basis`

**Mục tiêu:** Mở rộng $B$ ($n \times k$) có các cột trực chuẩn thành cơ sở trực chuẩn $n \times n$ đầy đủ.
**Đầu vào:** Ma trận $B$ kích thước $n \times k$ với các cột trực chuẩn.
**Đầu ra:** Ma trận $Q$ kích thước $n \times n$ trực chuẩn.

**Bước 1:** Khởi tạo $Q$ kích thước $n \times n$, gán $Q[:, :k] = B$.
**Bước 2:** Với $j = k$ đến $n-1$:
   - **Bước 2.1:** Tạo vector ngẫu nhiên $v$.
   - **Bước 2.2:** Trực giao hóa $v$ với các cột $Q[:, :j]$.
   - **Bước 2.3:** Chuẩn hóa $v$ và gán vào $Q[:, j]$.
**Bước 3:** Trả về $Q$.
