# Phân rã Giá trị Kỳ dị dùng Lũy thừa & Xuống thang (SVD with Power Method & Deflation) — v2

## Công thức toán học

Phân rã SVD của ma trận $A$ kích thước $m \times n$:

$$A = U \cdot \Sigma \cdot V^T$$

Khác với v1, v2 **không dùng `np.linalg.eigh`** mà dùng phương pháp lũy thừa + xuống thang để tìm trị riêng.

Tính toán qua ma trận $M = A \cdot A^T$ (khi $m \le n$):

$$M = A \cdot A^T \quad (m \times m)$$

Dùng `compute_all_eigenpairs` để tìm tất cả trị riêng và vector riêng của $M$:

$$M \cdot u_i = \lambda_i \cdot u_i$$

Giá trị kỳ dị: $\sigma_i = \sqrt{\lambda_i}$

Vector suy biến phải:

$$v_i = \frac{A^T \cdot u_i}{\sigma_i}$$

### Xuống thang (Deflation)

Sau khi tìm được một cặp trị riêng-vector riêng $(\lambda, v)$:

$$A_{new} = A - \lambda \cdot v \cdot w^T, \quad \text{với } A^T w = \lambda w,\; w^T v = 1$$

---

## Thuật toán `compute_svd` (qua $A \cdot A^T$, cho $m \le n$)

**Mục tiêu:** Tính SVD của $A$ dùng Power Method + Deflation.
**Đầu vào:** Ma trận $A$ kích thước $m \times n$ ($m \le n$).
**Đầu ra:** $U$, $\Sigma$, $V^T$.

**Bước 1:** Tính $M = A \cdot A^T$.
**Bước 2:** Gọi `compute_all_eigenpairs(M)` để tìm tất cả $(\lambda_i, u_i)$ bằng lũy thừa + xuống thang.
**Bước 3:** Sắp xếp $\lambda_i$ giảm dần.
**Bước 4:** Tính giá trị kỳ dị: $\sigma_i = \sqrt{\max(\lambda_i, 0)}$.
**Bước 5:** Xây dựng ma trận $\Sigma$ với $\sigma_i$ trên đường chéo.
**Bước 6:** Tính $v_i = A^T \cdot u_i / \sigma_i$ (với $\sigma_i > 0$).
**Bước 7:** Hoàn thiện $V$ thành cơ sở trực chuẩn đầy đủ.
**Bước 8:** $U$ là ma trận các vector riêng của $M$.
**Bước 9:** Trả về $U$, $\Sigma$, $V^T$.

---

## Thuật toán `compute_svd_2` (qua $A^T \cdot A$, cho $m > n$)

**Mục tiêu:** Tính SVD của $A$ khi $m > n$.
**Đầu vào:** Ma trận $A$ kích thước $m \times n$ ($m > n$).
**Đầu ra:** $U$, $\Sigma$, $V^T$.

**Bước 1:** Tính $N = A^T \cdot A$.
**Bước 2:** Gọi `compute_all_eigenpairs(N)` để tìm tất cả $(\lambda_i, v_i)$.
**Bước 3:** Sắp xếp $\lambda_i$ giảm dần.
**Bước 4:** Tính $\sigma_i = \sqrt{\max(\lambda_i, 0)}$.
**Bước 5:** Xây dựng $\Sigma$.
**Bước 6:** Tính $u_i = A \cdot v_i / \sigma_i$.
**Bước 7:** Hoàn thiện $U$ thành cơ sở trực chuẩn đầy đủ.
**Bước 8:** $V$ là ma trận các vector riêng của $N$, $V^T$ là chuyển vị.
**Bước 9:** Trả về $U$, $\Sigma$, $V^T$.
