# Phương pháp Danilevsky (Danilevsky Method)

## Công thức toán học

Đưa ma trận $A$ về dạng Frobenius (dạng companion) bằng phép biến đổi đồng dạng:

$$P = S^{-1} A S$$

Ma trận Frobenius có dạng:

$$F = \begin{pmatrix}
0 & 0 & \cdots & 0 & c_0 \\
1 & 0 & \cdots & 0 & c_1 \\
0 & 1 & \cdots & 0 & c_2 \\
\vdots & \vdots & \ddots & \vdots & \vdots \\
0 & 0 & \cdots & 1 & c_{n-1}
\end{pmatrix}$$

Đa thức đặc trưng:

$$P(\lambda) = \lambda^n - c_{n-1}\lambda^{n-1} - \cdots - c_1\lambda - c_0 = 0$$

---

## Thuật toán

**Mục tiêu:** Tìm trị riêng và vector riêng của ma trận $A$ bằng phương pháp Danilevsky.
**Đầu vào:** Ma trận vuông $A$ kích thước $n \times n$, ngưỡng sai số $tol$.
**Đầu ra:** Trị riêng $\lambda_i$ và vector riêng $v_i$.

### Phần A: Biến đổi về dạng Frobenius

**Bước 1:** Khởi tạo
   - Sao chép $P = A$.
   - Khởi tạo danh sách $S_{list} = []$.

**Bước 2:** Với $k = n-1$ xuống $1$:

   **Bước 2.1:** Kiểm tra phần tử $P[k, k-1]$

   - **Trường hợp 1:** $|P[k, k-1]| \ge tol$ (phần tử dưới đường chéo khác 0)
     
     **Bước 2.1.1:** Tạo ma trận biến đổi $S$:
     - $S[k-1] = -\dfrac{P[k]}{P[k,k-1]}$.
     - $S[k-1,k-1] = \dfrac{1}{P[k,k-1]}$.
     
     **Bước 2.1.2:** Thêm $S$ vào $S_{list}$.
     
     **Bước 2.1.3:** Cập nhật $P = S^{-1} \cdot P \cdot S$.

   - **Trường hợp 2:** $|P[k, k-1]| < tol$ và tồn tại $j < k-1$ sao cho $|P[k,j]| \ge tol$

     **Bước 2.2.1:** Tạo ma trận hoán vị $C$.
     
     **Bước 2.2.2:** Thực hiện $P = C \cdot P \cdot C$.
     
     **Bước 2.2.3:** Áp dụng Trường hợp 1 cho ma trận đã hoán vị.

   - **Trường hợp 3:** $|P[k,j]| < tol$ với mọi $j \le k-1$
     * Không cần biến đổi, tiếp tục.

### Phần B: Trích hệ số đa thức đặc trưng

**Bước 3:** Trích hệ số đa thức đặc trưng
   - $coeffs = [1, -P[0,0], -P[0,1], \dots, -P[0,n-1]]$.

**Bước 4:** Tìm trị riêng
   - $\lambda = roots(coeffs)$ (giải phương trình đặc trưng).

### Phần C: Tìm vector riêng

**Bước 5:** Với mỗi $\lambda_i$, tạo vector riêng của ma trận Frobenius
   - $v = [1, \lambda_i, \lambda_i^2, \dots, \lambda_i^{n-1}]^T$.

**Bước 6:** Biến đổi ngược vector riêng
   - $v = S_1 \cdot S_2 \cdot \dots \cdot S_k \cdot v$.

**Bước 7:** Chuẩn hóa vector riêng.

**Bước 8:** Trả về các cặp $(\lambda_i, v_i)$.
