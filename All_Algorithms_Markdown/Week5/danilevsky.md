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

## Thuật toán

**Đầu vào:** Ma trận vuông $A$ kích thước $n \times n$, ngưỡng sai số $tol$

**Đầu ra:** Trị riêng $\lambda_i$ và vector riêng $v_i$

1. Sao chép $P = A$, khởi tạo danh sách $S_{list} = []$
2. Với $k = n-1$ xuống $1$:
   - **TH1:** $|P[k, k-1]| \ge tol$ (phần tử dưới đường chéo khác 0):
     - Tạo ma trận biến đổi $S$:
       $$S[k-1] = -\dfrac{P[k]}{P[k,k-1]},\quad S[k-1,k-1] = \dfrac{1}{P[k,k-1]}$$
     - $S_{list}.append(S)$
     - $P = S^{-1} \cdot P \cdot S$
   - **TH2:** $|P[k, k-1]| < tol$ và tồn tại $j < k-1$ sao cho $|P[k,j]| \ge tol$:
     - Tạo ma trận hoán vị $C$, thực hiện $P = C \cdot P \cdot C$
     - Áp dụng TH1 cho ma trận đã hoán vị
   - **TH3:** $|P[k,j]| < tol$ với mọi $j \le k-1$ — không cần biến đổi
3. Trích hệ số đa thức đặc trưng: $coeffs = [1, -P[0,0], -P[0,1], \dots, -P[0,n-1]]$
4. Tìm trị riêng: $\lambda = roots(coeffs)$
5. Với mỗi $\lambda_i$, tạo vector riêng của ma trận Frobenius:
   $$v = [1, \lambda_i, \lambda_i^2, \dots, \lambda_i^{n-1}]^T$$
6. Biến đổi ngược: $v = S_1 \cdot S_2 \cdot \dots \cdot S_k \cdot v$
7. Chuẩn hóa vector riêng
