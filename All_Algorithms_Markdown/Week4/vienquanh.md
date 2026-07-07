# Phương pháp nghịch đảo ma trận khối (Vienquanh / Banachiewicz)

## Công thức toán học

### 1. Ý tưởng chính

Phương pháp Vienquanh (Block Matrix Recursion) tính ma trận nghịch đảo bằng cách xây dựng đệ quy nghịch đảo của các ma trận con chính (leading principal minors) từ kích thước $1 \times 1$ lên $n \times n$.

### 2. Công thức đệ quy

Cho $A_k$ là ma trận con chính kích thước $k \times k$ của $A$. Giả sử đã biết $B_{k-1} = A_{k-1}^{-1}$, ta xây dựng $B_k = A_k^{-1}$ như sau:

**Đại lượng $\theta_k$:**

$$\theta_k = r^T \cdot B_{k-1} \cdot c - a_{kk}$$

Trong đó:
- $c = A_{1:k-1, k}$ (cột cuối của $A_k$ trừ $a_{kk}$)
- $r^T = A_{k, 1:k-1}$ (hàng cuối của $A_k$ trừ $a_{kk}$)
- $a_{kk} = A_{k,k}$

**Phần tử góc dưới phải:** $$b_{k,k} = -\frac{1}{\theta_k}$$

**Cột biên phải:** $$\beta_{1:k-1, k} = \frac{B_{k-1} \cdot c}{\theta_k}$$

**Hàng biên dưới:** $$\beta_{k, 1:k-1} = \frac{r^T \cdot B_{k-1}}{\theta_k}$$

**Khối góc trên trái cập nhật:** $$B_{k-1}^{(new)} = B_{k-1} \cdot (I - c \otimes \beta_{k, 1:k-1})$$

### 3. Nghịch đảo gián tiếp qua $(A^T A)$

Khi $A$ không vuông, tính $M = A^T A$, nghịch đảo $M$ bằng đệ quy khối, rồi:

$$A^{-1} = M^{-1} \cdot A^T$$

---

## Thuật toán

**Mục tiêu:** Tính ma trận nghịch đảo $A^{-1}$ bằng phương pháp đệ quy khối.
**Đầu vào:** Ma trận $A_{m \times n}$ đọc từ file `BLMT_input_A.txt`.
**Đầu ra:** Ma trận nghịch đảo $A^{-1}$ (nếu $A$ vuông) hoặc ma trận giả nghịch đảo (nếu $A$ không vuông).

### Phần A: Xử lý ma trận đầu vào

**Bước 1:** Đọc dữ liệu đầu vào
   - **Bước 1.1:** Mở file `BLMT_input_A.txt`, đọc ma trận $A$.
   - **Bước 1.2:** Xác định số hàng $m$ và số cột $n$ của $A$.
   - **Bước 1.3:** Kiểm tra $m > 0$ và $n > 0$.

**Bước 2:** Phân nhánh theo dạng ma trận
   - **Bước 2.1:** Nếu $m = n$ (ma trận vuông): chuyển đến Phần B.
   - **Bước 2.2:** Nếu $m \neq n$ (ma trận không vuông): chuyển đến Phần C.

### Phần B: Nghịch đảo trực tiếp ma trận vuông (đệ quy khối)

**Bước 3:** Khởi tạo ma trận con đầu tiên $B_1 = A_1^{-1}$
   - **Bước 3.1:** Lấy $A_1 = [A[0][0]]$.
   - **Bước 3.2:** Kiểm tra $|A[0][0]| < 10^{-10}$: nếu đúng, báo lỗi và kết thúc.
   - **Bước 3.3:** Tính $B_1 = [1 / A[0][0]]$.
   - **Bước 3.4:** Gán $B_{k-1} = B_1$.

**Bước 4:** Với mỗi $k$ từ $2$ đến $n$ (mở rộng lên kích thước $k \times k$)

   **Bước 4.1:** Trích xuất các thành phần của $A_k$
   - Lấy cột $c$: $c[i] = A[i][k-1]$ với $i = 0, \dots, k-2$.
   - Lấy hàng $r^T$: $r^T[j] = A[k-1][j]$ với $j = 0, \dots, k-2$.
   - Lấy $a_{kk} = A[k-1][k-1]$.

   **Bước 4.2:** Tính $\theta_k$
   - $temp1 = B_{k-1} \times c$: $temp1[i] = \sum_{j=0}^{k-2} B_{k-1}[i][j] \times c[j]$.
   - $theta = r^T \times temp1$: $theta = \sum_{j=0}^{k-2} r^T[j] \times temp1[j]$.
   - $\theta_k = theta - a_{kk}$.
   - Kiểm tra $|\theta_k| < 10^{-10}$: nếu đúng, báo lỗi và kết thúc.

   **Bước 4.3:** Tính cột biên phải $\beta_{col}$
   - Với mỗi $i$ từ $0$ đến $k-2$: $\beta_{col}[i] = temp1[i] / \theta_k$.

   **Bước 4.4:** Tính hàng biên dưới $\beta_{row}$
   - $temp2 = r^T \times B_{k-1}$: $temp2[j] = \sum_{i=0}^{k-2} r^T[i] \times B_{k-1}[i][j]$.
   - Với mỗi $j$ từ $0$ đến $k-2$: $\beta_{row}[j] = temp2[j] / \theta_k$.

   **Bước 4.5:** Tính phần tử góc dưới phải $b_{kk} = -1 / \theta_k$.

   **Bước 4.6:** Cập nhật khối góc trên trái $B_{tl}$
   - $correction[i][j] = c[i] \times \beta_{row}[j]$.
   - $I\_minus\_correction = I - correction$.
   - $B_{tl} = B_{k-1} \times I\_minus\_correction$.

   **Bước 4.7:** Lắp ráp ma trận $B_k$
   - Tạo $B_k$ kích thước $k \times k$.
   - Gán khối trên trái: $B_k[0:k-1][0:k-1] = B_{tl}$.
   - Gán cột biên phải: $B_k[0:k-1][k-1] = \beta_{col}$.
   - Gán hàng biên dưới: $B_k[k-1][0:k-1] = \beta_{row}$.
   - Gán góc dưới phải: $B_k[k-1][k-1] = b_{kk}$.

   **Bước 4.8:** Cập nhật $B_{k-1} = B_k$ cho vòng lặp tiếp theo.

**Bước 5:** Kết thúc đệ quy khối
   - Sau khi $k = n$, $B_n$ chính là $A^{-1}$.
   - Trả về $B_n$.

### Phần C: Nghịch đảo gián tiếp cho ma trận không vuông

**Bước 6:** Tính ma trận $M = A^T \times A$
   - **Bước 6.1:** Tính $A^T$ (kích thước $n \times m$).
   - **Bước 6.2:** $M[i][j] = \sum_{t=0}^{m-1} A^T[i][t] \times A[t][j]$.

**Bước 7:** Tính nghịch đảo của $M$ bằng đệ quy khối
   - Thực hiện các bước từ 3 đến 5 với ma trận $M$.
   - Kết quả là $M^{-1}$ (kích thước $n \times n$).

**Bước 8:** Tính ma trận giả nghịch đảo $A^{-1}$
   - $A^{-1}[i][j] = \sum_{t=0}^{n-1} M^{-1}[i][t] \times A^T[t][j]$.
   - Trả về $A^{-1}$.

### Phần D: Hoàn tất

**Bước 9:** Xuất kết quả
   - In ma trận nghịch đảo (hoặc giả nghịch đảo) ra màn hình.
   - Ghi kết quả ra file nếu yêu cầu.
