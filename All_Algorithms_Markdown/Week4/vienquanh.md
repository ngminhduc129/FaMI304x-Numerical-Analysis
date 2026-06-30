# Phương pháp nghịch đảo ma trận khối (Vienquanh / Banachiewicz)

## Công thức toán học

### 1. Ý tưởng chính

Phương pháp Vienquanh (Block Matrix Recursion) tính ma trận nghịch đảo bằng cách xây dựng đệ quy nghịch đảo của các ma trận con chính (leading principal minors) từ kích thước $1 \times 1$ lên $n \times n$.

### 2. Công thức đệ quy

Cho $A_k$ là ma trận con chính kích thước $k \times k$ của $A$. Giả sử đã biết $B_{k-1} = A_{k-1}^{-1}$, ta xây dựng $B_k = A_k^{-1}$ như sau:

**Đại lượng $\theta_k$:**

$$\theta_k = r^T \cdot B_{k-1} \cdot c - a_{kk}$$

Trong đó:
- $c = A_{1:k-1, k}$ (cột cuối của $A_k$ trừ phần tử $a_{kk}$)
- $r^T = A_{k, 1:k-1}$ (hàng cuối của $A_k$ trừ phần tử $a_{kk}$)
- $a_{kk} = A_{k,k}$

**Phần tử góc dưới phải:**

$$b_{k,k} = -\frac{1}{\theta_k}$$

**Cột biên phải** ($\beta_{1:k-1, k}$):

$$\beta_{1:k-1, k} = \frac{B_{k-1} \cdot c}{\theta_k}$$

**Hàng biên dưới** ($\beta_{k, 1:k-1}$):

$$\beta_{k, 1:k-1} = \frac{r^T \cdot B_{k-1}}{\theta_k}$$

**Khối góc trên trái cập nhật:**

$$B_{k-1}^{(new)} = B_{k-1} \cdot (I - c \otimes \beta_{k, 1:k-1})$$

**Ma trận $B_k$ hoàn chỉnh:**

$$B_k = \begin{bmatrix}
B_{k-1}^{(new)} & \beta_{1:k-1, k} \\
\beta_{k, 1:k-1} & b_{k,k}
\end{bmatrix}$$

### 3. Nghịch đảo gián tiếp qua $(A^T A)$

Khi $A$ không vuông, tính $M = A^T A$, nghịch đảo $M$ bằng đệ quy khối, rồi:

$$A^{-1} = M^{-1} \cdot A^T$$

## Thuật toán

**Đầu vào:** Ma trận $A_{m \times n}$ đọc từ file `BLMT_input_A.txt`.  
**Đầu ra:** Ma trận nghịch đảo $A^{-1}$ (nếu $A$ vuông) hoặc ma trận giả nghịch đảo (nếu $A$ không vuông).

### Phần A: Xử lý ma trận đầu vào

1. **Đọc dữ liệu đầu vào:**
   a. Mở file `BLMT_input_A.txt`, đọc ma trận $A$.
   b. Xác định số hàng $m$ và số cột $n$ của $A$.
   c. Kiểm tra tính hợp lệ: $m > 0$ và $n > 0$.

2. **Phân nhánh theo dạng ma trận:**
   a. **Nếu $m = n$ (ma trận vuông):**
      - Chuyển đến Phần B: tính nghịch đảo trực tiếp bằng đệ quy khối trên $A$.
   b. **Nếu $m \neq n$ (ma trận không vuông):**
      - Chuyển đến Phần C: tính nghịch đảo gián tiếp qua $A^T A$.

### Phần B: Nghịch đảo trực tiếp ma trận vuông (đệ quy khối)

3. **Khởi tạo ma trận con đầu tiên $B_1 = A_1^{-1}$:**
   a. Lấy phần tử $A[0][0]$ của ma trận gốc làm ma trận $1 \times 1$: $A_1 = [A[0][0]]$.
   b. Kiểm tra $|A[0][0]| < 10^{-10}$:
      - Nếu đúng, thông báo ma trận suy biến và kết thúc.
   c. Tính $B_1 = [1 / A[0][0]]$ (ma trận $1 \times 1$ chứa nghịch đảo của $A[0][0]$).
   d. Gán $B_{k-1} = B_1$ (khởi tạo cho vòng lặp).

4. **Với mỗi $k$ từ $2$ đến $n$ (mở rộng dần ma trận con lên kích thước $k \times k$):**

   a. **Trích xuất các thành phần của $A_k$ từ ma trận gốc $A$:**
      - Lấy cột $c$: vector cột gồm $k-1$ phần tử từ $A[0:k-1][k-1]$ (cột thứ $k-1$, các hàng $0$ đến $k-2$).
        $$c[i] = A[i][k-1] \quad \text{với } i = 0, 1, \dots, k-2$$
      - Lấy hàng $r^T$: vector hàng gồm $k-1$ phần tử từ $A[k-1][0:k-1]$ (hàng thứ $k-1$, các cột $0$ đến $k-2$).
        $$r^T[j] = A[k-1][j] \quad \text{với } j = 0, 1, \dots, k-2$$
      - Lấy phần tử góc $a_{kk} = A[k-1][k-1]$.

   b. **Tính $\theta_k$:**
      - Tính vector trung gian $temp1 = B_{k-1} \times c$ (nhân ma trận $B_{k-1}$ với vector cột $c$).
        * Với mỗi $i$ từ $0$ đến $k-2$:
          $$temp1[i] = \sum_{j=0}^{k-2} B_{k-1}[i][j] \times c[j]$$
      - Tính vô hướng $theta = r^T \times temp1$:
        $$theta = \sum_{j=0}^{k-2} r^T[j] \times temp1[j]$$
      - Tính $\theta_k = theta - a_{kk}$.
      - Kiểm tra $|\theta_k| < 10^{-10}$:
        * Nếu đúng, ma trận suy biến, thông báo lỗi và kết thúc.

   c. **Tính cột biên phải $\beta_{col}$:**
      - $\beta_{col}$ là vector cột có $k-1$ phần tử.
      - Với mỗi $i$ từ $0$ đến $k-2$:
        $$\beta_{col}[i] = temp1[i] / \theta_k$$

   d. **Tính hàng biên dưới $\beta_{row}$:**
      - $\beta_{row}$ là vector hàng có $k-1$ phần tử.
      - Tính vector trung gian $temp2 = r^T \times B_{k-1}$:
        * Với mỗi $j$ từ $0$ đến $k-2$:
          $$temp2[j] = \sum_{i=0}^{k-2} r^T[i] \times B_{k-1}[i][j]$$
      - Với mỗi $j$ từ $0$ đến $k-2$:
        $$\beta_{row}[j] = temp2[j] / \theta_k$$

   e. **Tính phần tử góc dưới phải $b_{kk}$:**
      $$b_{kk} = -1 / \theta_k$$

   f. **Cập nhật khối góc trên trái $B_{tl}$:**
      - Tính ma trận hiệu chỉnh: với mỗi $i$ từ $0$ đến $k-2$, mỗi $j$ từ $0$ đến $k-2$:
        $$correction[i][j] = c[i] \times \beta_{row}[j]$$
      - Tạo ma trận đơn vị $I$ kích thước $(k-1) \times (k-1)$.
      - Tính $I\_minus\_correction = I - correction$.
      - Tính $B_{tl} = B_{k-1} \times I\_minus\_correction$:
        * Với mỗi $i$ từ $0$ đến $k-2$, mỗi $j$ từ $0$ đến $k-2$:
          $$B_{tl}[i][j] = \sum_{t=0}^{k-2} B_{k-1}[i][t] \times I\_minus\_correction[t][j]$$

   g. **Lắp ráp ma trận $B_k$:**
      - Tạo ma trận $B_k$ kích thước $k \times k$.
      - Gán khối trên trái: $B_k[0:k-1][0:k-1] = B_{tl}$.
      - Gán cột biên phải: $B_k[0:k-1][k-1] = \beta_{col}$.
      - Gán hàng biên dưới: $B_k[k-1][0:k-1] = \beta_{row}$.
      - Gán góc dưới phải: $B_k[k-1][k-1] = b_{kk}$.

   h. **Cập nhật cho vòng lặp tiếp theo:**
      - Gán $B_{k-1} = B_k$.

5. **Kết thúc đệ quy khối:**
   a. Sau khi $k = n$, $B_n$ chính là ma trận nghịch đảo $A^{-1}$.
   b. Trả về $B_n$.

### Phần C: Nghịch đảo gián tiếp cho ma trận không vuông

6. **Tính ma trận $M = A^T \times A$:**
   a. Tính $A^T$ (chuyển vị của $A$): $A^T$ có kích thước $n \times m$.
   b. Tạo ma trận $M$ kích thước $n \times n$, khởi tạo bằng $0$.
   c. Với mỗi $i$ từ $0$ đến $n-1$:
      - Với mỗi $j$ từ $0$ đến $n-1$:
        $$M[i][j] = \sum_{t=0}^{m-1} A^T[i][t] \times A[t][j]$$

7. **Tính nghịch đảo của $M$ bằng đệ quy khối:**
   a. Thực hiện các bước từ 3 đến 5 với ma trận $M$ thay vì $A$.
   b. Kết quả là $M^{-1}$ (ma trận kích thước $n \times n$).

8. **Tính ma trận giả nghịch đảo $A^{-1}$:**
   a. Tạo ma trận kết quả kích thước $n \times m$.
   b. Với mỗi $i$ từ $0$ đến $n-1$:
      - Với mỗi $j$ từ $0$ đến $m-1$:
        $$A^{-1}[i][j] = \sum_{t=0}^{n-1} M^{-1}[i][t] \times A^T[t][j]$$
   c. Trả về $A^{-1}$.

### Phần D: Hoàn tất

9. **Xuất kết quả:**
   a. In ma trận nghịch đảo (hoặc giả nghịch đảo) ra màn hình.
   b. Ghi kết quả ra file nếu yêu cầu.
   c. Kết thúc thuật toán.
