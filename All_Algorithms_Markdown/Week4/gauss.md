# Phương pháp khử Gauss (Gaussian Elimination)

## Công thức toán học

### 1. Khử tiến (Forward Elimination)

Cho hệ phương trình tuyến tính $Ax = B$, ma trận mở rộng $[A|B]$. Với mỗi cột pivot $i$, phần tử khử tại hàng $j$ ($j > i$):

$$A_{j,:} = A_{j,:} - \frac{A_{j, i}}{A_{i, i}} \cdot A_{i,:}$$

Trong đó $A_{i,i}$ là phần tử trụ (pivot) tại hàng $i$, cột $i$.

### 2. Sắp xếp lại hàng (Row Sorting)

Để đảm bảo cấu trúc bậc thang, các hàng được sắp xếp theo chỉ số cột chứa phần tử khác không đầu tiên:

$$\text{index}(hàng) = \min\{col : |A_{row, col}| > 10^{-10}\}$$

### 3. Thế ngược (Back Substitution)

Sau khi đưa ma trận về dạng bậc thang, nghiệm được tính từ dưới lên:

$$x_i = \frac{B_i - \sum_{k=i+1}^{n} A_{i,k} \cdot x_k}{A_{i,i}}$$

### 4. Xử lý biến tự do

Khi một cột không có pivot, cột đó được chuyển sang vế phải:

$$B = [B \; | \; -A_{:,i}]$$

---

## Thuật toán

**Mục tiêu:** Giải hệ phương trình tuyến tính $Ax = B$ bằng phương pháp khử Gauss.
**Đầu vào:** Ma trận $A_{m \times n}$, vector $B_{m \times 1}$ đọc từ file đầu vào.
**Đầu ra:** Ma trận nghiệm $X$ (dạng cột hoặc ma trận nếu hệ có vô số nghiệm).

### Phần A: Tiền xử lý và xây dựng ma trận mở rộng

**Bước 1:** Đọc dữ liệu đầu vào
   - **Bước 1.1:** Đọc ma trận $A$ từ file. Xác định số hàng $m$ và số cột $n$.
   - **Bước 1.2:** Đọc vector $B$ từ file. Kiểm tra $B$ có $m$ phần tử.
   - **Bước 1.3:** Nếu kích thước không tương thích, thông báo lỗi và kết thúc.

**Bước 2:** Xây dựng ma trận mở rộng
   - **Bước 2.1:** Tạo ma trận $AB$ kích thước $m \times (n + 1)$.
   - **Bước 2.2:** Gán $AB[:, :n] = A$: sao chép $A$ vào các cột đầu của $AB$.
   - **Bước 2.3:** Gán $AB[:, n] = B$: đặt $B$ vào cột cuối của $AB$.

### Phần B: Khử tiến (Forward Elimination)

**Bước 3:** Tạo mảng chỉ số pivot ban đầu
   - **Bước 3.1:** Với mỗi hàng $i$ từ $0$ đến $m-1$, tìm chỉ số cột $pivot\_idx[i]$ là cột đầu tiên mà $|AB[i][c]| > 10^{-10}$.
   - **Bước 3.2:** Nếu toàn bộ hàng là số 0, gán $pivot\_idx[i] = \infty$.

**Bước 4:** Sắp xếp hàng theo chỉ số pivot
   - **Bước 4.1:** Sắp xếp các hàng của $AB$ theo thứ tự tăng dần của $pivot\_idx$.
   - **Bước 4.2:** Cập nhật lại mảng $pivot\_idx$ sau khi sắp xếp.

**Bước 5:** Với mỗi hàng $i$ từ $0$ đến $m-1$ (hàng chủ)

   **Bước 5.1:** Nếu $pivot\_idx[i] \ge n$ (hàng toàn số 0), bỏ qua hàng này.
   
   **Bước 5.2:** Đặt $pivot\_col = pivot\_idx[i]$.
   
   **Bước 5.3:** Với mỗi hàng $j$ từ $i+1$ đến $m-1$ (hàng bên dưới):
   - Nếu $pivot\_idx[j] == pivot\_col$:
     * Tính hệ số khử: $factor = AB[j][pivot\_col] / AB[i][pivot\_col]$.
     * Với mỗi cột $k$ từ $pivot\_col$ đến $n$:
       $$AB[j][k] = AB[j][k] - factor \times AB[i][k]$$
     * Cập nhật lại $pivot\_idx[j]$: tìm cột khác không đầu tiên mới.
   
   **Bước 5.4:** Sau khi khử xong, sắp xếp lại các hàng từ $i+1$ đến $m-1$ theo $pivot\_idx$.
   
   **Bước 5.5:** Cập nhật lại mảng $pivot\_idx$ cho các hàng từ $i+1$ đến $m-1$.

**Bước 6:** Kiểm tra hệ vô nghiệm
   - **Bước 6.1:** Duyệt từng hàng $i$ từ $0$ đến $m-1$.
   - **Bước 6.2:** Nếu $pivot\_idx[i] == n$ (pivot tại cột $B$), tức hàng $[0\;0\;\dots\;0\;|\;b \neq 0]$:
     * Kết luận hệ vô nghiệm, kết thúc thuật toán.

### Phần C: Xác định biến tự do và thế ngược

**Bước 7:** Xác định biến cơ sở và biến tự do
   - **Bước 7.1:** Tạo danh sách $pivot\_cols$ chứa các $pivot\_idx[i]$ hợp lệ.
   - **Bước 7.2:** Tạo danh sách $free\_cols$ là các cột từ $0$ đến $n-1$ không có trong $pivot\_cols$.
   - **Bước 7.3:** Ghi nhận $r = |pivot\_cols|$ (số biến cơ sở) và $f = n - r$ (số biến tự do).

**Bước 8:** Xử lý biến tự do (nếu $f > 0$)
   - **Bước 8.1:** Với mỗi biến tự do tại cột $free\_cols[t]$:
     * Thêm một cột mới vào bên phải $AB$.
     * Cột mới là $-AB[:, free\_cols[t]]$ (chuyển cột sang vế phải).
   - **Bước 8.2:** Cập nhật số cột của ma trận nghiệm.

**Bước 9:** Thế ngược (Back Substitution)

   **Bước 9.1:** Tạo ma trận nghiệm $X$ kích thước $n \times (1 + f)$.
   
   **Bước 9.2:** Khởi tạo $X$ với toàn bộ phần tử bằng $0$.
   
   **Bước 9.3:** Với mỗi biến tự do tại cột $free\_cols[t]$:
   - Gán $X[free\_cols[t]][0] = 0$ (nghiệm tổng quát).
   - Gán $X[free\_cols[t]][t+1] = 1$ (hệ số tham số).
   
   **Bước 9.4:** Với mỗi hàng $i$ từ $r-1$ xuống $0$:
   - Đặt $row$ là hàng thực tế thứ $i$.
   - Đặt $pc = pivot\_cols[i]$.
   - **Tính nghiệm cột tổng quát (cột 0):**
     * $sum\_const = AB[row][n]$.
     * Với mỗi cột $k$ từ $pc+1$ đến $n-1$: nếu $AB[row][k] \neq 0$, $sum\_const = sum\_const - AB[row][k] \times X[k][0]$.
     * $X[pc][0] = sum\_const / AB[row][pc]$.
   - **Tính nghiệm các cột tham số (cột $t+1$):**
     * Với mỗi tham số $t$ từ $0$ đến $f-1$:
       - $sum\_param = AB[row][n + 1 + t]$ (nếu cột mở rộng tồn tại).
       - Với mỗi cột $k$ từ $pc+1$ đến $n-1$: nếu $AB[row][k] \neq 0$, $sum\_param = sum\_param - AB[row][k] \times X[k][t+1]$.
       - $X[pc][t+1] = sum\_param / AB[row][pc]$.

### Phần D: Hoàn tất

**Bước 10:** Trả về kết quả
   - **Bước 10.1:** Nếu $f = 0$: hệ có nghiệm duy nhất, trả về vector $X[:, 0]$.
   - **Bước 10.2:** Nếu $f > 0$: hệ có vô số nghiệm, trả về ma trận $X$ với $f$ tham số tự do.
   - **Bước 10.3:** In kết quả ra màn hình hoặc ghi file.
