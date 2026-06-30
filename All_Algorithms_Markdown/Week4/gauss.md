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

Khi một cột không có pivot (phần tử trên đường chéo chính bằng 0), cột đó được chuyển sang vế phải:

$$B = [B \; | \; -A_{:,i}]$$

## Thuật toán

**Đầu vào:** Ma trận $A_{m \times n}$, vector $B_{m \times 1}$ đọc từ file đầu vào.  
**Đầu ra:** Ma trận nghiệm $X$ (dạng cột hoặc ma trận nếu hệ có vô số nghiệm có tham số).

### Phần A: Tiền xử lý và xây dựng ma trận mở rộng

1. **Đọc dữ liệu đầu vào:**
   a. Đọc ma trận $A$ từ file. Xác định số hàng $m$ và số cột $n$ của $A$.
   b. Đọc vector $B$ từ file. Xác định số phần tử của $B$, phải bằng $m$.
   c. Kiểm tra tính hợp lệ: nếu kích thước không tương thích, thông báo lỗi và kết thúc.

2. **Xây dựng ma trận mở rộng:**
   a. Tạo ma trận $AB$ có kích thước $m \times (n + 1)$.
   b. Gán $AB[:, :n] = A$: sao chép toàn bộ ma trận $A$ vào các cột đầu của $AB$.
   c. Gán $AB[:, n] = B$: đặt vector $B$ vào cột cuối cùng của $AB$.

### Phần B: Khử tiến (Forward Elimination)

3. **Tạo mảng chỉ số pivot ban đầu:**
   a. Với mỗi hàng $i$ từ $0$ đến $m-1$, tìm chỉ số cột $pivot\_idx[i]$ là cột đầu tiên mà $|AB[i][c]| > 10^{-10}$.
   b. Nếu toàn bộ hàng là số 0, gán $pivot\_idx[i] = \infty$ (hoặc một số lớn hơn $n$).

4. **Sắp xếp hàng theo chỉ số pivot:**
   a. Sắp xếp các hàng của $AB$ theo thứ tự tăng dần của $pivot\_idx$ (hàng nào có pivot ở cột nhỏ hơn thì xếp lên trên).
   b. Cập nhật lại mảng $pivot\_idx$ sau khi sắp xếp.

5. **Với mỗi hàng $i$ từ $0$ đến $m-1$ (hàng chủ):**
   a. Nếu $pivot\_idx[i] \ge n$ (hàng toàn số 0), bỏ qua hàng này và tiếp tục.
   b. Đặt $pivot\_col = pivot\_idx[i]$.
   c. **Với mỗi hàng $j$ từ $i+1$ đến $m-1$ (hàng bên dưới):**
      - Nếu $pivot\_idx[j] == pivot\_col$ (hai hàng có cùng cột pivot):
        * Tính hệ số khử: $factor = AB[j][pivot\_col] / AB[i][pivot\_col]$.
        * Với mỗi cột $k$ từ $pivot\_col$ đến $n$:
          $$AB[j][k] = AB[j][k] - factor \times AB[i][k]$$
        * Cập nhật lại $pivot\_idx[j]$: tìm cột khác không đầu tiên mới của hàng $j$.
      - Nếu $pivot\_idx[j] \neq pivot\_col$, giữ nguyên hàng $j$.
   d. Sau khi khử xong tất cả các hàng $j > i$, sắp xếp lại toàn bộ các hàng từ $i+1$ đến $m-1$ theo $pivot\_idx$ để đảm bảo dạng bậc thang.
   e. Cập nhật lại mảng $pivot\_idx$ cho các hàng từ $i+1$ đến $m-1$.

6. **Kiểm tra hệ vô nghiệm:**
   a. Duyệt từng hàng $i$ từ $0$ đến $m-1$.
   b. Nếu $pivot\_idx[i] == n$ (pivot nằm ngay tại cột của $B$), tức là hàng có dạng $[0\;0\;\dots\;0\;|\;b \neq 0]$:
      - Kết luận hệ phương trình vô nghiệm.
      - Trả về kết quả thông báo vô nghiệm và kết thúc thuật toán.

### Phần C: Xác định biến tự do và thế ngược

7. **Xác định biến cơ sở và biến tự do:**
   a. Tạo danh sách $pivot\_cols$ chứa các $pivot\_idx[i]$ hợp lệ (khác $\infty$ và khác $n$).
   b. Tạo danh sách $free\_cols$ là các cột từ $0$ đến $n-1$ không có trong $pivot\_cols$.
   c. Ghi nhận số lượng biến cơ sở $r = |pivot\_cols|$ và số biến tự do $f = n - r$.

8. **Xử lý biến tự do (nếu có):**
   a. Nếu $f > 0$:
      - Với mỗi biến tự do tại cột $free\_cols[t]$:
        * Thêm một cột mới vào bên phải ma trận $AB$ (mở rộng $B$).
        * Cột mới này là vector $-AB[:, free\_cols[t]]$ (chuyển cột tương ứng với biến tự do sang vế phải).
        * Đánh dấu biến tự do $t$ tương ứng với tham số tự do.
   b. Cập nhật số cột của ma trận nghiệm.

9. **Thế ngược (Back Substitution):**
   a. Tạo ma trận nghiệm $X$ có kích thước $n \times (1 + f)$ (một cột cho nghiệm tổng quát, $f$ cột cho các tham số).
   b. Khởi tạo $X$ với toàn bộ phần tử bằng $0$.
   c. **Với mỗi biến tự do tại cột $free\_cols[t]$:**
      - Gán $X[free\_cols[t]][0] = 0$ (nghiệm tổng quát).
      - Gán $X[free\_cols[t]][t+1] = 1$ (hệ số của tham số).
   d. **Với mỗi hàng $i$ từ $r-1$ xuống $0$ (từ dưới lên trên):**
      - Đặt $row$ là hàng thực tế thứ $i$ trong ma trận bậc thang (đã sắp xếp).
      - Đặt $pc = pivot\_cols[i]$.
      - **Tính nghiệm cho cột tổng quát (cột 0):**
        * $sum\_const = AB[row][n]$ (giá trị ở cột $B$).
        * Với mỗi cột $k$ từ $pc+1$ đến $n-1$:
          - Nếu $AB[row][k] \neq 0$:
            $sum\_const = sum\_const - AB[row][k] \times X[k][0]$
        * $X[pc][0] = sum\_const / AB[row][pc]$
      - **Tính nghiệm cho các cột tham số (cột $t+1$ với $t = 0 \dots f-1$):**
        * Với mỗi tham số $t$ từ $0$ đến $f-1$:
          - $sum\_param = 0$ (nếu biến tự do nằm bên phải thì có thể có giá trị từ cột mở rộng)
          - Nếu cột mở rộng tương ứng tồn tại:
            $sum\_param = AB[row][n + 1 + t]$ (giá trị từ cột đã chuyển sang)
          - Với mỗi cột $k$ từ $pc+1$ đến $n-1$:
            Nếu $AB[row][k] \neq 0$:
              $sum\_param = sum\_param - AB[row][k] \times X[k][t+1]$
          - $X[pc][t+1] = sum\_param / AB[row][pc]$

### Phần D: Hoàn tất

10. **Trả về kết quả:**
    a. Nếu $f = 0$: hệ có nghiệm duy nhất, trả về vector nghiệm $X[:, 0]$.
    b. Nếu $f > 0$: hệ có vô số nghiệm, trả về ma trận nghiệm $X$ với $f$ tham số tự do.
    c. In ma trận nghiệm ra màn hình hoặc ghi ra file theo yêu cầu.
