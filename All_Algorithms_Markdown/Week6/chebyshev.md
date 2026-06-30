# Nút Chebyshev (Chebyshev Nodes)

## Mô tả
Tạo các điểm nội suy Chebyshev trên đoạn $[a, b]$ để tránh hiện tượng Runge khi nội suy đa thức bậc cao.

## Công thức toán học

Nút Chebyshev được tính theo công thức:

$$t_i = \cos\left(\frac{(2i+1)\pi}{2n+2}\right), \quad i = 0, 1, \dots, n$$

Nếu $a \neq -1$ hoặc $b \neq 1$, chuyển đổi sang đoạn $[a, b]$:

$$x_i = \frac{1}{2}\left[(b-a)t_i + (b+a)\right]$$

## Thuật toán

**Đầu vào:** Số nút $n$, đoạn $[a, b]$.

**Đầu ra:** Danh sách các nút Chebyshev $x_i$ trên đoạn $[a, b]$.

1. **Tính số thứ tự các nút:** Với $i$ chạy từ 0 đến $n$, gán $i$ là chỉ số của nút thứ $i+1$.
   a. $i = 0$ ứng với nút đầu tiên.
   b. $i = n$ ứng với nút cuối cùng.
   c. Có tổng cộng $n+1$ nút Chebyshev.

2. **Tính giá trị cosin trung gian $t_i$:** Với mỗi $i$ từ 0 đến $n$, tính:
   $$t_i = \cos\left(\frac{(2i+1)\pi}{2n+2}\right)$$
   a. Tử số góc: $(2i+1)\pi$.
   b. Mẫu số góc: $2n+2$.
   c. Lấy cosin của tỷ số trên để được $t_i \in [-1, 1]$.

3. **Kiểm tra đoạn chuẩn hóa:** Xét hai trường hợp:
   a. **Nếu $a = -1$ và $b = 1$:** Đặt $x_i = t_i$, không cần biến đổi.
   b. **Nếu $a \neq -1$ hoặc $b \neq 1$:** Biến đổi tuyến tính từ $[-1, 1]$ sang $[a, b]$ theo công thức:
      $$x_i = \frac{1}{2}\left[(b-a)t_i + (b+a)\right]$$
      - Bước 1: Tính độ dài đoạn $b-a$.
      - Bước 2: Tính $b+a$.
      - Bước 3: Nhân $t_i$ với $(b-a)$.
      - Bước 4: Cộng với $(b+a)$.
      - Bước 5: Chia cho 2 để được $x_i$.

4. **Hiển thị kết quả:** Trả về các cột $i$, $t_i$, $x_i$ dưới dạng bảng dữ liệu.
