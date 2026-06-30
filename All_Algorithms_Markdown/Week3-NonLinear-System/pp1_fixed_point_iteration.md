# Phương Pháp Lặp Điểm Cố Định Cho Hệ Phi Tuyến

## Cơ Sở Lý Thuyết

Cho hệ phương trình phi tuyến dạng:

$$
\begin{cases}
x = g_1(x, y, z, \dots) \\
y = g_2(x, y, z, \dots) \\
z = g_3(x, y, z, \dots) \\
\vdots
\end{cases}
$$

Phương pháp lặp điểm cố định xây dựng dãy lặp:

$$
\mathbf{x}^{(k+1)} = \mathbf{G}(\mathbf{x}^{(k)})
$$

với $\mathbf{x}^{(k)} = (x_1^{(k)}, x_2^{(k)}, \dots, x_n^{(k)})^T$ và $\mathbf{G} = (g_1, g_2, \dots, g_n)^T$.

## Điều Kiện Hội Tụ

Nếu $\mathbf{G}$ là ánh xạ co với hệ số co $q$ ($0 \le q < 1$):

$$
\|\mathbf{G}(\mathbf{x}) - \mathbf{G}(\mathbf{y})\| \le q \|\mathbf{x} - \mathbf{y}\|
$$

thì dãy lặp hội tụ đến nghiệm duy nhất.

### Sai số tiên nghiệm:

$$
\|\mathbf{x}^{(k)} - \mathbf{x}^*\| \le \frac{q^k}{1-q} \|\mathbf{x}^{(1)} - \mathbf{x}^{(0)}\|
$$

### Sai số hậu nghiệm:

$$
\|\mathbf{x}^{(k)} - \mathbf{x}^*\| \le \frac{q}{1-q} \|\mathbf{x}^{(k)} - \mathbf{x}^{(k-1)}\|
$$

## Tiêu Chuẩn Dừng

### Tiêu chuẩn sai số tuyệt đối:

$$
\|\mathbf{x}^{(k)} - \mathbf{x}^{(k-1)}\| < \varepsilon \cdot \frac{1-q}{q}
$$

với $\varepsilon$ là sai số cho phép.

### Tiêu chuẩn sai số tương đối:

$$
\frac{\|\mathbf{x}^{(k)} - \mathbf{x}^{(k-1)}\|}{\|\mathbf{x}^{(k)}\|} < \eta \cdot \frac{1-q}{q}
$$

với $\eta$ là sai số tương đối cho phép.

## Các Chuẩn Thường Dùng

- **Chuẩn 1 (norm=1):** $\|\mathbf{x}\|_1 = \sum_{i=1}^{n} |x_i|$
- **Chuẩn vô cùng (norm=-1):** $\|\mathbf{x}\|_\infty = \max_{1 \le i \le n} |x_i|$

## Thuật Toán

### fixed_point_recursion_absolute (tiêu chuẩn sai số tuyệt đối)

**Đầu vào:**
- `initial_values`: vector khởi tạo $\mathbf{x}^{(0)} = (x_1^{(0)}, x_2^{(0)}, \dots, x_n^{(0)})^T$
- `q`: hệ số co, thỏa mãn $0 \le q < 1$
- `eps`: sai số cho phép $\varepsilon > 0$
- `*funcs`: danh sách các hàm $g_1, g_2, \dots, g_n$
- `norm`: loại chuẩn (1 cho chuẩn 1, -1 cho chuẩn vô cùng)

**Đầu ra:** Nghiệm xấp xỉ $\mathbf{x}^*$ của hệ phương trình.

**Các bước:**

1. **Tính ngưỡng dừng từ sai số cho phép:**
   - Đặt $\text{new\_eps} = \text{eps} \times \dfrac{1 - q}{q}$
   - Giá trị $\text{new\_eps}$ là ngưỡng để so sánh trực tiếp với $\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\|$

2. **Khởi tạo giá trị lặp ban đầu và bảng kết quả:**
   - Gán $\mathbf{x} = \text{initial\_values}$, đây là $\mathbf{x}^{(0)}$
   - Tạo một bảng (danh sách) để lưu kết quả từng bước lặp, bao gồm: số thứ tự lần lặp, vector $\mathbf{x}$ hiện tại, vector chênh lệch, và tổng chênh lệch

3. **Thực hiện vòng lặp (số lần lặp không giới hạn, dừng khi thỏa mãn tiêu chuẩn):**
   - **Bước 3.1:** Tính vector xấp xỉ mới:
     - Với mỗi thành phần $i = 1, 2, \dots, n$, tính $x_{\text{new}, i} = g_i(\mathbf{x})$
     - Kết quả là $\mathbf{x}_\text{new} = (x_{\text{new}, 1}, x_{\text{new}, 2}, \dots, x_{\text{new}, n})$
   - **Bước 3.2:** Tính vector chênh lệch tuyệt đối giữa hai lần lặp liên tiếp:
     - Với mỗi thành phần $i = 1, 2, \dots, n$, tính $\text{diff}_i = |x_{\text{new}, i} - x_i|$
   - **Bước 3.3:** Tính tổng chênh lệch theo chuẩn đã chọn:
     - Nếu `norm = 1`: $\text{total\_diff} = \sum_{i=1}^{n} \text{diff}_i$ (chuẩn 1)
     - Nếu `norm = -1`: $\text{total\_diff} = \max_{1 \le i \le n} \text{diff}_i$ (chuẩn vô cùng)
   - **Bước 3.4:** Lưu kết quả của lần lặp hiện tại vào bảng (số thứ tự, $\mathbf{x}_\text{new}$, $\text{diff}$, $\text{total\_diff}$)
   - **Bước 3.5:** Cập nhật giá trị lặp: $\mathbf{x} = \mathbf{x}_\text{new}$
   - **Bước 3.6:** Kiểm tra điều kiện dừng:
     - Nếu $\text{total\_diff} < \text{new\_eps}$ thì **kết thúc vòng lặp** (nghiệm đã đạt độ chính xác yêu cầu)
     - Ngược lại, quay lại Bước 3.1 để tiếp tục lặp

4. **In bảng kết quả:** Hiển thị bảng kết quả lưu được trong quá trình lặp (gồm tất cả các lần lặp đã thực hiện)

5. **Trả về nghiệm xấp xỉ:** Giá trị $\mathbf{x}$ hiện tại chính là nghiệm xấp xỉ $\mathbf{x}^*$ của hệ phương trình

---

### fixed_point_recursion_relative (tiêu chuẩn sai số tương đối)

**Đầu vào:**
- `initial_values`: vector khởi tạo $\mathbf{x}^{(0)} = (x_1^{(0)}, x_2^{(0)}, \dots, x_n^{(0)})^T$
- `q`: hệ số co, thỏa mãn $0 \le q < 1$
- `eta`: sai số tương đối cho phép $\eta > 0$
- `*funcs`: danh sách các hàm $g_1, g_2, \dots, g_n$
- `norm`: loại chuẩn (1 cho chuẩn 1, -1 cho chuẩn vô cùng)

**Đầu ra:** Nghiệm xấp xỉ $\mathbf{x}^*$ của hệ phương trình.

**Các bước:**

1. **Tính ngưỡng dừng từ sai số tương đối cho phép:**
   - Đặt $\text{new\_eta} = \text{eta} \times \dfrac{1 - q}{q}$
   - Giá trị $\text{new\_eta}$ là ngưỡng để so sánh với sai số tương đối $\dfrac{\|\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}\|}{\|\mathbf{x}^{(k+1)}\|}$

2. **Khởi tạo giá trị lặp ban đầu và bảng kết quả:**
   - Gán $\mathbf{x} = \text{initial\_values}$, đây là $\mathbf{x}^{(0)}$
   - Tạo một bảng (danh sách) để lưu kết quả từng bước lặp

3. **Thực hiện vòng lặp với chỉ số $i = 0, 1, 2, \dots$ (dừng khi thỏa mãn tiêu chuẩn):**
   - **Bước 3.1:** Tính vector xấp xỉ mới:
     - Với mỗi thành phần $i = 1, 2, \dots, n$, tính $x_{\text{new}, i} = g_i(\mathbf{x})$
     - Kết quả là $\mathbf{x}_\text{new} = (x_{\text{new}, 1}, x_{\text{new}, 2}, \dots, x_{\text{new}, n})$
   - **Bước 3.2:** Tính vector chênh lệch tuyệt đối:
     - Với mỗi thành phần $i = 1, 2, \dots, n$, tính $\text{diff}_i = |x_{\text{new}, i} - x_i|$
   - **Bước 3.3:** Tính vector trị tuyệt đối của xấp xỉ mới:
     - Với mỗi thành phần $i = 1, 2, \dots, n$, tính $\text{itself}_i = |x_{\text{new}, i}|$
   - **Bước 3.4:** Tính tổng chênh lệch và tổng trị tuyệt đối theo chuẩn đã chọn:
     - Nếu `norm = 1`:
       - $\text{total\_diff} = \sum_{i=1}^{n} \text{diff}_i$ (chuẩn 1)
       - $\text{total\_itself} = \sum_{i=1}^{n} \text{itself}_i$ (chuẩn 1)
     - Nếu `norm = -1`:
       - $\text{total\_diff} = \max_{1 \le i \le n} \text{diff}_i$ (chuẩn vô cùng)
       - $\text{total\_itself} = \max_{1 \le i \le n} \text{itself}_i$ (chuẩn vô cùng)
   - **Bước 3.5:** Nếu $i \ne 0$ (không phải lần lặp đầu tiên), tính sai số tương đối:
     - $\text{relative\_diff} = \dfrac{\text{total\_diff}}{\text{total\_itself}}$
   - **Bước 3.6:** Lưu kết quả của lần lặp hiện tại vào bảng (số thứ tự, $\mathbf{x}_\text{new}$, $\text{diff}$, $\text{total\_diff}$, $\text{total\_itself}$, $\text{relative\_diff}$)
   - **Bước 3.7:** Cập nhật giá trị lặp: $\mathbf{x} = \mathbf{x}_\text{new}$
   - **Bước 3.8:** Kiểm tra điều kiện dừng:
     - Nếu $i \ne 0$ **và** $\text{relative\_diff} < \text{new\_eta}$ thì **kết thúc vòng lặp**
     - Ngược lại, tăng $i$ lên 1 và quay lại Bước 3.1 để tiếp tục lặp

4. **In bảng kết quả:** Hiển thị bảng kết quả lưu được trong quá trình lặp (gồm tất cả các lần lặp đã thực hiện)

5. **Trả về nghiệm xấp xỉ:** Giá trị $\mathbf{x}$ hiện tại chính là nghiệm xấp xỉ $\mathbf{x}^*$ của hệ phương trình
