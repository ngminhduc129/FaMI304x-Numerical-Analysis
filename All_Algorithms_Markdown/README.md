# All Algorithms — Numerical Analysis

Thư mục này chứa toàn bộ mã nguồn và thuật toán của các phương pháp giải tích số trong môn **FaMI304x - Numerical Analysis**.

## Cấu trúc thư mục

```
All_Algorithms_Markdown/
├── README.md
├── Week1/               — Sai số (Error)
├── Week2/               — Tìm nghiệm hàm số (Roots of Equations)
├── Week3-NonLinear-System/  — Hệ phương trình phi tuyến
├── Week4/               — Hệ phương trình tuyến tính (Chính xác + Lặp)
├── Week5/               — Vector riêng & Giá trị riêng (Eigen)
├── Week6/               — Nội suy (Interpolation)
├── Week7/               — Đạo hàm & Tích phân số (Integration)
└── Week8/               — Phương trình vi phân (ODE)
```

## Danh sách các phương pháp theo tuần

### Week1 — Sai số
| # | Phương pháp | File |
|---|------------|------|
| 1 | Sai số (Error) | `error.md` |

### Week2 — Tìm nghiệm hàm số
| # | Phương pháp | File |
|---|------------|------|
| 1 | Phương pháp chia đôi (Bisection) | `bisection_method.md` |
| 2 | Phương pháp cát tuyến (Secant) | `secant_method.md` |
| 3 | Phương pháp Newton-Raphson | `newton_raphson.md` |
| 4 | Phương pháp lặp điểm cố định (Fixed-point) | `fixed_point_iteration.md` |
| 5 | Đa thức (Poly utils) | `poly_utils.md` |

### Week3 — Hệ phương trình phi tuyến
| # | Phương pháp | File |
|---|------------|------|
| 1 | Lặp điểm cố định cho hệ phi tuyến | `pp1_fixed_point_iteration.md` |
| 2 | Newton cho hệ phi tuyến | `pp2_newton.md` |

### Week4 — Hệ phương trình tuyến tính
#### Phần 1 — Phương pháp chính xác
| # | Phương pháp | File |
|---|------------|------|
| 1 | Khử Gauss | `gauss.md` |
| 2 | Gauss-Jordan | `gauss_jordan.md` |
| 3 | Phân rã LU | `lu_phan_ra.md` |
| 4 | Cholesky | `cholesky.md` |
| 5 | Viền quanh (Block matrix inversion) | `vienquanh.md` |

#### Phần 2 — Phương pháp lặp
| # | Phương pháp | File |
|---|------------|------|
| 1 | Lặp điểm cố định tuyến tính | `fixed_point_iteration_linear.md` |
| 2 | Jacobi | `jacobi.md` |
| 3 | Gauss-Seidel | `gauss_seidel.md` |

### Week5 — Vector riêng & Giá trị riêng
| # | Phương pháp | File |
|---|------------|------|
| 1 | Danilevsky | `danilevsky.md` |
| 2 | Lũy thừa + Khử (Power + Deflation) | `power_method_deflation.md` |
| 3 | SVD (v1) | `svd_v1.md` |
| 4 | SVD (v2 - dùng Power+Deflation) | `svd_v2.md` |

### Week6 — Nội suy
| # | Phương pháp | File |
|---|------------|------|
| 1 | Nút Chebyshev | `chebyshev.md` |
| 2 | Sơ đồ Horner | `horner.md` |
| 3 | Nội suy Lagrange | `lagrange.md` |
| 4 | Nhân tử Lagrange tiến | `nhan_tu_lagrange_tinh_tien.md` |
| 5 | Newton tiến/lùi (tỉ sai phân) | `noidungsuy_Newtontien_Newtonlui.md` |
| 6 | Newton cách đều tiến | `noi_suy_newton_cach_deu_tien.md` |
| 7 | Newton cách đều lùi | `noi_suy_newton_cach_deu_lui.md` |
| 8 | Gauss I | `noi_suy_gauss_I.md` |
| 9 | Gauss II | `noi_suy_gauss_II.md` |
| 10 | Stirling | `noi_suy_stirling.md` |
| 11 | Bessel | `noi_suy_bessel.md` |
| 12 | Nội suy ngược (dạng hàm) | `noi_suy_nguoc_1.md` |
| 13 | Nội suy ngược (dạng lặp) | `noi_suy_nguoc_2.md` |
| 14 | Spline tuyến tính | `noi_suy_tuyen_tinh_tung_khuc.md` |
| 15 | Spline bậc 2 | `noi_suy_bac2_tung_khuc.md` |
| 16 | Bình phương tối thiểu (Least squares) | `noi_suy_binhphuong_toithieu.md` |

### Week7 — Đạo hàm & Tích phân số
| # | Phương pháp | File |
|---|------------|------|
| 1 | Đạo hàm Lagrange | `lagrange_dao_ham.md` |
| 2 | Tích phân Lagrange | `lagrange_tich_phan.md` |
| 3 | Công thức hình thang (Trapezoid) | `hinh_thang.md` |
| 4 | Công thức Simpson 1/3 | `simpson.md` |
| 5 | Công thức Newton-Cotes tổng quát | `newton_cotes.md` |

### Week8 — Phương trình vi phân (ODE)
| # | Phương pháp | File |
|---|------------|------|
| 1 | Euler tiến (Forward) | `euler_forward.md` |
| 2 | Euler lùi (Backward) | `euler_backward.md` |
| 3 | Heun (cải tiến Euler) | `heun.md` |
| 4 | Runge-Kutta bậc 2 | `runge_kutta_2.md` |
| 5 | Runge-Kutta bậc 3 | `runge_kutta_3.md` |
| 6 | Runge-Kutta bậc 4 | `runge_kutta_4.md` |
| 7 | Adams-Bashforth bậc 2 | `adams_bashforth_2.md` |
| 8 | Adams-Bashforth bậc 3 | `adams_bashforth_3.md` |
| 9 | Adams-Bashforth bậc 4 | `adams_bashforth_4.md` |
| 10 | Adams-Moulton (ẩn) | `adams_moulton.md` |

## Quy ước
- Mỗi file `.md` bao gồm: tiêu đề, công thức toán (LaTeX), thuật toán từng bước, và mã nguồn Python.
- Công thức toán được viết bằng ký hiệu `$$...$$` (LaTeX).
