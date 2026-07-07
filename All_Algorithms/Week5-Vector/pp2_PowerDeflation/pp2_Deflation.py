# =============================================================================
# pp2_Deflation.py - Phương pháp xuống thang (Deflation)
#
# Chức năng: Loại bỏ ảnh hưởng của một cặp trị riêng đã tìm được khỏi ma trận,
#            cho phép tìm các trị riêng còn lại bằng phương pháp lũy thừa.
#            Kết hợp với pp2_PowerMethod.py để tìm tất cả trị riêng.
#
# Các hàm chính:
#   deflate_once(A, lam, v, ...)            - xuống thang một bước
#   compute_all_eigenpairs(A, x0)           - tìm tất cả trị riêng
#   run_all_cases(A, x0)                    - chạy cả 3 case
#
# Input: Đọc từ PWDF_input_A.txt, PWDF_input_A3.txt
# Cách dùng: python pp2_Deflation.py
# =============================================================================
import contextlib
from pathlib import Path
from fractions import Fraction
from typing import List, Tuple, Union
import numpy as np
import pandas as pd
from pp2_PowerMethod import power_method, power_method_case2, power_method_case3, \
    input_matrix, output_matrix, verify_eigenpair

__dir__ = Path(__file__).parent.resolve()

pd.set_option('display.precision', 12)
pd.set_option('display.width', 300)
pd.set_option('display.max_columns', None)


def deflate_once(A: np.ndarray, lam: float, v: np.ndarray, x0_left: np.ndarray,
                 tol: float, max_iter: int, norm_ord: Union[int, float],
                 precision: int) -> np.ndarray:
    """
    Deflate matrix A by removing the effect of eigenpair (lam, v) for non-symmetric A:
    1. Compute left eigenvector w of A: solve A^T w = lam w via power_method
    2. Normalize so that w^T v = 1
    3. A_new = A - lam * np.outer(v, w)
    """
    # Step 1: compute left eigenvector of A
    _, w = power_method(
        A.T, x0_left, tol=tol, max_iter=max_iter,
        norm_ord=norm_ord, precision=precision,
        display=True
    )
    # Step 2: normalize w so that w^T v = 1
    scale = w.dot(v)
    if abs(scale) < 1e-12:
        raise ValueError("Vector riêng trái và phải gần như trực giao; không thể xuống thang.")
    w = w / scale
    # Step 3: rank-1 deflation
    return A - lam * np.outer(v, w)


def compute_all_eigenpairs(
    A: np.ndarray,
    x0: np.ndarray,
    tol: float = 1e-6,
    max_iter: int = 200,
    norm_ord: Union[int, float] = 2,
    precision: int = 8
) -> Tuple[List[float], List[np.ndarray]]:
    """
    Compute all real eigenvalues and eigenvectors of A via deflation.

    Returns:
      eigs: list of eigenvalues (sorted by absolute value descending)
      vecs: list of corresponding eigenvectors
    """
    n = A.shape[0]
    M = A.copy().astype(float)
    y = x0.astype(float).flatten()
    eigs = []
    vecs = []

    for k in range(n):
        print(f"Đang tìm cặp trị riêng {k+1}: kích thước ma trận {M.shape}")
        lam, v = power_method(
            M, y, tol=tol, max_iter=max_iter,
            norm_ord=norm_ord, precision=precision, display=True
        )
        eigs.append(lam)
        vecs.append(v)
        M = deflate_once(M, lam, v, y, tol, max_iter, norm_ord, precision)
        y = x0.astype(float).flatten()

    pairs = sorted(zip(eigs, vecs), key=lambda ev: -abs(ev[0]))
    eigs_sorted, vecs_sorted = zip(*pairs)

    print("Tất cả trị riêng và vector riêng (sắp xếp):")
    for i, (lam, vec) in enumerate(zip(eigs_sorted, vecs_sorted), start=1):
        vec_str = ", ".join(f"{x:.{precision}f}" for x in vec)
        print(f"  λ{i} = {lam:.{precision}f}, v = [{vec_str}]")

    return list(eigs_sorted), list(vecs_sorted)


def run_all_cases(A: np.ndarray, x0: np.ndarray,
                  tol: float = 1e-6, max_iter: int = 200,
                  norm_ord: Union[int, float] = 2, precision: int = 8):
    n = A.shape[0]
    threshold = 1e-4

    print("=" * 70)
    print("CHẠY CẢ 3 CASE TRÊN CÙNG MA TRẬN")
    print("=" * 70)

    results = {}

    # ===== CASE 1 =====
    print("\n" + "=" * 70)
    print("CASE 1 - Lũy thừa thường (trị riêng trội thực, bội đơn)")
    print("=" * 70)
    try:
        M = A.copy().astype(float)
        y = x0.copy().astype(float).flatten()
        eigs1, vecs1 = [], []
        for k in range(n):
            lam, v = power_method(M, y, tol=tol, max_iter=max_iter,
                                  norm_ord=norm_ord, precision=precision,
                                  display=True)
            eigs1.append(lam)
            vecs1.append(v)
            M = deflate_once(M, lam, v, y, tol, max_iter, norm_ord, precision)
            y = x0.copy().astype(float).flatten()

        pairs = sorted(zip(eigs1, vecs1), key=lambda ev: -abs(ev[0]))
        eigs1, vecs1 = zip(*pairs)

        print(f"✅ Thành công - {n} trị riêng:")
        for i, (lam, v) in enumerate(zip(eigs1, vecs1), 1):
            r = verify_eigenpair(A, lam, v)
            stt = "✅" if r < threshold else "⚠️"
            print(f"  λ{i} = {lam:.{precision}f}  residual = {r:.2e}  {stt}")
        results['case1'] = {'status': '✅', 'eigs': list(eigs1), 'vecs': list(vecs1)}
    except Exception as e:
        print(f"❌ Thất bại: {e}")
        results['case1'] = {'status': '❌', 'error': str(e)}

    # ===== CASE 2 =====
    print("\n" + "=" * 70)
    print("CASE 2 - Hai trị trội thực, bằng trị tuyệt đối, trái dấu (λ₂ = -λ₁)")
    print("=" * 70)
    try:
        lam2a, v2a, v2b = power_method_case2(
            A, x0.copy().flatten(), tol=tol, max_iter=max_iter,
            norm_ord=norm_ord, precision=precision, display=False)
        r2a = verify_eigenpair(A, lam2a, v2a)
        r2b = verify_eigenpair(A, -lam2a, v2b)

        if r2a < threshold and r2b < threshold:
            print(f"✅ Phù hợp: λ₁ = {lam2a:.{precision}f}, λ₂ = {-lam2a:.{precision}f}")
            print(f"   residual₁ = {r2a:.2e}, residual₂ = {r2b:.2e}")

            M = A.copy().astype(float)
            y = x0.copy().astype(float).flatten()
            M = deflate_once(M, lam2a, v2a, y, tol, max_iter, norm_ord, precision)
            M = deflate_once(M, -lam2a, v2b, y, tol, max_iter, norm_ord, precision)

            eigs2 = [lam2a, -lam2a]
            vecs2 = [v2a, v2b]
            y = x0.copy().astype(float).flatten()

            for k in range(n - 2):
                try:
                    lam, v = power_method(M, y, tol=tol, max_iter=max_iter,
                                          norm_ord=norm_ord, precision=precision,
                                          display=True)
                    eigs2.append(lam)
                    vecs2.append(v)
                    M = deflate_once(M, lam, v, y, tol, max_iter, norm_ord, precision)
                    y = x0.copy().astype(float).flatten()
                except Exception:
                    break

            pairs = sorted(zip(eigs2, vecs2), key=lambda ev: -abs(ev[0]))
            eigs2, vecs2 = zip(*pairs)

            print(f"   → {len(eigs2)}/{n} trị riêng (deflation):")
            for i, (lam, v) in enumerate(zip(eigs2, vecs2), 1):
                r = verify_eigenpair(A, lam, v)
                stt = "✅" if r < threshold else "⚠️"
                print(f"     λ{i} = {lam:.{precision}f}  residual = {r:.2e}  {stt}")
            results['case2'] = {'status': '✅', 'eigs': list(eigs2), 'vecs': list(vecs2)}
        else:
            print(f"❌ Không phù hợp (residual: {r2a:.2e}, {r2b:.2e})")
            results['case2'] = {'status': '❌', 'error': f'residual lớn: {r2a:.2e}, {r2b:.2e}'}
    except Exception as e:
        print(f"❌ Thất bại: {e}")
        results['case2'] = {'status': '❌', 'error': str(e)}

    # ===== CASE 3 =====
    print("\n" + "=" * 70)
    print("CASE 3 - Hai trị trội là phức liên hợp (λ₂ = λ̅₁)")
    print("=" * 70)
    try:
        lam3a, lam3b, df3 = power_method_case3(
            A, x0.copy().astype(float).flatten(), tol=tol,
            max_iter=max_iter // 2, norm_ord=norm_ord,
            precision=precision, display=False)
        v3a = df3['v1'].values
        v3b = df3['v2'].values
        r3a = verify_eigenpair(A, lam3a, v3a)
        r3b = verify_eigenpair(A, lam3b, v3b)

        if r3a < threshold and r3b < threshold:
            print(f"✅ Phù hợp: λ₁ = {lam3a:.{precision}f}, λ₂ = {lam3b:.{precision}f}")
            print(f"   residual₁ = {r3a:.2e}, residual₂ = {r3b:.2e}")

            M = A.copy().astype(complex)
            y = x0.copy().astype(complex).flatten()
            M = deflate_once(M, lam3a, v3a, y, tol, max_iter, norm_ord, precision)
            M = deflate_once(M, lam3b, v3b, y, tol, max_iter, norm_ord, precision)

            eigs3 = [lam3a, lam3b]
            vecs3 = [v3a, v3b]
            y = x0.copy().astype(complex).flatten()

            for k in range(n - 2):
                try:
                    lam, v = power_method(M, y, tol=tol, max_iter=max_iter,
                                          norm_ord=norm_ord, precision=precision,
                                          display=True)
                    eigs3.append(lam)
                    vecs3.append(v)
                    M = deflate_once(M, lam, v, y, tol, max_iter, norm_ord, precision)
                    y = x0.copy().astype(complex).flatten()
                except Exception:
                    break

            pairs = sorted(zip(eigs3, vecs3), key=lambda ev: -abs(ev[0]))
            eigs3, vecs3 = zip(*pairs)

            print(f"   → {len(eigs3)}/{n} trị riêng (deflation):")
            for i, (lam, v) in enumerate(zip(eigs3, vecs3), 1):
                r = verify_eigenpair(A, lam, v)
                stt = "✅" if r < threshold else "⚠️"
                if np.isreal(lam):
                    fstr = f"{lam.real:.{precision}f}"
                else:
                    fstr = f"{lam.real:.{precision}f}{lam.imag:+.{precision}f}j"
                print(f"     λ{i} = {fstr}  residual = {r:.2e}  {stt}")
            results['case3'] = {'status': '✅', 'eigs': list(eigs3), 'vecs': list(vecs3)}
        else:
            print(f"❌ Không phù hợp (residual: {r3a:.2e}, {r3b:.2e})")
            results['case3'] = {'status': '❌', 'error': f'residual lớn: {r3a:.2e}, {r3b:.2e}'}
    except Exception as e:
        print(f"❌ Thất bại: {e}")
        results['case3'] = {'status': '❌', 'error': str(e)}

    print("\n" + "=" * 70)
    print("KẾT THÚC")
    print("=" * 70)
    return results


# ==== DEMO ====
if __name__ == "__main__":
    output_path = str(__dir__ / "Deflation_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        A = input_matrix('PWDF_input_A.txt', convert_fractions=False)
        x0 = np.array([1., 1., 1., 1., 1.])
        print("\nMa trận A (5x5 - trị riêng thực):")
        output_matrix(A, precision=4)
        run_all_cases(A, x0, precision=7)

        A3 = input_matrix('PWDF_input_A3.txt', convert_fractions=False)
        x0_3 = np.array([-1., 1., 0., 0.])
        print("\n\nMa trận A3 (4x4 - trị riêng phức liên hợp):")
        output_matrix(A3, precision=4)
        run_all_cases(A3, x0_3, precision=7)

    print(f"Đã ghi kết quả vào {output_path}")
