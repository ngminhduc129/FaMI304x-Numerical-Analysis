# =============================================================================
# pp3_SVD.py - Phân rã giá trị kỳ dị (Singular Value Decomposition)
#
# Chức năng: Phân tích ma trận A thành A = U * Σ * V^T thông qua
#            trị riêng của A^T A.
#
# Các hàm chính:
#   read_matrix(path)            - đọc ma trận
#   compute_svd(A_df, tol)       -> (U, Sigma, Vt)
#   compute_svd_2(A_df, tol)     -> (U, Sigma, Vt) (biến thể)
#
# Input: Đọc từ file SVD_input_A.txt
# Cách dùng: python pp3_SVD.py
# =============================================================================
from pathlib import Path
import contextlib
import numpy as np
import pandas as pd
from typing import Tuple

__dir__ = Path(__file__).parent.resolve()

pd.set_option('display.precision', 12)  # Increase decimal precision
pd.set_option('display.width', 300)     # Wider display
pd.set_option('display.max_columns', None)  # Show all column

def read_matrix(path: str) -> pd.DataFrame:
    """
    Reads a whitespace-delimited matrix from a text file (allowing decimals
    or integer fractions) and returns it as a pandas DataFrame of floats.
    """
    data = []
    with open(path, 'r') as f:
        for lineno, line in enumerate(f, start=1):
            tokens = line.strip().split()
            if not tokens:
                continue
            row = []
            for tok in tokens:
                try:
                    if tok.count('/') == 1:
                        num, den = tok.split('/')
                        val = int(num) / int(den)
                    else:
                        val = float(tok)
                except Exception as e:
                    raise ValueError(
                        f"Lỗi phân tích token '{tok}' tại dòng {lineno} trong {path}: {e}"
                    )
                row.append(val)
            data.append(row)

    if not data:
        # empty file → empty DataFrame
        return pd.DataFrame()

    # ensure all rows have same length
    ncol = len(data[0])
    if any(len(r) != ncol for r in data):
        lengths = [len(r) for r in data]
        raise ValueError(
            f"Độ dài hàng không khớp trong {path}: mong đợi {ncol} cột, nhận được {lengths}"
        )

    return pd.DataFrame(data)


def print_matrix(df: pd.DataFrame, float_format: str = '{:,.9f}') -> None:
    """
    Nicely prints the DataFrame without row/col labels, formatting floats.
    """
    print(
        df.to_string(
            index=False,
            header=False,
            float_format=float_format.format
        )
    )
def complete_orthonormal_basis(B: np.ndarray, tol: float = 1e-12) -> np.ndarray:
    """
    Given B (n×k) with orthonormal columns, extend to a full n×n orthonormal basis.
    """
    n, k = B.shape
    Q = np.zeros((n, n))
    Q[:, :k] = B
    for j in range(k, n):
        # start with a random vector
        v = np.random.rand(n)
        # orthonormalize against existing columns
        for i in range(j):
            v -= np.dot(Q[:, i], v) * Q[:, i]
        norm = np.linalg.norm(v)
        # if vector is nearly zero, retry
        while norm < tol:
            v = np.random.rand(n)
            for i in range(j):
                v -= np.dot(Q[:, i], v) * Q[:, i]
            norm = np.linalg.norm(v)
        Q[:, j] = v / norm
    return Q
def compute_svd(A_df: pd.DataFrame, tol: float = 1e-10) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute the SVD of A (m×n) via the smaller of A·Aᵀ or Aᵀ·A.
    Prints intermediate steps using pandas formatting.

    Returns:
        U_df: m×m DataFrame of left singular vectors
        Sigma_df: m×n DataFrame of singular values on the diagonal
        Vt_df: n×n DataFrame of right singular vectors transposed
    """
    # Convert to NumPy
    A = A_df.values.astype(float)
    m, n = A.shape

    print("--- Ma trận đầu vào A ---")
    print_matrix(A_df)

    # Choose smaller route
    if True:
        print("Cách: m <= n, dùng M = A · Aᵀ ({}×{})".format(m, m))
        M = A.dot(A.T)
        M_df = pd.DataFrame(M)
        print("Ma trận M:")
        print_matrix(M_df)

        # Eigen-decomposition of M
        eigvals, eigvecs = np.linalg.eigh(M)
        # sort descending
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

        # Print eigenpairs
        for i, lam in enumerate(eigvals, start=1):
            print(f"Trị riêng λ{i} = {lam}")
            print(f"Vector riêng u{i}:")
            print_matrix(pd.DataFrame(eigvecs[:, i-1]))

        # Singular values
        sigmas = np.sqrt(np.clip(eigvals, 0, None))
        print("Giá trị kỳ dị:")
        for i, s in enumerate(sigmas, start=1):
            print(f"σ{i} = {s}")

        # Build Σ
        Sigma = np.zeros((m, n))
        for i in range(min(m, n)):
            Sigma[i, i] = sigmas[i]
        Sigma_df = pd.DataFrame(Sigma)
        print("Ma trận Σ:")
        print_matrix(Sigma_df)

        # Compute right singular vectors V
        V = np.zeros((n, n))
        for i in range(min(m, n)):
            if sigmas[i] > tol:
                v = A.T.dot(eigvecs[:, i]) / sigmas[i]
            else:
                v = np.zeros(n)
            V[:, i] = v
            print(f"Vector kỳ dị phải v{i+1}:")
            print_matrix(pd.DataFrame(v))

        # Complete V to full orthonormal basis
        V_complete = complete_orthonormal_basis(V[:, :m], tol)
        Vt_df = pd.DataFrame(V_complete.T)

        # U is full eigenvector matrix of M
        U_complete = eigvecs
        U_df = pd.DataFrame(U_complete)

    return U_df, Sigma_df, Vt_df

def compute_svd_2(A_df: pd.DataFrame, tol: float = 1e-10) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Compute the SVD of A (m×n) via the smaller of Aᵀ·A.
    Prints intermediate steps using pandas formatting.

    Returns:
        U_df: m×m DataFrame of left singular vectors
        Sigma_df: m×n DataFrame of singular values on the diagonal
        Vt_df: n×n DataFrame of right singular vectors transposed
    """
    # Convert to NumPy
    A = A_df.values.astype(float)
    m, n = A.shape

    print("--- Ma trận đầu vào A ---")
    print_matrix(A_df)

    # Choose smaller route
    if True:
        print("Cách: m > n, dùng N = Aᵀ · A ({}×{})".format(n, n))
        N = A.T.dot(A)
        N_df = pd.DataFrame(N)
        print("Ma trận N:")
        print_matrix(N_df)

        # Eigen-decomposition of N
        eigvals, eigvecs = np.linalg.eigh(N)
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

        for j, lam in enumerate(eigvals, start=1):
            print(f"Trị riêng λ{j} = {lam}")
            print(f"Vector riêng v{j}:")
            print_matrix(pd.DataFrame(eigvecs[:, j-1]))

        sigmas = np.sqrt(np.clip(eigvals, 0, None))
        print("Giá trị kỳ dị:")
        for j, s in enumerate(sigmas, start=1):
            print(f"σ{j} = {s}")

        Sigma = np.zeros((m, n))
        for j in range(min(m, n)):
            Sigma[j, j] = sigmas[j]
        Sigma_df = pd.DataFrame(Sigma)
        print("Ma trận Σ:")
        print_matrix(Sigma_df)

        # Compute left singular vectors U
        U = np.zeros((m, m))
        for j in range(min(m, n)):
            if sigmas[j] > tol:
                u = A.dot(eigvecs[:, j]) / sigmas[j]
            else:
                u = np.zeros(m)
            U[:, j] = u
            print(f"Vector kỳ dị trái u{j+1}:")
            print_matrix(pd.DataFrame(u))

        # Complete U to full orthonormal basis
        U_complete = complete_orthonormal_basis(U[:, :n], tol)
        U_df = pd.DataFrame(U_complete)

        # V is full eigenvector matrix of N
        V_complete = eigvecs
        Vt_df = pd.DataFrame(V_complete.T)

    return U_df, Sigma_df, Vt_df

if __name__ == "__main__":
    output_path = str(__dir__ / "pp3_SVD_result.txt")
    with open(output_path, "w", encoding="utf-8") as f, contextlib.redirect_stdout(f):
        A_df = read_matrix(str(__dir__ / "SVD_input_A.txt"))
        print_matrix(A_df)
        U, S, Vt = compute_svd(A_df, tol=1e-9)
        print("--- U cuối cùng ---")
        print_matrix(U)
        print("--- Σ cuối cùng ---")
        print_matrix(S)
        print("--- V^T cuối cùng ---")
        print_matrix(Vt)
        U, S, Vt = compute_svd_2(A_df, tol=1e-9)
        print("--- U cuối cùng ---")
        print_matrix(U)
        print("--- Σ cuối cùng ---")
        print_matrix(S)
        print("--- V^T cuối cùng ---")
        print_matrix(Vt)
    print(f"Đã ghi kết quả vào {output_path}")




