from scipy.linalg import eigh
import numpy as np

def solve_double_well(lam=1.0, gamma=1.4, N=4000, L=6.0):
    x = np.linspace(-L, L, N)
    dx = x[1] - x[0]

    # Kinetic matrix: tridiagonal finite difference for -1/2 * d2/dx2
    diag     =  np.ones(N)   / dx**2          # main diagonal
    off_diag = -np.ones(N-1) / (2.0 * dx**2)  # one off from main
    T = np.diag(diag) + np.diag(off_diag, 1) + np.diag(off_diag, -1)

    # Potential matrix: diagonal
    V = np.diag(lam * (x**2 - gamma**2)**2)

    # Full Hamiltonian
    H = T + V

    # Diagonalise — only ask for the lowest 4 eigenvalues/vectors
    # eigh is for symmetric/Hermitian matrices (H is real symmetric here)
    vals, vecs = eigh(H, subset_by_index=[0, 3])

    return x, vals, vecs

x, vals, vecs = solve_double_well()
print(f"E0 = {vals[0]:.8f}")
print(f"E1 = {vals[1]:.8f}")
print(f"ΔE = {vals[1]-vals[0]:.8f}")