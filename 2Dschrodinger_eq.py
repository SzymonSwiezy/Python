from scipy import sparse
import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse.linalg import eigs
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import PillowWriter
# plt.style.use(['science', 'notebook'])

N = 150
X, Y = np.meshgrid(np.linspace(0, 1, N, dtype=float),
                   np.linspace(0, 1, N, dtype=float))


def get_potential(x, y):
    return 0*X


V = get_potential(X, Y)
diag = np.ones([N])
diags = np.array([diag, -2*diag, diag])
D = sparse.spdiags(diags, np.array([-1, 0, 1]), N, N)
T = -1/2*sparse.kronsum(D, D)
U = sparse.diags(V.reshape(N**2), (0))
H = T+U
eigenvalues, eigenvectors = eigsh(H, k=10, which='SM')


def get_e(n):
    return eigenvectors.T[n].reshape((N, N))


plt.figure(figsize=(9, 9))
plt.contourf(X, Y, get_e(7)**2, 20)
plt.show()
