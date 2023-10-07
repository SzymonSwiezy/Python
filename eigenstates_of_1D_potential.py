import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal

N = 2000
dy = 1/N
y = np.linspace(0, 1, N+1)


# def mL2V(y):
#     return 1000*(y-1/2)**2
def mL2V(y):
    return 1000*np.exp(-(y-0.7)**2 / (2*0.05**2)) #gaussian 


V = mL2V(y)
main_diagonal = 1/dy**2+mL2V(y)[1:-1]
off_diagonal = -1/(2*dy**2)*np.ones(len(main_diagonal)-1)
w, v = eigh_tridiagonal(main_diagonal, off_diagonal)

plt.plot(v.T[0]**2)
plt.plot(v.T[1]**2)
plt.plot(v.T[2]**2)
plt.plot(v.T[3]**2)


plt.show()
