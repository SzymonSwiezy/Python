import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import PillowWriter
from scipy.integrate import odeint
from scipy.integrate import solve_ivp


m3 = 4
m2 = 1
m1 = 1

v1 = 0.9911981217
v2 = 0.7119472124

x1_0 = -1
y1_0 = 0
x2_0 = 1
y2_0 = 0
x3_0 = 0
y3_0 = 0

vx1_0 = v1
vy1_0 = v2  # v = sqrt(a*r) with a=m2 (since G=1, r12=1)
vx2_0 = v1
vy2_0 = v2
vx3_0 = -2*v1/m3
vy3_0 = -2*v2/m3


def dSdt(t, S):
    x1, y1, x2, y2, x3, y3, vx1, vy1, vx2, vy2, vx3, vy3 = S
    r12 = np.sqrt((x2-x1)**2 + (y2-y1)**2)
    r13 = np.sqrt((x3-x1)**2 + (y3-y1)**2)
    r23 = np.sqrt((x2-x3)**2 + (y2-y3)**2)

    return [vx1, vy1, vx2, vy2, vx3, vy3,
            m2/r12**3 * (x2-x1) + m3/r13**3 * (x3-x1),
            m2/r12**3 * (y2-y1) + m3/r13**3 * (y3-y1),
            m1/r12**3 * (x1-x2) + m3/r23**3 * (x3-x2),
            m1/r12**3 * (y1-y2) + m3/r23**3 * (y3-y2),
            m1/r13**3 * (x1-x3) + m2/r23**3 * (x2-x3),
            m1/r13**3 * (y1-y3) + m2/r23**3 * (y2-y3),
            ]


t = np.linspace(0, 20, 1000)
sol = solve_ivp(dSdt, (0, 20), y0=[x1_0, y1_0, x2_0, y2_0, x3_0, y3_0,
                                   vx1_0, vy1_0, vx2_0, vy2_0, vx3_0, vy3_0], method='DOP853', t_eval=t, rtol=1e-10, atol=1e-13)

t = sol.t

x1 = sol.y[0]
y1 = sol.y[1]
x2 = sol.y[2]
y2 = sol.y[3]
x3 = sol.y[4]
y3 = sol.y[5]
plt.plot(t, x1)
plt.show()


tt = 1/np.sqrt(6.67e-11*1.99e30/(1.5e11)**3)  # to secound
tt = tt/(3600*24*365.25)*np.diff(t)[0]  # to years


def animate(i):
    ln1.set_data([x1[i], x2[i], x3[i]], [y1[i], y2[i], y3[i]])
    text.set_text('Time = {:.1f} Years'.format(i*tt))


fig, ax = plt.subplots(1, 1, figsize=(8, 8))
ax.grid()
ln1, = plt.plot([], [], 'ro', lw=3, markersize=6)
text = plt.text(0, 1.75, 'asdasd', fontsize=20,
                backgroundcolor='white', ha='center')
ax.set_ylim(-2, 2)
ax.set_xlim(-2, 2)
ani = animation.FuncAnimation(fig, animate, frames=1000, interval=50)
ani.save('plan.gif', writer='pillow', fps=30)
