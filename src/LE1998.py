import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from lorenz96_func import l96_func

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

N = 40
F = 8.
DT = 0.05

# initial condition
x = np.zeros(N)
for i in range(N):
    x[i] = F
x[int(N/2 - 1)] = F * 1.001

# time sequence
t = np.arange(0.0, 0.45, DT)

# update state variables
x = solve_ivp(l96_func, [0, 0.4], x, method='RK45', args=(F, ), dense_output=True, t_eval=t, max_step=DT, rtol=1, atol=1)

# plot Fig. 1
fig = plt.figure(constrained_layout=True, figsize=(4.5, 4.5), dpi=300)
ax = fig.add_subplot(111)

n_site = np.arange(1, N+1)
nx, ny = x.y.shape
state = np.zeros((ny, N), dtype=float)
for i in range(ny):
    state[i] = x.y[:, i] - F - i * 0.05
    ax.plot(n_site, state[i], 'k-', lw=0.5)

ax.set_xlim(0, 40)
ax.set_ylim(-0.44, 0.01)

ax.set_yticks([0.0, -0.2, -0.4])
ax.set_yticklabels(['0', '1', '2'])

ax.set_xlabel('site')
ax.set_ylabel('time (days)')

plt.savefig('./Lorenz_and_Emauel_1998_Fig1.png')
plt.show()
