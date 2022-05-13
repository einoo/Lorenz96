import numpy as np
import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from lorenz96_func import l96_func
from rmse import rmse

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

N = 36
F = 8.
DT = 0.05

# initial condition
x = np.zeros(N)
for i in range(N):
    x[i] = F
x[int(N/2 - 1)] = F * 1.001

# 14400 steps for spin up (10 years)
t = np.arange(0.0, 14401*DT, DT)
x = solve_ivp(l96_func, [0, 14400*DT], x, method='RK45', args=(F, ), t_eval=[14400*DT], max_step=DT, atol=1, rtol=1)
x0 = x.y.reshape(N, )

# investigate the erro propogation
NS = 200 # Steps for each run
M = 250  # Total run numbers based on Lorenz 1996

# generate random normal distributed values
mu, sigma = 0.0, 0.001
noise = np.random.normal(mu, sigma, M)

t = np.arange(0.0, (NS + 1)*DT, DT)
sigma = np.zeros((M, NS))

for m in range(M):
    x0n = x0 + noise[m]
    # update state variables
    x = solve_ivp(l96_func, [0, NS * DT], x0, method='RK45', args=(F, ), dense_output=True, t_eval=t, max_step=DT, atol=1, rtol=1)
    xn = solve_ivp(l96_func, [0, NS * DT], x0n, method='RK45', args=(F, ), dense_outpu=True, t_eval=t, max_step=DT, atol=1, rtol=1)
    # calculate the rmse
    sigma[m] = rmse(xn.y.T, x.y.T)
    # succeed the initial value for the next run
    x0 = x.y.T[-1]
    print(m)

# Obtain the average value
sigma_average = np.mean(sigma, axis=0)

# Calculate the log10 value
log10sigma = np.zeros(NS, dtype=float)
for i in range(NS):
    log10sigma[i] = math.log10(sigma_average[i])

#  print(sigma_average, len(sigma_average))
# One-column figure
fig = plt.figure(constrained_layout=True, figsize=(4.5, 4.5), dpi=300)
ax = fig.add_subplot(111)


ax.plot(t[1:]/DT/4, log10sigma, 'k-', linewidth=0.5)

ax.set_xlim(0, 50)
ax.set_ylim(-5, 1.)

ax.set(xlabel='Days', ylabel='log$_{10} E$')

ax2 = ax.twinx()
color = 'tab:blue'
ax2.set_ylabel('Variations of average prediction error ($E$)', color=color)
ax2.plot(t[1:]/DT/4, sigma_average, 'tab:blue', linewidth=0.5)
ax2.set_ylim(-2, 10)
ax2.tick_params(axis='y', labelcolor=color)

plt.savefig('./Lorenz_1996_Fig2.png')
plt.show()
