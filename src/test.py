import numpy as np
import copy
import argparse
from scipy.integrate import solve_ivp, odeint
from lorenz96_func import l96_func, l96_func_ode
from rk4 import rk4

parser = argparse.ArgumentParser(description='Methods for solving ODE functions')
parser.add_argument('--method', '-m', type=str, dest='method', required=True, help='ODE method')
args=parser.parse_args()

N = 40
F = 8.
DT = 0.05
method = args.method

# initial condition
x = np.zeros(N, dtype=float)
for i in range(N):
    x[i] = F
x[int(N/2 - 1)] = F * 1.001
x0 = copy.deepcopy(x)

NS = 40

for i in range(NS):
    x = rk4(x, N, F, DT)
print(x)

t = np.arange(0.0, (NS+1)*DT, DT)
xa = solve_ivp(l96_func, [0, NS*DT], x0, method='RK45', args=(F, ), dense_output=True, max_step=DT)
print(xa.y.T[-1])

xc = solve_ivp(l96_func, [0, NS*DT], x0, method='RK45', args=(F, ), dense_output=True, max_step=DT, atol=10, rtol=10)
print(xc.y.T[-1])

xb = odeint(l96_func_ode, x0, t, args=(F, ))
print(xb[-1])
breakpoint()
