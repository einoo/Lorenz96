import numpy as np
import copy
import argparse
from scipy.integrate import solve_ivp, odeint
from rk4 import rk4
from lorenz96_func import l96_func, l96_func_ode

parser = argparse.ArgumentParser(description='Methods for solving ODE functions')
parser.add_argument('--method', '-m', type=str, dest='method', required=True, help='ODE method')
args=parser.parse_args()

N = 40
F = 8.
DT = 0.05
method = args.method

# initial condition
x = np.zeros(N)
for i in range(N):
    x[i] = F
x[int(N/2 - 1)] = F * 1.001
xi = copy.deepcopy(x)

if method == 'ivp':
    # 1440 steps for spin up (1 year)
    NS = 1 * 365 * 4
    x = solve_ivp(l96_func, [0, NS*DT], xi, method='RK45', args=(F, ), t_eval=[NS * DT], max_step=DT, atol=1, rtol=1)
    # 1440 steps for true data (1 year)
    x0 = x.y.reshape(N, )
    t = np.arange(0.0, (NS+1)*DT, DT)
    x = solve_ivp(l96_func, [0, NS*DT], x0, method='RK45', args=(F, ), dense_output=True, t_eval=t, max_step=DT, atol=1, rtol=1)
    np.save('./true_%s.npy' % (method), x.y.reshape(NS+1, N))

if method == 'ode':
    # 1440 steps for spin up (1 year)
    NS = 1 * 365 * 4
    t = np.arange(0.0, (NS+1)*DT, DT)
    x = odeint(l96_func_ode, xi, t, args=(F, ))
    # 1440 steps for true data (1 year)
    x0 = x[-1]
    t = np.arange(0.0, (NS+1)*DT, DT)
    x = odeint(l96_func_ode, x0, t, args=(F, ))
    np.save('./true_%s.npy' % (method), x)

if method == 'rk4':
    # 1440 steps for spin up (1 year)
    NS = 1 * 365 * 4
    xr= xi
    for i in range(NS):
        xr = rk4(xr, N, F, DT)
    # 1440 steps for true data (1 year)
    xt = np.zeros((NS+1, N), dtype=float)
    xt[0] = xr
    for i in range(NS):
        xr = rk4(xr, N, F, DT)
        xt[i+1] = xr
    np.save('./true_%s.npy' % (method), xt)
