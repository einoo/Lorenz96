'''
Alogrithm x.x (Asch et al., 2016, Data assimilation)
'''

import numpy as np
import argparse
from rk4 import rk4
from rmse import rmse

parser = argparse.ArgumentParser(description='Extensive chaos in the Lorenz-96 model')
parser.add_argument('--inflation', '-i', type=float, dest='infl', required=True, help='Inflation factor')
args=parser.parse_args()

N = 40     # number of state variables
F = 8.0    # external force
P = 40     # number of observations
DT = 0.05  # time interval
NS = 1460  # total time steps (one year)

# requirments
R = np.identity(P)   # observation error covariance matrix
# Q = 0.  ignore the model error covariance matrix
H = np.identity(P)   # tangent linear observation model
# M = lorenz96_func  forward model

# data
xtrue = np.load('./Data/true.npy')
xobs = np.load('./Data/obs.npy')

# (1) initial condition
xf = xobs[0]
Pb = np.identity(N) * args.infl
t = np.array([0.0, DT])

# (2) loop for target time periods
err_xa = []
for i in range(NS):
    # (3) compute the gain
    K = Pb @ H.T @ np.linalg.inv(H @ Pb @ H.T + R)
    # (4) compute the state analysis
    xa = xf + K @ (xobs[i+1] - H @ xf)
    # (5) compute the forecast state
    #  x = solve_ivp(l96_func, [0, DT], xa, method='RK45', args=(F, ), t_eval=[DT])
    #  xf = x.y.reshape(N, )
    xf = rk4(xa, N, F, DT)

    err_xa.append(rmse(xa.reshape(1, N), xtrue[i+1].reshape(1, N)))

print(args.infl, np.mean(err_xa[-500:]))
with open('./3dvar.txt', 'a') as f:
    f.write('%0.2f %0.4f \n' % (args.infl, np.mean(err_xa[-500:])))
