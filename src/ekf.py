'''
Alogrithm 6.1 (Asch et al., 2016, Data assimilation)
'''

import numpy as np
import copy
import argparse
from lorenz96_tlm import l96_tlm
import matplotlib.pyplot as plt
from rmse import rmse
from rk4 import rk4, jacobian_l96

parser = argparse.ArgumentParser(description='Extensive chaos in the Lorenz-96 model')
parser.add_argument('--inflation', '-i', type=float, dest='infl', required=True, help='Inflation factor')
parser.add_argument('--method', '-m', type=str, dest='meth', required=True, help='M construction method')
parser.add_argument('--split', '-s', type=int, dest='splt', required=False, help='M seperating number')
#  parser.add_argument('--force', '-f', type=int, dest='force', required=True, help='External force')
args=parser.parse_args()

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

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
# DM = lorenz96_tlm  tangent linear forward model

# data
xtrue = np.load('./Data/true.npy')
xobs = np.load('./Data/obs.npy')

# (1) initial condition
xf = xobs[0]
Pb = np.identity(N) * 25.0
t = np.array([0.0, DT])

# test of tangent linear model
#  Mk = l96_tlm(xf, DT)
#  print(Mk[0])
#  # build by numerical way
#  e = 1.e-10
#  I = np.identity(N)
#  JM = (odeint(l96_func, xf + e * I[0], t, args=(F, ))[-1] - odeint(l96_func, xf, t, args=(F, ))[-1])/e
#  print(JM)
#  breakpoint()

# (2) loop for target time periods
delta = args.infl
E = 1.e-5
I = np.identity(N)
tr = []
err_xa = []
err_obs = []
x_true = []
x_obs = []
x_anal = []
x_bg = []
for i in range(NS):
    x_bg.append(xf[0])
    # (3) compute the gain
    K = Pb @ H.T @ np.linalg.inv(H @ Pb @ H.T + R)
    #  print(np.max(K), np.min(K))
    # (4) compute the state analysis
    xa = xf + K @ (xobs[i+1] - H @ xf)
    # (5) compute the analysis error covariance matrix
    Pa = (np.identity(N) - K @ H) @ Pb
    # (6) compute the forecast state
    xf = rk4(xa, N, F, DT)
    # (7) comput the forecast/background error covariance matrix
    if args.meth == 'm1':
        #  Method 1
        JM = np.zeros((N, N))
        for j in range(N):
            JM[:, j] = jacobian_l96(xa, N, F, DT, E, I[:, j])
        Pb = (1.0 + delta) * JM @ Pa @ JM.T
    if args.meth == 'm2':
        #  Method 2
        Mk = l96_tlm(xa, DT)
        Pb = (1.0 + delta) * Mk @ Pa @ Mk.T
    if args.meth == 'm3':
        #  Method 3 (4D VAR)
        xav = copy.deepcopy(xa)
        SPLIT = args.splt
        MV = np.zeros((SPLIT, N, N))
        for k in range(SPLIT):
            MV[k] = l96_tlm(xav, DT/SPLIT)
            xav = rk4(xav, N, F, DT/SPLIT)
        if SPLIT == 1:
            M4V = MV[0]
        else:
            M4V = MV[SPLIT-1]
            for s in range(SPLIT-2, -1, -1):
                M4V = M4V @ MV[s]
        #  M4V = MV[4] @ MV[3] @ MV[2] @ MV[1] @ MV[0]
        Pb = (1.0 + delta) * M4V @ Pa @ M4V.T

    #  print(xf, xa, xobs[i], xtrue[i], Pa, Pb, K)
    #  breakpoint()
    tr.append(np.sqrt(np.trace(Pa)/N))
    x_true.append(xtrue[i+1][0])
    x_obs.append(xobs[i+1][0])
    x_anal.append(xa[0])

    err_xa.append(rmse(xa.reshape(1, N), xtrue[i+1].reshape(1, N)))
    err_obs.append(rmse(xobs[i+1].reshape(1, N), xtrue[i+1].reshape(1, N)))

np.save('./ekf_%s_%d.npy' % (args.meth, delta * 100), err_xa)

#  plt.plot(tr, label='Trace')
#  plt.plot(err_xa, label='Analysis RMSE, Delta = %0.2f' % (delta))
#  plt.plot(err_obs, label='Observation RMSE')
print(delta, np.mean(err_xa[500:]))
#  plt.plot(x_anal[:], label='Analysis')
#  plt.plot(x_obs[:], label='Observation')
#  plt.plot(x_bg[:], label='Forecast')

#  plt.plot(x_true, label='True')
#  plt.legend(loc=0)
#  plt.savefig('./rmse_%0.2f.png' % (delta))
#  plt.show()
#  plt.close()
