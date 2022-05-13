'''
Alogrithm 6.3 (Asch et al., 2016, Data assimilation)
'''

import numpy as np
import copy
import argparse
import matplotlib.pyplot as plt
from rmse import rmse
from rk4 import rk4

parser = argparse.ArgumentParser(description='Stochastic EnKF')
#  parser.add_argument('--inflation', '-i', type=float, dest='infl', required=True, help='Inflation factor')
parser.add_argument('--method', '-t', type=str, dest='met', required=True, help='Method of perturbe observation')
parser.add_argument('--member', '-m', type=int, dest='mem', required=True, help='Ensemble member size')
#  parser.add_argument('--split', '-s', type=int, dest='splt', required=False, help='M seperating number')
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
M = args.mem         # ensemble memeber size

# data
xtrue = np.load('./Data/true.npy')
xobs = np.load('./Data/obs.npy')

# (1) initial condition
xf = np.array([xtrue[-1] for i in range(M)]) + np.random.normal(0, 1.0, (M, N))

# (2) loop for target time periods
tr = []
err_xa = []
err_obs = []
x_true = []
x_obs = []
x_anal = []
x_bg = []
for i in range(NS):
    # (3) draw a statistically consistent boservation set
    y = np.array([xobs[i+1] for r in range(M)])
    np.random.seed()
    u = np.random.normal(0, 1.0, (M, P))
    yk = y + u
    # (4) compute the ensemble means
    xf_mean = np.mean(xf, axis=0)
    u_mean = np.mean(u, axis=0)
    yf_mean = np.mean(H @ xf.T, axis=1)
    #  and the normalized anomalies (perturbation)
    Xf = np.array([r - xf_mean for r in xf]) / np.sqrt(M - 1.)
    #  Yf = (np.array([r - yf_mean for r in (H @ xf.T).T]) - np.array([r - u_mean for r in u])) / np.sqrt(M - 1)
    Yf = np.zeros((M, P))
    for s in range(M):
        Yf[s] = (xf[s] - u[s] - yf_mean + u_mean) / np.sqrt(M - 1.)
    #  Yf = (np.array([r - xf_mean for r in xf]) - np.array([r - u_mean for r in u])) / np.sqrt(M - 1)
    # (5) compute the gain
    PHT = Xf.T @ (H @ Xf.T).T
    HPHT = H @ Xf.T @ (H @ Xf.T).T
    if args.met == 'm1':
        K = PHT @ np.linalg.inv(HPHT + R)
    if args.met == 'm2':
        RE = np.array([r - u_mean for r in u]) / np.sqrt(M - 1.)
        K = PHT @ np.linalg.pinv(HPHT + RE.T @ RE)
    if args.met == 'm3':
        K = Xf.T @ Yf @ np.linalg.pinv(Yf.T @ Yf)
    # (6) update the ensemble
    xa = xf + (K @ (yk.T - H @ xf.T)).T
    # (7) compute the ensemble forecast
    xf = np.array([rk4(r, N, F, DT) for r in xa])

    #  print(xf, xa, xobs[i], xtrue[i], Pa, Pb, K)
    #  tr.append(np.sqrt(np.trace(Pa)/N))
    #  x_true.append(xtrue[i+1][0])
    #  x_obs.append(xobs[i+1][0])
    #  x_anal.append(xa[0])

    err_xa.append(rmse(np.mean(xa, axis=0).reshape(1, N), xtrue[i+1].reshape(1, N)))
    err_obs.append(rmse(xobs[i+1].reshape(1, N), xtrue[i+1].reshape(1, N)))

#  np.save('./enkf_%s_%d.npy' % (args.meth, delta * 100), err_xa)

#  plt.plot(tr, label='Trace')
delta = 0.00
plt.plot(err_xa, label='Analysis RMSE, Delta = %0.2f' % (delta))
plt.plot(err_obs, label='Observation RMSE')
print(M, delta, np.mean(err_xa[500:]))
#  plt.plot(x_anal[:], label='Analysis')
#  plt.plot(x_obs[:], label='Observation')
#  plt.plot(x_bg[:], label='Forecast')

#  plt.plot(x_true, label='True')
plt.legend(loc=0)
#  plt.savefig('./rmse_%0.2f.png' % (delta))
plt.show()
#  plt.close()
