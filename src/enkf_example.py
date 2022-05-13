import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

N = 2
M = 10
P = 1
TN = 2500
TS = 25    # observation time step
H = np.array([1, 0])

w = 3.5e-2
v = 3.0e-4

def harmonic_oscillator(x):
    MP = np.array([[2.0 + w**2 - v**2 * x[0]**2, -1.], [1.0, 0.0]])
    x = MP @ x
    return x

def jacobian(x, e, u):
    jx = (harmonic_oscillator(x+e*u) - harmonic_oscillator(x))/e
    return jx

def tml(x):
    TML = np.array([[2.0 + w**2 - 3. * v**2 * x[0]**2, -1.], [1.0, 0.0]])
    return TML

# Initial condition
x = np.array([1., 0.])

# Generate the truth / nature run
xtrue = []
t = []
xtrue.append(x[0])
t.append(0)
# Generate the observation / add noise (0, 10)
xobs = np.zeros(int(TN/TS)+1)
tobs = np.zeros(int(TN/TS)+1)
noise = np.random.normal(0, 10, int(TN/TS)+1)
for i in range(TN + 1):
    t.append(i+1)
    x = harmonic_oscillator(x)
    xtrue.append(x[0])
    if (i % TS == 0):
        tobs[int(i/TS)] = i
        xobs[int(i/TS)] = x[0] + noise[int(i/TS)]

# data assimilation by EKF
E = 0.5
I = np.identity(N)
xf = np.array([xtrue[-1], xtrue[-2]])
Pb = np.identity(N) * 25.5
R = np.identity(P) * 10.
xekf = []
xekf.append(xf[0])
for i in range(TN + 1):
    # Having observation, conduct the data assimilation
    xf = harmonic_oscillator(xf)
    if i > 0 and i % TS == 0:
        # (3) compute the gain
        H = H.reshape(P, N)
        K = Pb @ H.T @ np.linalg.inv(H @ Pb @ H.T + R)
        # (4) compute the state analysis
        xa = xf + K @ (xobs[int(i/TS)] - H @ xf)
        # (5) compute the analysis error covariance matrix
        Pa = (np.identity(N) - K @ H) @ Pb
        # (6) compute the forecast state
        xf = harmonic_oscillator(xa)
        #  # (7) comput the forecast/background error covariance matrix
        Mk = tml(xa)
        Pb = Mk @ Pa @ Mk.T
        #  JM = np.zeros((N, N))
        #  for j in range(N):
        #      JM[:, j] = jacobian(xa, E, I[:, j])
        #  Pb = JM @ Pa @ JM.T
        print(xobs[int(i/TS)], xtrue[i], xf, Pb, K)
    xekf.append(xf[0])

fig = plt.figure(constrained_layout=True, figsize=(4.5, 3.0), dpi=300)
ax = fig.add_subplot(111)

ax.plot(t, xtrue, 'k-', linewidth=0.5, label='True')
ax.plot(tobs, xobs, 'ko', fillstyle='none', ms=3, markeredgewidth=0.8, label='Observations')
ax.plot(t, xekf, 'b--', linewidth=0.5, label='EKF')
#  ax.plot(t, xtrue, 'mo--', fillstyle='none', ms=5, markeredgewidth=0.8, linewidth=0.5, label='RMSE (analysis-true)')
ax.set_xlim(500, 2500)
ax.set_ylim(-200, 200)
ax.set(xlabel='Time', ylabel='')
plt.legend(loc=0)
plt.grid(b=True, which='both', color='gray', linestyle=':', lw=0.5)

#  plt.savefig('./3DAVR.png')
plt.show()
