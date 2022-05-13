import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description='Plot EKF figures')
parser.add_argument('--method', '-m', type=str, dest='meth', required=True, help='M construction method')
args=parser.parse_args()

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

# Plot the evolution of RMSE with different inflation factors
fig = plt.figure(constrained_layout=True, figsize=(4.5, 3.0), dpi=300)
ax = fig.add_subplot(111)

for i in [0.0, 0.03, 0.05, 0.10, 0.20]:
    data = np.load('./ekf_%s_%d.npy' % (args.meth, i * 100))
    x = []
    for j in range(len(data)):
        x.append(j/4.0)
    ax.plot(x, data, '-', linewidth=1.0, label='delta=%0.2f' % (i))
ax.set_xlim(250, 300)
ax.set_ylim(0.0, 1.0)
ax.set(xlabel='Time (day)', ylabel='RMSE (analysis-true)')
plt.grid(b=True, which='both', color='gray', linestyle=':', lw=0.5)
plt.legend(loc=0)

plt.savefig('./ekf_%s.png' % (args.meth))
plt.show()
