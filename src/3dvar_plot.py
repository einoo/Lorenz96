import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'

data = pd.read_csv('./3dvar.txt', index_col=None, header=None, sep=r'\s+')

# Plot the sensitivity of the Bii to RMSE
fig = plt.figure(constrained_layout=True, figsize=(4.5, 3.0), dpi=300)
ax = fig.add_subplot(111)

ax.plot(data.loc[:, 0], data.loc[:, 1], 'mo--', fillstyle='none', ms=5, markeredgewidth=0.8, linewidth=0.5, label='RMSE (analysis-true)')
ax.set_xlim(0, 0.6)
ax.set_ylim(0, 2.5)
ax.set(xlabel='Average B', ylabel='RMSE (forecast-true)')
plt.legend(loc=0)

plt.savefig('./3DAVR.png')
plt.show()
