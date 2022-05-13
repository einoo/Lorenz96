import numpy as np

def l96_tlm(x, dt):
    n = len(x)
    tlm = np.zeros((n, n))
    for i in range(n):
        tlm[i, i] = 1. - dt
        tlm[i, (i+1) % n] = x[(i-1+n) % n] * dt
        tlm[i, (i-1+n) % n] = (x[(i+1) % n] - x[(i-2+n) % n]) * dt
        tlm[i, (i-2+n) % n] = -x[(i-1+n) % n] * dt
    return tlm
