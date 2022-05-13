import numpy as np


def l96_func(t, x, exf):
    '''
    Define the differential equations

    Args:
        x (array):  variables
        t (-):      deviative

    Return:
        dx (array): deriative equations
    '''
    N = len(x)
    dx = np.zeros(N)

    for i in range(N):
        dx[i] = (x[(i+1) % N] - x[(i-2+N) % N]) * x[(i-1+N) % N] - x[i] + exf

    return dx


def l96_func_ode(x, t, exf):
    '''
    Define the differential equations

    Args:
        x (array):  variables
        t (-):      deviative

    Return:
        dx (array): deriative equations
    '''
    N = len(x)
    dx = np.zeros(N)

    for i in range(N):
        dx[i] = (x[(i+1) % N] - x[(i-2+N) % N]) * x[(i-1+N) % N] - x[i] + exf

    return dx
