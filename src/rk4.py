import numpy as np


def slope(x, n, f):
    '''
    Obtain the slope based on the partial difference equation

    args:
        x (array): Input value at each location
        n  (int):   System size
        f (float): Constant value for the force term

    return:
        dx (array): Slope at each location
    '''
    dx = np.zeros(n)

    for i in range(n):
        dx[i] = (x[(i+1) % n] - x[(i-2+n) % n]) * x[(i-1+n) % n] - x[i] + f

    return dx


def rk4(x, n, f, dt):
    '''
    Obtain the integrated value based on 4th Runge Kutta method

    args:
        x  (array): Input value at each location
        n  (int):   System size
        f  (float): Magnitude of the external forcing
        dt (float): Time step

    returns:
        rx (array): Integrated values based on 4th Runge Kutta method
    '''
    k = np.zeros((5, n))  # weight factors for the slope

    k[1] = slope(x, n, f)

    k[2] = slope(x + dt/2.*k[1], n, f)

    k[3] = slope(x + dt/2.*k[2], n, f)

    k[4] = slope(x + dt * k[3], n, f)

    k[0] = k[1] + 2. * k[2] + 2. * k[3] + k[4]

    rx = x + 1./6. * dt * k[0]

    return rx


def jacobian_l96(x, n, f, dt, e, u):
    '''
    Return the Jacobian matrix of Lorenz 96 model

    args:
        x (array): State variables
        f (float): Constant external force
        dt (float): Time step
        e (float);  Taylor series limit
        I (array): Identify matrix

    returns:
        jx (array): Tagent linear approximation
    '''
    jx = (rk4(x+e*u, n, f, dt) - rk4(x, n, f, dt))/e
    return jx
