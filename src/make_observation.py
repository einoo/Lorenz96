import numpy as np

true = np.load('./Data/true.npy')

nx, ny = true.shape

# generate random noise
mu, sigma = 0.0, 1.0
noise = np.random.normal(mu, sigma, nx * ny).reshape(nx, ny)

true += noise

np.save('./obs.npy', true)
