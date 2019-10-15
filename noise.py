import numpy as np

def white_noise(mu=0, rms=1.0, t_size):

	return sigma*np.random.randn(t_size)+mu

def one_f_noise(mu,rms,t_size):

    return 1/(white_noise(mu,rms,t_size).sort())

	