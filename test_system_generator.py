""" This file will contain the necessary functions to generate the rate/markov matrices
for some simple test systems and also to generate trajectories (double well, triple well,
szabo-berezhovski 2-D potential) """

from msm_scripts import *
import numpy as np
from scipy.linalg import expm
import ipdb

# i want to make a simulation object that contains all the information of the chosen simulation
# and then have a function i can call on it to create a data object to pass to my msm scripts

class simulation:
    def __init__(self, potential, sim_length, dt, biasing_protocol, num_of_states, T, initial_state):
        self.pot = potential
        self.length = sim_length
        self.bias = biasing_protocol
        self.states = num_of_states
        self.temp = T
        self.initial = initial_state
        self.dt = dt
        if not biasing_protocol[1] in range(num_of_states):
            raise ValueError('Biasing position is not a valid state index. Make sure these values are integers from 1 to num_of_states')

    def generate_data(self, num_of_traj=1):

        """ next i would need some functions to a) create a selection of different potentials of interest b) define a
        markov chain monte carlo transition step c) generate the specified number of transitions according to a
        particular biasing protocol"""
        x ,y = potential(self.pot, self.states)
        traj = trajectory(x, y, self.length, self.dt, self.bias, self.temp,self.initial)
        simulation_data = data_set(traj, bias_pos=self.bias[1], force=self.bias[0], T=self.temp)
        self.x = x
        self.y = y
        return simulation_data


def potential(pot_type,num_of_states):

    if pot_type.lower() == 'single_well':
        x = np.linspace(-np.pi, np.pi, num_of_states)
        y = -1*np.sin(0.5*x + 0.5*np.pi)

    elif pot_type.lower() == 'double_well':
        x = np.linspace(-np.pi, np.pi, num_of_states)
        y = 1*np.sin(2*(x - 3*np.pi / 4))

    elif pot_type.lower() == 'triple_well':
        x = np.linspace(-np.pi, np.pi, num_of_states)
        y = 1*np.sin(3*(x - 7*np.pi / 2))

    elif pot_type.lower() == 'flat':
        x = np.linspace(-np.pi, np.pi, num_of_states)
        y = 0*x


    else:
        raise ValueError('invalid input potential type')

    y = y-np.min(y)

    return x, y

def trajectory(x, y, sim_length, dt, biasing_protocol, T, initial_state):

    # first construct the Markov matrix for the process
    v_bias = 0.5*biasing_protocol[0]*(x - x[biasing_protocol[1]])**2
    y = y + v_bias
    A = 10 # Arrhenius factor
    KbT = 0.001987 * T

    K = np.zeros([len(x),len(x)])
    for i in range(len(x)-1):
        K[i,i+1] = A* np.exp(-(y[i+1]-y[i])/KbT/2)
        K[i+1, i] = A * np.exp(-(y[i] - y[i+1]) / KbT/2)

    for i in range(len(x)):
        K[i, i] = -np.sum(K[i,:])

    mm = expm(K*dt)

    G, eq = free_energy_profile(mm, T)

    # next we want to run a simulation of the Markov model produced above
    if initial_state != 0:
        traj = [initial_state]
    else:
        traj = [np.random.choice(np.arange(1, len(x)+1), p=np.real(eq))]

    count = 0
    #ipdb.set_trace()
    while count < (sim_length-1):
        rand = np.random.uniform(0, 1)
        for j in range(len(x)):
            if np.sum(mm[traj[count]-1, 0:j+1]) > rand:
                traj.append(j+1)
                break

        count += 1
    #ipdb.set_trace()
    #traj += 1
    return traj



