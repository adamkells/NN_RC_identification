""" Python scripts to compute Markov model from cleaned data files.
 Single column files of numerical values"""

""" Primary outstanding issues with code:
A) Adapt to handle two dimensional data set
B) Parallelize the computing of the data
"""

""" Packages to import """
import numpy as np
from tqdm import tqdm  # package for making a progress bar
import ipdb
import matplotlib.pyplot as plt
import time
from multiprocessing import Pool
#from contextlib import closing
from functools import partial

""" So new plan, I'm going to build my code around two objects. One object contains all of the information about our 
simulation and the other all the information about the Markov model which we want to create """

""" Here is the first object, my data set. It will take 4 inputs. A list of data files, a file containing bias
 positions (all zeros for unbiased MD), a file containing the force constants (all zeros for unbiased MD) and 
 the simulation temperature (default to 300) """

class data_set:
    def __init__(self, data, bias_pos ,force ,T=300):
        self.data = data
        self.bias_pos = bias_pos
        self.force = force
        self.T = T
        if isinstance(self.data, list) == 0:
            raise ValueError('Data does not appear to be of correct dimensions, please check inputs.')

        if isinstance(self.data[0], str):
            print("(Initialized a data set object containing {} simulations)".format(len(self.data)))
        elif isinstance(self.data[0], int):
            print("(Initialized a data set object containing 1 simulation of length {})".format(len(self.data)))

    def plot_data(self):
        plt.figure()
        plt.plot(self.data)
        plt.xlabel('Time')
        plt.ylabel('Position')
        plt.title('Simulation Trajectory: Position vs Time')

    def hist_data(self):
        plt.figure()
        plt.hist(self.data, bins=(max(self.data)))
        plt.xlabel('Position')
        plt.ylabel('Counts')
        #plt.title('Simulation Trajectory: Position vs Time')

""" The second object which I will require is one which contains all of the information of my model. I want my model
object to take the data set as an input. Then by calling some functions on my model object, I can output all of the
nice things that we might want to (free energies, relaxation times, mean first passage times etc etc) """
class cont_markov_model:
    def __init__(self, data_set: list, lag_time: int = 1, num_bins: int = 100, bin_edges: list = [-40, 40]):
        self.lag_time = lag_time
        self.num_bins = num_bins
        self.data_set = data_set
        self.bin_edges = bin_edges
        self.MSM = self.markov_analysis()

        print("(Created a Markov model object at lag time {} consisting of {} trajctories)".format(self.lag_time,
                                                                                                   len(self.data_set)))

    def markov_analysis(self):
        # need to parallelize this

        # define my bin edges
        bin_edges = np.linspace(self.bin_edges[0], self.bin_edges[1], self.num_bins + 1)
        bin_centers = bin_edges[:-1] + 0.5 * (bin_edges[1] - bin_edges[0])

        # sort through my input files and obtain state and transition counts for each simulation
        #sim_counts = np.empty([len(self.data_set.data), self.num_bins])
        #sim_ntr = np.empty([self.num_bins, self.num_bins])


        """ This bit is to do the parallelization. I don't expect this will work the first time round but the main 
        issues that I see potentially occuring are related to the passing of arguments and also the way the output of
        the looped function is formatted. """
        """ I'm also not completely confident on the way I'm outputing the result of the parallel loop"""
        def main():
            iterable = self.data_set.data
            pool = Pool()
            func = partial(loading_data, bin_edges, self.num_bins, self.lag_time)
            counts, ntr = pool.map(func, iterable)
            pool.close()
            pool.join()
            return counts, ntr

        counts, ntr = loading_data(self.data_set.data, bin_edges, self.num_bins, self.lag_time)
        #counts, ntr = main()

        sim_counts = counts
        if len(self.data_set.data) > 1 and isinstance(self.data_set.data[0], str):
            sim_ntr = np.sum(ntr)
        else:
            sim_ntr = ntr

        unvisited = [i for i in range(self.num_bins) if np.max(sim_counts[0][i]) < 1]
        sim_counts = np.delete(sim_counts, unvisited, 1)
        sim_ntr = np.delete(sim_ntr, unvisited, 0)
        sim_ntr = np.delete(sim_ntr, unvisited, 1)
        bin_centers = np.delete(bin_centers, unvisited, 0)
        self.num_bins = len(bin_centers)

        MSM = dham_msm(self.num_bins, bin_centers, sim_ntr, sim_counts, self.data_set.force, self.data_set.bias_pos)
        #need to rewrite this function and then update this line

        return MSM


    def free_energy_profile(self):
        free_energy_profile(self.MSM)

class discrete_markov_model:
    def __init__(self, data_set: list, simulation, lag_time: int = 1):
        self.lag_time = lag_time
        self.data_set = data_set
        self.x = simulation.x
        self.y = simulation.y
        self.num_bins = simulation.states
        self.temp = simulation.temp
        self.MSM = self.markov_analysis()

        print("(Created a Markov model object at lag time {} consisting of {} trajectories)".format(self.lag_time,
                                                                                                   len(self.data_set)))
    def markov_analysis(self):
        # define my bin edges
        bin_centers = self.x

        counts = []
        ntr = []
        bias_pos = []
        force = []
        for data_inst in self.data_set:
            counts_tmp, ntr_tmp = dham_prep(data_inst.data, self.num_bins, self.lag_time)

            counts.append(counts_tmp)
            ntr.append(ntr_tmp)
            bias_pos.append(self.x[data_inst.bias_pos])
            force.append(data_inst.force)

        if len(self.data_set) > 1:
            #ipdb.set_trace()
            ntr = np.sum(ntr,0)

        kbT = 0.001987 *self.temp

        MSM = dham_msm(self.num_bins, bin_centers, ntr, counts, force, bias_pos, kbT)

        return MSM


    def free_energy_profile(self):
        free_energy_profile(self.MSM, self.temp)

def loading_data(input_file, bin_edges, num_bins, lag_time):
    time_series = bin_coordinate(input_file, bin_edges, num_bins)
    counts, ntr = dham_prep(time_series, num_bins, lag_time)
    return counts, ntr

def bin_coordinate(input_name, bin_edges, num_bins):

    """
    :param input_name: this can come in two forms, a list of numbers (a trajectory) or the name of a text file in which
    case the text file is loaded and the trajectory is loaded from there
    :param bin_edges: this is a list of linearly increasing numbers defining the bin edges. The data is sorted in to the
    bins defined
    :return: this will give the input time series back as a time series of discrete integers (the original trajectory
    sorted in to the defined bins).
    """
    # check what kind of input we have and import it
    if isinstance(input_name[0], str) == 1:
        time_series = []
        # open the colvar file and attempt to read it line by line
        with open(input_name[0], 'r') as f:
            for line in f:
                try:
                    time_series.append(float(line.strip('\r\n')))  #time_series is a list with my trajectory
                except:
                    pass
    else:
        time_series = input_name

    values_range = max(bin_edges) - min(bin_edges)

    # cycle through the time series obtained from colvar and place in the appropriate bin
    # also uses tqdm to include a progress bar
    for i in tqdm(range(len(time_series))):
        time_series[i] = np.ceil(((time_series[i]-bin_edges[0])*(num_bins)/values_range))
        if time_series[i] > num_bins:
            time_series[i] = num_bins
        if time_series[i] == 0:
            time_series[i] = 1

    time_series=[int(x) for x in time_series]

    return time_series

""" function to construct a MSM from simulation data using DHAM """
def dham_prep(data, num_bins, lag_time):

    """
    :param data: This is a single trajectory
    :param num_bins: the number of bins
    :param lag_time: the lag time of the model
    :return:
    """
    bins_sim = range(num_bins+1)
    bins_sim = [ (x+0.5) for x in bins_sim]

    counts=[]
    #ipdb.set_trace()
    hist, bin_edges = np.histogram(data,bins=bins_sim)
    counts.append(hist)

    # next we need to count transition numbers
    ntr = np.zeros([num_bins, num_bins])
    for j in range(len(data)-lag_time-1):
        ntr[data[j]-1,data[j+lag_time]-1]+=1

    return counts, ntr


""" This function needs a massive overhaul, particularly to account for multiple simulation trajectories"""


def dham_msm(num_bins, bin_centers, ntr, counts, force, bias_pos, kbT):
    #ipdb.set_trace()
    MSM = np.zeros([num_bins, num_bins])
    for i in range(num_bins):
        for j in range(num_bins):
            if ntr[i][j] > 0:
                msm_tmp = 0
                for sim_num, counts_k in enumerate(counts):
                    u=0.5*force[sim_num]*((bias_pos[sim_num]-bin_centers)**2)
                    if counts_k[0][i]>0:
                         msm_tmp=msm_tmp+counts_k[0][i]*np.exp(-(u[j]-u[i])/2/kbT)
                MSM[i][j]=ntr[i][j]/msm_tmp

    row_sums = MSM.sum(axis=1)
    MSM = MSM / row_sums[:, np.newaxis]

    return MSM


def spectral_analysis(matrix):

    """
    :param matrix: Takes a square matrix as an input
    :return: Gives back the ordered eigenvectors and eigenvalues (sorted by ascending order of eigenvalue)
    """

    eigenValues, eigenVectors = np.linalg.eig(matrix)
    idx = eigenValues.argsort()  # idx is the way we order it
    eigenValues = eigenValues[idx]  # Order eigenvalues in ascending order
    eigenVectors = eigenVectors[:, idx]  # Order eigenvectors
    return eigenValues, eigenVectors
    
def iterative_hummer_szabo_method(MSM,N):

    I_N=np.identity(N)
    eq = free_energy_profile(MSM)
    # for loop running over
    trial_clusterings=[[np.linspace(0,N-1)]]
    while trial_clusterings[-1][0]<num_bins-N:

        for i in range(N):
            a=sum(eq[1:5])
    """ This code is unfinished"""

def free_energy_profile(MSM,temp,image=0):

    """
    :param MSM: Takes a Markov state model as an input (an N by N matrix describing transition probabilities)
    :return: G (the free energy profile of the MSM) and equilibrium (the equilibrium probabilities of the MSM states)
    """
    kbT = 0.001987 * temp
    val, vec = spectral_analysis(MSM.T)
    val = [x for x in val]
    maxpos = val.index(max(val))
    equilibrium=vec[:,maxpos]
    equilibrium=equilibrium/sum(equilibrium)
    G = -kbT*np.log(equilibrium)
    G = G - np.min(G)

    if image == 1:
        fig1, ax1 = plt.subplots()
        ax1.plot(G)
        fig1.savefig('test.png')

    return G , equilibrium


############################################################################################################
"""“Persons attempting to find a motive in this narrative will be prosecuted; persons attempting to find a 
moral in it will be banished; persons attempting to find a plot in it will be shot.
BY ORDER OF THE AUTHOR
G.G., CHIEF OF ORDNANCE”
― Mark Twain, The Adventures of Huck Finn """

""" In other words, the code below this line is not general, so be careful using it. The generalisation of these 
functions will be a longer term project."""


""" FUNCTIONS TO GET MEAN FIRST PASSAGE TIMES"""

data_length = 100000 # number of time steps to use in computing MFPT by direct counting

""" This function uses a method by Jensen et al. to calculate the MFPT between regions"""
def mfpt_jensen(MM,state1,state2):
    state1.append(state2)
    #MM = MM[state1,state1]
    M_til = MM
    M_til = np.delete(M_til,state2, axis=0)
    M_til = np.delete(M_til, state2, axis=1)


    L = MM[state2,:]
    L = np.delete(L, state2, axis=0)

    O = np.ones([1,1])

    G, eq = free_energy_profile(MM)

    V = np.zeros([num_bins,1])
    for i in state1:
        V[i]=eq[i]/np.sum(eq[state1])
    V = np.delete(V,state2,axis=0)

    eigval, eigvec = spectral_analysis(M_til)
    a = np.matmul(np.linalg.inv(eigvec),np.reshape(V,[num_bins-1,1]))

    t = 0

    for ii in range(len(a)):
        t += (a[ii]/(1 - eigval[ii]) ** 2) * O * np.matmul(L, eigvec[:, ii])

    return t

"""This function aims to calculate the mfpt by directly counting from the trajectory"""

#state 1 is the starting region, state 2 is the region being entered and state 3 is the region we want to not enter
# so this will give the mean first passage time to enter 2 from 1 without first entering 3.
def mfpt_state_two_periodic_regions(state1, state2, state3,data_length):
    if not isinstance(state1, list):
        state1 = [state1]
    if not isinstance(state2, list):
        state2 = [state2]

    # initialise lists
    mfpt = []
    num_of_mfpt = []
    tmp_data = data
    for traj in tmp_data:

        traj = traj[0:data_length] #take subset of total data

        tmp = 0
        # tmp will be zero as long as there is another instance of state 1 in the future
        while tmp == 0:
            #ipdb.set_trace()
            t = time.time()
            # look for the index of the next occurrence of each of the states
            next_1 = next((i for i, v in enumerate(traj) if v in list(state1)),-1)
            next_2 = next((i for i, v in enumerate(traj) if v in list(state2) and i > next_1),-1)
            next_3 = next((i for i, v in enumerate(traj[:next_2]) if v in list(state3) and i > next_1 and i < next_2),-1)
            elapsed_1 = time.time() - t

            """ compares the indices from above to add to a list of mfpt events and then deletes the first however many
            frames of the trajectory since these have been examined"""
            if next_3 == -1:
                if next_2 != -1:
                    mfpt.append(0.5 * (next_2 - next_1)*(next_2 - next_1 + 1)) # result of sum from i = 1 to N of i
                    num_of_mfpt.append(next_2 - next_1)
                    del traj[:next_2+1]
            else:
                del traj[:next_3+1]

            if next_1 == -1:
                tmp = 1
            if next_2 == -1:
                tmp = 1
            elapsed_2 = time.time() - elapsed_1 - t

    # get the average by deleting the sum of all mfpt events by the number of events counted.
    mfpt_value = np.sum(mfpt)/np.sum(num_of_mfpt)

    return mfpt_value

def meyer_mfpt(MM,state1,state2):

    A = np.identity(num_bins) - MM
    A[state2,:]=0
    A[state2,state2] = 1
    b = np.ones([num_bins,1])
    b[state2] = 0

    MFP = np.linalg.solve(A, b)

    G, eq = free_energy_profile(MM.T)
    eq = np.real(eq)
    output = np.dot(eq[state1],MFP[state1])/np.sum(eq[state1])

    return output