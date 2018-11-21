""" This file will contain python scripts which can take in a dataset object and return the bias position and
force for a new iteration of umbrella potentials.
The functions should have a range of loss functions to test out.
Future directions should include more flexible biasing potentials"""

from msm_scripts import *
from test_system_generator import *

def linear_regression_enhanced_sampling(dataset):

    data = []
    forces = []
    bias_pos = []

    # want to get the existing data and associated biases
    for item in data_set:
        data.append(item.data)
        forces.append(item.force)
        bias_pos.append(item.bias_pos)

    

