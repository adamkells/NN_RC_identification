from unittest import TestCase
from msm_scripts import *
#####
class TestMarkov_model(TestCase):
    a=['colvar_1b']
    b=[0]
    c=[0]
    d=300
    ala5_data = data_set(a, b, c, d)
    ala5_model = markov_model(ala5_data, bin_edges=[-180,180])
    ala5_model.markov_analysis()
    pass
