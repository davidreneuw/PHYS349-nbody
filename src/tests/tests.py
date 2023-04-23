import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
sys.path.append('./src/v1')
from model1 import Model1
import assertpy as ap
import unittest

tempSystem = [
    {'time': 0, 'data': [
        {'mass': 0, 'ps': [0, 0, 0, 0, 0, 0]},
        {'mass': 0, 'ps': [0, 0, 0, 0, 0, 0]}
    ]}
]

class TestAcceleration(unittest.TestCase):

    def test_magnitude(self):
        system = copy.deepcopy(Model1.seriesTemp)
        Me = copy.deepcopy(Model1.phaseTemp)
        Me['mass'] = 65
        Me['ps'] = [6_378_000, 0, 0, 0, 0, 0]
        Earth = copy.deepcopy(Model1.phaseTemp)
        Earth['mass'] = 5.92e24
        Earth['ps'] = [0, 0, 0, 0, 0, 0]
        system['data'] = [Me, Earth]
        model = Model1([system])
        test = model.dpdt([system])
        self.assertAlmostEqual(np.abs(test[0][3]), 9.8, delta=0.2)
        self.assertAlmostEqual(np.abs(test[1][3]), 0, delta=0.01)
        self.assertAlmostEqual(np.abs(test[0][4]), 0, delta=0.01)
        self.assertAlmostEqual(np.abs(test[1][4]), 0, delta=0.01)
        self.assertAlmostEqual(np.abs(test[0][5]), 0, delta=0.01)
        self.assertAlmostEqual(np.abs(test[1][5]), 0, delta=0.01)

    def test_orientation(self):
        system = copy.deepcopy(tempSystem)
        system[0]['data'][0]['ps'] = [100, 0, 0, 0, 0, 0]
        system[0]['data'][1]['ps'] = [-100, 0, 0, 0, 0, 0]
        system[0]['data'][0]['mass'] = 1e10
        system[0]['data'][1]['mass'] = 1e10
        model = Model1(system)
        test = model.dpdt(system)
        self.assertLessEqual(test[0][3], 0)
        self.assertGreaterEqual(test[1][3], 0)
        self.assertAlmostEqual(test[0][4], 0)
        self.assertAlmostEqual(test[1][4], 0)
        self.assertAlmostEqual(test[0][5], 0)
        self.assertAlmostEqual(test[1][5], 0)

class TestIntegration(unittest.TestCase):
    
        def test_integration(self):
            system = copy.deepcopy(tempSystem)
            system[0]['data'][0]['ps'] = [100, 0, 0, 0, 0, 0]
            system[0]['data'][1]['ps'] = [-100, 0, 0, 0, 0, 0]
            system[0]['data'][0]['mass'] = 1e10
            system[0]['data'][1]['mass'] = 1e10
            model = Model1(system)
            test = model.integrate(system, 1, 1)
            self.assertAlmostEqual(test[0][0]['ps'][0], 100, delta=0.01)
            self.assertAlmostEqual(test[0][1]['ps'][0], -100, delta=0.01)

if __name__ == '__main__':
    unittest.main()
