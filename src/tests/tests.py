import numpy as np
import copy
import matplotlib.pyplot as plt
import sys
sys.path.append('./src/v1')
from model1 import Model1
import assertpy as ap
import unittest
from datetime import datetime

# Constants
tempSystem = [
    {'time': 0, 'data': [
        {'mass': 0, 'ps': [0, 0, 0, 0, 0, 0]},
        {'mass': 0, 'ps': [0, 0, 0, 0, 0, 0]}
    ]}
]

Me = copy.deepcopy(Model1.phaseTemp)
Me['mass'] = 65
Me['ps'] = [6_378_000, 0, 0, 0, 0, 0]
Earth = copy.deepcopy(Model1.phaseTemp)
Earth['mass'] = 5.92e24
Earth['ps'] = [0, 0, 0, 0, 0, 0]
OtherEarth = copy.deepcopy(Model1.phaseTemp)
OtherEarth['mass'] = 5.92e24
OtherEarth['ps'] = [2*6_378_000, 0, 0, 0, 0, 0]

class TestAcceleration(unittest.TestCase):
    """Tests relating to the acceleration function and the derivative function  """

    def test_magnitude(self):
        """Tests the magnitude of the acceleration of the me on the Earth
        """
        system = copy.deepcopy(Model1.seriesTemp)
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
        """ Tests the orientation of the acceleration of 2 particles aligned in the x axis"""
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
    """Tests relating to the Leapfrog integration method"""

    def test_no_acc(self):
        """Tests the integration of a particle with no acceleration"""
        system = copy.deepcopy(tempSystem)
        system[0]['data'][0]['ps'] = [0, 0, 0, 1, 0, 0]
        system[0]['data'].pop()
        system[0]['data'][0]['mass'] = 10
        model = Model1(system)
        t = np.arange(0, 10.1, 0.1)
        solved = model.integrate(model.dpdt, t)
        endPos = solved[-1]['data'][0]['ps'][0]
        self.assertAlmostEqual(endPos, 10, delta=0.1)

    def test_acc(self):
        """Tests the integration of a particle with acceleration"""
        system = copy.deepcopy(Model1.seriesTemp)
        system['data'] = [Me, Earth]
        model = Model1([system])
        t = np.arange(0, 10.1, 0.1)
        solved = model.integrate(model.dpdt, t)
        startPos = solved[0]['data'][0]['ps'][0]
        endPos = solved[-1]['data'][0]['ps'][0]
        a = model.a(Me, Earth)[0]
        kinematicAnswer = 0.5 * a * 10**2
        self.assertAlmostEqual(endPos-startPos, kinematicAnswer, delta=0.5)

    def test_conservation_energy(self):
        """Tests the conservation of energy in the system"""
        system = copy.deepcopy(Model1.seriesTemp)
        system['data'] = [Me, Earth]
        h_init = Me['ps'][0]
        model = Model1([system])
        a1 = model.a(Me, Earth)[0]
        a2 = model.a(Earth, Me)[0]
        initEnergy = Me['mass']*np.abs(a1)*h_init
        t = np.arange(0, 30, 0.1)
        solved = model.integrate(model.dpdt, t)
        v_Me = np.sqrt(solved[-1]['data'][0]['ps'][3]**2 + solved[-1]['data'][0]['ps'][4]**2 + solved[-1]['data'][0]['ps'][5]**2)
        v_Earth = np.sqrt(solved[-1]['data'][1]['ps'][3]**2 + solved[-1]['data'][1]['ps'][4]**2 + solved[-1]['data'][1]['ps'][5]**2)
        h_endMe = np.sqrt(solved[-1]['data'][0]['ps'][0]**2 + solved[-1]['data'][0]['ps'][1]**2 + solved[-1]['data'][0]['ps'][2]**2)
        h_endEarth = np.sqrt(solved[-1]['data'][1]['ps'][0]**2 + solved[-1]['data'][1]['ps'][1]**2 + solved[-1]['data'][1]['ps'][2]**2)
        finalEnergy = Me['mass']*np.abs(a1)*h_endMe + Earth['mass']*np.abs(a2)*h_endEarth + 0.5*Me['mass']*v_Me**2 + 0.5*Earth['mass']*v_Earth**2
        self.assertAlmostEqual(initEnergy, finalEnergy, delta=0.001*initEnergy)

    def test_3body(self):
        """Tests the integration of 3 particles"""
        system = copy.deepcopy(Model1.seriesTemp)
        system['data'] = [Earth, Me, OtherEarth]
        model = Model1([system])
        t = np.arange(0, 5, 0.1)
        solved = model.integrate(model.dpdt, t)
        dx_Me, dy_Me, dz_Me = solved[-1]['data'][1]['ps'][0]-solved[0]['data'][1]['ps'][0], solved[-1]['data'][1]['ps'][1]-solved[0]['data'][1]['ps'][1], solved[-1]['data'][1]['ps'][2]-solved[0]['data'][1]['ps'][2]
        dx_Earth, dy_Earth, dz_Earth = solved[-1]['data'][0]['ps'][0]-solved[0]['data'][0]['ps'][0], solved[-1]['data'][0]['ps'][1]-solved[0]['data'][0]['ps'][1], solved[-1]['data'][0]['ps'][2]-solved[0]['data'][0]['ps'][2]
        dx_OtherEarth, dy_OtherEarth, dz_OtherEarth = solved[-1]['data'][2]['ps'][0]-solved[0]['data'][2]['ps'][0], solved[-1]['data'][2]['ps'][1]-solved[0]['data'][2]['ps'][1], solved[-1]['data'][2]['ps'][2]-solved[0]['data'][2]['ps'][2]
        disp_Me = np.sqrt(dx_Me**2 + dy_Me**2 + dz_Me**2)
        disp_Earth = np.sqrt(dx_Earth**2 + dy_Earth**2 + dz_Earth**2)
        disp_OtherEarth = np.sqrt(dx_OtherEarth**2 + dy_OtherEarth**2 + dz_OtherEarth**2)
        self.assertAlmostEqual(disp_Me, 0, delta=0.1)
        self.assertGreater(disp_Earth, 10)
        self.assertGreater(disp_OtherEarth, 10)



if __name__ == '__main__':
    now = str(datetime.now()).replace(" ", "_").replace(":", "_").split(".")[0]
    log_file = f'./logs/{str(now)}.txt'
    with open(log_file, "w") as f:
       runner = unittest.TextTestRunner(f, verbosity=2)
       unittest.main(testRunner=runner, verbosity=2)
