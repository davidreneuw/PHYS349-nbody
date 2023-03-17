import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

G = 6.67e-11

class Model(object):
    unitX = [1, 0, 0]
    unitY = [0, 1, 0]
    unitZ = [0, 0, 1]
    seriesTemp = {'time': 0, 'data': []}
    phaseTemp = {'mass': 0, 'ps': []}

    def __init__(self, system) -> None:
        self.system = system
        
    def __str__(self) -> str:
        __out = ""
        for series in self.system:
            __out += "="*100+"\n"
            __out += "Time: "+str(series['time'])+"\n"
            __out += "Data:\n"
            for p in series['data']:
                __out += f"---Mass: {p['mass']}kg   Phase Space: {p['ps']}\n"
        return __out
    
    def dpdt(self, system, t=0):
        derivativeSpace = []
        for i, p in enumerate(self.system['data'][t]):
            otherP:list = self.system['data'][t].copy()
            otherP.pop(i)
            ax, ay, az = 0, 0, 0
            for other in otherP:
                a = self.a(p, other)
                ax += a[0]
                ay += a[1]
                az += a[2]
            derivativeSpace.append([p[3], p[4], p[5], ax, ay, az])

    def unit_vector(self, vector):
        """ Returns the unit vector of the vector.  """
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
                
    def a(self, p1, p2):
        dx = p2['ps'][0]-p1['ps'][0]
        dy = p2['ps'][1]-p1['ps'][1]
        dz = p2['ps'][2]-p1['ps'][2]
        r = np.sqrt(dx**2+dy**2+dz**2)
        return G*p2['mass']/r**2

# nSystem = Model.seriesTemp.copy()
# p1 = Model.phaseTemp.copy()
# p1['mass'] = 10
# p1['ps'] = [0, 0, 0, 1, 0, 0]
# p2 = Model.phaseTemp.copy()
# p2['mass'] = 20
# p2['ps'] = [2, 1, 1, 1, 0, 0]
# nSystem['data'] = [p1, p2]
# model = Model([nSystem])
# print(model)
