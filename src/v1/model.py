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
    
    # Test 2
                
    def a(self, p1, p2):
        dx = p2['ps'][0]-p1['ps'][0]
        dy = p2['ps'][1]-p1['ps'][1]
        dz = p2['ps'][2]-p1['ps'][2]
        vec = [dx, dy, dz]
        alpha = self.angle_between(vec, Model.unitX)
        beta = self.angle_between(vec, Model.unitY)
        gamma = self.angle_between(vec, Model.unitZ)
        r = np.sqrt(dx**2+dy**2+dz**2)
        norm = G*p2['mass']/r**2
        acc = [norm*np.cos(alpha), norm*np.cos(beta), norm*np.cos(gamma)]
        return acc
    
    

    def integrate(f, init, t, h, method="Leapfrog"):
        if method == 'Leapfrog':
            # For Leapfrog, we start by finding the initial values that will help us integrate the system
            h = t[1]-t[0]
            v0 = solved[0][1]
            x0 = solved[0][0]
            # We find the first Euler half-step
            vh2 = v0 + (h/2)*self.a(th=x0)
            xh2 = x0 + h*vh2
            solved.append([xh2, vh2]) # We add that half-step to our system
            # Then from then on we find each new step. The big issue here is that we take half-steps between each time steps, and as such,
            # we must find values for each half-steps, totalling at 2*len(t)-1 half-steps. We start at 2 since we already have the initial and first half-step
            for i in range(2,2*len(t)-1):
                # In our MS Teams, Prof. Hudson mentioned that position and velocity evaluated at different times could be stored in the same row.
                # We find v and x, then append them to our system as being at the same time.
                v = solved[i-2][1]+h*self.a(solved[i-1][0])
                x = solved[i-2][0]+h*solved[i-1][1]
                solved.append([x,v])
            del solved[1::2] # We remove all odd index elements, which correspond to all the half-steps
            self.state = solved[-1]


nSystem = Model.seriesTemp.copy()
p1 = Model.phaseTemp.copy()
p1['mass'] = 10
p1['ps'] = [0, 0, 0, 1, 0, 0]
p2 = Model.phaseTemp.copy()
p2['mass'] = 6e24
p2['ps'] = [2, 1, 1, 1, 0, 0]
nSystem['data'] = [p1, p2]
model = Model([nSystem])
model.a(p1,p2)

