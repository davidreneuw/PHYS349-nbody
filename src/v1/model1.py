from datetime import datetime
from matplotlib.animation import FuncAnimation
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import copy

class Model1(object):
    """Model 1: N-Body Simulation"""

    # Class constants
    axisLabel = ['x', 'y', 'z', 'vx', 'vy', 'vz']
    G = 6.67e-11
    M_E = 5.972e24
    M_S = 1.989e30
    AU = 1.5041e11
    V_E = 3.0287e4
    L2 = [1.01*AU, 0, 0]
    L4 = [0.5e11, 0.866e11, 0]
    unitX = [1, 0, 0]
    unitY = [0, 1, 0]
    unitZ = [0, 0, 1]
    seriesTemp = {'time': 0, 'data': []}
    phaseTemp = {'mass': 0, 'ps': [], 'name': ''}

    def __init__(self, system):
        """Initializes the model

        Args:
            system (list | str): List of time series, or string of system name: "Lagrange L2", "Lagrange L4", "Galaxy"
        """
        self.setSystem([{}])
        if type(system) is list:
            self.setSystem(system)
        elif system == "LagrangeL2":
            self.setSystem(self.getLagrangeSystems()[0])
        elif system == "LagrangeL4":
            self.setSystem(self.getLagrangeSystems()[1])
        elif system == "Galaxy":
            self.setSystem(self.getGalaxySystem())
        

    def setSystem(self, system):
        """Sets the system

        Args:
            system (list): List of time series

        Raises:
            TypeError: System must be a list of time series
            TypeError: System must be a list of time series
        """
        if type(system) is not list:
            raise TypeError("System must be a list of time series")
        if type(system[0]) is not dict:
            raise TypeError("System must be a list of time series")
        self.system = system
    
    def getLagrangeSystems(self):
        """Gets the Lagrange L2 and L4 systems

        Returns:
            list: List of time series for L2 and L4
        """
        Earth = copy.deepcopy(self.phaseTemp)
        Earth['mass'], Earth['ps'], Earth['name'] = self.M_E, [self.AU, 0, 0, 0, self.V_E, 0], "Earth"
        Sun = copy.deepcopy(self.phaseTemp)
        Sun['mass'], Sun['ps'], Sun['name'] = self.M_S, [0, 0, 0, 0, 0, 0], "Sun"
        SatelliteL2 = copy.deepcopy(self.phaseTemp)
        SatelliteL2['mass'], SatelliteL2['ps'], SatelliteL2['name'] = 0.001, [*self.L2, 0, 0, 0], "L2"
        SatelliteL4 = copy.deepcopy(self.phaseTemp)
        SatelliteL4['mass'], SatelliteL4['ps'], SatelliteL4['name'] = 0.001, [*self.L4, 0, 0, 0], "L4"
        timeSeriesL2 = copy.deepcopy(self.seriesTemp)
        timeSeriesL2['data'] = [Earth, Sun, SatelliteL2]
        timeSeriesL4 = copy.deepcopy(self.seriesTemp)
        timeSeriesL4['data'] = [Earth, Sun, SatelliteL4]
        return [[timeSeriesL2], [timeSeriesL4]]

    def getGalaxySystem(self):
        """Creates a galaxy system

        Returns:
            list: List of time series for galaxy system
        """
        mainGalaxy = copy.deepcopy(self.phaseTemp)
        mainGalaxy['mass'], mainGalaxy['ps'], mainGalaxy['name'] = 1e11*self.M_S, [0, 0, 0, 0, 0, 0], "Main Galaxy" #TODO
        otherGalaxy = copy.deepcopy(self.phaseTemp)
        otherGalaxy['mass'], otherGalaxy['ps'], otherGalaxy['name'] = 1e10*self.M_S, [3.086e21, 0, 0, 0, 0, 0], "Perturbor Galaxy" #TODO
        starRing1, starRing2, starRing3 = [], [], []
        angles = np.arange(0, 2*np.pi, 2*np.pi/6)
        r1 = 1.543e20
        r2 = 2*r1
        r3 = 3*r1
        for angle in angles:
            star1 = copy.deepcopy(self.phaseTemp)
            star2 = copy.deepcopy(self.phaseTemp)
            star3 = copy.deepcopy(self.phaseTemp)
            star1['mass'], star1['ps'] = 0.001, [r1*np.cos(angle), r1*np.sin(angle), 0, 0, 0, 0] #TODO
            starRing1.append(star1)
            star2['mass'], star2['ps'] = 0.001, [r2*np.cos(angle), r2*np.sin(angle), 0, 0, 0, 0] #TODO
            starRing2.append(star2)
            star3['mass'], star3['ps'] = 0.001, [r3*np.cos(angle), r3*np.sin(angle), 0, 0, 0, 0] #TODO
            starRing3.append(star3)
        timeSeries = copy.deepcopy(self.seriesTemp)
        timeSeries['data'] = [mainGalaxy, otherGalaxy, *starRing1, *starRing2, *starRing3]
        return [timeSeries]
        
        
    def __str__(self):
        """Creates a string representation of the object

        Returns:
            str: String representation of the object
        """
        __out = ""
        for series in self.system:
            __out += "="*100+"\n"
            __out += "Time: "+str(series['time'])+"\n"
            __out += "Data:\n"
            for p in series['data']:
                __out += f"---Mass: {p['mass']}kg   Phase Space: {p['ps']}\n"
        return __out
    
    def dpdt(self, system, t=0):
        """Calculates the derivative of the phase space

        Args:
            system (list): List of time series
            t (int, optional): Time at which to find the derivative. Defaults to 0.

        Returns:
            list: list of phase space derivatives
        """
        derivativeSpace = []
        for i, p in enumerate(system[t]['data']):
            otherP:list = copy.deepcopy(system[t]['data'])
            otherP.pop(i)
            ax, ay, az = 0, 0, 0
            for other in otherP:
                a = self.a(p, other)
                ax += a[0]
                ay += a[1]
                az += a[2]
            derivativeSpace.append([p['ps'][3], p['ps'][4], p['ps'][5], ax, ay, az])
        return derivativeSpace

    def unit_vector(self, vector):
        """ Returns the unit vector of the vector.  
        
        Returns:
            list: Unit vector
        """
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        """ Returns the angle in radians between vectors 'v1' and 'v2'.
        
        Returns:
            float: Angle in radians
        """
        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
                
    def a(self, p1, p2, softening=0.01):
        """Calculates the acceleration between two particles

        Args:
            p1 (dict): Phase space dictionary for particle 1
            p2 (dict): Phase space dictionary for particle 2
            softening (float, optional): Softening parameter. Defaults to 0.01.

        Returns:
            list: Acceleration vector
        """
        dx = p2['ps'][0]-p1['ps'][0]
        dy = p2['ps'][1]-p1['ps'][1]
        dz = p2['ps'][2]-p1['ps'][2]
        vec = [dx, dy, dz]
        alpha, beta, gamma = self.angle_between(vec, Model1.unitX),self.angle_between(vec, Model1.unitY), self.angle_between(vec, Model1.unitZ)
        r = np.sqrt(dx**2+dy**2+dz**2+softening**2)
        norm = self.G*p2['mass']/r**2
        acc = [norm*np.cos(alpha), norm*np.cos(beta), norm*np.cos(gamma)]
        for i, val in enumerate(acc):
            if np.abs(val) < 1e-8:
                acc[i] = 0
        return acc

    def integrate(self, f, t):
        """Integrates the system using the Leapfrog method

        Args:
            f (function): dpdt function
            t (list): List of time values

        Returns:
            list: List of time series for the solved system
        """
        # For Leapfrog, we start by finding the initial values that will help us integrate the system
        h = t[1]-t[0]
        solved = []

        initTime = copy.deepcopy(self.seriesTemp)
        initTime['time'] = t[0]
        h2Time = copy.deepcopy(self.seriesTemp)
        h2Time['time'] = t[0]+h/2
        for j, p in enumerate(self.system[0]['data']):
            initSpace = copy.deepcopy(self.phaseTemp)
            h2Space = copy.deepcopy(self.phaseTemp)
            initSpace['mass'] = p['mass']
            initSpace['ps'] = p['ps']
            initSpace['name'] = p['name']

            x0, y0, z0 = initSpace['ps'][0], initSpace['ps'][1], initSpace['ps'][2]
            vx0, vy0, vz0 = initSpace['ps'][3], initSpace['ps'][4], initSpace['ps'][5]
            initDpdt = f(self.system, t=0)
            ax0, ay0, az0 = initDpdt[j][3], initDpdt[j][4], initDpdt[j][5]
            vxh2, vyh2, vzh2 = vx0 + (h/2)*ax0, vy0 + (h/2)*ay0, vz0 + (h/2)*az0
            xh2, yh2, zh2 = x0 + h*vx0, y0 + h*vy0, z0 + h*vz0

            h2Space['mass'] = p['mass']
            h2Space['ps'] = [xh2, yh2, zh2, vxh2, vyh2, vzh2]
            h2Space['name'] = p['name']
            initTime['data'].append(initSpace)
            h2Time['data'].append(h2Space)
            
        solved.append(initTime)
        solved.append(h2Time)

        # Then from then on we find each new step. The big issue here is that we take half-steps between each time steps, and as such,
        # we must find values for each half-steps, totalling at 2*len(t)-1 half-steps. We start at 2 since we already have the initial and first half-step
        for i in range(2,2*len(t)-1):
            timeSeries = copy.deepcopy(self.seriesTemp)
            timeSeries['time'] = t[i//2]
            for j, p in enumerate(solved[i-1]['data']):
                phaseSpace = copy.deepcopy(self.phaseTemp)
                t2 = solved[i-2]['data'][j]['ps']
                t1 = solved[i-1]['data'][j]['ps']
                a = f(solved, t=i-1)[j][3:]
                vx, vy, vz = t2[3]+h*a[0], t2[4]+h*a[1], t2[5]+h*a[2]
                x, y, z = t2[0]+h*t1[3], t2[1]+h*t1[4], t2[2]+h*t1[5]
                phaseSpace['mass'] = p['mass']
                phaseSpace['ps'] = [x, y, z, vx, vy, vz]
                timeSeries['data'].append(phaseSpace)
            solved.append(timeSeries)
        del solved[1::2] # We remove all odd index elements, which corresponds to all the half-steps
        self.setSystem(solved)
        return solved
    
    def getPos(self, timeSeries, axis=0):
        """Get position series of each atom given an axis.

        Args:
            timeSeries (list): List of time steps dictionaries generated by leapfrog
            axis (int, optional): Axis to use. Defaults to 0.

        Returns:
            list: List of position series for each atom of the given axis.
        """
        pos = np.array([[s['data'][i]['ps'][axis] for i in range(len(s['data']))] for s in timeSeries]).T.tolist()
        return pos

    def spacePlot(self, system, ax1=0, ax2=1, markevery=100):
        """Plots the position of each atom in a 2D space
        
        Args:
            system (list): List of time series for the solved system
            ax1 (int, optional): Axis to use for the x-axis. Defaults to 0.
            ax2 (int, optional): Axis to use for the y-axis. Defaults to 1.
            markevery (int, optional): Number of points to skip between each plotted point. Defaults to 100.
        
        Returns:
            matplotlib.figure.Figure: Figure object
        """
        fig, ax = plt.subplots()
        pos1= self.getPos(system, axis=ax1)
        pos2= self.getPos(system, axis=ax2)
        for i, x in enumerate(pos1):
            ax.plot(x, pos2[i], label=system[0]['data'][i]['name'], markevery=markevery, marker='o', ms=4, ls='-')
        ax.legend()
        ax.set_ylabel(f"{self.axisLabel[ax2]} position in natural units")
        ax.set_xlabel(f"{self.axisLabel[ax1]} position in natural units")
        title = f"{self.axisLabel[ax2]} as a function of {self.axisLabel[ax1]}"
        ax.set_title(title)
        return fig
    
    # Animation function that uses the data generated by the integrate() function
    def animate(self, system, ax1=0, ax2=1, interval=100, save=False, filename=None, markevery=10):
        """Animates the position of each atom in a 2D space
        
        Args:
            system (list): List of time series for the solved system
            ax1 (int, optional): Axis to use for the x-axis. Defaults to 0.
            ax2 (int, optional): Axis to use for the y-axis. Defaults to 1.
            interval (int, optional): Time interval between each frame in milliseconds. Defaults to 100.
            save (bool, optional): True to save the animation as a mp4 file. Defaults to False.
            filename (str, optional): Name of the file to save the animation. Defaults to 'animation.mp4'.
        
        Returns:
            matplotlib.animation.FuncAnimation: Animation object
        """
        if filename is None:
            now = datetime.now().strftime("%d-%m-%Y_%H-%M-%S")
            filename = f'./assets/{now}.mp4'
        else:
            filename = f'./assets/{filename}.mp4'
        pos1= self.getPos(system, axis=ax1)
        pos2= self.getPos(system, axis=ax2)
        fig, ax = plt.subplots()
        ax.set_xlim(np.min(pos1), np.max(pos1))
        ax.set_ylim(np.min(pos2), np.max(pos2))
        ax.set_ylabel(f"{self.axisLabel[ax2]} position")
        ax.set_xlabel(f"{self.axisLabel[ax1]} position")
        title = f"{self.axisLabel[ax2]} as a function of {self.axisLabel[ax1]}"
        ax.set_title(title)
        lines = []
        for i, x in enumerate(pos1):
            line, = ax.plot([], [], label=system[0]['data'][i]['name'], marker='o', ms=4, ls='-', markevery=markevery)
            lines.append(line)
        ax.legend()
        def init():
            for line in lines:
                line.set_data([], [])
            return lines
        def animate(i):
            for j, line in enumerate(lines):
                line.set_data(pos1[j][:i], pos2[j][:i])
            return lines
        anim = FuncAnimation(fig, animate, init_func=init, frames=len(system), interval=interval, blit=True)
        if save:
            anim.save(filename, writer='ffmpeg', fps=30)
        return anim
