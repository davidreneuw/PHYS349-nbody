import tkinter as tk
from tkinter import ttk
import os
import sys
sys.path.append('./src/tests')
from test_v1 import Tests
import copy

class View(object):

    sceneTemp = {
        "title": "N-Body Gravitational Simulation - David René",
        "data": [],
        "size": (600, 200)
    }

    def __init__(self):
        self.root = tk.Tk()
        self.history = []
        self.width = 600
        self.height = 200
        self.sceneLst = []
        self.root.geometry(f"{self.width}x{self.height}")

        btn_back = ttk.Button(self.root, text="Back", width=20, command=lambda: self.back())

        lbl_select = ttk.Label(self.root, text="Please select a case study:")
        btn_lagrange = ttk.Button(self.root, text="Lagrange points", width=20, command=lambda: self.Lagrange())
        btn_galaxy = ttk.Button(self.root, text="Galaxy collision", width=20, command=self.Galaxy)
        btn_test = ttk.Button(self.root, text="Run tests", width=20, command=self.Tests)
        menu_scene = copy.deepcopy(self.sceneTemp)
        menu_scene['data'] = [lbl_select, btn_lagrange, btn_galaxy, btn_test]
        self.sceneLst.append(menu_scene)

        frame_lagrange = tk.Frame(self.root)
        lbl_lagrange = ttk.Label(self.root, text="Which Lagrange point would you like to observe?")
        ttk.Button(frame_lagrange, text="L2", width=20, command=self.L2).grid(row=0, column=0)
        ttk.Button(frame_lagrange, text="L4", width=20, command=self.L4).grid(row=0, column=1)
        lagrange_scene1 = copy.deepcopy(self.sceneTemp)
        lagrange_scene1['data'] = [lbl_lagrange, frame_lagrange, btn_back]
        lagrange_scene1['title'] = "Lagrange points - David René"
        self.sceneLst.append(lagrange_scene1)

        frame_L2 = tk.Frame(self.root)
        tk.Canvas(frame_L2, width=600, height=400).pack(anchor=tk.CENTER, expand=True)
        L2_scene = copy.deepcopy(self.sceneTemp)
        L2_scene['data'] = [frame_L2, btn_back]
        L2_scene['title'] = "L2 Observation - David René"
        L2_scene['size'] = (600, 500)
        self.sceneLst.append(L2_scene)

        frame_L4 = tk.Frame(self.root)
        tk.Canvas(frame_L4, width=600, height=400).pack(anchor=tk.CENTER, expand=True)
        L4_scene = copy.deepcopy(self.sceneTemp)
        L4_scene['data'] = [frame_L4, btn_back]
        L4_scene['title'] = "L4 Observation - David René"
        L4_scene['size'] = (600, 500)
        self.sceneLst.append(L4_scene)

        frame_galaxy = tk.Frame(self.root)
        tk.Canvas(frame_galaxy, width=600, height=400).pack(anchor=tk.CENTER, expand=True)
        galaxy_scene = copy.deepcopy(self.sceneTemp)
        galaxy_scene['data'] = [frame_galaxy, btn_back]
        galaxy_scene['title'] = "Galaxy Collision - David René"
        galaxy_scene['size'] = (600, 500)
        self.sceneLst.append(galaxy_scene)

        self.build(menu_scene)
        self.root.mainloop()

    def Lagrange(self):
        self.build(self.sceneLst[1])

    def L2(self):
        self.build(self.sceneLst[2])

    def L4(self):
        self.build(self.sceneLst[3])
    
    def Galaxy(self):
        self.build(self.sceneLst[4])
    
    def Tests(self):
        os.system('python3 -m unittest ./src/tests/test_v1.py')

    def setGeo(self, size):
        self.width = size[0]
        self.height = size[1]
        self.root.geometry(f"{self.width}x{self.height}")

    def build(self, endScene, back=False):
        if len(self.history) == 0:
            startScene = None
        else:
            startScene = self.history[-1]
        if not back:
            self.history.append(endScene)
        if startScene is not None:
            for widget in startScene['data']:
                widget.pack_forget()
        self.setGeo(endScene['size'])
        self.root.title(endScene['title'])
        for widget in endScene['data']:
            widget.pack(anchor=tk.CENTER, expand=True)
        
        return endScene

    def back(self):
        prev_scene = self.history[-2]
        self.build(prev_scene, back=True)
        self.history.pop(-1)

    
if __name__ == "__main__":
    view = View()





