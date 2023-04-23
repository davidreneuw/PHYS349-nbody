import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('./src/v1')
from model1 import Model1 as m1
import assertpy as ap


# generate test for model 1
def test1():
    system = [
        {'time': 0, 'data': [
            {'mass': 1, 'ps': [1, 0, 0, 0, 0, 0]},
            {'mass': 1, 'ps': [-1, 0, 0, 0, 0, 0]}
        ]}
    ]
    model = m1(system)
    test = model.dpdt(system)
    ap.assert_that(test).is_equal_to([[0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0]])
test1()