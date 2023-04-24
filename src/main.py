import sys
import os
sys.path.append('./src/view')
from view_v1 import View

if __name__ == "__main__":
    os.system('python3 -m unittest ./src/tests/test_v1.py')
    View() 