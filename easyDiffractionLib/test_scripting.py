__author__ = 'github.com/wardsimon'
__version__ = '0.0.1'

from plotter import winPlotter
from easyDiffractionLib import Phase

p = Phase('test')
p.add_atom('Cu', 'Cu2+', 1, 0, 0, 0)
w = winPlotter(p)
w.plot()
