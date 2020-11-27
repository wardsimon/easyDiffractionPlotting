__author__ = 'github.com/wardsimon'
__version__ = '0.0.1'

from abc import abstractmethod


class BasePlotter:
    def __init__(self, plotter, obj):
        self.plotter = plotter
        self.obj = obj
        self.actors = []

    @abstractmethod
    def plot(self, *args, **kwargs):
        pass
